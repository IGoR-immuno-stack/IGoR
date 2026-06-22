#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Core/Aligner.h>
#ifdef IGOR_WITH_SEQAN2
#  include <igor/Core/SeqAn2Aligner.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Catch::Matchers::WithinRel;

namespace {

struct ExpectedAlignment
{
    std::string gene_name;
    std::string cigar;
    double score;
};

struct ActualAlignment
{
    std::string gene_name;
    std::string cigar;
    double score;
};

Matrix<double> build_test_score_matrix(double match_score = 2.0, double mismatch_score = -1.0)
{
    std::vector<double> values(15 * 15, mismatch_score);
    for (int i = 0; i < 15; ++i) {
        values[i + 15 * i] = match_score;
    }
    return Matrix<double>(15, 15, values);
}

Aligner make_legacy_aligner(const Matrix<double> &matrix, int gap_penalty, Gene_class gene_class,
                            const std::vector<std::pair<std::string, std::string>> &genomic_templates)
{
    Aligner aligner(matrix, gap_penalty, gene_class);
    aligner.set_genomic_sequences(genomic_templates);
    return aligner;
}

#ifdef IGOR_WITH_SEQAN2
std::unique_ptr<SeqAn2Aligner>
make_seqan2_aligner(const Matrix<double> &matrix, int gap_penalty, Gene_class gene_class,
                    const std::vector<std::pair<std::string, std::string>> &genomic_templates)
{
    auto aligner = std::make_unique<SeqAn2Aligner>(matrix, gap_penalty, gene_class, AlignmentPreset::with_band(20));
    aligner->set_genomic_sequences(genomic_templates);
    return aligner;
}
#endif

std::string find_germline_sequence(const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const std::string &gene_name)
{
    for (const auto &entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry.second;
        }
    }
    return std::string();
}

std::pair<std::string, std::string>
find_genomic_template(const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                      const std::string &gene_name)
{
    for (const auto &entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry;
        }
    }
    return std::pair<std::string, std::string>{ };
}

std::string cigar_visual(const std::string &cigar, const std::string &read, const std::string &germline)
{
    std::string read_gapped;
    std::string match_line;
    std::string germ_gapped;
    std::size_t read_i = 0;
    std::size_t germ_i = 0;

    auto push_pair = [&](char rc, char gc, char op) {
        read_gapped.push_back(rc);
        germ_gapped.push_back(gc);
        if (rc == '-' || gc == '-') {
            match_line.push_back(' ');
        } else if (op == '=') {
            match_line.push_back('|');
        } else if (op == 'X') {
            match_line.push_back('*');
        } else {
            match_line.push_back((rc == gc) ? '|' : '*');
        }
    };

    for (const auto &entry : parse_cigar(cigar)) {
        const int count = entry.first;
        const char op = entry.second;
        for (int i = 0; i < count; ++i) {
            switch (op) {
            case '=':
            case 'X':
            case 'M': {
                const char rc = (read_i < read.size()) ? read[read_i] : '?';
                const char gc = (germ_i < germline.size()) ? germline[germ_i] : '?';
                push_pair(rc, gc, op);
                ++read_i;
                ++germ_i;
                break;
            }
            case 'I':
            case 'S': {
                const char rc = (read_i < read.size()) ? read[read_i] : '?';
                push_pair(rc, '-', op);
                ++read_i;
                break;
            }
            case 'D':
            case 'N': {
                const char gc = (germ_i < germline.size()) ? germline[germ_i] : '?';
                push_pair('-', gc, op);
                ++germ_i;
                break;
            }
            case 'H':
            case 'P':
            default:
                break;
            }
        }
    }

    std::ostringstream out;
    const std::size_t width = 100;
    for (std::size_t i = 0; i < read_gapped.size(); i += width) {
        out << "\n      read " << read_gapped.substr(i, width) << "\n           " << match_line.substr(i, width)
            << "\n      germ " << germ_gapped.substr(i, width);
    }
    return out.str();
}

std::vector<ActualAlignment> normalize_actuals(const std::forward_list<Alignment_data> &alignments)
{
    std::vector<ActualAlignment> normalized;
    for (const auto &aln : alignments) {
        normalized.push_back({ aln.gene_name, alignment_data_to_cigar(aln), aln.score });
    }
    std::sort(normalized.begin(), normalized.end(), [](const ActualAlignment &a, const ActualAlignment &b) {
        return std::tie(a.gene_name, a.cigar, a.score) < std::tie(b.gene_name, b.cigar, b.score);
    });
    return normalized;
}

std::vector<ExpectedAlignment> normalize_expected(std::vector<ExpectedAlignment> expected)
{
    std::sort(expected.begin(), expected.end(), [](const ExpectedAlignment &a, const ExpectedAlignment &b) {
        return std::tie(a.gene_name, a.cigar, a.score) < std::tie(b.gene_name, b.cigar, b.score);
    });
    return expected;
}

ActualAlignment best_alignment(const std::forward_list<Alignment_data> &alignments)
{
    auto best_it = std::max_element(alignments.begin(), alignments.end(),
                                    [](const Alignment_data &a, const Alignment_data &b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());
    return { best_it->gene_name, alignment_data_to_cigar(*best_it), best_it->score };
}

std::string add_cigar_visual_info(const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  const std::string &label, const std::string &gene_name, const std::string &cigar)
{
    const std::string germline = find_germline_sequence(genomic_templates, gene_name);
    if (germline.empty()) {
        std::ostringstream msg;
        msg << label << " gene=" << gene_name << " has no matching germline sequence for visual explanation";
        return msg.str();
    }
    std::ostringstream msg;
    msg << label << " gene=" << gene_name << " cigar=" << cigar << cigar_visual(cigar, query, germline);
    return msg.str();
}

void assert_alignment_matches(const Alignment_data &alignment, const std::string &query,
                              const std::pair<std::string, std::string> &genomic_template,
                              const ExpectedAlignment &expected)
{
    const std::string actual_cigar =
            alignment_data_to_cigar_full_span(alignment, query.length(), genomic_template.second.length());

    INFO("expected: gene=" << expected.gene_name << " cigar=" << expected.cigar << " score=" << expected.score);
    INFO("actual: gene=" << alignment.gene_name << " cigar=" << actual_cigar << " score=" << alignment.score);

    const std::vector<std::pair<std::string, std::string>> templates = { genomic_template };
    INFO(add_cigar_visual_info(query, templates, "expected", expected.gene_name, expected.cigar));
    INFO(add_cigar_visual_info(query, templates, "actual", alignment.gene_name, actual_cigar));

    REQUIRE(alignment.gene_name == expected.gene_name);
    REQUIRE(actual_cigar == expected.cigar);
    REQUIRE_THAT(alignment.score, WithinRel(expected.score));
}

void assert_best_alignment_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                   const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const ExpectedAlignment &expected)
{
    REQUIRE(!alignments.empty());
    auto best_it = std::max_element(alignments.begin(), alignments.end(),
                                    [](const Alignment_data &a, const Alignment_data &b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());

    auto genomic_template = find_genomic_template(genomic_templates, expected.gene_name);
    if (genomic_template.first.empty()) {
        genomic_template = find_genomic_template(genomic_templates, best_it->gene_name);
    }
    assert_alignment_matches(*best_it, query, genomic_template, expected);
}

void assert_alignment_set_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  std::vector<ExpectedAlignment> expected)
{
    std::vector<const Alignment_data *> actual;
    for (const auto &aln : alignments) {
        actual.push_back(&aln);
    }
    std::sort(actual.begin(), actual.end(), [](const Alignment_data *a, const Alignment_data *b) {
        return std::make_tuple(a->gene_name, alignment_data_to_cigar(*a), a->score)
                < std::make_tuple(b->gene_name, alignment_data_to_cigar(*b), b->score);
    });
    expected = normalize_expected(std::move(expected));

    INFO("actual count=" << actual.size() << " expected count=" << expected.size());
    REQUIRE(actual.size() == expected.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
        INFO("alignment index=" << i);
        auto genomic_template = find_genomic_template(genomic_templates, expected[i].gene_name);
        if (genomic_template.first.empty()) {
            genomic_template = find_genomic_template(genomic_templates, actual[i]->gene_name);
        }
        assert_alignment_matches(*actual[i], query, genomic_template, expected[i]);
    }
}

} // namespace

TEST_CASE("Aligner V gene multiplex dummy alignments.", "[aligner][V_gene][dummy]")
{

    /*
    * A nucleotides for NT of germline of interest
    * T nucleotides for other NT in the read, or to make a germline variation allowing to constrain alignment.
    * G nucleotides for mismatches
    * C nucleotides for penalized gaps
    * 
    * Score costs are set with prime numbers such that total score mismatch is easier to debug. 
    */
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const int gap_penalty = 13;
    std::string query_read;
    std::string germline_ref;
    std::string expected_cigar;
    double expected_score;

    SECTION("complete alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_cigar = "20=";
        expected_score = 140;
    }

    SECTION("short read, free leading reference deletion")
    {
        query_read = std::string(5, 'A') + std::string(1, 'T');
        germline_ref = std::string(9, 'A') + std::string(1, 'T');
        expected_cigar = "4D6=";
        expected_score = 42;
    }

    SECTION("partial overlap read, no 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T');
        expected_cigar = "13D7=13I";
        expected_score = 49;
    }

    SECTION("partial overlap read, single free 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T') + std::string(1, 'A');
        expected_cigar = "13D7=1D13I";
        expected_score = 49;
    }

    SECTION("partial overlap read, three free 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T') + std::string(3, 'A');
        expected_cigar = "13D7=3D13I";
        expected_score = 49;
    }

    SECTION("free trailing germline treated as deletion")
    {
        query_read = std::string(14, 'T');
        germline_ref = std::string(6, 'A') + std::string(14, 'T') + std::string(2, 'A');
        expected_cigar = "6D14=2D";
        expected_score = 98;
    }

    SECTION("mismatch in long match")
    {
        query_read = std::string(20, 'A');
        query_read[4] = 'G';
        germline_ref = std::string(20, 'A');
        expected_cigar = "4=1X15=";
        expected_score = 122;
    }

    SECTION("mismatch in first NT")
    {
        query_read = std::string(20, 'A');
        query_read[0] = 'G';
        germline_ref = std::string(20, 'A');
        expected_cigar = "1X19=";
        expected_score = 122;
    }

    SECTION("mismatch in last NT seen as free gap")
    {
        query_read = std::string(20, 'A');
        query_read[19] = 'G';
        germline_ref = std::string(20, 'A');
        expected_cigar = "19=1D1I";
        expected_score = 133;
    }

    SECTION("penalized insertion")
    {
        query_read = std::string(5, 'A') + "C" + std::string(15, 'A');
        germline_ref = std::string(20, 'A');
        expected_cigar = "5=1I15=";
        expected_score = 127;
    }

    SECTION("penalized deletion")
    {
        query_read = "T" + std::string(17, 'A') + "T";
        germline_ref = "T" + std::string(7, 'A') + "T" + std::string(10, 'A') + "T";
        expected_cigar = "8=1D11=";
        expected_score = 120;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_cigar, expected_score });
}
 
TEST_CASE("Aligner D gene dummy local alignments.", "[aligner][D_gene][dummy]")
{

    /*
    * A nucleotides for matching core of germline/read.
    * T nucleotides for flanking context to shape local overlap.
    * G nucleotides for mismatches.
    * C nucleotides for penalized indels in the core.
    *
    * Scores are prime-valued for easier manual debugging.
    */
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const int gap_penalty = 13;
    std::string query_read;
    std::string germline_ref;
    std::string expected_cigar;
    double expected_score;

    SECTION("complete local alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_cigar = "20=";
        expected_score = 140;
    }

    SECTION("short local core in long read")
    {
        query_read = std::string(4, 'T') + std::string(6, 'A') + std::string(4, 'T');
        germline_ref = std::string(6, 'A');
        expected_cigar = "4I6=4I";
        expected_score = 42;
    }

    SECTION("long local core in long read")
    {
        query_read = std::string(2, 'T') + std::string(16, 'A') + std::string(2, 'T');
        germline_ref = std::string(16, 'A');
        expected_cigar = "2I16=2I";
        expected_score = 112;
    }

    SECTION("mismatch in several places")
    {
        query_read = "AAAGAAAAGAAAAGAA";
        germline_ref = std::string(16, 'A');
        expected_cigar = "3=1X4=1X4=1X2=";
        expected_score = 58;
    }

    SECTION("left mismatch shortens alignment")
    {
        query_read = "AGAGAAAAGAAAAGAA";
        germline_ref = std::string(16, 'A');
        expected_cigar = "4D4I4=1X4=1X2=";
        expected_score = 48;
    }

    SECTION("right mismatch shortens alignment")
    {
        query_read = "AAAGAAAAGAAAGGAA";
        germline_ref = std::string(16, 'A');
        expected_cigar = "3=1X4=1X3=4D4I";
        expected_score = 48;
    }

    SECTION("local overlap with one trailing germline base")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(1, 'T');
        expected_cigar = "8=1D";
        expected_score = 56;
    }

    SECTION("local overlap with three trailing germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(3, 'T');
        expected_cigar = "8=3D";
        expected_score = 56;
    }

    SECTION("local overlap with five trailing germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(5, 'T');
        expected_cigar = "8=5D";
        expected_score = 56;
    }

    SECTION("local overlap with one leading germline base")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(1, 'T') + std::string(8, 'A');
        expected_cigar = "1D8=";
        expected_score = 56;
    }

    SECTION("local overlap with three leading germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(3, 'T') + std::string(8, 'A');
        expected_cigar = "3D8=";
        expected_score = 56;
    }

    SECTION("local overlap with five leading germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(5, 'T') + std::string(8, 'A');
        expected_cigar = "5D8=";
        expected_score = 56;
    }

    SECTION("penalized insertion in core")
    {
        query_read = std::string(5, 'T') + std::string(5, 'A') + "C" + std::string(5, 'T');
        germline_ref = std::string(5, 'T') + std::string(5, 'A') + std::string(5, 'T');
        expected_cigar = "10=1I5=";
        expected_score = 92;
    }

    SECTION("penalized deletion in core")
    {
        query_read = std::string(5, 'T') + std::string(5, 'A') + std::string(5, 'T');
        germline_ref = std::string(5, 'T') + std::string(5, 'A') + "C" + std::string(5, 'T');
        expected_cigar = "10=1D5=";
        expected_score = 92;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_cigar, expected_score });
}

TEST_CASE("Aligner J gene multiplex dummy alignments.", "[aligner][J_gene][dummy]")
{

    /*
    * A nucleotides for NT of germline of interest
    * T nucleotides for other NT in the read, or to make a germline variation allowing to constrain alignment.
    * G nucleotides for mismatches
    * C nucleotides for penalized gaps
    * 
    * Score costs are set with prime numbers such that total score mismatch is easier to debug. 
    */
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const int gap_penalty = 13;
    std::string query_read;
    std::string germline_ref;
    std::string expected_cigar;
    double expected_score;

    SECTION("complete alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_cigar = "20=";
        expected_score = 140;
    }

    SECTION("short read, free trailing reference deletion")
    {
        query_read = std::string(1, 'T') + std::string(5, 'A');
        germline_ref = std::string(1, 'T') + std::string(9, 'A');
        expected_cigar = "6=4D";
        expected_score = 42;
    }

    SECTION("partial overlap read, no 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(1, 'T') + std::string(19, 'A');
        expected_cigar = "13I7=13D";
        expected_score = 49;
    }

    SECTION("partial overlap read, single free 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(1, 'A') + std::string(1, 'T') + std::string(19, 'A');
        expected_cigar = "1D13I7=13D";
        expected_score = 49;
    }

    SECTION("partial overlap read, three free 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(3, 'A') + std::string(1, 'T') + std::string(19, 'A');
        expected_cigar = "3D13I7=13D";
        expected_score = 49;
    }

    SECTION("free leading germline treated as deletion")
    {
        query_read = std::string(14, 'T');
        germline_ref = std::string(2, 'A') + std::string(14, 'T') + std::string(6, 'A');
        expected_cigar = "2D14=6D";
        expected_score = 98;
    }

    SECTION("mismatch in long match")
    {
        query_read = std::string(20, 'A');
        query_read[4] = 'G';
        germline_ref = std::string(20, 'A');
        expected_cigar = "4=1X15=";
        expected_score = 122;
    }

    SECTION("mismatch in last NT")
    {
        query_read = std::string(20, 'A');
        query_read[19] = 'G';
        germline_ref = std::string(20, 'A');
        expected_cigar = "19=1X";
        expected_score = 122;
    }

    SECTION("mismatch in first NT seen as free gap")
    {
        query_read = std::string(19, 'A') + std::string(1, 'T');
        query_read[0] = 'G';
        germline_ref = std::string(19, 'A') + std::string(1, 'T');
        expected_cigar = "1D1I19=";
        expected_score = 133;
    }

    SECTION("penalized insertion")
    {
        query_read = std::string(5, 'A') + "C" + std::string(15, 'A');
        germline_ref = std::string(20, 'A');
        expected_cigar = "5=1I15=";
        expected_score = 127;
    }

    SECTION("penalized deletion")
    {
        query_read = "T" + std::string(17, 'A') + "T";
        germline_ref = "T" + std::string(7, 'A') + "T" + std::string(10, 'A') + "T";
        expected_cigar = "8=1D11=";
        expected_score = 120;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, J_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_cigar, expected_score });
}

TEST_CASE("Aligner V gene best alignment matches expected CIGAR and score on realistic data.",
          "[aligner][V_gene][realistic]")
{
    // Align a single query to a single reference, check that best alignment matches
    const Matrix<double> matrix = build_test_score_matrix(5.0, -14.0);
    const int gap_penalty = 30;
    std::string query_read;
    std::string germline_ref;
    std::string expected_cigar;
    double expected_score;

    SECTION("Long match without V deletions")
    {
        query_read = "ACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGATATAGGACTAGATTCACAGATACGCA";
        germline_ref = "GATACTGGAATTACCCAGACACCAAAATACCTGGTCACAGCAATGGGGAGTAAAAGGACAATGAAACGTGAGCATCTGGGACATGATTCTATGTA"
                       "TTGGTACAGACAGAAAGCTAAGAAATCCCTGGAGTTCATGTTTTACTACAACTGTAAGGAATTCATTGAAAACAAGACTGTGCCAAATCACTTCA"
                       "CACCTGAATGCCCTGACAGCTCTCGCTTATACCTTCATGTGGTCGCACTGCAGCAAGAAGACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGA";
        expected_cigar = "250D34=26I";
        expected_score = 170;
    }

    SECTION("Long match with single V deletion")
    {
        query_read = "ACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGTGTGTCCCGGACAGACGACTATGGCTA";
        germline_ref =
                "GACACAGCTGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAACGACAAGTCCATTAAATGTGAACAAAATCTGGGCCATGATACTATGTATTGG"
                "TATAAACAGGACTCTAAGAAATTTCTGAAGATAATGTTTAGCTACAATAACAAGGAGATCATTATAAATGAAACAGTTCCAAATCGATTCTCACCTAAATCT"
                "CCAGACAAAGCTAAATTAAATCTTCACATCAATTCCCTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA";
        expected_cigar = "253D33=1D27I";
        expected_score = 165;
    }

    SECTION("Long match with long V deletions")
    {
        query_read = "CCAACCAGACAGCTCTTTACTTCTGTGCCACCCTACGAACAGGGAAAGGAACACTGAAGC";
        germline_ref =
                "GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGG"
                "TATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGA"
                "CAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_cigar = "247D32=9D28I";
        expected_score = 160;
    }

    SECTION("Random NT mismatch")
    {
        query_read = "CAGACAGCTCTTTACTTCTGTGTCACCAGTGATTTGCACTGGACAGGGGGAAGAGACCCA";
        germline_ref =
                "GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGG"
                "TATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGA"
                "CAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_cigar = "252D22=1X13=24I";
        expected_score = 161;
    }

    SECTION("First NT mismatch")
    {
        query_read = "GCTGTGTACTTCTGTGCCAGCAGTTCGGGACTAGCGGGGAATGCCAGCCGCAGATACGCA";
        germline_ref =
                "AATGCTGGTGTCACTCAGACCCCAAAATTCCGCATCCTGAAGATAGGACAGAGCATGACACTGCAGTGTGCCCAGGATATGAACCATAACTACATGTACTGG"
                "TATCGACAAGACCCAGGCATGGGGCTGAAGCTGATTTATTATTCAGTTGGTGCTGGTATCACTGACAAAGGAGAAGTCCCGAATGGCTACAACGTCTCCAGA"
                "TCAACCACAGAGGATTTCCCGCTCAGGCTGGAGTTGGCTGCTCCCTCCCAGACATCTGTGTACTTCTGTGCCAGCAGTTACTC";
        expected_cigar = "258D1X24=4D35I";
        expected_score = 106;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_cigar, expected_score });
}

TEST_CASE("Legacy Aligner strict set matching when best_only is false", "[aligner][sw][cigar]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::string query = "ACGT";
    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", "ACGT" }, { "g2", "ACGT" } };

    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query, -1000.0, false, false, INT16_MIN, INT16_MAX);

    assert_alignment_set_matches(alignments, query, genomic_templates, { { "g1", "4=", 8.0 }, { "g2", "4=", 8.0 } });
}


