#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>



#include <igor/Core/Aligner.h>
#ifdef IGOR_WITH_SEQAN2
#include <igor/Core/SeqAn2Aligner.h>
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

struct ExpectedAlignment {
    std::string gene_name;
    std::string cigar;
    double score;
};

struct ActualAlignment {
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

Aligner make_legacy_aligner(const Matrix<double>& matrix, int gap_penalty, Gene_class gene_class,
                            const std::vector<std::pair<std::string, std::string>>& genomic_templates)
{
    Aligner aligner(matrix, gap_penalty, gene_class);
    aligner.set_genomic_sequences(genomic_templates);
    return aligner;
}

#ifdef IGOR_WITH_SEQAN2
std::unique_ptr<SeqAn2Aligner> make_seqan2_aligner(const Matrix<double>& matrix, int gap_penalty, Gene_class gene_class,
                                                   const std::vector<std::pair<std::string, std::string>>& genomic_templates)
{
    auto aligner = std::make_unique<SeqAn2Aligner>(matrix, gap_penalty, gene_class, AlignmentPreset::with_band(20));
    aligner->set_genomic_sequences(genomic_templates);
    return aligner;
}
#endif

std::string find_germline_sequence(const std::vector<std::pair<std::string, std::string>>& genomic_templates,
                                   const std::string& gene_name)
{
    for (const auto& entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry.second;
        }
    }
    return std::string();
}

std::pair<std::string, std::string> find_genomic_template(
    const std::vector<std::pair<std::string, std::string>>& genomic_templates, const std::string& gene_name)
{
    for (const auto& entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry;
        }
    }
    return std::pair<std::string, std::string>{};
}

std::string cigar_visual(const std::string& cigar, const std::string& read, const std::string& germline)
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

    for (const auto& entry : parse_cigar(cigar)) {
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
        out << "\n      read " << read_gapped.substr(i, width)
            << "\n           " << match_line.substr(i, width)
            << "\n      germ " << germ_gapped.substr(i, width);
    }
    return out.str();
}

std::vector<ActualAlignment> normalize_actuals(const std::forward_list<Alignment_data>& alignments)
{
    std::vector<ActualAlignment> normalized;
    for (const auto& aln : alignments) {
        normalized.push_back({aln.gene_name, alignment_data_to_cigar(aln), aln.score});
    }
    std::sort(normalized.begin(), normalized.end(), [](const ActualAlignment& a, const ActualAlignment& b) {
        return std::tie(a.gene_name, a.cigar, a.score) < std::tie(b.gene_name, b.cigar, b.score);
    });
    return normalized;
}

std::vector<ExpectedAlignment> normalize_expected(std::vector<ExpectedAlignment> expected)
{
    std::sort(expected.begin(), expected.end(), [](const ExpectedAlignment& a, const ExpectedAlignment& b) {
        return std::tie(a.gene_name, a.cigar, a.score) < std::tie(b.gene_name, b.cigar, b.score);
    });
    return expected;
}

ActualAlignment best_alignment(const std::forward_list<Alignment_data>& alignments)
{
    auto best_it = std::max_element(
        alignments.begin(), alignments.end(),
        [](const Alignment_data& a, const Alignment_data& b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());
    return {best_it->gene_name, alignment_data_to_cigar(*best_it), best_it->score};
}

void add_cigar_visual_info(const std::string& query,
                           const std::vector<std::pair<std::string, std::string>>& genomic_templates,
                           const std::string& label, const std::string& gene_name, const std::string& cigar)
{
    const std::string germline = find_germline_sequence(genomic_templates, gene_name);
    if (germline.empty()) {
        INFO(label << " gene=" << gene_name << " has no matching germline sequence for visual explanation");
        return;
    }
    INFO(label << " gene=" << gene_name << " cigar=" << cigar << cigar_visual(cigar, query, germline));
}

void assert_alignment_matches(const Alignment_data& alignment, const std::string& query,
                              const std::pair<std::string, std::string>& genomic_template,
                              const ExpectedAlignment& expected)
{
    const std::string actual_cigar = alignment_data_to_cigar_full_span(alignment, query.length(), genomic_template.second.length());

    INFO("expected: gene=" << expected.gene_name << " cigar=" << expected.cigar << " score=" << expected.score);
    INFO("actual: gene=" << alignment.gene_name << " cigar=" << actual_cigar << " score=" << alignment.score);

    if (actual_cigar != expected.cigar) {
        const std::vector<std::pair<std::string, std::string>> templates = {genomic_template};
        add_cigar_visual_info(query, templates, "expected", expected.gene_name, expected.cigar);
        add_cigar_visual_info(query, templates, "actual", alignment.gene_name, actual_cigar);
    }

    REQUIRE(alignment.gene_name == expected.gene_name);
    REQUIRE(actual_cigar == expected.cigar);
    REQUIRE_THAT(alignment.score, WithinRel(expected.score));
}

void assert_best_alignment_matches(const std::forward_list<Alignment_data>& alignments, const std::string& query,
                                   const std::vector<std::pair<std::string, std::string>>& genomic_templates,
                                   const ExpectedAlignment& expected)
{
    REQUIRE(!alignments.empty());
    auto best_it = std::max_element(
        alignments.begin(), alignments.end(),
        [](const Alignment_data& a, const Alignment_data& b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());

    auto genomic_template = find_genomic_template(genomic_templates, expected.gene_name);
    if (genomic_template.first.empty()) {
        genomic_template = find_genomic_template(genomic_templates, best_it->gene_name);
    }
    assert_alignment_matches(*best_it, query, genomic_template, expected);
}

void assert_alignment_set_matches(const std::forward_list<Alignment_data>& alignments, const std::string& query,
                                  const std::vector<std::pair<std::string, std::string>>& genomic_templates,
                                  std::vector<ExpectedAlignment> expected)
{
    std::vector<const Alignment_data*> actual;
    for (const auto& aln : alignments) {
        actual.push_back(&aln);
    }
    std::sort(actual.begin(), actual.end(), [](const Alignment_data* a, const Alignment_data* b) {
         return std::make_tuple(a->gene_name, alignment_data_to_cigar(*a), a->score) <
             std::make_tuple(b->gene_name, alignment_data_to_cigar(*b), b->score);
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

TEST_CASE("Aligner V gene best alignment matches expected CIGAR and score on realistic data.", "[aligner][V_gene][realistic]")
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
        query_read="ACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGATATAGGACTAGATTCACAGATACGCA";
        germline_ref="GATACTGGAATTACCCAGACACCAAAATACCTGGTCACAGCAATGGGGAGTAAAAGGACAATGAAACGTGAGCATCTGGGACATGATTCTATGTATTGGTACAGACAGAAAGCTAAGAAATCCCTGGAGTTCATGTTTTACTACAACTGTAAGGAATTCATTGAAAACAAGACTGTGCCAAATCACTTCACACCTGAATGCCCTGACAGCTCTCGCTTATACCTTCATGTGGTCGCACTGCAGCAAGAAGACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGA";
        expected_cigar="250D34=26I";
        expected_score=170;
    }

    SECTION("Long match with single V deletion")
    {
        query_read="ACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGTGTGTCCCGGACAGACGACTATGGCTA";
        germline_ref="GACACAGCTGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAACGACAAGTCCATTAAATGTGAACAAAATCTGGGCCATGATACTATGTATTGGTATAAACAGGACTCTAAGAAATTTCTGAAGATAATGTTTAGCTACAATAACAAGGAGATCATTATAAATGAAACAGTTCCAAATCGATTCTCACCTAAATCTCCAGACAAAGCTAAATTAAATCTTCACATCAATTCCCTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA";
        expected_cigar="253D33=1D27I";
        expected_score=165;
    }

    SECTION("Long match with long V deletions")
    {
        query_read="CCAACCAGACAGCTCTTTACTTCTGTGCCACCCTACGAACAGGGAAAGGAACACTGAAGC";
        germline_ref="GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGGTATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGACAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_cigar="247D32=9D28I";
        expected_score=160;
    }

    SECTION("Random NT mismatch")
    {
        query_read="CAGACAGCTCTTTACTTCTGTGTCACCAGTGATTTGCACTGGACAGGGGGAAGAGACCCA";
        germline_ref="GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGGTATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGACAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_cigar="252D22=1X13=24I";
        expected_score=161;
    }

    SECTION("First NT mismatch")
    {
        query_read="GCTGTGTACTTCTGTGCCAGCAGTTCGGGACTAGCGGGGAATGCCAGCCGCAGATACGCA";
        germline_ref="AATGCTGGTGTCACTCAGACCCCAAAATTCCGCATCCTGAAGATAGGACAGAGCATGACACTGCAGTGTGCCCAGGATATGAACCATAACTACATGTACTGGTATCGACAAGACCCAGGCATGGGGCTGAAGCTGATTTATTATTCAGTTGGTGCTGGTATCACTGACAAAGGAGAAGTCCCGAATGGCTACAACGTCTCCAGATCAACCACAGAGGATTTCCCGCTCAGGCTGGAGTTGGCTGCTCCCTCCCAGACATCTGTGTACTTCTGTGCCAGCAGTTACTC";
        expected_cigar="258D1X24=4D35I";
        expected_score=106;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", germline_ref}};
    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, {"g1", expected_cigar, expected_score});

}

TEST_CASE("Legacy Aligner best alignment supports mismatch and indel CIGAR assertions", "[aligner][sw][cigar]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "LALALA"}};

    SECTION("single mismatch")
    {
        auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
        const auto alignments = aligner.align_seq("ACGTTCGT", -1000.0, true, true, INT16_MIN, INT16_MAX);
        assert_best_alignment_matches(alignments, "ACGTTCGT", genomic_templates, {"g1", "4=1X3=", 13.0});
    }

    SECTION("single insertion in read")
    {
        const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "ACGTACGT"}};
        auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
        const auto alignments = aligner.align_seq("ACGTTACGT", -1000.0, true, true, INT16_MIN, INT16_MAX);
        assert_best_alignment_matches(alignments, "ACGTTACGT", genomic_templates, {"g1", "4=1I4=", 14.0});
    }

    SECTION("single deletion in read")
    {
        const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "ACGTACGT"}};
        auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
        const auto alignments = aligner.align_seq("ACGACGT", -1000.0, true, true, INT16_MIN, INT16_MAX);
        assert_best_alignment_matches(alignments, "ACGACGT", genomic_templates, {"g1", "3=1D4=", 12.0});
    }
}

TEST_CASE("Legacy Aligner strict set matching when best_only is false", "[aligner][sw][cigar]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::string query = "ACGT";
    const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "ACGT"}, {"g2", "ACGT"}};

    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query, -1000.0, false, false, INT16_MIN, INT16_MAX);

    assert_alignment_set_matches(alignments, query, genomic_templates,
                                 {{"g1", "4=", 8.0}, {"g2", "4=", 8.0}});
}

TEST_CASE("Legacy Aligner align_seqs supports indexed CIGAR assertions", "[aligner][sw][cigar]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "ACGTACGT"}};

    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);

    const std::vector<std::pair<const int, const std::string>> reads = {
        {10, "ACGTACGT"},
        {11, "ACGT"},
    };

    const auto indexed = aligner.align_seqs(reads, -1000.0, true, true);

    REQUIRE(indexed.count(10) == 1);
    REQUIRE(indexed.count(11) == 1);
    assert_best_alignment_matches(indexed.at(10), "ACGTACGT", genomic_templates, {"g1", "8=", 16.0});
    assert_best_alignment_matches(indexed.at(11), "ACGT", genomic_templates, {"g1", "4=", 8.0});
}

#ifdef IGOR_WITH_SEQAN2
TEST_CASE("Legacy and SeqAn2 best alignments share CIGAR and score on deterministic cases", "[aligner][sw][cigar][seqan2]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::vector<std::pair<std::string, std::string>> genomic_templates = {{"g1", "ACGTACGT"}};

    auto legacy = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    auto seqan2 = make_seqan2_aligner(matrix, gap_penalty, D_gene, genomic_templates);

    SECTION("exact match")
    {
        const std::string query = "ACGTACGT";
        const auto legacy_alns = legacy.align_seq(query, -1000.0, true, true, INT16_MIN, INT16_MAX);
        const auto seqan2_alns = seqan2->align_seq(query, -1000.0, true, true);

        const ActualAlignment legacy_best = best_alignment(legacy_alns);
        const ActualAlignment seqan2_best = best_alignment(seqan2_alns);

        INFO("legacy: gene=" << legacy_best.gene_name << " cigar=" << legacy_best.cigar << " score=" << legacy_best.score);
        INFO("seqan2: gene=" << seqan2_best.gene_name << " cigar=" << seqan2_best.cigar << " score=" << seqan2_best.score);

        if (legacy_best.cigar != seqan2_best.cigar) {
            add_cigar_visual_info(query, genomic_templates, "legacy", legacy_best.gene_name, legacy_best.cigar);
            add_cigar_visual_info(query, genomic_templates, "seqan2", seqan2_best.gene_name, seqan2_best.cigar);
        }

        REQUIRE(legacy_best.gene_name == seqan2_best.gene_name);
        REQUIRE_THAT(legacy_best.score, WithinRel(seqan2_best.score));
        REQUIRE(legacy_best.cigar == seqan2_best.cigar);
    }

    SECTION("single insertion")
    {
        const std::string query = "ACGTTACGT";
        const auto legacy_alns = legacy.align_seq(query, -1000.0, true, true, INT16_MIN, INT16_MAX);
        const auto seqan2_alns = seqan2->align_seq(query, -1000.0, true, true);

        const ActualAlignment legacy_best = best_alignment(legacy_alns);
        const ActualAlignment seqan2_best = best_alignment(seqan2_alns);

        INFO("legacy: gene=" << legacy_best.gene_name << " cigar=" << legacy_best.cigar << " score=" << legacy_best.score);
        INFO("seqan2: gene=" << seqan2_best.gene_name << " cigar=" << seqan2_best.cigar << " score=" << seqan2_best.score);

        if (legacy_best.cigar != seqan2_best.cigar) {
            add_cigar_visual_info(query, genomic_templates, "legacy", legacy_best.gene_name, legacy_best.cigar);
            add_cigar_visual_info(query, genomic_templates, "seqan2", seqan2_best.gene_name, seqan2_best.cigar);
        }

        REQUIRE(legacy_best.gene_name == seqan2_best.gene_name);
        REQUIRE_THAT(legacy_best.score, WithinRel(seqan2_best.score));
        REQUIRE(legacy_best.cigar == seqan2_best.cigar);
    }
}
#endif
