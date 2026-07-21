#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <igor/Core/Aligner.h>
#include <igor/Core/AlignerInternal.h>
#include "AlignerTestUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Catch::Matchers::WithinRel;
using namespace igor::test::align;

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
    std::string expected_core_cigar;
    std::string expected_extended_cigar;
    double expected_score;

    SECTION("complete alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "20=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 140;
    }

    SECTION("short read, free leading reference deletion")
    {
        query_read = std::string(5, 'A') + std::string(1, 'T');
        germline_ref = std::string(9, 'A') + std::string(1, 'T');
        expected_core_cigar = "4N6=";
        expected_extended_cigar = expected_core_cigar;  // Extended: N for reference-only gaps
        expected_score = 42;
    }

    SECTION("partial overlap read, no 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T');
        expected_core_cigar = "13N7=13S";
        expected_extended_cigar = expected_core_cigar;  // Extended: N for ref gaps, S for query gaps
        expected_score = 49;
    }

    SECTION("partial overlap read, single free 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T') + std::string(1, 'A');
        expected_core_cigar = "13N7=1N13S";
        expected_extended_cigar = "13N7=1X12S";  // Extended: N for ref gaps, S for query gaps
        expected_score = 49;
    }

    SECTION("partial overlap read, three free 3' deletion")
    {
        query_read = std::string(6, 'A') + std::string(14, 'T');
        germline_ref = std::string(19, 'A') + std::string(1, 'T') + std::string(3, 'A');
        expected_core_cigar = "13N7=3N13S";
        expected_extended_cigar = "13N7=3X10S";  // Extended: N for ref gaps, S for query gaps
        expected_score = 49;
    }

    SECTION("free trailing germline treated as deletion")
    {
        query_read = std::string(14, 'T');
        germline_ref = std::string(6, 'A') + std::string(14, 'T') + std::string(2, 'A');
        expected_core_cigar = "6N14=2N";
        expected_extended_cigar = expected_core_cigar;  // Extended: N for reference gaps
        expected_score = 98;
    }

    SECTION("mismatch in long match")
    {
        query_read = std::string(20, 'A');
        query_read[4] = 'G';
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "4=1X15=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 122;
    }

    SECTION("mismatch in first NT")
    {
        query_read = std::string(20, 'A');
        query_read[0] = 'G';
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "1X19=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 122;
    }

    SECTION("mismatch in last NT seen as free gap")
    {
        query_read = std::string(20, 'A');
        query_read[19] = 'G';
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "19=1N1S";
        expected_extended_cigar = "19=1X";  // Extended: trailing query gap S
        expected_score = 133;
    }

    SECTION("penalized insertion")
    {
        query_read = std::string(5, 'A') + "C" + std::string(15, 'A');
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "5=1I15=";  // Core insertion - keep I
        expected_extended_cigar = expected_core_cigar;
        expected_score = 127;
    }

    SECTION("penalized deletion")
    {
        query_read = "T" + std::string(17, 'A') + "T";
        germline_ref = "T" + std::string(7, 'A') + "T" + std::string(10, 'A') + "T";
        expected_core_cigar = "8=1D11=";  // Core deletion - keep D
        expected_extended_cigar = expected_core_cigar;
        expected_score = 120;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_core_cigar, expected_extended_cigar, expected_score });
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
    std::string expected_core_cigar;
    std::string expected_extended_cigar;
    double expected_score;

    SECTION("complete local alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "20=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 140;
    }

    SECTION("short local core in long read")
    {
        query_read = std::string(4, 'T') + std::string(6, 'A') + std::string(4, 'T');
        germline_ref = std::string(6, 'A');
        expected_core_cigar = "4S6=4S";  // Updated to AIRR: S for query gaps at ends
        expected_extended_cigar = expected_core_cigar;
        expected_score = 42;
    }

    SECTION("long local core in long read")
    {
        query_read = std::string(2, 'T') + std::string(16, 'A') + std::string(2, 'T');
        germline_ref = std::string(16, 'A');
        expected_core_cigar = "2S16=2S";  // Updated to AIRR: S for query gaps at ends
        expected_extended_cigar = expected_core_cigar;
        expected_score = 112;
    }

    SECTION("mismatch in several places")
    {
        query_read = "AAAGAAAAGAAAAGAA";
        germline_ref = std::string(16, 'A');
        expected_core_cigar = "3=1X4=1X4=1X2=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 58;
    }

    SECTION("left mismatch shortens alignment")
    {
        query_read = "AGAGAAAAGAAAAGAA";
        germline_ref = std::string(16, 'A');
        expected_core_cigar = "4N4S4=1X4=1X2=";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "1=1X1=1X4=1X4=1X2=";
        expected_score = 48;
    }

    SECTION("right mismatch shortens alignment")
    {
        query_read = "AAAGAAAAGAAAGGAA";
        germline_ref = std::string(16, 'A');
        expected_core_cigar = "3=1X4=1X3=4N4S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "3=1X4=1X3=2X2=";
        expected_score = 48;
    }

    SECTION("local overlap with one trailing germline base")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(1, 'T');
        expected_core_cigar = "8=1N";  // Updated to AIRR: N for trailing reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("local overlap with three trailing germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(3, 'T');
        expected_core_cigar = "8=3N";  // Updated to AIRR: N for trailing reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("local overlap with five trailing germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(8, 'A') + std::string(5, 'T');
        expected_core_cigar = "8=5N";  // Updated to AIRR: N for trailing reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("local overlap with one leading germline base")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(1, 'T') + std::string(8, 'A');
        expected_core_cigar = "1N8=";  // Updated to AIRR: N for leading reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("local overlap with three leading germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(3, 'T') + std::string(8, 'A');
        expected_core_cigar = "3N8=";  // Updated to AIRR: N for leading reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("local overlap with five leading germline bases")
    {
        query_read = std::string(8, 'A');
        germline_ref = std::string(5, 'T') + std::string(8, 'A');
        expected_core_cigar = "5N8=";  // Updated to AIRR: N for leading reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 56;
    }

    SECTION("penalized insertion in core")
    {
        query_read = std::string(5, 'T') + std::string(5, 'A') + "C" + std::string(5, 'T');
        germline_ref = std::string(5, 'T') + std::string(5, 'A') + std::string(5, 'T');
        expected_core_cigar = "10=1I5=";  // Core insertion - keep I
        expected_extended_cigar = expected_core_cigar;
        expected_score = 92;
    }

    SECTION("penalized deletion in core")
    {
        query_read = std::string(5, 'T') + std::string(5, 'A') + std::string(5, 'T');
        germline_ref = std::string(5, 'T') + std::string(5, 'A') + "C" + std::string(5, 'T');
        expected_core_cigar = "10=1D5=";  // Core deletion - keep D
        expected_extended_cigar = expected_core_cigar;
        expected_score = 92;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_core_cigar, expected_extended_cigar, expected_score });
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
    std::string expected_core_cigar;
    std::string expected_extended_cigar;
    double expected_score;

    SECTION("complete alignment")
    {
        query_read = std::string(20, 'A');
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "20=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 140;
    }

    SECTION("short read, free trailing reference deletion")
    {
        query_read = std::string(1, 'T') + std::string(5, 'A');
        germline_ref = std::string(1, 'T') + std::string(9, 'A');
        expected_core_cigar = "6=4N";  // Updated to AIRR: N for trailing reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 42;
    }

    SECTION("partial overlap read, no 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(1, 'T') + std::string(19, 'A');
        expected_core_cigar = "13S7=13N";  // Updated to AIRR: S for query gaps, N for ref gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 49;
    }

    SECTION("partial overlap read, single free 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(1, 'A') + std::string(1, 'T') + std::string(19, 'A');
        expected_core_cigar = "1N13S7=13N";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "12S1X7=13N";
        expected_score = 49;
    }

    SECTION("partial overlap read, three free 5' deletion")
    {
        query_read = std::string(14, 'T') + std::string(6, 'A');
        germline_ref = std::string(3, 'A') + std::string(1, 'T') + std::string(19, 'A');
        expected_core_cigar = "3N13S7=13N";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "10S3X7=13N";
        expected_score = 49;
    }

    SECTION("free leading germline treated as deletion")
    {
        query_read = std::string(14, 'T');
        germline_ref = std::string(2, 'A') + std::string(14, 'T') + std::string(6, 'A');
        expected_core_cigar = "2N14=6N";  // Updated to AIRR: N for reference gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 98;
    }

    SECTION("mismatch in long match")
    {
        query_read = std::string(20, 'A');
        query_read[4] = 'G';
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "4=1X15=";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 122;
    }

    SECTION("mismatch in last NT")
    {
        query_read = std::string(20, 'A');
        query_read[19] = 'G';
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "19=1X";
        expected_extended_cigar = expected_core_cigar;
        expected_score = 122;
    }

    SECTION("mismatch in first NT seen as free gap")
    {
        query_read = std::string(19, 'A') + std::string(1, 'T');
        query_read[0] = 'G';
        germline_ref = std::string(19, 'A') + std::string(1, 'T');
        expected_core_cigar = "1N1S19=";  // Updated to AIRR: N for ref gap, S for query gap
        expected_extended_cigar = "1X19=";
        expected_score = 133;
    }

    SECTION("penalized insertion")
    {
        query_read = std::string(5, 'A') + "C" + std::string(15, 'A');
        germline_ref = std::string(20, 'A');
        expected_core_cigar = "5=1I15=";  // Core insertion - keep I
        expected_extended_cigar = expected_core_cigar;
        expected_score = 127;
    }

    SECTION("penalized deletion")
    {
        query_read = "T" + std::string(17, 'A') + "T";
        germline_ref = "T" + std::string(7, 'A') + "T" + std::string(10, 'A') + "T";
        expected_core_cigar = "8=1D11=";  // Core deletion - keep D
        expected_extended_cigar = expected_core_cigar;
        expected_score = 120;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, J_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_core_cigar, expected_extended_cigar, expected_score });
}

TEST_CASE("Aligner V gene best alignment matches expected CIGAR and score on realistic data.",
          "[aligner][V_gene][realistic]")
{
    // Align a single query to a single reference, check that best alignment matches
    const Matrix<double> matrix = build_test_score_matrix(5.0, -14.0);
    const int gap_penalty = 30;
    std::string query_read;
    std::string germline_ref;
    std::string expected_core_cigar;
    std::string expected_extended_cigar;
    double expected_score;

    SECTION("Long match without V deletions")
    {
        query_read = "ACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGATATAGGACTAGATTCACAGATACGCA";
        germline_ref = "GATACTGGAATTACCCAGACACCAAAATACCTGGTCACAGCAATGGGGAGTAAAAGGACAATGAAACGTGAGCATCTGGGACATGATTCTATGTA"
                       "TTGGTACAGACAGAAAGCTAAGAAATCCCTGGAGTTCATGTTTTACTACAACTGTAAGGAATTCATTGAAAACAAGACTGTGCCAAATCACTTCA"
                       "CACCTGAATGCCCTGACAGCTCTCGCTTATACCTTCATGTGGTCGCACTGCAGCAAGAAGACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGA";
        expected_core_cigar = "250N34=26S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 170;
    }

    SECTION("Long match with single V deletion")
    {
        query_read = "ACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGTGTGTCCCGGACAGACGACTATGGCTA";
        germline_ref =
                "GACACAGCTGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAACGACAAGTCCATTAAATGTGAACAAAATCTGGGCCATGATACTATGTATTGG"
                "TATAAACAGGACTCTAAGAAATTTCTGAAGATAATGTTTAGCTACAATAACAAGGAGATCATTATAAATGAAACAGTTCCAAATCGATTCTCACCTAAATCT"
                "CCAGACAAAGCTAAATTAAATCTTCACATCAATTCCCTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA";
        expected_core_cigar = "253N33=1N27S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "253N33=1X26S";
        expected_score = 165;
    }

    SECTION("Long match with long V deletions")
    {
        query_read = "CCAACCAGACAGCTCTTTACTTCTGTGCCACCCTACGAACAGGGAAAGGAACACTGAAGC";
        germline_ref =
                "GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGG"
                "TATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGA"
                "CAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_core_cigar = "247N32=9N28S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "247N32=9X19S";
        expected_score = 160;
    }

    SECTION("Random NT mismatch")
    {
        query_read = "CAGACAGCTCTTTACTTCTGTGTCACCAGTGATTTGCACTGGACAGGGGGAAGAGACCCA";
        germline_ref =
                "GATGCTGATGTTACCCAGACCCCAAGGAATAGGATCACAAAGACAGGAAAGAGGATTATGCTGGAATGTTCTCAGACTAAGGGTCATGATAGAATGTACTGG"
                "TATCGACAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGA"
                "CAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACCAGTGATTTG";
        expected_core_cigar = "252N22=1X13=24S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = expected_core_cigar;
        expected_score = 161;
    }

    SECTION("First NT mismatch")
    {
        query_read = "GCTGTGTACTTCTGTGCCAGCAGTTCGGGACTAGCGGGGAATGCCAGCCGCAGATACGCA";
        germline_ref =
                "AATGCTGGTGTCACTCAGACCCCAAAATTCCGCATCCTGAAGATAGGACAGAGCATGACACTGCAGTGTGCCCAGGATATGAACCATAACTACATGTACTGG"
                "TATCGACAAGACCCAGGCATGGGGCTGAAGCTGATTTATTATTCAGTTGGTGCTGGTATCACTGACAAAGGAGAAGTCCCGAATGGCTACAACGTCTCCAGA"
                "TCAACCACAGAGGATTTCCCGCTCAGGCTGGAGTTGGCTGCTCCCTCCCAGACATCTGTGTACTTCTGTGCCAGCAGTTACTC";
        expected_core_cigar = "258N1X24=4N35S";  // Updated to AIRR: N for ref gaps, S for query gaps
        expected_extended_cigar = "258N1X24=4X31S";
        expected_score = 106;
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query_read, -1000.0, true, true, INT16_MIN, INT16_MAX);
    assert_best_alignment_matches(alignments, query_read, genomic_templates, { "g1", expected_core_cigar, expected_extended_cigar, expected_score });
}

TEST_CASE("Legacy Aligner strict set matching when best_only is false", "[aligner][sw][cigar]")
{
    const Matrix<double> matrix = build_test_score_matrix();
    const int gap_penalty = 2;
    const std::string query = "ACGT";
    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", "ACGT" }, { "g2", "ACGT" } };

    auto aligner = make_legacy_aligner(matrix, gap_penalty, D_gene, genomic_templates);
    const auto alignments = aligner.align_seq(query, -1000.0, false, false, INT16_MIN, INT16_MAX);

    assert_alignment_set_matches(alignments, query, genomic_templates, { { "g1", "4=", "4=", 8.0 }, { "g2", "4=", "4=", 8.0 } });
}

// ============================================================================
//  Smith-Waterman DP-level tests
// ============================================================================

TEST_CASE("fill_sw_score_matrix produces the expected DP matrices", "[aligner][sw][dp_matrix]")
{
    // Simple fixed example with a hand-derived reference.
    const std::string query = "AAAAAT"; // data sequence -> DP rows
    const std::string reference = "AAAAAAAAAT"; // genomic sequence -> DP cols

    // TODO(user): must match the values used when deriving the reference matrices below.
    const double match_score = 7;
    const double mismatch_score = -11;
    const int gap_penalty = 13;

    SwDPConfig config;
    config.score_threshold = -1000.0; // unused by fill_sw_score_matrix itself
    config.best_only = false; // unused by fill_sw_score_matrix itself
    config.min_offset = INT16_MIN; // unused by fill_sw_score_matrix itself
    config.max_offset = INT16_MAX; // unused by fill_sw_score_matrix itself
    config.substitution_matrix = build_test_score_matrix(match_score, mismatch_score);
    config.gap_penalty = gap_penalty;
    config.alignment_mode = SwAlignmentMode{ false, false, true, false, false }; // full semi global (V-gene-like)

    const Int_Str int_query = nt2int(query);
    const Int_Str int_reference = nt2int(reference);

    // --- Setup mirroring sw_align()'s pipeline up to (not including) fill_sw_score_matrix ---
    const swalign::SwPreparedInputs prepared = swalign::prepare_sw_inputs(int_query, int_reference, config);
    const int n_rows = static_cast<int>(prepared.data_sequence.size()) + 1; // 7
    const int n_cols = static_cast<int>(prepared.genomic_sequence.size()) + 1; // 11

    swalign::SwDPState dp(n_rows, n_cols);
    swalign::initialize_sw_matrices(dp, config);

    // --- Call under test ---
    swalign::fill_sw_score_matrix(prepared.data_sequence, prepared.genomic_sequence, dp, config);

    // --- Expected matrices: row i = query position (0 = init boundary), col j = ref position (0 = init boundary) ---
    // row_memory(i,j): 1 = predecessor is (i-1,*) [diagonal or up move], 0 = predecessor is (i,j-1) [left move]
    // col_memory(i,j): 1 = predecessor is (*,j-1) [diagonal or left move], 0 = predecessor is (i-1,j) [up move]
    std::vector<double> scores_arr = {
        0,0,0,0,0,0,0,0,0,0,0,
        -13,7,7,7,7,7,7,7,7,7,-6,
        -26,-6,14,14,14,14,14,14,14,14,1,
        -39,-19,1,21,21,21,21,21,21,21,8,
        -52,-32,-12,8,28,28,28,28,28,28,15,
        -65,-45,-25,-5,15,35,35,35,35,35,22,
        -78,-58,-38,-18,2,22,24,24,24,24,42,
    };
    std::vector<int> row_mem_arr = {
        0,0,0,0,0,0,0,0,0,0,0,
        0,1,1,1,1,1,1,1,1,1,0,
        0,1,1,1,1,1,1,1,1,1,0,
        0,1,1,1,1,1,1,1,1,1,0,
        0,1,1,1,1,1,1,1,1,1,0,
        0,1,1,1,1,1,1,1,1,1,0,
        0,1,1,1,1,1,1,1,1,1,1,
    };
    std::vector<int> col_mem_arr = {
        0,0,0,0,0,0,0,0,0,0,0,
        0,1,1,1,1,1,1,1,1,1,1,
        0,1,1,1,1,1,1,1,1,1,1,
        0,1,1,1,1,1,1,1,1,1,1,
        0,1,1,1,1,1,1,1,1,1,1,
        0,1,1,1,1,1,1,1,1,1,1,
        0,0,0,0,0,0,1,1,1,1,1,
    };
    Matrix<double> expected_score = Matrix(n_cols, n_rows, scores_arr.data()).transpose();
    Matrix<int> expected_row_memory= Matrix(n_cols, n_rows, row_mem_arr.data()).transpose();
    Matrix<int> expected_col_memory = Matrix(n_cols, n_rows, col_mem_arr.data()).transpose();

    assert_matrix_equals(dp.score_matrix, expected_score);
    assert_matrix_equals(dp.row_memory_matrix, expected_row_memory);
    assert_matrix_equals(dp.col_memory_matrix, expected_col_memory);
}

TEST_CASE("Aligner emits all candidate local alignments without filtering, and filters correctly with thresholds",
          "[aligner][sw][align_set]")
{
    // Same query/reference pair as the DP-matrix test above, so the expected candidate list here
    // can be derived directly from that test's hand-worked matrices/dp.candidates.
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const int gap_penalty = 13;
    const std::string query_read = "AAAAAT";
    const std::string germline_ref = "AAAAAAAAAT";
    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };

    auto aligner = make_legacy_aligner(matrix, gap_penalty, V_gene, genomic_templates); // D_gene -> full local mode

    SECTION("No score threshold, no offset bounds: every candidate is emitted")
    {
        const auto alignments = aligner.align_seq(query_read, -1000.0, /*best_align_only=*/false,
                                                  /*best_gene_only=*/false, INT16_MIN, INT16_MAX);
        assert_alignment_set_matches(
                alignments, query_read, genomic_templates,
                {
                        // TODO(user): one ExpectedAlignment{gene, core_cigar, extended_cigar, score}
                        // per candidate derived from dp.candidates in the DP-matrix test above.
                        ExpectedAlignment{ "g1", "4N6=", "4N6=", 42 },
                        ExpectedAlignment{ "g1", "3N5=2N1S", "3N5=1X1N", 35 },
                        ExpectedAlignment{ "g1", "2N5=3N1S", "2N5=1X2N", 35 },
                        ExpectedAlignment{ "g1", "1N5=4N1S", "1N5=1X3N", 35 },
                        ExpectedAlignment{ "g1", "5=5N1S", "5=1X4N", 35 },
                        ExpectedAlignment{ "g1", "5N4=1N2S", "5N4=1X1S", 28 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "6N3=1N3S", "6N3=1X2S", 21 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "7N2=1N4S", "7N2=1X3S", 14 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "1S4=6N1S", "1S4=1X5N", 15 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "8N1=1N5S", "8N1=1X4S", 7 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "2S3=7N1S", "2S3=1X6N", -5 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "3S2=8N1S", "3S2=1X7N", -25 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "4S1=9N1S", "4S1=1X8N", -45 }, // Watterman-Eggert issue?
                });
    }

    SECTION("Score threshold removes low-scoring candidates")
    {
        const double score_threshold = 20;
        const auto alignments = aligner.align_seq(query_read, score_threshold, false, false, INT16_MIN, INT16_MAX);
        assert_alignment_set_matches(
                alignments, query_read, genomic_templates,
                {
                        ExpectedAlignment{ "g1", "4N6=", "4N6=", 42 },
                        ExpectedAlignment{ "g1", "3N5=2N1S", "3N5=1X1N", 35 },
                        ExpectedAlignment{ "g1", "2N5=3N1S", "2N5=1X2N", 35 },
                        ExpectedAlignment{ "g1", "1N5=4N1S", "1N5=1X3N", 35 },
                        ExpectedAlignment{ "g1", "5=5N1S", "5=1X4N", 35 },
                        ExpectedAlignment{ "g1", "5N4=1N2S", "5N4=1X1S", 28 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "6N3=1N3S", "6N3=1X2S", 21 }, // Watterman-Eggert issue?
                });
    }

    SECTION("Offset bounds remove out-of-range candidates")
    {
        const int min_offset = -5;
        const int max_offset = 1;
        const auto alignments = aligner.align_seq(query_read, -1000.0, false, false, min_offset, max_offset);
        assert_alignment_set_matches(
                alignments, query_read, genomic_templates,
                {
                        // TODO(user): one ExpectedAlignment{gene, core_cigar, extended_cigar, score}
                        // per candidate derived from dp.candidates in the DP-matrix test above.
                        ExpectedAlignment{ "g1", "4N6=", "4N6=", 42 },
                        ExpectedAlignment{ "g1", "3N5=2N1S", "3N5=1X1N", 35 },
                        ExpectedAlignment{ "g1", "2N5=3N1S", "2N5=1X2N", 35 },
                        ExpectedAlignment{ "g1", "1N5=4N1S", "1N5=1X3N", 35 },
                        ExpectedAlignment{ "g1", "5=5N1S", "5=1X4N", 35 },
                        ExpectedAlignment{ "g1", "5N4=1N2S", "5N4=1X1S", 28 }, // Watterman-Eggert issue?
                        ExpectedAlignment{ "g1", "1S4=6N1S", "1S4=1X5N", 15 }, // Watterman-Eggert issue?
                });
    }
}

TEST_CASE("Reversed local alignment matches forward local alignment", "[aligner][sw][reverse]")
{
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const int gap_penalty = 13;
    const std::string query = "AAAAAT";
    const std::string reference = "AAAAAAAAAT";

    const Int_Str int_query = nt2int(query);
    const Int_Str int_reference = nt2int(reference);

    SwDPConfig forward_config;
    forward_config.score_threshold = -1000.0;
    forward_config.min_offset = INT16_MIN;
    forward_config.max_offset = INT16_MAX;
    forward_config.substitution_matrix = matrix;
    forward_config.gap_penalty = gap_penalty;
    forward_config.alignment_mode = SwAlignmentMode{ true, true, true, true, false };

    SwDPConfig reverse_config = forward_config;
    reverse_config.alignment_mode.reverse_sequences = true; // mirroring true/true/true/true is a no-op

    auto forward_alignments = sw_align(int_query, int_reference, /*best_only=*/false, forward_config);
    auto reverse_alignments = sw_align(int_query, int_reference, /*best_only=*/false, reverse_config);

    // Candidate order isn't part of the contract being tested here, so compare the two candidate
    // sets order-independently.
    assert_alignment_set_matches(forward_alignments, reverse_alignments, "forward", "reverse");
}

TEST_CASE("Dropping extended gaps must trigger failure of Alignment data comparison.",
          "[aligner][V_gene][legacy_csv][!shouldfail]")
{
    const Matrix<double> matrix = build_test_score_matrix(5.0, -14.0);
    const int gap_penalty = 30;
    std::string query_read;
    std::string germline_ref;
    std::string expected_csv_line;
    std::string obtained_csv_line;

    SECTION("Long match without V deletions")
    {
        query_read = "TCAGAACCCAGGGACTCAGCTGTGTATTTTTGTGCTAGTGGTTTGGTACAATCAGCCCCA";
        germline_ref =
                "gaagctggagttgcccagtctcccagatataagattatagagaaaaggcagagtgtggctttttggtgcaatcctatatctggccatgctaccctttactggtaccagcagatcctgggacagggcccaaagcttctgattcagtttcagaataacggtgtagtggatgattcacagttgcctaaggatcgattttctgcagagaggctcaaaggagtagactccactctcaagatccagcctgcaaagcttgaggactcggccgtgtatctctgtgccagcagcttaga";
        expected_csv_line = "1;g1;125;-238;{};{};{1,10,13,19,23,44,45,48};44;0;43";
        obtained_csv_line = "1;g1;125;-238;{};{};{1,10,13,19,23};44;0;43";
    }

    const std::vector<std::pair<std::string, std::string>> genomic_templates = { { "g1", germline_ref } };
    auto align = parse_single_alignment_csv_line(obtained_csv_line);
    assert_alignment_data_matches(align.second, expected_csv_line, query_read, genomic_templates);
}

// ============================================================================
//  Enhanced Alignment Data Representation Tests
// ============================================================================

TEST_CASE("Alignment_data getters", "[aligner][alignment_data][accessors]")
{
    SECTION("Default construction with lengths")
    {
        Alignment_data aln("test_gene", 10, 50, 100);
        REQUIRE(aln.gene_name == "test_gene");
        REQUIRE(aln.offset == 10);
        REQUIRE(aln.query_length == 50);
        REQUIRE(aln.germline_length == 100);
        REQUIRE(aln.query_align_start() == aln.five_p_offset);
        REQUIRE(aln.query_align_end() == aln.three_p_offset);
    }

    SECTION("Bounds getters with negative offset no in-dels")
    {
        //  qqQQqq
        // rrrRRrr
        Alignment_data aln("test", -1, 2, 3, 2, { }, { }, { }, 100.0, 6, 7);

        // Core alignment 
        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 3);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 4);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 5);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 6);
    }

    SECTION("Bounds getters with positive offset no in-dels")
    {
        // qqqqqQQqq
        //    rrRRrr
        Alignment_data aln("test", 3, 5, 6, 2, { }, { }, { }, 100.0, 9, 6);

        // Core alignment 
        REQUIRE(aln.query_align_start() == 5);
        REQUIRE(aln.query_align_end() == 6);
        REQUIRE(aln.reference_align_start() == 2);
        REQUIRE(aln.reference_align_end() == 3);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 3);
        REQUIRE(aln.extended_query_align_end() == 8);
        REQUIRE(aln.extended_reference_align_start() == 0);
        REQUIRE(aln.extended_reference_align_end() == 5);
    }

    SECTION("Bounds getters with negative offset, del in core")
    {
        //  qqQ-QQqq
        // rrrRRRRrr
        Alignment_data aln("test", -1, 2, 4, 4, { }, { 4 }, { }, 100.0, 7, 9);

        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 4);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 6);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 6);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 8);
    }

    SECTION("Bounds getters with positive offset, del in core")
    {
        // qqqqQQ-Qqq
        //   rrRRRRrr
        Alignment_data aln("test", 2, 4, 6, 4, { }, { 4 }, { }, 100.0, 9, 8);

        // Core alignment
        REQUIRE(aln.query_align_start() == 4);
        REQUIRE(aln.query_align_end() == 6);
        REQUIRE(aln.reference_align_start() == 2);
        REQUIRE(aln.reference_align_end() == 5);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 2);
        REQUIRE(aln.extended_query_align_end() == 8);
        REQUIRE(aln.extended_reference_align_start() == 0);
        REQUIRE(aln.extended_reference_align_end() == 7);
    }

    SECTION("Bounds getters with negative offset, del in 5p extension")
    {
        //  q-QQQQq
        // rrrRRRRr
        Alignment_data aln("test", -1, 1, 4, 4, { }, { 2 }, { }, 100.0, 6, 8);

        // Core alignment
        REQUIRE(aln.query_align_start() == 1);
        REQUIRE(aln.query_align_end() == 4);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 6);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 5);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 7);
    }

    SECTION("Bounds getters with positive offset, del in 5p extension")
    {
        // qqq-QQQQqqq
        //   rrRRRRrrr
        Alignment_data aln("test", 2, 3, 6, 4, { }, { 1 }, { }, 100.0, 10, 9);

        // Core alignment
        REQUIRE(aln.query_align_start() == 3);
        REQUIRE(aln.query_align_end() == 6);
        REQUIRE(aln.reference_align_start() == 2);
        REQUIRE(aln.reference_align_end() == 5);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 2);
        REQUIRE(aln.extended_query_align_end() == 9);
        REQUIRE(aln.extended_reference_align_start() == 0);
        REQUIRE(aln.extended_reference_align_end() == 8);
    }

    SECTION("Bounds getters with negative offset, insertion in core")
    {
        //  qqQQQQ
        // rrrR-RR
        Alignment_data aln("test", -1, 2, 5, 4, { 3 }, { }, { }, 100.0, 6, 6);

        // Core alignment
        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 5);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 5);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 5);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 5);
    }

    SECTION("Bounds getters with positive offset, insertion in core")
    {
        // qqqqQQQQqqq
        //   rrRR-Rrrr
        Alignment_data aln("test", 2, 4, 7, 4, { 6 }, { }, { }, 100.0, 11, 8);

        // Core alignment
        REQUIRE(aln.query_align_start() == 4);
        REQUIRE(aln.query_align_end() == 7);
        REQUIRE(aln.reference_align_start() == 2);
        REQUIRE(aln.reference_align_end() == 4);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 2);
        REQUIRE(aln.extended_query_align_end() == 10);
        REQUIRE(aln.extended_reference_align_start() == 0);
        REQUIRE(aln.extended_reference_align_end() == 7);
    }

    SECTION("Bounds getters with negative offset, insertion in 5p extension")
    {
        //  qqQQQQqqqq
        // rr-RRRRqqqq
        Alignment_data aln("test", -1, 2, 5, 4, { 1 }, { }, { }, 100.0, 10, 10);

        // Core alignment
        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 5);
        REQUIRE(aln.reference_align_start() == 2);
        REQUIRE(aln.reference_align_end() == 5);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 9);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 9);
    }

    SECTION("Bounds getters with positive offset, insertion in 5p extension")
    {
        // qqqqQQQQqq
        //   r-RRRRqq
        Alignment_data aln("test", 2, 4, 7, 4, { 1 }, { }, { }, 100.0, 10, 7);

        // Core alignment
        REQUIRE(aln.query_align_start() == 4);
        REQUIRE(aln.query_align_end() == 7);
        REQUIRE(aln.reference_align_start() == 1);
        REQUIRE(aln.reference_align_end() == 4);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 2);
        REQUIRE(aln.extended_query_align_end() == 9);
        REQUIRE(aln.extended_reference_align_start() == 0);
        REQUIRE(aln.extended_reference_align_end() == 6);
    }

    SECTION("Bounds getters with extended clipped by reference length")
    {
        //  qqQQqqqq
        // rrrRRrr
        Alignment_data aln("test", -1, 2, 3, 2, { }, { }, { }, 100.0, 8, 7);

        // Core alignment 
        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 3);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 4);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 5);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 6);
    }

    SECTION("Bounds getters with extended clipped by query length")
    {
        //  qqQQqq
        // rrrRRrrrr
        Alignment_data aln("test", -1, 2, 3, 2, { }, { }, { }, 100.0, 6, 9);

        // Core alignment 
        REQUIRE(aln.query_align_start() == 2);
        REQUIRE(aln.query_align_end() == 3);
        REQUIRE(aln.reference_align_start() == 3);
        REQUIRE(aln.reference_align_end() == 4);

        // Extended alignment
        REQUIRE(aln.extended_query_align_start() == 0);
        REQUIRE(aln.extended_query_align_end() == 5);
        REQUIRE(aln.extended_reference_align_start() == 1);
        REQUIRE(aln.extended_reference_align_end() == 6);
    }

}

// ============================================================================
//  Deletion Categorization Tests
// ============================================================================

TEST_CASE("Alignment_data deletion categorization", "[aligner][alignment_data][deletions]")
{
    SECTION("Single 5p extended deletion")
    {
        //  q-qQQQQ
        // rrrrRRRR
        // Deletion at ref_pos = 2 is before ref_start, so it's 5p extended
        Alignment_data aln("test", -1, 2, 5, 4, { }, { 2 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 2 }));
        REQUIRE(core_del.empty());
        REQUIRE(three_p_ext_del.empty());
    }

    SECTION("Single core deletion")
    {
        //  qqQ-QQ
        // rrrRRRR
        // Deletion at ref_pos = 4 is within [3, 5], so it's core
        Alignment_data aln("test", -1, 2, 4, 4, { }, { 4 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE(ext_del.empty());
        REQUIRE(core_del.size() == 1);
        REQUIRE(core_del[0] == 4);
        REQUIRE(three_p_ext_del.empty());
    }

    SECTION("Single 3p extended deletion")
    {
        //  qqQQQQq-qq
        // rrrRRRRrrrr
        // Deletion at ref_pos = 6 is after ref_end, so it's 3p extended
        Alignment_data aln("test", -1, 2, 4, 4, { }, { 8 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE(ext_del.empty());
        REQUIRE(core_del.empty());
        REQUIRE_THAT(three_p_ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 8 }));
    }

    SECTION("Multiple 5p extended deletions")
    {
        // q----qqqqQQQQQQqqq
        // rrrrrrrrrRRRRRRrrr
        Alignment_data aln("test", 0, 5, 10, 10, { }, { 1, 2, 3, 4 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 1, 2, 3, 4 }));
        REQUIRE(core_del.empty());
        REQUIRE(three_p_ext_del.empty());
    }

    SECTION("Mixed deletions: 5p extended, core, 3p extended")
    {
        // q--qq-q-qQ-Q-QQQQ-qqqq-qq
        // rrrrrrrrrRRRRRRRRrrrrrrrr
        Alignment_data aln("test", 0, 5, 10, 10, { }, { 1, 2, 5, 7, 10, 12, 17, 22 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 1, 2, 5, 7 }));
        REQUIRE_THAT(core_del, Catch::Matchers::Equals(std::vector<size_t>{ 10, 12 }));
        REQUIRE_THAT(three_p_ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 17, 22 }));
    }

    SECTION("Deletion at boundary between 5p and core")
    {
        // qq-QQQQ
        // rrrRRRR
        Alignment_data aln("test", 0, 2, 5, 4, { }, { 2 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 2 }));
        REQUIRE(core_del.empty());
    }

    SECTION("Deletion just after core end")
    {
        // qqqQQQQ-q
        // rrrRRRRrr
        Alignment_data aln("test", 0, 3, 6, 4, { }, { 7 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE(ext_del.empty());
        REQUIRE(core_del.empty());
        REQUIRE_THAT(three_p_ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 7 }));
    }

    SECTION("Positive offset with 5p extended deletion")
    {
        // qq-qQQQQ
        //   rrRRRR
        Alignment_data aln("test", 2, 3, 6, 4, { }, { 0 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 0 }));
        REQUIRE(core_del.empty());
    }

    SECTION("Negative offset with multiple deletions")
    {
        //    qq-q-qQ-Q-QQQQ-qqqq-qq
        // rrrrrrrrrRRRRRRRRrrrrrrrr
        Alignment_data aln("test", -3, 4, 9, 10, { }, { 5, 7, 10, 12, 17, 22 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 5, 7 }));
        REQUIRE_THAT(core_del, Catch::Matchers::Equals(std::vector<size_t>{ 10, 12 }));
        REQUIRE_THAT(three_p_ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 17, 22 }));
    }

    SECTION("Negative offset with multiple deletions and insertions")
    {
        //    qq-qqq-qQ-Q-QQQQQ-qqqqq-qq
        // rrrrrrr--rrRRRRRR-RRrrr-rrrrr
        Alignment_data aln("test", -3, 6, 12, 9, { 3, 4, 10, 15 }, { 5, 7, 10, 12, 17, 22 }, { }, 100.0);

        auto ext_del = aln.get_5p_extended_deletions();
        auto core_del = aln.get_core_deletions();
        auto three_p_ext_del = aln.get_3p_extended_deletions();

        REQUIRE_THAT(ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 5, 7 }));
        REQUIRE_THAT(core_del, Catch::Matchers::Equals(std::vector<size_t>{ 10, 12 }));
        REQUIRE_THAT(three_p_ext_del, Catch::Matchers::Equals(std::vector<size_t>{ 17, 22 }));
    }
}

// ============================================================================
//  Insertion Categorization Tests
// ============================================================================

TEST_CASE("Alignment_data insertion categorization", "[aligner][alignment_data][insertions]")
{
    SECTION("Single 5p extended insertion")
    {
        //  qqQQQQ
        // r-rRRRR
        // Deletion at ref_pos = 2 is before ref_start, so it's 5p extended
        Alignment_data aln("test", -1, 2, 5, 4, { 0 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE_THAT(five_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 0 }));
        REQUIRE(core_ins.empty());
        REQUIRE(three_p_ext_ins.empty());
    }

    SECTION("Single core deletion")
    {
        //  qqQQQQ
        // rrrR-RR
        // Deletion at ref_pos = 4 is within [3, 5], so it's core
        Alignment_data aln("test", -1, 2, 5, 4, { 3 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE(five_p_ext_ins.empty());
        REQUIRE_THAT(core_ins, Catch::Matchers::Equals(std::vector<size_t>{ 3 }));
        REQUIRE(three_p_ext_ins.empty());
    }

    SECTION("Single 3p extended deletion")
    {
        //  qqQQQQqqqq
        // rrrRRRRrr-r
        // Deletion at ref_pos = 6 is after ref_end, so it's 3p extended
        Alignment_data aln("test", -1, 2, 5, 4, { 8 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE(five_p_ext_ins.empty());
        REQUIRE(core_ins.empty());
        REQUIRE_THAT(three_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 8 }));
    }

    SECTION("Multiple 5p extended deletions")
    {
        // qqqqqqqqqQQQQQQqqq
        // r----rrrrRRRRRRrrr
        Alignment_data aln("test", 0, 5, 10, 10, { 1, 2, 3, 4 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE_THAT(five_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 1, 2, 3, 4 }));
        REQUIRE(core_ins.empty());
        REQUIRE(three_p_ext_ins.empty());
    }

    SECTION("Mixed deletions: 5p extended, core, 3p extended")
    {
        // qqqqqqqqqQQQQQQQQqqqqqqqq
        // r--rr-r-rR-R-RRRR-rrrr-rr
        Alignment_data aln("test", 0, 9, 16, 8, { 1, 2, 5, 7, 10, 12, 17, 22 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE_THAT(five_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 1, 2, 5, 7 }));
        REQUIRE_THAT(core_ins, Catch::Matchers::Equals(std::vector<size_t>{ 10, 12 }));
        REQUIRE_THAT(three_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 17, 22 }));
    }

    SECTION("Insertion just before core start")
    {
        // qqqqQQQQ
        //  rr-RRRR
        Alignment_data aln("test", 1, 4, 7, 4, { 3 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE_THAT(five_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 3 }));
        REQUIRE(core_ins.empty());
        REQUIRE(three_p_ext_ins.empty());
    }

    SECTION("Insertion just after core end")
    {
        // qqqqQQQQqqq
        //  rrrRRRR-rr
        Alignment_data aln("test", 0, 4, 7, 10, { 8 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE(five_p_ext_ins.empty());
        REQUIRE(core_ins.empty());
        REQUIRE_THAT(three_p_ext_ins, Catch::Matchers::Equals(std::vector<size_t>{ 8 }));
    }

    SECTION("Insertion just after core start")
    {
        // qqqqQQQQ
        //  rrr-RRR
        Alignment_data aln("test", 1, 4, 7, 4, { 4 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE(five_p_ext_ins.empty());
        REQUIRE_THAT(core_ins, Catch::Matchers::Equals(std::vector<size_t>{ 4 }));
        REQUIRE(three_p_ext_ins.empty());
    }

    SECTION("Insertion just before core end")
    {
        // qqqqQQQQqqq
        //  rrrRRR-rrr
        Alignment_data aln("test", 0, 5, 10, 10, { 7 }, { }, { }, 100.0);

        auto five_p_ext_ins = aln.get_5p_extended_insertions();
        auto core_ins = aln.get_core_insertions();
        auto three_p_ext_ins = aln.get_3p_extended_insertions();

        REQUIRE(five_p_ext_ins.empty());
        REQUIRE_THAT(core_ins, Catch::Matchers::Equals(std::vector<size_t>{ 7 }));
        REQUIRE(three_p_ext_ins.empty());
    }
}