#include <catch2/catch_test_macros.hpp>

#include <igor/Core/Aligner.h>
#include "AlignerTestUtils.h"

#include <algorithm>
#include <cstdio>
#include <forward_list>
#include <string>
#include <vector>

using namespace igor::test::align;

TEST_CASE("parse extended CIGAR", "[cigar]")
{
    auto ops = parse_cigar("3=1X2I4D5M");
    REQUIRE(ops.size() == 5);
    REQUIRE(ops[0] == std::make_pair(3, '='));
    REQUIRE(ops[4] == std::make_pair(5, 'M'));
    REQUIRE_THROWS_AS(parse_cigar("0="), std::invalid_argument);
    REQUIRE_THROWS_AS(parse_cigar("=10"), std::invalid_argument);
    REQUIRE_THROWS_AS(parse_cigar("10Z"), std::invalid_argument);
}

TEST_CASE("alignment data to CIGAR basic cases", "[cigar]")
{
    REQUIRE(alignment_data_to_core_cigar(Alignment_data("g", 0, 0, 9, 10, { }, { }, { }, 10)) == "10=");
    REQUIRE(alignment_data_to_core_cigar(
                    Alignment_data("g", 0, 0, 9, 10, { }, { }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, 10))
            == "10X");
    REQUIRE(alignment_data_to_core_cigar(Alignment_data("g", 0, 0, 12, 13, { 7, 6, 5 }, { }, { }, 10)) == "5=3I5=");
    REQUIRE(alignment_data_to_core_cigar(Alignment_data("g", 0, 0, 9, 12, { }, { 5, 4 }, { }, 10)) == "4=2D6=");
}

TEST_CASE("alignment data full-span CIGAR includes terminal gaps", "[cigar]")
{
    Alignment_data aln("gene", -237, 0, 52, 53, { 50, 51 }, { }, { }, 246.0);
    REQUIRE(alignment_data_to_core_cigar(aln) == "50=2I1=");
    REQUIRE(alignment_data_to_core_cigar(aln, 60, 288) == "237N50=2I1=7S");
}

TEST_CASE("alignment data CIGAR mixed round trip", "[cigar]")
{
    Alignment_data aln("gene", -2, 0, 9, 12, { 6 }, { 3, 4 }, { 2, 8 }, 42.0);
    std::string cigar = alignment_data_to_core_cigar(aln);
    REQUIRE(cigar == "1=2D1=1X3=1I1=1X1=");

    Alignment_data round = alignment_data_from_cigar(
            aln.gene_name, cigar, alignment_data_sequence_start(aln), alignment_data_sequence_end(aln),
            alignment_data_germline_start(aln), alignment_data_germline_end(aln), aln.score);
    REQUIRE(round.gene_name == aln.gene_name);
    REQUIRE(round.offset == aln.offset);
    REQUIRE(round.five_p_offset == aln.five_p_offset);
    REQUIRE(round.three_p_offset == aln.three_p_offset);
    REQUIRE(round.align_length == aln.align_length);
    REQUIRE(sorted_list(round.insertions) == sorted_list(aln.insertions));
    REQUIRE(sorted_list(round.deletions) == sorted_list(aln.deletions));
    std::sort(round.mismatches.begin(), round.mismatches.end());
    auto expected_mismatches = aln.mismatches;
    std::sort(expected_mismatches.begin(), expected_mismatches.end());
    REQUIRE(round.mismatches == expected_mismatches);
}

TEST_CASE("CIGAR M fallback loses mismatch detail by design", "[cigar]")
{
    Alignment_data aln = alignment_data_from_cigar("gene", "5M", 1, 5, 1, 5, 10);
    REQUIRE(aln.align_length == 5);
    REQUIRE(aln.mismatches.empty());
    REQUIRE(aln.five_p_offset == 0);
    REQUIRE(aln.three_p_offset == 4);
}

TEST_CASE("CSV reader preserves alignment length and offsets", "[cigar]")
{
    std::string path = "alignment_cigar_tmp.csv";
    Aligner aligner;
    std::unordered_map<int, std::forward_list<Alignment_data>> alignments;
    Alignment_data aln("gene", 1, 1, 10, 11, { 6 }, { }, { 4 }, 12.0);
    alignments[7].push_front(aln);
    aligner.write_alignments_seq_csv(path, alignments);

    auto loaded = aligner.read_alignments_seq_csv(path, 0.0, true);
    std::remove(path.c_str());
    REQUIRE(loaded.count(7) == 1);
    REQUIRE(!loaded[7].empty());
    const Alignment_data &got = loaded[7].front();
    REQUIRE(got.align_length == aln.align_length);
    REQUIRE(got.five_p_offset == aln.five_p_offset);
    REQUIRE(got.three_p_offset == aln.three_p_offset);
    REQUIRE(alignment_data_to_core_cigar(got) == alignment_data_to_core_cigar(aln));
}

TEST_CASE("alignment data to extended CIGAR basic cases", "[cigar]")
{
    // Basic cases without leading or trailing gaps, and no extended mismatches
    // i.e. extended and core CIGAR are the same

    // Test case 1: Complete alignment with no gaps
    {
        Alignment_data aln("g", 0, 0, 9, 10, { }, { }, { }, 10);
        std::string expected_core = "10=";
        std::string expected_extended = expected_core; // No end gaps, so same as core
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }

    // Test case 2: Core alignment with mismatches
    {
        Alignment_data aln("g", 0, 0, 9, 10, { }, { }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, 10);
        std::string expected_core = "10X";
        std::string expected_extended = expected_core; // No end gaps, so same as core
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }

    // Test case 3: Core alignment with insertions
    {
        Alignment_data aln("g", 0, 0, 12, 13, { 7, 6, 5 }, { }, { }, 10);
        std::string expected_core = "5=3I5=";
        std::string expected_extended = expected_core; // No end gaps, so same as core
        REQUIRE(alignment_data_to_core_cigar(aln, 13, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 13, 10) == expected_extended);
    }

    // Test case 4: Core alignment with deletions
    {
        Alignment_data aln("g", 0, 0, 9, 12, { }, { 5, 4 }, { }, 10);
        std::string expected_core = "4=2D6=";
        std::string expected_extended = expected_core; // No end gaps, so same as core
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 12) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 12) == expected_extended);
    }
}

TEST_CASE("alignment data to extended CIGAR with end gaps", "[cigar]")
{
    // Basic cases with non extendable leading or trailing gaps, and no extended mismatches
    // i.e. extended and core CIGAR are the same
    
    // Test case 1: Leading reference gaps (N operators)
    {
        Alignment_data aln("gene", -5, 0, 10, 11, { }, { }, { }, 10);
        // offset = -5 means germline starts 5 positions before query
        // five_p_offset = 0, three_p_offset = 10, so alignment covers query positions 0-10 (11 positions)
        // germline starts at position 0 - (-5) = 5, so there are 5 leading reference gaps
        std::string expected_core = "5N11="; // 5 leading reference gaps + 11 matches
        std::string expected_extended = expected_core; 
        REQUIRE(alignment_data_to_core_cigar(aln, 11, 16) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 11, 16) == expected_extended);
    }

    // Test case 2: Trailing reference gaps (N operators)
    {
        Alignment_data aln("gene", 0, 0, 10, 11, { }, { }, { }, 10);
        // offset = 0, five_p_offset = 0, three_p_offset = 10
        // query length = 11, germline length = 16
        // alignment covers query positions 0-10 and germline positions 0-10
        // there should be 5 trailing reference gaps (16-11=5)
        std::string expected_core = "11=5N";  // 11 matches + 5 trailing reference gaps
        std::string expected_extended = expected_core; 
        REQUIRE(alignment_data_to_core_cigar(aln, 11, 16) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 11, 16) == expected_extended);
    }

    // Test case 3: Leading query gaps (S operators)
    {
        Alignment_data aln("gene", 5, 5, 15, 11, { }, { }, { }, 10);
        // five_p_offset = 5, so there are 5 leading query gaps
        std::string expected_core = "5S11="; // 5 leading query gaps + 11 matches
        std::string expected_extended = expected_core;
        REQUIRE(alignment_data_to_core_cigar(aln, 16, 11) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 16, 11) == expected_extended);
    }

    // Test case 4: Trailing query gaps (S operators)
    {
        Alignment_data aln("gene", 0, 0, 10, 11, { }, { }, { }, 10);
        // query length = 16, alignment ends at position 10
        // so there are 5 trailing query gaps
        std::string expected_core = "11=5S"; // 11 matches + 5 trailing query gaps
        std::string expected_extended = expected_core;
        REQUIRE(alignment_data_to_core_cigar(aln, 16, 11) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 16, 11) == expected_extended);
    }
}

TEST_CASE("alignment data to extended CIGAR with extended mismatches", "[cigar]")
{
    // Only matches in leading extended region
    {
        Alignment_data aln("gene", 0, 2, 9, 8, { }, { }, { }, 10);
        // No mismatches 
        // offset = 0 & five_p_offset = 2, implies 2 leading five prime deletions
        std::string expected_core = "2N2S8=";
        std::string expected_extended = "10=";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }
    //  Mismatches in leading extended region
    {
        Alignment_data aln("gene", 0, 2, 9, 8, { }, { }, { 0, 1 }, 10);
        // Mismatches at positions 0,1 which are masked by leading deletions
        // offset = 0 & five_p_offset = 2, implies 2 leading five prime deletions
        std::string expected_core = "2N2S8="; // 11 matches + 5 trailing query gaps
        std::string expected_extended = "2X8=";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }
    //  Mixed case with negative offset implying remaining leading deletion
    {
        Alignment_data aln("gene", -2, 2, 9, 8, { }, { }, { 1 }, 10);
        // Mismatches at positions 0,1 which are masked by leading deletions
        // offset = -2 & five_p_offset = 2, implies 4 leading five prime deletions
        std::string expected_core = "4N2S8="; // 11 matches + 5 trailing query gaps
        std::string expected_extended = "2N1=1X8=";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 12) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 12) == expected_extended);
    }

    // Only matches in trailing extended region
    {
        Alignment_data aln("gene", 0, 0, 7, 8, { }, { }, { }, 10);
        // No mismatches 
        // offset = 0 & five_p_offset = 2, implies 2 leading five prime deletions
        std::string expected_core = "8=2N2S";
        std::string expected_extended = "10=";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }
    //  Mismatches in trailing extended region
    {
        Alignment_data aln("gene", 0, 0, 7, 8, { }, { }, { 8,9 }, 10);
        // Mismatches at positions 0,1 which are masked by leading deletions
        // offset = 0 & five_p_offset = 2, implies 2 leading five prime deletions
        std::string expected_core = "8=2N2S"; // 11 matches + 5 trailing query gaps
        std::string expected_extended = "8=2X";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 10) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 10) == expected_extended);
    }
    //  Mixed case with longer reference implying remaining trailing deletion
    {
        Alignment_data aln("gene", 0, 0, 7, 8, { }, { }, { 8 }, 10);
        // Mismatches at positions 0,1 which are masked by leading deletions
        // offset = -2 & five_p_offset = 2, implies 4 leading five prime deletions
        std::string expected_core = "8=4N2S"; // 11 matches + 5 trailing query gaps
        std::string expected_extended = "8=1X1=2N";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 12) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 12) == expected_extended);
    }
    // Both mixed cases at once
    {
        Alignment_data aln("gene", -2, 2, 7, 6, { }, { }, { 1, 8 }, 10);
        // Mismatches at positions 0,1 which are masked by leading deletions
        // offset = -2 & five_p_offset = 2, implies 4 leading five prime deletions
        std::string expected_core = "4N2S6=4N2S"; // 11 matches + 5 trailing query gaps
        std::string expected_extended = "2N1=1X6=1X1=2N";
        REQUIRE(alignment_data_to_core_cigar(aln, 10, 14) == expected_core);
        REQUIRE(alignment_data_to_extended_cigar(aln, 10, 14) == expected_extended);
    }
}

TEST_CASE("alignment data to extended CIGAR edge cases", "[cigar]")
{
    // Test case 1: No alignment (edge case)
    SECTION("No alignment (edge case)")
    {
        Alignment_data aln("gene", 0, 0, 0, 0, { }, { }, { }, 0);
        std::string extended_cigar = alignment_data_to_extended_cigar(aln, 10, 10);
        INFO("No alignment - Extended CIGAR: " << extended_cigar);
        // Should produce CIGAR with all gaps
        REQUIRE(!extended_cigar.empty());
    }

    // Test case 3: Very short sequences
    SECTION("Very short sequences")
    {
        Alignment_data aln("gene", 0, 0, 0, 1, { }, { }, { }, 5);
        std::string extended_cigar = alignment_data_to_extended_cigar(aln, 1, 1);
        INFO("Short sequences - Extended CIGAR: " << extended_cigar);
        REQUIRE(!extended_cigar.empty());
    }

    // Test case for D gene alignment with extended mismatches
    SECTION("D gene alignment with extended mismatches")
    {
        // offset=0, five_p_offset=4, three_p_offset=14, align_length=11
        // empty insertions and deletions, mismatches at positions 1, 3, 8, 13
        // query and germline lengths are both 16
        Alignment_data aln("g", 0, 4, 14, 11, { }, { }, { 1, 3, 8, 13 }, 0.0);
        std::string extended_cigar = alignment_data_to_extended_cigar(aln, 16, 16);
        std::string expected = "1=1X1=1X4=1X4=1X2=";
        INFO("D gene extended CIGAR: " << extended_cigar);
        REQUIRE(extended_cigar == expected);
    }

    // Test case for D gene alignment without extended mismatches (local alignment)
    SECTION("D gene local alignment without extended mismatches")
    {
        // Query: TTTTAAAAAATTTT (14 chars)
        // Germline: AAAAAA (6 chars)
        // offset=4, five_p_offset=4, three_p_offset=9, align_length=6
        // empty insertions, deletions, and mismatches
        // query_length=14, germline_length=6
        // Expected: 4S6=4S (same as core CIGAR)
        Alignment_data aln("g", 4, 4, 9, 6, { }, { }, { }, 0.0);
        std::string extended_cigar = alignment_data_to_extended_cigar(aln, 14, 6);
        std::string expected = "4S6=4S";
        INFO("D gene local alignment extended CIGAR: " << extended_cigar);
        REQUIRE(extended_cigar == expected);
    }
}

TEST_CASE("CIGAR to alignment data basic cases", "[cigar]")
{
    Alignment_data expected_aln("dummy", 0, 0, 0, 0, { }, { }, { }, 0);
    std::string core_cigar;
    std::string extended_cigar;
    std::string gene_name;
    double score;

    SECTION("Complete alignment with no gaps")
    {
        expected_aln = Alignment_data("g", 0, 0, 9, 10, { }, { }, { }, 10);
        core_cigar = "10=";
        extended_cigar = core_cigar;
        gene_name = "g";
        score = 10;
    }

    SECTION("Core alignment with mismatches")
    {
        expected_aln = Alignment_data("g", 0, 0, 9, 10, { }, { }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, 10);
        core_cigar = "10X";
        extended_cigar = core_cigar;
        gene_name = "g";
        score = 10;
    }

    SECTION("Core alignment with insertions")
    {
        expected_aln = Alignment_data("g", 0, 0, 12, 13, { 7, 6, 5 }, { }, { }, 10);
        core_cigar = "5=3I5=";
        extended_cigar = core_cigar;
        gene_name = "g";
        score = 10;
    }

    SECTION("Core alignment with deletions")
    {
        expected_aln = Alignment_data("g", 0, 0, 9, 12, { }, { 5, 4 }, { }, 10);
        core_cigar = "4=2D6=";
        extended_cigar = core_cigar;
        gene_name = "g";
        score = 10;
    }

    auto aln = alignment_data_from_cigar_and_extended(gene_name, core_cigar, extended_cigar, score);
    check_alignment_data_equal(expected_aln, aln);
}

TEST_CASE("CIGAR to alignment data with end gaps", "[cigar]")
{
    Alignment_data expected_aln("dummy", 0, 0, 0, 0, { }, { }, { }, 0);
    std::string core_cigar;
    std::string extended_cigar;
    std::string gene_name;
    double score;

    SECTION("Leading reference gaps")
    {
        expected_aln = Alignment_data("gene", -5, 0, 10, 11, { }, { }, { }, 10);
        core_cigar = "5N11=";
        extended_cigar = core_cigar;
        gene_name = "gene";
        score = 10;
    }

    SECTION("Trailing reference gaps")
    {
        expected_aln = Alignment_data("gene", 0, 0, 10, 11, { }, { }, { }, 10);
        core_cigar = "11=5N";
        extended_cigar = core_cigar;
        gene_name = "gene";
        score = 10;
    }

    SECTION("Leading query gaps")
    {
        expected_aln = Alignment_data("gene", 5, 5, 15, 11, { }, { }, { }, 10);
        core_cigar = "5S11=";
        extended_cigar = core_cigar;
        gene_name = "gene";
        score = 10;
    }

    SECTION("Trailing query gaps")
    {
        expected_aln = Alignment_data("gene", 0, 0, 10, 11, { }, { }, { }, 10);
        core_cigar = "11=5S";
        extended_cigar = core_cigar;
        gene_name = "gene";
        score = 10;
    }

    auto aln = alignment_data_from_cigar_and_extended(gene_name, core_cigar, extended_cigar, score);
    check_alignment_data_equal(expected_aln, aln);
}

TEST_CASE("CIGAR to alignment data with extended mismatches", "[cigar]")
{
    Alignment_data expected_aln("dummy", 0, 0, 0, 0, { }, { }, { }, 0);
    std::string core_cigar;
    std::string extended_cigar;
    std::string gene_name;
    double score;

    SECTION("Only matches in leading extended region")
    {
        expected_aln = Alignment_data("gene", 0, 2, 9, 8, { }, { }, { }, 10);
        core_cigar = "2N2S8=";
        extended_cigar = "10=";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Mismatches in leading extended region")
    {
        expected_aln = Alignment_data("gene", 0, 2, 9, 8, { }, { }, { 0, 1 }, 10);
        core_cigar = "2N2S8=";
        extended_cigar = "2X8=";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Mixed case with negative offset implying remaining leading deletion")
    {
        expected_aln = Alignment_data("gene", -2, 2, 9, 8, { }, { }, { 1 }, 10);
        core_cigar = "4N2S8=";
        extended_cigar = "2N1=1X8=";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Only matches in trailing extended region")
    {
        expected_aln = Alignment_data("gene", 0, 0, 7, 8, { }, { }, { }, 10);
        core_cigar = "8=2N2S";
        extended_cigar = "10=";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Mismatches in trailing extended region")
    {
        expected_aln = Alignment_data("gene", 0, 0, 7, 8, { }, { }, { 8, 9 }, 10);
        core_cigar = "8=2N2S";
        extended_cigar = "8=2X";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Mixed case with longer reference implying remaining trailing deletion")
    {
        expected_aln = Alignment_data("gene", 0, 0, 7, 8, { }, { }, { 8 }, 10);
        core_cigar = "8=4N2S";
        extended_cigar = "8=1X1=2N";
        gene_name = "gene";
        score = 10;
    }

    SECTION("Both mixed cases at once")
    {
        expected_aln = Alignment_data("gene", -2, 2, 7, 6, { }, { }, { 1, 8 }, 10);
        core_cigar = "4N2S6=4N2S";
        extended_cigar = "2N1=1X6=1X1=2N";
        gene_name = "gene";
        score = 10;
    }

    
    REQUIRE_THROWS_AS(alignment_data_from_cigar_and_extended(gene_name, core_cigar, core_cigar, score),  std::invalid_argument);
    auto extend_only_aln = alignment_data_from_cigar_and_extended(gene_name, extended_cigar, extended_cigar, score);
    auto aln = alignment_data_from_cigar_and_extended(gene_name, core_cigar, extended_cigar, score);
    REQUIRE(!alignment_data_equal(extend_only_aln, aln));
    check_alignment_data_equal(expected_aln, aln);
}
