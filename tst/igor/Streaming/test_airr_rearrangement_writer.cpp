/*
 * test_airr_rearrangement_writer.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  Unit tests for AIRR writer functions
 *
 *  Uses Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures
 *  - SECTION for grouping related tests
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "StreamingTestUtils.h"
#include <igor/Streaming/AIRRRearrangementWriter.h>
#include <igor/Streaming/AIRRRearrangementReader.h>

#include <fstream>
#include <sstream>

using namespace igor;
using namespace igor::airr::rearrangement;
using namespace igor::test;
using namespace Catch::Matchers;

// Import shared utility functions from parent namespace
using igor::airr::delimiter_char;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for AIRR writer tests with automatic cleanup
 */
struct AIRRWriterFixture
{
    TestDirectory test_dir{"igor_airr_writer_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }

    /**
     * @brief Create test sequences with V/D/J alignments
     */
    std::vector<SequenceData> create_test_sequences()
    {
        std::vector<SequenceData> sequences;

        // Sequence 1: Full V/D/J
        {
            std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;

            std::forward_list<int> insertions;
            std::forward_list<int> deletions;
            std::vector<int> mismatches;

            aligns[V_gene].emplace_back(
                "IGHV1-2*01", 0, 0, 0, 50, insertions, deletions, mismatches, 150.5);
            aligns[D_gene].emplace_back(
                "IGHD1-1*01", 50, 0, 0, 20, insertions, deletions, mismatches, 45.2);
            aligns[J_gene].emplace_back(
                "IGHJ3*01", 70, 0, 0, 30, insertions, deletions, mismatches, 120.0);

            sequences.emplace_back(1, "ATCGATCGATCGATCGATCGATCGATCGATCG", aligns);
        }

        // Sequence 2: V and J only (no D)
        {
            std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;

            std::forward_list<int> insertions;
            std::forward_list<int> deletions;
            std::vector<int> mismatches;

            aligns[V_gene].emplace_back(
                "IGHV3-9*01", 0, 0, 0, 45, insertions, deletions, mismatches, 140.0);
            aligns[J_gene].emplace_back(
                "IGHJ4*02", 45, 0, 0, 25, insertions, deletions, mismatches, 115.0);

            sequences.emplace_back(2, "GCTAGCTAGCTAGCTAGCTAGCTA", aligns);
        }

        // Sequence 3: No alignments
        {
            std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;
            sequences.emplace_back(3, "TTAATTAATTAATTAA", aligns);
        }

        return sequences;
    }

    /**
     * @brief Create sequences with insertions/deletions for CIGAR tests
     */
    std::vector<SequenceData> create_sequences_with_indels()
    {
        std::vector<SequenceData> sequences;

        std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;

        std::forward_list<int> insertions = {5, 10};  // 2 insertions
        std::forward_list<int> deletions = {15};       // 1 deletion
        std::vector<int> mismatches = {3, 8, 20};

        aligns[V_gene].emplace_back(
            "IGHV1-2*01", 0, 0, 0, 50, insertions, deletions, mismatches, 140.0);

        sequences.emplace_back(1, "ATCGATCGATCGATCG", aligns);

        return sequences;
    }

    /**
     * @brief Read file content as string
     */
    std::string read_file_content(const std::string& filepath)
    {
        std::ifstream file(filepath);
        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }
};

//==============================================================================
// write_tsv tests
//==============================================================================

TEST_CASE_METHOD(AIRRWriterFixture, "write_tsv", "[airr][writer]")
{
    SECTION("writes valid AIRR TSV")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("test.tsv");

        write_tsv(filepath, sequences);

        // Verify file exists and has content
        auto content = read_file_content(filepath);
        REQUIRE_FALSE(content.empty());

        // Check header
        REQUIRE_THAT(content, ContainsSubstring("sequence_id"));
        REQUIRE_THAT(content, ContainsSubstring("sequence"));
        REQUIRE_THAT(content, ContainsSubstring("v_call"));
        REQUIRE_THAT(content, ContainsSubstring("d_call"));
        REQUIRE_THAT(content, ContainsSubstring("j_call"));
    }

    SECTION("writes correct number of rows")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("test_rows.tsv");

        write_tsv(filepath, sequences);

        // Read back with AIRRReader
        auto info = get_file_info(filepath);
        REQUIRE(info.num_rows == 3);
    }

    SECTION("handles empty sequence list")
    {
        std::vector<SequenceData> empty_sequences;
        auto filepath = path("empty.tsv");

        write_tsv(filepath, empty_sequences);

        auto info = get_file_info(filepath);
        REQUIRE(info.num_rows == 0);
        REQUIRE(info.has_sequence_id);
        REQUIRE(info.has_sequence);
    }
}

//==============================================================================
// write_csv tests
//==============================================================================

TEST_CASE_METHOD(AIRRWriterFixture, "write_csv", "[airr][writer]")
{
    SECTION("writes valid AIRR CSV")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("test.csv");

        write_csv(filepath, sequences);

        // Check delimiter
        auto content = read_file_content(filepath);
        REQUIRE_THAT(content, ContainsSubstring(","));
        REQUIRE_THAT(content, !ContainsSubstring("\t"));
    }

    SECTION("can be read back")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("test_read.csv");

        write_csv(filepath, sequences);

        // Read back
        auto info = get_file_info(filepath, Delimiter::COMMA);
        REQUIRE(info.num_rows == 3);
        REQUIRE(info.delimiter == Delimiter::COMMA);
    }
}

//==============================================================================
// Round-trip tests
//==============================================================================

TEST_CASE_METHOD(AIRRWriterFixture, "round-trip", "[airr][writer][roundtrip]")
{
    SECTION("TSV round-trip preserves data")
    {
        auto original = create_test_sequences();
        auto filepath = path("roundtrip.tsv");

        // Write
        write_tsv(filepath, original);

        // Read back
        auto read_back = read_sequences(filepath);

        REQUIRE(read_back.size() == original.size());

        // Check first sequence
        REQUIRE(read_back[0].index == original[0].index);
        REQUIRE(read_back[0].sequence == original[0].sequence);

        // Check V alignment
        REQUIRE(read_back[0].alignments.count(V_gene) == 1);
        const auto& v_orig = original[0].alignments.at(V_gene)[0];
        const auto& v_read = read_back[0].alignments.at(V_gene)[0];
        REQUIRE(v_read.gene_name == v_orig.gene_name);
        REQUIRE(v_read.score == v_orig.score);
        // offset is converted: 0-based -> 1-based -> 0-based
        REQUIRE(v_read.offset == v_orig.offset);
        REQUIRE(v_read.align_length == v_orig.align_length);
    }

    SECTION("CSV round-trip preserves data")
    {
        auto original = create_test_sequences();
        auto filepath = path("roundtrip.csv");

        write_csv(filepath, original);
        auto read_back = read_sequences(filepath);

        REQUIRE(read_back.size() == original.size());
        REQUIRE(read_back[1].sequence == original[1].sequence);
    }

    SECTION("handles missing D gene correctly")
    {
        auto original = create_test_sequences();
        auto filepath = path("no_d.tsv");

        write_tsv(filepath, original);
        auto read_back = read_sequences(filepath);

        // Second sequence has no D gene
        REQUIRE(read_back[1].alignments.count(D_gene) == 0);
        // But has V and J
        REQUIRE(read_back[1].alignments.count(V_gene) == 1);
        REQUIRE(read_back[1].alignments.count(J_gene) == 1);
    }
}

//==============================================================================
// make_cigar tests
//==============================================================================

TEST_CASE("make_cigar", "[airr][writer][cigar]")
{
    SECTION("simple match only")
    {
        std::forward_list<int> insertions;
        std::forward_list<int> deletions;
        std::vector<int> mismatches;

        Alignment_data align("IGHV1", 0, 0, 0, 50, insertions, deletions, mismatches, 100.0);

        REQUIRE(make_cigar(align) == "50M");
    }

    SECTION("with insertions")
    {
        std::forward_list<int> insertions = {5, 10};
        std::forward_list<int> deletions;
        std::vector<int> mismatches;

        Alignment_data align("IGHV1", 0, 0, 0, 50, insertions, deletions, mismatches, 100.0);

        // With corrected logic: align_length is used directly (not reduced by insertions)
        REQUIRE(make_cigar(align) == "50M2I");
    }

    SECTION("with deletions")
    {
        std::forward_list<int> insertions;
        std::forward_list<int> deletions = {15, 20, 25};
        std::vector<int> mismatches;

        Alignment_data align("IGHV1", 0, 0, 0, 50, insertions, deletions, mismatches, 100.0);

        REQUIRE(make_cigar(align) == "50M3D");
    }

    SECTION("with insertions and deletions")
    {
        std::forward_list<int> insertions = {5};
        std::forward_list<int> deletions = {10, 15};
        std::vector<int> mismatches;

        Alignment_data align("IGHV1", 0, 0, 0, 50, insertions, deletions, mismatches, 100.0);

        REQUIRE(make_cigar(align) == "50M1I2D");
    }

    SECTION("empty alignment")
    {
        std::forward_list<int> insertions;
        std::forward_list<int> deletions;
        std::vector<int> mismatches;

        Alignment_data align("IGHV1", 0, 0, 0, 0, insertions, deletions, mismatches, 0.0);

        REQUIRE(make_cigar(align) == "");
    }
}

//==============================================================================
// get_airr_columns tests
//==============================================================================

TEST_CASE("get_airr_columns", "[airr][writer]")
{
    SECTION("with alignment details")
    {
        auto columns = get_airr_columns(true);

        REQUIRE_THAT(columns, SizeIs(17));  // 2 base + 5*3 alignment
        REQUIRE(columns[0] == "sequence_id");
        REQUIRE(columns[1] == "sequence");
        REQUIRE(columns[2] == "v_call");
    }

    SECTION("without alignment details")
    {
        auto columns = get_airr_columns(false);

        REQUIRE_THAT(columns, SizeIs(2));
        REQUIRE(columns[0] == "sequence_id");
        REQUIRE(columns[1] == "sequence");
    }
}

//==============================================================================
// write_legacy_tsv tests
//==============================================================================

TEST_CASE_METHOD(AIRRWriterFixture, "write_legacy_tsv", "[airr][writer][legacy]")
{
    // Create legacy format sequences
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>> legacy_seqs;

    std::forward_list<int> insertions;
    std::forward_list<int> deletions;
    std::vector<int> mismatches;

    std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;
    aligns[V_gene].emplace_back("IGHV1-2*01", 0, 0, 0, 50, insertions, deletions, mismatches, 150.0);

    legacy_seqs.emplace_back(1, "ATCGATCGATCG", aligns);

    auto filepath = path("legacy.tsv");
    write_legacy_tsv(filepath, legacy_seqs);

    // Verify
    auto info = get_file_info(filepath);
    REQUIRE(info.num_rows == 1);
}

//==============================================================================
// Error handling tests
//==============================================================================

TEST_CASE_METHOD(AIRRWriterFixture, "write error handling", "[airr][writer][error]")
{
    SECTION("throws on invalid path")
    {
        std::vector<SequenceData> sequences;

        REQUIRE_THROWS_AS(
            write_tsv("/nonexistent/directory/file.tsv", sequences),
            std::runtime_error
        );
    }
}
