/*
 * test_airr_alignment.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  Unit tests for AIRR Alignment schema reader and writer
 *
 *  Uses Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures
 *  - SECTION for grouping related tests
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "StreamingTestUtils.h"
#include <igor/Streaming/AIRRAlignmentReader.h>
#include <igor/Streaming/AIRRAlignmentWriter.h>

#include <fstream>
#include <sstream>

using namespace igor;
using namespace igor::airr::alignment;
using namespace igor::test;
using namespace Catch::Matchers;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for AIRR Alignment tests with automatic cleanup
 */
struct AIRRAlignmentFixture
{
    TestDirectory test_dir{"igor_airr_alignment_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }

    /**
     * @brief Create a minimal AIRR Alignment TSV file
     */
    std::string create_minimal_alignment_tsv()
    {
        auto filepath = path("minimal_alignment.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsegment\tcall\tscore\tsequence_start\tsequence_end\tgermline_start\tcigar\trank\n";
        file << "1\tV\tIGHV1-2*01\t150.5\t1\t50\t1\t50M\t1\n";
        file << "1\tD\tIGHD1-1*01\t45.2\t51\t70\t1\t20M\t1\n";
        file << "1\tJ\tIGHJ3*01\t120.0\t71\t100\t1\t30M\t1\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an AIRR Alignment file with multiple alignments per gene
     */
    std::string create_with_ranked_alignments()
    {
        auto filepath = path("ranked_alignments.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsegment\tcall\tscore\tsequence_start\tsequence_end\tgermline_start\tcigar\trank\n";
        // Sequence 1: primary and secondary V alignments
        file << "1\tV\tIGHV1-2*01\t150.5\t1\t50\t1\t50M\t1\n";
        file << "1\tV\tIGHV1-3*01\t145.0\t1\t48\t1\t48M\t2\n";
        file << "1\tJ\tIGHJ3*01\t120.0\t71\t100\t1\t30M\t1\n";
        // Sequence 2: V and J only
        file << "2\tV\tIGHV3-9*01\t140.0\t1\t45\t1\t45M\t1\n";
        file << "2\tJ\tIGHJ4*02\t115.0\t46\t70\t1\t25M\t1\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an AIRR Alignment CSV file
     */
    std::string create_csv()
    {
        auto filepath = path("alignment.csv");
        std::ofstream file(filepath);
        file << "sequence_id,segment,call,score,rank\n";
        file << "1,V,IGHV1-2*01,150.5,1\n";
        file << "1,J,IGHJ3*01,120.0,1\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an empty alignment file (header only)
     */
    std::string create_empty_file()
    {
        auto filepath = path("empty.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsegment\tcall\tscore\trank\n";
        file.close();
        return filepath;
    }

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

        // Sequence 2: V with secondary alignment and J
        {
            std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns;

            std::forward_list<int> insertions;
            std::forward_list<int> deletions;
            std::vector<int> mismatches;

            // Primary V
            aligns[V_gene].emplace_back(
                "IGHV3-9*01", 0, 0, 0, 45, insertions, deletions, mismatches, 140.0);
            // Secondary V
            aligns[V_gene].emplace_back(
                "IGHV3-11*01", 0, 0, 0, 43, insertions, deletions, mismatches, 135.0);
            // J
            aligns[J_gene].emplace_back(
                "IGHJ4*02", 45, 0, 0, 25, insertions, deletions, mismatches, 115.0);

            sequences.emplace_back(2, "GCTAGCTAGCTAGCTAGCTAGCTA", aligns);
        }

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
// get_file_info tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment::get_file_info", "[airr][alignment][reader]")
{
    SECTION("basic alignment file")
    {
        auto filepath = create_minimal_alignment_tsv();
        auto info = get_file_info(filepath);

        REQUIRE(info.filepath == filepath);
        REQUIRE(info.num_rows == 3);
        REQUIRE(info.delimiter == Delimiter::TAB);
        REQUIRE(info.has_sequence_id);
        REQUIRE(info.has_segment);
        REQUIRE(info.has_call);
        REQUIRE(info.has_cigar);
        REQUIRE(info.has_rank);
    }

    SECTION("CSV file auto-detection")
    {
        auto filepath = create_csv();
        auto info = get_file_info(filepath, Delimiter::AUTO);

        REQUIRE(info.delimiter == Delimiter::COMMA);
        REQUIRE(info.num_rows == 2);
    }

    SECTION("empty file (header only)")
    {
        auto filepath = create_empty_file();
        auto info = get_file_info(filepath);

        REQUIRE(info.num_rows == 0);
        REQUIRE(info.has_sequence_id);
        REQUIRE(info.has_segment);
        REQUIRE(info.has_call);
    }
}

//==============================================================================
// validate_schema tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment::validate_schema", "[airr][alignment][reader]")
{
    SECTION("valid alignment file")
    {
        auto filepath = create_minimal_alignment_tsv();
        REQUIRE(validate_schema(filepath));
    }

    SECTION("file without segment column is invalid")
    {
        auto filepath = path("no_segment.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tcall\tscore\n";
        file << "1\tIGHV1-2*01\t150.5\n";
        file.close();

        REQUIRE_FALSE(validate_schema(filepath));
    }

    SECTION("file without call column is invalid")
    {
        auto filepath = path("no_call.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsegment\tscore\n";
        file << "1\tV\t150.5\n";
        file.close();

        REQUIRE_FALSE(validate_schema(filepath));
    }
}

//==============================================================================
// read_batch tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment::read_batch", "[airr][alignment][reader]")
{
    SECTION("reads minimal alignment TSV")
    {
        auto filepath = create_minimal_alignment_tsv();
        auto batch = read_batch(filepath);

        REQUIRE(batch.nb_rows() == 3);
        REQUIRE(batch.nb_columns() == 9);
    }

    SECTION("reads CSV file")
    {
        auto filepath = create_csv();
        auto batch = read_batch(filepath, Delimiter::AUTO);

        REQUIRE(batch.nb_rows() == 2);
        REQUIRE(batch.nb_columns() == 5);
    }

    SECTION("handles empty file")
    {
        auto filepath = create_empty_file();
        auto batch = read_batch(filepath);

        REQUIRE(batch.nb_rows() == 0);
        REQUIRE(batch.nb_columns() == 5);
    }
}

//==============================================================================
// read_sequences tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment::read_sequences", "[airr][alignment][reader]")
{
    SECTION("reads and groups by sequence_id")
    {
        auto filepath = create_minimal_alignment_tsv();
        auto sequences = read_sequences(filepath);

        REQUIRE(sequences.size() == 1);
        REQUIRE(sequences[0].alignments.count(V_gene) == 1);
        REQUIRE(sequences[0].alignments.count(D_gene) == 1);
        REQUIRE(sequences[0].alignments.count(J_gene) == 1);
    }

    SECTION("handles multiple alignments per gene")
    {
        auto filepath = create_with_ranked_alignments();
        auto sequences = read_sequences(filepath);

        REQUIRE(sequences.size() == 2);

        // Find sequence 1 (has 2 V alignments)
        auto it = std::find_if(sequences.begin(), sequences.end(),
            [](const SequenceData& s) { return s.index == 1; });
        REQUIRE(it != sequences.end());
        REQUIRE(it->alignments.at(V_gene).size() == 2);
    }

    SECTION("extracts alignment scores")
    {
        auto filepath = create_minimal_alignment_tsv();
        auto sequences = read_sequences(filepath);

        const auto& v_align = sequences[0].alignments.at(V_gene)[0];
        REQUIRE(v_align.gene_name == "IGHV1-2*01");
        REQUIRE(v_align.score == 150.5);
    }

    SECTION("handles empty file")
    {
        auto filepath = create_empty_file();
        auto sequences = read_sequences(filepath);
        REQUIRE_THAT(sequences, IsEmpty());
    }
}

//==============================================================================
// parse_cigar tests
//==============================================================================

TEST_CASE("alignment::parse_cigar", "[airr][alignment][reader][cigar]")
{
    std::forward_list<int> insertions;
    std::forward_list<int> deletions;
    size_t align_length = 0;

    SECTION("simple match only")
    {
        REQUIRE(parse_cigar("50M", insertions, deletions, align_length));
        REQUIRE(align_length == 50);
        REQUIRE(insertions.empty());
        REQUIRE(deletions.empty());
    }

    SECTION("match with insertions")
    {
        REQUIRE(parse_cigar("50M2I", insertions, deletions, align_length));
        REQUIRE(align_length == 50);  // I does not add to align_length
        REQUIRE(std::distance(insertions.begin(), insertions.end()) == 2);
    }

    SECTION("match with deletions")
    {
        REQUIRE(parse_cigar("50M3D", insertions, deletions, align_length));
        REQUIRE(align_length == 53);  // D adds to align_length
        REQUIRE(std::distance(deletions.begin(), deletions.end()) == 3);
    }

    SECTION("complex CIGAR")
    {
        REQUIRE(parse_cigar("30M2I10M3D", insertions, deletions, align_length));
        REQUIRE(align_length == 43);  // 30 + 10 + 3
        REQUIRE(std::distance(insertions.begin(), insertions.end()) == 2);
        REQUIRE(std::distance(deletions.begin(), deletions.end()) == 3);
    }

    SECTION("empty CIGAR returns false")
    {
        REQUIRE_FALSE(parse_cigar("", insertions, deletions, align_length));
    }

    SECTION("invalid CIGAR returns false")
    {
        REQUIRE_FALSE(parse_cigar("X", insertions, deletions, align_length));
    }
}

//==============================================================================
// write_sequences tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment::write_sequences", "[airr][alignment][writer]")
{
    SECTION("writes valid alignment TSV")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("output.tsv");

        write_sequences(filepath, sequences);

        auto content = read_file_content(filepath);
        REQUIRE_FALSE(content.empty());

        // Check header
        REQUIRE_THAT(content, ContainsSubstring("sequence_id"));
        REQUIRE_THAT(content, ContainsSubstring("segment"));
        REQUIRE_THAT(content, ContainsSubstring("call"));
        REQUIRE_THAT(content, ContainsSubstring("score"));
        REQUIRE_THAT(content, ContainsSubstring("rank"));
    }

    SECTION("writes one row per alignment")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("rows.tsv");

        write_sequences(filepath, sequences);

        auto info = get_file_info(filepath);
        // Sequence 1: 3 alignments (V, D, J)
        // Sequence 2: 3 alignments (V primary, V secondary, J)
        REQUIRE(info.num_rows == 6);
    }

    SECTION("handles empty sequence list")
    {
        std::vector<SequenceData> empty_sequences;
        auto filepath = path("empty_output.tsv");

        write_sequences(filepath, empty_sequences);

        auto info = get_file_info(filepath);
        REQUIRE(info.num_rows == 0);
    }

    SECTION("writes CSV format")
    {
        auto sequences = create_test_sequences();
        auto filepath = path("output.csv");

        write_sequences(filepath, sequences, Delimiter::COMMA);

        auto content = read_file_content(filepath);
        REQUIRE_THAT(content, ContainsSubstring(","));
        REQUIRE_THAT(content, !ContainsSubstring("\t"));
    }
}

//==============================================================================
// make_cigar tests
//==============================================================================

TEST_CASE("alignment::make_cigar", "[airr][alignment][writer][cigar]")
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
// Round-trip tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment round-trip", "[airr][alignment][roundtrip]")
{
    SECTION("TSV round-trip preserves data")
    {
        auto original = create_test_sequences();
        auto filepath = path("roundtrip.tsv");

        // Write
        write_sequences(filepath, original);

        // Read back
        auto read_back = read_sequences(filepath);

        REQUIRE(read_back.size() == original.size());

        // Find sequence 1
        auto it = std::find_if(read_back.begin(), read_back.end(),
            [](const SequenceData& s) { return s.index == 1; });
        REQUIRE(it != read_back.end());

        // Check V alignment
        REQUIRE(it->alignments.count(V_gene) == 1);
        const auto& v_read = it->alignments.at(V_gene)[0];
        REQUIRE(v_read.gene_name == "IGHV1-2*01");
        REQUIRE(v_read.score == 150.5);
    }

    SECTION("CSV round-trip preserves data")
    {
        auto original = create_test_sequences();
        auto filepath = path("roundtrip.csv");

        write_sequences(filepath, original, Delimiter::COMMA);
        auto read_back = read_sequences(filepath);

        REQUIRE(read_back.size() == original.size());
    }

    SECTION("preserves multiple alignments per gene")
    {
        auto original = create_test_sequences();
        auto filepath = path("multi_align.tsv");

        write_sequences(filepath, original);
        auto read_back = read_sequences(filepath);

        // Find sequence 2 (has 2 V alignments)
        auto it = std::find_if(read_back.begin(), read_back.end(),
            [](const SequenceData& s) { return s.index == 2; });
        REQUIRE(it != read_back.end());
        REQUIRE(it->alignments.at(V_gene).size() == 2);
    }
}

//==============================================================================
// Error handling tests
//==============================================================================

TEST_CASE_METHOD(AIRRAlignmentFixture, "alignment error handling", "[airr][alignment][error]")
{
    SECTION("throws on non-existent file for read")
    {
        std::string bad_path = "/nonexistent/path/file.tsv";

        REQUIRE_THROWS_AS(read_batch(bad_path), std::runtime_error);
        REQUIRE_THROWS_AS(read_sequences(bad_path), std::runtime_error);
        REQUIRE_THROWS_AS(get_file_info(bad_path), std::runtime_error);
    }

    SECTION("throws on invalid path for write")
    {
        std::vector<SequenceData> sequences;

        REQUIRE_THROWS_AS(
            write_sequences("/nonexistent/directory/file.tsv", sequences),
            std::runtime_error
        );
    }
}
