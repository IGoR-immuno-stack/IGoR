/*
 * test_airr_rearrangement_reader.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  Unit tests for AIRR reader functions
 *
 *  Uses Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures
 *  - GENERATE for parametric tests
 *  - SECTION for grouping related tests
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "StreamingTestUtils.h"
#include <igor/Streaming/AIRRRearrangementReader.h>

#include <fstream>

using namespace igor;
using namespace igor::airr::rearrangement;
using namespace igor::test;
using namespace Catch::Matchers;

// Import shared utility functions from parent namespace
using igor::airr::detect_delimiter;
using igor::airr::delimiter_char;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for AIRR reader tests with automatic cleanup
 *
 * Creates test AIRR files in a temporary directory that is
 * automatically cleaned up when the fixture is destroyed.
 */
struct AIRRReaderFixture
{
    TestDirectory test_dir{"igor_airr_reader_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }

    /**
     * @brief Create a minimal AIRR TSV file
     */
    std::string create_minimal_tsv()
    {
        auto filepath = path("minimal.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsequence\n";
        file << "seq_001\tATCGATCGATCG\n";
        file << "seq_002\tGCTAGCTAGCTA\n";
        file << "seq_003\tTTAATTAATTAA\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an AIRR TSV file with V/D/J calls
     */
    std::string create_with_vdj_calls()
    {
        auto filepath = path("with_vdj.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsequence\tv_call\td_call\tj_call\tv_score\td_score\tj_score\n";
        file << "1\tATCGATCGATCGATCG\tIGHV1-2*01\tIGHD1-1*01\tIGHJ3*01\t150.5\t50.2\t120.3\n";
        file << "2\tGCTAGCTAGCTAGCTA\tIGHV3-9*01\t\tIGHJ4*02\t145.0\t\t115.8\n";
        file << "3\tTTAATTAATTAATTAA\tIGHV4-4*01\tIGHD2-2*01\tIGHJ5*01\t160.0\t55.5\t130.0\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an AIRR CSV file
     */
    std::string create_csv()
    {
        auto filepath = path("data.csv");
        std::ofstream file(filepath);
        file << "sequence_id,sequence,v_call,productive\n";
        file << "1,ATCGATCG,IGHV1-2*01,T\n";
        file << "2,GCTAGCTA,IGHV3-9*01,F\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an AIRR TSV with full alignment info
     */
    std::string create_with_alignment_positions()
    {
        auto filepath = path("with_positions.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsequence\tv_call\tv_sequence_start\tv_sequence_end\t"
             << "v_germline_start\tv_germline_end\tv_alignment_length\tv_score\n";
        file << "1\tATCGATCGATCGATCG\tIGHV1-2*01\t1\t50\t1\t50\t50\t150.5\n";
        file << "2\tGCTAGCTAGCTAGCTA\tIGHV3-9*01\t5\t45\t10\t50\t40\t140.0\n";
        file.close();
        return filepath;
    }

    /**
     * @brief Create an empty file (header only)
     */
    std::string create_empty_file()
    {
        auto filepath = path("empty.tsv");
        std::ofstream file(filepath);
        file << "sequence_id\tsequence\n";
        file.close();
        return filepath;
    }
};

//==============================================================================
// Delimiter detection tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "detect_delimiter", "[airr][reader]")
{
    SECTION("detects TSV delimiter")
    {
        auto filepath = create_minimal_tsv();
        REQUIRE(detect_delimiter(filepath) == Delimiter::TAB);
    }

    SECTION("detects CSV delimiter")
    {
        auto filepath = create_csv();
        REQUIRE(detect_delimiter(filepath) == Delimiter::COMMA);
    }

    SECTION("throws on non-existent file")
    {
        REQUIRE_THROWS_AS(
            detect_delimiter("/nonexistent/file.tsv"),
            std::runtime_error
        );
    }
}

TEST_CASE("delimiter_char", "[airr][reader][utility]")
{
    REQUIRE(delimiter_char(Delimiter::TAB) == '\t');
    REQUIRE(delimiter_char(Delimiter::COMMA) == ',');
}

//==============================================================================
// get_file_info tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "get_file_info", "[airr][reader]")
{
    SECTION("basic TSV file")
    {
        auto filepath = create_minimal_tsv();
        auto info = get_file_info(filepath);

        REQUIRE(info.filepath == filepath);
        REQUIRE(info.num_rows == 3);
        REQUIRE(info.delimiter == Delimiter::TAB);
        REQUIRE(info.has_sequence_id);
        REQUIRE(info.has_sequence);
        REQUIRE_FALSE(info.has_v_call);
        REQUIRE_THAT(info.column_names, SizeIs(2));
    }

    SECTION("file with V/D/J calls")
    {
        auto filepath = create_with_vdj_calls();
        auto info = get_file_info(filepath);

        REQUIRE(info.num_rows == 3);
        REQUIRE(info.has_v_call);
        REQUIRE(info.has_d_call);
        REQUIRE(info.has_j_call);
        REQUIRE_THAT(info.column_names, SizeIs(8));
    }

    SECTION("empty file (header only)")
    {
        auto filepath = create_empty_file();
        auto info = get_file_info(filepath);

        REQUIRE(info.num_rows == 0);
        REQUIRE(info.has_sequence_id);
        REQUIRE(info.has_sequence);
    }

    SECTION("CSV file auto-detection")
    {
        auto filepath = create_csv();
        auto info = get_file_info(filepath, Delimiter::AUTO);

        REQUIRE(info.delimiter == Delimiter::COMMA);
        REQUIRE(info.num_rows == 2);
    }
}

//==============================================================================
// validate_schema tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "validate_schema", "[airr][reader]")
{
    SECTION("valid AIRR file")
    {
        auto filepath = create_minimal_tsv();
        REQUIRE(validate_schema(filepath));
    }

    SECTION("file with V/D/J calls is valid")
    {
        auto filepath = create_with_vdj_calls();
        REQUIRE(validate_schema(filepath));
    }

    SECTION("file without sequence_id is invalid")
    {
        auto filepath = path("no_seq_id.tsv");
        std::ofstream file(filepath);
        file << "sequence\tv_call\n";
        file << "ATCGATCG\tIGHV1-2*01\n";
        file.close();

        REQUIRE_FALSE(validate_schema(filepath));
    }
}

//==============================================================================
// read_batch tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "read_batch", "[airr][reader]")
{
    SECTION("reads minimal TSV")
    {
        auto filepath = create_minimal_tsv();
        auto batch = read_batch(filepath);

        REQUIRE(batch.nb_rows() == 3);
        REQUIRE(batch.nb_columns() == 2);
    }

    SECTION("reads file with V/D/J calls")
    {
        auto filepath = create_with_vdj_calls();
        auto batch = read_batch(filepath);

        REQUIRE(batch.nb_rows() == 3);
        REQUIRE(batch.nb_columns() == 8);
    }

    SECTION("reads CSV file")
    {
        auto filepath = create_csv();
        auto batch = read_batch(filepath, Delimiter::AUTO);

        REQUIRE(batch.nb_rows() == 2);
        REQUIRE(batch.nb_columns() == 4);
    }

    SECTION("handles empty file")
    {
        auto filepath = create_empty_file();
        auto batch = read_batch(filepath);

        REQUIRE(batch.nb_rows() == 0);
        REQUIRE(batch.nb_columns() == 2);
    }
}

//==============================================================================
// read_sequences tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "read_sequences", "[airr][reader]")
{
    SECTION("reads minimal TSV")
    {
        auto filepath = create_minimal_tsv();
        auto sequences = read_sequences(filepath);

        REQUIRE(sequences.size() == 3);
        REQUIRE(sequences[0].sequence == "ATCGATCGATCG");
        REQUIRE(sequences[1].sequence == "GCTAGCTAGCTA");
        REQUIRE(sequences[2].sequence == "TTAATTAATTAA");
    }

    SECTION("extracts V/D/J alignments")
    {
        auto filepath = create_with_vdj_calls();
        auto sequences = read_sequences(filepath);

        REQUIRE(sequences.size() == 3);

        // First sequence has V, D, J alignments
        const auto& aligns0 = sequences[0].alignments;
        REQUIRE(aligns0.count(V_gene) == 1);
        REQUIRE(aligns0.count(D_gene) == 1);
        REQUIRE(aligns0.count(J_gene) == 1);

        REQUIRE(aligns0.at(V_gene)[0].gene_name == "IGHV1-2*01");
        REQUIRE(aligns0.at(V_gene)[0].score == 150.5);
        REQUIRE(aligns0.at(D_gene)[0].gene_name == "IGHD1-1*01");
        REQUIRE(aligns0.at(J_gene)[0].gene_name == "IGHJ3*01");

        // Second sequence has only V and J (D is empty)
        const auto& aligns1 = sequences[1].alignments;
        REQUIRE(aligns1.count(V_gene) == 1);
        REQUIRE(aligns1.count(D_gene) == 0); // No D call in file
        REQUIRE(aligns1.count(J_gene) == 1);
    }

    SECTION("extracts alignment positions")
    {
        auto filepath = create_with_alignment_positions();
        auto sequences = read_sequences(filepath);

        REQUIRE(sequences.size() == 2);

        const auto& v_align = sequences[0].alignments.at(V_gene)[0];
        REQUIRE(v_align.gene_name == "IGHV1-2*01");
        REQUIRE(v_align.offset == 0);  // 1-based converted to 0-based
        REQUIRE(v_align.align_length == 50);
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
// read_legacy tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "read_legacy", "[airr][reader]")
{
    auto filepath = create_with_vdj_calls();
    auto sequences = read_legacy(filepath);

    REQUIRE(sequences.size() == 3);

    // Check tuple structure
    auto& [id, seq, aligns] = sequences[0];
    REQUIRE(id == 1);
    REQUIRE(seq == "ATCGATCGATCGATCG");
    REQUIRE(aligns.count(V_gene) == 1);
}

//==============================================================================
// read_columns tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "read_columns", "[airr][reader]")
{
    SECTION("reads specific columns")
    {
        auto filepath = create_with_vdj_calls();
        auto batch = read_columns(filepath, {"sequence_id", "sequence"});

        REQUIRE(batch.nb_rows() == 3);
        REQUIRE(batch.nb_columns() == 2);
    }

    SECTION("throws on missing column")
    {
        auto filepath = create_minimal_tsv();
        REQUIRE_THROWS_AS(
            read_columns(filepath, {"nonexistent_column"}),
            std::runtime_error
        );
    }
}

//==============================================================================
// Error handling tests
//==============================================================================

TEST_CASE_METHOD(AIRRReaderFixture, "error handling", "[airr][reader][error]")
{
    SECTION("throws on non-existent file")
    {
        std::string bad_path = "/nonexistent/path/file.tsv";

        REQUIRE_THROWS_AS(read_batch(bad_path), std::runtime_error);
        REQUIRE_THROWS_AS(read_sequences(bad_path), std::runtime_error);
        REQUIRE_THROWS_AS(get_file_info(bad_path), std::runtime_error);
    }

    SECTION("throws on empty file content")
    {
        auto filepath = path("truly_empty.tsv");
        std::ofstream file(filepath);
        file.close();  // File with no content

        REQUIRE_THROWS_AS(read_batch(filepath), std::runtime_error);
    }
}
