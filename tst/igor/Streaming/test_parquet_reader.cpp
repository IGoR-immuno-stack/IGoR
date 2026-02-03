/*
 * test_parquet_reader.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ParquetReader module
 *
 *  Refactored to use Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures with RAII cleanup
 *  - GENERATE for parametric tests
 *  - SECTION for grouping related tests
 *  - Matchers for cleaner assertions
 *  - Proper BENCHMARK macro (hidden by default)
 */

// Include Arrow compatibility workarounds before any other headers
#include <igor/Streaming/ArrowCompatibility.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "StreamingTestUtils.h"
#include <igor/Streaming/ParquetReader.h>
#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

using namespace igor;
using namespace igor::test;
using namespace Catch::Matchers;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for ParquetReader tests
 *
 * Creates a test directory and provides helper methods for creating
 * test files. Automatically cleans up on destruction.
 */
struct ParquetReaderFixture
{
    TestDirectory test_dir{"igor_parquet_reader_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }

    /// Create a test file with the given number of sequences
    std::string create_test_file(const std::string& name, size_t count)
    {
        auto sequences = create_test_sequences(count);
        auto filepath = path(name);
        ParquetWriter::write_sequences(filepath, sequences);
        return filepath;
    }
};

//==============================================================================
// Basic read tests
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "read simple parquet file", "[parquet][reader]")
{
    auto file = create_test_file("simple.parquet", 10);

    auto batch = ParquetReader::read_batch(file);

    REQUIRE(batch.nb_rows() == 10);
    REQUIRE(batch.nb_columns() >= 2);  // At least id and sequence
}

TEST_CASE_METHOD(ParquetReaderFixture, "read_batch returns correct data", "[parquet][reader]")
{
    auto file = create_test_file("read_data.parquet", 5);
    auto batch = ParquetReader::read_batch(file);

    REQUIRE(batch.nb_rows() == 5);

    for (size_t i = 0; i < batch.nb_rows(); ++i) {
        auto seq_data = row_to_sequence_data(batch, i);
        REQUIRE(seq_data.index == static_cast<int>(i));
        REQUIRE_THAT(seq_data.sequence, !IsEmpty());
    }
}

TEST_CASE_METHOD(ParquetReaderFixture, "read_sequences converts to legacy format", "[parquet][reader]")
{
    auto original = create_test_sequences(10);
    ParquetWriter::write_sequences(path("legacy.parquet"), original);

    auto read_sequences = ParquetReader::read_sequences(path("legacy.parquet"));

    REQUIRE(read_sequences.size() == 10);

    for (size_t i = 0; i < read_sequences.size(); ++i) {
        REQUIRE(std::get<0>(read_sequences[i]) == std::get<0>(original[i]));
        REQUIRE(std::get<1>(read_sequences[i]) == std::get<1>(original[i]));
    }
}

//==============================================================================
// Round-trip tests
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "round-trip preserves all data", "[parquet][reader][writer]")
{
    auto count = GENERATE(10, 50, 100);
    CAPTURE(count);

    auto original = create_test_sequences(static_cast<size_t>(count));
    ParquetWriter::write_sequences(path("roundtrip.parquet"), original);

    auto read_back = ParquetReader::read_sequences(path("roundtrip.parquet"));

    REQUIRE(read_back.size() == original.size());

    for (size_t i = 0; i < original.size(); ++i) {
        REQUIRE(std::get<0>(read_back[i]) == std::get<0>(original[i]));
        REQUIRE(std::get<1>(read_back[i]) == std::get<1>(original[i]));
    }
}

TEST_CASE_METHOD(ParquetReaderFixture, "read-write-read preserves data exactly", "[parquet][reader][writer][roundtrip]")
{
    auto original = create_test_sequences(100);

    // Write original data
    ParquetWriter::write_sequences(path("file1.parquet"), original);

    // Read as batch and write again
    auto batch = ParquetReader::read_batch(path("file1.parquet"));
    ParquetWriter::write_batch(path("file2.parquet"), batch);

    // Read both and compare
    auto from_file1 = ParquetReader::read_sequences(path("file1.parquet"));
    auto from_file2 = ParquetReader::read_sequences(path("file2.parquet"));

    REQUIRE(from_file1.size() == from_file2.size());
    REQUIRE(from_file1.size() == original.size());

    for (size_t i = 0; i < from_file1.size(); ++i) {
        REQUIRE(std::get<0>(from_file1[i]) == std::get<0>(from_file2[i]));
        REQUIRE(std::get<1>(from_file1[i]) == std::get<1>(from_file2[i]));
    }

    // Verify file metadata consistency
    auto info1 = ParquetReader::get_file_info(path("file1.parquet"));
    auto info2 = ParquetReader::get_file_info(path("file2.parquet"));

    REQUIRE(info1.num_rows == info2.num_rows);
    REQUIRE(info1.num_columns == info2.num_columns);
    REQUIRE(info1.column_names == info2.column_names);
}

//==============================================================================
// Metadata tests
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "get_file_info returns correct metadata", "[parquet][reader][metadata]")
{
    auto sequences = create_test_sequences(100);
    ParquetWriter::write_sequences(path("metadata.parquet"), sequences, CompressionType::SNAPPY);

    auto info = ParquetReader::get_file_info(path("metadata.parquet"));

    REQUIRE(info.num_rows == 100);
    REQUIRE(info.num_columns >= 2);
    REQUIRE(info.compression == "SNAPPY");
    REQUIRE(info.column_names.size() >= 2);

    // Check expected columns exist
    REQUIRE_THAT(info.column_names, Contains(std::string("sequence_id")));
    REQUIRE_THAT(info.column_names, Contains(std::string("sequence")));
}

//==============================================================================
// Column selection tests
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "read_columns", "[parquet][reader]")
{
    create_test_file("columns.parquet", 10);

    SECTION("reads single column")
    {
        auto batch = ParquetReader::read_columns(path("columns.parquet"), {"sequence"});

        REQUIRE(batch.nb_rows() == 10);
        REQUIRE(batch.nb_columns() == 1);

        auto column_names = batch.names();
        REQUIRE(column_names.size() == 1);
        REQUIRE(column_names[0] == "sequence");
    }

    SECTION("reads multiple columns")
    {
        auto batch = ParquetReader::read_columns(path("columns.parquet"), {"sequence_id", "sequence"});

        REQUIRE(batch.nb_rows() == 10);
        REQUIRE(batch.nb_columns() == 2);
    }

    SECTION("throws on invalid column")
    {
        REQUIRE_THROWS_AS(
            ParquetReader::read_columns(path("columns.parquet"), {"non_existent_column"}),
            std::runtime_error
        );
    }
}

//==============================================================================
// Compression reading tests
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "read different compression types", "[parquet][reader][compression]")
{
    auto compression = GENERATE(
        CompressionType::NONE,
        CompressionType::SNAPPY,
        CompressionType::GZIP,
        CompressionType::ZSTD,
        CompressionType::LZ4
    );

    auto name = ParquetWriter::compression_name(compression);
    CAPTURE(name);

    auto sequences = create_test_sequences(20);
    ParquetWriter::write_sequences(path("compression_" + name + ".parquet"), sequences, compression);

    auto info = ParquetReader::get_file_info(path("compression_" + name + ".parquet"));
    REQUIRE(info.compression == name);

    auto read_back = ParquetReader::read_sequences(path("compression_" + name + ".parquet"));
    REQUIRE(read_back.size() == sequences.size());
}

//==============================================================================
// Edge cases
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "edge cases", "[parquet][reader][edge]")
{
    SECTION("read empty file")
    {
        std::vector<SequenceTuple> empty;
        ParquetWriter::write_sequences(path("empty.parquet"), empty);

        auto batch = ParquetReader::read_batch(path("empty.parquet"));
        REQUIRE(batch.nb_rows() == 0);

        auto sequences = ParquetReader::read_sequences(path("empty.parquet"));
        REQUIRE_THAT(sequences, IsEmpty());
    }

    SECTION("non-existent file throws")
    {
        std::string non_existent = path("does_not_exist.parquet");

        REQUIRE_THROWS_AS(ParquetReader::read_batch(non_existent), std::runtime_error);
        REQUIRE_THROWS_AS(ParquetReader::read_sequences(non_existent), std::runtime_error);
        REQUIRE_THROWS_AS(ParquetReader::get_file_info(non_existent), std::runtime_error);
    }
}

//==============================================================================
// Real data integration tests (hidden by default)
//==============================================================================

TEST_CASE_METHOD(ParquetReaderFixture, "Murugan dataset round-trip", "[parquet][reader][.integration]")
{
    // Load real test data (300 sequences with 73K+ alignments)
    auto real_sequences = load_murugan_dataset();
    REQUIRE(real_sequences.size() == 300);

    INFO("Loaded " << real_sequences.size() << " sequences from Murugan dataset");

    // Write to Parquet
    ParquetWriter::write_sequences(path("murugan.parquet"), real_sequences);
    REQUIRE(std::filesystem::exists(path("murugan.parquet")));

    // Check file info
    auto file_info = ParquetReader::get_file_info(path("murugan.parquet"));
    REQUIRE(file_info.num_rows == 300);

    constexpr int EXPECTED_COLUMNS = 29;  // 2 basic + 27 alignment fields
    REQUIRE(file_info.num_columns == EXPECTED_COLUMNS);

    // Read back
    auto read_sequences = ParquetReader::read_sequences(path("murugan.parquet"));
    REQUIRE(read_sequences.size() == 300);

    // Verify a sample of sequences with alignments
    for (size_t i = 0; i < std::min(size_t{10}, read_sequences.size()); ++i) {
        const auto& original = real_sequences[i];
        const auto& recovered = read_sequences[i];

        // Check basic fields
        REQUIRE(std::get<0>(recovered) == std::get<0>(original));
        REQUIRE(std::get<1>(recovered) == std::get<1>(original));

        // Check alignments for each gene class
        const auto& orig_aligns = std::get<2>(original);
        const auto& recov_aligns = std::get<2>(recovered);

        for (auto gene_class : {V_gene, D_gene, J_gene}) {
            if (orig_aligns.count(gene_class) > 0 && !orig_aligns.at(gene_class).empty()) {
                CAPTURE(i, static_cast<int>(gene_class));
                REQUIRE(recov_aligns.count(gene_class) > 0);

                const auto& orig_align = orig_aligns.at(gene_class)[0];
                const auto& recov_align = recov_aligns.at(gene_class)[0];

                CHECK(recov_align.gene_name == orig_align.gene_name);
                CHECK(recov_align.score == orig_align.score);
                CHECK(recov_align.offset == orig_align.offset);
                CHECK(recov_align.five_p_offset == orig_align.five_p_offset);
                CHECK(recov_align.three_p_offset == orig_align.three_p_offset);
                CHECK(recov_align.align_length == orig_align.align_length);

                // Check list fields
                std::vector<int> orig_insertions(orig_align.insertions.begin(), orig_align.insertions.end());
                std::vector<int> recov_insertions(recov_align.insertions.begin(), recov_align.insertions.end());
                CHECK(recov_insertions == orig_insertions);

                std::vector<int> orig_deletions(orig_align.deletions.begin(), orig_align.deletions.end());
                std::vector<int> recov_deletions(recov_align.deletions.begin(), recov_align.deletions.end());
                CHECK(recov_deletions == orig_deletions);

                CHECK(recov_align.mismatches == orig_align.mismatches);
            }
        }
    }
}

//==============================================================================
// Benchmarks (hidden by default)
//==============================================================================

TEST_CASE("Read benchmarks", "[.benchmark]")
{
    TestDirectory dir{"igor_parquet_read_benchmarks"};

    // Create test files of various sizes
    auto small = create_test_sequences(100);
    auto medium = create_test_sequences(1000);
    auto large = create_test_sequences(10000);

    ParquetWriter::write_sequences(dir.path("small.parquet"), small);
    ParquetWriter::write_sequences(dir.path("medium.parquet"), medium);
    ParquetWriter::write_sequences(dir.path("large.parquet"), large);

    BENCHMARK("read 100 sequences")
    {
        return ParquetReader::read_sequences(dir.path("small.parquet"));
    };

    BENCHMARK("read 1000 sequences")
    {
        return ParquetReader::read_sequences(dir.path("medium.parquet"));
    };

    BENCHMARK("read 10000 sequences")
    {
        return ParquetReader::read_sequences(dir.path("large.parquet"));
    };

    BENCHMARK("read_batch 10000 sequences")
    {
        return ParquetReader::read_batch(dir.path("large.parquet"));
    };
}
