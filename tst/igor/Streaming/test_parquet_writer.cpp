/*
 * test_parquet_writer.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ParquetWriter module
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
#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

using namespace igor;
using namespace igor::test;
using namespace Catch::Matchers;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for ParquetWriter tests with automatic directory cleanup
 */
struct ParquetWriterFixture
{
    TestDirectory test_dir{"igor_parquet_writer_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }
};

//==============================================================================
// Basic write tests
//==============================================================================

TEST_CASE_METHOD(ParquetWriterFixture, "write simple parquet file", "[parquet][writer]")
{
    auto sequences = create_test_sequences(10);
    auto output = path("simple.parquet");

    REQUIRE_NOTHROW(ParquetWriter::write_sequences(output, sequences));
    REQUIRE(std::filesystem::exists(output));

    auto file_size = std::filesystem::file_size(output);
    REQUIRE(file_size > 0);
    REQUIRE(file_size < 10000);  // Small file for 10 sequences
}

TEST_CASE_METHOD(ParquetWriterFixture, "write from record_batch", "[parquet][writer]")
{
    // Create batch using builder
    std::vector<int32_t> ids = {0, 1, 2, 3, 4};
    std::vector<std::string> sequences = {"ATCG", "GCTA", "TTAA", "GGCC", "ATAT"};

    std::vector<std::string> names = {"sequence_id", "sequence"};
    std::vector<sparrow::array> arrays;

    arrays.emplace_back(sparrow::build(ids));
    arrays.emplace_back(sparrow::build(sequences));

    sparrow::record_batch batch(std::move(names), std::move(arrays));

    // Write batch
    REQUIRE_NOTHROW(ParquetWriter::write_batch(path("batch.parquet"), batch));
    REQUIRE(std::filesystem::exists(path("batch.parquet")));
}

//==============================================================================
// Compression tests (using GENERATE)
//==============================================================================

TEST_CASE_METHOD(ParquetWriterFixture, "write with compression", "[parquet][writer][compression]")
{
    auto compression = GENERATE(
        CompressionType::NONE,
        CompressionType::SNAPPY,
        CompressionType::GZIP,
        CompressionType::ZSTD,
        CompressionType::LZ4
    );

    CAPTURE(ParquetWriter::compression_name(compression));

    auto sequences = create_test_sequences(100);
    auto output = path("test_" + ParquetWriter::compression_name(compression) + ".parquet");

    REQUIRE_NOTHROW(ParquetWriter::write_sequences(output, sequences, compression));
    REQUIRE(std::filesystem::exists(output));
}

TEST_CASE_METHOD(ParquetWriterFixture, "compression reduces file size", "[parquet][writer][compression]")
{
    auto sequences = create_test_sequences(1000);

    ParquetWriter::write_sequences(path("uncompressed.parquet"), sequences, CompressionType::NONE);
    ParquetWriter::write_sequences(path("snappy.parquet"), sequences, CompressionType::SNAPPY);

    auto uncompressed_size = std::filesystem::file_size(path("uncompressed.parquet"));
    auto snappy_size = std::filesystem::file_size(path("snappy.parquet"));

    CAPTURE(uncompressed_size, snappy_size);
    INFO("Compression ratio: " << static_cast<double>(uncompressed_size) / snappy_size << "x");

    REQUIRE(snappy_size < uncompressed_size);
}

//==============================================================================
// Edge cases
//==============================================================================

TEST_CASE_METHOD(ParquetWriterFixture, "edge cases", "[parquet][writer][edge]")
{
    SECTION("empty sequences")
    {
        std::vector<SequenceTuple> empty;
        REQUIRE_NOTHROW(ParquetWriter::write_sequences(path("empty.parquet"), empty));
        REQUIRE(std::filesystem::exists(path("empty.parquet")));

        // Verify it's readable
        auto input = arrow::io::ReadableFile::Open(path("empty.parquet"));
        REQUIRE(input.ok());

        std::unique_ptr<parquet::arrow::FileReader> reader;
        REQUIRE(parquet::arrow::OpenFile(*input, arrow::default_memory_pool(), &reader).ok());

        auto metadata = reader->parquet_reader()->metadata();
        REQUIRE(metadata->num_rows() == 0);
    }

    SECTION("invalid path throws")
    {
        auto sequences = create_test_sequences(10);
        REQUIRE_THROWS_AS(
            ParquetWriter::write_sequences("/nonexistent_directory/invalid/test.parquet", sequences),
            std::runtime_error
        );
    }
}

//==============================================================================
// Validation tests
//==============================================================================

TEST_CASE_METHOD(ParquetWriterFixture, "parquet file is readable by Arrow", "[parquet][writer][validation]")
{
    auto sequences = create_test_sequences(50);
    ParquetWriter::write_sequences(path("arrow_readable.parquet"), sequences);

    // Read with Arrow Parquet reader
    auto input = arrow::io::ReadableFile::Open(path("arrow_readable.parquet"));
    REQUIRE(input.ok());

    std::unique_ptr<parquet::arrow::FileReader> reader;
    REQUIRE(parquet::arrow::OpenFile(*input, arrow::default_memory_pool(), &reader).ok());

    // Verify metadata
    auto metadata = reader->parquet_reader()->metadata();
    REQUIRE(metadata->num_rows() == 50);

    // 2 basic columns + 9 fields × 3 gene classes = 29 columns
    constexpr int EXPECTED_COLUMNS = 29;
    REQUIRE(metadata->num_columns() == EXPECTED_COLUMNS);

    // Read full table
    std::shared_ptr<arrow::Table> table;
    REQUIRE(reader->ReadTable(&table).ok());
    REQUIRE(table->num_rows() == 50);
}

//==============================================================================
// Utility function tests
//==============================================================================

TEST_CASE("compression_name", "[parquet][writer][utility]")
{
    auto [type, expected] = GENERATE(table<CompressionType, std::string>({
        {CompressionType::NONE, "NONE"},
        {CompressionType::SNAPPY, "SNAPPY"},
        {CompressionType::GZIP, "GZIP"},
        {CompressionType::ZSTD, "ZSTD"},
        {CompressionType::LZ4, "LZ4"}
    }));

    REQUIRE(ParquetWriter::compression_name(type) == expected);
}

//==============================================================================
// Real data integration tests (hidden by default)
//==============================================================================

TEST_CASE_METHOD(ParquetWriterFixture, "write Murugan dataset", "[parquet][writer][.integration]")
{
    auto real_sequences = load_murugan_dataset();
    REQUIRE(real_sequences.size() == 300);

    INFO("Loaded " << real_sequences.size() << " sequences from Murugan dataset");

    REQUIRE_NOTHROW(ParquetWriter::write_sequences(path("murugan.parquet"), real_sequences));
    REQUIRE(std::filesystem::exists(path("murugan.parquet")));

    auto file_size = std::filesystem::file_size(path("murugan.parquet"));
    CAPTURE(file_size);
    REQUIRE(file_size > 10000);  // Should have substantial data

    // Verify with Arrow
    auto input = arrow::io::ReadableFile::Open(path("murugan.parquet"));
    REQUIRE(input.ok());

    std::unique_ptr<parquet::arrow::FileReader> reader;
    REQUIRE(parquet::arrow::OpenFile(*input, arrow::default_memory_pool(), &reader).ok());

    std::shared_ptr<arrow::Table> table;
    REQUIRE(reader->ReadTable(&table).ok());

    constexpr int EXPECTED_COLUMNS = 29;  // 2 basic + 27 alignment fields
    REQUIRE(table->num_rows() == 300);
    REQUIRE(table->num_columns() == EXPECTED_COLUMNS);

    // Verify alignment column presence
    auto schema = table->schema();
    bool has_v_insertions = false;
    bool has_d_gene_name = false;
    bool has_j_score = false;

    for (int i = 0; i < schema->num_fields(); ++i) {
        std::string name = schema->field(i)->name();
        if (name == "v_gene_insertions") has_v_insertions = true;
        if (name == "d_gene_name") has_d_gene_name = true;
        if (name == "j_gene_score") has_j_score = true;
    }

    REQUIRE(has_v_insertions);
    REQUIRE(has_d_gene_name);
    REQUIRE(has_j_score);
}

//==============================================================================
// Benchmarks (hidden by default)
//==============================================================================

TEST_CASE("Write benchmarks", "[.benchmark]")
{
    TestDirectory dir{"igor_parquet_write_benchmarks"};

    BENCHMARK("write 100 sequences")
    {
        auto sequences = create_test_sequences(100);
        ParquetWriter::write_sequences(dir.path("bench_100.parquet"), sequences, CompressionType::SNAPPY);
    };

    BENCHMARK("write 1000 sequences")
    {
        auto sequences = create_test_sequences(1000);
        ParquetWriter::write_sequences(dir.path("bench_1000.parquet"), sequences, CompressionType::SNAPPY);
    };

    BENCHMARK("write 10000 sequences")
    {
        auto sequences = create_test_sequences(10000);
        ParquetWriter::write_sequences(dir.path("bench_10000.parquet"), sequences, CompressionType::SNAPPY);
    };
}
