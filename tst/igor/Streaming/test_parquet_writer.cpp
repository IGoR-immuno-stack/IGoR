/*
 * test_parquet_writer.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ParquetWriter module
 */

// Include Arrow compatibility workarounds before any other headers
#include <igor/Streaming/ArrowCompatibility.h>

#include <catch2/catch_test_macros.hpp>
#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

#include <sparrow/record_batch.hpp>
#include <sparrow/builder.hpp>

#include <filesystem>
#include <fstream>

using namespace igor;

// Test fixture for Parquet writer tests
static const std::string TEST_OUTPUT_DIR = "/tmp/igor_parquet_tests";

// Helper to ensure test directory exists
static void setup_test_directory()
{
    std::filesystem::create_directories(TEST_OUTPUT_DIR);
}

// Helper to clean up test files
static void cleanup_test_file(const std::string &filename)
{
    std::filesystem::path filepath = std::filesystem::path(TEST_OUTPUT_DIR) / filename;
    if (std::filesystem::exists(filepath)) {
        std::filesystem::remove(filepath);
    }
}

// Helper to get full test path
static std::string test_path(const std::string &filename)
{
    return (std::filesystem::path(TEST_OUTPUT_DIR) / filename).string();
}

// Helper to create test sequences
static std::vector<std::tuple<int, std::string,
                              std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
create_test_sequences(size_t count)
{
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            sequences;
    sequences.reserve(count);

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;

    for (size_t i = 0; i < count; ++i) {
        std::string seq = "ATCG";
        for (size_t j = 0; j < i % 10; ++j) {
            seq += "ACGT";
        }
        sequences.emplace_back(static_cast<int>(i), seq, empty_alignments);
    }

    return sequences;
}

// Test: Write simple Parquet file
TEST_CASE("write simple parquet file", "[parquet][writer]")
{
    setup_test_directory();
    std::string output_file = test_path("simple_test.parquet");
    cleanup_test_file("simple_test.parquet");

    // Create test sequences
    auto sequences = create_test_sequences(10);

    // Write to Parquet
    REQUIRE_NOTHROW(ParquetWriter::write_sequences(output_file, sequences));

    // Verify file exists
    REQUIRE(std::filesystem::exists(output_file));

    // Verify file size is reasonable (should be small)
    auto file_size = std::filesystem::file_size(output_file);
    REQUIRE(file_size > 0);
    REQUIRE(file_size < 10000); // Should be small for 10 sequences

    cleanup_test_file("simple_test.parquet");
}

// Test: Write with different compression types
TEST_CASE("write with compression types", "[parquet][writer][compression]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(100);

    std::vector<CompressionType> compressions = {
        CompressionType::NONE, CompressionType::SNAPPY, CompressionType::GZIP,
        CompressionType::ZSTD, CompressionType::LZ4
    };

    for (auto compression : compressions) {
        std::string filename = "test_" + ParquetWriter::compression_name(compression) + ".parquet";
        std::string output_file = test_path(filename);
        cleanup_test_file(filename);

        INFO("Testing compression: " << ParquetWriter::compression_name(compression));

        REQUIRE_NOTHROW(ParquetWriter::write_sequences(output_file, sequences, compression));
        REQUIRE(std::filesystem::exists(output_file));

        cleanup_test_file(filename);
    }
}

// Test: Compression reduces file size
TEST_CASE("compression reduces file size", "[parquet][writer][compression]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(1000);

    std::string uncompressed_file = test_path("uncompressed.parquet");
    std::string snappy_file = test_path("snappy.parquet");

    cleanup_test_file("uncompressed.parquet");
    cleanup_test_file("snappy.parquet");

    // Write without compression
    ParquetWriter::write_sequences(uncompressed_file, sequences, CompressionType::NONE);

    // Write with Snappy compression
    ParquetWriter::write_sequences(snappy_file, sequences, CompressionType::SNAPPY);

    auto uncompressed_size = std::filesystem::file_size(uncompressed_file);
    auto snappy_size = std::filesystem::file_size(snappy_file);

    INFO("Uncompressed: " << uncompressed_size << " bytes");
    INFO("Snappy: " << snappy_size << " bytes");
    INFO("Compression ratio: " << (double)uncompressed_size / snappy_size << "x");

    // Snappy should provide some compression
    REQUIRE(snappy_size < uncompressed_size);

    cleanup_test_file("uncompressed.parquet");
    cleanup_test_file("snappy.parquet");
}

// Test: Write from record_batch
TEST_CASE("write from record_batch", "[parquet][writer]")
{
    setup_test_directory();
    std::string output_file = test_path("batch_test.parquet");
    cleanup_test_file("batch_test.parquet");

    // Create batch using builder
    std::vector<int32_t> ids = {0, 1, 2, 3, 4};
    std::vector<std::string> sequences = {"ATCG", "GCTA", "TTAA", "GGCC", "ATAT"};

    std::vector<std::string> names = {"sequence_id", "sequence"};
    std::vector<sparrow::array> arrays;

    arrays.emplace_back(sparrow::build(ids));
    arrays.emplace_back(sparrow::build(sequences));

    sparrow::record_batch batch(std::move(names), std::move(arrays));

    // Write batch
    REQUIRE_NOTHROW(ParquetWriter::write_batch(output_file, batch));

    // Verify file exists
    REQUIRE(std::filesystem::exists(output_file));

    cleanup_test_file("batch_test.parquet");
}

// Test: Verify Parquet file is readable by Arrow
TEST_CASE("parquet file is readable by arrow", "[parquet][writer][validation]")
{
    setup_test_directory();
    std::string output_file = test_path("arrow_readable.parquet");
    cleanup_test_file("arrow_readable.parquet");

    // Create and write test data
    auto sequences = create_test_sequences(50);
    ParquetWriter::write_sequences(output_file, sequences);

    // Try to read with Arrow Parquet reader
    auto input_result = arrow::io::ReadableFile::Open(output_file);
    REQUIRE(input_result.ok());

    auto input_stream = input_result.ValueOrDie();

    std::unique_ptr<parquet::arrow::FileReader> reader;
    auto reader_status = parquet::arrow::OpenFile(input_stream, arrow::default_memory_pool(), &reader);
    REQUIRE(reader_status.ok());

    // Read metadata
    auto metadata = reader->parquet_reader()->metadata();
    REQUIRE(metadata->num_rows() == 50);
    REQUIRE(metadata->num_columns() == 2); // sequence_id, sequence

    // Read table
    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    REQUIRE(read_status.ok());
    REQUIRE(table->num_rows() == 50);

    cleanup_test_file("arrow_readable.parquet");
}

// Test: Empty sequences
TEST_CASE("write empty sequences", "[parquet][writer][edge-case]")
{
    setup_test_directory();
    std::string output_file = test_path("empty.parquet");
    cleanup_test_file("empty.parquet");

    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            empty_sequences;

    REQUIRE_NOTHROW(ParquetWriter::write_sequences(output_file, empty_sequences));
    REQUIRE(std::filesystem::exists(output_file));

    // Verify it's readable
    auto input_result = arrow::io::ReadableFile::Open(output_file);
    REQUIRE(input_result.ok());

    std::unique_ptr<parquet::arrow::FileReader> reader;
    auto reader_status = parquet::arrow::OpenFile(input_result.ValueOrDie(),
                                                  arrow::default_memory_pool(), &reader);
    REQUIRE(reader_status.ok());

    auto metadata = reader->parquet_reader()->metadata();
    REQUIRE(metadata->num_rows() == 0);

    cleanup_test_file("empty.parquet");
}

// Test: Large dataset (performance test)
TEST_CASE("write large dataset", "[parquet][writer][performance][!benchmark]")
{
    setup_test_directory();
    std::string output_file = test_path("large.parquet");
    cleanup_test_file("large.parquet");

    const size_t num_sequences = 10000;
    auto sequences = create_test_sequences(num_sequences);

    // Measure write time
    auto start = std::chrono::high_resolution_clock::now();

    ParquetWriter::write_sequences(output_file, sequences, CompressionType::SNAPPY);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    double sequences_per_second = (num_sequences * 1000.0) / duration.count();

    INFO("Wrote " << num_sequences << " sequences in " << duration.count() << " ms");
    INFO("Write speed: " << sequences_per_second << " sequences/second");

    // Target: ≥100K sequences/second
    // We'll use a relaxed threshold for CI
    REQUIRE(sequences_per_second > 5000);

    auto file_size = std::filesystem::file_size(output_file);
    INFO("File size: " << file_size << " bytes");

    cleanup_test_file("large.parquet");
}

// Test: Compression name function
TEST_CASE("compression name function", "[parquet][writer][utility]")
{
    REQUIRE(ParquetWriter::compression_name(CompressionType::NONE) == "NONE");
    REQUIRE(ParquetWriter::compression_name(CompressionType::SNAPPY) == "SNAPPY");
    REQUIRE(ParquetWriter::compression_name(CompressionType::GZIP) == "GZIP");
    REQUIRE(ParquetWriter::compression_name(CompressionType::ZSTD) == "ZSTD");
    REQUIRE(ParquetWriter::compression_name(CompressionType::LZ4) == "LZ4");
}

// Test: Invalid output path
TEST_CASE("write to invalid path fails gracefully", "[parquet][writer][error]")
{
    std::string invalid_path = "/nonexistent_directory/invalid/test.parquet";
    auto sequences = create_test_sequences(10);

    REQUIRE_THROWS_AS(ParquetWriter::write_sequences(invalid_path, sequences),
                     std::runtime_error);
}
