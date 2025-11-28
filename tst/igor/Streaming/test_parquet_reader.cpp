/*
 * test_parquet_reader.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ParquetReader module
 */

// Include Arrow compatibility workarounds before any other headers
#include <igor/Streaming/ArrowCompatibility.h>

#include <catch2/catch_test_macros.hpp>
#include <igor/Streaming/ParquetReader.h>
#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <sparrow/record_batch.hpp>
#include <sparrow/builder.hpp>

#include <filesystem>
#include <fstream>

using namespace igor;

// Test fixture for Parquet reader tests
static const std::string TEST_OUTPUT_DIR = "/tmp/igor_parquet_reader_tests";

// Helper to ensure test directory exists
static void setup_test_directory()
{
    std::filesystem::create_directories(TEST_OUTPUT_DIR);
}

// Helper to create test sequences
static std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
create_test_sequences(size_t count)
{
    std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>> sequences;
    sequences.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        int id = static_cast<int>(i);
        std::string seq = "ACGTACGT" + std::to_string(i);
        std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

        sequences.emplace_back(id, seq, alignments);
    }

    return sequences;
}

TEST_CASE("read simple parquet file", "[parquet][reader]")
{
    setup_test_directory();

    // Create a test file using the writer
    auto sequences = create_test_sequences(10);
    std::string test_file = TEST_OUTPUT_DIR + "/read_simple.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    // Now read it back
    REQUIRE_NOTHROW([&]() {
        auto batch = ParquetReader::read_batch(test_file);
        REQUIRE(batch.nb_rows() == 10);
        REQUIRE(batch.nb_columns() >= 2); // At least id and sequence
    }());
}

TEST_CASE("read_batch returns correct data", "[parquet][reader]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(5);
    std::string test_file = TEST_OUTPUT_DIR + "/read_batch_data.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    auto batch = ParquetReader::read_batch(test_file);

    REQUIRE(batch.nb_rows() == 5);

    // Verify we can access the data
    for (size_t i = 0; i < batch.nb_rows(); ++i) {
        auto seq_data = row_to_sequence_data(batch, i);
        REQUIRE(seq_data.index == static_cast<int>(i));
        REQUIRE(seq_data.sequence == "ACGTACGT" + std::to_string(i));
    }
}

TEST_CASE("read_sequences converts to legacy format", "[parquet][reader]")
{
    setup_test_directory();

    auto original_sequences = create_test_sequences(10);
    std::string test_file = TEST_OUTPUT_DIR + "/read_sequences.parquet";

    ParquetWriter::write_sequences(test_file, original_sequences);

    auto read_sequences = ParquetReader::read_sequences(test_file);

    REQUIRE(read_sequences.size() == 10);

    // Verify data matches
    for (size_t i = 0; i < read_sequences.size(); ++i) {
        REQUIRE(std::get<0>(read_sequences[i]) == std::get<0>(original_sequences[i]));
        REQUIRE(std::get<1>(read_sequences[i]) == std::get<1>(original_sequences[i]));
    }
}

TEST_CASE("round-trip preserves all data", "[parquet][reader][writer]")
{
    setup_test_directory();

    auto original = create_test_sequences(50);
    std::string test_file = TEST_OUTPUT_DIR + "/roundtrip.parquet";

    // Write
    ParquetWriter::write_sequences(test_file, original);

    // Read
    auto read_back = ParquetReader::read_sequences(test_file);

    // Verify
    REQUIRE(read_back.size() == original.size());

    for (size_t i = 0; i < original.size(); ++i) {
        REQUIRE(std::get<0>(read_back[i]) == std::get<0>(original[i]));
        REQUIRE(std::get<1>(read_back[i]) == std::get<1>(original[i]));
    }
}

TEST_CASE("get_file_info returns correct metadata", "[parquet][reader][metadata]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(100);
    std::string test_file = TEST_OUTPUT_DIR + "/metadata_test.parquet";

    ParquetWriter::write_sequences(test_file, sequences, CompressionType::SNAPPY);

    auto info = ParquetReader::get_file_info(test_file);

    REQUIRE(info.num_rows == 100);
    REQUIRE(info.num_columns >= 2);
    REQUIRE(info.compression == "SNAPPY");
    REQUIRE(info.column_names.size() >= 2);

    // Check that we have the expected columns
    bool has_id = false;
    bool has_sequence = false;
    for (const auto &name : info.column_names) {
        if (name == "sequence_id") has_id = true;
        if (name == "sequence") has_sequence = true;
    }
    REQUIRE(has_id);
    REQUIRE(has_sequence);
}

TEST_CASE("read empty parquet file", "[parquet][reader][edge-case]")
{
    setup_test_directory();

    std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>> empty;
    std::string test_file = TEST_OUTPUT_DIR + "/empty.parquet";

    ParquetWriter::write_sequences(test_file, empty);

    auto batch = ParquetReader::read_batch(test_file);
    REQUIRE(batch.nb_rows() == 0);

    auto sequences = ParquetReader::read_sequences(test_file);
    REQUIRE(sequences.empty());
}

TEST_CASE("read non-existent file throws", "[parquet][reader][error]")
{
    std::string non_existent = TEST_OUTPUT_DIR + "/does_not_exist.parquet";

    REQUIRE_THROWS_AS(ParquetReader::read_batch(non_existent), std::runtime_error);
    REQUIRE_THROWS_AS(ParquetReader::read_sequences(non_existent), std::runtime_error);
    REQUIRE_THROWS_AS(ParquetReader::get_file_info(non_existent), std::runtime_error);
}

TEST_CASE("read_columns reads subset of columns", "[parquet][reader]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(10);
    std::string test_file = TEST_OUTPUT_DIR + "/columns_test.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    // Read only the sequence column
    auto batch = ParquetReader::read_columns(test_file, {"sequence"});

    REQUIRE(batch.nb_rows() == 10);
    REQUIRE(batch.nb_columns() == 1);

    auto column_names = batch.names();
    REQUIRE(column_names.size() == 1);
    REQUIRE(column_names[0] == "sequence");
}

TEST_CASE("read_columns with multiple columns", "[parquet][reader]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(10);
    std::string test_file = TEST_OUTPUT_DIR + "/multi_columns.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    auto batch = ParquetReader::read_columns(test_file, {"sequence_id", "sequence"});

    REQUIRE(batch.nb_rows() == 10);
    REQUIRE(batch.nb_columns() == 2);
}

TEST_CASE("read_columns with invalid column throws", "[parquet][reader][error]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(5);
    std::string test_file = TEST_OUTPUT_DIR + "/invalid_column.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    REQUIRE_THROWS_AS(
        ParquetReader::read_columns(test_file, {"non_existent_column"}),
        std::runtime_error
    );
}

TEST_CASE("read different compression types", "[parquet][reader][compression]")
{
    setup_test_directory();

    auto sequences = create_test_sequences(20);

    std::vector<std::pair<CompressionType, std::string>> compressions = {
        {CompressionType::NONE, "NONE"},
        {CompressionType::SNAPPY, "SNAPPY"},
        {CompressionType::GZIP, "GZIP"},
        {CompressionType::ZSTD, "ZSTD"},
        {CompressionType::LZ4, "LZ4"}
    };

    for (const auto &[type, name] : compressions) {
        std::string test_file = TEST_OUTPUT_DIR + "/compression_" + name + ".parquet";

        INFO("Testing compression: " << name);

        ParquetWriter::write_sequences(test_file, sequences, type);

        auto info = ParquetReader::get_file_info(test_file);
        REQUIRE(info.compression == name);

        auto read_back = ParquetReader::read_sequences(test_file);
        REQUIRE(read_back.size() == sequences.size());
    }
}

TEST_CASE("read large dataset efficiently", "[parquet][reader][performance]")
{
    setup_test_directory();

    const size_t num_sequences = 10000;
    auto sequences = create_test_sequences(num_sequences);
    std::string test_file = TEST_OUTPUT_DIR + "/large_dataset.parquet";

    ParquetWriter::write_sequences(test_file, sequences);

    auto start = std::chrono::high_resolution_clock::now();
    auto read_sequences = ParquetReader::read_sequences(test_file);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    REQUIRE(read_sequences.size() == num_sequences);

    // Should be able to read 10K sequences in under 1 second
    INFO("Read " << num_sequences << " sequences in " << duration_ms << " ms");
    REQUIRE(duration_ms < 1000);

    double sequences_per_second = (num_sequences * 1000.0) / std::max(duration_ms, 1LL);
    INFO("Read speed: " << sequences_per_second << " sequences/second");
    REQUIRE(sequences_per_second > 5000); // Should read at least 5K seqs/sec
}
TEST_CASE("read-write-read preserves data exactly", "[parquet][reader][writer][roundtrip]")
{
    setup_test_directory();

    // Create original test data
    auto original_sequences = create_test_sequences(100);
    std::string file1 = TEST_OUTPUT_DIR + "/read_write_read_1.parquet";
    std::string file2 = TEST_OUTPUT_DIR + "/read_write_read_2.parquet";
    
    // Write original data
    ParquetWriter::write_sequences(file1, original_sequences);
    
    // Read it back as a batch
    auto batch = ParquetReader::read_batch(file1);
    
    // Write the batch to a new file
    ParquetWriter::write_batch(file2, batch);
    
    // Read both files and compare
    auto sequences_from_file1 = ParquetReader::read_sequences(file1);
    auto sequences_from_file2 = ParquetReader::read_sequences(file2);
    
    // Verify same number of sequences
    REQUIRE(sequences_from_file1.size() == sequences_from_file2.size());
    REQUIRE(sequences_from_file1.size() == original_sequences.size());
    
    // Verify all data matches exactly
    for (size_t i = 0; i < sequences_from_file1.size(); ++i) {
        REQUIRE(std::get<0>(sequences_from_file1[i]) == std::get<0>(sequences_from_file2[i]));
        REQUIRE(std::get<0>(sequences_from_file1[i]) == std::get<0>(original_sequences[i]));
        
        REQUIRE(std::get<1>(sequences_from_file1[i]) == std::get<1>(sequences_from_file2[i]));
        REQUIRE(std::get<1>(sequences_from_file1[i]) == std::get<1>(original_sequences[i]));
    }
    
    // Verify file metadata is consistent
    auto info1 = ParquetReader::get_file_info(file1);
    auto info2 = ParquetReader::get_file_info(file2);
    
    REQUIRE(info1.num_rows == info2.num_rows);
    REQUIRE(info1.num_columns == info2.num_columns);
    REQUIRE(info1.column_names == info2.column_names);
}
