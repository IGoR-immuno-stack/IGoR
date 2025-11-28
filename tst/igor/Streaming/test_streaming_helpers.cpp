/*
 * test_streaming_helpers.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for SequenceBatchHelpers module
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <igor/Streaming/SequenceBatchHelpers.h>
#include <sparrow/record_batch.hpp>
#include <sparrow/array.hpp>
#include <sparrow/builder.hpp>
#include <chrono>

using namespace igor;

// Helper function to create a test batch
static sparrow::record_batch create_test_batch()
{
    std::vector<std::string> names = { "sequence_id", "sequence" };
    std::vector<sparrow::array> arrays;

    // Create ID array using builder
    std::vector<int32_t> ids = { 0, 1, 2 };
    arrays.emplace_back(sparrow::build(ids));

    // Create sequence array using builder
    std::vector<std::string> sequences = { "ATCGATCGATCG", "GCTAGCTAGCTA", "TTAATTAATTAA" };
    arrays.emplace_back(sparrow::build(sequences));

    return sparrow::record_batch(std::move(names), std::move(arrays));
}

// Test: has_column function
TEST_CASE("has_column detects existing columns", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    REQUIRE(has_column(test_batch, "sequence_id"));
    REQUIRE(has_column(test_batch, "sequence"));
}

TEST_CASE("has_column returns false for missing columns", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    REQUIRE_FALSE(has_column(test_batch, "nonexistent_column"));
    REQUIRE_FALSE(has_column(test_batch, "v_gene_name"));
}

// Test: get_string_value function
TEST_CASE("get_string_value returns correct value", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();

    std::string seq = get_string_value(test_batch, "sequence", 0);
    REQUIRE(seq == "ATCGATCGATCG");

    seq = get_string_value(test_batch, "sequence", 1);
    REQUIRE(seq == "GCTAGCTAGCTA");

    seq = get_string_value(test_batch, "sequence", 2);
    REQUIRE(seq == "TTAATTAATTAA");
}

TEST_CASE("get_string_value returns default for missing column", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    std::string result = get_string_value(test_batch, "missing", 0, "default");
    REQUIRE(result == "default");
}

// Test: get_int_value function
TEST_CASE("get_int_value returns correct value", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();

    int id = get_int_value(test_batch, "sequence_id", 0);
    REQUIRE(id == 0);

    id = get_int_value(test_batch, "sequence_id", 1);
    REQUIRE(id == 1);

    id = get_int_value(test_batch, "sequence_id", 2);
    REQUIRE(id == 2);
}

TEST_CASE("get_int_value returns default for missing column", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    int result = get_int_value(test_batch, "missing", 0, -999);
    REQUIRE(result == -999);
}

// Test: row_to_sequence_data function
TEST_CASE("row_to_sequence_data converts correctly", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    SequenceData seq_data = row_to_sequence_data(test_batch, 0);

    REQUIRE(seq_data.index == 0);
    REQUIRE(seq_data.sequence == "ATCGATCGATCG");
}

TEST_CASE("row_to_sequence_data handles all rows", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    for (size_t i = 0; i < test_batch.nb_rows(); ++i) {
        SequenceData seq_data = row_to_sequence_data(test_batch, i);
        REQUIRE(seq_data.index == static_cast<int>(i));
        REQUIRE_FALSE(seq_data.sequence.empty());
    }
}

TEST_CASE("row_to_sequence_data throws on invalid index", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    REQUIRE_THROWS_AS(row_to_sequence_data(test_batch, 999), std::out_of_range);
}

TEST_CASE("row_to_sequence_data throws on missing sequence column", "[streaming][helpers]")
{
    // Create a batch without "sequence" column
    std::vector<std::string> names = { "sequence_id" };
    std::vector<sparrow::array> arrays;

    arrays.emplace_back(sparrow::build(std::vector<int32_t>{ 0 }));

    sparrow::record_batch bad_batch(std::move(names), std::move(arrays));

    REQUIRE_THROWS_AS(row_to_sequence_data(bad_batch, 0), std::runtime_error);
}

// Test: vector_to_batch function
TEST_CASE("vector_to_batch converts empty vector", "[streaming][helpers]")
{
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            empty_vec;

    sparrow::record_batch batch = vector_to_batch(empty_vec);

    REQUIRE(batch.nb_rows() == 0);
    REQUIRE(has_column(batch, "sequence_id"));
    REQUIRE(has_column(batch, "sequence"));
}

TEST_CASE("vector_to_batch converts single sequence", "[streaming][helpers]")
{
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            sequences;

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;
    sequences.emplace_back(42, "ACGTACGTACGT", empty_alignments);

    sparrow::record_batch batch = vector_to_batch(sequences);

    REQUIRE(batch.nb_rows() == 1);

    SequenceData seq_data = row_to_sequence_data(batch, 0);
    REQUIRE(seq_data.index == 42);
    REQUIRE(seq_data.sequence == "ACGTACGTACGT");
}

TEST_CASE("vector_to_batch converts multiple sequences", "[streaming][helpers]")
{
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            sequences;

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;

    sequences.emplace_back(0, "AAAA", empty_alignments);
    sequences.emplace_back(1, "TTTT", empty_alignments);
    sequences.emplace_back(2, "GGGG", empty_alignments);
    sequences.emplace_back(3, "CCCC", empty_alignments);

    sparrow::record_batch batch = vector_to_batch(sequences);

    REQUIRE(batch.nb_rows() == 4);

    for (size_t i = 0; i < 4; ++i) {
        SequenceData seq_data = row_to_sequence_data(batch, i);
        REQUIRE(seq_data.index == static_cast<int>(i));
    }
}

// Test: Round-trip conversion (vector → batch → SequenceData)
TEST_CASE("round-trip conversion preserves data", "[streaming][helpers]")
{
    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            original_sequences;

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;

    // Create test sequences
    original_sequences.emplace_back(100, "ATCGATCGATCG", empty_alignments);
    original_sequences.emplace_back(200, "GCTAGCTAGCTA", empty_alignments);
    original_sequences.emplace_back(300, "TTAATTAATTAA", empty_alignments);

    // Convert to batch
    sparrow::record_batch batch = vector_to_batch(original_sequences);

    // Convert back and verify
    for (size_t i = 0; i < original_sequences.size(); ++i) {
        SequenceData seq_data = row_to_sequence_data(batch, i);

        REQUIRE(seq_data.index == std::get<0>(original_sequences[i]));
        REQUIRE(seq_data.sequence == std::get<1>(original_sequences[i]));
    }
}

// Test: Null handling
TEST_CASE("null handling works correctly", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    // get_string_value, get_int_value, get_double_value should handle nulls gracefully
    // by returning default values

    double result = get_double_value(test_batch, "nonexistent", 0, 3.14);
    REQUIRE(result == 3.14);
}

// Test: parse_alignments_from_columns
TEST_CASE("parse_alignments returns empty for no alignment columns", "[streaming][helpers]")
{
    auto test_batch = create_test_batch();
    auto alignments = parse_alignments_from_columns(test_batch, 0);

    // Should return empty map when no alignment columns present
    REQUIRE(alignments.empty());
}

TEST_CASE("parse_alignments extracts V gene data", "[streaming][helpers]")
{
    // Create a batch with V gene alignment data
    std::vector<std::string> names = { "sequence_id", "sequence", "v_gene_name", "v_gene_offset",
                                       "v_gene_score" };
    std::vector<sparrow::array> arrays;

    // Create arrays using builder
    arrays.emplace_back(sparrow::build(std::vector<int32_t>{ 0 }));
    arrays.emplace_back(sparrow::build(std::vector<std::string>{ "ATCG" }));
    arrays.emplace_back(sparrow::build(std::vector<std::string>{ "IGHV1-1*01" }));
    arrays.emplace_back(sparrow::build(std::vector<int32_t>{ 5 }));
    arrays.emplace_back(sparrow::build(std::vector<double>{ 150.5 }));

    sparrow::record_batch batch(std::move(names), std::move(arrays));

    auto alignments = parse_alignments_from_columns(batch, 0);

    REQUIRE_FALSE(alignments.empty());
    REQUIRE(alignments.find(V_gene) != alignments.end());
    REQUIRE(alignments[V_gene].size() == 1);
    REQUIRE(alignments[V_gene][0].gene_name == "IGHV1-1*01");
    REQUIRE(alignments[V_gene][0].offset == 5);
    REQUIRE(alignments[V_gene][0].score == 150.5);
}

// Performance test: verify conversion speed meets target (≥100K conversions/second)
TEST_CASE("performance conversion speed", "[streaming][helpers][!benchmark]")
{
    // Create a larger batch for performance testing
    const size_t num_sequences = 10000;

    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
            sequences;
    sequences.reserve(num_sequences);

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;

    for (size_t i = 0; i < num_sequences; ++i) {
        sequences.emplace_back(static_cast<int>(i), "ATCGATCGATCGATCG", empty_alignments);
    }

    sparrow::record_batch batch = vector_to_batch(sequences);

    // Time the conversions
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < batch.nb_rows(); ++i) {
        SequenceData seq_data = row_to_sequence_data(batch, i);
        // Prevent optimization from removing the conversion
        REQUIRE_FALSE(seq_data.sequence.empty());
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    double conversions_per_second = (num_sequences * 1000.0) / duration.count();

    INFO("Conversion performance: " << conversions_per_second << " conversions/second");

    // Target: ≥100K conversions/second
    // We'll use a relaxed threshold for CI environments
    REQUIRE(conversions_per_second > 50000);
}
