/*
 * test_streaming_helpers.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for SequenceBatchHelpers module
 *
 *  Refactored to use Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures
 *  - GENERATE for parametric tests
 *  - SECTION for grouping related tests
 *  - Matchers for cleaner assertions
 *  - Proper BENCHMARK macro (hidden by default)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "StreamingTestUtils.h"
#include <igor/Streaming/SequenceBatchHelpers.h>

using namespace igor;
using namespace igor::test;
using namespace Catch::Matchers;

//==============================================================================
// Fixtures
//==============================================================================

struct HelpersFixture
{
    sparrow::record_batch batch = create_test_batch();
};

//==============================================================================
// has_column tests
//==============================================================================

TEST_CASE_METHOD(HelpersFixture, "has_column", "[streaming][helpers]")
{
    SECTION("detects existing columns")
    {
        REQUIRE(has_column(batch, "sequence_id"));
        REQUIRE(has_column(batch, "sequence"));
    }

    SECTION("returns false for missing columns")
    {
        auto col = GENERATE("nonexistent_column", "v_gene_name", "");
        CAPTURE(col);
        REQUIRE_FALSE(has_column(batch, col));
    }
}

//==============================================================================
// get_*_value tests
//==============================================================================

TEST_CASE_METHOD(HelpersFixture, "get_string_value", "[streaming][helpers]")
{
    SECTION("returns correct values")
    {
        auto [row, expected] = GENERATE(table<size_t, std::string>({
            {0, "ATCGATCGATCG"},
            {1, "GCTAGCTAGCTA"},
            {2, "TTAATTAATTAA"}
        }));

        CAPTURE(row);
        REQUIRE(get_string_value(batch, "sequence", row) == expected);
    }

    SECTION("returns default for missing column")
    {
        REQUIRE(get_string_value(batch, "missing", 0, "default") == "default");
    }
}

TEST_CASE_METHOD(HelpersFixture, "get_int_value", "[streaming][helpers]")
{
    SECTION("returns correct values")
    {
        auto row = GENERATE(0, 1, 2);
        CAPTURE(row);
        REQUIRE(get_int_value(batch, "sequence_id", row) == row);
    }

    SECTION("returns default for missing column")
    {
        REQUIRE(get_int_value(batch, "missing", 0, -999) == -999);
    }
}

TEST_CASE_METHOD(HelpersFixture, "get_double_value", "[streaming][helpers]")
{
    SECTION("returns default for missing column")
    {
        REQUIRE(get_double_value(batch, "nonexistent", 0, 3.14) == 3.14);
    }
}

//==============================================================================
// row_to_sequence_data tests
//==============================================================================

TEST_CASE_METHOD(HelpersFixture, "row_to_sequence_data", "[streaming][helpers]")
{
    SECTION("converts first row correctly")
    {
        auto seq_data = row_to_sequence_data(batch, 0);
        REQUIRE(seq_data.index == 0);
        REQUIRE(seq_data.sequence == "ATCGATCGATCG");
    }

    SECTION("handles all rows")
    {
        auto row = GENERATE(size_t{0}, size_t{1}, size_t{2});
        CAPTURE(row);

        auto seq_data = row_to_sequence_data(batch, row);
        REQUIRE(seq_data.index == static_cast<int>(row));
        REQUIRE_THAT(seq_data.sequence, !IsEmpty());
    }

    SECTION("throws on invalid index")
    {
        REQUIRE_THROWS_AS(row_to_sequence_data(batch, 999), std::out_of_range);
    }

    SECTION("throws on missing sequence column")
    {
        // Create batch without "sequence" column
        std::vector<std::string> names = {"sequence_id"};
        std::vector<sparrow::array> arrays;
        arrays.emplace_back(sparrow::build(std::vector<int32_t>{0}));
        sparrow::record_batch bad_batch(std::move(names), std::move(arrays));

        REQUIRE_THROWS_AS(row_to_sequence_data(bad_batch, 0), std::runtime_error);
    }
}

//==============================================================================
// vector_to_batch tests
//==============================================================================

TEST_CASE("vector_to_batch", "[streaming][helpers]")
{
    SECTION("converts empty vector")
    {
        std::vector<SequenceTuple> empty_vec;
        auto batch = vector_to_batch(empty_vec);

        REQUIRE(batch.nb_rows() == 0);
        REQUIRE(has_column(batch, "sequence_id"));
        REQUIRE(has_column(batch, "sequence"));
    }

    SECTION("converts various sizes")
    {
        auto count = GENERATE(1, 4, 50);
        CAPTURE(count);

        auto sequences = create_test_sequences(static_cast<size_t>(count));
        auto batch = vector_to_batch(sequences);

        REQUIRE(batch.nb_rows() == static_cast<size_t>(count));

        for (size_t i = 0; i < batch.nb_rows(); ++i) {
            auto seq_data = row_to_sequence_data(batch, i);
            REQUIRE(seq_data.index == static_cast<int>(i));
        }
    }
}

//==============================================================================
// Round-trip tests
//==============================================================================

TEST_CASE("round-trip conversion preserves data", "[streaming][helpers]")
{
    auto count = GENERATE(1, 10, 100);
    CAPTURE(count);

    auto original = create_test_sequences(static_cast<size_t>(count));
    auto batch = vector_to_batch(original);

    REQUIRE(batch.nb_rows() == original.size());

    for (size_t i = 0; i < original.size(); ++i) {
        auto recovered = row_to_sequence_data(batch, i);
        REQUIRE(recovered.index == std::get<0>(original[i]));
        REQUIRE(recovered.sequence == std::get<1>(original[i]));
    }
}

TEST_CASE("round-trip with populated alignments", "[streaming][helpers][regression]")
{
    // Regression test: ensures alignment list fields survive round-trip
    auto [id, seq, alignments] = create_sequence_with_v_alignment(1, "ACGT");

    std::vector<SequenceTuple> sequences = {{id, seq, alignments}};
    auto batch = vector_to_batch(sequences);
    auto recovered = row_to_sequence_data(batch, 0);

    REQUIRE(recovered.alignments.count(V_gene) == 1);
    REQUIRE(recovered.alignments.at(V_gene).size() == 1);

    const auto& orig = alignments.at(V_gene)[0];
    const auto& res = recovered.alignments.at(V_gene)[0];

    REQUIRE(res.gene_name == orig.gene_name);
    REQUIRE(res.offset == orig.offset);
    REQUIRE(res.score == orig.score);

    // Compare forward_list fields
    std::vector<int> orig_ins(orig.insertions.begin(), orig.insertions.end());
    std::vector<int> res_ins(res.insertions.begin(), res.insertions.end());
    REQUIRE_THAT(res_ins, RangeEquals(orig_ins));

    std::vector<int> orig_del(orig.deletions.begin(), orig.deletions.end());
    std::vector<int> res_del(res.deletions.begin(), res.deletions.end());
    REQUIRE_THAT(res_del, RangeEquals(orig_del));

    REQUIRE(res.mismatches == orig.mismatches);
}

//==============================================================================
// parse_alignments tests
//==============================================================================

TEST_CASE_METHOD(HelpersFixture, "parse_alignments_from_columns", "[streaming][helpers]")
{
    SECTION("returns empty for batch without alignment columns")
    {
        auto alignments = parse_alignments_from_columns(batch, 0);
        REQUIRE_THAT(alignments, IsEmpty());
    }

    SECTION("extracts V gene data")
    {
        // Create batch with V gene alignment data
        std::vector<std::string> names = {
            "sequence_id", "sequence",
            "v_gene_name", "v_gene_offset", "v_gene_score"
        };
        std::vector<sparrow::array> arrays;

        arrays.emplace_back(sparrow::build(std::vector<int32_t>{0}));
        arrays.emplace_back(sparrow::build(std::vector<std::string>{"ATCG"}));
        arrays.emplace_back(sparrow::build(std::vector<std::string>{"IGHV1-1*01"}));
        arrays.emplace_back(sparrow::build(std::vector<int32_t>{5}));
        arrays.emplace_back(sparrow::build(std::vector<double>{150.5}));

        sparrow::record_batch v_batch(std::move(names), std::move(arrays));
        auto alignments = parse_alignments_from_columns(v_batch, 0);

        REQUIRE_THAT(alignments, !IsEmpty());
        REQUIRE(alignments.count(V_gene) == 1);
        REQUIRE(alignments.at(V_gene).size() == 1);

        const auto& v = alignments.at(V_gene)[0];
        REQUIRE(v.gene_name == "IGHV1-1*01");
        REQUIRE(v.offset == 5);
        REQUIRE(v.score == 150.5);
    }
}

//==============================================================================
// Benchmarks (hidden by default, run with [.benchmark])
//==============================================================================

TEST_CASE("Conversion benchmarks", "[.benchmark]")
{
    constexpr size_t NUM_SEQUENCES = 10000;
    auto sequences = create_test_sequences(NUM_SEQUENCES);
    auto batch = vector_to_batch(sequences);

    BENCHMARK("row_to_sequence_data - single row")
    {
        return row_to_sequence_data(batch, 0);
    };

    BENCHMARK("row_to_sequence_data - all rows")
    {
        for (size_t i = 0; i < batch.nb_rows(); ++i) {
            auto seq = row_to_sequence_data(batch, i);
            (void)seq;  // Prevent optimization
        }
    };

    BENCHMARK("vector_to_batch - 100 sequences")
    {
        auto small = create_test_sequences(100);
        return vector_to_batch(small);
    };

    BENCHMARK("vector_to_batch - 1000 sequences")
    {
        auto medium = create_test_sequences(1000);
        return vector_to_batch(medium);
    };
}
