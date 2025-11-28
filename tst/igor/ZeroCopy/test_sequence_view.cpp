/*
 * test_sequence_view.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ZeroCopy SequenceView and AlignmentView
 */

#include <igor/ZeroCopy/SequenceView.h>

#include <catch2/catch_test_macros.hpp>
#include <sparrow/array.hpp>
#include <sparrow/builder.hpp>
#include <sparrow/record_batch.hpp>

using namespace igor;

// Helper to create a test batch with AIRR columns
static sparrow::record_batch create_airr_batch() {
    std::vector<std::string> names = {
        "sequence_id",    "sequence", "productive", "v_call",  "v_score", "v_sequence_start",
        "v_sequence_end", "v_cigar",  "d_call",     "d_score", "j_call",  "j_score"};
    std::vector<sparrow::array> arrays;

    // sequence_id
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"seq1", "seq2"}));
    // sequence
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"ACGT", "TGCA"}));
    // productive
    arrays.emplace_back(sparrow::build(std::vector<bool>{true, false}));

    // V gene
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"IGHV1-69", "IGHV3-23"}));
    arrays.emplace_back(sparrow::build(std::vector<double>{100.5, 90.0}));
    arrays.emplace_back(sparrow::build(std::vector<int32_t>{1, 2}));
    arrays.emplace_back(sparrow::build(std::vector<int32_t>{10, 12}));
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"10M", "11M"}));

    // D gene (partial data)
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"IGHD1-1", ""}));
    arrays.emplace_back(sparrow::build(std::vector<double>{50.0, 0.0}));

    // J gene
    arrays.emplace_back(sparrow::build(std::vector<std::string>{"IGHJ4", "IGHJ6"}));
    arrays.emplace_back(sparrow::build(std::vector<double>{80.0, 85.5}));

    return sparrow::record_batch(std::move(names), std::move(arrays));
}

TEST_CASE("SequenceView reads basic fields", "[zerocopy][sequence_view]") {
    auto batch = create_airr_batch();
    SequenceView view(batch, 0);

    REQUIRE(view.is_valid());
    REQUIRE(view.get_sequence_id() == "seq1");
    REQUIRE(view.get_sequence() == "ACGT");
    REQUIRE(view.productive() == true);

    SequenceView view2(batch, 1);
    REQUIRE(view2.get_sequence_id() == "seq2");
    REQUIRE(view2.get_sequence() == "TGCA");
    REQUIRE(view2.productive() == false);
}

TEST_CASE("AlignmentView reads V gene fields", "[zerocopy][alignment_view]") {
    auto batch = create_airr_batch();
    SequenceView view(batch, 0);

    auto v_align = view.get_alignment(V_gene);
    REQUIRE(v_align.is_valid());
    REQUIRE(v_align.get_gene_call() == "IGHV1-69");
    REQUIRE(v_align.get_score() == 100.5);
    REQUIRE(v_align.get_sequence_start() == 1);
    REQUIRE(v_align.get_sequence_end() == 10);
    REQUIRE(v_align.get_cigar() == "10M");

    // Legacy mapping
    REQUIRE(v_align.get_offset() == 0);         // 1-based start 1 -> 0-based offset 0
    REQUIRE(v_align.get_align_length() == 10);  // 10 - 1 + 1
}

TEST_CASE("AlignmentView handles missing/empty data", "[zerocopy][alignment_view]") {
    auto batch = create_airr_batch();
    SequenceView view(batch, 1);  // Row 2

    auto d_align = view.get_alignment(D_gene);
    REQUIRE_FALSE(d_align.is_valid());  // d_call is empty string
    REQUIRE(d_align.get_gene_call() == "");
    REQUIRE(d_align.get_score() == 0.0);
}
