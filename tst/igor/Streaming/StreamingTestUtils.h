/*
 * StreamingTestUtils.h
 *
 *  Created on: Feb 3, 2026
 *      Author: IGoR Development Team
 *
 *  Shared test utilities for Streaming module tests
 */

#pragma once

#include <igor/Streaming/SequenceBatchHelpers.h>
#include <igor/Core/Utils.h>
#include <igor/Core/Aligner.h>

#include <exception>
#include <sparrow/record_batch.hpp>
#include <sparrow/array.hpp>
#include <sparrow/builder.hpp>

#include <filesystem>
#include <forward_list>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace igor::test {

//==============================================================================
// Type Aliases (eliminates verbose tuple repetition)
//==============================================================================

/// Legacy sequence format type alias for readability
using SequenceTuple = std::tuple<
    int,
    std::string,
    std::unordered_map<Gene_class, std::vector<Alignment_data>>
>;

//==============================================================================
// RAII Test Directory Manager (replaces manual setup/cleanup)
//==============================================================================

/**
 * @brief RAII wrapper for test directory management
 *
 * Creates a unique temporary directory on construction,
 * removes it (and all contents) on destruction.
 * Ensures cleanup happens even if tests fail.
 */
class TestDirectory
{
public:
    explicit TestDirectory(const std::string& name = "igor_streaming_tests")
        : m_path(std::filesystem::temp_directory_path() / name)
    {
        std::filesystem::create_directories(m_path);
    }

    ~TestDirectory()
    {
        std::error_code ec;
        std::filesystem::remove_all(m_path, ec);
#ifndef NDEBUG
        if (ec) {
            std::cerr << "[TestDirectory] Warning: Failed to cleanup test directory '"
                      << m_path.string() << "': " << ec.message() << std::endl;
        }
#endif
    }

    // Non-copyable, non-movable
    TestDirectory(const TestDirectory&) = delete;
    TestDirectory& operator=(const TestDirectory&) = delete;
    TestDirectory(TestDirectory&&) = delete;
    TestDirectory& operator=(TestDirectory&&) = delete;

    /// Get full path to a file within the test directory
    [[nodiscard]] std::string path(const std::string& filename) const
    {
        return (m_path / filename).string();
    }

    /// Get the root test directory path
    [[nodiscard]] const std::filesystem::path& root() const { return m_path; }

private:
    std::filesystem::path m_path;
};

//==============================================================================
// Test Data Generators
//==============================================================================

/**
 * @brief Create a simple test batch with sequence_id and sequence columns
 *
 * Creates a batch with 3 rows containing:
 * - sequence_id: [0, 1, 2]
 * - sequence: ["ATCGATCGATCG", "GCTAGCTAGCTA", "TTAATTAATTAA"]
 *
 * @return A sparrow::record_batch with 2 columns and 3 rows
 */
inline sparrow::record_batch create_test_batch()
{
    std::vector<std::string> names = {"sequence_id", "sequence"};
    std::vector<sparrow::array> arrays;

    arrays.emplace_back(sparrow::build(std::vector<int32_t>{0, 1, 2}));
    arrays.emplace_back(sparrow::build(std::vector<std::string>{
        "ATCGATCGATCG", "GCTAGCTAGCTA", "TTAATTAATTAA"
    }));

    return sparrow::record_batch(std::move(names), std::move(arrays));
}

/**
 * @brief Create test sequences without alignments
 *
 * Generates sequences with IDs from 0 to count-1. Each sequence starts
 * with "ATCGATCGATCGATCG" (16 bp). When `varied` is true, each sequence i
 * has (i % 10) additional "ACGT" suffixes appended, creating varied lengths.
 *
 * @param count Number of sequences to generate
 * @param varied If true, vary sequence length; if false, use fixed 16 bp length
 * @return Vector of SequenceTuple with empty alignment maps
 */
inline std::vector<SequenceTuple> create_test_sequences(
    size_t count,
    bool varied = true)
{
    std::vector<SequenceTuple> sequences;
    sequences.reserve(count);

    std::unordered_map<Gene_class, std::vector<Alignment_data>> empty_alignments;

    for (size_t i = 0; i < count; ++i) {
        std::string seq = "ATCGATCGATCGATCG";
        if (varied) {
            for (size_t j = 0; j < i % 10; ++j) {
                seq += "ACGT";
            }
        }
        sequences.emplace_back(static_cast<int>(i), std::move(seq), empty_alignments);
    }

    return sequences;
}

/**
 * @brief Create a test sequence with a pre-populated V gene alignment
 *
 * Creates a sequence with one V_gene alignment containing:
 * - gene_name: "IGHV1-1*01"
 * - offset: 10, five_p_offset: 5, three_p_offset: 20, align_length: 100
 * - insertions: [1, 2], deletions: [3, 4], mismatches: [5, 6]
 * - score: 123.45
 *
 * Useful for testing alignment round-trip serialization.
 *
 * @param id Sequence identifier
 * @param seq Nucleotide sequence string
 * @return SequenceTuple with the specified id, sequence, and V alignment
 */
inline SequenceTuple create_sequence_with_v_alignment(int id, const std::string& seq)
{
    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

    Alignment_data v_align(
        "IGHV1-1*01",                     // gene_name
        10,                                // offset
        5,                                 // five_p_offset
        20,                                // three_p_offset
        100,                               // align_length
        std::forward_list<int>{1, 2},      // insertions
        std::forward_list<int>{3, 4},      // deletions
        std::vector<int>{5, 6},            // mismatches
        123.45                             // score
    );
    alignments[V_gene].push_back(v_align);

    return {id, seq, alignments};
}

/**
 * @brief Load Murugan test dataset from CSV files
 *
 * Loads 300 sequences with ~73K alignments from the test data directory.
 * This is a real-world dataset for integration testing.
 *
 * @note Requires IGOR_TEST_DATA_DIR to be defined
 */
std::vector<SequenceTuple> load_murugan_dataset();

//==============================================================================
// Comparison Helpers
//==============================================================================

/**
 * @brief Compare two Alignment_data structures for field-by-field equality
 *
 * Compares all fields including gene_name, offset, score, five_p_offset,
 * three_p_offset, align_length, insertions, deletions, and mismatches.
 * Forward lists are converted to vectors for comparison.
 *
 * @param a First alignment to compare
 * @param b Second alignment to compare
 * @return true if all fields are equal, false otherwise
 */
inline bool alignments_equal(
    const Alignment_data& a,
    const Alignment_data& b)
{
    if (a.gene_name != b.gene_name) return false;
    if (a.offset != b.offset) return false;
    if (a.score != b.score) return false;
    if (a.five_p_offset != b.five_p_offset) return false;
    if (a.three_p_offset != b.three_p_offset) return false;
    if (a.align_length != b.align_length) return false;

    // Compare forward_lists
    std::vector<int> a_ins(a.insertions.begin(), a.insertions.end());
    std::vector<int> b_ins(b.insertions.begin(), b.insertions.end());
    if (a_ins != b_ins) return false;

    std::vector<int> a_del(a.deletions.begin(), a.deletions.end());
    std::vector<int> b_del(b.deletions.begin(), b.deletions.end());
    if (a_del != b_del) return false;

    if (a.mismatches != b.mismatches) return false;

    return true;
}

} // namespace igor::test
