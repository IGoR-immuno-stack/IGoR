/*
 * AlignerTestUtils.h
 *
 * Test utilities for Aligner module tests.
 * 
 * This source code is distributed as part of the IGoR software.
 * IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 * generation, selection, mutation and all other processes.
 *  Copyright (C) 2017  Quentin Marcou
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Core/Aligner.h>
#include <igor/Core/Utils.h>

#include <algorithm>
#include <forward_list>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace igor::test {
namespace align {

using Catch::Matchers::WithinRel;

// ============================================================================
// Test Data Structures
// ============================================================================

/**
 * Expected alignment data for test assertions
 */
struct ExpectedAlignment
{
    std::string gene_name;
    std::string core_cigar;
    std::string extended_cigar;
    double score;
};

/**
 * Actual alignment data extracted from Alignment_data for comparison
 */
struct ActualAlignment
{
    std::string gene_name;
    std::string cigar;
    double score;
};

// ============================================================================
// Test Fixtures
// ============================================================================

/**
 * Build a test score matrix with specified match and mismatch scores.
 * All diagonal elements (self-matches) get match_score, others get mismatch_score.
 * 
 * \param match_score Score for matching nucleotides (default: 2.0)
 * \param mismatch_score Score for mismatches (default: -1.0)
 * \return Matrix<double> 15x15 substitution matrix
 */
Matrix<double> build_test_score_matrix(double match_score = 2.0, double mismatch_score = -1.0);

/**
 * Create a legacy aligner configured for testing.
 * 
 * \param matrix Substitution score matrix
 * \param gap_penalty Gap penalty value
 * \param gene_class Gene class (V_gene, D_gene, J_gene)
 * \param genomic_templates Vector of gene name and sequence pairs
 * \return Configured Aligner instance
 */
Aligner make_legacy_aligner(const Matrix<double> &matrix, int gap_penalty, Gene_class gene_class,
                            const std::vector<std::pair<std::string, std::string>> &genomic_templates);

// ============================================================================
// Lookup Helpers
// ============================================================================

/**
 * Find germline sequence by gene name from genomic templates.
 * 
 * \param genomic_templates Vector of gene name and sequence pairs
 * \param gene_name Name of gene to find
 * \return The germline sequence string, or empty string if not found
 */
std::string find_germline_sequence(const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const std::string &gene_name);

/**
 * Find genomic template by gene name from genomic templates.
 * 
 * \param genomic_templates Vector of gene name and sequence pairs
 * \param gene_name Name of gene to find
 * \return The genomic template pair, or empty pair if not found
 */
std::pair<std::string, std::string> find_genomic_template(
    const std::vector<std::pair<std::string, std::string>> &genomic_templates,
    const std::string &gene_name);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Convert a forward_list<int> to a sorted vector<int>.
 * 
 * \param xs The forward_list to convert
 * \return Sorted vector containing the same elements
 */
std::vector<int> sorted_list(std::forward_list<int> xs);

/**
 * Compare two Alignment_data objects for equality with detailed diagnostics.
 * Two alignments are considered equal if they have:
 * - Same gene name
 * - Same offset
 * - Same five_p_offset and three_p_offset
 * - Same insertions, deletions, mismatches (as sets, order-independent)
 * - Same align_length
 * - Same score (within tolerance)
 * 
 * This function uses Catch2 REQUIRE statements to provide detailed failure information
 * for each field that doesn't match, making test debugging much easier.
 * 
 * \param a First alignment to compare
 * \param b Second alignment to compare
 * \param score_tolerance Tolerance for score comparison (default: 1e-9)
 * \return true if alignments are considered equal, false otherwise
 */
bool check_alignment_data_equal(const Alignment_data &a, const Alignment_data &b, double score_tolerance = 1e-9);

// ============================================================================
// Normalization and Extraction
// ============================================================================

/**
 * Convert a list of Alignment_data to a sorted vector of ActualAlignment structs.
 * Useful for comparing sets of alignments regardless of order.
 * 
 * \param alignments The forward_list of Alignment_data to normalize
 * \return Sorted vector of ActualAlignment structs
 */
std::vector<ActualAlignment> normalize_actuals(const std::forward_list<Alignment_data> &alignments);

/**
 * Sort a vector of ExpectedAlignment structs.
 * 
 * \param expected The vector of ExpectedAlignment to sort
 * \return Sorted vector of ExpectedAlignment
 */
std::vector<ExpectedAlignment> normalize_expected(std::vector<ExpectedAlignment> expected);

/**
 * Find the best (highest-scoring) alignment from a list.
 * 
 * \param alignments The forward_list of Alignment_data to search
 * \return ActualAlignment representing the best alignment
 */
ActualAlignment best_alignment(const std::forward_list<Alignment_data> &alignments);

// ============================================================================
// Visualization Helpers
// ============================================================================

/**
 * Create a visual representation of a CIGAR alignment.
 * Shows the read, match line, and germline sequences with alignment symbols.
 * 
 * \param cigar The CIGAR string to visualize
 * \param read The query/read sequence
 * \param germline The reference/germline sequence
 * \return String containing the visual alignment representation
 */
std::string cigar_visual(const std::string &cigar, const std::string &read, const std::string &germline);

/**
 * Create visual CIGAR information for a specific gene alignment.
 * 
 * \param query The query sequence
 * \param genomic_templates Vector of genomic templates (gene name, sequence)
 * \param label Label to prefix the output (e.g., "expected", "actual")
 * \param gene_name Name of the gene being aligned
 * \param cigar The CIGAR string to visualize
 * \return String containing labeled visual information
 */
std::string add_cigar_visual_info(const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  const std::string &label, const std::string &gene_name, const std::string &cigar);

/**
 * Create visual representation of Alignment_data.
 * 
 * \param aln The Alignment_data to visualize
 * \param query The query sequence
 * \param germline The germline sequence
 * \return String containing the visual alignment representation
 */
std::string alignment_data_visual(const Alignment_data &aln, const std::string &query, const std::string &germline);

// ============================================================================
// Test Assertions
// ============================================================================

/**
 * Assert that an alignment matches expected values with detailed diagnostics.
 * Performs the following checks:
 * - Gene name matches
 * - Core CIGAR matches
 * - Extended CIGAR matches
 * - Score matches (within relative tolerance)
 * 
 * Provides detailed visual output for debugging failures.
 * 
 * \param alignment The actual Alignment_data from the aligner
 * \param query The query sequence used for alignment
 * \param genomic_template The genomic template pair (gene_name, sequence)
 * \param expected The expected alignment values
 */
void assert_alignment_matches(const Alignment_data &alignment, const std::string &query,
                              const std::pair<std::string, std::string> &genomic_template,
                              const ExpectedAlignment &expected);

/**
 * Assert that the best alignment in a list matches expected values.
 * 
 * \param alignments The forward_list of Alignment_data from the aligner
 * \param query The query sequence used for alignment
 * \param genomic_templates Vector of genomic templates
 * \param expected The expected alignment values
 */
void assert_best_alignment_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                   const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const ExpectedAlignment &expected);

/**
 * Assert that a set of alignments matches expected values.
 * Alignments are sorted and compared in order.
 * 
 * \param alignments The forward_list of Alignment_data from the aligner
 * \param query The query sequence used for alignment
 * \param genomic_templates Vector of genomic templates
 * \param expected Vector of ExpectedAlignment values
 */
void assert_alignment_set_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  std::vector<ExpectedAlignment> expected);

/**
 * Assert that alignment data matches expected values from a CSV line.
 * 
 * \param actual The actual Alignment_data from the aligner
 * \param expected_csv_line The expected alignment in CSV format
 * \param query The query sequence used for alignment
 * \param genomic_templates Vector of genomic templates
 */
void assert_alignment_data_matches(const Alignment_data &actual,
                                   const std::string &expected_csv_line,
                                   const std::string &query,
                                   const std::vector<std::pair<std::string, std::string>> &genomic_templates);

} // namespace Aligner
} // namespace IgorTestUtils
