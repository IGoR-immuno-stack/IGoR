/*
 * AlignerTestUtils.cpp
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

#include "AlignerTestUtils.h"

#include <cmath>
#include <cstdio>
#include <limits>

namespace igor::test {
namespace align {

// ============================================================================
// Test Fixtures
// ============================================================================

Matrix<double> build_test_score_matrix(double match_score, double mismatch_score)
{
    std::vector<double> values(15 * 15, mismatch_score);
    for (int i = 0; i < 15; ++i) {
        values[i + 15 * i] = match_score;
    }
    return Matrix<double>(15, 15, values);
}

Aligner make_legacy_aligner(const Matrix<double> &matrix, int gap_penalty, Gene_class gene_class,
                            const std::vector<std::pair<std::string, std::string>> &genomic_templates)
{
    Aligner aligner(matrix, gap_penalty, gene_class);
    aligner.set_genomic_sequences(genomic_templates);
    return aligner;
}

// ============================================================================
// Lookup Helpers
// ============================================================================

std::string find_germline_sequence(const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const std::string &gene_name)
{
    for (const auto &entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry.second;
        }
    }
    return std::string();
}

std::pair<std::string, std::string> find_genomic_template(
    const std::vector<std::pair<std::string, std::string>> &genomic_templates,
    const std::string &gene_name)
{
    for (const auto &entry : genomic_templates) {
        if (entry.first == gene_name) {
            return entry;
        }
    }
    return std::pair<std::string, std::string>{ };
}

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<int> sorted_list(std::forward_list<int> xs)
{
    std::vector<int> out(xs.begin(), xs.end());
    std::sort(out.begin(), out.end());
    return out;
}

bool check_alignment_data_equal(const Alignment_data &a, const Alignment_data &b, double score_tolerance)
{
    using namespace std;
    // Check basic fields
    REQUIRE(a.gene_name == b.gene_name);
    REQUIRE(a.offset == b.offset);
    REQUIRE(a.five_p_offset == b.five_p_offset);
    REQUIRE(a.three_p_offset == b.three_p_offset);
    REQUIRE(a.align_length == b.align_length);
    REQUIRE(fabs(a.score - b.score) <= score_tolerance);
    
    // Check insertions (convert to sets for order-independent comparison)
    unordered_set<int> a_ins(a.insertions.begin(), a.insertions.end());
    unordered_set<int> b_ins(b.insertions.begin(), b.insertions.end());
    REQUIRE(a_ins == b_ins);
    
    // Check deletions
    unordered_set<int> a_del(a.deletions.begin(), a.deletions.end());
    unordered_set<int> b_del(b.deletions.begin(), b.deletions.end());
    REQUIRE(a_del == b_del);
    
    // Check mismatches (already sorted, but compare as sets to be safe)
    unordered_set<int> a_mis(a.mismatches.begin(), a.mismatches.end());
    unordered_set<int> b_mis(b.mismatches.begin(), b.mismatches.end());
    REQUIRE(a_mis == b_mis);
    
    return true;
}

// ============================================================================
// Normalization and Extraction
// ============================================================================

std::vector<ActualAlignment> normalize_actuals(const std::forward_list<Alignment_data> &alignments)
{
    std::vector<ActualAlignment> normalized;
    for (const auto &aln : alignments) {
        normalized.push_back({ aln.gene_name, alignment_data_to_core_cigar(aln), aln.score });
    }
    std::sort(normalized.begin(), normalized.end(), [](const ActualAlignment &a, const ActualAlignment &b) {
        return std::tie(a.gene_name, a.cigar, a.score) < std::tie(b.gene_name, b.cigar, b.score);
    });
    return normalized;
}

std::vector<ExpectedAlignment> normalize_expected(std::vector<ExpectedAlignment> expected)
{
    std::sort(expected.begin(), expected.end(), [](const ExpectedAlignment &a, const ExpectedAlignment &b) {
        return std::tie(a.gene_name, a.core_cigar, a.extended_cigar, a.score) < std::tie(b.gene_name, b.core_cigar, b.extended_cigar, b.score);
    });
    return expected;
}

ActualAlignment best_alignment(const std::forward_list<Alignment_data> &alignments)
{
    auto best_it = std::max_element(alignments.begin(), alignments.end(),
                                    [](const Alignment_data &a, const Alignment_data &b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());
    return { best_it->gene_name, alignment_data_to_core_cigar(*best_it), best_it->score };
}

// ============================================================================
// Visualization Helpers
// ============================================================================

std::string cigar_visual(const std::string &cigar, const std::string &read, const std::string &germline)
{
    std::string read_gapped;
    std::string match_line;
    std::string germ_gapped;
    std::size_t read_i = 0;
    std::size_t germ_i = 0;

    auto push_pair = [&](char rc, char gc, char op) {
        read_gapped.push_back(rc);
        germ_gapped.push_back(gc);
        if (rc == '-' || gc == '-') {
            match_line.push_back(' ');
        } else if (op == '=') {
            match_line.push_back('|');
        } else if (op == 'X') {
            match_line.push_back('*');
        } else {
            match_line.push_back((rc == gc) ? '|' : '*');
        }
    };

    for (const auto &entry : parse_cigar(cigar)) {
        const int count = entry.first;
        const char op = entry.second;
        for (int i = 0; i < count; ++i) {
            switch (op) {
            case '=':
            case 'X':
            case 'M': {
                const char rc = (read_i < read.size()) ? read[read_i] : '?';
                const char gc = (germ_i < germline.size()) ? germline[germ_i] : '?';
                push_pair(rc, gc, op);
                ++read_i;
                ++germ_i;
                break;
            }
            case 'I':
            case 'S': {
                const char rc = (read_i < read.size()) ? read[read_i] : '?';
                push_pair(rc, '-', op);
                ++read_i;
                break;
            }
            case 'D':
            case 'N': {
                const char gc = (germ_i < germline.size()) ? germline[germ_i] : '?';
                push_pair('-', gc, op);
                ++germ_i;
                break;
            }
            case 'H':
            case 'P':
            default:
                break;
            }
        }
    }

    std::ostringstream out;
    const std::size_t width = 100;
    for (std::size_t i = 0; i < read_gapped.size(); i += width) {
        out << "\n      read " << read_gapped.substr(i, width) << "\n           " << match_line.substr(i, width)
            << "\n      germ " << germ_gapped.substr(i, width);
    }
    return out.str();
}

std::string add_cigar_visual_info(const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  const std::string &label, const std::string &gene_name, const std::string &cigar)
{
    const std::string germline = find_germline_sequence(genomic_templates, gene_name);
    if (germline.empty()) {
        std::ostringstream msg;
        msg << label << " gene=" << gene_name << " has no matching germline sequence for visual explanation";
        return msg.str();
    }
    std::ostringstream msg;
    msg << label << " gene=" << gene_name << " cigar=" << cigar << cigar_visual(cigar, query, germline);
    return msg.str();
}

std::string alignment_data_visual(const Alignment_data &aln, const std::string &query, const std::string &germline)
{
    const std::string cigar = alignment_data_to_extended_cigar(aln, query.length(), germline.length());
    return cigar_visual(cigar, query, germline);
}

// ============================================================================
// Test Assertions
// ============================================================================

void assert_alignment_matches(const Alignment_data &alignment, const std::string &query,
                              const std::pair<std::string, std::string> &genomic_template,
                              const ExpectedAlignment &expected)
{
    const std::string actual_core_cigar = alignment_data_to_core_cigar(alignment, query.length(), genomic_template.second.length());
    const std::string actual_extended_cigar = alignment_data_to_extended_cigar(alignment, query.length(), genomic_template.second.length());

    INFO("expected: gene=" << expected.gene_name << " core_cigar=" << expected.core_cigar << " extended_cigar=" << expected.extended_cigar << " score=" << expected.score);
    INFO("actual: gene=" << alignment.gene_name << " core_cigar=" << actual_core_cigar << " extended_cigar=" << actual_extended_cigar << " score=" << alignment.score);


    const std::vector<std::pair<std::string, std::string>> templates = { genomic_template };
    INFO("=== CORE CIGAR COMPARISON ===");
    INFO(add_cigar_visual_info(query, templates, "expected", expected.gene_name, expected.core_cigar));
    INFO(add_cigar_visual_info(query, templates, "actual", alignment.gene_name, actual_core_cigar));
    INFO("=== EXTENDED CIGAR COMPARISON ===");
    INFO(add_cigar_visual_info(query, templates, "expected", expected.gene_name, expected.extended_cigar));
    INFO(add_cigar_visual_info(query, templates, "actual", alignment.gene_name, actual_extended_cigar));

    REQUIRE(alignment.gene_name == expected.gene_name);
    REQUIRE(actual_core_cigar == expected.core_cigar);
    REQUIRE(actual_extended_cigar == expected.extended_cigar);
    REQUIRE_THAT(alignment.score, WithinRel(expected.score));

    // Redundantly assert match in Alignment_data struct space
    Alignment_data expected_aln = alignment_data_from_cigar_and_extended(expected.gene_name, expected.core_cigar,
                                                                         expected.extended_cigar, expected.score);
    check_alignment_data_equal(expected_aln, alignment);
}

void assert_best_alignment_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                   const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                   const ExpectedAlignment &expected)
{
    REQUIRE(!alignments.empty());
    auto best_it = std::max_element(alignments.begin(), alignments.end(),
                                    [](const Alignment_data &a, const Alignment_data &b) { return a.score < b.score; });
    REQUIRE(best_it != alignments.end());

    auto genomic_template = find_genomic_template(genomic_templates, expected.gene_name);
    if (genomic_template.first.empty()) {
        genomic_template = find_genomic_template(genomic_templates, best_it->gene_name);
    }
    assert_alignment_matches(*best_it, query, genomic_template, expected);
}

void assert_alignment_set_matches(const std::forward_list<Alignment_data> &alignments, const std::string &query,
                                  const std::vector<std::pair<std::string, std::string>> &genomic_templates,
                                  std::vector<ExpectedAlignment> expected)
{
    std::vector<const Alignment_data *> actual;
    for (const auto &aln : alignments) {
        actual.push_back(&aln);
    }
    std::sort(actual.begin(), actual.end(), [](const Alignment_data *a, const Alignment_data *b) {
        return std::make_tuple(a->gene_name, alignment_data_to_core_cigar(*a), a->score)
                < std::make_tuple(b->gene_name, alignment_data_to_core_cigar(*b), b->score);
    });
    expected = normalize_expected(std::move(expected));

    INFO("actual count=" << actual.size() << " expected count=" << expected.size());
    REQUIRE(actual.size() == expected.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
        INFO("alignment index=" << i);
        auto genomic_template = find_genomic_template(genomic_templates, expected[i].gene_name);
        if (genomic_template.first.empty()) {
            genomic_template = find_genomic_template(genomic_templates, actual[i]->gene_name);
        }
        assert_alignment_matches(*actual[i], query, genomic_template, expected[i]);
    }
}

void assert_alignment_data_matches(const Alignment_data &actual,
                                   const std::string &expected_csv_line,
                                   const std::string &query,
                                   const std::vector<std::pair<std::string, std::string>> &genomic_templates)
{
    const auto parsed = parse_single_alignment_csv_line(expected_csv_line);
    const Alignment_data expected = parsed.second;
    auto genomic_template = find_genomic_template(genomic_templates, expected.gene_name);

    const std::string actual_extended_cigar = alignment_data_to_extended_cigar(actual, query.length(), genomic_template.second.length());
    const std::string expected_extended_cigar = alignment_data_to_extended_cigar(expected, query.length(), genomic_template.second.length());

    INFO("expected: gene=" << expected.gene_name << " extended_cigar=" << expected_extended_cigar << " score=" << expected.score);
    INFO("actual: gene=" << actual.gene_name << " extended_cigar=" << actual_extended_cigar << " score=" << actual.score);

    std::vector<std::pair<std::string, std::string>> templates = { genomic_template };
    INFO("expected: gene=" << expected.gene_name << " extended_cigar=" << expected_extended_cigar 
         << cigar_visual(expected_extended_cigar, query, genomic_template.second));
    INFO("actual: gene=" << actual.gene_name << " extended_cigar=" << actual_extended_cigar 
         << cigar_visual(actual_extended_cigar, query, genomic_template.second));

    REQUIRE(actual.gene_name == expected.gene_name);
    REQUIRE(actual_extended_cigar == expected_extended_cigar);
    REQUIRE(align_compare(actual, expected));
    REQUIRE_THAT(actual.score, WithinRel(expected.score));
}

} // namespace Aligner
} // namespace IgorTestUtils
