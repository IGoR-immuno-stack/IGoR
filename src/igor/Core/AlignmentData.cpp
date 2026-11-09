/*
 * AlignmentData.cpp
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <igor/Core/AlignmentData.h>

#include <unordered_set>

// Moved verbatim from Alignment/Legacy/Aligner.cpp (step 1c); the using-declarations below keep
// the bodies as they were written.
using std::invalid_argument;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::sort;
using std::string;
using std::to_string;
using std::unordered_set;
using std::vector;

namespace igor::core {

namespace {

void append_cigar_op(vector<pair<int, char>> &ops, char op)
{
    if (!ops.empty() && ops.back().second == op) {
        ++ops.back().first;
    } else {
        ops.push_back(make_pair(1, op));
    }
}

string cigar_ops_to_string(const vector<pair<int, char>> &ops)
{
    string cigar_out;
    for (const auto &op : ops) {
        cigar_out += to_string(op.first);
        cigar_out += op.second;
    }
    return cigar_out;
}

void append_cigar_run(vector<pair<int, char>> &ops, int count, char op)
{
    if (count <= 0)
        return;
    if (!ops.empty() && ops.back().second == op)
        ops.back().first += count;
    else
        ops.push_back(make_pair(count, op));
}

} // namespace

std::vector<std::pair<int, char>> parseCigar(const std::string &cigar)
{
    std::vector<std::pair<int, char>> ops;
    if (cigar.empty()) {
        throw invalid_argument("empty CIGAR string");
    }
    size_t i = 0;
    while (i < cigar.size()) {
        if (!isdigit(static_cast<unsigned char>(cigar[i]))) {
            throw invalid_argument("CIGAR count expected");
        }
        long count = 0;
        while (i < cigar.size() && isdigit(static_cast<unsigned char>(cigar[i]))) {
            count = count * 10 + (cigar[i] - '0');
            ++i;
        }
        if (count <= 0 || i == cigar.size()) {
            throw invalid_argument("invalid CIGAR count");
        }
        char op = cigar[i++];
        const string valid_ops = "=XIDNSHPM";
        if (valid_ops.find(op) == string::npos) {
            throw invalid_argument("invalid CIGAR operation");
        }
        ops.push_back(make_pair(static_cast<int>(count), op));
    }
    return ops;
}

/**
 * Convert AlignmentData to core CIGAR string (1-parameter overload).
 * This function produces a CIGAR string for just the core alignment region
 * between five_p_offset and three_p_offset, using =, X, D, I operators.
 */
std::string toCoreCigar(const AlignmentData &aln)
{
    vector<size_t> insertions = aln.get_all_insertions();
    vector<size_t> deletions = aln.get_all_deletions();
    sort(insertions.begin(), insertions.end());
    sort(deletions.begin(), deletions.end());
    unordered_set<int> mismatches(aln.get_all_mismatches().begin(), aln.get_all_mismatches().end());

    size_t ins_i = 0;
    size_t del_i = 0;
    int t = static_cast<int>(aln.five_p_offset);
    int g = static_cast<int>(aln.five_p_offset) - aln.offset;
    const int t_end = static_cast<int>(aln.three_p_offset);
    vector<pair<int, char>> ops;

    while (t <= t_end) {
        while (del_i < deletions.size() && deletions[del_i] == g) {
            append_cigar_op(ops, 'D');
            ++g;
            ++del_i;
        }
        if (ins_i < insertions.size() && insertions[ins_i] == t) {
            append_cigar_op(ops, 'I');
            ++t;
            ++ins_i;
        } else {
            append_cigar_op(ops, mismatches.count(t) ? 'X' : '=');
            ++t;
            ++g;
        }
    }
    while (del_i < deletions.size() && deletions[del_i] == g) {
        append_cigar_op(ops, 'D');
        ++g;
        ++del_i;
    }

    return cigar_ops_to_string(ops);
}

/**
 * Convert AlignmentData to core CIGAR string (3-parameter version).
 * This was renamed from alignment_data_to_cigar_full_span as requested.
 * This function produces a CIGAR string that includes the full span with
 * AIRR-compliant N/S operators for end gaps.
 */
std::string toCoreCigar(const AlignmentData &aln, size_t sequence_length, size_t germline_length)
{
    int t = static_cast<int>(aln.five_p_offset);
    int g = static_cast<int>(aln.five_p_offset) - aln.offset;
    vector<pair<int, char>> ops;
    // Use AIRR-standard operators: N for reference-only gaps, S for query-only gaps
    append_cigar_run(ops, g, 'N');  // Leading reference gaps
    append_cigar_run(ops, t, 'S');  // Leading query gaps
    
    // Reuse the core version for the core region
    for (const auto &entry : parseCigar(toCoreCigar(aln))) {
        append_cigar_run(ops, entry.first, entry.second);
        switch (entry.second) {
        case '=':
        case 'X':
        case 'M':
            t += entry.first;
            g += entry.first;
            break;
        case 'I':
        case 'S':
            t += entry.first;
            break;
        case 'D':
        case 'N':
            g += entry.first;
            break;
        default:
            break;
        }
    }
    // Use AIRR-standard operators: N for reference-only gaps, S for query-only gaps
    append_cigar_run(ops, static_cast<int>(germline_length) - g, 'N');  // Trailing reference gaps
    append_cigar_run(ops, static_cast<int>(sequence_length) - t, 'S');  // Trailing query gaps
    return cigar_ops_to_string(ops);
}

/**
 * Convert AlignmentData to extended CIGAR string including extended alignment regions.
 * This function creates a comprehensive CIGAR representation that includes:
 * - Leading reference gaps as N operators
 * - Leading query gaps as S operators  
 * - Core alignment with =, X, D, I operators (extended to include extended mismatches)
 * - Trailing reference gaps as N operators
 * - Trailing query gaps as S operators
 * 
 * Unlike alignment_data_to_core_cigar which only processes the core alignment region,
 * this function extends the alignment to include ALL mismatch positions from alignment_data.mismatches,
 * representing extended mismatches as X operators rather than as gaps.
 */
std::string toExtendedCigar(const AlignmentData &aln, size_t sequence_length, size_t germline_length)
{
    AlignmentData aln_copy = aln;
    // Update offset values to extend alignment bounds
    aln_copy.five_p_offset = max(0, aln.offset);
    size_t n_del = aln.get_all_deletions().size();
    size_t n_ins = aln.get_all_insertions().size();
    aln_copy.three_p_offset = min(aln_copy.offset + germline_length - 1 + n_ins - n_del, sequence_length - 1);
    aln_copy.align_length = aln_copy.three_p_offset - aln_copy.five_p_offset + n_del;
    return toCoreCigar(aln_copy, sequence_length, germline_length);
}

// ============================================================================
// AlignmentData computed getter method implementations
// ============================================================================

std::vector<size_t> AlignmentData::get_core_mismatches() const {
    std::vector<size_t> core;
    for (int pos : mismatches) {
        if (pos >= five_p_offset && pos <= three_p_offset) {
            core.push_back(pos);
        }
    }
    return core;
}

std::vector<size_t> AlignmentData::get_5p_extended_mismatches() const {
    std::vector<size_t> result;
    for (int pos : extended_mismatches) {
        if (pos < static_cast<int>(five_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<size_t> AlignmentData::get_3p_extended_mismatches() const {
    std::vector<size_t> result;
    for (int pos : extended_mismatches) {
        if (pos > static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<size_t> AlignmentData::get_core_insertions() const
{
    std::vector<size_t> result;
    for (int pos : get_all_insertions()) {
        if (pos >= static_cast<int>(five_p_offset) && pos <= static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<size_t> AlignmentData::get_core_deletions() const
{
    std::vector<size_t> result;
    size_t n_ins_5p = get_5p_extended_insertions().size();
    size_t n_ins_all = get_5p_extended_insertions().size() + get_core_insertions().size();
    int ref_start = static_cast<int>(five_p_offset) - offset - n_ins_5p;
    int ref_end = static_cast<int>(three_p_offset) - offset - n_ins_all;
    for (int pos : get_all_deletions()) {
        if (pos < ref_start) {
            ref_start++;
            ref_end++;
        } else if (pos > ref_start && pos <= ref_end) {
            result.push_back(pos);
            ref_end++; // Account for reference coordinate shift
        } else {
            // Beyond core alignement
            break;
        }
    }
    return result;
}

// Extended insertion accessors (query coordinates)
std::vector<size_t> AlignmentData::get_5p_extended_insertions() const {
    std::vector<size_t> result;
    for (int pos : get_all_insertions()) {
        if (pos < static_cast<int>(five_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<size_t> AlignmentData::get_3p_extended_insertions() const {
    std::vector<size_t> result;
    for (int pos : get_all_insertions()) {
        if (pos > static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

// Extended deletion accessors (reference coordinates with dynamic threshold adjustment)
std::vector<size_t> AlignmentData::get_5p_extended_deletions() const
{
    std::vector<size_t> result;
    size_t n_ins_5p = get_5p_extended_insertions().size();
    int ref_start = static_cast<int>(five_p_offset) - offset - n_ins_5p;
    for (int pos : get_all_deletions()) {
        if (pos <= ref_start) {
            result.push_back(pos);
            ref_start++; // Account for reference coordinate shift
        }
    }
    return result;
}

std::vector<size_t> AlignmentData::get_3p_extended_deletions() const {
    std::vector<size_t> result;
    size_t n_ins_all = get_5p_extended_insertions().size() + get_core_insertions().size();
    int ref_end = static_cast<int>(three_p_offset) - offset - n_ins_all;
    
    // Second pass: collect 3' extended deletions (after adjusted ref_end)
    for (int pos : get_all_deletions()) {
        if (pos > ref_end) {
            result.push_back(pos);
        }
        else{
            ref_end++;
        }
    }
    return result;
}

std::string AlignmentData::core_cigar() const {
    if (query_length == 0 || germline_length == 0) {
        throw std::logic_error("query_length and germline_length required for CIGAR generation");
    }
    return toCoreCigar(*this, query_length, germline_length);
}

std::string AlignmentData::extended_cigar() const {
    if (query_length == 0 || germline_length == 0) {
        throw std::logic_error("query_length and germline_length required for CIGAR generation");
    }
    return toExtendedCigar(*this, query_length, germline_length);
}

bool AlignmentData::validate() const {
    // Check that mismatches are sorted
    for (size_t i = 1; i < mismatches.size(); ++i) {
        if (mismatches[i-1] > mismatches[i]) {
            return false; // Not sorted
        }
    }

    // Check that insertions/deletions are sorted (relied upon by
    // get_core_deletions/get_5p_extended_deletions/get_3p_extended_deletions)
    for (size_t i = 1; i < insertions.size(); ++i) {
        if (insertions[i-1] > insertions[i]) {
            return false; // Not sorted
        }
    }
    for (size_t i = 1; i < deletions.size(); ++i) {
        if (deletions[i-1] > deletions[i]) {
            return false; // Not sorted
        }
    }

    // Check that extended_mismatches are sorted
    for (size_t i = 1; i < extended_mismatches.size(); ++i) {
        if (extended_mismatches[i-1] > extended_mismatches[i]) {
            return false; // Not sorted
        }
    }
    
    // Check that extended_mismatches are actually outside core alignment
    for (int pos : extended_mismatches) {
        if (pos >= static_cast<int>(five_p_offset) && pos <= static_cast<int>(three_p_offset)) {
            return false; // Extended mismatch is within core alignment
        }
    }
    
    // Check that all mismatches in extended_mismatches are also in mismatches
    std::unordered_set<int> all_mismatches_set(mismatches.begin(), mismatches.end());
    for (int pos : extended_mismatches) {
        if (all_mismatches_set.count(pos) == 0) {
            return false; // Extended mismatch not in main mismatches
        }
    }
    
    // Check that core mismatches + extended mismatches = all mismatches
    std::vector<size_t> core = get_core_mismatches();
    std::vector<size_t> combined = core;
    combined.insert(combined.end(), extended_mismatches.begin(), extended_mismatches.end());
    std::sort(combined.begin(), combined.end());
    if (combined != mismatches) {
        return false; // Inconsistent mismatch categorization
    }
    
    return true;
}

} // namespace igor::core
