/*
 * AlignmentData.h
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

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <igor/Core/Export.h>

namespace igor::core {

/**
 * \class AlignmentData Aligner.h
 * \brief Stores information on the alignment of one genomic template against the target
 * \author Q.Marcou
 * \version 1.0
 *
 * Stores information on the alignment of one genomic template against the target.
 * It contains:
 * - the gene name
 * - the offset of the alignment (index on the target sequence on which the first letter of the FULL genomic template aligns (can be negative or lie outside the target)
 * - 5' and 3' offset positions of the best alignment first and last aligned nucleotide
 * - insertions : indices on the TARGET of inserted nucleotides
 * - deletions : indices on the GENOMIC TEMPLATE of deleted nucleotides
 * - alignment length
 * - list of mismatches (that lie even outside the best alignment to allow IGoR to know mismatch positions in advance while exploring different deletions numbers)
 *   NOTE: The mismatches vector may contain mismatches beyond the five_p_offset to three_p_offset range,
 *   representing mismatches in the extended alignment regions (e.g., in deleted V/J nucleotides).
 *   The vector is SORTED and contains ALL mismatches (both within the core alignment and in extended regions).
 * - the alignment score
 * - query_length and germline_length: lengths of the query/target and reference/genomic sequences (0 = unknown)
 * - extended_mismatches: mismatches outside the core alignment region [five_p_offset, three_p_offset]
 * - extended_insertions/deletions: indels in extended regions (reserved for future use, currently empty)
 *
 * Computed properties (via getter methods):
 * - query_align_start/end: core alignment bounds in query sequence
 * - reference_align_start/end: core alignment bounds in reference sequence
 * - extended_query/reference_align_start/end: extended alignment bounds
 * - core_cigar/extended_cigar: CIGAR string representations
 * - get_all/mismatches: all mismatches (core + extended)
 * - get_core_mismatches: mismatches within the core alignment region
 * - get_5p/3p_extended_mismatches: extended mismatches by side
 * - get_extended_mismatches: all extended mismatches
 * - get_core/extended/all_insertions: insertion accessors
 * - get_core/extended/all_deletions: deletion accessors
 */
struct CORE_EXPORT AlignmentData
{
    std::string gene_name;
    int offset;
    size_t five_p_offset;
    size_t three_p_offset;
    std::vector<size_t> insertions; //gap in the genomic sequence
    std::vector<size_t> deletions; //gap in the data sequence
    size_t align_length;
    mutable std::vector<size_t> mismatches;
    double score;
    
    // Sequence lengths (0 = unknown)
    size_t query_length = 0;
    size_t germline_length = 0;
    
    // Extended region tracking
    std::vector<size_t> extended_mismatches;  // Mismatches outside [five_p_offset, three_p_offset]
    std::vector<size_t> extended_insertions;  // Insertions in extended regions (future use)
    std::vector<size_t> extended_deletions;  // Deletions in extended regions (future use)

    // Computed property getters
    
    // Core alignment bounds (query = target sequence, reference = genomic template)
    size_t query_align_start() const { return five_p_offset; }
    size_t query_align_end() const { return three_p_offset; }
    size_t reference_align_start() const
    {
        size_t n_ins_5p = get_5p_extended_insertions().size();
        size_t n_del_5p = get_5p_extended_deletions().size();
        return static_cast<size_t>(five_p_offset - offset - n_ins_5p + n_del_5p);
    }
    size_t reference_align_end() const
    {
        size_t n_ins_all = get_5p_extended_insertions().size() + get_core_insertions().size();
        size_t n_del_all = get_5p_extended_deletions().size() + get_core_deletions().size();
        return static_cast<size_t>(three_p_offset - offset - n_ins_all + n_del_all);
    }

    // Extended alignment bounds (clamped to sequence lengths)
    size_t extended_query_align_start() const
    {
        if (query_length == 0)
            throw std::logic_error("query_length required for extended bounds");
        return static_cast<size_t>((std::max)(offset, 0));
    }
    size_t extended_query_align_end() const
    {
        if (query_length == 0 || germline_length == 0)
            throw std::logic_error("query_length and germline_length required for extended bounds");
        size_t n_del = get_all_deletions().size();
        size_t n_ins = get_all_insertions().size();
        return static_cast<size_t>((std::min)(offset + germline_length - 1 + n_ins - n_del, query_length - 1));
    }
    size_t extended_reference_align_start() const
    {
        if (germline_length == 0)
            throw std::logic_error("germline_length required for extended bounds");
        return static_cast<size_t>((std::max)(-offset, 0)); // Full germline reference
    }
    size_t extended_reference_align_end() const
    {
        if (germline_length == 0)
            throw std::logic_error("germline_length required for extended bounds");
        size_t n_del = get_all_deletions().size();
        size_t n_ins = get_all_insertions().size();
        return static_cast<size_t>((std::min)(query_length - 1 - n_ins + n_del - offset,
                                              germline_length - 1)); // Full germline reference
    }

    // Mismatch and indel access
    const std::vector<size_t>& get_all_mismatches() const { return mismatches; }
    std::vector<size_t> get_core_mismatches() const;
    std::vector<size_t> get_5p_extended_mismatches() const;
    std::vector<size_t> get_3p_extended_mismatches() const;
    std::vector<size_t> get_extended_mismatches() const { return extended_mismatches; }

    const std::vector<size_t>& get_all_insertions() const { return insertions; }
    std::vector<size_t> get_core_insertions() const;
    std::vector<size_t> get_5p_extended_insertions() const;
    std::vector<size_t> get_3p_extended_insertions() const;
    std::vector<size_t> get_extended_insertions() const { return extended_insertions; }

    const std::vector<size_t>& get_all_deletions() const { return deletions; }
    std::vector<size_t> get_core_deletions() const;
    std::vector<size_t> get_5p_extended_deletions() const;
    std::vector<size_t> get_3p_extended_deletions() const;
    std::vector<size_t> get_extended_deletions() const { return extended_deletions; }

    
    // CIGAR string access
    std::string core_cigar() const;
    std::string extended_cigar() const;
    
    // Validation
    bool validate() const;

    AlignmentData(std::string gene, int off, size_t query_len = 0, size_t germline_len = 0)
        : gene_name(gene),
          offset(off),
          score(0),
          query_length(query_len),
          germline_length(germline_len)
    {
    }
    AlignmentData(int off, size_t five_p_off, size_t three_p_off, size_t align_len, std::vector<size_t> ins,
                   std::vector<size_t> del, std::vector<size_t> mis, double alignment_score,
                   size_t query_len = 0, size_t germline_len = 0)
        : gene_name(std::string()),
          offset(off),
          five_p_offset(five_p_off),
          three_p_offset(three_p_off),
          insertions(std::move(ins)),
          deletions(std::move(del)),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score),
          query_length(query_len),
          germline_length(germline_len)
    {
        // Insertions/deletions are stored sorted ascending (invariant relied upon by
        // get_core_deletions/get_5p_extended_deletions/get_3p_extended_deletions)
        std::sort(insertions.begin(), insertions.end());
        std::sort(deletions.begin(), deletions.end());
        // Populate extended_mismatches from mismatches
        for (int pos : mismatches) {
            if (pos < five_p_offset || pos > three_p_offset) {
                extended_mismatches.push_back(pos);
            }
        }
    }
    AlignmentData(std::string gene, int off, size_t align_len, std::vector<size_t> ins, std::vector<size_t> del,
                   std::vector<size_t> mis, double alignment_score,
                   size_t query_len = 0, size_t germline_len = 0)
        : gene_name(gene),
          offset(off),
          insertions(std::move(ins)),
          deletions(std::move(del)),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score),
          query_length(query_len),
          germline_length(germline_len)
    {
        std::sort(insertions.begin(), insertions.end());
        std::sort(deletions.begin(), deletions.end());
    }
    AlignmentData(std::string gene, int off, size_t five_p_off, size_t three_p_off, size_t align_len,
                   std::vector<size_t> ins, std::vector<size_t> del, std::vector<size_t> mis, double alignment_score,
                   size_t query_len = 0, size_t germline_len = 0)
        : gene_name(gene),
          offset(off),
          five_p_offset(five_p_off),
          three_p_offset(three_p_off),
          insertions(std::move(ins)),
          deletions(std::move(del)),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score),
          query_length(query_len),
          germline_length(germline_len)
    {
        std::sort(insertions.begin(), insertions.end());
        std::sort(deletions.begin(), deletions.end());
        // Populate extended_mismatches from mismatches
        for (int pos : mismatches) {
            if (pos < five_p_offset || pos > three_p_offset) {
                extended_mismatches.push_back(pos);
            }
        }
    }

    /*	bool operator<(const AlignmentData& align){
		//Hardcode to get the alignments in descending order using sort()
		return this->score > align.score;
	}*/
};

/// Splits a CIGAR string into (count, operator) pairs.
CORE_EXPORT std::vector<std::pair<int, char>> parseCigar(const std::string &cigar);

/// CIGAR of the core alignment region only, with =, X, D, I operators.
CORE_EXPORT std::string toCoreCigar(const AlignmentData &aln);

/// CIGAR of the core region spanning the full sequences, with AIRR N/S operators for end gaps.
CORE_EXPORT std::string toCoreCigar(const AlignmentData &aln, std::size_t sequence_length, std::size_t germline_length);

/// As toCoreCigar over the full span, extended to every mismatch position of the alignment.
CORE_EXPORT std::string toExtendedCigar(const AlignmentData &aln, std::size_t sequence_length,
                                        std::size_t germline_length);

} // namespace igor::core
