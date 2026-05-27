/*
 * AIRRAlignmentWriter.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 *  analyze and model immune receptors generation, selection, mutation and all
 *  other processes.
 *   Copyright (C) 2026  IGoR Development Team
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

#include <igor/Streaming/AIRRAlignmentWriter.h>

#include <sparrow/array.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>

namespace igor::airr::alignment {

//==============================================================================
// Internal helpers
//==============================================================================

namespace {

// Convert Gene_class to segment string
std::string gene_class_to_segment(Gene_class gc)
{
    switch (gc) {
        case Gene_class::V_gene: return "V";
        case Gene_class::D_gene: return "D";
        case Gene_class::J_gene: return "J";
        default: return "";
    }
}

} // anonymous namespace

//==============================================================================
// Public API
//==============================================================================

// Use shared utility function from parent namespace
using airr::delimiter_char;

std::string make_cigar(const Alignment_data& align)
{
    std::ostringstream cigar;

    // Count insertions and deletions
    size_t num_insertions = std::distance(align.insertions.begin(), align.insertions.end());
    size_t num_deletions = std::distance(align.deletions.begin(), align.deletions.end());

    // Generate simplified CIGAR: matches first, then insertions, then deletions
    if (align.align_length > 0) {
        cigar << align.align_length << "M";
    }
    if (num_insertions > 0) {
        cigar << num_insertions << "I";
    }
    if (num_deletions > 0) {
        cigar << num_deletions << "D";
    }

    return cigar.str();
}

void write_sequences(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences,
    Delimiter delimiter)
{
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::alignment::write_sequences: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Write header
    file << "sequence_id" << delim
         << "segment" << delim
         << "call" << delim
         << "score" << delim
         << "sequence_start" << delim
         << "sequence_end" << delim
         << "germline_start" << delim
         << "cigar" << delim
         << "rank\n";

    // Write alignments (one row per alignment)
    for (const auto& seq : sequences) {
        // Iterate over all gene classes
        for (const auto& [gene_class, align_vec] : seq.alignments) {
            // Skip if no alignments for this gene class
            if (align_vec.empty()) continue;

            // Get segment string
            std::string segment = gene_class_to_segment(gene_class);
            if (segment.empty()) continue; // Skip unsupported gene classes

            // Write each alignment with rank
            size_t rank = 1;
            for (const auto& align : align_vec) {
                file << seq.index << delim
                     << segment << delim
                     << align.gene_name << delim
                     << align.score << delim
                     << (align.offset + 1) << delim  // 0-based → 1-based
                     << (align.offset + static_cast<int>(align.align_length)) << delim
                     << (align.five_p_offset + 1) << delim  // 0-based → 1-based
                     << make_cigar(align) << delim
                     << rank << "\n";
                ++rank;
            }
        }
    }
}

void write_batch(
    const std::string& filepath,
    const sparrow::record_batch& batch,
    Delimiter delimiter)
{
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::alignment::write_batch: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Get column names from batch using proper API
    const auto& names = batch.names();
    std::vector<std::string> column_names(names.begin(), names.end());

    // Write header
    for (size_t i = 0; i < column_names.size(); ++i) {
        if (i > 0) file << delim;
        file << column_names[i];
    }
    file << "\n";

    // Write data rows
    size_t num_rows = batch.nb_rows();
    for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < column_names.size(); ++col) {
            if (col > 0) file << delim;
            file << get_string_value(batch, column_names[col], row);
        }
        file << "\n";
    }
}

} // namespace igor::airr::alignment
