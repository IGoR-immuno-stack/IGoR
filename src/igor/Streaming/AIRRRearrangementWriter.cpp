/*
 * AIRRRearrangementWriter.cpp
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

#include <igor/Streaming/AIRRRearrangementWriter.h>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iterator>

namespace igor::airr::rearrangement {

//==============================================================================
// Internal helpers
//==============================================================================

namespace {

/**
 * @brief Get alignment for a specific gene class, or nullptr if not present
 */
const Alignment_data* get_alignment(
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>& alignments,
    Gene_class gene_class)
{
    auto it = alignments.find(gene_class);
    if (it == alignments.end() || it->second.empty()) {
        return nullptr;
    }
    return &it->second[0];  // Return first alignment
}

/**
 * @brief Write alignment fields for a gene class
 */
void write_alignment_fields(
    std::ostream& out,
    const Alignment_data* align,
    char delim)
{
    if (align) {
        // gene_call
        out << delim << align->gene_name;
        // score
        out << delim << align->score;
        // sequence_start (convert 0-based to 1-based)
        out << delim << (align->offset + 1);
        // alignment_length
        out << delim << align->align_length;
        // cigar
        out << delim << make_cigar(*align);
    } else {
        // Empty fields for missing alignment
        out << delim << ""   // call
            << delim << ""   // score
            << delim << ""   // start
            << delim << ""   // length
            << delim << "";  // cigar
    }
}

/**
 * @brief Write header row
 */
void write_header(std::ostream& out, char delim)
{
    out << "sequence_id" << delim << "sequence";

    // V gene columns
    out << delim << "v_call" << delim << "v_score"
        << delim << "v_sequence_start" << delim << "v_alignment_length"
        << delim << "v_cigar";

    // D gene columns
    out << delim << "d_call" << delim << "d_score"
        << delim << "d_sequence_start" << delim << "d_alignment_length"
        << delim << "d_cigar";

    // J gene columns
    out << delim << "j_call" << delim << "j_score"
        << delim << "j_sequence_start" << delim << "j_alignment_length"
        << delim << "j_cigar";

    out << "\n";
}

/**
 * @brief Write a single sequence row
 */
void write_sequence_row(
    std::ostream& out,
    const SequenceData& seq,
    char delim)
{
    // sequence_id and sequence
    out << seq.index << delim << seq.sequence;

    // V, D, J alignments
    write_alignment_fields(out, get_alignment(seq.alignments, V_gene), delim);
    write_alignment_fields(out, get_alignment(seq.alignments, D_gene), delim);
    write_alignment_fields(out, get_alignment(seq.alignments, J_gene), delim);

    out << "\n";
}

} // anonymous namespace

//==============================================================================
// Public API
//==============================================================================

std::string make_cigar(const Alignment_data& alignment)
{
    if (alignment.align_length == 0) {
        return "";
    }

    std::ostringstream cigar;

    // Count insertions and deletions using std::distance
    size_t num_insertions = static_cast<size_t>(std::distance(
        alignment.insertions.begin(), alignment.insertions.end()));
    size_t num_deletions = static_cast<size_t>(std::distance(
        alignment.deletions.begin(), alignment.deletions.end()));

    // Build simplified CIGAR
    // Format: {align_length}M[{insertions}I][{deletions}D]
    cigar << alignment.align_length << "M";

    if (num_insertions > 0) {
        cigar << num_insertions << "I";
    }

    if (num_deletions > 0) {
        cigar << num_deletions << "D";
    }

    return cigar.str();
}

std::vector<std::string> get_airr_columns(bool include_alignment_details)
{
    std::vector<std::string> columns = {"sequence_id", "sequence"};

    if (include_alignment_details) {
        // V gene
        columns.insert(columns.end(), {
            "v_call", "v_score", "v_sequence_start", "v_alignment_length", "v_cigar"
        });

        // D gene
        columns.insert(columns.end(), {
            "d_call", "d_score", "d_sequence_start", "d_alignment_length", "d_cigar"
        });

        // J gene
        columns.insert(columns.end(), {
            "j_call", "j_score", "j_sequence_start", "j_alignment_length", "j_cigar"
        });
    }

    return columns;
}

void write(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences,
    Delimiter delimiter)
{
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::write: Cannot open file for writing: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Write header
    write_header(file, delim);

    // Write data rows
    for (const auto& seq : sequences) {
        write_sequence_row(file, seq, delim);
    }
}

void write_tsv(const std::string& filepath, const std::vector<SequenceData>& sequences)
{
    write(filepath, sequences, Delimiter::TAB);
}

void write_csv(const std::string& filepath, const std::vector<SequenceData>& sequences)
{
    write(filepath, sequences, Delimiter::COMMA);
}

void write_legacy_tsv(
    const std::string& filepath,
    const std::vector<std::tuple<int, std::string,
                                 std::unordered_map<Gene_class, std::vector<Alignment_data>>>>& sequences)
{
    // Convert legacy tuples to SequenceData
    std::vector<SequenceData> seq_data;
    seq_data.reserve(sequences.size());

    for (const auto& [index, sequence, alignments] : sequences) {
        seq_data.emplace_back(index, sequence, alignments);
    }

    write_tsv(filepath, seq_data);
}

void write_batch(
    const std::string& filepath,
    const sparrow::record_batch& batch,
    Delimiter delimiter)
{
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::write_batch: Cannot open file for writing: " + filepath);
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
    for (size_t row = 0; row < batch.nb_rows(); ++row) {
        for (size_t col = 0; col < column_names.size(); ++col) {
            if (col > 0) file << delim;

            // Get value using helper function
            file << get_string_value(batch, column_names[col], row);
        }
        file << "\n";
    }
}

} // namespace igor::airr::rearrangement

