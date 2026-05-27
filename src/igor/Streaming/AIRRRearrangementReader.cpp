/*
 * AIRRRearrangementReader.cpp
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

#include <igor/Streaming/AIRRRearrangementReader.h>

#include <sparrow/array.hpp>
#include <sparrow/builder.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <optional>

namespace igor::airr::rearrangement {

//==============================================================================
// Internal helpers (anonymous namespace)
//==============================================================================

namespace {

std::vector<std::string> parse_line(const std::string& line, char delim)
{
    std::vector<std::string> fields;
    std::string field;
    std::istringstream stream(line);

    while (std::getline(stream, field, delim)) {
        // Trim whitespace
        size_t start = field.find_first_not_of(" \t\r\n");
        size_t end = field.find_last_not_of(" \t\r\n");
        if (start != std::string::npos && end != std::string::npos) {
            fields.push_back(field.substr(start, end - start + 1));
        } else {
            fields.push_back("");
        }
    }

    // Handle trailing delimiter (empty last field)
    if (!line.empty() && line.back() == delim) {
        fields.push_back("");
    }

    return fields;
}

std::unordered_map<std::string, size_t> build_column_index(
    const std::vector<std::string>& header)
{
    std::unordered_map<std::string, size_t> index;
    for (size_t i = 0; i < header.size(); ++i) {
        index[header[i]] = i;
    }
    return index;
}

std::string get_field(
    const std::vector<std::string>& row,
    const std::unordered_map<std::string, size_t>& col_idx,
    const std::string& name)
{
    auto it = col_idx.find(name);
    if (it == col_idx.end() || it->second >= row.size()) {
        return "";
    }
    return row[it->second];
}

int get_int_field(
    const std::vector<std::string>& row,
    const std::unordered_map<std::string, size_t>& col_idx,
    const std::string& name,
    int default_value = 0)
{
    std::string value = get_field(row, col_idx, name);
    if (value.empty()) {
        return default_value;
    }
    try {
        return std::stoi(value);
    } catch (...) {
        return default_value;
    }
}

double get_double_field(
    const std::vector<std::string>& row,
    const std::unordered_map<std::string, size_t>& col_idx,
    const std::string& name,
    double default_value = 0.0)
{
    std::string value = get_field(row, col_idx, name);
    if (value.empty()) {
        return default_value;
    }
    try {
        return std::stod(value);
    } catch (...) {
        return default_value;
    }
}

std::optional<Alignment_data> extract_alignment(
    const std::vector<std::string>& row,
    const std::unordered_map<std::string, size_t>& col_idx,
    Gene_class gene_class)
{
    // Determine prefix based on gene class
    std::string prefix;
    switch (gene_class) {
        case V_gene: prefix = "v_"; break;
        case D_gene: prefix = "d_"; break;
        case J_gene: prefix = "j_"; break;
        default: return std::nullopt;
    }

    // Get gene call (e.g., "IGHV1-2*01")
    std::string gene_name = get_field(row, col_idx, prefix + "call");
    if (gene_name.empty()) {
        return std::nullopt; // No alignment for this gene class
    }

    // Get alignment parameters
    int offset = get_int_field(row, col_idx, prefix + "sequence_start", 0);
    double score = get_double_field(row, col_idx, prefix + "score", 0.0);
    size_t align_length = static_cast<size_t>(
        get_int_field(row, col_idx, prefix + "alignment_length", 0));

    // AIRR uses 1-based positions, convert to 0-based
    if (offset > 0) {
        offset -= 1;
    }

    // Get germline start/end for calculating offsets
    int germline_start = get_int_field(row, col_idx, prefix + "germline_start", 1);
    int germline_end = get_int_field(row, col_idx, prefix + "germline_end", 0);

    // Calculate five_p_offset and three_p_offset
    size_t five_p_offset = germline_start > 0 ? static_cast<size_t>(germline_start - 1) : 0;
    size_t three_p_offset = 0;
    if (germline_end > 0 && align_length > 0) {
        three_p_offset = germline_end > 0 ? static_cast<size_t>(germline_end) : 0;
    }

    // Create alignment data
    // Note: insertions, deletions, and mismatches would require CIGAR parsing
    std::forward_list<int> insertions;
    std::forward_list<int> deletions;
    std::vector<int> mismatches;

    return Alignment_data(
        gene_name,
        offset,
        five_p_offset,
        three_p_offset,
        align_length,
        insertions,
        deletions,
        mismatches,
        score
    );
}

} // anonymous namespace

//==============================================================================
// Public API
//==============================================================================

// Use shared utility functions from parent namespace
using airr::detect_delimiter;
using airr::delimiter_char;

FileInfo get_file_info(const std::string& filepath, Delimiter delimiter)
{
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::get_file_info: Cannot open file: " + filepath);
    }

    FileInfo info;
    info.filepath = filepath;
    info.delimiter = delimiter;

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::get_file_info: File is empty: " + filepath);
    }

    info.column_names = parse_line(header_line, delim);

    // Check for required/common columns
    for (const auto& col : info.column_names) {
        if (col == "sequence_id") info.has_sequence_id = true;
        else if (col == "sequence") info.has_sequence = true;
        else if (col == "v_call") info.has_v_call = true;
        else if (col == "d_call") info.has_d_call = true;
        else if (col == "j_call") info.has_j_call = true;
    }

    // Count data rows
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            ++info.num_rows;
        }
    }

    return info;
}

bool validate_schema(const std::string& filepath, Delimiter delimiter)
{
    auto info = get_file_info(filepath, delimiter);

    // AIRR minimal requirements: sequence_id and sequence
    return info.has_sequence_id && info.has_sequence;
}

sparrow::record_batch read_batch(const std::string& filepath, Delimiter delimiter)
{
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::read_batch: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::read_batch: File is empty: " + filepath);
    }

    auto column_names = parse_line(header_line, delim);
    size_t num_cols = column_names.size();

    // Prepare column data storage (all as strings initially)
    std::vector<std::vector<std::string>> column_data(num_cols);

    // Read data rows
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = parse_line(line, delim);

        // Pad with empty strings if row has fewer fields
        while (fields.size() < num_cols) {
            fields.push_back("");
        }

        for (size_t i = 0; i < num_cols; ++i) {
            column_data[i].push_back(i < fields.size() ? fields[i] : "");
        }
    }

    // Build record batch
    std::vector<std::string> names = column_names;
    std::vector<sparrow::array> arrays;

    for (size_t i = 0; i < num_cols; ++i) {
        arrays.emplace_back(sparrow::build(column_data[i]));
    }

    return sparrow::record_batch(std::move(names), std::move(arrays));
}

std::vector<SequenceData> read_sequences(const std::string& filepath, Delimiter delimiter)
{
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::read_sequences: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::read_sequences: File is empty: " + filepath);
    }

    auto header = parse_line(header_line, delim);
    auto col_idx = build_column_index(header);

    std::vector<SequenceData> sequences;
    std::string line;
    int row_num = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = parse_line(line, delim);

        // Get sequence_id (use row number as fallback)
        int seq_id = row_num;
        std::string seq_id_str = get_field(fields, col_idx, "sequence_id");
        if (!seq_id_str.empty()) {
            try {
                seq_id = std::stoi(seq_id_str);
            } catch (...) {
                // Keep row_num as fallback if not a valid integer
            }
        }

        // Get sequence
        std::string sequence = get_field(fields, col_idx, "sequence");

        // Extract alignments for V, D, J genes
        std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

        auto v_align = extract_alignment(fields, col_idx, V_gene);
        if (v_align) {
            alignments[V_gene].push_back(*v_align);
        }

        auto d_align = extract_alignment(fields, col_idx, D_gene);
        if (d_align) {
            alignments[D_gene].push_back(*d_align);
        }

        auto j_align = extract_alignment(fields, col_idx, J_gene);
        if (j_align) {
            alignments[J_gene].push_back(*j_align);
        }

        sequences.emplace_back(seq_id, sequence, alignments);
        ++row_num;
    }

    return sequences;
}

std::vector<std::tuple<int, std::string,
                       std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_legacy(const std::string& filepath, Delimiter delimiter)
{
    auto sequences = read_sequences(filepath, delimiter);

    std::vector<std::tuple<int, std::string,
                           std::unordered_map<Gene_class, std::vector<Alignment_data>>>> result;
    result.reserve(sequences.size());

    for (const auto& seq : sequences) {
        result.emplace_back(seq.index, seq.sequence, seq.alignments);
    }

    return result;
}

sparrow::record_batch read_columns(
    const std::string& filepath,
    const std::vector<std::string>& columns,
    Delimiter delimiter)
{
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::read_columns: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::read_columns: File is empty: " + filepath);
    }

    auto all_columns = parse_line(header_line, delim);
    auto col_idx = build_column_index(all_columns);

    // Validate requested columns exist
    std::vector<size_t> column_indices;
    for (const auto& col : columns) {
        auto it = col_idx.find(col);
        if (it == col_idx.end()) {
            throw std::runtime_error("airr::read_columns: Column not found: " + col);
        }
        column_indices.push_back(it->second);
    }

    // Prepare storage for selected columns
    std::vector<std::vector<std::string>> column_data(columns.size());

    // Read data rows
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = parse_line(line, delim);

        for (size_t i = 0; i < column_indices.size(); ++i) {
            size_t idx = column_indices[i];
            column_data[i].push_back(idx < fields.size() ? fields[idx] : "");
        }
    }

    // Build record batch with selected columns
    std::vector<std::string> names = columns;
    std::vector<sparrow::array> arrays;

    for (size_t i = 0; i < columns.size(); ++i) {
        arrays.emplace_back(sparrow::build(column_data[i]));
    }

    return sparrow::record_batch(std::move(names), std::move(arrays));
}

} // namespace igor::airr::rearrangement

