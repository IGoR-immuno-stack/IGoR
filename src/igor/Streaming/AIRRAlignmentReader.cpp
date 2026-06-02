/*
 * AIRRAlignmentReader.cpp
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

#include <igor/Streaming/AIRRAlignmentReader.h>

#include <sparrow/array.hpp>
#include <sparrow/builder.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <optional>
#include <cctype>

namespace igor::airr::alignment {

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

std::optional<std::string> get_field(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& col_index,
    const std::string& col_name)
{
    auto it = col_index.find(col_name);
    if (it == col_index.end() || it->second >= fields.size()) {
        return std::nullopt;
    }
    const auto& value = fields[it->second];
    return value.empty() ? std::nullopt : std::optional<std::string>(value);
}

template<typename T>
std::optional<T> parse_numeric(const std::optional<std::string>& str_val)
{
    if (!str_val) {
        return std::nullopt;
    }
    try {
        if constexpr (std::is_same_v<T, int>) {
            return std::stoi(*str_val);
        } else if constexpr (std::is_same_v<T, size_t>) {
            return std::stoull(*str_val);
        } else if constexpr (std::is_same_v<T, double>) {
            return std::stod(*str_val);
        }
    } catch (...) {
        return std::nullopt;
    }
    return std::nullopt;
}

// Map segment string to Gene_class
std::optional<Gene_class> parse_segment(const std::string& segment)
{
    if (segment.empty()) {
        return std::nullopt;
    }

    char first = std::toupper(segment[0]);
    switch (first) {
        case 'V': return Gene_class::V_gene;
        case 'D': return Gene_class::D_gene;
        case 'J': return Gene_class::J_gene;
        case 'C': return std::nullopt; // C genes not supported in IGoR
        default: return std::nullopt;
    }
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
    // Auto-detect delimiter if needed
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    FileInfo info;
    info.filepath = filepath;
    info.delimiter = delimiter;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::alignment::get_file_info: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::alignment::get_file_info: File is empty: " + filepath);
    }

    info.column_names = parse_line(header_line, delim);

    // Check for key columns
    for (const auto& col : info.column_names) {
        if (col == "sequence_id") info.has_sequence_id = true;
        else if (col == "segment") info.has_segment = true;
        else if (col == "call") info.has_call = true;
        else if (col == "cigar") info.has_cigar = true;
        else if (col == "rank") info.has_rank = true;
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

    // Required fields for AIRR Alignment schema
    return info.has_sequence_id && info.has_segment && info.has_call;
}

sparrow::record_batch read_batch(const std::string& filepath, Delimiter delimiter)
{
    // Auto-detect delimiter if needed
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::alignment::read_batch: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::alignment::read_batch: File is empty: " + filepath);
    }

    auto header = parse_line(header_line, delim);
    auto col_index = build_column_index(header);

    // Collect all rows first
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = parse_line(line, delim);

        // Ensure we have enough fields
        if (fields.size() < header.size()) {
            fields.resize(header.size(), "");
        }
        rows.push_back(std::move(fields));
    }

    // Build arrays for each column
    std::vector<sparrow::array> arrays;
    arrays.reserve(header.size());
    for (size_t col = 0; col < header.size(); ++col) {
        std::vector<std::string> column_data;
        column_data.reserve(rows.size());
        for (const auto& row : rows) {
            column_data.push_back(row[col]);
        }

        sparrow::array arr(sparrow::build(std::move(column_data)));
        arr.set_name(header[col]);
        arrays.push_back(std::move(arr));
    }

    return sparrow::record_batch(std::move(arrays));
}

std::vector<SequenceData> read_sequences(const std::string& filepath, Delimiter delimiter)
{
    // Auto-detect delimiter if needed
    if (delimiter == Delimiter::AUTO) {
        delimiter = detect_delimiter(filepath);
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::alignment::read_sequences: Cannot open file: " + filepath);
    }

    char delim = delimiter_char(delimiter);

    // Read header
    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("airr::alignment::read_sequences: File is empty: " + filepath);
    }

    auto header = parse_line(header_line, delim);
    auto col_index = build_column_index(header);

    // Validate required columns
    if (col_index.find("sequence_id") == col_index.end() ||
        col_index.find("segment") == col_index.end() ||
        col_index.find("call") == col_index.end()) {
        throw std::runtime_error(
            "airr::alignment::read_sequences: Missing required columns (sequence_id, segment, call)");
    }

    // Group alignments by sequence_id
    std::unordered_map<std::string, std::vector<std::tuple<Gene_class, Alignment_data>>>
        alignment_map;

    // Read all alignment rows
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto fields = parse_line(line, delim);
        if (fields.size() < header.size()) {
            continue; // Skip malformed rows
        }

        // Extract required fields
        auto sequence_id = get_field(fields, col_index, "sequence_id");
        auto segment = get_field(fields, col_index, "segment");
        auto call = get_field(fields, col_index, "call");

        if (!sequence_id || !segment || !call) {
            continue; // Skip rows with missing required fields
        }

        // Parse segment to gene class
        auto gene_class = parse_segment(*segment);
        if (!gene_class) {
            continue; // Skip unsupported segments (C genes)
        }

        // Parse offset (required for constructor)
        int offset = 0;
        if (auto seq_start = get_field(fields, col_index, "sequence_start")) {
            // Convert 1-based to 0-based
            auto offset_val = parse_numeric<int>(*seq_start);
            if (offset_val && *offset_val > 0) {
                offset = *offset_val - 1;
            }
        }

        // Create Alignment_data using simple constructor
        Alignment_data align(*call, offset);

        // Parse optional fields
        if (auto score = get_field(fields, col_index, "score")) {
            align.score = parse_numeric<double>(score).value_or(0.0);
        }

        if (auto seq_end = get_field(fields, col_index, "sequence_end");
            seq_end && offset >= 0) {
            auto end = parse_numeric<int>(*seq_end);
            if (end && *end > offset) {
                align.align_length = static_cast<size_t>(*end - offset);
            }
        }

        if (auto germ_start = get_field(fields, col_index, "germline_start")) {
            auto offset_val = parse_numeric<size_t>(*germ_start);
            if (offset_val && *offset_val > 0) {
                align.five_p_offset = *offset_val - 1;
            }
        }

        // TODO: Parse CIGAR string to extract insertions/deletions/mismatches
        // For now, leave them empty

        // Add to alignment map
        alignment_map[*sequence_id].emplace_back(*gene_class, align);
    }

    // Convert to SequenceData vector
    std::vector<SequenceData> sequences;
    sequences.reserve(alignment_map.size());

    for (auto& [seq_id, align_list] : alignment_map) {
        SequenceData seq;

        // Try to parse sequence_id as integer index
        try {
            seq.index = std::stoi(seq_id);
        } catch (...) {
            // If not a number, use hash
            seq.index = std::hash<std::string>{}(seq_id);
        }

        // Note: sequence field is not in alignment schema, would need to be provided separately
        seq.sequence = "";

        // Group alignments by gene class
        for (auto& [gene_class, align] : align_list) {
            seq.alignments[gene_class].push_back(std::move(align));
        }

        sequences.push_back(std::move(seq));
    }

    return sequences;
}

bool parse_cigar(
    const std::string& cigar,
    std::forward_list<int>& insertions,
    std::forward_list<int>& deletions,
    size_t& align_length)
{
    if (cigar.empty()) {
        return false;
    }

    insertions.clear();
    deletions.clear();
    align_length = 0;

    size_t i = 0;
    while (i < cigar.size()) {
        // Parse number
        size_t num_start = i;
        while (i < cigar.size() && std::isdigit(cigar[i])) {
            ++i;
        }

        if (i == num_start || i >= cigar.size()) {
            return false; // No number or no operation
        }

        int count = std::stoi(cigar.substr(num_start, i - num_start));
        char op = cigar[i++];

        switch (op) {
            case 'M': // Match/mismatch
                align_length += count;
                break;
            case 'I': // Insertion
                // Simplified: just record count
                for (int j = 0; j < count; ++j) {
                    insertions.push_front(static_cast<int>(align_length));
                }
                break;
            case 'D': // Deletion
                // Simplified: just record count
                for (int j = 0; j < count; ++j) {
                    deletions.push_front(static_cast<int>(align_length));
                }
                align_length += count;
                break;
            case 'N': // Skipped region
            case 'S': // Soft clipping
            case 'H': // Hard clipping
            case 'P': // Padding
            case '=': // Sequence match
            case 'X': // Sequence mismatch
                // Not yet supported in simplified parsing
                break;
            default:
                return false; // Unknown operation
        }
    }

    return true;
}

} // namespace igor::airr::alignment
