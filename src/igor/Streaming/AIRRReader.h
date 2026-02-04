/*
 * AIRRReader.h
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

/**
 * @file AIRRReader.h
 * @brief Reader for AIRR-Community standard TSV/CSV rearrangement files
 *
 * This module provides functionality to read immune receptor rearrangement data
 * in the AIRR (Adaptive Immune Receptor Repertoire) standard format.
 *
 * **AIRR Rearrangement Format (v1.4):**
 * The AIRR Community defines a standardized format for sharing immune repertoire
 * sequencing data. Key characteristics:
 * - Tab-delimited (TSV) or comma-delimited (CSV) text files
 * - Required header row with column names
 * - All fields are nullable (empty string = no data)
 * - UTF-8 or ASCII encoding
 *
 * **Field Mapping to IGoR:**
 * ```
 * AIRR Field          → IGoR Field
 * -----------         → ----------
 * sequence_id         → SequenceData::index (parsed as int if possible)
 * sequence            → SequenceData::sequence
 * v_call              → Alignment_data::gene_name (V_gene)
 * d_call              → Alignment_data::gene_name (D_gene)
 * j_call              → Alignment_data::gene_name (J_gene)
 * v_score             → Alignment_data::score (V_gene)
 * d_score             → Alignment_data::score (D_gene)
 * j_score             → Alignment_data::score (J_gene)
 * v_sequence_start    → Alignment_data::offset (V_gene)
 * d_sequence_start    → Alignment_data::offset (D_gene)
 * j_sequence_start    → Alignment_data::offset (J_gene)
 * v_alignment_length  → Alignment_data::align_length (V_gene)
 * d_alignment_length  → Alignment_data::align_length (D_gene)
 * j_alignment_length  → Alignment_data::align_length (J_gene)
 * ```
 *
 * @see https://docs.airr-community.org/en/stable/datarep/rearrangements.html
 */

#pragma once

#include <igor/Streaming/Export.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>
#include <unordered_map>

namespace igor::airr {

/**
 * @brief Delimiter type for AIRR files
 */
enum class Delimiter
{
    TAB,   ///< Tab-separated values (TSV) - AIRR standard
    COMMA, ///< Comma-separated values (CSV) - Alternative
    AUTO   ///< Auto-detect based on file content
};

/**
 * @brief Metadata about an AIRR file
 */
struct STREAMING_EXPORT FileInfo
{
    std::string filepath;                     ///< Path to the file
    size_t num_rows{0};                       ///< Number of data rows (excluding header)
    std::vector<std::string> column_names;    ///< Column names from header
    Delimiter delimiter{Delimiter::TAB};      ///< Detected delimiter
    bool has_sequence_id{false};              ///< Has sequence_id column
    bool has_sequence{false};                 ///< Has sequence column
    bool has_v_call{false};                   ///< Has v_call column
    bool has_d_call{false};                   ///< Has d_call column
    bool has_j_call{false};                   ///< Has j_call column
};

/**
 * @brief Get metadata about an AIRR file without reading all data
 *
 * Reads the header and counts rows to provide file information.
 * Auto-detects delimiter if not specified.
 *
 * @param filepath Path to the AIRR TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return FileInfo with file metadata
 * @throws std::runtime_error if file cannot be opened or parsed
 *
 * Example usage:
 * @code
 *   using namespace igor::airr;
 *
 *   auto info = get_file_info("data.tsv");
 *   std::cout << "Found " << info.num_rows << " sequences\n";
 * @endcode
 */
STREAMING_EXPORT
FileInfo get_file_info(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read AIRR file as Sparrow record_batch
 *
 * Converts AIRR columnar data to Arrow-compatible batch format.
 * This is useful for analytical processing.
 *
 * @param filepath Path to the AIRR TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return sparrow::record_batch with AIRR columns
 * @throws std::runtime_error if file cannot be opened or parsed
 */
STREAMING_EXPORT
sparrow::record_batch read_batch(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read AIRR file as vector of SequenceData
 *
 * Converts AIRR data to IGoR's internal SequenceData format,
 * mapping V/D/J calls to alignment structures.
 *
 * @param filepath Path to the AIRR TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return Vector of SequenceData with sequences and alignments
 * @throws std::runtime_error if file cannot be opened or parsed
 */
STREAMING_EXPORT
std::vector<SequenceData> read_sequences(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read AIRR file as legacy tuple format
 *
 * Provides backward compatibility with legacy code that expects
 * the tuple format (index, sequence, alignments).
 *
 * @param filepath Path to the AIRR TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return Vector of tuples in legacy format
 * @throws std::runtime_error if file cannot be opened or parsed
 */
STREAMING_EXPORT
std::vector<std::tuple<int, std::string,
                       std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_legacy(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read specific columns from AIRR file
 *
 * Reads only the specified columns for efficiency when not all
 * data is needed.
 *
 * @param filepath Path to the AIRR TSV/CSV file
 * @param columns List of column names to read
 * @param delimiter Delimiter type (default: AUTO)
 * @return sparrow::record_batch with only specified columns
 * @throws std::runtime_error if file cannot be opened or columns not found
 */
STREAMING_EXPORT
sparrow::record_batch read_columns(
    const std::string& filepath,
    const std::vector<std::string>& columns,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Auto-detect delimiter from file content
 *
 * Examines the first line to determine if file is TSV or CSV.
 *
 * @param filepath Path to the file
 * @return Detected delimiter (TAB or COMMA)
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
Delimiter detect_delimiter(const std::string& filepath);

/**
 * @brief Validate AIRR schema compliance
 *
 * Checks if the file has the required AIRR columns.
 * Required: sequence_id, sequence
 * Optional but common: v_call, d_call, j_call, productive, junction, etc.
 *
 * @param filepath Path to the AIRR file
 * @param delimiter Delimiter type (default: AUTO)
 * @return true if file has required columns, false otherwise
 */
STREAMING_EXPORT
bool validate_schema(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Get delimiter character
 *
 * @param delimiter Delimiter enum
 * @return Character to use as delimiter
 */
STREAMING_EXPORT
char delimiter_char(Delimiter delimiter);

} // namespace igor::airr
