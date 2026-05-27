/*
 * AIRRAlignmentReader.h
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
 * @file AIRRAlignmentReader.h
 * @brief Reader for AIRR-Community Alignment schema TSV/CSV files
 *
 * This module provides functionality to read immune receptor alignment data
 * in the AIRR Alignment schema format (experimental).
 *
 * **AIRR Alignment Schema vs Rearrangement Schema:**
 * - **Rearrangement Schema**: One row per sequence, V/D/J calls inlined as columns
 * - **Alignment Schema**: One row per alignment (multiple rows per sequence, one per V/D/J/C)
 *
 * The Alignment schema enables:
 * - Multiple candidate alignments per gene class
 * - Ranked alignments (primary, secondary, tertiary)
 * - More detailed alignment metadata
 *
 * **Key Fields:**
 * ```
 * AIRR Field          → IGoR Field
 * -----------         → ----------
 * sequence_id         → SequenceData::index
 * segment             → Gene class (V/D/J/C)
 * call                → Alignment_data::gene_name
 * score               → Alignment_data::score
 * sequence_start      → Alignment_data::offset (1-based → 0-based)
 * sequence_end        → offset + align_length
 * germline_start      → Alignment_data::five_p_offset (1-based → 0-based)
 * cigar               → Parse to insertions/deletions/mismatches
 * rank                → Alignment priority (1 = best)
 * ```
 *
 * **Status:** Experimental - AIRR Alignment schema is not yet finalized
 *
 * @see https://docs.airr-community.org/en/stable/datarep/alignments.html
 */

#pragma once

#include <igor/Streaming/Export.h>
#include <igor/Streaming/SequenceBatchHelpers.h>
#include <igor/Streaming/AIRRCommon.h>

#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>
#include <forward_list>

namespace igor::airr::alignment {

// Use shared types from parent namespace
using airr::Delimiter;
using airr::FileInfo;

/**
 * @brief Get metadata about an AIRR Alignment file
 *
 * @param filepath Path to the AIRR Alignment TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return FileInfo with file metadata
 * @throws std::runtime_error if file cannot be opened or parsed
 *
 * Example usage:
 * @code
 *   using namespace igor::airr::alignment;
 *
 *   // Read file info
 *   auto info = get_file_info("alignments.tsv");
 *   std::cout << "Found " << info.num_rows << " alignments\n";
 *
 *   // Read as SequenceData (groups by sequence_id)
 *   auto sequences = read_sequences("alignments.tsv");
 * @endcode
 */
STREAMING_EXPORT
FileInfo get_file_info(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Validate AIRR Alignment schema compliance
 *
 * Checks if the file has the required AIRR Alignment columns.
 * Required: sequence_id, segment, call
 *
 * @param filepath Path to the AIRR Alignment file
 * @param delimiter Delimiter type (default: AUTO)
 * @return true if file has required columns, false otherwise
 */
STREAMING_EXPORT
bool validate_schema(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read AIRR Alignment file as vector of SequenceData
 *
 * Groups alignments by sequence_id and converts to IGoR's internal format.
 * Multiple alignments per gene class are stored in the vector.
 *
 * @param filepath Path to the AIRR Alignment TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return Vector of SequenceData with grouped alignments
 * @throws std::runtime_error if file cannot be opened or parsed
 */
STREAMING_EXPORT
std::vector<SequenceData> read_sequences(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Read AIRR Alignment file as Sparrow record_batch
 *
 * Returns the raw alignment data without grouping.
 *
 * @param filepath Path to the AIRR Alignment TSV/CSV file
 * @param delimiter Delimiter type (default: AUTO)
 * @return sparrow::record_batch with alignment rows
 * @throws std::runtime_error if file cannot be opened or parsed
 */
STREAMING_EXPORT
sparrow::record_batch read_batch(
    const std::string& filepath,
    Delimiter delimiter = Delimiter::AUTO);

/**
 * @brief Parse CIGAR string to extract insertions, deletions, and mismatches
 *
 * Parses a CIGAR string (e.g., "50M2I3D") and extracts the operation counts.
 * This is the inverse of make_cigar() from AIRRAlignmentWriter.
 *
 * **CIGAR Semantics:**
 * - `M` (match/mismatch): Adds to alignment length (reference span)
 * - `I` (insertion): Does NOT add to alignment length (insertion in query)
 * - `D` (deletion): ADDS to alignment length (deletion from query = gap in query)
 *
 * **Note:** The returned `align_length` represents the reference span covered
 * by the alignment, which includes matches and deletions but not insertions.
 * This follows standard CIGAR conventions.
 *
 * **Limitation:** This only extracts counts from simplified CIGAR strings.
 * Full CIGAR parsing with exact positions is not yet implemented.
 *
 * @param cigar CIGAR string
 * @param insertions Output: list of insertion positions (currently just counts)
 * @param deletions Output: list of deletion positions (currently just counts)
 * @param align_length Output: reference span (M + D operations)
 * @return true if parsed successfully, false if CIGAR is invalid or empty
 */
STREAMING_EXPORT
bool parse_cigar(
    const std::string& cigar,
    std::forward_list<int>& insertions,
    std::forward_list<int>& deletions,
    size_t& align_length);

} // namespace igor::airr::alignment
