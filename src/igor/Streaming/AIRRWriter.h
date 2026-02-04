/*
 * AIRRWriter.h
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
 * @file AIRRWriter.h
 * @brief Writer for AIRR-Community standard TSV/CSV rearrangement files
 *
 * This module provides functionality to write immune receptor rearrangement data
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
 * **Field Mapping from IGoR:**
 * ```
 * IGoR Field                      → AIRR Field
 * ----------                      → -----------
 * SequenceData::index             → sequence_id
 * SequenceData::sequence          → sequence
 * Alignment_data::gene_name       → v_call, d_call, j_call
 * Alignment_data::score           → v_score, d_score, j_score
 * Alignment_data::offset          → v_sequence_start, d_sequence_start, j_sequence_start
 * Alignment_data::align_length    → v_alignment_length, d_alignment_length, j_alignment_length
 * insertions/deletions/mismatches → v_cigar, d_cigar, j_cigar (simplified)
 * ```
 *
 * @see https://docs.airr-community.org/en/stable/datarep/rearrangements.html
 */

#pragma once

#include <igor/Streaming/Export.h>
#include <igor/Streaming/SequenceBatchHelpers.h>
#include <igor/Streaming/AIRRReader.h>

#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>
#include <unordered_map>

namespace igor::airr {

/**
 * @brief Write sequences to an AIRR TSV file
 *
 * Writes sequence data in AIRR v1.4 TSV format with all available
 * alignment information.
 *
 * @param filepath Path to output file
 * @param sequences Vector of SequenceData to write
 * @throws std::runtime_error if file cannot be opened
 *
 * Example:
 * @code
 *   using namespace igor::airr;
 *   std::vector<SequenceData> seqs = ...;
 *   write_tsv("output.tsv", seqs);
 * @endcode
 */
STREAMING_EXPORT
void write_tsv(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences);

/**
 * @brief Write sequences to an AIRR CSV file
 *
 * Writes sequence data in AIRR format with comma delimiters.
 *
 * @param filepath Path to output file
 * @param sequences Vector of SequenceData to write
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
void write_csv(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences);

/**
 * @brief Write sequences to an AIRR file with specified delimiter
 *
 * Low-level write function that allows specifying the delimiter.
 *
 * @param filepath Path to output file
 * @param sequences Vector of SequenceData to write
 * @param delimiter Delimiter to use (TAB or COMMA)
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
void write(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences,
    Delimiter delimiter);

/**
 * @brief Write legacy tuple format to an AIRR TSV file
 *
 * Provides backward compatibility with legacy code that uses
 * the tuple format (index, sequence, alignments).
 *
 * @param filepath Path to output file
 * @param sequences Vector of tuples in legacy format
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
void write_legacy_tsv(
    const std::string& filepath,
    const std::vector<std::tuple<int, std::string,
                                 std::unordered_map<Gene_class, std::vector<Alignment_data>>>>& sequences);

/**
 * @brief Write a Sparrow record_batch to an AIRR TSV file
 *
 * Writes batch data directly to AIRR format.
 *
 * **Note:** This function expects all columns in the batch to be
 * string-typed. For numeric or list columns, convert to strings first
 * or use write() with SequenceData instead.
 *
 * @param filepath Path to output file
 * @param batch Record batch with sequence data (string columns only)
 * @param delimiter Delimiter to use (default: TAB)
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
void write_batch(
    const std::string& filepath,
    const sparrow::record_batch& batch,
    Delimiter delimiter = Delimiter::TAB);

/**
 * @brief Generate a CIGAR string from alignment data
 *
 * Creates a simplified CIGAR string representation of the alignment.
 * Format: {align_length}M[{num_insertions}I][{num_deletions}D]
 *
 * **Limitation:** This creates a simplified CIGAR that summarizes
 * the total counts (e.g., "50M2I3D") but does NOT preserve the
 * exact positions of each insertion/deletion/mismatch within the
 * alignment. Full CIGAR reconstruction would require access to
 * the aligned sequences.
 *
 * @param alignment Alignment data
 * @return CIGAR string (e.g., "50M2I3D"), or empty string if align_length is 0
 */
STREAMING_EXPORT
std::string make_cigar(const Alignment_data& alignment);

/**
 * @brief Get AIRR column headers
 *
 * Returns the list of column names for AIRR output.
 *
 * @param include_alignment_details If true, includes detailed alignment columns
 * @return Vector of column names
 */
STREAMING_EXPORT
std::vector<std::string> get_airr_columns(bool include_alignment_details = true);

} // namespace igor::airr
