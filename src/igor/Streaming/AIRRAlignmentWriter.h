/*
 * AIRRAlignmentWriter.h
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
 * @file AIRRAlignmentWriter.h
 * @brief Writer for AIRR-Community Alignment schema TSV/CSV files
 *
 * This module provides functionality to write immune receptor alignment data
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
 * **Field Mapping from IGoR:**
 * ```
 * IGoR Field                      → AIRR Field
 * ----------                      → -----------
 * SequenceData::index             → sequence_id
 * Gene class (V/D/J)              → segment
 * Alignment_data::gene_name       → call
 * Alignment_data::score           → score
 * Alignment_data::offset          → sequence_start (0-based → 1-based)
 * offset + align_length           → sequence_end
 * Alignment_data::five_p_offset   → germline_start (0-based → 1-based)
 * insertions/deletions            → cigar
 * (alignment order in vector)     → rank
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

#include <exception>
#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>

namespace igor::airr::alignment {

// Use shared types from parent namespace
using airr::Delimiter;

/**
 * @brief Write sequences to an AIRR Alignment TSV file
 *
 * Writes alignment data in AIRR Alignment schema format (one row per alignment).
 * Each V/D/J alignment becomes a separate row in the file.
 *
 * @param filepath Path to output file
 * @param sequences Vector of SequenceData to write
 * @param delimiter Delimiter to use (default: TAB)
 * @throws std::runtime_error if file cannot be opened
 *
 * Example usage:
 * @code
 *   using namespace igor::airr::alignment;
 *
 *   std::vector<SequenceData> seqs = ...;
 *   write_sequences("alignments.tsv", seqs);
 * @endcode
 */
STREAMING_EXPORT
void write_sequences(
    const std::string& filepath,
    const std::vector<SequenceData>& sequences,
    Delimiter delimiter = Delimiter::TAB);

/**
 * @brief Write record_batch to an AIRR Alignment file
 *
 * Writes a Sparrow record_batch directly to AIRR Alignment format.
 * The batch must contain string columns with AIRR Alignment field names.
 *
 * @param filepath Path to output file
 * @param batch Record batch to write
 * @param delimiter Delimiter to use (default: TAB)
 * @throws std::runtime_error if file cannot be opened
 */
STREAMING_EXPORT
void write_batch(
    const std::string& filepath,
    const sparrow::record_batch& batch,
    Delimiter delimiter = Delimiter::TAB);

/**
 * @brief Generate CIGAR string from Alignment_data
 *
 * Creates a simplified CIGAR string from alignment data.
 * Format: operations in order (e.g., "50M2I3D" = 50 matches, 2 insertions, 3 deletions)
 *
 * **Limitation:** This generates a simplified CIGAR string that only includes
 * total counts of insertions and deletions, not their exact positions within
 * the alignment. For full CIGAR generation, per-base alignment information
 * would be required.
 *
 * @param align Alignment data
 * @return CIGAR string (e.g., "50M2I3D")
 */
STREAMING_EXPORT
std::string make_cigar(const Alignment_data& align);

} // namespace igor::airr::alignment
