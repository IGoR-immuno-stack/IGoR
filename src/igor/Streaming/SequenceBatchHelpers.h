/*
 * SequenceBatchHelpers.h
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 *  analyze and model immune receptors generation, selection, mutation and all
 *  other processes.
 *   Copyright (C) 2025  IGoR Development Team
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

#include <igor/Streaming/Export.h>

#include <igor/Core/Utils.h>
#include <igor/Core/Aligner.h>

#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace igor {

/**
 * @brief Represents sequence data in IGoR's internal format
 *
 * This is a simple struct that packages together the essential components
 * of a sequence for processing in IGoR.
 */
struct STREAMING_EXPORT SequenceData
{
    int index;
    std::string sequence;
    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

    SequenceData() : index(-1), sequence(""), alignments() { }

    SequenceData(int idx, const std::string &seq) : index(idx), sequence(seq), alignments() { }

    SequenceData(int idx, const std::string &seq,
                 const std::unordered_map<Gene_class, std::vector<Alignment_data>> &aligns)
        : index(idx), sequence(seq), alignments(aligns)
    {
    }
};

/**
 * @brief Convert a Sparrow record_batch row to IGoR SequenceData
 *
 * Extracts sequence information from a Sparrow record_batch row and converts
 * it to IGoR's internal SequenceData format. Handles missing values gracefully.
 *
 * @param batch The Sparrow record_batch containing sequence data
 * @param row_index Zero-based row index (must be < batch.nb_rows())
 * @return SequenceData object with parsed sequence and alignments
 * @throws std::out_of_range if row_index is invalid
 * @throws std::runtime_error if required columns are missing
 *
 * Expected columns in batch:
 * - "sequence_id" or "seq_index" (optional, uses row_index as fallback)
 * - "sequence" (required)
 * - Alignment columns (optional, parsed if present)
 */
STREAMING_EXPORT
SequenceData row_to_sequence_data(const sparrow::record_batch &batch, size_t row_index);

/**
 * @brief Convert legacy vector format to Sparrow record_batch
 *
 * Provides backward compatibility by converting IGoR's legacy sequence format
 * (vector of tuples) to a Sparrow record_batch. This allows existing code and
 * tests to work with the new streaming architecture without modification.
 *
 * @param sequences Vector of tuples containing (index, sequence, alignments)
 * @return Sparrow record_batch with columns for id, sequence, and alignments
 *
 * The resulting batch will have columns:
 * - "sequence_id" (int32)
 * - "sequence" (string)
 * - Alignment-related columns (if alignments present)
 */
STREAMING_EXPORT
sparrow::record_batch vector_to_batch(
        const std::vector<std::tuple<int, std::string,
                                     std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
                &sequences);

/**
 * @brief Parse alignment data from batch columns for a specific row
 *
 * Internal helper function to extract alignment information from columnar format.
 * Looks for columns matching alignment patterns (e.g., "v_gene_0_name", "v_gene_0_offset").
 *
 * @param batch The record_batch containing alignment columns
 * @param row_index The row to extract alignments from
 * @return Map of gene class to vector of alignment data
 */
std::unordered_map<Gene_class, std::vector<Alignment_data>>
parse_alignments_from_columns(const sparrow::record_batch &batch, size_t row_index);

/**
 * @brief Check if a column exists in the batch
 *
 * @param batch The record_batch to check
 * @param column_name Name of the column
 * @return true if column exists, false otherwise
 */
bool has_column(const sparrow::record_batch &batch, const std::string &column_name);

/**
 * @brief Get string value from batch, handling nulls
 *
 * @param batch The record_batch
 * @param column_name Name of the column
 * @param row_index Row index
 * @param default_value Value to return if null or missing
 * @return String value or default
 */
std::string get_string_value(const sparrow::record_batch &batch, const std::string &column_name,
                             size_t row_index, const std::string &default_value = "");

/**
 * @brief Get integer value from batch, handling nulls
 *
 * @param batch The record_batch
 * @param column_name Name of the column
 * @param row_index Row index
 * @param default_value Value to return if null or missing
 * @return Integer value or default
 */
int get_int_value(const sparrow::record_batch &batch, const std::string &column_name,
                  size_t row_index, int default_value = -1);

/**
 * @brief Get double value from batch, handling nulls
 *
 * @param batch The record_batch
 * @param column_name Name of the column
 * @param row_index Row index
 * @param default_value Value to return if null or missing
 * @return Double value or default
 */
double get_double_value(const sparrow::record_batch &batch, const std::string &column_name,
                        size_t row_index, double default_value = 0.0);

} // namespace igor