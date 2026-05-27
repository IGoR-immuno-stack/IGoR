/*
 * AIRRCommon.h
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
 * @file AIRRCommon.h
 * @brief Common types and utilities for AIRR format handling
 *
 * This header provides shared types used across AIRR Rearrangement and
 * Alignment schemas.
 */

#pragma once

#include <igor/Streaming/Export.h>

#include <string>
#include <vector>

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

    // Common fields present in both schemas
    bool has_sequence_id{false};              ///< Has sequence_id column
    bool has_sequence{false};                 ///< Has sequence column

    // Rearrangement-specific fields
    bool has_v_call{false};                   ///< Has v_call column
    bool has_d_call{false};                   ///< Has d_call column
    bool has_j_call{false};                   ///< Has j_call column

    // Alignment-specific fields
    bool has_segment{false};                  ///< Has segment column (alignment schema)
    bool has_call{false};                     ///< Has call column (alignment schema)
    bool has_cigar{false};                    ///< Has cigar column
    bool has_rank{false};                     ///< Has rank column (alignment schema)
};

/**
 * @brief Get delimiter character
 *
 * @param delimiter Delimiter enum
 * @return Character to use as delimiter
 */
STREAMING_EXPORT
char delimiter_char(Delimiter delimiter);

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

} // namespace igor::airr
