/*
 * ParquetReader.h
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

namespace igor {

/**
 * @brief Reader for Parquet format files
 *
 * This class provides functionality to read IGoR sequence data from Parquet format
 * files, enabling fast loading of large datasets with efficient memory usage.
 *
 * Key features:
 * - Reads Parquet files written by ParquetWriter
 * - Returns data as sparrow::record_batch for efficient processing
 * - Can convert to legacy vector format for backward compatibility
 * - Supports reading metadata (row count, column info)
 *
 * Example usage:
 * @code
 * // Read entire file as record_batch
 * auto batch = ParquetReader::read_batch("sequences.parquet");
 *
 * // Read and convert to legacy format
 * auto sequences = ParquetReader::read_sequences("sequences.parquet");
 *
 * // Get metadata without reading data
 * auto info = ParquetReader::get_file_info("sequences.parquet");
 * std::cout << "Rows: " << info.num_rows << std::endl;
 * @endcode
 */
class STREAMING_EXPORT ParquetReader
{
public:
    /**
     * @brief Metadata about a Parquet file
     */
    struct FileInfo
    {
        int64_t num_rows;           ///< Number of rows in the file
        int num_columns;            ///< Number of columns in the file
        std::vector<std::string> column_names;  ///< Names of all columns
        std::string compression;    ///< Compression type used
    };

    /**
     * @brief Read entire Parquet file into a record_batch
     *
     * @param input_path Path to input Parquet file
     * @return sparrow::record_batch containing all data
     * @throws std::runtime_error if file cannot be read or is invalid
     */
    static sparrow::record_batch read_batch(const std::string &input_path);

    /**
     * @brief Read Parquet file and convert to legacy vector format
     *
     * Reads a Parquet file and converts it to the legacy IGoR sequence format.
     * This is useful for backward compatibility with existing code.
     *
     * @param input_path Path to input Parquet file
     * @return Vector of tuples (id, sequence, alignments)
     * @throws std::runtime_error if file cannot be read or conversion fails
     */
    static std::vector<std::tuple<int, std::string,
                                  std::unordered_map<Gene_class,
                                                     std::vector<Alignment_data>>>>
    read_sequences(const std::string &input_path);

    /**
     * @brief Get metadata about a Parquet file without reading all data
     *
     * @param input_path Path to input Parquet file
     * @return FileInfo struct with metadata
     * @throws std::runtime_error if file cannot be opened
     */
    static FileInfo get_file_info(const std::string &input_path);

    /**
     * @brief Read only specific columns from Parquet file
     *
     * @param input_path Path to input Parquet file
     * @param column_names Names of columns to read
     * @return sparrow::record_batch containing only requested columns
     * @throws std::runtime_error if file cannot be read or columns don't exist
     */
    static sparrow::record_batch read_columns(const std::string &input_path,
                                              const std::vector<std::string> &column_names);
};

} // namespace igor
