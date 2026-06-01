/*
 * ParquetWriter.h
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

#include <exception>
#include <sparrow/record_batch.hpp>

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>

namespace igor {

/**
 * @brief Compression types supported for Parquet files
 */
enum class STREAMING_EXPORT CompressionType
{
    NONE,    ///< No compression
    SNAPPY,  ///< Snappy compression (default - fast, ~3-4x)
    GZIP,    ///< GZIP compression (slower, better ratio)
    ZSTD,    ///< Zstandard compression (good balance, ~5-6x)
    LZ4      ///< LZ4 compression (very fast, moderate ratio)
};

/**
 * @brief Writer for Parquet format files
 *
 * This class provides functionality to write IGoR sequence data to Parquet format,
 * which offers significant compression and faster loading compared to CSV.
 *
 * Key features:
 * - Multiple compression algorithms (Snappy, GZIP, ZSTD, LZ4)
 * - Preserves all sequence and alignment data
 * - Compatible with Apache Arrow ecosystem
 *
 * Example usage:
 * @code
 * // Convert legacy vector format to Parquet
 * std::vector<std::tuple<int, std::string, AlignmentMap>> sequences = ...;
 * ParquetWriter::write_sequences("output.parquet", sequences);
 *
 * // Write from record_batch with custom compression
 * sparrow::record_batch batch = ...;
 * ParquetWriter::write_batch("output.parquet", batch, CompressionType::ZSTD);
 * @endcode
 */
class STREAMING_EXPORT ParquetWriter
{
public:
    /**
     * @brief Write sequences from legacy vector format to Parquet file
     *
     * @param output_path Path to output Parquet file
     * @param sequences Vector of tuples (index, sequence, alignments)
     * @param compression Compression algorithm to use (default: SNAPPY)
     * @throws std::runtime_error if file cannot be written
     */
    static void write_sequences(
            const std::string &output_path,
            const std::vector<std::tuple<int, std::string,
                                         std::unordered_map<Gene_class,
                                                            std::vector<Alignment_data>>>>
                    &sequences,
            CompressionType compression = CompressionType::SNAPPY);

    /**
     * @brief Write record_batch to Parquet file
     *
     * @param output_path Path to output Parquet file
     * @param batch Sparrow record_batch to write
     * @param compression Compression algorithm to use (default: SNAPPY)
     * @throws std::runtime_error if file cannot be written
     */
    static void write_batch(const std::string &output_path,
                           const sparrow::record_batch &batch,
                           CompressionType compression = CompressionType::SNAPPY);

    /**
     * @brief Get compression type name as string
     *
     * @param compression Compression type
     * @return String name of compression (e.g., "SNAPPY", "ZSTD")
     */
    static std::string compression_name(CompressionType compression);

private:
    /**
     * @brief Convert IGoR compression type to Arrow compression
     */
    static int get_arrow_compression(CompressionType compression);
};

} // namespace igor
