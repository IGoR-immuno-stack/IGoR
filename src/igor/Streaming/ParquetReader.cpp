/*
 * ParquetReader.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

// Include Arrow compatibility workarounds before any Arrow headers
#include <igor/Streaming/ArrowCompatibility.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/c/bridge.h>
#include <parquet/arrow/reader.h>

#include <igor/Streaming/ParquetReader.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <sstream>
#include <stdexcept>

namespace igor {

sparrow::record_batch ParquetReader::read_batch(const std::string &input_path)
{
    try {
        // Open the Parquet file
        auto input_result = arrow::io::ReadableFile::Open(input_path);
        if (!input_result.ok()) {
            throw std::runtime_error("Failed to open input file '" + input_path
                                     + "': " + input_result.status().ToString());
        }
        auto input_stream = input_result.ValueOrDie();

        // Create Parquet reader
        std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
        auto reader_status = parquet::arrow::OpenFile(input_stream,
                                                       arrow::default_memory_pool(),
                                                       &arrow_reader);
        if (!reader_status.ok()) {
            throw std::runtime_error("Failed to create Parquet reader: "
                                     + reader_status.ToString());
        }

        // Read entire file as Arrow Table
        std::shared_ptr<arrow::Table> table;
        auto read_status = arrow_reader->ReadTable(&table);
        if (!read_status.ok()) {
            throw std::runtime_error("Failed to read Parquet table: "
                                     + read_status.ToString());
        }

        // Convert Arrow Table to RecordBatch
        // Arrow Table can have multiple chunks, we'll combine them into one RecordBatch
        auto combine_result = table->CombineChunks(arrow::default_memory_pool());
        if (!combine_result.ok()) {
            throw std::runtime_error("Failed to combine table chunks: "
                                     + combine_result.status().ToString());
        }
        auto combined_table = combine_result.ValueOrDie();

        // Convert the first (and only) RecordBatch from the table
        arrow::TableBatchReader batch_reader(*combined_table);
        std::shared_ptr<arrow::RecordBatch> record_batch;
        auto batch_status = batch_reader.ReadNext(&record_batch);
        if (!batch_status.ok()) {
            throw std::runtime_error("Failed to read record batch: "
                                     + batch_status.ToString());
        }

        if (!record_batch) {
            // Empty file - return empty batch
            std::vector<std::string> empty_names;
            std::vector<sparrow::array> empty_arrays;
            return sparrow::record_batch(std::move(empty_names), std::move(empty_arrays));
        }

        // Export Arrow RecordBatch to C structures
        ArrowSchema arrow_schema;
        ArrowArray arrow_array;

        auto export_status = arrow::ExportRecordBatch(*record_batch, &arrow_array, &arrow_schema);
        if (!export_status.ok()) {
            throw std::runtime_error("Failed to export record batch: "
                                     + export_status.ToString());
        }

        // Import into Sparrow using the Arrow C Data Interface
        auto sparrow_batch = sparrow::record_batch(&arrow_array, &arrow_schema);

        return sparrow_batch;

    } catch (const std::exception &e) {
        std::ostringstream oss;
        oss << "Error reading Parquet file '" << input_path << "': " << e.what();
        throw std::runtime_error(oss.str());
    }
}

std::vector<std::tuple<int, std::string,
                       std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
ParquetReader::read_sequences(const std::string &input_path)
{
    // Read the entire file as a record_batch
    sparrow::record_batch batch = read_batch(input_path);

    // Convert record_batch to vector format
    std::vector<std::tuple<int, std::string,
                          std::unordered_map<Gene_class, std::vector<Alignment_data>>>> sequences;

    size_t num_rows = batch.nb_rows();
    sequences.reserve(num_rows);

    for (size_t i = 0; i < num_rows; ++i) {
        // Use helper function to convert each row
        SequenceData seq_data = row_to_sequence_data(batch, i);
        sequences.emplace_back(seq_data.index, seq_data.sequence, seq_data.alignments);
    }

    return sequences;
}

ParquetReader::FileInfo ParquetReader::get_file_info(const std::string &input_path)
{
    try {
        // Open the Parquet file
        auto input_result = arrow::io::ReadableFile::Open(input_path);
        if (!input_result.ok()) {
            throw std::runtime_error("Failed to open input file '" + input_path
                                     + "': " + input_result.status().ToString());
        }
        auto input_stream = input_result.ValueOrDie();

        // Create Parquet reader
        std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
        auto reader_status = parquet::arrow::OpenFile(input_stream,
                                                       arrow::default_memory_pool(),
                                                       &arrow_reader);
        if (!reader_status.ok()) {
            throw std::runtime_error("Failed to create Parquet reader: "
                                     + reader_status.ToString());
        }

        // Get file metadata
        auto parquet_reader = arrow_reader->parquet_reader();
        auto metadata = parquet_reader->metadata();

        FileInfo info;
        info.num_rows = metadata->num_rows();
        info.num_columns = metadata->num_columns();

        // Get schema to extract column names
        std::shared_ptr<arrow::Schema> schema;
        auto schema_status = arrow_reader->GetSchema(&schema);
        if (schema_status.ok()) {
            for (const auto &field : schema->fields()) {
                info.column_names.push_back(field->name());
            }
        }

        // Get compression info from first row group (if exists)
        if (metadata->num_row_groups() > 0) {
            auto row_group = metadata->RowGroup(0);
            if (row_group->num_columns() > 0) {
                auto column = row_group->ColumnChunk(0);
                auto compression = column->compression();

                switch (compression) {
                case parquet::Compression::UNCOMPRESSED:
                    info.compression = "NONE";
                    break;
                case parquet::Compression::SNAPPY:
                    info.compression = "SNAPPY";
                    break;
                case parquet::Compression::GZIP:
                    info.compression = "GZIP";
                    break;
                case parquet::Compression::ZSTD:
                    info.compression = "ZSTD";
                    break;
                case parquet::Compression::LZ4:
                    info.compression = "LZ4";
                    break;
                default:
                    info.compression = "UNKNOWN";
                }
            }
        }

        return info;

    } catch (const std::exception &e) {
        std::ostringstream oss;
        oss << "Error reading Parquet file info '" << input_path << "': " << e.what();
        throw std::runtime_error(oss.str());
    }
}

sparrow::record_batch ParquetReader::read_columns(const std::string &input_path,
                                                  const std::vector<std::string> &column_names)
{
    try {
        // Open the Parquet file
        auto input_result = arrow::io::ReadableFile::Open(input_path);
        if (!input_result.ok()) {
            throw std::runtime_error("Failed to open input file '" + input_path
                                     + "': " + input_result.status().ToString());
        }
        auto input_stream = input_result.ValueOrDie();

        // Create Parquet reader
        std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
        auto reader_status = parquet::arrow::OpenFile(input_stream,
                                                       arrow::default_memory_pool(),
                                                       &arrow_reader);
        if (!reader_status.ok()) {
            throw std::runtime_error("Failed to create Parquet reader: "
                                     + reader_status.ToString());
        }

        // Convert column names to indices
        std::shared_ptr<arrow::Schema> schema;
        auto schema_status = arrow_reader->GetSchema(&schema);
        if (!schema_status.ok()) {
            throw std::runtime_error("Failed to get schema: " + schema_status.ToString());
        }

        std::vector<int> column_indices;
        for (const auto &name : column_names) {
            int index = schema->GetFieldIndex(name);
            if (index == -1) {
                throw std::runtime_error("Column '" + name + "' not found in Parquet file");
            }
            column_indices.push_back(index);
        }

        // Read only the specified columns
        std::shared_ptr<arrow::Table> table;
        auto read_status = arrow_reader->ReadTable(column_indices, &table);
        if (!read_status.ok()) {
            throw std::runtime_error("Failed to read columns: " + read_status.ToString());
        }

        // Convert to RecordBatch
        auto combine_result = table->CombineChunks(arrow::default_memory_pool());
        if (!combine_result.ok()) {
            throw std::runtime_error("Failed to combine chunks: "
                                     + combine_result.status().ToString());
        }
        auto combined_table = combine_result.ValueOrDie();

        arrow::TableBatchReader batch_reader(*combined_table);
        std::shared_ptr<arrow::RecordBatch> record_batch;
        auto batch_status = batch_reader.ReadNext(&record_batch);
        if (!batch_status.ok()) {
            throw std::runtime_error("Failed to read batch: " + batch_status.ToString());
        }

        if (!record_batch) {
            std::vector<std::string> empty_names;
            std::vector<sparrow::array> empty_arrays;
            return sparrow::record_batch(std::move(empty_names), std::move(empty_arrays));
        }

        // Export to Sparrow
        ArrowSchema arrow_schema;
        ArrowArray arrow_array;

        auto export_status = arrow::ExportRecordBatch(*record_batch, &arrow_array, &arrow_schema);
        if (!export_status.ok()) {
            throw std::runtime_error("Failed to export batch: " + export_status.ToString());
        }

        return sparrow::record_batch(&arrow_array, &arrow_schema);

    } catch (const std::exception &e) {
        std::ostringstream oss;
        oss << "Error reading columns from Parquet file '" << input_path << "': " << e.what();
        throw std::runtime_error(oss.str());
    }
}

} // namespace igor
