/*
 * ParquetWriter.cpp
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
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>

#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <sstream>
#include <stdexcept>

namespace igor {

int ParquetWriter::get_arrow_compression(CompressionType compression)
{
    switch (compression) {
    case CompressionType::NONE:
        return arrow::Compression::UNCOMPRESSED;
    case CompressionType::SNAPPY:
        return arrow::Compression::SNAPPY;
    case CompressionType::GZIP:
        return arrow::Compression::GZIP;
    case CompressionType::ZSTD:
        return arrow::Compression::ZSTD;
    case CompressionType::LZ4:
        return arrow::Compression::LZ4;
    default:
        return arrow::Compression::SNAPPY;
    }
}

std::string ParquetWriter::compression_name(CompressionType compression)
{
    switch (compression) {
    case CompressionType::NONE:
        return "NONE";
    case CompressionType::SNAPPY:
        return "SNAPPY";
    case CompressionType::GZIP:
        return "GZIP";
    case CompressionType::ZSTD:
        return "ZSTD";
    case CompressionType::LZ4:
        return "LZ4";
    default:
        return "UNKNOWN";
    }
}

void ParquetWriter::write_sequences(
        const std::string &output_path,
        const std::vector<std::tuple<int, std::string,
                                     std::unordered_map<Gene_class,
                                                        std::vector<Alignment_data>>>>
                &sequences,
        CompressionType compression)
{
    // Convert legacy format to record_batch
    sparrow::record_batch batch = vector_to_batch(sequences);

    // Write the batch
    write_batch(output_path, batch, compression);
}

void ParquetWriter::write_batch(const std::string &output_path,
                                const sparrow::record_batch &batch,
                                CompressionType compression)
{
    try {
        // Convert sparrow record_batch to Arrow RecordBatch
        // We need to make a copy since extract_arrow_structures takes ownership
        sparrow::record_batch batch_copy = batch;

        // Extract Arrow C structures from Sparrow
        auto [arrow_array, arrow_schema] = sparrow::extract_arrow_structures(std::move(batch_copy));

        // Import into Arrow C++ RecordBatch
        auto record_batch_result = arrow::ImportRecordBatch(&arrow_array, &arrow_schema);
        if (!record_batch_result.ok()) {
            throw std::runtime_error("Failed to import Arrow record batch: "
                                     + record_batch_result.status().ToString());
        }

        auto record_batch = record_batch_result.ValueOrDie();

        // Open output file
        auto output_result = arrow::io::FileOutputStream::Open(output_path);
        if (!output_result.ok()) {
            throw std::runtime_error("Failed to open output file '" + output_path
                                     + "': " + output_result.status().ToString());
        }
        auto output_stream = output_result.ValueOrDie();

        // Configure Parquet writer properties
        parquet::WriterProperties::Builder props_builder;
        props_builder.compression(static_cast<arrow::Compression::type>(
                get_arrow_compression(compression)));

        // Use dictionary encoding for string columns (better compression)
        props_builder.enable_dictionary();

        auto writer_properties = props_builder.build();

        // Configure Arrow writer properties
        auto arrow_props = parquet::ArrowWriterProperties::Builder().build();

        // Create Parquet file writer using Result-returning API (modern)
        auto writer_result = parquet::arrow::FileWriter::Open(
                *record_batch->schema(),
                arrow::default_memory_pool(),
                output_stream,
                writer_properties,
                arrow_props);

        if (!writer_result.ok()) {
            throw std::runtime_error("Failed to create Parquet writer: "
                                     + writer_result.status().ToString());
        }

        auto writer = std::move(writer_result).ValueOrDie();

        // Write the record batch
        auto write_status = writer->WriteRecordBatch(*record_batch);
        if (!write_status.ok()) {
            throw std::runtime_error("Failed to write record batch: "
                                     + write_status.ToString());
        }

        // Close the writer
        auto close_status = writer->Close();
        if (!close_status.ok()) {
            throw std::runtime_error("Failed to close Parquet writer: "
                                     + close_status.ToString());
        }

        // Close the output stream
        auto stream_close_status = output_stream->Close();
        if (!stream_close_status.ok()) {
            throw std::runtime_error("Failed to close output stream: "
                                     + stream_close_status.ToString());
        }

    } catch (const std::exception &e) {
        std::ostringstream oss;
        oss << "Error writing Parquet file '" << output_path << "': " << e.what();
        throw std::runtime_error(oss.str());
    }
}

} // namespace igor
