# IGoR Streaming Module

This module provides high-performance streaming data I/O for IGoR using Apache Sparrow and Apache Parquet.

## Overview

The Streaming module enables IGoR to read and write immunoglobulin/TCR sequence data in efficient columnar formats. It provides:

- **Parquet I/O**: Read and write sequence data in Apache Parquet format with compression
- **Batch Processing**: Work with sequences in columnar batches for better performance
- **Legacy Compatibility**: Convert between traditional IGoR formats and modern columnar layouts
- **Zero-Copy Operations**: Efficient data transfer using Arrow C Data Interface

## Components

### ParquetWriter

Write IGoR sequence data to Parquet files with optional compression.

**Key Features:**
- Multiple compression formats (SNAPPY, GZIP, ZSTD, LZ4, NONE)
- Dictionary encoding for string columns
- Static API for simple usage
- Performance: >10M sequences/second

**Example:**
```cpp
#include <Streaming/ParquetWriter.h>

// Write sequences to Parquet file
std::vector<SequenceData> sequences = {
    {1, "ATCG", {}},
    {2, "GCTA", {}}
};
sparrow::record_batch batch = vector_to_batch(sequences);
ParquetWriter::write("output.parquet", batch, ParquetWriter::Compression::SNAPPY);
```

### ParquetReader

Read Parquet files into IGoR sequence batches.

**Key Features:**
- Zero-copy reading via Arrow C Data Interface
- Selective column reading
- File metadata inspection
- Legacy format conversion
- Performance: >5K sequences/second

**Example:**
```cpp
#include <Streaming/ParquetReader.h>

// Read entire file
sparrow::record_batch batch = ParquetReader::read_batch("input.parquet");

// Read specific columns
std::vector<std::string> columns = {"sequence_index", "sequence"};
sparrow::record_batch subset = ParquetReader::read_columns("input.parquet", columns);

// Get file metadata
ParquetReader::FileInfo info = ParquetReader::get_file_info("input.parquet");
std::cout << "Sequences: " << info.num_rows << std::endl;

// Convert to legacy format
std::vector<SequenceData> sequences = ParquetReader::read_sequences("input.parquet");
```

### SequenceBatchHelpers

Utilities for converting between IGoR's internal data structures and Sparrow's columnar format.

**Key Functions:**
- `row_to_sequence_data()` - Extract single sequence from a batch row
- `vector_to_batch()` - Convert legacy vector to columnar batch
- `parse_alignments_from_columns()` - Extract alignment data from columns
- Null-safe column accessors

**Example:**
```cpp
#include <Streaming/SequenceBatchHelpers.h>

// Convert legacy data to batch
std::vector<SequenceData> sequences = load_legacy_data();
sparrow::record_batch batch = vector_to_batch(sequences);

// Process individual rows
for (size_t i = 0; i < batch.nb_rows(); ++i) {
    SequenceData seq = row_to_sequence_data(batch, i);
    process_sequence(seq);
}
```

## Building

The streaming module is built automatically when you build IGoR:

```bash
cd build
cmake ..
make
```

The module requires:
- Apache Arrow 13.0.0 or later
- Apache Parquet 13.0.0 or later
- Sparrow (header-only)

## Testing

Run the streaming tests:

```bash
cd build
ctest -R streaming_tests -V
```

Or run the test executable directly:

```bash
./bin/streaming_tests
```

The test suite includes:
- Parquet write/read round-trip validation
- Compression format verification
- Performance benchmarks
- Legacy format compatibility tests
- Error handling validation

## Performance

Benchmarks on typical hardware (Apple M-series):

| Operation | Throughput | Notes |
|-----------|-----------|-------|
| Parquet Write | >10M seq/s | SNAPPY compression |
| Parquet Read | >5K seq/s | Full batch loading |
| Batch Conversion | >100K seq/s | Legacy ↔ Sparrow |

## Architecture

The module uses a layered architecture:

1. **Public API** (`ParquetWriter`, `ParquetReader`, helpers) - Static methods, simple interface
2. **Conversion Layer** (`SequenceBatchHelpers`) - Format transformations
3. **Arrow Integration** - Zero-copy via C Data Interface (private implementation)
4. **Storage** - Apache Parquet columnar format

Dependencies on Apache Arrow and Parquet are private implementation details - users only need to work with the public API and Sparrow types.
