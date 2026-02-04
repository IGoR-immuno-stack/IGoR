# IGoR Streaming Module

This module provides high-performance streaming data I/O for IGoR supporting both Apache Parquet and AIRR Community standards.

## Overview

The Streaming module enables IGoR to read and write immunoglobulin/TCR sequence data in multiple formats for maximum interoperability and performance:

- **Parquet I/O**: Read and write sequence data in Apache Parquet format with compression
- **AIRR I/O**: Full support for AIRR Rearrangement and Alignment schemas (TSV/CSV)
- **Batch Processing**: Work with sequences in columnar batches for better performance
- **Legacy Compatibility**: Convert between traditional IGoR formats and modern columnar layouts
- **Zero-Copy Operations**: Efficient data transfer using Arrow C Data Interface

## Components

### Parquet I/O

#### ParquetWriter

Write IGoR sequence data to Parquet files with optional compression.

**Key Features:**
- Multiple compression formats (SNAPPY, GZIP, ZSTD, LZ4, NONE)
- Dictionary encoding for string columns
- Native Arrow list arrays for insertions/deletions/mismatches
- 29 columns: 2 basic + 9 fields × 3 gene classes (V, D, J)
- Performance: ~1.5M sequences/second (10K batch)

**Example:**
```cpp
#include <igor/Streaming/ParquetWriter.h>

// Write sequences to Parquet file
std::vector<SequenceData> sequences = {
    {1, "ATCG", {}},
    {2, "GCTA", {}}
};
ParquetWriter::write_sequences("output.parquet", sequences, CompressionType::SNAPPY);
```

#### ParquetReader

Read Parquet files into IGoR sequence batches.

**Key Features:**
- Zero-copy reading via Arrow C Data Interface
- Selective column reading
- File metadata inspection
- Legacy format conversion
- Performance: ~1.4M sequences/second (full conversion), ~7M seq/s (batch only)

**Example:**
```cpp
#include <igor/Streaming/ParquetReader.h>

// Read entire file
sparrow::record_batch batch = ParquetReader::read_batch("input.parquet");

// Read specific columns
std::vector<std::string> columns = {"sequence_id", "sequence"};
sparrow::record_batch subset = ParquetReader::read_columns("input.parquet", columns);

// Get file metadata
auto info = ParquetReader::get_file_info("input.parquet");
std::cout << "Sequences: " << info.num_rows << std::endl;

// Convert to SequenceData
std::vector<SequenceData> sequences = ParquetReader::read_sequences("input.parquet");
```

### AIRR I/O

Full support for the AIRR Community standards for data sharing.

#### AIRR Rearrangement Schema

**AIRRRearrangementReader** and **AIRRRearrangementWriter** (`namespace igor::airr::rearrangement`)

- One row per sequence
- V/D/J calls inlined as columns
- TSV and CSV support with auto-detection
- Coordinate conversion (1-based AIRR ↔ 0-based IGoR)

**Example:**
```cpp
#include <igor/Streaming/AIRRRearrangementReader.h>
#include <igor/Streaming/AIRRRearrangementWriter.h>

using namespace igor::airr::rearrangement;

// Read AIRR Rearrangement file
auto sequences = read_sequences("input.tsv");

// Write to AIRR format
write_sequences("output.tsv", sequences, Delimiter::TAB);
```

#### AIRR Alignment Schema

**AIRRAlignmentReader** and **AIRRAlignmentWriter** (`namespace igor::airr::alignment`)

- One row per alignment (multiple rows per sequence)
- Supports ranked alignments (primary, secondary, etc.)
- CIGAR string generation and parsing
- Groups alignments by sequence_id

**Example:**
```cpp
#include <igor/Streaming/AIRRAlignmentReader.h>
#include <igor/Streaming/AIRRAlignmentWriter.h>

using namespace igor::airr::alignment;

// Read AIRR Alignment file
auto sequences = read_sequences("alignments.tsv");

// Write to AIRR Alignment format
write_sequences("output_alignments.tsv", sequences);
```

#### AIRR Common Types

Shared utilities in `AIRRCommon.h`:
- `Delimiter` enum (TAB, COMMA, AUTO)
- `FileInfo` struct (unified for both schemas)
- `delimiter_char()` and `detect_delimiter()` functions

### SequenceBatchHelpers

Utilities for converting between IGoR's internal data structures and Sparrow's columnar format.

**Key Functions:**
- `row_to_sequence_data()` - Extract single sequence from a batch row
- `vector_to_batch()` - Convert legacy vector to columnar batch
- `parse_alignments_from_columns()` - Extract alignment data from columns
- Null-safe column accessors

**Example:**
```cpp
#include <igor/Streaming/SequenceBatchHelpers.h>

// Convert legacy data to batch
std::vector<SequenceData> sequences = load_legacy_data();
sparrow::record_batch batch = vector_to_batch(sequences);

// Process individual rows
for (size_t i = 0; i < batch.nb_rows(); ++i) {
    SequenceData seq = row_to_sequence_data(batch, i);
    process_sequence(seq);
}
```

## Format Comparison

| Format | Use Case | Performance | Interoperability |
|--------|----------|-------------|------------------|
| **Parquet** | Large datasets, internal storage | ★★★★★ Fast | Arrow/pandas/R |
| **AIRR Rearrangement** | Data sharing, one sequence per row | ★★★☆☆ Medium | immunarch, scirpy |
| **AIRR Alignment** | Multiple alignments, ranked results | ★★★☆☆ Medium | AIRR-compliant tools |

## Building

The streaming module is built automatically when you build IGoR:

```bash
pixi run cmake --build build
```

The module requires:
- Apache Arrow 13.0.0 or later
- Apache Parquet 13.0.0 or later
- Sparrow (header-only)

## Testing

Run the streaming tests:

```bash
pixi run ./build/bin/streaming_tests
```

Or run specific test categories:

```bash
# Parquet tests only
pixi run ./build/bin/streaming_tests "[parquet]"

# AIRR rearrangement tests
pixi run ./build/bin/streaming_tests "[airr][rearrangement]"

# AIRR alignment tests
pixi run ./build/bin/streaming_tests "[airr][alignment]"

# Integration tests (real Murugan dataset)
pixi run ./build/bin/streaming_tests "[.integration]"

# Benchmarks
pixi run ./build/bin/streaming_tests "[.benchmark]"
```

The test suite includes:
- **50 test cases**, **1131 assertions**
- Parquet write/read round-trip validation
- AIRR Rearrangement schema tests
- AIRR Alignment schema tests
- Compression format verification
- Performance benchmarks
- Legacy format compatibility tests
- Error handling validation

## Performance

Benchmarks on typical hardware (Apple M-series):

| Operation | Throughput | Notes |
|-----------|-----------|-------|
| Parquet Write | ~1.5M seq/s | SNAPPY compression, 10K batch |
| Parquet Read (batch) | ~7M seq/s | Batch only, no conversion |
| Parquet Read (full) | ~1.4M seq/s | Full SequenceData conversion |
| AIRR Rearrangement I/O | ~100K seq/s | TSV format |
| AIRR Alignment I/O | ~100K seq/s | TSV format, grouped by sequence |

## Architecture

The module uses a layered architecture:

1. **Public API Layer**
   - `ParquetWriter`, `ParquetReader` - Parquet I/O
   - `AIRRRearrangementReader`, `AIRRRearrangementWriter` - AIRR Rearrangement schema
   - `AIRRAlignmentReader`, `AIRRAlignmentWriter` - AIRR Alignment schema
   - `SequenceBatchHelpers` - Conversion utilities

2. **Common Layer**
   - `AIRRCommon` - Shared AIRR utilities
   - Delimiter detection and handling

3. **Integration Layer**
   - Arrow C Data Interface (zero-copy)
   - Sparrow columnar types

4. **Storage Layer**
   - Apache Parquet (binary columnar)
   - AIRR TSV/CSV (text-based)

Dependencies on Apache Arrow and Parquet are private implementation details - users only need to work with the public API and Sparrow types.

## Known Limitations

### AIRR Alignment Schema

1. **CIGAR Simplification**: CIGAR strings are simplified to counts (e.g., `50M2I3D`). Per-base positions are not preserved. Use Parquet for lossless storage.

2. **Sequence Field Not Stored**: The AIRR Alignment schema doesn't include a `sequence` column. Sequence data must be joined from a Rearrangement file.

3. **Optional Fields**: `identity`, `support`, `germline_end` fields are not implemented (not in IGoR's `Alignment_data` structure).

4. **C Gene Segment**: Constant (C) gene segments are not supported (not in IGoR's `Gene_class` enum).

5. **Sequence ID Parsing**: Non-numeric sequence IDs are hashed to integers; original strings are not preserved.

See `STREAMING_TASK_PLAN.md` for detailed limitation descriptions and workarounds.

## References

- [AIRR Community Standards](https://docs.airr-community.org/en/stable/)
- [Apache Arrow](https://arrow.apache.org/)
- [Apache Parquet](https://parquet.apache.org/)
- [Sparrow](https://github.com/man-group/sparrow)
