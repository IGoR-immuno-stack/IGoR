# Streaming Layer Task Plan

**Branch:** feature/streaming
**Last Updated:** 4 February 2026

---

## Objective

Enable IGoR to read and write receptor sequence data in standard formats for interoperability with the immunology research ecosystem.

**Target Formats:**
- **Parquet** — High-performance binary format for large datasets, interoperable with Arrow/pandas/R
- **AIRR** — AIRR-Community standard (TSV/CSV) for sharing data between immune repertoire tools

**Why This Matters:**
- **Data Exchange** — Share sequence/alignment data with tools like immunarch (R), scirpy (Python), and other AIRR-compliant software
- **Performance** — Parquet enables fast I/O for large datasets (millions of sequences)
- **Preservation** — Full round-trip of all 9 alignment fields (gene name, offset, score, insertions, deletions, mismatches, etc.)

---

## Status Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Parquet I/O | ✅ Complete | Production-ready with native Arrow list arrays |
| AIRR Rearrangement I/O | ✅ Complete | Reader & Writer with full test coverage |
| AIRR Common Types | ✅ Complete | Shared utilities for all AIRR formats |
| AIRR Alignment I/O | ✅ Complete | Reader & Writer with full test coverage |
| Test Suite | ✅ Complete | 50 tests, 1131 assertions |

**Overall Completion: ~95%** (Core functionality complete, minor optional features pending)

---

## Completed Work

### ✅ Parquet Writer (`ParquetWriter.{h,cpp}`)

Production-ready with:
- 5 compression codecs (NONE, SNAPPY, GZIP, ZSTD, LZ4)
- Zero-copy via Arrow C Data Interface
- Native Arrow `list<int32>` arrays for insertions/deletions/mismatches
- 29 columns: 2 basic + 9 fields × 3 gene classes (V, D, J)
- Performance: ~1.5M sequences/second (10K batch)

### ✅ Parquet Reader (`ParquetReader.{h,cpp}`)

Production-ready with:
- Metadata inspection (`get_file_info`)
- Selective column reading (`read_columns`)
- Legacy format conversion (`read_sequences`)
- Native list array deserialization
- Performance: ~1.4M sequences/second (full conversion), ~7M seq/s (batch only)

### ✅ Batch Helpers (`SequenceBatchHelpers.{h,cpp}`)

All 9 `Alignment_data` fields supported:

| Field | Arrow Type | Notes |
|-------|-----------|-------|
| gene_name | string | |
| offset | int32 | |
| score | double | |
| five_p_offset | uint64 | |
| three_p_offset | uint64 | |
| align_length | uint64 | |
| insertions | list\<int32\> | Native Arrow list array |
| deletions | list\<int32\> | Native Arrow list array |
| mismatches | list\<int32\> | Native Arrow list array |

### ✅ Test Suite

**Architecture:**
- `test_streaming_helpers.cpp` — SequenceBatchHelpers tests
- `test_parquet_reader.cpp` — ParquetReader tests
- `test_parquet_writer.cpp` — ParquetWriter tests
- `StreamingTestUtils.{h,cpp}` — Shared utilities

**Catch2 v3 Features Used:**
- `TEST_CASE_METHOD` with RAII fixtures
- `GENERATE` for parametric tests
- `SECTION` for grouping related tests
- Matchers (`IsEmpty()`, `Contains()`, `RangeEquals()`)
- `CAPTURE()` for failure diagnostics
- Hidden tests: `[.integration]`, `[.benchmark]`

**Shared Utilities:**
- `TestDirectory` — RAII wrapper with automatic cleanup
- `SequenceTuple` — Type alias for verbose legacy tuple
- `create_test_batch()` — Inline test batch creation
- `create_test_sequences()` — Synthetic sequence generation
- `create_sequence_with_v_alignment()` — Alignment test data
- `load_murugan_dataset()` — Real data loader (300 seqs, 73K alignments)

**Test Data:**
- Location: `tst/igor/Streaming/test_data/`
- Files:
  - `murugan_naive1_noncoding_demo_seqs_indexed_seq.csv` (300 sequences)
  - `murugan_naive1_noncoding_demo_seqs_alignments_V.csv` (6.5K alignments)
  - `murugan_naive1_noncoding_demo_seqs_alignments_D.csv` (59K alignments)
  - `murugan_naive1_noncoding_demo_seqs_alignments_J.csv` (7.6K alignments)

**Running Tests:**
```bash
# All tests (excluding hidden)
pixi run ./build/bin/streaming_tests

# Integration tests (real Murugan data)
pixi run ./build/bin/streaming_tests "[.integration]"

# Benchmarks
pixi run ./build/bin/streaming_tests "[.benchmark]"
```

---

## Remaining Work

### ✅ Task 2: AIRR TSV/CSV Reading (COMPLETE)

**Actual Effort:** ~700 LOC (header + implementation + tests)
**Completion Date:** February 4, 2026

**Files Created:**
- `src/igor/Streaming/AIRRReader.{h,cpp}` (290 + 450 lines)
- `tst/igor/Streaming/test_airr_reader.cpp` (340 lines)

**Features Implemented:**
1. ✅ Parse AIRR v1.4 TSV format
2. ✅ Parse AIRR v1.4 CSV format
3. ✅ Auto-detect delimiter (tabs vs commas)
4. ✅ Map AIRR fields to IGoR format:
   - sequence_id → SequenceData::index
   - sequence → SequenceData::sequence
   - v_call, d_call, j_call → Alignment_data::gene_name
   - v_score, d_score, j_score → Alignment_data::score
   - v_sequence_start, etc. → Alignment_data::offset (1-based → 0-based conversion)
   - v_alignment_length, etc. → Alignment_data::align_length
5. ✅ Schema validation (checks required columns)
6. ✅ Multiple output formats: batch, SequenceData, legacy tuples, selective columns

**Test Coverage:**
- 9 test cases, 65 assertions
- Tests for: delimiter detection, file info, schema validation, batch reading,
  sequence reading with V/D/J alignments, legacy format, column selection, error handling

**Acceptance Criteria:**
- [x] Can read AIRR v1.4 TSV files
- [x] Can read AIRR v1.4 CSV files
- [x] Correctly maps all required fields
- [x] Error handling for malformed files

---

### ✅ Task 3: AIRR TSV/CSV Writing (COMPLETE)

**Actual Effort:** ~480 LOC (header + implementation + tests)
**Completion Date:** February 4, 2026

**Files Created:**
- `src/igor/Streaming/AIRRWriter.{h,cpp}` (190 + 270 lines)
- `tst/igor/Streaming/test_airr_writer.cpp` (400 lines)

**Features Implemented:**
1. ✅ Write AIRR v1.4 TSV format (`write_tsv()`)
2. ✅ Write AIRR v1.4 CSV format (`write_csv()`)
3. ✅ Generate simplified CIGAR strings from Alignment_data (`make_cigar()`)
4. ✅ Write column headers matching AIRR schema
5. ✅ Round-trip validation with AIRRReader
6. ✅ Support for legacy tuple format (`write_legacy_tsv()`)
7. ✅ Direct batch writing (`write_batch()`)

**Test Coverage:**
- 7 test cases, 41 assertions (optimized from 45)
- Tests for: TSV writing, CSV writing, round-trip validation,
  CIGAR generation (match, insertions, deletions, combined), column headers,
  legacy format support, error handling

**Code Review Fixes Applied:**
1. ✅ Removed unused `escape_csv()` function
2. ✅ Fixed CIGAR logic: `align_length` used directly (biologically correct)
3. ✅ Added limitation note for simplified CIGAR in docstring
4. ✅ Added string-column requirement note for `write_batch()` in docstring
5. ✅ Replaced manual loop with `std::distance()` for counting indels

**Acceptance Criteria:**
- [x] Can write AIRR v1.4 TSV files
- [x] Can write AIRR v1.4 CSV files
- [x] Output validates against AIRR schema
- [x] Correctly generates CIGAR strings

---

### � Task 3.5: AIRR Alignment Schema Support

**Estimated Effort:** ~500-600 LOC
**Priority:** Medium (before Parquet support)

**Background:**
The AIRR Alignment Schema (experimental) differs from the Rearrangement Schema:
- **Rearrangement**: One row per sequence, V/D/J calls inlined as columns
- **Alignment**: One row per alignment (multiple rows per sequence, one per V/D/J/C)

This enables storing multiple candidate alignments per gene class with ranking.

**✅ Part A: Rename Existing Files & Refactor Common Types (COMPLETE)**
**Completion Date:** February 4, 2026

- ✅ Renamed `AIRRReader.{h,cpp}` → `AIRRRearrangementReader.{h,cpp}`
- ✅ Renamed `AIRRWriter.{h,cpp}` → `AIRRRearrangementWriter.{h,cpp}`
- ✅ Updated namespace from `igor::airr` → `igor::airr::rearrangement`
- ✅ Updated `test_airr_reader.cpp` → `test_airr_rearrangement_reader.cpp`
- ✅ Updated `test_airr_writer.cpp` → `test_airr_rearrangement_writer.cpp`
- ✅ Updated CMakeLists.txt references
- ✅ Created `AIRRCommon.{h,cpp}` with shared types:
  - `enum class Delimiter { TAB, COMMA, AUTO }`
  - `struct FileInfo` (unified for both schemas)
  - `char delimiter_char(Delimiter)`
  - `Delimiter detect_delimiter(const std::string&)`
- ✅ Refactored existing files to use shared types from `igor::airr` namespace
- ✅ All tests pass (106 assertions in 16 test cases)

**Namespace Organization:**
- `igor::airr` — Shared types and utilities
- `igor::airr::rearrangement` — Rearrangement schema (one row per sequence)
- `igor::airr::alignment` — Alignment schema (one row per alignment)

**✅ Part B: Implement Alignment Schema (COMPLETE)**
**Completion Date:** February 4, 2026

**Files Created:**
- ✅ `src/igor/Streaming/AIRRAlignmentReader.{h,cpp}` (~185 + 444 lines)
- ✅ `src/igor/Streaming/AIRRAlignmentWriter.{h,cpp}` (~139 + 174 lines)
- ✅ `tst/igor/Streaming/test_airr_alignment.cpp` (~480 lines, 9 test cases, 76 assertions)

**AIRR Alignment Schema Fields:**
| Field | Type | Required | IGoR Mapping | Status |
|-------|------|----------|--------------|--------|
| `sequence_id` | string | ✅ | `SequenceData::index` | ✅ |
| `segment` | string | ✅ | Gene class (V/D/J) | ✅ |
| `call` | string | ✅ | `Alignment_data::gene_name` | ✅ |
| `score` | number | ✅ | `Alignment_data::score` | ✅ |
| `cigar` | string | ✅ | Generated/parsed | ✅ |
| `sequence_start` | integer | | `offset + 1` (1-based) | ✅ |
| `sequence_end` | integer | | `offset + align_length` | ✅ |
| `germline_start` | integer | | `five_p_offset + 1` | ✅ |
| `rank` | integer | | Vector position | ✅ |
| `identity` | number | | Fractional identity | ⬜ |
| `support` | number | | E-value/p-value | ⬜ |
| `germline_end` | integer | | Calculated | ⬜ |

**Implementation Status:**
- ✅ Created `igor::airr::alignment` namespace
- ✅ Implemented alignment reader with segment parsing
- ✅ Implemented alignment writer (one row per alignment)
- ✅ Handles ranked alignments (multiple per gene class)
- ✅ CIGAR string generation and parsing
- ✅ Uses shared types from `AIRRCommon.h`
- ✅ Field mapping and coordinate conversion (0-based ↔ 1-based)
- ✅ Updated CMakeLists.txt
- ✅ Library compiles successfully
- ✅ All tests pass (1131 assertions in 50 test cases)
- ✅ Test suite for alignment schema (9 test cases, 76 assertions)
- ✅ Round-trip validation tests

**Key Features:**
- Groups multiple alignments per gene class by sequence_id
- Generates rank field based on alignment order in vector
- Simplified CIGAR support (counts only, not per-base positions)
- Validates required columns (sequence_id, segment, call)
- Supports both TSV and CSV formats

---

### ⚠️ Known Limitations

The following limitations exist in the current implementation and may be addressed in future work:

#### 1. CIGAR String Simplification
**Affected:** `AIRRAlignmentReader::parse_cigar()`, `AIRRAlignmentWriter::make_cigar()`

- **Limitation:** CIGAR strings are simplified to aggregate counts only (e.g., `50M2I3D`)
- **Issue:** Per-base insertion/deletion positions are lost during round-trip
- **Impact:** The exact positions of insertions and deletions cannot be reconstructed from the CIGAR string
- **Workaround:** Use Parquet format for lossless storage of all alignment fields

#### 2. Alignment Schema: Sequence Field Not Stored
**Affected:** `AIRRAlignmentReader::read_sequences()`

- **Limitation:** The AIRR Alignment schema does not include a `sequence` column
- **Issue:** `SequenceData::sequence` is set to empty string when reading alignment files
- **Impact:** Sequence data must be joined from a separate source (e.g., Rearrangement file)
- **Workaround:** Use alignment files alongside rearrangement files or Parquet

#### 3. Optional AIRR Fields Not Implemented
**Affected:** Both reader and writer

- **Not implemented:** `identity`, `support`, `germline_end`
- **Reason:** These fields are not used in IGoR's internal `Alignment_data` structure
- **Impact:** Cannot round-trip these fields from external AIRR files

#### 4. C Gene Segment Not Supported
**Affected:** `AIRRAlignmentReader::parse_segment()`

- **Limitation:** Constant (C) gene segments are ignored
- **Reason:** IGoR's `Gene_class` enum does not include C genes
- **Impact:** Alignment rows with `segment=C` are silently skipped

#### 5. Sequence ID Parsing
**Affected:** `AIRRAlignmentReader::read_sequences()`

- **Limitation:** Non-numeric sequence IDs are hashed to integers
- **Issue:** Original string IDs are not preserved in `SequenceData::index`
- **Impact:** May require changes if string-based sequence IDs are needed

#### 6. write_batch() Limited String Support
**Affected:** `AIRRAlignmentWriter::write_batch()`

- **Limitation:** Assumes all columns are string-convertible via `get_string_value()`
- **Impact:** Numeric types may have formatting differences on round-trip

---

### �🟢 Task 4: AIRR Parquet Support (Low Priority)

**Estimated Effort:** ~200-300 LOC

**Subtasks:**
1. Add AIRR-compliant column names to ParquetWriter
2. Add AIRR schema metadata to Parquet files
3. Auto-detect AIRR Parquet files in ParquetReader
4. Test interoperability with immunarch (R) and scirpy (Python)

**Acceptance Criteria:**
- [ ] Can write AIRR-schema Parquet files
- [ ] Can read AIRR-schema Parquet files
- [ ] Interoperable with external tools

---

## Technical Notes

### Native Arrow List Arrays

List fields use native Arrow `list<int32>` arrays with targeted `detail::array_access` for child naming:

```cpp
// Create native list array
sparrow::array arr(sparrow::build(std::move(int_vectors)));
arr.set_name(column_name);

// Name child "item" (Arrow C interface requirement)
auto& proxy = sparrow::detail::array_access::get_arrow_proxy(arr);
if (!proxy.children().empty()) {
    proxy.children()[0].set_name("item");
}
```

Benefits:
- 20-40% smaller files vs string serialization
- Proper Arrow data model
- Minimal internal API usage (only for child naming)

---

## Dependencies

| Dependency | Version | Status |
|------------|---------|--------|
| Apache Arrow | 13.0.0 | ✅ Installed via pixi |
| Apache Parquet | 13.0.0 | ✅ Installed via pixi |
| Sparrow | latest | ✅ Header-only |
| Catch2 | 3.x | ✅ Installed via pixi |

---

## References

- GitHub Issue: #12 "Random Access and Batch Reading of Receptor Sequences"
- AIRR Standards: https://docs.airr-community.org/en/stable/
- Branch: feature/streaming
