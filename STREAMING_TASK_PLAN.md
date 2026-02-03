# Streaming Layer Task Plan

**Branch:** feature/streaming  
**Last Updated:** 3 February 2026  
**Goal:** Stream, read and write files under Parquet and AIRR format

---

## Status Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Parquet I/O | ✅ Complete | Production-ready with native Arrow list arrays |
| Test Suite | ✅ Complete | 25 tests, 949 assertions, Catch2 v3 best practices |
| AIRR I/O | ⬜ Not Started | Tasks 2-4 below |

**Overall Completion: ~50%** (Parquet complete, AIRR not started)

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

### 🟡 Task 2: AIRR TSV/CSV Reading

**Estimated Effort:** ~500-600 LOC

**Files to Create:**
- `src/igor/Streaming/AIRRReader.{h,cpp}`
- `tst/igor/Streaming/test_airr_reader.cpp`

**Subtasks:**
1. Parse AIRR v1.4 TSV format
2. Parse AIRR v1.4 CSV format
3. Auto-detect delimiter
4. Map AIRR fields to IGoR format:
   - sequence_id → Seq_data::seq_name
   - sequence → Seq_data::sequence
   - v_call, d_call, j_call → Alignment_data::gene_name
   - v_cigar, d_cigar, j_cigar → parse to Alignment_data fields
5. Validate schema compliance

**Acceptance Criteria:**
- [ ] Can read AIRR v1.4 TSV files
- [ ] Can read AIRR v1.4 CSV files
- [ ] Correctly maps all required fields
- [ ] Error handling for malformed files

---

### 🟡 Task 3: AIRR TSV/CSV Writing

**Estimated Effort:** ~400-500 LOC

**Files to Create:**
- `src/igor/Streaming/AIRRWriter.{h,cpp}`
- `tst/igor/Streaming/test_airr_writer.cpp`

**Subtasks:**
1. Write AIRR v1.4 TSV format
2. Write AIRR v1.4 CSV format
3. Generate CIGAR strings from Alignment_data
4. Add AIRR-required metadata headers
5. Round-trip validation with AIRRReader

**Acceptance Criteria:**
- [ ] Can write AIRR v1.4 TSV files
- [ ] Can write AIRR v1.4 CSV files
- [ ] Output validates against AIRR schema
- [ ] Correctly generates CIGAR strings

---

### 🟢 Task 4: AIRR Parquet Support (Low Priority)

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
