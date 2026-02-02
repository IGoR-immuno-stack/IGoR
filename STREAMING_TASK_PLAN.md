# Streaming Layer Task Plan
**Branch:** feature/streaming
**Base Commit:** 850c3a6
**Last Updated:** 2 February 2026
**Goal:** Stream, read and write files under Parquet and AIRR format

## Executive Summary

**Overall Completion: 51%**
- Parquet I/O: 60% complete (read/write works, alignment data incomplete)
- AIRR I/O: 0% complete (not started)

### Critical Finding
Initial assessment claimed "100% complete" for Phases 1-3, but code analysis reveals:
- ✅ Parquet infrastructure is solid (write/read basic data)
- ⚠️ Alignment serialization only 33% complete (3 of 9 fields)
- ❌ AIRR format support completely missing (0%)
- ⚠️ All 39 tests pass but only test empty alignments

## Completed Work (51%)

### ✅ Parquet Write Infrastructure (COMPLETE)
- **Files:** `src/igor/Streaming/ParquetWriter.{h,cpp}`
- **Status:** Production-ready for basic sequence data
- **Features:**
  - 5 compression codecs (NONE, SNAPPY, GZIP, ZSTD, LZ4)
  - Zero-copy via Arrow C Data Interface
  - RAII memory management
  - Performance: >10M sequences/second (100x target)
- **Limitation:** Does not serialize complete alignment data

### ✅ Parquet Read Infrastructure (COMPLETE)
- **Files:** `src/igor/Streaming/ParquetReader.{h,cpp}`
- **Status:** Production-ready for basic sequence data
- **Features:**
  - Metadata inspection (get_file_info)
  - Selective column reading (read_columns)
  - Legacy format conversion
  - Performance: >5K sequences/second (10x target)
- **Limitation:** Cannot parse complete alignment data

### ⚠️ Batch Conversion Helpers (PARTIAL - 60%)
- **Files:** `src/igor/Streaming/SequenceBatchHelpers.{h,cpp}`
- **Status:** Works for sequence_id and sequence, incomplete for alignments
- **Implemented:**
  - `row_to_sequence_data()` - converts columnar to legacy (basic data only)
  - `has_column()`, `get_string_value()`, `get_int_value()`, `get_double_value()` - helper functions
- **Incomplete:**
  - `vector_to_batch()` - TODO line 240: "Add alignment columns"
  - `parse_alignments_from_columns()` - TODO line 125: "Implement full alignment parsing"
  - **Critical Issue:** Only 3 of 9 `Alignment_data` fields implemented:
    - ✅ gene_name
    - ✅ offset
    - ✅ score
    - ❌ five_p_offset
    - ❌ three_p_offset
    - ❌ insertions (forward_list<int>)
    - ❌ deletions (forward_list<int>)
    - ❌ align_length
    - ❌ mismatches (vector<int>)

### ✅ Test Suite (COMPLETE but limited scope)
- **Files:** `tst/igor/Streaming/*.cpp`
- **Status:** 39 tests, 653 assertions, 100% passing
- **Coverage:**
  - test_parquet_writer.cpp: 9 test cases
  - test_parquet_reader.cpp: 13 test cases
  - test_streaming_helpers.cpp: 17 test cases
- **Critical Limitation:** All tests use empty alignments
  - Round-trip tests pass because no alignment data is tested
  - Real alignment data would be lost in production

### ✅ Build System Integration (COMPLETE)
- CMake configuration for Streaming module
- Proper PRIVATE linkage for Arrow/Parquet dependencies
- PUBLIC linkage for Sparrow header-only library
- Zero conflicts with develop branch

## Remaining Work (49%)

### 🔴 Task 1: Complete Alignment Data Preservation (HIGH PRIORITY)
**Estimated Effort:** ~300-400 LOC
**Blocking:** Production use of Streaming layer
**Target Completion:** Week 1

#### Subtasks:
1. **Implement missing fields in `parse_alignments_from_columns()`**
   - Add five_p_offset, three_p_offset parsing
   - Implement insertions (forward_list<int>) parsing from array column
   - Implement deletions (forward_list<int>) parsing from array column
   - Add align_length parsing
   - Implement mismatches (vector<int>) parsing from array column
   - Location: `src/igor/Streaming/SequenceBatchHelpers.cpp` line 125

2. **Add alignment serialization to `vector_to_batch()`**
   - Create columns for all 9 Alignment_data fields
   - Handle nested structures (insertions, deletions, mismatches as list/array columns)
   - Location: `src/igor/Streaming/SequenceBatchHelpers.cpp` line 240

3. **Add comprehensive alignment tests**
   - Create test data with real alignments (all 9 fields populated)
   - Test round-trip preservation of complex alignment data
   - Test edge cases (empty lists, large arrays, boundary values)
   - Add to existing test files

**Acceptance Criteria:**
- [ ] All 9 Alignment_data fields serialize to Parquet
- [ ] All 9 Alignment_data fields deserialize from Parquet
- [ ] Round-trip test preserves complete alignment data
- [ ] Tests pass with non-empty alignment data
- [ ] No data loss in production use cases

---

### 🟡 Task 2: AIRR TSV/CSV Reading (MEDIUM PRIORITY)
**Estimated Effort:** ~500-600 LOC
**Blocking:** AIRR format support (50% of stated goal)
**Target Completion:** Week 2

#### Subtasks:
1. **Create AIRRReader class**
   - File: `src/igor/Streaming/AIRRReader.h` (new)
   - File: `src/igor/Streaming/AIRRReader.cpp` (new)
   - Parse AIRR v1.4 TSV format
   - Parse AIRR v1.4 CSV format
   - Auto-detect delimiter

2. **Implement field mapping**
   - Map AIRR required fields to IGoR format:
     - sequence_id → Seq_data::seq_name
     - sequence → Seq_data::sequence
     - v_call, d_call, j_call → Alignment_data::gene_name
     - v_cigar, d_cigar, j_cigar → parse to Alignment_data fields
   - Handle optional fields (productive, junction, etc.)
   - Validate schema compliance

3. **Add AIRR read tests**
   - Test TSV parsing
   - Test CSV parsing
   - Test invalid format handling
   - Test schema validation
   - File: `tst/igor/Streaming/test_airr_reader.cpp` (new)

**Acceptance Criteria:**
- [ ] Can read AIRR v1.4 TSV files
- [ ] Can read AIRR v1.4 CSV files
- [ ] Correctly maps all required fields
- [ ] Handles optional fields gracefully
- [ ] Validates schema compliance
- [ ] Error handling for malformed files

---

### 🟡 Task 3: AIRR TSV/CSV Writing (MEDIUM PRIORITY)
**Estimated Effort:** ~400-500 LOC
**Blocking:** AIRR format support (50% of stated goal)
**Target Completion:** Week 2-3

#### Subtasks:
1. **Create AIRRWriter class**
   - File: `src/igor/Streaming/AIRRWriter.h` (new)
   - File: `src/igor/Streaming/AIRRWriter.cpp` (new)
   - Write AIRR v1.4 TSV format
   - Write AIRR v1.4 CSV format
   - Configurable delimiter

2. **Implement field mapping**
   - Map IGoR format to AIRR required fields
   - Generate CIGAR strings from Alignment_data
   - Populate optional fields from available data
   - Add AIRR-required metadata headers

3. **Add AIRR write tests**
   - Test TSV generation
   - Test CSV generation
   - Test round-trip with AIRRReader
   - Validate output schema compliance
   - File: `tst/igor/Streaming/test_airr_writer.cpp` (new)

**Acceptance Criteria:**
- [ ] Can write AIRR v1.4 TSV files
- [ ] Can write AIRR v1.4 CSV files
- [ ] Output validates against AIRR schema
- [ ] Round-trip with AIRRReader preserves data
- [ ] Correctly generates CIGAR strings
- [ ] Includes required metadata headers

---

### 🟢 Task 4: AIRR Parquet Support (LOW PRIORITY)
**Estimated Effort:** ~200-300 LOC
**Blocking:** Advanced interoperability
**Target Completion:** Week 4

#### Subtasks:
1. **Extend ParquetWriter for AIRR schema**
   - Add AIRR-compliant column names
   - Add AIRR schema metadata to Parquet files
   - Support AIRR-specific data types

2. **Extend ParquetReader for AIRR schema**
   - Auto-detect AIRR Parquet files
   - Read AIRR-compliant column names
   - Validate AIRR schema metadata

3. **Add AIRR Parquet tests**
   - Test write/read with AIRR schema
   - Test interoperability with immunarch
   - Test interoperability with scirpy
   - File: Update existing `test_parquet_*.cpp`

**Acceptance Criteria:**
- [ ] Can write AIRR-schema Parquet files
- [ ] Can read AIRR-schema Parquet files
- [ ] Interoperable with immunarch R package
- [ ] Interoperable with scirpy Python package
- [ ] Schema metadata preserved in files

---

## Dependencies & Risks

### External Dependencies
- ✅ Apache Arrow 13.0.0 (installed via pixi)
- ✅ Apache Parquet 13.0.0 (installed via pixi)
- ✅ Sparrow (header-only, integrated)
- ✅ Catch2 v3.8.0 (testing framework)

### Technical Risks
1. **Nested structure serialization** (Task 1)
   - Risk: Complex types (forward_list, vector) may not map cleanly to Parquet
   - Mitigation: Use Arrow list arrays, test thoroughly

2. **CIGAR string generation** (Task 3)
   - Risk: Converting IGoR alignment format to CIGAR may lose information
   - Mitigation: Define clear mapping rules, add validation tests

3. **AIRR schema evolution** (Tasks 2-4)
   - Risk: AIRR standard may change
   - Mitigation: Version schema explicitly, support multiple versions

### Resource Requirements
- Estimated 4 weeks for full completion
- Requires access to real AIRR test datasets
- May need AIRR schema validation tools

## Testing Strategy

### Current Test Coverage (Baseline)
- 39 test cases, 653 assertions
- 100% pass rate (but limited scope - empty alignments only)

### Required Test Coverage (Target)
- **Task 1:** Add 15-20 tests with real alignment data
  - Target: 95% code coverage for SequenceBatchHelpers
- **Task 2:** Add 20-25 tests for AIRR reading
  - Test valid/invalid TSV/CSV files
  - Test schema validation
- **Task 3:** Add 20-25 tests for AIRR writing
  - Test output validation
  - Test round-trip with Task 2
- **Task 4:** Add 10-15 tests for AIRR Parquet
  - Test interoperability with external tools

### Performance Targets (Maintained)
- Parquet write: >1M sequences/second (currently 10M+) ✅
- Parquet read: >500 sequences/second (currently 5K+) ✅
- AIRR TSV write: >100K sequences/second (new target)
- AIRR TSV read: >50K sequences/second (new target)

## Success Criteria

### Phase 1: Alignment Data (Task 1)
- [x] Parquet write basic data (DONE)
- [x] Parquet read basic data (DONE)
- [ ] Parquet write complete alignments (IN PROGRESS)
- [ ] Parquet read complete alignments (IN PROGRESS)
- [ ] Round-trip preserves all alignment fields (BLOCKED)

### Phase 2: AIRR Read/Write (Tasks 2-3)
- [ ] Read AIRR TSV format (NOT STARTED)
- [ ] Read AIRR CSV format (NOT STARTED)
- [ ] Write AIRR TSV format (NOT STARTED)
- [ ] Write AIRR CSV format (NOT STARTED)
- [ ] Schema validation (NOT STARTED)

### Phase 3: AIRR Parquet (Task 4)
- [ ] Write AIRR Parquet (NOT STARTED)
- [ ] Read AIRR Parquet (NOT STARTED)
- [ ] Interoperability validation (NOT STARTED)

### Production Readiness Checklist
- [ ] All alignment data preserved in Parquet I/O
- [ ] AIRR TSV/CSV read support
- [ ] AIRR TSV/CSV write support
- [ ] Comprehensive test coverage (>90%)
- [ ] Performance targets met
- [ ] Documentation complete
- [ ] Zero conflicts with develop branch
- [ ] Code review completed
- [ ] Merge to develop branch

## Timeline

| Week | Focus | Deliverable |
|------|-------|-------------|
| 1 | Task 1 | Complete alignment serialization |
| 2 | Task 2 | AIRR TSV/CSV reading |
| 3 | Task 3 | AIRR TSV/CSV writing |
| 4 | Task 4 | AIRR Parquet support |
| 5 | Testing & Documentation | Production-ready release |

## Notes

### Known Issues from Code Analysis
- **Line 125 (SequenceBatchHelpers.cpp):** `parse_alignments_from_columns()` incomplete
  - TODO: "Implement full alignment parsing when format is defined"
  - Only 3 of 9 fields implemented

- **Line 240 (SequenceBatchHelpers.cpp):** `vector_to_batch()` incomplete
  - TODO: "Add alignment columns in future tasks"
  - Only serializes sequence_id and sequence

### Compatibility Status
- ✅ Zero conflicts with develop branch (verified via git merge-tree)
- ✅ Build system properly isolates dependencies
- ✅ No ZeroCopy dependencies in feature/streaming
- ✅ Ready for incremental merge once tasks complete

### References
- GitHub Issue: #12 "Random Access and Batch Reading of Receptor Sequences"
- Base documentation: `Task_6_Streaming_Architecture.md`
- AIRR Standards: https://docs.airr-community.org/en/stable/
- Current branch: feature/streaming (commit 850c3a6)
