# CIGAR / SeqAn2 Implementation Conclusion

## Summary

The `CIGAR_plan.md` implementation has been completed on branch `feature/CIGAR` and pushed to `origin/feature/CIGAR`.

The implementation covers all planned phases:

- Phase 0: CSV reader fix
- Phase 1: CIGAR parsing and conversion helpers
- Phase 2: Optional SeqAn2 alignment backend
- Phase 3: `AlignmentPreset` and band calibration helpers
- Phase 4: Benchmark harness with legacy, banded SeqAn2, and unbanded SeqAn2 modes

`CIGAR_plan.md` was intentionally left untracked, per instruction.

## Verification

The SeqAn2-enabled test suite was run successfully:

```bash
./build-seqan2/bin/igor_tests "[cigar],[preset],[seqan2],[calibration]"
```

Result:

```text
All tests passed (86 assertions in 12 test cases)
```

The benchmark harness also builds and runs:

```bash
./build-seqan2/bin/bench_seqan2_vs_igor \
  demo/murugan_naive1_noncoding_demo_seqs.txt \
  1 \
  demo/genomicVs_with_primers.fasta \
  models/heavy_pen_substitution_matrix.csv
```

Example output from the improved real-germline benchmark:

```json
{
  "input_sequences": 300,
  "repeat": 1,
  "benchmark_sequences": 300,
  "germline_templates": 89,
  "avg_read_length": 60,
  "avg_germline_length": 286.798,
  "germline_fasta": "demo/genomicVs_with_primers.fasta",
  "legacy_ms": 1451,
  "legacy_sequences": 300,
  "legacy_alignments": 27903,
  "legacy_best_score": 274,
  "seqan2_banded_ms": 261,
  "seqan2_banded_sequences": 300,
  "seqan2_banded_alignments": 26700,
  "seqan2_banded_band_lower_diag": -400,
  "seqan2_banded_band_upper_diag": 400,
  "seqan2_banded_best_score": 8,
  "seqan2_unbanded_ms": 199,
  "seqan2_unbanded_sequences": 300,
  "seqan2_unbanded_alignments": 26700,
  "seqan2_unbanded_best_score": 8,
  "score_delta_seqan2_unbanded_minus_legacy": -266
}
```

## Implemented components

### CSV reader fix

`read_alignments_seq_csv` now reads or reconstructs:

- `align_length`
- `five_p_offset`
- `three_p_offset`

This enables reliable CIGAR reconstruction after CSV round-trips.

### CIGAR conversion API

Added free functions:

```cpp
parse_cigar(...)
alignment_data_to_cigar(...)
alignment_data_from_cigar(...)
alignment_data_sequence_start(...)
alignment_data_sequence_end(...)
alignment_data_germline_start(...)
alignment_data_germline_end(...)
```

These are covered by dedicated Catch2 tests.

### AlignmentPreset

Added `AlignmentPreset.h` with protocol-aware factory methods for:

- V gene, 5' RACE
- V gene, multiplex PCR
- J gene, C-region primer
- J gene, J-gene primer
- D gene
- symmetric custom band

Also added `BandCalibrationParams` and calibration helper APIs.

### SeqAn2 backend

Added guarded SeqAn2 backend files:

- `src/igor/Core/SeqAn2Aligner.h`
- `src/igor/Core/SeqAn2Aligner.cpp`

The backend is enabled only with:

```cmake
-DIGOR_WITH_SEQAN2=ON
```

The backend now uses SeqAn2 `Iupac` sequences and a full 15×15 IGoR substitution matrix loaded into SeqAn2's `ScoreMatrix<Iupac>`.

D-gene/local alignment uses SeqAn2 `LocalAlignmentEnumerator`:

- banded Waterman-Eggert enumeration when a valid band is configured;
- unbanded Waterman-Eggert enumeration when `band_lower_diag > band_upper_diag`.

Pixi/CMake integration was added for SeqAn2. `bioconda::seqan` was available as `noarch`; the `sbl` channel was checked but did not provide `seqan`.

### Band calibration

Calibration helpers now include:

- diagonal estimation from full-matrix seed alignments;
- percentile-based initial band construction;
- score-equivalence verification against full-matrix scores;
- progressive symmetric band widening when the band clips optimal alignments.

The calibration API is covered by dedicated tests.

### Benchmark harness

Added:

```text
tst/bench_seqan2_vs_igor.cpp
```

This builds as a separate target and emits JSON-style timing and comparison output.

The benchmark now supports:

```bash
bench_seqan2_vs_igor [reads] [repeat] [germline_fasta] [substitution_matrix]
```

Default benchmark data now uses real longer sequences:

- reads: `demo/murugan_naive1_noncoding_demo_seqs.txt`
- germlines: `demo/genomicVs_with_primers.fasta`
- substitution matrix: `models/heavy_pen_substitution_matrix.csv`

The benchmark reports:

- input sequence count;
- repeat factor;
- total benchmark sequence count;
- number of germline templates;
- average read length;
- average germline length;
- legacy runtime and alignment count;
- SeqAn2 banded runtime, alignment count, and band bounds;
- SeqAn2 unbanded runtime and alignment count;
- best-score summaries and SeqAn2-vs-legacy score delta.

## Remaining caveats

The structural implementation is complete, but the benchmark output currently shows a large score discrepancy between legacy `Aligner` and SeqAn2:

```text
legacy_best_score = 274
seqan2_unbanded_best_score = 8
```

This means the benchmark is useful for exercising the code paths and measuring wall-time behavior, but the scoring/coordinate equivalence between the legacy custom aligner and SeqAn2 still needs scientific validation before SeqAn2 can be considered a drop-in replacement.

The likely causes to investigate are:

1. differences in semi-global boundary conditions between legacy `Aligner` and the current benchmark preset;
2. differences in gap-penalty sign/convention;
3. differences between IGoR's custom modified Smith-Waterman scoring and SeqAn2's global/local alignment modes;
4. read/germline orientation or primer-region assumptions in the benchmark setup.

The SIMD two-pass acceleration path remains a future optimization. The CMake option and compile flag wiring exist, but score-only SIMD filtering plus traceback reconstruction has not yet been implemented.

## Git status

The implementation was committed and pushed on:

```text
feature/CIGAR
```

Representative commits include:

```text
Implement CIGAR conversion and SeqAn2 alignment backend
Add unbanded SeqAn2 benchmark mode
Scale SeqAn2 benchmark workload
Complete SeqAn2 fidelity and benchmark metrics
Use real germline data in SeqAn2 benchmark
```

Remote branch:

```text
origin/feature/CIGAR
```

## Final conclusion

The repository now contains a working, tested implementation of the CIGAR conversion infrastructure, protocol-aware alignment presets, optional SeqAn2 backend integration, Waterman-Eggert local alignment enumeration, band calibration, and an improved benchmark harness using real germline data.

The core code paths are verified by unit tests and the feature branch has been pushed for review. The remaining work is scientific validation of SeqAn2 scoring/boundary equivalence against IGoR's legacy aligner, plus optional SIMD acceleration.
