# Plan: CIGAR Conversion + SeqAn2 Banded Alignment Integration (rev7)

## Implementation Task List
- [x] Phase 0: Fix `read_alignments_seq_csv` to load/reconstruct `align_length`, `five_p_offset`, and `three_p_offset`.
- [x] Phase 1: Add CIGAR parsing/conversion/accessor free functions and CIGAR unit tests.
- [x] Phase 2: Add optional SeqAn2 build flags/dependency wiring, `SeqAn2Aligner` backend, and guarded tests.
- [x] Phase 3: Add `AlignmentPreset` factories, calibration APIs, and tests.
- [x] Phase 4: Add benchmark harness target and JSON timing output.
- [x] Verification: Run `./build-seqan2/bin/igor_tests "[cigar],[preset],[seqan2],[calibration]"` and `bench_seqan2_vs_igor`.

## TL;DR
Four sequential phases: (0) fix CSV reader, (1) CIGAR conversion free functions, (2) `SeqAn2Aligner` basic backend with inline config struct, (3) `AlignmentPreset` technology-aware factory with automatic band calibration (Step 3.5). Phase 4 benchmarks SeqAn2 against legacy `Aligner`. SeqAn2 from bioconda (`seqan` v2.5.3), header-only.

---

## Background and Glossary

This section is written for software engineers who are not familiar with immunoinformatics or the IGoR codebase.

### What is IGoR?
IGoR (Inference and Generation Of Repertoires) is a C++ tool that models the stochastic process by which immune cells generate unique antibody and T-cell receptor gene sequences. The core computation is a probabilistic inference over the alignments between an observed sequencing read and a set of *germline gene templates*. Before any probabilistic inference can happen, IGoR must first identify which parts of the observed read correspond to which germline templates — this alignment step is what the `Aligner` class does.

### V, D, and J Genes
The immune receptor locus is divided into three types of gene segments:
- **V (Variable) gene**: encodes the 5' (left) portion of the receptor. The read is expected to start at a known primer site in the V gene.
- **D (Diversity) gene**: a short (~10–30 nt) segment in the middle of the receptor, present only in heavy-chain loci (IgH, TRB, TRD). It is flanked by random non-templated nucleotides and can appear at varying positions within the read.
- **J (Joining) gene**: encodes the 3' (right) portion of the receptor.

The alignment problem for each gene class has different boundary conditions (i.e., which ends of the read/germline can be left unmatched). This directly drives the alignment configuration.

### Smith-Waterman (SW) Local Alignment
Smith-Waterman (1981) is a dynamic programming algorithm that finds the highest-scoring contiguous substring alignment between two sequences. It fills an $(m+1) \times (n+1)$ scoring matrix where each cell $(i, j)$ holds the best alignment score ending at position $i$ of the query and $j$ of the reference, flooring at zero so alignments restart rather than extend with negative scores. The best local alignment is found by backtracking from the cell with the maximum value. Time and space complexity: $O(mn)$.

IGoR's `Aligner` is documented as a *modified* Smith-Waterman: for V and J genes the DP matrix boundary conditions are altered to implement a *semi-global* alignment (some ends are forced to align) while keeping local behavior at other ends. The `local_align` flag in `Aligner` controls whether the zero-floor is active (local) or disabled (global/semi-global).

### Waterman-Eggert Multi-hit Alignment
The **Waterman-Eggert** algorithm (1987) is an extension of Smith-Waterman that finds the $k$ best **non-overlapping** local alignments in decreasing score order. The procedure is:
1. Run SW to find the best alignment.
2. Zero out all DP matrix cells that were on the traceback path of that alignment.
3. Re-run SW (which now ignores the used cells) to find the next best alignment.
4. Repeat until no more alignments exceed the score threshold.

This is the correct approach for **D-gene alignment**, where the D segment may appear multiple times (or multiple candidate D templates may match at different positions within the CDR3 region of the read), and each candidate is a genuinely non-overlapping short alignment.

### Does the Current IGoR Aligner Implement Waterman-Eggert?

**No. The current `Aligner::sw_align` does not implement Waterman-Eggert.**

The IGoR class comment says "modified Smith-Waterman" and the code is indeed a custom multi-hit strategy, but it works by a fundamentally different mechanism than Waterman-Eggert:

**IGoR's actual mechanism — single-pass connected-region tracking:**

The DP matrix is filled exactly once. A secondary integer matrix called `alignment_numb_tracker` is filled in parallel with the score matrix. Each cell in `alignment_numb_tracker` records which "alignment group" owns that cell:
- When a match/mismatch cell $(i, j)$ is computed and its diagonal predecessor $(i-1, j-1)$ has tracker value `-1` (meaning that predecessor was never claimed by any group, i.e., it has score 0), a **new group ID** is created and assigned to $(i, j)$.
- Otherwise $(i, j)$ **inherits** the group ID of its predecessor (diagonal, up, or left, depending on which direction the traceback pointer follows).
- For each group, the tracker records the cell $(r, c)$ that achieved the maximum score seen so far within that group.

After the matrix is fully filled, all groups whose maximum score exceeds `score_threshold` are traced back from their recorded maximum-score cell.

**Key differences from Waterman-Eggert:**

| Property | Waterman-Eggert | IGoR's tracker |
|---|---|---|
| Number of DP passes | $k$ passes (one per alignment) | 1 pass |
| Cell masking | Yes — cells from prior alignments are zeroed | Never |
| Alignment independence | Guaranteed non-overlapping in query coordinates | **Not guaranteed** — groups can cover overlapping query ranges |
| Score ordering | Strict descending order | Within a group: max score; between groups: arbitrary |
| Missing alignments | None (exhaustive) | May miss alignments that share cells with a dominant group |
| Time complexity | $O(k \cdot mn)$ | $O(mn)$ |

**Practical consequence for D-gene alignment**: because IGoR never masks previously used cells, two candidate D alignments that partially share DP cells may be merged into a single group rather than reported as two separate hits. Conversely, a high-scoring D segment at one position may "shadow" a second valid D hit at a different position within the same group. SeqAn2's `LocalAlignmentEnumerator` implements true Waterman-Eggert and should therefore return a strictly larger or equal number of non-overlapping sub-threshold hits.

**Note on the `best_only` flag**: when `best_only=true`, IGoR post-filters and returns only the groups whose maximum score equals the global maximum across all groups. This is independent of the multi-hit mechanism.

---

## Why SeqAn2
- Provides `LocalAlignmentEnumerator<Score, Banded>` — a proper **Waterman-Eggert** implementation with banding for $O(k \cdot (2b+1) \cdot n)$ complexity per sequence pair, where $b$ is the band half-width
- The `Iupac` alphabet (16-letter) directly covers IGoR's 15-letter IUPAC substitution matrix without re-encoding
- SeqAn3 (the successor library) lacks Waterman-Eggert as of version 3.4.x (it is on the roadmap but unimplemented)
- `<seqan/align_parallel.h>` provides SIMD-vectorized batch alignment for high-throughput score filtering

---

## Coordinate System and `Alignment_data` Fields

Understanding the coordinate system is essential before implementing CIGAR conversion. All positions below are 0-indexed unless noted.

`Alignment_data` is a struct in `src/igor/Core/Aligner.h` that holds a single alignment result. Its fields:

| Field | Type | Meaning |
|---|---|---|
| `gene_name` | `string` | Name of the germline template |
| `offset` | `int` | Target position where **position 0 of the full germline template** would align. Formula: `target_pos = germline_pos + offset`. Can be negative if the germline starts before the read. |
| `five_p_offset` | `size_t` | 0-indexed **target** position of the first aligned nucleotide (the 5' end of the aligned region) |
| `three_p_offset` | `size_t` | 0-indexed **target** position of the last aligned nucleotide (the 3' end of the aligned region) |
| `insertions` | `forward_list<int>` | **Target** positions where the germline has a gap (target has an extra nucleotide → CIGAR `I`) |
| `deletions` | `forward_list<int>` | **Germline** positions where the target has a gap (germline has an extra nucleotide → CIGAR `D`) |
| `mismatches` | `mutable vector<int>` | **Target** positions of mismatches (→ CIGAR `X`) |
| `align_length` | `size_t` | Total number of traceback steps (matches + mismatches + insertions + deletions) |
| `score` | `double` | Raw alignment score |

Coordinate relationship:
$$\text{germline\_pos} = \text{target\_pos} - \text{offset}$$
$$\text{five\_p\_offset} = \max(0,\ \text{offset})$$

---

## CIGAR Format

A CIGAR string (Compact Idiosyncratic Gapped Alignment Report) encodes an alignment as a run-length encoded sequence of operations. Example: `3=1X2I4D` means 3 matches, 1 mismatch, 2 insertions (in the target), 4 deletions (in the reference). Operations:

| Op | Meaning | Advances target | Advances germline |
|---|---|---|---|
| `=` | Sequence match | yes | yes |
| `X` | Sequence mismatch | yes | yes |
| `I` | Insertion in target (gap in germline) | yes | no |
| `D` | Deletion from target (gap in target, present in germline) | no | yes |
| `M` | Match or mismatch (legacy, no distinction) | yes | yes |

CIGAR is the standard interchange format used by SAM/BAM files, AIRR-seq formats, and external alignment libraries including SeqAn. IGoR currently stores insertions/deletions/mismatches in separate lists; CIGAR conversion bridges IGoR's internal representation to the outside world.

---

## Phase 0: Fix CSV Reader

### Step 0.1 — Fix `read_alignments_seq_csv` in `Aligner.cpp`

**Context**: When alignments are saved to disk and then reloaded, IGoR uses a semicolon-delimited CSV format. The CSV *writer* (`write_alignments_seq_csv`) already serializes all ten fields including `align_length` (column 8), `five_p_offset` (column 9), and `three_p_offset` (column 10). However, the CSV *reader* (`read_alignments_seq_csv`) has a comment `//FIXME read alignment length` and silently skips these three columns, leaving the fields at their zero-initialized defaults.

This bug has no impact on IGoR's core probabilistic inference because it does not use `five_p_offset`/`three_p_offset` during inference — it recomputes what it needs from `offset`. But it breaks the new CIGAR conversion functions (Phase 1), which require valid `five_p_offset` and `three_p_offset` to walk the alignment region correctly.

**What to change**: In `read_alignments_seq_csv` in `Aligner.cpp`, after parsing the existing 7 semicolon-delimited fields, also parse:
- Column 8 → `align_length`
- Column 9 → `five_p_offset`
- Column 10 → `three_p_offset`

**Backward compatibility**: Old CSV files written by a version of IGoR that did not output these three columns have only 8 fields per row (columns 0–7). Detect this by checking whether columns 9–10 are present after parsing. If absent, reconstruct:
```
five_p_offset  = max(0, offset)
three_p_offset = five_p_offset + align_length - 1 - count(deletions)
```
If `align_length` is also absent (7-column format), set `five_p_offset = 0` as a safe fallback; CIGAR conversion will produce incorrect results on such data but will not crash.

---

## Phase 1: CIGAR Conversion Free Functions

**Why this phase exists**: SeqAn2 (and most external alignment libraries) express alignment results as CIGAR strings paired with 1-based sequence/reference coordinates (the AIRR-seq standard). IGoR internally represents alignments in its own struct (`Alignment_data`). These two functions bridge the gap, enabling:
- Export of IGoR alignments to AIRR/SAM format.
- Import of alignments produced by SeqAn2 into IGoR's data structures for use in inference.

All functions are declared as `CORE_EXPORT` free functions in `Aligner.h` (not methods of any class) because they operate purely on `Alignment_data` values with no dependency on aligner configuration.

### Step 1.1 — `parse_cigar`
```cpp
CORE_EXPORT std::vector<std::pair<int,char>> parse_cigar(const std::string& cigar);
```
Parses a CIGAR string such as `"3=1X2I4D"` into a list of `(count, op_char)` pairs. Must:
- Handle extended CIGAR ops: `=`, `X`, `I`, `D`, `N`, `S`, `H`, `P`
- Accept legacy `M` (match-or-mismatch) as a valid op — needed to import alignments from tools that do not emit extended CIGAR
- Throw `std::invalid_argument` for malformed input (non-digit followed by non-op character, count of zero, etc.)

### Step 1.2 — `alignment_data_to_cigar`
```cpp
CORE_EXPORT std::string alignment_data_to_cigar(const Alignment_data& aln);
```
Converts an `Alignment_data` to an extended CIGAR string. Requires `five_p_offset` and `three_p_offset` to be valid (see Phase 0).

Algorithm:
1. Initialize two cursors: `t = five_p_offset` (target position), `g = five_p_offset - offset` (germline position).
2. Sort `aln.insertions` and `aln.deletions` into sorted containers (they may arrive in arbitrary order from the IGoR traceback).
3. Build an `unordered_set` from `aln.mismatches` for O(1) lookup.
4. Step from `t = five_p_offset` to `t = three_p_offset` (inclusive), emitting one op per step:
   - If `t` is in the insertions set → emit `I`, advance `t` only.
   - If `g` is in the deletions set → emit `D`, advance `g` only.
   - If `t` is in the mismatches set → emit `X`, advance both.
   - Otherwise → emit `=`, advance both.
5. RLE-compress consecutive identical ops into `count + op` format.

### Step 1.3 — `alignment_data_from_cigar`
```cpp
CORE_EXPORT Alignment_data alignment_data_from_cigar(
    const std::string& gene_name,
    const std::string& cigar,
    int seq_start_1based,   // 1-based first aligned position in the target read
    int seq_end_1based,     // 1-based last  aligned position in the target read
    int ref_start_1based,   // 1-based first aligned position in the germline
    int ref_end_1based,     // 1-based last  aligned position in the germline
    double score);
```
Constructs an `Alignment_data` from AIRR-style 1-based coordinates and a CIGAR string. The coordinate mapping:
```
offset        = seq_start_1based - ref_start_1based   (may be negative)
five_p_offset = seq_start_1based - 1                  (convert to 0-based target)
three_p_offset= seq_end_1based   - 1
```
Walk CIGAR ops to populate `insertions`, `deletions`, `mismatches`:
- `=`: advance `t` and `g`; no record
- `X`: record `t` in `mismatches`; advance `t` and `g`
- `I`: record `t` in `insertions`; advance `t` only
- `D`: record `g` in `deletions`; advance `g` only
- `M`: advance `t` and `g`; **no mismatch info** (positions are not added to `mismatches`)
- `align_length = sum of all op counts`

### Step 1.4 — AIRR accessor helpers
Thin inline free functions that translate `Alignment_data` fields to AIRR-seq 1-based column values. Declared in `Aligner.h` as `CORE_EXPORT`:
```cpp
CORE_EXPORT int alignment_data_sequence_start (const Alignment_data&); // five_p_offset + 1
CORE_EXPORT int alignment_data_sequence_end   (const Alignment_data&); // three_p_offset + 1
CORE_EXPORT int alignment_data_germline_start (const Alignment_data&); // five_p_offset - offset + 1
CORE_EXPORT int alignment_data_germline_end   (const Alignment_data&); // three_p_offset - offset + 1
```

### Step 1.5 — Tests: `tst/test_alignment_cigar.cpp`
New Catch2 test file. The build system's `file(GLOB TEST_SOURCES)` will pick it up automatically.

Test cases:
- **Pure match**: 10-nt perfect alignment → CIGAR `10=`
- **Pure mismatch**: 10-nt alignment with all mismatches → CIGAR `10X`
- **Insertion**: 3-nt target insertion at a known position → CIGAR `5=3I5=`
- **Deletion**: 2-nt germline deletion → CIGAR `4=2D6=`
- **Mixed**: combination of matches, mismatches, insertions, deletions
- **Round-trip**: `alignment_data_from_cigar(alignment_data_to_cigar(aln))` must reproduce identical `offset`, `five_p_offset`, `three_p_offset`, `insertions`, `deletions`, `mismatches`, `align_length` (note: round-tripping through `M` ops loses mismatch positions by design)
- **Negative offset**: V-gene alignment where the germline starts before the read (`offset < 0`)
- **`M` fallback**: `parse_cigar("5M")` succeeds; `from_cigar` with `M` ops produces valid `Alignment_data` with empty `mismatches`
- **CSV-loaded round-trip**: construct an `Alignment_data`, write it to CSV, read it back (exercises Phase 0), then call `to_cigar` and verify the result matches the pre-CSV value. This is gated on Phase 0 being done first.

---

## Phase 2: SeqAn2 Basic Aligner Backend

**What this phase builds**: A new class `SeqAn2Aligner` that provides the same public API as `Aligner` (`align_seq` / `align_seqs` returning `forward_list<Alignment_data>`) but uses SeqAn2 under the hood. It is entirely `#ifdef IGOR_WITH_SEQAN2`-guarded so that existing IGoR builds without SeqAn2 are unaffected.

**Banded alignment**: In banded Smith-Waterman/Needleman-Wunsch, only cells within a diagonal band of width $2b+1$ around the main diagonal are computed. This reduces complexity from $O(mn)$ to $O((2b+1) \cdot n)$ per sequence pair. The band is expressed as `[band_lower_diag, band_upper_diag]` — the range of diagonals to compute. For sequences of similar length, a narrow band (e.g., ±20) is sufficient to capture alignments with few indels while dramatically reducing runtime.

**Semi-global alignment and `AlignConfig`**: For V and J genes, the alignment is *semi-global* — some ends are free (no penalty for overhangs) while others are pinned (must be fully aligned). SeqAn2's `AlignConfig<bool,bool,bool,bool>` template controls which ends are free. A `true` value for an end means that end is free (no gap-open penalty for overhanging). See Phase 3 for per-gene-class settings.

**Package management**: IGoR uses `pixi` (a conda-based package manager using `conda-forge` channels) for reproducible builds. SeqAn2 is available as `bioconda::seqan` v2.5.3 (header-only, no compiled library). A new pixi *feature* (`seqan2`) is added so that the default pixi environment is not affected; only environments that include the `seqan2` feature will pull in SeqAn headers.

### Step 2.1 — `pixi.toml` SeqAn2 feature
```toml
[feature.seqan2.dependencies]
seqan = { version = ">=2.4.0", channel = "bioconda" }

[environments]
seqan2 = { features = ["seqan2"], solve-group = "default" }
```
`bioconda::seqan` is `noarch` (header-only, latest 2.5.3).

### Step 2.2 — Root `CMakeLists.txt`
```cmake
option(IGOR_WITH_SEQAN2 "Enable SeqAn2 banded alignment backend" OFF)
if(IGOR_WITH_SEQAN2)
  # FindSeqAn.cmake ships inside the bioconda seqan package.
  # When building inside a pixi/conda environment, CMAKE_PREFIX_PATH already
  # points to the conda prefix, so find_package should locate it automatically.
  # If not, set CMAKE_MODULE_PATH to ${CONDA_PREFIX}/share/seqan/cmake.
  find_package(SeqAn REQUIRED)
endif()
```
`SEQAN_DEFINITIONS` (provided by `FindSeqAn.cmake`) contains macro definitions that SeqAn headers require (e.g., to enable ZLIB / BZ2 / OpenMP support). These must be passed to `target_compile_definitions`. OpenMP is already linked `PRIVATE` to the `Core` target by the existing build system.

### Step 2.3 — `src/igor/Core/CMakeLists.txt`
```cmake
if(IGOR_WITH_SEQAN2)
  target_sources(Core PRIVATE SeqAn2Aligner.cpp)
  target_include_directories(Core PRIVATE ${SEQAN_INCLUDE_DIRS})
  # IGOR_WITH_SEQAN2 is PUBLIC so downstream targets (app, tests) also see it
  target_compile_definitions(Core PUBLIC IGOR_WITH_SEQAN2 ${SEQAN_DEFINITIONS})
  target_link_libraries(Core PRIVATE ${SEQAN_LIBRARIES})

  option(IGOR_SEQAN2_SIMD "Compile SeqAn2 with SIMD vectorization (requires -msse4 or -mavx2)" OFF)
  if(IGOR_SEQAN2_SIMD)
    # Change -msse4 to -mavx2 for AVX2 support on capable hardware
    target_compile_options(Core PRIVATE -msse4)
  endif()
endif()
```

### Step 2.4 — `SeqAn2AlignConfig` inline struct (in `SeqAn2Aligner.h`, no separate file)
```cpp
struct SeqAn2AlignConfig {
    int  band_lower_diag    = -20;  // most negative diagonal to compute
    int  band_upper_diag    =  20;  // most positive diagonal to compute
    bool seq2_leading_free  = false;  // germline 5' end: free overhang allowed?
    bool seq1_leading_free  = false;  // read     5' end: free overhang allowed?
    bool seq1_trailing_free = true;   // read     3' end: free overhang allowed?
    bool seq2_trailing_free = true;   // germline 3' end: free overhang allowed?
    bool local_alignment    = false;  // true = use LocalAlignmentEnumerator (D gene)
};
```

**SeqAn2 `AlignConfig` semantics** (confirmed from SeqAn2 2.5.3 source and docs):
`AlignConfig<seq2_leading, seq1_leading, seq1_trailing, seq2_trailing>`
where `seq1 = row(ali,0)` = the first sequence placed in the alignment object = the **read** (horizontal axis), `seq2 = row(ali,1)` = the **germline** (vertical axis).
A `true` argument makes that end free (no penalty for unaligned overhangs at that end).

This struct is a private implementation detail of `SeqAn2Aligner`. It will be replaced by `AlignmentPreset` in Phase 3; the fields map 1:1.

### Step 2.5 — `SeqAn2Aligner.h` (guarded by `#ifdef IGOR_WITH_SEQAN2`)

The header must not expose any SeqAn2 types in its public interface (they are too heavy to include transitively). All SeqAn2 types are used only in `SeqAn2Aligner.cpp`.

```cpp
#ifdef IGOR_WITH_SEQAN2
#pragma once
#include "Aligner.h"          // Alignment_data, Gene_class
#include "igor_utils.h"       // Matrix<double>

struct SeqAn2AlignConfig { /* ... see Step 2.4 */ };

class CORE_EXPORT SeqAn2Aligner {
public:
    // subst: 15×15 IUPAC substitution score matrix (same layout as Aligner)
    // gap_penalty: linear gap penalty (subtracted for each gap column)
    // gene: V_gene, D_gene, or J_gene — used to validate config at construction
    SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene,
                  SeqAn2AlignConfig config = {});
    ~SeqAn2Aligner();  // defined in .cpp to avoid incomplete-type errors

    // Align one read against all stored genomic sequences.
    // Returns alignments above threshold; if best_only=true, returns only
    // the highest-scoring alignment per genomic template.
    std::forward_list<Alignment_data>
    align_seq(const std::string& seq, double threshold,
              bool allow_in_dels, bool best_only = false) const;

    // Align multiple reads in parallel (OpenMP).
    // seqs: vector of (sequence_index, sequence_string) pairs
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(const std::vector<std::pair<int, std::string>>& seqs,
               double threshold, bool allow_in_dels) const;

    // Supply the germline templates to align against.
    // Each pair is (gene_name, nucleotide_sequence).
    void set_genomic_sequences(std::vector<std::pair<std::string, std::string>>);

private:
    struct Impl;           // pimpl — hides SeqAn2 types from the header
    std::unique_ptr<Impl> impl_;
};
#endif // IGOR_WITH_SEQAN2
```

The pimpl (`Impl`) pattern keeps SeqAn2 types (e.g., `Score<int, ScoreMatrix<Iupac, Default>>`) out of the header. Downstream translation units that include `SeqAn2Aligner.h` do not need to include any SeqAn2 headers themselves.

### Step 2.6 — `SeqAn2Aligner.cpp` implementation

**Alphabet conversion**: IGoR internally encodes nucleotides as small integers using the `Nucleotides` enum (0–14 for A, C, G, T, and IUPAC ambiguity codes). SeqAn2's `Iupac` alphabet covers 16 letters (0–15). At construction, build a `char igor_to_iupac_char[15]` lookup table that maps each IGoR `Nucleotides` integer to the corresponding `Iupac` character, using `Nucl_arr` (IGoR's nucleotide-to-char table). This allows direct use of `Iupac(c)` with no per-alignment overhead.

**Scoring scheme construction**:
```cpp
typedef seqan2::Score<int, seqan2::ScoreMatrix<seqan2::Iupac, seqan2::Default>> TScore;
TScore scoring_scheme(gap_penalty, gap_penalty);  // (gap-extend, gap-open) for linear gaps
for (int i = 0; i < 15; ++i)
    for (int j = 0; j < 15; ++j)
        seqan2::setScore(scoring_scheme,
                         seqan2::Iupac(igor_to_iupac_char[i]),
                         seqan2::Iupac(igor_to_iupac_char[j]),
                         static_cast<int>(subst(i, j)));
```
Note: SeqAn2 version 2.5.3 uses the `seqan2::` namespace (older versions used `seqan::`).

**V/J semi-global banded alignment** (one genomic template at a time):
```cpp
using TSeq  = seqan2::String<seqan2::Iupac>;
using TAlign = seqan2::Align<TSeq, seqan2::ArrayGaps>;
TAlign ali;
seqan2::resize(seqan2::rows(ali), 2);
seqan2::assignSource(seqan2::row(ali, 0), read_iupac);    // query = read
seqan2::assignSource(seqan2::row(ali, 1), germline_iupac); // subject = germline
int sc = seqan2::globalAlignment(
    ali,
    scoring_scheme_,
    make_align_config(config_),      // converts SeqAn2AlignConfig → AlignConfig<>
    config_.band_lower_diag,
    config_.band_upper_diag);
if (sc < threshold) continue;
results.push_front(from_seqan2_align(ali, gene_name, sc));
```
`make_align_config` is a small helper that maps the four `bool` fields of `SeqAn2AlignConfig` to the corresponding `seqan2::AlignConfig<bool,bool,bool,bool>` template instantiation using `if constexpr` or a runtime dispatch table.

**D-gene multi-hit via Waterman-Eggert** (true Waterman-Eggert in SeqAn2):
```cpp
seqan2::LocalAlignmentEnumerator<TScore, seqan2::Banded> enumerator(
    scoring_scheme_,
    config_.band_lower_diag,
    config_.band_upper_diag,
    static_cast<int>(threshold));   // minimum score — alignments below this are not returned
TAlign ali;
seqan2::resize(seqan2::rows(ali), 2);
seqan2::assignSource(seqan2::row(ali, 0), read_iupac);
seqan2::assignSource(seqan2::row(ali, 1), germline_iupac);
while (seqan2::nextLocalAlignment(ali, enumerator)) {
    if (!allow_in_dels && has_indels(ali)) continue;
    results.push_front(from_seqan2_align(ali, gene_name,
                                          seqan2::getScore(enumerator)));
    if (best_only) break;
}
```
Each call to `nextLocalAlignment` masks the cells used by the previous alignment (true Waterman-Eggert) and returns the next best non-overlapping hit. Hits are returned in strictly decreasing score order.

`has_indels(ali)` inspects `row(ali,0)` and `row(ali,1)` for gap characters using `seqan2::isGap`.

**CIGAR extraction in `from_seqan2_align`**:
After alignment, SeqAn2 stores the result as a gapped sequence pair inside the `Align` object. To extract IGoR's position lists, walk the clipped aligned region column by column:
```cpp
auto& r0 = seqan2::row(ali, 0);  // read row (with gaps inserted)
auto& r1 = seqan2::row(ali, 1);  // germline row
int t = seqan2::clippedBeginPosition(r0);  // 0-based start in read
int g = seqan2::clippedBeginPosition(r1);  // 0-based start in germline
for (int col = seqan2::clippedBeginPosition(r0);
         col < seqan2::clippedEndPosition(r0);
         ++col) {
    if      (seqan2::isGap(r0, col)) { deletions.push_front(g++); }
    else if (seqan2::isGap(r1, col)) { insertions.push_front(t++); }
    else {
        if (r0[col] != r1[col]) mismatches.push_back(t);
        ++t; ++g;
    }
}
```
Then call `alignment_data_from_cigar()` (Phase 1) with 1-based coordinates derived from the clipped positions. This ensures the coordinate conventions are consistent between IGoR's legacy aligner and SeqAn2.

**Multithreading**: `align_seqs` uses an OpenMP parallel-for loop identical in structure to `Aligner::align_seqs`. Each loop iteration constructs its own `TAlign` object on the stack — SeqAn2 alignment objects are not shared between threads. The `scoring_scheme_` stored in `Impl` is read-only after construction and is therefore thread-safe.

### Step 2.7 — Tests: `tst/test_seqan2_aligner.cpp` (guarded by `#ifdef IGOR_WITH_SEQAN2`)
```cpp
#ifdef IGOR_WITH_SEQAN2
// ... test cases ...
#endif
```
Test cases:
- **V-gene basic**: align a synthetic read with a known offset against a synthetic germline; verify `five_p_offset`, `three_p_offset`, and score match expected values
- **V-gene negative offset**: read starts before the germline template
- **J-gene with/without primer**: supply the same read and germline but switch `AlignConfig` bools to simulate the two protocol variants; verify that the free-end variant produces a higher or equal score when there is a primer mismatch region
- **D-gene single hit**: one short D segment buried in random flanking sequence; verify the enumerator finds it
- **D-gene multiple hits**: two non-overlapping D segment hits above threshold; verify at least two hits are returned in decreasing score order
- **`allow_in_dels=false`**: construct an alignment that requires an indel; verify the result is absent from the output
- **Qualitative vs `Aligner`**: run both aligners on at least one `demo/` sequence; print a side-by-side score comparison; assert scores are within ±5 (the exact tolerance is TBD after initial testing)

---

## Phase 3: AlignmentPreset — Technology-Aware Config

*Depends on Phase 2. Replaces `SeqAn2AlignConfig` as the config source.*

**Why this phase exists**: Phase 2 introduces `SeqAn2AlignConfig` as an internal struct with raw Boolean fields whose meaning is not obvious without reading the SeqAn2 documentation. Phase 3 replaces this with `AlignmentPreset`, a higher-level, named configuration object whose factory methods express intent ("V gene in 5' RACE protocol" or "J gene with J-gene reverse primer") rather than mechanism ("seq1_leading_free=false"). This makes the configuration reusable across different aligner backends and testable in isolation.

`AlignmentPreset` is completely independent of SeqAn2 — it holds only primitive fields — so it can be included in any translation unit without pulling in SeqAn2 headers.

### Sequencing protocols and the AlignConfig mapping

This section explains the two most common immune-receptor sequencing protocols and how each determines which alignment boundaries are **free** vs **anchored**. Understanding this is essential to choose or author the correct `AlignmentPreset`.

#### The free-end concept for software engineers

Each observed read is a raw DNA sequence produced by the sequencer. Its start and end are set by where PCR primers anneal and how many bases the sequencer reads. When IGoR aligns a read against a germline template, only one (read, template) pair is aligned at a time. That alignment has four boundary points:

```
         read 5'                                             read 3'
         |                                                       |
Read:    |== prefix (unaligned) ==|=== matched region ===|= suffix =|
                                  |                       |
                          germline 5'             germline 3'
```

For each boundary:
- **Free** (`true` in `AlignConfig`): unmatched nucleotides on that side carry **no gap penalty**. Equivalent to a SAM soft-clip.
- **Anchored** (`false`): unmatched nucleotides on that side are **penalised as gaps** in the DP score.

`AlignConfig<seq2_leading, seq1_leading, seq1_trailing, seq2_trailing>` where `seq1` = read (horizontal) and `seq2` = germline (vertical).

#### Protocol 1 — 5' RACE (Rapid Amplification of cDNA Ends)

**In one sentence**: a single pair of primers captures all receptor variants from RNA without gene-family-specific sequences. Reads start *before* the V gene and end in the C region.

**Mechanism**: mRNA is reverse-transcribed into cDNA. The reverse transcriptase adds untemplated C nucleotides to the 3' end of the fresh cDNA (the 5' end of the original mRNA). A **Template-Switch Oligonucleotide (TSO)** with a poly-G tail anneals to these C nucleotides and causes the polymerase to switch templates, appending a universal adaptor sequence to the cDNA's 5' end. A UMI (Unique Molecular Identifier — a random sequence) is often embedded in the adaptor to tag individual molecules for PCR deduplication. One forward primer targeting the adaptor and one reverse primer targeting the C region then amplify all receptor variants.

**Resulting read structure**:
```
5'-[TSO adaptor]-[UMI?]-[5'UTR / leader ~50-100 nt]-[V gene]-[CDR3]-[J gene]-[C region]-3'
                                                                                         ^
                                                                         C-region primer
```

**V gene alignment — 5' RACE**:
- Read has adaptor + 5'UTR/leader (~50–100 nt) **before** the V gene starts → `seq1_leading_free = true`.
- V germline is fully covered from position 0 → `seq2_leading_free = false`.
- Read continues into CDR3 after V gene → `seq1_trailing_free = true`.
- V germline ends at V–CDR3 junction → `seq2_trailing_free = false`.
- Expected DP diagonal: `+L` (leader length). Band: `[L - max_indel, L + max_indel]`.

**J gene alignment — 5' RACE** (read covers CDR3 + J gene + C region):
- Read has CDR3 prefix (up to C nt) before J germline → `seq1_leading_free = true`.
- J germline starts at position 0 → `seq2_leading_free = false`.
- Read continues into C gene after J → `seq1_trailing_free = true`.
- J germline ends at J–C junction → `seq2_trailing_free = false`.
- Expected DP diagonal: `0` to `+C`. Band: `[0, C + max_indel]`.

#### Protocol 2 — Multiplex PCR

**In one sentence**: a panel of V-gene-family-specific primers amplifies all variants and works on both RNA and genomic DNA (gDNA), but introduces per-family amplification bias. Reads start **inside** the V gene at the primer binding site.

**Mechanism**: a set of forward primers — one per V gene family — is pooled together. Each primer anneals at a **known position inside** its target V gene (not at V gene position 0, but typically 50–250 nt in). A reverse primer targets either the J gene (gDNA, unprocessed) or the C region (cDNA, spliced). Multiple PCR rounds amplify the pooled product.

**gDNA vs cDNA**: on genomic DNA the VDJ segments are separated by introns, so J-region primers are needed (shorter products, fewer primer variants required). On cDNA the introns are already spliced out, so a single C-region primer covers all J genes.

**V gene alignment — Multiplex PCR**:
- V germline has a prefix (positions 0 to `primer_pos_in_gene`) absent from the read → `seq2_leading_free = true`.
- Read starts at the primer, no read-side prefix → `seq1_leading_free = false`.
- Read continues into CDR3 → `seq1_trailing_free = true`.
- V germline ends at V–CDR3 junction → `seq2_trailing_free = false`.
- Expected DP diagonal: `-P` (primer position). Band: `[-(P + max_indel), -(P - max_indel)]`.

**J gene alignment — C-region primer** (cDNA Multiplex PCR or any 5' RACE-like J setup):
- Identical to the 5' RACE J case: `seq1_leading=true`, `seq2_leading=false`, `seq1_trailing=true`, `seq2_trailing=false`. Band: `[0, C + max_indel]`.

**J gene alignment — J-gene primer** (gDNA Multiplex PCR):
- A J-gene-specific reverse primer anneals **within** the J gene; read ends at the primer site before the J gene's 3' end.
- Read arrives from CDR3 → `seq1_leading_free = true`. J germline starts at position 0 → `seq2_leading_free = false`.
- Read ends at J primer → `seq1_trailing_free = false`. J germline extends past the read end → `seq2_trailing_free = true`.
- Band: `[0, C + max_indel]` (same starting-position range as C-primer case).

#### D gene — both protocols

The D gene is a short (~10–30 nt) segment surrounded by non-templated nucleotides (N/P additions) at an unpredictable position within the CDR3. No semi-global configuration is appropriate — its read position is unknown. The correct approach is **local alignment** (Waterman-Eggert via SeqAn2's `LocalAlignmentEnumerator`). When `local_alignment = true`, the `AlignConfig` booleans are ignored.

#### Band diagonal reference

In the DP matrix, read is horizontal (index `i`) and germline is vertical (index `j`). Diagonal `k = i - j`. Positive k: read starts before germline (e.g. 5' RACE). Negative k: germline starts before read (e.g. Multiplex PCR V gene).

| Protocol + gene | Expected diagonal | Band |
|---|---|---|
| 5' RACE, V gene (leader length `L`) | `+L` | `[L - max_indel, L + max_indel]` |
| Multiplex PCR, V gene (primer at pos `P`) | `-P` | `[-(P+max_indel), -(P-max_indel)]` |
| J gene, C primer (max CDR3 `C`) | `0` to `+C` | `[0, C + max_indel]` |
| J gene, J primer (max CDR3 `C`) | `0` to `+C` | `[0, C + max_indel]` |
| D gene | unpredictable | `[-max_indel, +max_indel]` |

### Step 3.1 — `AlignmentPreset.h` (no SeqAn dependency, standalone public header)
```cpp
/// Technology-agnostic alignment boundary configuration.
/// Maps directly to SeqAn2's AlignConfig and banded-alignment parameters.
/// Does NOT depend on SeqAn2 headers — safe to include anywhere.
struct AlignmentPreset {
    // Banded DP: only cells where
    //   band_lower_diag <= (read_pos - germline_pos) <= band_upper_diag
    // are computed. Narrower band = faster; must cover expected offset +/- max_indel.
    int  band_lower_diag;
    int  band_upper_diag;

    // Maps to SeqAn2: AlignConfig<seq2_leading, seq1_leading, seq1_trailing, seq2_trailing>
    //   seq1 = read     (DP row 0, horizontal axis, index i)
    //   seq2 = germline (DP row 1, vertical axis,   index j)
    //   true  = FREE:     unmatched nucleotides at that side carry no gap penalty
    //   false = ANCHORED: unmatched nucleotides at that side are penalised as gaps
    bool seq2_leading_free;   // germline 5' end: true -> germline has prefix before read starts
    bool seq1_leading_free;   // read     5' end: true -> read    has prefix before germline starts
    bool seq1_trailing_free;  // read     3' end: true -> read continues past germline end
    bool seq2_trailing_free;  // germline 3' end: true -> germline extends past read end

    // If true:  use SeqAn2 LocalAlignmentEnumerator (Waterman-Eggert multi-hit).
    //           The AlignConfig booleans above are ignored.
    // If false: use globalAlignment with the AlignConfig booleans above.
    bool local_alignment;

    // --- Factory methods ---

    // 5' RACE, V gene.
    //   Read = [TSO adaptor + 5'UTR/leader, ~max_leader_nt nt] + [V gene] + [CDR3 ...]
    //   seq1_leading=true  (adaptor+leader prefix in read before germline)
    //   seq2_leading=false (germline starts at position 0)
    //   seq1_trailing=true (read continues past V gene into CDR3)
    //   seq2_trailing=false(germline ends at V-CDR3 junction)
    //   band = [max_leader_nt - max_indel, max_leader_nt + max_indel]
    static AlignmentPreset v_gene_5p_race(int max_leader_nt = 80, int max_indel = 15);

    // Multiplex PCR, V gene.
    //   V primer anneals primer_pos_in_gene nt inside the V gene; read starts there.
    //   seq1_leading=false (read starts at primer, no prefix)
    //   seq2_leading=true  (germline has prefix of length primer_pos_in_gene)
    //   seq1_trailing=true (read continues past V gene into CDR3)
    //   seq2_trailing=false(germline ends at V-CDR3 junction)
    //   band = [-(primer_pos_in_gene + max_indel), -(primer_pos_in_gene - max_indel)]
    static AlignmentPreset v_gene_multiplex_pcr(int primer_pos_in_gene = 0, int max_indel = 15);

    // J gene, C-region reverse primer (5' RACE or Multiplex PCR on cDNA).
    //   Read = [CDR3, up to max_cdr3_nt nt] + [J gene] + [C region up to primer]
    //   seq1_leading=true  (CDR3 prefix before J germline)
    //   seq2_leading=false (J germline starts at position 0)
    //   seq1_trailing=true (read continues into C region past J end)
    //   seq2_trailing=false(J germline ends at J-C junction)
    //   band = [0, max_cdr3_nt + max_indel]
    static AlignmentPreset j_gene_c_primer(int max_cdr3_nt = 60, int max_indel = 5);

    // J gene, J-gene reverse primer (Multiplex PCR on gDNA).
    //   Read = [CDR3, up to max_cdr3_nt nt] + [J gene up to J primer site]
    //   Read ends inside J gene; J germline extends past the read's 3' end.
    //   seq1_leading=true  (CDR3 prefix before J germline)
    //   seq2_leading=false (J germline starts at position 0)
    //   seq1_trailing=false(read ends at J primer; no trailing read sequence)
    //   seq2_trailing=true (J germline extends past the read end)
    //   band = [0, max_cdr3_nt + max_indel]
    static AlignmentPreset j_gene_j_primer(int max_cdr3_nt = 60, int max_indel = 5);

    // D gene (both protocols).
    //   Short (~10-30 nt) segment at unknown position within CDR3.
    //   local_alignment = true -- uses Waterman-Eggert; AlignConfig booleans ignored.
    //   band = [-max_indel, +max_indel]
    static AlignmentPreset d_gene(int max_indel = 3);

    // Symmetric band, all ends anchored -- utility preset for custom use and unit tests.
    static AlignmentPreset with_band(int half_width);
};
```

Preset boundary configuration summary:

| Preset | seq2_leading (germ 5') | seq1_leading (read 5') | seq1_trailing (read 3') | seq2_trailing (germ 3') | local |
|---|:---:|:---:|:---:|:---:|:---:|
| `v_gene_5p_race` | false | **true** | **true** | false | false |
| `v_gene_multiplex_pcr` | **true** | false | **true** | false | false |
| `j_gene_c_primer` | false | **true** | **true** | false | false |
| `j_gene_j_primer` | false | **true** | false | **true** | false |
| `d_gene` | — | — | — | — | **true** |
| `with_band` | false | false | false | false | false |

`v_gene_5p_race` and `j_gene_c_primer` have identical boolean fields because in both cases the read has a prefix (leader or CDR3) before the germline and continues past it. The difference is entirely in the band: a positive offset for the V-gene leader, a positive-from-zero offset for the J-gene CDR3.

Band values:

| Preset | `band_lower_diag` | `band_upper_diag` |
|---|---|---|
| `v_gene_5p_race(L, indel)` | `L - indel` | `L + indel` |
| `v_gene_multiplex_pcr(P, indel)` | `-(P + indel)` | `-(P - indel)` |
| `j_gene_c_primer(C, indel)` | `0` | `C + indel` |
| `j_gene_j_primer(C, indel)` | `0` | `C + indel` |
| `d_gene(indel)` | `-indel` | `+indel` |
| `with_band(w)` | `-w` | `+w` |

### Step 3.2 — Refactor `SeqAn2Aligner` to accept `AlignmentPreset`
Add a constructor overload:
```cpp
SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene,
              AlignmentPreset preset);
```
Since `AlignmentPreset` and `SeqAn2AlignConfig` share the same fields with identical semantics, the constructor body converts `AlignmentPreset` → `SeqAn2AlignConfig` trivially (or replaces `SeqAn2AlignConfig` entirely; prefer the latter to avoid dead code).

### Step 3.3 — Install `AlignmentPreset.h` unconditionally
`AlignmentPreset.h` is added to the `Core` target's unconditional `FILE_SET HEADERS` — it carries no SeqAn2 dependency and is part of the public IGoR API surface.

### Step 3.4 — Tests: `tst/test_alignment_preset.cpp`
- **Band `v_gene_5p_race`**: `v_gene_5p_race(80, 15)` → `band_lower_diag=65`, `band_upper_diag=95`
- **Band `v_gene_multiplex_pcr`**: `v_gene_multiplex_pcr(50, 15)` → `band_lower_diag=-65`, `band_upper_diag=-35`
- **Band `j_gene_c_primer`**: `j_gene_c_primer(60, 5)` → `band_lower_diag=0`, `band_upper_diag=65`
- **Band `j_gene_j_primer`**: `j_gene_j_primer(60, 5)` → `band_lower_diag=0`, `band_upper_diag=65`
- **Zero indel**: `v_gene_5p_race(80, 0)` → `band_lower_diag=80`, `band_upper_diag=80` (zero-width band, only the exact leader-length diagonal is explored)
- **Zero primer offset**: `v_gene_multiplex_pcr(0, 10)` → `band_lower_diag=-10`, `band_upper_diag=10` (primer at gene start, equivalent to `with_band(10)` but with seq2_leading=true)
- **Boolean fields**: verify every row of the preset table above for all six factory methods
- **D-gene band**: `d_gene(3)` → `band_lower_diag=-3`, `band_upper_diag=3`; `d_gene().local_alignment == true`
- **`with_band` symmetry**: `with_band(20)` → `band_lower_diag=-20`, `band_upper_diag=20`; all four free-end bools are `false`

### Step 3.5 — Automatic band calibration via full-matrix warm-up subset

**Why this step exists**: The `AlignmentPreset` factory methods require the caller to supply parameters like `max_leader_nt` (5' RACE) or `primer_pos_in_gene` (Multiplex PCR) that encode domain knowledge about the sequencing protocol. When these parameters are unknown — for example when processing data from an undocumented protocol, or when fine-tuning the band to a specific germline database — requiring the user to guess them by hand undermines the goal of a safe, reproducible alignment. Step 3.5 adds a calibration procedure that infers the correct band from the data itself.

**Algorithm — two-stage bootstrap**:

*Stage 1 — full-matrix seed alignment (unbanded).*
Run standard `globalAlignment` (or `localAlignment` for D gene) without any band restriction on a small random sample of $N_1$ sequences (default 50). For each (read, germline) pair in the sample, record the **optimal alignment diagonal** — the dominant diagonal of the traceback path — defined as:
$$d_i = \text{clippedBeginPosition}(\text{row}(ali, 0)) - \text{clippedBeginPosition}(\text{row}(ali, 1))$$
This is the diagonal index $k = i - j$ at the point where the alignment begins, which equals the read start offset minus the germline start offset. Collect the set $\{d_i\}$ across all sample pairs and compute:
- $\hat{d}$: median of $\{d_i\}$
- $\delta_-$: $\hat{d} - P_5(\{d_i\})$ (downward spread to the 5th percentile)
- $\delta_+$: $P_{95}(\{d_i\}) - \hat{d}$ (upward spread to the 95th percentile)

Initial band estimate:
$$[\hat{d} - \delta_- - b_{\text{init}},\quad \hat{d} + \delta_+ + b_{\text{init}}]$$
where $b_{\text{init}}$ is a caller-supplied minimum slack (default 15 nt, corresponding to `max_indel`).

*Stage 2 — progressive band refinement (optional, iterative).*
Using the Stage 1 band estimate, run *banded* `globalAlignment` on a second independent sample of $N_2$ sequences (default 100). For each pair compare the banded score against the previously recorded full-matrix score:
- If **all** banded scores match the full-matrix scores (within `score_tolerance`, default 0): the band is verified; proceed to Stage 3 with the current band.
- If **any** banded score falls below the full-matrix score: the band clipped the optimal path for that pair. Widen the band symmetrically by $\Delta b$ (default 10 nt) and repeat Stage 2 with the widened band. Stop when all scores agree or after `max_iterations` (default 3) expansions.

*Stage 3 — commit.*
Update the `AlignmentPreset`'s `band_lower_diag` and `band_upper_diag` to the verified values. The calibrated preset is ready for use in `SeqAn2Aligner`.

**Diagonal computation note**: $d_i$ as defined above equals the `offset` field that IGoR would record in `Alignment_data` for the same alignment start, so if Stage 1 is run after a legacy `Aligner` pass the offsets can be recycled directly without a second full-matrix run.

**New API on `AlignmentPreset`**:
```cpp
struct BandCalibrationParams {
    int  n_seed_seqs      = 50;    // Stage 1: number of sequences for full-matrix seed
    int  n_verify_seqs    = 100;   // Stage 2: number of sequences for banded verification
    int  min_band_slack   = 15;    // minimum slack added to percentile spread (nt)
    int  band_expand_step = 10;    // per-iteration band widening in Stage 2 (nt)
    int  max_iterations   = 3;     // maximum Stage 2 expansion iterations
    int  score_tolerance  = 0;     // allowed score loss vs full-matrix (0 = exact match)
    int  seed             = 42;    // RNG seed for reproducible sampling
};

/// Calibrate the band of this preset by running full-matrix alignment on a
/// random sample of (read, germline) pairs, then verifying with banded alignment.
/// The preset's band_lower_diag and band_upper_diag are updated in-place.
/// Requires IGOR_WITH_SEQAN2 at compile time; no-op otherwise.
///
/// @param reads          pool of candidate read sequences to sample from
/// @param germlines      all germline templates to sample against (names ignored)
/// @param scoring        scoring scheme to use for calibration alignments
/// @param params         tuning parameters (defaults are reasonable)
/// @return               reference to *this (for method chaining)
AlignmentPreset& calibrate(
    const std::vector<std::string>& reads,
    const std::vector<std::string>& germlines,
    const seqan2::Score<int>& scoring,  // forward-declared; defined only if IGOR_WITH_SEQAN2
    BandCalibrationParams params = {});
```

Because `calibrate()` depends on SeqAn2 types in its signature (the `Score` parameter), it is declared inside an `#ifdef IGOR_WITH_SEQAN2` block in `AlignmentPreset.h`. This is the single exception to the rule that `AlignmentPreset.h` is SeqAn2-free; all other fields and factory methods remain unconditional.

**Alternative signature without SeqAn2 in the header** (preferred if `AlignmentPreset.h` must stay fully unconditional): move `calibrate` to a free function in `SeqAn2Aligner.h`:
```cpp
#ifdef IGOR_WITH_SEQAN2
CORE_EXPORT AlignmentPreset calibrate_preset(
    AlignmentPreset preset,
    const std::vector<std::string>& reads,
    const std::vector<std::string>& germlines,
    const seqan2::Score<int>& scoring,
    BandCalibrationParams params = {});
#endif
```
This keeps `AlignmentPreset.h` fully SeqAn2-free. **Prefer this form.** `BandCalibrationParams` lives in `AlignmentPreset.h` (unconditional, no SeqAn2 dependency).

**Integration with `SeqAn2Aligner`**: add an optional `auto_calibrate` flag to the `SeqAn2Aligner` constructor:
```cpp
// If calibration_reads is non-empty, runs calibrate_preset() automatically
// before the first align_seqs call, using the first min(n_seed_seqs, reads.size())
// sequences from calibration_reads as Stage 1 input.
SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene,
              AlignmentPreset preset,
              std::vector<std::string> calibration_reads = {},
              BandCalibrationParams calib_params = {});
```
When `calibration_reads` is non-empty the constructor calls `calibrate_preset` and stores the refined preset before any `align_seq` call. The germline sequences come from those already loaded via `set_genomic_sequences`. The calibration is run at construction time (not lazily) so that the preset is fixed and reproducible across all subsequent calls.

**Performance note**: Stage 1 runs at most $N_1 \times |\text{germlines}|$ full-matrix alignments. For typical IGoR V-gene databases (~60 templates) with $N_1 = 50$, this is 3 000 unbanded DP passes, each $O(mn)$ with $m \approx 300$, $n \approx 300$: approximately $270 \times 10^6$ cell evaluations, well under 1 second on a modern core. Stage 2 adds at most $3 \times 100 \times 60 = 18 000$ banded passes, which are faster by a factor of $mn / ((2b+1)n) \approx 300 / 41 \approx 7$. Total calibration overhead is expected to be $\lesssim 2$ seconds regardless of dataset size.

### Step 3.6 — Tests: `tst/test_band_calibration.cpp` (guarded by `#ifdef IGOR_WITH_SEQAN2`)

- **Synthetic 5' RACE calibration**: generate 200 synthetic reads with a random leader of length drawn from $\mathcal{U}[60, 100]$ nt; calibrate with 50 seed + 100 verify; assert the resulting band contains the true interval $[60 - 15, 100 + 15] = [45, 115]$ and has width $\leq 80$ nt.
- **Synthetic Multiplex PCR calibration**: same construction with a primer at fixed offset 70 nt inside the germline; assert band is centred near $-70$ with width $\leq 40$.
- **Stage 2 expansion triggered**: construct a dataset where Stage 1 gives a deliberately narrow estimate; assert Stage 2 expands the band until all banded scores equal the full-matrix scores.
- **`score_tolerance > 0`**: with `score_tolerance = 2`, assert the calibration accepts a band that misses a few alignments by up to 2 points, producing a narrower band than `score_tolerance = 0`.
- **Reproducibility**: run calibration twice with the same `seed`; assert identical `band_lower_diag` and `band_upper_diag`.
- **Empty calibration reads**: if `calibration_reads` is empty, `SeqAn2Aligner` must use the preset band unchanged without crashing.
- **Calibrated preset feeds `SeqAn2Aligner`**: run `SeqAn2Aligner` with a calibrated preset on the full synthetic dataset; assert ≥99% of returned best scores match those from the unbanded `SeqAn2Aligner` on the same data.

---

## Phase 4: Performance and Results Comparison

*Depends on Phase 2+3. Independent of Phase 1.*

**Goal**: Establish quantitative evidence about (a) how closely SeqAn2's scores and traceback positions match IGoR's existing aligner, and (b) how much faster SeqAn2 can run, both with OpenMP and with SIMD vectorization.

### Step 4.1 — Benchmark harness
The harness runs both `Aligner` and `SeqAn2Aligner` on the same input sequence set and collects structured JSON/CSV results that can be processed offline.

Suggested implementation: a standalone binary `tst/bench_seqan2_vs_igor.cpp` compiled as a separate CMake target (not part of `igor_tests`). This avoids polluting the unit-test run time with benchmark execution.

Inputs:
- `demo/murugan_naive1_noncoding_demo_seqs.txt` (built-in, always available)
- Optionally a larger external AIRR dataset (passed as CLI argument) for the scaling tests

Metrics collected per run:
- Wall time (using `std::chrono::high_resolution_clock`) for alignment of the entire sequence set
- Per-sequence: number of alignments above threshold, best score, gene name of best hit
- Memory usage (peak RSS, if available on the platform)

### Step 4.2 — Results correctness comparison

Because IGoR's multi-hit mechanism differs from Waterman-Eggert (see "Background" section), the two aligners are **not expected to produce identical results**. The comparison verifies that SeqAn2 is a valid replacement rather than a regression.

**V and J genes** (single best alignment per template):
- Score difference histogram: for each (read, template) pair, compute `score_seqan2 - score_igor`. Expect the distribution to be centered near 0 with spread ≤ 5. Large discrepancies indicate a misconfigured `AlignConfig` or a scoring scheme mismatch.
- `offset` agreement rate: expect ≥95% of pairs to agree on `offset` (same alignment start).
- `five_p_offset` / `three_p_offset` agreement: same threshold.

**D gene** (multi-hit):
- SeqAn2 should return **at least as many** hits above threshold as IGoR, because Waterman-Eggert is exhaustive while IGoR's region tracker may merge hits.
- Compare the set of `(five_p_offset, three_p_offset)` intervals between the two aligners. Compute intersection-over-union (IoU) per sequence to quantify overlap.
- Investigate any case where IGoR returns a hit that SeqAn2 does not: this would imply a D-gene alignment that IGoR's region tracker finds but Waterman-Eggert misses (possible if a hit lies within the band of a prior high-scoring hit).

### Step 4.3 — Multithreading scaling
```bash
for T in 1 2 4 8 16; do
  OMP_NUM_THREADS=$T ./bench_seqan2_vs_igor \
    --aligner seqan2 --gene V --seqs demo/... --threads $T
done
```
Report wall time vs thread count. IGoR's `Aligner` already uses OpenMP; include it in the same scaling test for a fair comparison. Expected outcome: both aligners scale linearly to the number of cores up to memory-bandwidth saturation.

### Step 4.4 — SIMD hardware acceleration (`IGOR_SEQAN2_SIMD=ON` only)

**Background**: modern CPUs can process multiple integer lanes in parallel using SIMD (Single Instruction Multiple Data) registers. SSE4 provides 128-bit registers (8 × int16 per register); AVX2 provides 256-bit registers (16 × int16). SeqAn2's `<seqan/align_parallel.h>` exploits this by packing multiple sequence pairs into a single SIMD register and executing the DP fill for all of them simultaneously (inter-sequence parallelism).

SeqAn2's `<seqan/align_parallel.h>` execution modes:

| `ExecutionPolicy` | Parallelism | SIMD | Throughput (single core) |
|---|---|---|---|
| `<Parallel, Serial>` | OpenMP thread chunking | none | N× on N cores |
| `<Parallel, Vectorial>` | OpenMP + inter-sequence SIMD | SSE4 / AVX2 | ×8 per core (SSE4+int16), ×16 (AVX2+int16) |
| `<WavefrontAlignment<BlockOffsetOptimization>, Vectorial>` | intra-sequence wave-front + SIMD | SSE4 / AVX2 | best for sequences >1000 nt |

Compile flags: `-msse4` (SSE4) or `-mavx2` (AVX2). `SEQAN_DEFINITIONS` macro required.

**Critical constraint**: all batch `ExecutionPolicy` modes return **scores only** — no traceback, no positions. The `globalAlignmentScore` / `localAlignmentScore` functions with an `ExecutionPolicy` argument do not produce an `Align` object.

To produce `Alignment_data` (CIGAR, `five_p_offset`, `three_p_offset`), use a **two-pass strategy**:

**Pass 1 — Fast score filter** (SIMD-vectorized, all pairs):
```cpp
// Pack all (read, germline) pairs and run score-only DP
using TExecPolicy = seqan2::ExecutionPolicy<seqan2::Parallel, seqan2::Vectorial>;
TExecPolicy exec;
seqan2::setNumThreads(exec, n_threads);
// seqs_query[i] and seqs_germline[i] are aligned against each other
auto scores = seqan2::globalAlignmentScore(exec, seqs_query, seqs_germline, scoring_scheme);
```

**Pass 2 — Traceback for passing pairs only** (single-pair, standard `globalAlignment`):
```cpp
for (int i = 0; i < scores.size(); ++i) {
    if (scores[i] < threshold) continue;
    // Reconstruct full alignment for this pair only
    results[i] = align_single_pair(seqs_query[i], seqs_germline[i], scoring_scheme);
}
```

This two-pass strategy means the expensive traceback is only performed for the fraction of pairs that pass the score threshold, while the score filter benefits from SIMD throughput.

The two-pass mode is optional, gated on `IGOR_SEQAN2_SIMD`. For typical IGoR sequence lengths (200–400 nt reads, ~300 nt germlines) with a moderate number of sequences, the OpenMP loop in Phase 2's `align_seqs` often already saturates available cores; SIMD two-pass gives the largest gains for large datasets (tens of thousands of sequences).

Benchmark targets: sequential vs OpenMP vs SIMD two-pass wall time on 10k sequences. Expected: ≥4× speedup with `-mavx2` vs sequential (accounting for the fraction of pairs that require traceback).

API reference for completeness:
```cpp
// Wavefront + SIMD for long sequences (score-only pass):
using TExecPolicy2 = seqan2::ExecutionPolicy<
    seqan2::WavefrontAlignment<seqan2::BlockOffsetOptimization>,
    seqan2::Vectorial>;
TExecPolicy2 exec2;
seqan2::setNumThreads(exec2, n_threads);
seqan2::setBlockSize(exec2, 500);          // DP matrix block size (tune per hardware)
seqan2::setParallelAlignments(exec2, 10);  // max concurrent wavefront alignments
```

---

## Relevant Files
| File | Change |
|---|---|
| `src/igor/Core/Aligner.h` | add `CORE_EXPORT` declarations for Phase 1 free functions |
| `src/igor/Core/Aligner.cpp` | Phase 0 CSV fix + Phase 1 CIGAR functions |
| `src/igor/Core/AlignmentPreset.h` | **new** (Phase 3; no SeqAn dependency, unconditional) |
| `src/igor/Core/SeqAn2Aligner.h` | **new** (Phase 2–3; `#ifdef IGOR_WITH_SEQAN2`) |
| `src/igor/Core/SeqAn2Aligner.cpp` | **new** |
| `src/igor/Core/CMakeLists.txt` | conditional SeqAn2 sources + `IGOR_SEQAN2_SIMD` option |
| `CMakeLists.txt` (root) | `IGOR_WITH_SEQAN2` + `IGOR_SEQAN2_SIMD` options |
| `pixi.toml` | `[feature.seqan2]` with `bioconda::seqan` |
| `tst/test_alignment_cigar.cpp` | **new** (Phase 1) |
| `tst/test_seqan2_aligner.cpp` | **new** (Phase 2; `#ifdef IGOR_WITH_SEQAN2`) |
| `tst/test_alignment_preset.cpp` | **new** (Phase 3) |
| `tst/test_band_calibration.cpp` | **new** (Phase 3, Step 3.6; `#ifdef IGOR_WITH_SEQAN2`) |
| `tst/bench_seqan2_vs_igor.cpp` | **new** (Phase 4; separate CMake target, not part of `igor_tests`) |

---

## Verification

Each verification step corresponds to a concrete passing criterion, not just "run and hope":

1. **Phase 0 + Phase 1 (no SeqAn2)**:
   ```bash
   cmake -DIGOR_WITH_SEQAN2=OFF .. && ninja igor_tests && ./tst/igor_tests "[cigar]"
   ```
   All CIGAR tests pass. Existing tests must not regress.

2. **CSV round-trip** (Phase 0 prerequisite for Phase 1 CSV test):
   Load a V-gene CSV file → call `alignment_data_to_cigar()` → assert the result is non-empty and parses without error → call `alignment_data_from_cigar` → assert `offset`, `five_p_offset`, `three_p_offset` match the originally loaded values.

3. **CIGAR round-trip** (Phase 1 standalone):
   For any `Alignment_data` with known fields: `alignment_data_from_cigar(alignment_data_to_cigar(aln))` must reproduce `offset`, `five_p_offset`, `three_p_offset`, `insertions`, `deletions`, `mismatches`, `align_length` identically. Note: round-tripping through `M` ops legitimately loses mismatch positions.

4. **SeqAn2 build** (Phase 2):
   ```bash
   pixi run -e seqan2 cmake -DIGOR_WITH_SEQAN2=ON .. && ninja igor_tests
   ./tst/igor_tests "[seqan2]"
   ```
   All SeqAn2 tests pass.

5. **AlignmentPreset factory values** (Phase 3, Steps 3.1–3.4):
   ```bash
   ./tst/igor_tests "[preset]"
   ```
   Every factory method assertion passes.

6. **Band calibration** (Phase 3, Steps 3.5–3.6):
   ```bash
   ./tst/igor_tests "[calibration]"
   ```
   - Synthetic 5' RACE: calibrated band contains $[45, 115]$ and width $\leq 80$.
   - Synthetic Multiplex PCR: band centred within 5 nt of $-70$.
   - Stage 2 expansion triggered test passes.
   - Reproducibility: two runs with same seed return identical band bounds.
   - Full-dataset accuracy: ≥99% of best scores match unbanded baseline.

7. **Correctness on demo data** (Phase 4):
   - V-gene: SeqAn2 scores within ±5 of IGoR scores on `demo/murugan_naive1_noncoding_demo_seqs.txt`; offset agreement ≥90%.
   - D-gene: SeqAn2 returns ≥ IGoR hit count per sequence at the same score threshold; IoU ≥0.7 for matching intervals.
   - SIMD two-pass (if `IGOR_SEQAN2_SIMD=ON`): ≥4× wall-time speedup vs sequential baseline on 10k sequences with `-mavx2`.

---

## Decisions and Scope

### What is NOT in scope
- **No AIRR TSV I/O**: Phase 1 adds CIGAR conversion functions and AIRR coordinate accessors but does not add a full AIRR TSV writer/reader.
- **No changes to `Aligner` class**: the existing `Aligner` class and its CSV format are unchanged (except the Phase 0 reader bug fix).
- **No CSV format changes**: the on-disk CSV format is unchanged; Phase 0 only fixes the reader to parse what the writer already writes.
- **`SeqAn2Aligner` is an addition, not a replacement**: the existing `Aligner` class remains in the codebase and continues to be the default. `SeqAn2Aligner` is an opt-in alternative.

### Key technical decisions
- **SeqAn2 not SeqAn3**: SeqAn3 does not have Waterman-Eggert. SeqAn2 v2.5.3 from bioconda is the target.
- **`seqan2::` namespace**: SeqAn2 version 2.5.3 uses `seqan2::` (older versions used `seqan::`). All code must use `seqan2::` to avoid namespace collisions.
- **pimpl in `SeqAn2Aligner.h`**: SeqAn2 headers are large and introduce heavy template machinery. Using a pimpl prevents them from leaking into every translation unit that includes `SeqAn2Aligner.h`.
- **`AlignConfig` order** (confirmed from SeqAn2 2.5.3 source): `AlignConfig<seq2_leading, seq1_leading, seq1_trailing, seq2_trailing>`. seq1 = read (horizontal, `row(ali,0)`), seq2 = germline (vertical, `row(ali,1)`).
  - V gene: `<false, false, true, true>`
  - J no primer: `<true, true, true, false>`
  - J 3' primer: `<true, true, false, false>`
- **Phase 2 uses `SeqAn2AlignConfig`**; Phase 3 replaces it 1:1 with `AlignmentPreset` — the field names and semantics are identical, so refactoring is mechanical.
- **`AlignmentPreset.h` is unconditional**: it carries no SeqAn2 dependency and should be usable with any future aligner backend (e.g., a potential SeqAn3 backend once Waterman-Eggert is available).
- **SIMD requires compile flags** (`-msse4`/`-mavx2`) and `SEQAN_DEFINITIONS`; gated on `IGOR_SEQAN2_SIMD` CMake option to avoid forcing a specific CPU ISA on all builds.
- **SIMD batch modes are score-only**: the two-pass strategy is mandatory to produce `Alignment_data` when SIMD is enabled.
- **`FindSeqAn.cmake` location**: the bioconda `seqan` package installs `FindSeqAn.cmake` under `${CONDA_PREFIX}/share/seqan/cmake/`. If `cmake` does not find it automatically via `CMAKE_PREFIX_PATH`, set `CMAKE_MODULE_PATH` explicitly.
- **Phase ordering**: 0 → 1 → 2 → 3 → 4 sequential. Phase 1 can be done partially in parallel with Phase 0 (the CSV round-trip test requires Phase 0, but the other CIGAR tests do not).
- **`calibrate_preset` as free function, not method**: keeps `AlignmentPreset.h` fully SeqAn2-free; `BandCalibrationParams` lives in `AlignmentPreset.h` (unconditional) since it contains only plain integers.
- **Calibration runs at `SeqAn2Aligner` construction, not lazily**: ensures the band is deterministic and fixed for the lifetime of the aligner object. Thread-safe by construction.
- **Stage 1 uses full-matrix (unbanded) alignment on a sample, not the legacy `Aligner`**: avoids any dependency on `Aligner` being available and keeps the calibration self-contained within the SeqAn2 backend. If a prior `Aligner` pass has already been run, its `offset` values can be fed directly as pre-computed diagonals (API accepts `std::vector<int> known_diagonals` as an optional shortcut; documented in `calibrate_preset` overload).
