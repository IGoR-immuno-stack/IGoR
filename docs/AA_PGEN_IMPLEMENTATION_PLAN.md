# AA-Level Pgen Implementation Plan

**Created**: May 2026  
**Status**: Partially implemented on `feature/AA_PGEN` — Phases 1–3 complete, Phases 4–7
landed but incomplete/incorrect in the places listed in **Part II**, Phase 8 largely missing.  
**Last status review**: 2026-08-27 (commits `35440f9`..`88a8dfa`)  
**Relates to**: `ITERATE_METHOD_ANALYSIS.md`, context refactoring (Phases 1–5)

> **Read Part II first.** Sections 1–9 below are the *original design intent* and have been
> left unedited so that deviations can be diffed against them. Part II (Section 10 onwards)
> records what was actually built, where it departs from this design, what is left, and what
> must be watched when merging with the `iterate()` refactor branch.

---

## 1. Problem Statement

IGoR currently computes generation probabilities (`Pgen`) of immune receptor sequences at the
**nucleotide level**. Computing Pgen at the **amino acid level** — i.e., the probability of
generating any nucleotide sequence that encodes a given amino acid sequence — is a desirable
extension for several immunological analyses.

The naïve approach of enumerating all synonymous nucleotide sequences and summing their Pgen
values is **combinatorially intractable**: a receptor of length 100 aa with all Leu positions
would have 6^(number of Leu) synonymous sequences.

This document describes an efficient approach that integrates directly into the existing
`Rec_Event::iterate()` workflow.

---

## 2. The Encoding Obstacle: Leu, Arg, Ser

For 17 of the 20 amino acids, a single IUPAC codon triplet exactly represents the synonymous
codon set (e.g., Phe = TTY, Gly = GGN). For three amino acids, positional correlations between
codon positions make a single IUPAC triplet insufficient:

| Amino acid | Codons | Naïve IUPAC | Problem |
|-----------|--------|-------------|---------|
| **Leu** | TTA, TTG, CTT, CTC, CTA, CTG | YTN | Includes TTT=Phe, TTC=Phe |
| **Arg** | CGT, CGC, CGA, CGG, AGA, AGG | RGR | Includes wrong triplets |
| **Ser** | TCT, TCC, TCA, TCG, AGT, AGC | WSN | Includes many non-Ser |

The root cause is that the allowed set at position 2 of a Leu codon depends on position 0:
if position 0 = C then position 2 ∈ {A,C,G,T}, but if position 0 = T then position 2 ∈ {A,G}.
A single IUPAC code cannot encode this correlation.

These three amino acids can be represented as the union of two exact IUPAC triplets:

| Amino acid | Group 1 | Group 2 |
|-----------|---------|---------|
| Leu | `TTR` | `CTN` |
| Arg | `CGN` | `AGR` |
| Ser | `TCN` | `AGY` |

This motivates a per-codon allowed-set structure for exact leaf-node accounting, combined with
a single IUPAC reference sequence for alignment and conservative mismatch detection.

### 2.3 AA Motif Queries: Subsets and Wildcards

The same obstacle — and its solution — generalizes to motif queries where one or more positions
allow a *set* of amino acids rather than a single one. Two notations are supported:

- **Wildcard `X`**: any amino acid accepted at that position (all non-stop codons allowed).
- **Bracket group `[KSM]`**: exactly the listed amino acids are accepted. This is the notation
  used by OLGA (e.g., `CAV[KSM]DSNYQLI[WF]`). IUPAC parenthesis notation `(K,S,M)` conveys
  the same intent but is not standardised in practice; bracket notation is the canonical input
  format here.

For a position accepting amino acids $\{a_1, a_2, \ldots, a_k\}$, the allowed codon set is the
union of the synonymous codon sets for each listed amino acid. This maps directly onto the
`CodonMask` type (defined in Phase 1.3):

```cpp
// Single AA (exact query):
CodonMask mask_M   = codon_mask_for_aa('M');              // 1 bit (ATG)

// AA subset [KSM]:
CodonMask mask_KSM = codon_mask_for_aa('K')
                   | codon_mask_for_aa('S')
                   | codon_mask_for_aa('M');              // up to 18 bits

// Wildcard X:
CodonMask mask_X = 0;
for (char aa : "ACDEFGHIKLMNPQRSTVWY") mask_X |= codon_mask_for_aa(aa);  // 61 bits
```

The IUPAC reference nucleotide at each codon position is built as the union of all nucleotides
that appear at that position in any codon allowed by the mask. For large subsets this frequently
collapses to `int_N` (matches anything), which is correct: those positions impose no mismatch
penalty during alignment.

`CodonMask` is a **preprocessing-only** utility: it is the internal data structure used by the
factory function (Phase 2) to compute per-position IUPAC sets and patch alternatives. It does
not appear in the iterate loop, mismatch recording, or leaf-check code paths at runtime.

---

## 3. Existing Infrastructure

Before describing new code, the relevant existing infrastructure:

### 3.1 `Int_nt` encoding and IUPAC support

`Utils.h` defines all 15 IUPAC nucleotide codes as an enum:

```cpp
enum Int_nt {
    int_A = 0, int_C = 1, int_G = 2, int_T = 3,
    int_R = 4, int_Y = 5, int_K = 6, int_M = 7,
    int_S = 8, int_W = 9, int_B = 10, int_D = 11,
    int_H = 12, int_V = 13, int_N = 14
};
```

`Aligner.cpp` already provides:
- `nt2int()`: parses all 15 IUPAC characters into `Int_Str`
- `comp_nt_int(nt_1, nt_2)`: returns `true` if the two codes are compatible (at least one
  shared nucleotide). Used by `sw_align()` for mismatch detection during alignment.
- `get_ambiguous_nt_list()`: expands any IUPAC code to its constituent A/C/G/T list.

The Smith-Waterman aligner already uses `comp_nt_int` when computing mismatch lists in
`Alignment_data`. **This means the aligner can already align IUPAC-encoded query sequences.**

### 3.2 The mismatch list's current contract

`Mismatch_vectors_map` (typedef in `Utils.h`) maps `Seq_type` → `vector<int>*` where the
integers are positions in the query where `comp_nt_int(gene_nt, query_nt) == false`.

This list is consumed by two independent clients:
1. **`downstream_proba_map`** (in `ExplorationContext`): uses the mismatch count to compute an
   upper bound on the error-rate contribution for pruning.
2. **`Error_rate::compute_scenario_error_probability()`** (at leaf nodes): interprets the
   positions to compute the actual error-weighted probability.

Currently both clients see the same list, implicitly assuming that `comp_nt_int` semantics are
appropriate for both. This is an approximation: for IUPAC queries, `comp_nt_int` returns
**no mismatch** whenever at least one branch of the query is compatible — which is correct for
pruning (conservative) but incorrect for leaf accumulation (under-counts mismatches for
incompatible branches).

### 3.3 `get_err_rate_upper_bound(n_mismatch, n_genomic)`

`Error_rate` provides a precomputed matrix accessed via `get_err_rate_upper_bound(n, m)` that
returns the maximum error-weighted probability given `n` mismatches and `m` non-mismatch
positions. This is called in `Gene_choice` and `Deletion` to set the error component of the
downstream probability bound.

---

## 4. Core Design: Journaled Query Representation

### 4.1 Two-track mismatch semantics

During each iterate step in `Gene_choice` and `Deletion`, a position-by-position comparison
loop already runs to build pruning bounds. That **same loop** records two results
simultaneously, with no additional pass over the gene positions:

- **Floor** (`pruning_mismatch_floor`): `!comp_nt_int(gene_nt, query.iupac_union[p])` — position
  is incompatible with **every** alternative NT at that position — confirmed mismatch in all
  branches. Used for conservative downstream pruning (`downstream_proba_map`).
- **Upper bound** (`mismatches_lists`):
  `query.empty_isect[p] || !comp_nt_int(gene_nt, query.iupac_intersection[p])` — position is
  a mismatch in **at least one** branch. When `empty_isect[p]` is true (no IUPAC code
  represents the intersection, e.g. position 0 of Leu), the position is unconditionally
  flagged as an upper-bound mismatch.

The loop is **mode-agnostic**: there is no branch on `journaled_query` presence.

**Backward compatibility**: for exact NT queries, `iupac_union == iupac_intersection ==
reference` (all codes in [0..3]), `empty_isect` is all-false, and `comp_nt_int` reduces to
`gene_nt == query_nt`. Both tracks are identical. No behaviour change for existing workflows.

### 4.2 `JournaledQuery` structure

The `JournaledQuery` struct unifies exact NT queries, IUPAC NT queries, AA motifs, and
arbitrary NT-patch queries under one representation. A `Patch` describes a span of positions
where the query has multiple valid NT alternatives (not just the reference).

```cpp
/// One span of positions with multiple valid NT alternatives.
/// Covers arbitrary-length patches (e.g. one codon = 3 positions).
struct Patch {
    int start;                         // first position within the full receptor nt seq
    int length;                        // number of positions in this patch
    std::vector<Int_Str> alternatives; // valid NT subsequences (each of length `length`)
                                       // does NOT include the reference subsequence
};

/// Query representation for Pgen computation — exact NT, IUPAC NT, AA motif, or patched NT.
struct JournaledQuery {
    Int_Str reference;            // one concrete valid NT sequence (leaf check + display)
    Int_Str iupac_union;          // per-pos IUPAC union of all alternatives (alignment + floor)
    Int_Str iupac_intersection;   // per-pos IUPAC intersection of all alternatives (upper bound)
    std::vector<bool> empty_isect;// true where intersection is ∅ → always flag upper-bound
    std::vector<Patch> patches;   // empty for exact NT queries
    std::string display_string;   // original input (AA motif, NT, etc.) for diagnostics
};
```

**Field semantics**:
- `reference`: a single concrete valid NT sequence. Used for alignment and for leaf-node
  display. For AA motifs, this is an arbitrary representative encoding (e.g., first valid
  codon per position). For exact NT queries, it equals the query itself.
- `iupac_union`: per-position IUPAC union across all valid NT sequences (i.e., all NT values
  that appear at position `p` in any valid sequence). Stored in
  `QuerySequenceContext::int_sequence`. Used for alignment and for the floor track.
- `iupac_intersection`: per-position IUPAC intersection across all valid NT sequences (only
  NT values that appear at position `p` in **every** valid sequence). Used for the upper-bound
  track. For Leu/Arg/Ser codon position 0, the intersection can be empty.
- `empty_isect[p]`: set to `true` when the intersection at position `p` is empty (no IUPAC
  code can represent it). Position is unconditionally flagged as an upper-bound mismatch.
- `patches`: list of spans where the query has alternatives. Empty for exact NT queries (all
  fields are then trivially identical: `union == intersection == reference`, `empty_isect` all
  false, `patches` empty).
- `display_string`: the original human-readable input, for diagnostics and logging.

---

## 5. Context Placement

Consistent with the existing context architecture:

| Member | Context | Rationale |
|--------|---------|-----------|
| `mismatches_lists` | `ScenarioContext` | Upper-bound mismatch positions; consumed by error rate at leaf |
| `pruning_mismatch_floor` | `ExplorationContext` | NT-floor mismatch positions; consumed by `downstream_proba_map` |
| `journaled_query` | `QuerySequenceContext` | Read-only input; absent for exact NT inference (trivially equivalent when present with no patches) |

`int_sequence` holds `journaled_query.iupac_union` in patched/motif mode. For exact NT
queries `iupac_union == reference`, so `int_sequence` is unchanged from today.

Both `mismatches_lists` and `pruning_mismatch_floor` require the same memory-layer treatment:
each recursion level uses a distinct layer to allow backtracking without copying.

---

## 6. Implementation Phases

### Phase 1 — Utility Foundation

**Files**: `Aligner.h`, `Aligner.cpp`, `EventUtils.h`, `EventUtils.cpp`

#### 1.1 Fix D-gene exhaustive search path

`Genechoice.cpp` lines 625 and 733 use direct `!=` comparison in the D-gene exhaustive search
path, inconsistent with the aligner's `comp_nt_int` semantics:

```cpp
// BEFORE (lines 625, 733):
if (gene_seq[i] != query.int_sequence[d_5_off + i]) {
    no_d_mismatches.push_back(d_5_off + i);

// AFTER:
if (!comp_nt_int(gene_seq[i], query.int_sequence[d_5_off + i])) {
    no_d_mismatches.push_back(d_5_off + i);  // floor semantics, consistent with aligner
```

For exact NT queries this is a no-op. For IUPAC queries it is a correctness fix.

#### 1.2 Add genetic code utilities to `EventUtils`

```cpp
// EventUtils.h
namespace EventUtils {

// Codon index encoding: codon_index(n0,n1,n2) = n0*16 + n1*4 + n2
// where n0,n1,n2 ∈ {int_A=0, int_C=1, int_G=2, int_T=3}
constexpr int codon_index(int n0, int n1, int n2) { return n0*16 + n1*4 + n2; }

// 64-bit bitmask: bit k set iff codon k (in codon_index encoding) encodes the target AA
using CodonMask = uint64_t;

// Standard genetic code translation
CORE_EXPORT char translate_codon(int n0, int n1, int n2);

// Returns CodonMask for a given amino acid (single-letter code; '*' = stop)
CORE_EXPORT CodonMask codon_mask_for_aa(char aa);

// Translate a complete Int_Str; reading frame starts at frame_offset.
// Returns empty string if length past frame_offset is not a multiple of 3.
CORE_EXPORT std::string translate_int_seq(const Int_Str& seq, int frame_offset);

// For a single amino acid, return the IUPAC triplet that conservatively covers all
// synonymous codons (union of nucleotide sets per position).
// For Leu/Arg/Ser this over-approximates; exact handling uses CodonMask.
CORE_EXPORT std::array<int,3> aa_to_iupac_codon(char aa);

// For a CodonMask, return the IUPAC nucleotide triplet [pos0, pos1, pos2] that conservatively
// covers all codons whose bit is set (union of nucleotide sets per codon position).
// Generalises aa_to_iupac_codon to arbitrary codon subsets including motif groups.
CORE_EXPORT std::array<int,3> mask_to_iupac_codon(CodonMask mask);

// For a motif character, return the corresponding CodonMask:
// - single-letter AA code: identical to codon_mask_for_aa(ch)
// - 'X': OR of all non-stop amino acid masks (any of the 20 standard amino acids)
// Throws std::invalid_argument for unrecognised characters.
CORE_EXPORT CodonMask motif_char_to_mask(char ch);

// Parse an OLGA-style bracket-notation motif into a vector of CodonMasks, one per position.
// Supports: single-letter AA codes, 'X' wildcard, '[...]' subset groups.
// Example: "CAV[KSM]DSX" → {mask_C, mask_A, mask_V, mask_K|mask_S|mask_M, mask_D, mask_S, mask_X}
// Throws std::invalid_argument on malformed input (unclosed bracket, empty bracket, etc.).
CORE_EXPORT std::vector<CodonMask> parse_aa_motif(const std::string& motif);

} // namespace EventUtils
```

`translate_codon` and `codon_mask_for_aa` are implemented with a 64-entry `constexpr` lookup
table (standard genetic code, no ambiguity needed). `mask_to_iupac_codon` iterates over bits
set in the mask to collect allowed nucleotides per codon position.

**Dependency**: none. Can be implemented and merged independently.

> **Note on `CodonMask` in Phase 1.2**: the `codon_index(n0,n1,n2)` encoding is defined here
> because it is needed in Phase 5 (`Single_error_rate` motif-mode computation) and Phase 6
> (`Dinucl_markov` AA sum). The genetic code lookup table does not appear anywhere in the
> mismatch-building or pruning code paths.

---

### Phase 2 — `JournaledQuery` and Factory

**New file**: `JournaledQuery.h`  
**Modified file**: `EventUtils.h`, `EventUtils.cpp`

The `JournaledQuery` and `Patch` structs (see Section 4.2) are declared in a new header
included by `QuerySequenceContext.h`. The factory `motif_to_journaled_query()` is the primary
construction path; `CodonMask` is used internally only.

Factory function in `EventUtils`:

```cpp
// EventUtils.h
namespace EventUtils {

/// Build a JournaledQuery from an OLGA-style AA motif string.
/// @param motif        AA motif, e.g. "CAVX[KSM]DS" (single letters, 'X', or '[...]' groups)
/// @param frame_offset Position of the first codon start within the receptor nt sequence
/// @param receptor_len Total expected nt length of the receptor
CORE_EXPORT JournaledQuery motif_to_journaled_query(
    const std::string& motif,
    int frame_offset,
    int receptor_len
);

/// Convenience overload for a plain AA sequence (no brackets); equivalent to motif input.
inline JournaledQuery aa_to_journaled_query(
    const std::string& aa_seq,
    int frame_offset,
    int receptor_len)
{
    return motif_to_journaled_query(aa_seq, frame_offset, receptor_len);
}

} // namespace EventUtils
```

Implementation of `motif_to_journaled_query`:
1. Call `parse_aa_motif(motif)` to get `std::vector<CodonMask> masks` (one per codon position).
2. Allocate `reference`, `iupac_union`, `iupac_intersection` of `receptor_len` positions,
   initialized to `int_N` for non-codon positions. Allocate `empty_isect` of `receptor_len`
   initialized to `false`.
3. For each mask at position `i`:
   a. Enumerate all valid concrete codons (NT triplets whose `codon_index` bit is set in
      `masks[i]`). Build the list of alternatives as `std::vector<Int_Str>` (each of length 3).
   b. Pick the first valid concrete codon as `reference` for positions
      `[frame_offset + 3i, frame_offset + 3i + 3)`.
   c. For each of the three codon positions `j ∈ {0,1,2}`, collect the set of NT values that
      appear at position `j` across all valid codons. Encode the union as an IUPAC code →
      store in `iupac_union[frame_offset + 3i + j]`. Encode the intersection as an IUPAC
      code → store in `iupac_intersection[frame_offset + 3i + j]`. If the intersection is
      empty (no NT value common to all valid codons), set `empty_isect[frame_offset + 3i + j]
      = true` and store any placeholder (e.g. `int_N`) in `iupac_intersection` at that
      position (it will not be consulted at runtime).
   d. Construct a `Patch{frame_offset + 3i, 3, alternatives}` where `alternatives` contains
      all valid concrete codon `Int_Str` values **except** the one chosen as `reference`.
      Append to `patches`. (Omit the patch if the mask has exactly one valid codon —
      `alternatives` would be empty, making the patch trivial.)
4. Store `motif` as `display_string` in the result.

There is no special-casing for Leu/Arg/Ser: their multi-group structure is encoded exactly in
their `CodonMask`, and the per-position set enumeration in step 3c handles the empty-
intersection case uniformly.

**Dependency**: Phase 1 (`parse_aa_motif`, `codon_mask_for_aa`, `codon_index`).

---

### Phase 3 — `QuerySequenceContext` Extension

**Modified file**: `QuerySequenceContext.h`

```cpp
#include <igor/Core/JournaledQuery.h>
#include <optional>

struct QuerySequenceContext {
    const std::string& sequence;
    const Int_Str& int_sequence;   // iupac_union in patched/motif mode; exact NT otherwise
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>& gene_alignments;

    // Present only in patched/motif mode; std::nullopt for standard NT inference
    const std::optional<JournaledQuery> journaled_query;

    // NT mode constructor (unchanged behaviour)
    QuerySequenceContext(
        const std::string& sequence_,
        const Int_Str& int_sequence_,
        const std::unordered_map<Gene_class, std::vector<Alignment_data>>& gene_alignments_
    ) : sequence(sequence_), int_sequence(int_sequence_),
        gene_alignments(gene_alignments_), journaled_query(std::nullopt) {}

    // Patched/motif mode constructor
    QuerySequenceContext(
        const std::string& sequence_,
        const Int_Str& int_sequence_,
        const std::unordered_map<Gene_class, std::vector<Alignment_data>>& gene_alignments_,
        JournaledQuery jq
    ) : sequence(sequence_), int_sequence(int_sequence_),
        gene_alignments(gene_alignments_), journaled_query(std::move(jq)) {}

    QuerySequenceContext(const QuerySequenceContext&) = delete;
    QuerySequenceContext& operator=(const QuerySequenceContext&) = delete;
    QuerySequenceContext(QuerySequenceContext&&) = delete;
    QuerySequenceContext& operator=(QuerySequenceContext&&) = delete;
};
```

`int_sequence` holds `journaled_query->iupac_union` in patched/motif mode. All existing
iterate code (gene comparisons, `comp_nt_int` calls, offset-based substr operations) works
transparently because `iupac_union` is a valid `Int_Str` of the same length as the receptor.

**Dependency**: Phase 2.

---

### Phase 4 — Two-track Mismatch Recording

**Modified files**: `Utils.h`, `ExplorationContext.h`, `Genechoice.cpp`, `Deletion.cpp`

#### 4.1 New typedef

```cpp
// Utils.h (alongside existing typedefs)
typedef Enum_fast_memory_map<Seq_type, std::vector<int>*> Pruning_mismatch_floor_map;
```

Same underlying type as `Mismatch_vectors_map`.

#### 4.2 Add `pruning_mismatch_floor` to `ExplorationContext`

```cpp
struct ExplorationContext {
    // ... existing members ...

    /// NT-floor mismatch positions per Seq_type.
    /// Position p is in floor[seg] iff !comp_nt_int(gene_nt, query.int_sequence[p]).
    /// Invariant: floor[seg] ⊆ mismatches_lists[seg] (floor is subset of upper bound).
    /// Consumed exclusively by downstream_proba_map (conservative pruning bound).
    /// Memory-layer management: identical to mismatches_lists in ScenarioContext.
    Pruning_mismatch_floor_map& pruning_mismatch_floor;
};
```

#### 4.3 Single-pass two-track computation in `Gene_choice` and `Deletion`

The existing comparison loop (one pass per position) is extended to write both tracks
simultaneously from the same evaluation. The loop is **fully mode-agnostic**: no branch on
`journaled_query` presence.

```cpp
// Inside Gene_choice/Deletion position loop (single pass, no extra iteration):
bool floor_mismatch = !comp_nt_int(gene_nt, query.int_sequence[qpos]);
                   // int_sequence == iupac_union (patched mode) or exact NT (NT mode)

bool upper_mismatch;
if (query.journaled_query.has_value()) {
    // Patched/motif mode: use iupac_intersection + empty_isect
    upper_mismatch = query.journaled_query->empty_isect[qpos]
                  || !comp_nt_int(gene_nt, query.journaled_query->iupac_intersection[qpos]);
} else {
    // Exact NT mode: is_exact_match upper bound
    upper_mismatch = (gene_nt != query.int_sequence[qpos]);
}

if (floor_mismatch) floor_vec.push_back(qpos);
if (upper_mismatch) upper_vec.push_back(qpos);
```

**Why the NT-mode branch is retained**: `iupac_intersection` and `empty_isect` are fields of
`JournaledQuery`, which is absent in NT mode — they cannot be accessed without the branch.
For **exact** NT queries (all codes in [0..3]), `comp_nt_int(a, b)` returns `true` iff
`a == b`, so `!comp_nt_int` reduces to `!=` with no behavioural difference. The branch is
structural necessity for that case.

> **Warning — IUPAC NT queries without `journaled_query`**: the existing code already accepts
> IUPAC-encoded NT queries in the non-`journaled_query` path. In that case the `!=` fallback
> is **semantically incorrect** for the upper-bound track: `gene_nt = int_A` against
> `query_nt = int_R` (purine) satisfies `comp_nt_int` (compatible) but fails `!=` (different
> integer values), so the position is over-counted as a mismatch. Any caller that passes an
> IUPAC NT query without a `JournaledQuery` will get an inflated `mismatches_lists` and a
> spuriously low accumulated probability. The correct fix is to always wrap IUPAC NT queries
> in a `JournaledQuery` (with `iupac_union == iupac_intersection == the query itself`, empty
> `patches`, `empty_isect` all false); this routes them through the `comp_nt_int`-based
> upper-bound path and gives exact results.

`downstream_proba_map` entries use `floor_vec.size()` (conservative lower bound).
`mismatches_lists` receives `upper_vec`, which `Single_error_rate` reads at the leaf
(Phase 5).

`Deletion::iterate()` trims both vectors identically as deletions grow (same existing trim
logic applied to each vector in lock-step).

Memory-layer management: `pruning_mismatch_floor` and `mismatches_lists` request and
release layers in lock-step.

**D-gene exhaustive path**: the two-track loop replaces the direct `!=` check on lines 625
and 733 of `Genechoice.cpp`, resolving the Phase 1.1 D-gene fix in the same change.

**Dependency**: Phase 3 not required; Phase 4 applies correctly for both NT and motif queries.

---

### Phase 5 — Patch-Aware Leaf Check in `Single_error_rate`

**Modified files**: `Singleerrorrate.h`, `Singleerrorrate.cpp`

`Single_error_rate::compute_scenario_error_probability` at the leaf reuses the pre-computed
`mismatches_lists` (upper-bound track from Phase 4). The `comp_nt_int` comparisons already
happened during traversal and are **not re-run at the leaf**.

**NT mode** (no `journaled_query`): `mismatches_lists[seg].size()` is the full NT-level
mismatch count. Apply error model directly. No sequence reconstruction, no change.

**Patch mode** (`journaled_query` present with non-empty `patches`): the leaf check verifies
whether the assembled concrete sequence lies within the set of valid sequences encoded by
`journaled_query`. The unit of mismatch is one **patch** (an `alternative` span), not one
nucleotide.

```cpp
double Single_error_rate::compute_scenario_error_probability(
    const QuerySequenceContext& query,
    const ModelContext& model,
    ScenarioContext& scenario,
    ExplorationContext& exploration)
{
    size_t mismatch_count = 0;
    size_t total_positions = 0;

    if (query.journaled_query.has_value() && !query.journaled_query->patches.empty()) {
        // ── Patch/motif mode ───────────────────────────────────────────────────
        const JournaledQuery& jq = *query.journaled_query;
        total_positions = jq.patches.size();

        // Collect upper-bound mismatch positions (pre-computed during traversal).
        std::unordered_set<int> upper_mismatches;
        for (const auto& [seg, vec_ptr] : scenario.mismatches_lists)
            for (int p : *vec_ptr)
                upper_mismatches.insert(p);

        // Assemble full sequence for patch-level verification.
        // Dinucl_markov has already stored a compatible representative insertion.
        Int_Str full_seq = EventUtils::build_scenario_sequence(
            scenario.constructed_sequences);

        for (const auto& patch : jq.patches) {
            // Fast path: any position in the patch is a confirmed upper-bound mismatch.
            bool fast_miss = false;
            for (int k = patch.start; k < patch.start + patch.length && !fast_miss; ++k)
                if (upper_mismatches.count(k)) fast_miss = true;
            if (fast_miss) { ++mismatch_count; continue; }

            // Slow path: extract the assembled concrete subsequence for this patch span.
            Int_Str assembled = full_seq.substr(patch.start, patch.length);
            // Check against reference span and each alternative.
            Int_Str ref_span = jq.reference.substr(patch.start, patch.length);
            if (assembled == ref_span) continue;  // matches reference → no mismatch
            bool found = false;
            for (const auto& alt : patch.alternatives)
                if (assembled == alt) { found = true; break; }
            if (!found) ++mismatch_count;
        }
    } else {
        // ── NT mode (unchanged) ────────────────────────────────────────────────
        for (const auto& [seg, vec_ptr] : scenario.mismatches_lists)
            mismatch_count += vec_ptr->size();
        total_positions = compute_total_genomic_positions(scenario, model);
    }

    // Apply the error-rate formula — identical for both modes.
    //   error_rate == 0 → mismatch_count > 0 gives 0.0 (exact Pgen)
    //   error_rate  > 0 → patch mismatches penalised (fuzzy / hypermutated Pgen)
    return apply_error_model(mismatch_count, total_positions,
                             scenario.scenario_proba, exploration);
}
```

**Properties**:
- `comp_nt_int` is not re-run at the leaf — it is reused from `mismatches_lists`.
- Patches where any position has an upper-bound mismatch take the fast path.
- The slow path compares the assembled concrete subsequence against `reference` span and
  `patch.alternatives`. Full sequence assembly is still required when patches span segment
  boundaries (e.g., the last V-gene codon continues into the VD junction); insertions are
  handled by the representative stored by `Dinucl_markov` (Phase 6).
- The mechanism is general: `patches` of arbitrary length cover codon-level AA queries,
  dinucleotide patches, or any other NT-level local constraint.
- **No dedicated `AA_pgen_error_rate` subclass** needed.
- `error_rate == 0` → exact motif Pgen; `error_rate > 0` → fuzzy/hypermutated Pgen.
- `build_upper_bound_matrix` should be called with `patches.size()` as the maximum.

**Dependency**: Phase 3 (`QuerySequenceContext` extension) + Phase 4 (`mismatches_lists`
upper-bound semantics).

---

### Phase 6 — `Dinucl_markov` AA Pgen Mode

**Modified files**: `Dinuclmarkov.h`, `Dinuclmarkov.cpp`

For **gene-derived** positions (V, D, J sequence segments), the constructed nucleotides are
fixed at the time `Gene_choice` and `Deletion` run. The translation check at the leaf node
handles these positions correctly.

For **inserted** positions (VD, DJ, VJ junctions), the `Dinucl_markov::iterate()` currently
uses a greedy heuristic: pick the single most probable nucleotide sequence for the insertion.
For AA Pgen, this must be replaced by a **Markov forward sum** over all nucleotide sequences
that encode the correct amino acid at each codon position touching the insertion.

#### 6.1 Markov forward sum

For an insertion of length `L` governed by a first-order Markov chain with transition matrix
`T` (the `dinuc_proba_matrix`), and starting from previous nucleotide context `prev`:

$$P(\text{insertion encodes target codon sequence}) = \sum_{\substack{(n_0, n_1, \ldots, n_{L-1})\\ \text{encodes target AAs}}} T_{\text{prev}, n_0} \cdot \prod_{i=1}^{L-1} T_{n_{i-1}, n_i}$$

For a single codon wholly within the insertion (positions `k`, `k+1`, `k+2`):

$$\sum_{(n_k, n_{k+1}, n_{k+2}) \in \text{allowed\_codons}} T_{\text{prev}^{(k)}, n_k} \cdot T_{n_k, n_{k+1}} \cdot T_{n_{k+1}, n_{k+2}}$$

This is a sum over at most 6 terms (Leu, the worst case). For a full junction of length `L`
spanning multiple codons, the sum factorizes across codons if the junction length is a
multiple of 3. If not, terms span codon boundaries and the full forward sum must be computed
position by position using dynamic programming over the 4-state Markov chain restricted to
compatible nucleotides.

For codons straddling the gene/insertion boundary, one or two positions are fixed by the gene
template. Only the free positions within the insertion are summed:

```
V gene:  ... v_{n-1}  v_n   | VD insertion: i_0  i_1  i_2 | D gene: d_0 ...
                              ^-- codon boundary may be here
```

If the codon boundary falls at position `v_n | i_0`, only `i_0`, `i_1`, `i_2` are free;
`v_{n-1}` and `v_n` are fixed from the V gene template.

#### 6.2 New private helper

```cpp
// Dinuclmarkov.h (private section)

/// Compute the probability sum over all nucleotide sequences of length `ins_len`
/// that are compatible with the target codon structure, given:
/// @param prev_nt   Context nucleotide before the insertion (Int_nt value, 0..3)
/// @param ins_len   Length of the insertion
/// @param frame_start_in_ins  Which codon boundary the insertion starts on
///                  (0 = first base of a codon, 1 = second, 2 = third)
/// @param codon_masks  Vector of CodonMasks, one per full codon covered
/// @param model_params Marginal array for reading dinucleotide probabilities
/// @return Total probability weight for this insertion in the AA Pgen context
double compute_markov_aa_sum(
    int prev_nt,
    int ins_len,
    int frame_start_in_ins,
    const std::vector<EventUtils::CodonMask>& codon_masks,
    const Marginal_array_p& model_params
) const;
```

For each (prev_nt, target_aa, junction_length) combination, the result can be cached in a
precomputed table during `initialize()` to avoid repeated summation during inference.

#### 6.3 Mode switch in `iterate()`

```cpp
void Dinucl_markov::iterate(
    QuerySequenceContext& query, const ModelContext& model,
    ScenarioContext& scenario, ExplorationContext& exploration,
    AccumulationContext& accumulation)
{
    if (query.journaled_query.has_value() && !query.journaled_query->patches.empty()) {
        iterate_patched_pgen(query, model, scenario, exploration, accumulation);
        return;
    }
    // ... existing NT mode logic unchanged ...
}
```

`iterate_patched_pgen()` is a private method that:
1. Determines the reading frame position for the inserted junction from
   `journaled_query->patches` (which patches touch the insertion span) and current V/D/J
   offsets in `scenario.seq_offsets`.
2. Collects patch alternatives for each patch span overlapping the insertion region.
3. Calls `compute_markov_aa_sum()` (passing patch alternatives rather than `CodonMask`
   directly) to compute the total probability weight.
4. Sets `scenario.scenario_proba *= weight` (replaces the greedy heuristic contribution).
5. Calls `Rec_Event::iterate_wrap_up()`.

**Why this works with Phase 5**: `iterate_aa_pgen()` folds the Markov sum weight over all
compatible insertions into `scenario.scenario_proba` and stores a representative compatible
insertion in `constructed_sequences`. `Single_error_rate` at the leaf node assembles the
full sequence, finds zero codon mismatches (compatible insertion was stored), and returns
`scenario.scenario_proba` — the accumulated sum is therefore exactly
$P_\text{gen}^\text{AA}(\text{motif})$. For `error_rate > 0`, gene-derived codon mismatches
are penalised naturally by the same code path.

**Dependency**: Phase 3.

---

### Phase 7 — Input Preprocessing and Public API

**Modified files**: `GenModel.h`, `GenModel.cpp`

#### 7.1 Alignment with IUPAC reference

The aligner already handles IUPAC query sequences via `comp_nt_int` in `sw_align()`. No
aligner changes are needed. The caller passes `syn_query.reference` (the `Int_Str`) as the
query sequence.

Alignment offset bounds may need widening for IUPAC queries since the IUPAC reference has
`int_N` at non-codon positions, which matches every gene nucleotide. Standard bounds derived
from CDR3 annotation should still constrain V and J alignment correctly.

#### 7.2 Result type

```cpp
struct AAPgenResult {
    double pgen;        // P_gen(aa_seq) = sum over all contributing scenarios
    size_t n_scenarios; // number of scenarios encoding the target AA sequence
};
```

#### 7.3 Entry point

```cpp
// GenModel.h

/// Compute AA-level Pgen for a plain AA sequence or OLGA-style motif.
/// @param motif        AA sequence or bracket-notation motif (e.g. "CAVX[KSM]DS")
/// @param frame_offset Position of the first codon start within the receptor nt sequence
/// @param receptor_nt_len Expected full receptor nt length
CORE_EXPORT AAPgenResult compute_aa_pgen(
    const std::string& motif,
    int frame_offset,
    int receptor_nt_len,
    GenModel& model
);
```

Implementation:

```cpp
AAPgenResult compute_aa_pgen(const std::string& aa_seq, int frame_offset,
                              int receptor_nt_len, GenModel& model)
{
    // 1. Build JournaledQuery (parses motif string; plain AA sequences work unchanged)
    JournaledQuery jq = EventUtils::motif_to_journaled_query(
        motif, frame_offset, receptor_nt_len);

    // 2. Align IUPAC union reference against all gene templates
    std::string iupac_nt_str = EventUtils::int_str_to_nt(jq.iupac_union);
    auto gene_alignments = model.align_all_genes(iupac_nt_str, ...);

    // 3. Construct QuerySequenceContext in patched/motif mode
    QuerySequenceContext query(iupac_nt_str, jq.iupac_union, gene_alignments, jq);

    // 4. Build AccumulationContext with Single_error_rate.
    //    rate=0  → exact AA/motif Pgen (reject any patch-level mismatch)
    //    rate>0  → fuzzy Pgen (penalise patch mismatches, useful for hypermutation)
    auto pgen_error_rate = std::make_shared<Single_error_rate>(0.0);
    int n_patches = static_cast<int>(jq.patches.size());
    pgen_error_rate->build_upper_bound_matrix(n_patches + 1, n_patches + 1);
    AccumulationContext accumulation(updated_marginals, counters, pgen_error_rate);

    // 5. Build ExplorationContext with pruning_mismatch_floor (Phase 4)
    Pruning_mismatch_floor_map floor_map(6);
    ExplorationContext exploration(downstream_proba_map, seq_max_prob, threshold,
                                   index_map, next_event_ptr_arr, safety_set, floor_map);

    // 6. Run iterate chain (same as standard NT inference)
    model_queue[0]->iterate(query, model_ctx, scenario_ctx, exploration, accumulation);

    return { pgen_error_rate->get_seq_likelihood(), pgen_error_rate->debug_number_scenarios };
}
```

**Dependency**: Phases 3, 4, 5, 6.

---

### Phase 8 — Testing

#### 8.1 Unit tests: new utilities

| Test | Assertion |
|------|-----------|
| `translate_codon(int_A, int_T, int_G)` | `'M'` (Met) |
| `translate_codon(int_T, int_A, int_A)` | `'*'` (Stop) |
| `codon_mask_for_aa('L')` | bits set for TTA,TTG,CTT,CTC,CTA,CTG only (6 bits) |
| `codon_mask_for_aa('M')` | exactly one bit (ATG) |
| `motif_char_to_mask('X')` | 61 bits set (all non-stop codons) |
| `motif_char_to_mask('!')` | throws `std::invalid_argument` |
| `parse_aa_motif("ML")` | vector length 2; masks[0]=ATG only; masks[1]=6 Leu codons |
| `parse_aa_motif("CAV[KSM]DS")` | vector length 6; masks[3] = mask_K\|mask_S\|mask_M |
| `parse_aa_motif("CAVX")` | vector length 4; masks[3] = 61-bit wildcard mask |
| `parse_aa_motif("[")` | throws `std::invalid_argument` (unclosed bracket) |
| `parse_aa_motif("[]")` | throws `std::invalid_argument` (empty bracket) |
| `translate_int_seq(int_str("ATGTTA"), 0)` | `"ML"` (still useful for diagnostics) |

**`JournaledQuery` construction tests**:

| Test | Assertion |
|------|-----------|
| `motif_to_journaled_query("M", 0, 3)` | `patches.size() == 0` (single codon, one alternative → trivial) |
| `motif_to_journaled_query("L", 0, 3)` | `patches.size() == 1`; `patch.alternatives.size() == 5`; `reference` ∈ {TTA,TTG,CTT,CTC,CTA,CTG} |
| Leu `iupac_union[0]` | `int_Y` (T or C = union of {T} and {C}) |
| Leu `iupac_union[1]` | `int_T` (only T at position 1 of all Leu codons) |
| Leu `iupac_union[2]` | `int_N` (A,G,T,C all appear at position 2) |
| Leu `iupac_intersection[0]` | placeholder (don't care); `empty_isect[0] == true` |
| Leu `iupac_intersection[1]` | `int_T`; `empty_isect[1] == false` |
| Met `iupac_union == iupac_intersection == reference` | all three fields equal `{int_A, int_T, int_G}` |
| Met `empty_isect` | all false |
| `motif_to_journaled_query("[KS]", 0, 3)` | `patches.size() == 1`; alternatives covers all K+S codons minus reference; `iupac_union` and `iupac_intersection` computed from K∪S codon set |

#### 8.2 Mismatch semantics regression tests

For exact NT queries (all codes in [0..3]):
- `pruning_mismatch_floor[seg]` must equal `mismatches_lists[seg]` at every iterate step
  (floor = upper bound when query has no IUPAC codes)
- Full inference run on existing test dataset must produce numerically identical marginals
  before and after Phase 4 changes (no behaviour change for exact NT queries)

For IUPAC NT queries:
- At any position with an IUPAC query code: floor must not flag the position if
  `comp_nt_int` succeeds; `mismatches_lists` may flag it if `is_exact_match` fails
- `floor_vec.size() <= upper_vec.size()` at every node in the iterate tree
- No scenario that would contribute to NT Pgen is incorrectly pruned

For motif (AA) queries:
- `pruning_mismatch_floor[seg]` ⊆ `mismatches_lists[seg]` (floor uses `iupac_union`;
  upper bound uses `iupac_intersection` + `empty_isect` — so upper bound is at least as
  large as floor, as required)
- No scenario that would contribute to AA Pgen is incorrectly pruned by the floor

#### 8.3 `Single_error_rate` patch-mode unit tests

**rate = 0 (exact Pgen):**
- Single scenario, no insertions, assembled sequence matches `reference` for all patches:
  `pgen == scenario_proba`
- Single scenario, no insertions, one assembled codon not in `{reference_span} ∪ alternatives`:
  `pgen == 0.0`
- Single scenario, assembled codon equals an `alternative` (not `reference`): `pgen == scenario_proba`
  (alternative is a valid encoding → no mismatch)
- Two scenarios encoding the same motif (different V genes): `pgen == proba_1 + proba_2`
- Pruning threshold: scenario with `proba < seq_max * threshold` → not accumulated

**rate > 0 (fuzzy Pgen):**
- Single scenario, one patch mismatch: `pgen == scenario_proba * err_factor(1, n_patches)`
- Additivity: `pgen(motif, rate=r)` equals sum of `pgen(single_aa, rate=r)` for all AAs in
  `[KSM]` — holds because patch mismatch count is at codon level and AA sets are disjoint
- Verify `build_upper_bound_matrix` called with `patches.size()`, not nt length

#### 8.4 `Dinucl_markov` AA Pgen mode unit tests

- Insertion of length 3 fully within VD junction, target AA = Met (ATG only):
  `proba_contribution == T[prev][A] * T[A][T] * T[T][G]`
- Insertion of length 3, target AA = Leu:
  `proba_contribution == sum over {TTA,TTG,CTT,CTC,CTA,CTG} of T[prev][n0]*T[n0][n1]*T[n1][n2]`
- Codon straddling V/VD boundary (first 2 positions from V gene, third from insertion):
  verify only the free insertion position is summed

#### 8.5 Integration correctness test

For any AA sequence where all codons are unambiguous (Met = ATG, Trp = TGG only):

$$P_\text{gen}^\text{AA}(\text{seq}) = P_\text{gen}^\text{NT}(\text{unique encoding})$$

Both values are computable; assert numeric equality to tolerance $10^{-12}$.

For a sequence with one Leu position:

$$P_\text{gen}^\text{AA}(\text{...L...}) = \sum_{c \in \{TTA,TTG,CTT,CTC,CTA,CTG\}} P_\text{gen}^\text{NT}(\text{...}c\text{...})$$

Run NT inference six times (once per codon), sum results, compare to single AA Pgen call.

For a motif with one `[KS]` position (subset additivity):

$$P_\text{gen}^\text{AA}(\text{...[KS]...}) = P_\text{gen}^\text{AA}(\text{...K...}) + P_\text{gen}^\text{AA}(\text{...S...})$$

Run two single-AA Pgen calls, sum results, compare to the single bracket-group call. The
additivity holds because the codon sets for K and S are disjoint.

For a wildcard `X` at one position:

$$P_\text{gen}^\text{AA}(\text{...X...}) = \sum_{\text{aa} \in \text{all 20 AAs}} P_\text{gen}^\text{AA}(\text{...aa...})$$

Run 20 single-AA Pgen calls, sum results, compare to wildcard call. This also validates that
stop codons are excluded from the wildcard mask.

---

## 7. Dependency Graph

```
Phase 1
  (D-gene comp_nt_int fix; genetic code utilities:
   translate_codon, codon_mask_for_aa, [mask_to_iupac_codon],
   parse_aa_motif)
         │
         ▼
Phase 2 (JournaledQuery + Patch + motif_to_journaled_query factory;
         CodonMask used internally only — not exposed at runtime)
         │
         ▼
Phase 3 (QuerySequenceContext extension — adds journaled_query;
         int_sequence = iupac_union)
         │
         ├────────────────────────┐
         ▼                        ▼
Phase 4                     Phase 6
(mode-agnostic two-track     (Dinucl_markov patched mode —
 mismatch recording;          Markov forward sum over
 iupac_union for floor;       patch alternatives)
 iupac_intersection +
 empty_isect for upper bound)
         │                        │
         ▼                        │
Phase 5                           │
(Single_error_rate                │
 patch-based leaf check)          │
         │                        │
         └────────────────────────┘
                        │
                        ▼
                  Phase 7 (API)
                        │
                        ▼
                  Phase 8 (testing)
```

**Phase 1 is independent** and can be merged as a standalone correctness improvement
(D-gene `comp_nt_int` fix — though Phase 4 also handles the D-gene fix in the same
two-track loop). Phases 2, 3, 4, 5, 6 constitute the core patched/motif Pgen feature and
should be merged together. Phase 7 provides the user-facing entry point.

---

## 8. Unchanged Components

The following are explicitly not modified:

| Component | Reason |
|-----------|--------|
| `Insertion::iterate()` | Insertion count is determined by offsets; no nucleotide comparison |
| `Safety_bool_map` / overlap checks | Gene placement logic unchanged |
| `ScenarioContext` structure (except `mismatches_lists` upper-bound semantics) | No new counter or sequence fields |
| Legacy iterate interface | Passes through to context-based interface unchanged |
| `HypermutationfullNmererrorrate` | Not relevant for AA Pgen mode |
| Marginal update logic | Pgen computation does not update model marginals |

---

## 9. Notes on Insertion Probability in Patched Pgen Mode

The `Dinucl_markov` sum (Phase 6) produces a total probability weight $W$ over all insertion
nucleotide sequences that satisfy the patch alternatives touching the insertion span. $W$ is
multiplied into `scenario.scenario_proba`, and a representative compatible insertion is stored
in `constructed_sequences`.

`Single_error_rate` (Phase 5) then:
1. Reads `mismatches_lists` — the upper-bound mismatch positions pre-computed during traversal.
   These positions are not re-evaluated at the leaf.
2. For each patch where no position in `mismatches_lists` is set (fast path passes), assembles
   the concrete NT subsequence for that patch span and checks whether it matches `reference`
   or any entry in `patch.alternatives`.
3. For the representative insertion stored by `Dinucl_markov`: because the insertion was
   chosen to be patch-compatible, the leaf check for insertion-spanning patches passes (zero
   additional patch mismatches from insertions).
4. Applies the error model to the total patch mismatch count.

This is correct because:
1. The Markov model generates nucleotide sequences, not amino acid sequences.
2. $W = \sum_{\text{nt seq in patch set}} P(\text{nt seq} | \text{Markov model})$.
3. Upper-bound mismatch positions (step 1) capture all patch spans containing an
   NT incompatible with `iupac_intersection` (or flagged by `empty_isect`).
4. The patch-level check (step 2) catches the remaining cases where all NTs are
   intersection-compatible but the assembled NT subsequence is not in
   `{reference_span} ∪ alternatives`. This is the generalization of the Leu/Arg/Ser
   IUPAC over-approximation to arbitrary patch lengths. Full sequence assembly is
   unavoidable for patches that span segment boundaries.
5. The total accumulated probability is $P_\text{gen}^\text{patched}(\text{motif})$ for
   `error_rate = 0`, and the Pgen of sequences close to the motif for `error_rate > 0`.

For non-zero error rates, each mismatched patch contributes one error event regardless of
how many nucleotides within it differ. This is a deliberate simplification; finer
nucleotide-level weighting within mismatched patches is a possible future extension.

---
---

# Part II — Implementation Status Review (2026-08-27)

Scope of review: all non-merge commits on `feature/AA_PGEN` after `f393974f251`, i.e.

```
35440f9 first step                                   (Phases 1 + 2 + test scaffold)
f8e655a Phase 4 — two-track mismatch recording
83ac46b Phase 5 — patch-aware leaf check
01c89e3 Phase 6 — Dinucl_markov AA Pgen mode
2a7ecc6 override keyword sweep
818879d Fix three bugs in Dinuclmarkov AA-level Pgen
6e853ff Phase 7 — public API (compute_aa_pgen)
6d31346 fix regression
9bb37df simplify code
11c42d0 fix(Deletion): two-track mismatch for pruning bounds
059a27e fix(Genechoice): OpenMP data races
8474394 fix types
a1ac166 fix(Utils): heap-buffer-overflow in Enum_fast_memory_map::set_value
c49cb1f fix(Deletion): two-track mismatch for V_gene, D_gene, J_gene
e8c33a0 revert vj_length_d_position_proba and fix no_d_mismatches
819ff51 fix regression
88cf6d3 fix Windows build: popcountll
f4def2e fix test on osx
88a8dfa gene alignment using Aligner in compute_aa_pgen + more tests
```

Merge commits (`54c6198`, `c464d80`, `d517131`, `b322f60`, `286c94d`) pull in `develop` /
`feature/override` / pixi-env work and are excluded, except where a revert on this branch
(`e8c33a0`) undid part of a fix that had been made here.

Net footprint: **25 files, ~3900 insertions**, of which ~1270 lines are the two new test files.

---

## 10. Phase-by-phase status

| Phase | Status | Summary |
|---|---|---|
| 1 — Utility foundation | ✅ **Complete** | All declared functions implemented; D-gene `comp_nt_int` fix folded into Phase 4 as anticipated |
| 2 — `JournaledQuery` + factory | ✅ **Complete** (minor deviations) | Factory lives in a new `JournaledQuery.cpp`, not `EventUtils.cpp`; no bounds checking |
| 3 — `QuerySequenceContext` | ✅ **Complete** | Implemented exactly as specified |
| 4 — Two-track mismatch | ⚠️ **Structurally landed, functionally incomplete** | Floor track is written but never read; not populated on the aligner-derived paths |
| 5 — Patch-aware leaf check | ⚠️ **Landed, coordinate bug** | Logic matches design; patch coordinates are not mapped into assembled-sequence coordinates |
| 6 — `Dinucl_markov` AA mode | ❌ **Landed but incorrect** | Over-counts (Leu/Arg/Ser and all straddling codons); DJ junction direction wrong |
| 7 — Public API | ❌ **Landed but non-functional** | Uninitialised `Aligner`, IUPAC degraded to `N`, pruning threshold set so that everything is pruned |
| 8 — Testing | ⚠️ **§8.1 done, §8.2–§8.5 effectively absent** | Phases 5/6/7 are tested against test-local re-implementations, not production code |

---

### 10.1 Phase 1 — Utility foundation ✅

**Files**: `EventUtils.h` (+132), `EventUtils.cpp` (+247).

Implemented as designed: `codon_index`, `CodonMask`, `translate_codon`,
`codon_mask_for_aa`, `translate_int_seq`, `aa_to_iupac_codon`, `mask_to_iupac_codon`,
`motif_char_to_mask`, `parse_aa_motif`.

**Deviations**
- **Extra exported symbol**: `EventUtils::iupac_from_bits(int bits)` was added to the public
  header. The plan treated bit→IUPAC encoding as a private detail of the factory. It is now
  part of the DLL surface (`CORE_EXPORT`).
- `motif_char_to_mask('*')` returns the three stop codons rather than throwing. The plan only
  specified single-letter AAs and `'X'`, with everything else throwing. Harmless, but it means
  a motif containing `*` is silently accepted and will produce a `JournaledQuery` whose codon
  set is disjoint from any productive receptor.
- **Phase 1.1 (D-gene `!=` → `comp_nt_int`)** was *not* done as a standalone change; it was
  absorbed into the Phase 4 two-track loop in `Genechoice.cpp` (`no_d_mismatches`), as the
  plan's own note anticipated. It was then broken and re-fixed by `e8c33a0`
  ("fix `no_d_mismatches`") — worth a second look, see §10.9.

---

### 10.2 Phase 2 — `JournaledQuery` and factory ✅

**Files**: `JournaledQuery.h` (new, 106), `JournaledQuery.cpp` (new, 153),
`CMakeLists.txt` (+2).

`Patch` and `JournaledQuery` match Section 4.2 field-for-field.
`motif_to_journaled_query()` implements steps 1–4 of Section 2 faithfully, including the
uniform empty-intersection handling for Leu/Arg/Ser and the "omit the patch when the codon has
a single encoding" rule.

**Deviations**
- **Factory placement**: the plan put `motif_to_journaled_query` in `EventUtils.cpp`. It is in
  a dedicated `JournaledQuery.cpp` instead, still in `namespace EventUtils` and still declared
  in `JournaledQuery.h`. Cosmetic, but note that `JournaledQuery.h` now declares functions
  (the factory) as well as types, so it is no longer a pure data header.
- **No bounds checking against `receptor_len`.** The loop writes at
  `frame_offset + 3*codon_idx + j` with no check that this is `< receptor_len`. A motif longer
  than `(receptor_len - frame_offset)/3` writes past the end of `reference`, `iupac_union`,
  `iupac_intersection` and `empty_isect`. `compute_aa_pgen` does not validate the relation
  either, so this is reachable from the public API.
- **`reference` is not a valid NT sequence outside the motif span.** Non-codon positions stay
  `int_N`. Section 4.2 specifies `reference` as "one concrete valid NT sequence". Two
  downstream consumers assume this: `Single_error_rate`'s slow-path `ref_span` comparison
  (§10.5) and `Dinucl_markov`'s representative-insertion copy (§10.6) — both read
  `jq.reference` and both will read `int_N` if the motif does not tile the whole receptor.

---

### 10.3 Phase 3 — `QuerySequenceContext` extension ✅

**File**: `QuerySequenceContext.h` (+43). Implemented exactly as specified: `optional`
member, two constructors, copy/move still deleted.

One structural note for the refactor branch (see §13): `journaled_query` is the **only owning
member** of a struct that is otherwise entirely `const&`. It is held by value and moved in.

---

### 10.4 Phase 4 — Two-track mismatch recording ⚠️

**Files**: `Utils.h` (+typedef), `ExplorationContext.h` (+12), `Rec_Event.h`/`.cpp` (+14),
`Genechoice.cpp`/`.h`, `Deletion.cpp` (+148)/`.h`, `GenModel.cpp`, `test_Contexts.cpp`.

**What landed**
- `Pruning_mismatch_floor_map` typedef (`Utils.h:488`).
- `ExplorationContext::pruning_mismatch_floor` reference member + 7th constructor parameter.
- Threading through the **legacy** `Rec_Event::iterate()` / `iterate_wrap_up()` overloads
  (`Rec_Event.h:329`, `Rec_Event.cpp:139/174/221/245`) — see §13.1.
- The two-track loop, exactly in the shape of Section 4.3, in:
  - `Genechoice.cpp:622–668` and `748–835` (the two D-gene exhaustive-search branches);
  - `Deletion.cpp` negative-deletion (palindrome) loops for V (`~372–415`),
    D 5′ (`~630–662`), D 3′ (`~885–915`), J (`~1163–1203`);
  - lock-step trimming of `floor_mismatches_vector` alongside `mismatches_vector` in every
    positive-deletion branch (added incrementally by `11c42d0` and `c49cb1f`).
- New scratch members `Deletion::floor_mismatches_vector`, `Gene_choice::no_d_floor_mismatches`.

**Deviations and gaps — this is the most consequential group**

1. **The floor track is write-only. Nothing consumes it.**
   Section 4.3 states "`downstream_proba_map` entries use `floor_vec.size()`". They do not.
   Every `get_err_rate_upper_bound(...)` call in `Genechoice.cpp` (lines 340, 507, 661, 829,
   1007) and `Deletion.cpp` (lines 435, 681, 932, 1223) still passes the **upper-bound**
   count (`mismatches_vector.size()` / `endogeneous_mismatches` derived from it).
   A grep for `pruning_mismatch_floor` outside `Deletion.cpp`/`Genechoice.cpp` finds only the
   typedef, the context member, and the plumbing.

   Consequence: in motif mode the pruning bound is computed from the *larger* of the two
   counts, so the bound is **anti-conservative in the pruning direction** — scenarios that do
   contribute to AA Pgen can be discarded. This is precisely the failure mode §8.2 says must
   be tested against ("No scenario that would contribute to AA Pgen is incorrectly pruned by
   the floor"), and no such test exists.

2. **`Gene_choice` never populates the floor on the aligner-derived paths.**
   `set_mismatches(V_gene_seq, &(*iter).mismatches, …)` at `Genechoice.cpp:312`, and the
   equivalents at `:467` (D, aligned) and `:979` (J), install the aligner's mismatch list.
   No matching `pruning_mismatch_floor.set_value(...)` accompanies them. Only the two D-gene
   *exhaustive* branches call `set_value` (`:667`, `:834`).

   `Deletion::iterate` then does, whenever `journaled_query.has_value()`:
   ```cpp
   exploration.pruning_mismatch_floor.at(V_gene_seq, memory_layer_mismatches - 1)   // :247
   ```
   (and `:517`, `:775`, `:1053` for D 5′, D 3′, J). In motif mode this reads a slot that was
   never written — an uninitialised `std::vector<int>*` — and dereferences it. **Latent
   use-of-uninitialised-pointer / crash on the first real motif run.** It has not been
   observed because no test drives `compute_aa_pgen` end to end.

3. **Aligner-derived mismatch lists carry *floor* semantics but are stored in the *upper-bound*
   track.** `sw_align()` builds `Alignment_data::mismatches` with `comp_nt_int`, i.e. exactly
   the floor rule. Those lists are what `Gene_choice` puts into `scenario.mismatches_lists`.
   In motif mode, therefore, the upper bound on all aligner-derived gene positions is
   **under-counted** — the very approximation Section 3.2 identifies as "incorrect for leaf
   accumulation". Section 4 assumed every gene position would go through the two-track loop;
   in reality only the palindromic (negative-deletion) nucleotides and the D exhaustive path do.

   The plan has **no phase** covering this. Either the aligner must produce both tracks, or
   `Gene_choice` must re-run the two-track comparison over the aligned span, or the leaf check
   must stop relying on `mismatches_lists` for the fast path (see §10.5). This is new work.

4. **NT-mode floor is a per-realization copy on the hot path.**
   ```cpp
   vector<int> v_floor_mismatch_storage;
   vector<int>* v_floor_mismatch_ptr = query.journaled_query.has_value()
       ? exploration.pruning_mismatch_floor.at(...)
       : (v_floor_mismatch_storage = v_mismatch_list, std::addressof(v_floor_mismatch_storage));
   ```
   (`Deletion.cpp:243–249`, repeated four times.) In NT mode this copies the mismatch vector
   once per `Deletion::iterate` call, plus `floor_mismatches_vector` is assigned/trimmed in
   lock-step throughout. Section 4.1 promised "no behaviour change for existing workflows";
   that is true numerically but not for cost. Worth measuring before merging into `develop`.

5. **`Rec_Event` legacy interface changed.** Section 8 lists "Legacy iterate interface" under
   *Unchanged Components*. In fact `Pruning_mismatch_floor_map&` was appended to both legacy
   overloads. See §13.1.

---

### 10.5 Phase 5 — Patch-aware leaf check ⚠️

**File**: `Singleerrorrate.cpp:131–…`.

Implemented in the shape of Section 5: patch/motif branch collects upper-bound positions into
an `unordered_set`, assembles the full sequence, then fast-path / slow-path per patch. NT
branch preserved verbatim (re-indented into an `else`). `number_errors` = patch mismatch
count, `genomic_nucl` = `jq.patches.size()`, and the existing
`(rate/3)^n_err * (1-rate)^(n - n_err)` formula is reused unchanged — consistent with the
plan's "no dedicated `AA_pgen_error_rate` subclass".

**Deviations**

1. **Coordinate-space bug (blocking).** `patch.start` is a position in *receptor / query*
   coordinates (`frame_offset + 3i`). `EventUtils::build_scenario_sequence()` returns a plain
   concatenation `V ++ [VD] ++ D ++ [DJ] ++ J`, whose index 0 is the V gene's **5′ offset in
   the query**, not query position 0. The code indexes `full_seq` directly with `patch.start`
   and `jq.reference` with the same index:
   ```cpp
   Int_Str assembled(full_seq.begin() + patch.start, full_seq.begin() + patch.start + patch.length);
   ```
   No `scenario.get_offset(V_gene_seq, Five_prime)` correction is applied. Correct only when
   that offset happens to be 0. Also unguarded against `patch.start + patch.length >
   full_seq.size()` → out-of-range iterator.

2. **The fast path inherits gap (3) from §10.4.** It tests membership in `mismatches_lists`,
   which on aligner-derived positions holds floor-semantics data. So the fast path
   under-triggers, and correctness then rests entirely on the slow path — which is exactly the
   full-assembly comparison, so the *result* may still be right, but the "fast path" is not
   doing the work the design assigned it and the `upper_mismatches` set is not the object
   Section 4.1 defines.

3. **Segment presence is probed through `mismatches_lists`.**
   `has_vd_ins = scenario.mismatches_lists.exist(VD_ins_seq)` etc. Insertion segments only
   appear in that map if some event registered them. This couples "does the scenario have a
   VD junction" to "did anything record mismatches for VD". Fragile, and it is the same probe
   used to decide the layout of the assembled sequence in (1).

4. `exploration` is now an unused parameter in the patch branch except for the final threshold
   check — fine, but note the threshold check still runs in motif mode and interacts with the
   `proba_threshold_factor` problem described in §10.7.

---

### 10.6 Phase 6 — `Dinucl_markov` AA Pgen mode ❌

**Files**: `Dinuclmarkov.h` (+20), `Dinuclmarkov.cpp` (+351).

Landed: `compute_markov_aa_sum()` (forward DP, `:637`), `iterate_patched_pgen()` (`:737`),
dispatch from `iterate()` (`:94–97`). Three bugs were found and fixed in `818879d`
(DP base-case swap, codon-mask decoding, patch-index confusion) with narrative regression
tests. The following are **still open**.

1. **The DP loses intra-codon correlation — it reproduces the YTN obstacle.**
   `compute_markov_aa_sum` reduces each codon mask to three independent per-position
   nucleotide sets (`pos_allowed[pos]`), then runs a position-wise 4-state DP. For Leu that
   gives `pos0 ∈ {C,T}`, `pos1 = {T}`, `pos2 ∈ {A,C,G,T}` — i.e. **YTN**, which includes
   TTT and TTC (Phe). The sum therefore over-counts Leu, Arg and Ser insertions.
   Section 2 exists specifically to rule this out. The DP state must be widened to
   `(codon-position, partial codon prefix)` — 1 + 4 + 16 states per codon — or the codon set
   must be enumerated explicitly (≤ 6 terms, as Section 6.1 originally proposed).

2. **Gene-fixed positions of straddling codons are not constrained.** Section 6.1 requires
   that for a codon spanning the gene/insertion boundary, the positions fixed by the template
   pin the codon and only the free positions are summed. The implementation builds
   `codon_masks` purely from the patch/reference and never intersects with the actual
   V/D/J nucleotides in `scenario.constructed_sequences`. Every junction has at least one
   straddling codon in general, so this over-counts almost always.

3. **DJ junction direction is wrong.** NT mode (`Dinuclmarkov.cpp:135–150`) takes
   `previous_seq = J_gene_seq`, `previous_nt_str = j_seq.front()`, reverses the data substring,
   runs `iterate_common`, then reverses `dj_seq` back — the DJ chain runs 3′→5′ from the J side.
   `iterate_patched_pgen` instead uses
   ```cpp
   int prev_nt = d_gene_seq.front();          // D gene's FIRST nucleotide
   int d_end_pos = j_start - dj_len;
   ```
   with no reversal and no reversal of `rep_ins`. Both the context nucleotide and the chain
   direction are wrong; the comment `// Use D gene end as context` contradicts the code.

4. **Frame arithmetic ignores `frame_offset`.**
   `frame_start_in_ins = v_end_pos % 3` and the codon lookup searches for
   `patch.start == ci*3` where `ci = v_end_pos/3`. Patches actually live at
   `frame_offset + 3i`. Both are correct only when `frame_offset % 3 == 0`. The API accepts an
   arbitrary `frame_offset`.

5. **Representative insertion may be incompatible.** `rep_ins` is copied straight out of
   `jq.reference` over the insertion span. Section 9 step 3 assumes the stored insertion is
   patch-compatible so the leaf check sees zero extra mismatches. Since `reference` is chosen
   independently of the gene-fixed neighbours (and is `int_N` outside the motif span,
   §10.2), a straddling patch can fail the leaf check — which at `rate = 0` zeroes the entire
   scenario. Silent, total loss of probability mass.

6. **`model_params` is accepted and ignored** — `compute_markov_aa_sum` reads the member
   `dinuc_proba_matrix`. No caching table was built, contrary to the §6.2 suggestion.

7. **Three near-identical ~90-line blocks** (VD / DJ / VJ) do the same codon-mask lookup.
   `818879d` had to apply the same fix three times, which is exactly the maintenance hazard
   that matters when the refactor branch lands (§13.5).

8. **`iterate_patched_pgen` does not populate `*_realizations_indices` / the write index
   list**, so marginals and counters are not driven from this path. Acceptable for Pgen (§8
   "Marginal update logic … unchanged"), but it means the patched path cannot be reused for
   inference, and any refactor that assumes every `iterate()` fills the index list will trip.

---

### 10.7 Phase 7 — Public API ❌

**Files**: `GenModel.h` (+30), `GenModel.cpp` (+~230).

Shipped as a **`GenModel` member** `GenModel::compute_aa_pgen(motif, frame_offset,
receptor_nt_len)` plus a free-function wrapper `compute_aa_pgen(..., GenModel&)` declared in
`GenModel.h` (the plan specified only the free function). `AAPgenResult { double pgen; size_t
n_scenarios; }` matches, with `debug_number_scenarios` widened to `size_t` in `Errorrate.h` to
suit.

The eight-step body follows the plan's skeleton. It is nonetheless **not usable as written**:

1. **Alignment uses default-constructed `Aligner` objects.** `88a8dfa` added
   ```cpp
   Aligner v_aligner;                       // Aligner::Aligner() is an empty stub
   v_aligner.set_genomic_sequences(get_templates(v_gene_choice));
   auto v_alignments = v_aligner.align_seq(iupac_nt_str, -10.0, false, 0, iupac_nt_str.size(), false);
   ```
   `Aligner::Aligner()` in `Aligner.cpp` has an empty body (`// TODO Auto-generated
   constructor stub`) — the substitution matrix, gap penalty and `gene_class` are
   uninitialised. `gene_class` in particular decides the anchoring semantics (V anchored 3′,
   J anchored 5′). Section 7.1 assumed the *model's configured* aligner would be reused.
   The accompanying tests (`test_aa_pgen_bugs.cpp` #4–#6) only assert `num_alignments > 0`,
   which an uninitialised aligner can still satisfy.

2. **IUPAC is degraded to `N` before alignment.** The `Int_Str → std::string` conversion maps
   only `A/C/G/T/N` and sends every other code (`R Y K M S W B D H V`) to `'N'`:
   ```cpp
   case 14: iupac_nt_str += 'N'; break;
   default: iupac_nt_str += 'N'; break;
   ```
   So the string that is actually aligned has no degeneracy information, and every mismatch
   list produced by the aligner is empty at exactly the interesting positions.
   `query.int_sequence` still holds the true `iupac_union`, so `query.sequence` and
   `query.int_sequence` **disagree**. Section 7.3 step 2 called for a real `int_str_to_nt`
   helper; it was never written. `Aligner::nt2int` already parses all 15 codes, so the inverse
   is a 15-entry table.

3. **Pruning threshold makes the result identically zero.**
   ```cpp
   double max_proba_scenario = 1.0;
   double proba_threshold_factor = 0.001;
   ```
   `seq_max_prob_scenario` is never updated from a first pass, so `should_prune` compares every
   scenario upper bound against an absolute `1e-3`. Real scenario probabilities are many orders
   of magnitude below that, so essentially everything is pruned and `get_seq_likelihood()`
   returns 0. For exact Pgen the threshold must be 0, or a proper max-probability pre-pass must
   run (as normal inference does).

4. **`Marginal_array_p updated_marginals(new long double[1024], …)`** — fixed size, never
   initialised. Any event that writes marginals overflows it for a real model.

5. **`if (num_events == 0 || num_events > 100) return {0.0, 0ULL};`** — an arbitrary guard that
   returns a *valid-looking* zero Pgen.

6. **`catch (const std::exception&) { return {0.0, 0ULL}; }`** — swallows every failure, so a
   crash-free bug (e.g. the uninitialised floor pointer of §10.4.2, if it throws) is
   indistinguishable from a genuine zero. This should at minimum log, and probably rethrow.

7. **No `frame_offset` / `receptor_nt_len` validation** (see §10.2 OOB write), no CDR3-derived
   alignment offset bounds (Section 7.1 flagged this explicitly), and no CLI / `ExtractFeatures`
   exposure — `compute_aa_pgen` is not reachable from any command line.

---

### 10.8 Phase 8 — Testing ⚠️

Two new files: `tst/igor/Core/test_aa_pgen.cpp` (1119 lines, 38 `TEST_CASE`s) and
`tst/igor/Core/test_aa_pgen_bugs.cpp` (151 lines, 3 `TEST_CASE`s).

| Plan section | Status |
|---|---|
| §8.1 unit tests (genetic code + `JournaledQuery`) | ✅ Covered thoroughly, and beyond the plan (Arg, Ser, multi-codon, mid-sequence `frame_offset`, mixed bracket+wildcard motifs) |
| §8.2 mismatch semantics regression | ❌ **Absent** |
| §8.3 `Single_error_rate` patch mode | ❌ **Tests a test-local copy, not the production function** |
| §8.4 `Dinucl_markov` AA mode | ❌ **Does not exercise the code** |
| §8.5 integration correctness | ❌ **Absent** |

Detail:

- **§8.2.** The three `[phase4]` cases (`test_aa_pgen.cpp:577–670`) build expected floor/upper
  sets by hand from a `JournaledQuery`. They never call `Gene_choice::iterate` or
  `Deletion::iterate`. There is no `floor ⊆ upper` invariant check at iterate nodes, and no
  "marginals numerically identical before/after Phase 4" regression run on an existing dataset —
  which is the one test that would have justified the "no behaviour change" claim.

- **§8.3.** `test_aa_pgen.cpp` defines a file-local `check_patch(full_seq, reference, patch,
  upper_mismatches)` that **re-implements** the leaf logic, and the four `[phase5]` cases
  exercise that copy. `Single_error_rate::compute_scenario_error_probability` is never called
  from any test. Missing outright: all `rate > 0` cases, the `[KSM]` additivity case, and the
  `build_upper_bound_matrix(patches.size())` sizing assertion.

- **§8.4.** `compute_markov_aa_sum` is `private`, so the `[phase6]` cases construct a
  `Dinucl_markov` and then only assert properties of `CodonMask` bit patterns
  (`REQUIRE(mask == (1ULL << 58))`, `REQUIRE(bit_count == 6)`). One test even says so:
  *"We can't directly call compute_markov_aa_sum (private) … For now, just verify the codon
  mask is correct."* None of the plan's numeric assertions
  (`T[prev][A]·T[A][T]·T[T][G]`, the 6-term Leu sum, the straddling-boundary case) exist —
  which is why §10.6.1 and §10.6.2 went unnoticed.

- **§8.5.** `compute_aa_pgen` is **never invoked** anywhere in the test suite. None of the four
  integration identities (unambiguous-codon AA == NT Pgen; Leu == sum of 6 NT Pgens;
  `[KS]` additivity; `X` == sum over 20 AAs) is implemented. These are the tests that would
  catch §10.5.1, §10.6, and every item in §10.7.

- **Note on the suite's health.** Commit messages record a pre-existing inference test failure
  (SIGABRT/SIGSEGV) at `01c89e3`, subsequently addressed by `059a27e` / `a1ac166` /
  `e8c33a0` / `819ff51`. Re-confirm the full suite is green on this branch before merging.

---

### 10.9 Content added that the plan does not describe

Beyond the feature itself, this branch carries several changes that are independent of AA Pgen
and should be reviewed (and possibly split out) on their own merits:

| Change | Commit | Note |
|---|---|---|
| `popcountll()` portability shim (`Utils.h:109`, MSVC `__popcnt64`) | `88cf6d3` | Windows build fix; general-purpose utility |
| `Enum_fast_memory_map::set_value` buffer growth + relaxed first-call assertion | `a1ac166` | Fixes an ASan heap-buffer-overflow reached from `Deletion::iterate`. **Changes a memory-layer invariant** — see §13.8 |
| `Errorrate::debug_number_scenarios` `int` → `size_t` | `8474394` | Public member type change |
| `override` keyword sweep | `2a7ecc6` | From `feature/override` |
| `Gene_choice` OpenMP `#pragma omp single` around init; `D_position_info` storing `std::string` instead of a dangling `const char*`; incremental `vj_length_d_position_range` tracking | `059a27e` | **Largely reverted by `e8c33a0`** ("revert vj_length_d_position_proba"). Confirm this was intentional — the dangling-pointer and data-race fixes described in `059a27e` do not appear in the current `Genechoice.h`/`.cpp`. If the crash it fixed was real, it is back |
| `JournaledQuery.cpp` as a separate TU | `35440f9` | Plan put the factory in `EventUtils.cpp` |
| `EventUtils::iupac_from_bits` exported | `35440f9` | Plan treated this as private |
| "Documentation-style reproduction tests" | `818879d`, `88a8dfa` | A testing convention not in the plan: numbered `BUG #n` / `BUG_TEST #n` cases whose bodies narrate the defect. Useful as a record; note that #4–#6 test `Aligner` directly rather than `compute_aa_pgen`, so they do not actually cover the code that was changed |

---

## 11. Remaining work

Ordered by dependency. Items marked **(new)** are not in the original Sections 1–9.

### 11.1 Blocking correctness

1. **Define one coordinate convention and apply it everywhere (new).**
   Three coordinate spaces are currently in play and are silently mixed:
   *receptor coordinates* (`JournaledQuery`, `patch.start = frame_offset + 3i`);
   *query coordinates* (`scenario.seq_offsets`, `query.int_sequence` indices, all mismatch
   positions); *assembled-sequence coordinates* (`build_scenario_sequence`, index 0 = V 5′
   offset). Fix `Single_error_rate` (§10.5.1) and `Dinucl_markov` (§10.6.4) against that
   convention, and add the offset translation explicitly rather than by coincidence.
2. **Wire the floor track into the pruning bound** (§10.4.1) — Section 4.3 as originally
   written. Every `get_err_rate_upper_bound` call site in `Genechoice.cpp` / `Deletion.cpp`.
3. **Populate the floor on the aligner-derived V/D/J paths** (§10.4.2), or make the
   `Deletion` read tolerate an unset slot. Currently a latent uninitialised-pointer read on
   the first motif run.
4. **Give the upper-bound track correct semantics on aligned gene spans (new)** (§10.4.3) —
   either dual-track output from `sw_align`, or a re-scan in `Gene_choice`, or drop the
   fast path in `Single_error_rate` and rely on the slow path alone.
5. **Rewrite `compute_markov_aa_sum` to preserve intra-codon correlation** (§10.6.1) and to
   pin gene-fixed positions of straddling codons (§10.6.2).
6. **Fix the DJ junction direction and context nucleotide** (§10.6.3).
7. **Choose a patch-compatible representative insertion** rather than copying `jq.reference`
   (§10.6.5).
8. **Make `compute_aa_pgen` functional**: real aligner configuration (§10.7.1), lossless
   `int_str_to_nt` (§10.7.2), threshold 0 or a max-probability pre-pass (§10.7.3),
   correctly-sized marginal buffer (§10.7.4), remove the `num_events > 100` and blanket-catch
   guards (§10.7.5–6), validate `frame_offset`/`receptor_nt_len` (§10.2).

### 11.2 Tests that must exist before this is trustworthy

9. §8.5 integration identities — **all four**. These are the highest-value tests on the list;
   each one independently catches several of items 1–8.
10. §8.4 numeric Markov tests, including the straddling-codon case. Requires making
    `compute_markov_aa_sum` testable (friend, or move to a free function in `EventUtils`).
11. §8.3 against the real `Single_error_rate::compute_scenario_error_probability`; delete the
    test-local `check_patch` copy.
12. §8.2 regression: NT-mode marginals bit-identical pre/post Phase 4, plus a `floor ⊆ upper`
    invariant assertion at iterate nodes.
13. **(new)** A performance check on NT-mode `Deletion::iterate` for the copy introduced in
    §10.4.4.

### 11.3 Follow-ups

14. Deduplicate the three VD/DJ/VJ blocks in `iterate_patched_pgen` (§10.6.7).
15. CLI / `ExtractFeatures` exposure of `compute_aa_pgen`.
16. Resolve the status of the `059a27e` OpenMP / dangling-pointer fixes that `e8c33a0`
    reverted (§10.9).

---

## 12. How the deviations propagate into the unimplemented work

The parts that are still missing are *not* independent of the parts that landed. Concretely:

- **§8.5 integration tests cannot pass until items 1–8 are done.** They are currently written
  nowhere, which conceals that fact. Expect the first end-to-end run to fail for at least four
  distinct reasons at once (threshold-prunes-everything, uninitialised aligner, floor read of
  an unset slot, patch coordinate offset). Land items 3 and 8 first so the pipeline runs at
  all, then use the Met/Trp identity (unique encoding ⇒ AA Pgen == NT Pgen) as the first
  green light — it exercises the whole chain while keeping `patches` empty-ish and the Markov
  sum trivial, so it isolates §10.7 from §10.6.

- **The Leu identity test (§8.5) will fail on §10.6.1 even after everything else is fixed.**
  The per-position DP over YTN adds TTT/TTC weight whenever a Leu codon lands inside a
  junction. Any attempt to validate Phase 6 numerically before rewriting the DP will produce
  a confusing "close but wrong by a few percent" result. Fix the DP *before* writing the test,
  or the test will be tuned to the wrong number.

- **Fixing the floor→pruning wiring (item 2) will change NT-mode behaviour unless item 4 is
  done first.** In NT mode floor == upper, so switching the bound to the floor is a no-op *for
  positions that went through the two-track loop*. But if item 4 changes what
  `mismatches_lists` contains on aligned spans, NT-mode numbers move too. Sequence these two
  together and re-run the §8.2 regression across both.

- **`Single_error_rate`'s fast path depends on the semantics decided in item 4.** If the
  upper-bound track stops being aligner-derived, the fast path starts firing on positions where
  it previously did not, and `number_errors` (hence Pgen at `rate > 0`) changes. The §8.3
  tests must be written against the *post*-item-4 semantics, not the current ones.

- **The representative-insertion choice (item 7) is coupled to the leaf check (item 1).**
  Both read `jq.reference` at receptor coordinates. If item 1 introduces an offset translation
  in `Single_error_rate` but not in `Dinucl_markov`, the two will disagree and every
  insertion-spanning patch will read as a mismatch. Change them in the same commit.

- **`reference` containing `int_N` outside the motif span (§10.2) makes items 1 and 7 fragile.**
  Decide now whether `motif_to_journaled_query` must tile the whole receptor (then validate it)
  or whether `reference` is only meaningful within the motif span (then every consumer needs a
  span guard). The current code assumes the former and enforces neither.

- **`Dinucl_markov` counters/marginals (§10.6.8)** are unimplemented on the patched path. That
  is fine for Pgen but means the "fuzzy Pgen with `rate > 0`" use case in §5 — which the plan
  presents as a headline capability — cannot be combined with any counter
  (`Pgencounter`, `Coverageerrcounter`, `Bestscenarioscounter`). If that combination is wanted,
  it is additional unplanned work.

---

## 13. Merge watch-list against the `iterate()` refactor branch

Ordered by expected difficulty. Items 1–5 are where the two branches genuinely collide;
6–10 are smaller but easy to lose in a mechanical merge.

### 13.1 `Rec_Event` legacy `iterate()` / `iterate_wrap_up()` signatures — **highest risk**

Section 8 of this plan claimed the legacy interface was untouched. It is not:
`Pruning_mismatch_floor_map& pruning_mismatch_floor` was appended to the deprecated
overloads (`Rec_Event.h:329`, `Rec_Event.cpp:139`, `:174`, `:221`, `:245`), and to the legacy
`iterate()` declarations in `Dinuclmarkov.h`, `Genechoice.h`, `Deletion.h`, `Insertion.h`.

If the refactor branch **deletes** the legacy overloads (the likely direction), this parameter
should simply disappear: the floor map already lives in `ExplorationContext` and that is its
correct home. Resolve in favour of the refactor branch and drop the parameter — do **not**
mechanically re-add it. If the legacy overloads survive, the parameter must be threaded again.

### 13.2 `ExplorationContext` constructor arity

Gained a 7th parameter. Construction sites: `GenModel.cpp:459` (inference),
`GenModel.cpp:1151` (`compute_aa_pgen`), `tst/igor/Core/test_Contexts.cpp` (+18 lines of
churn). Any construction site the refactor branch adds or moves needs the same argument.
Consider making `pruning_mismatch_floor` a defaulted/optional member so the arity stops being
a merge tripwire.

### 13.3 `Deletion.cpp` — **highest textual conflict risk**

~148 changed lines, spread across all four gene/side sections, and interleaved *line by line*
with the existing mismatch-trimming logic rather than confined to new blocks. Four repeated
patterns to preserve:
- floor-list acquisition at `:243–249`, `:513–519`, `:771–777`, `:1049–1055`
  (each is a ternary on `journaled_query.has_value()`);
- `floor_mismatches_vector.clear()` / `.assign(...)` mirroring every
  `mismatches_vector.clear()` / `.assign(...)` in the positive-deletion branches;
- the two-track `floor_mismatch` / `upper_mismatch` block in each negative-deletion
  (palindrome) loop;
- `exploration.pruning_mismatch_floor.set_value(...)` paired with every
  `scenario.set_mismatches(...)` at `:415`, `:662`, `:915`, `:1203`.

If the refactor restructures the deletion loops, re-derive these from the pattern rather than
accepting hunks — a partially-applied merge here produces exactly the "floor set on some
paths, not others" bug already present in `Gene_choice` (§10.4.2), and it will be silent.

### 13.4 `Genechoice.cpp` D-gene exhaustive-search branches

Two two-track blocks at `:622–668` and `:748–835`. These sit inside the `no_d_align`
exhaustive path, which `e8c33a0` and `6d31346` already churned once on this branch (and where
`059a27e`'s OpenMP fixes were reverted). Expect conflicts on both sides. Note also that the
*aligned* V/D/J paths at `:312`, `:467`, `:979` deliberately have **no** floor call today
(§10.4.2) — if the merge adds one, that is a fix, not a regression, but it changes behaviour.

### 13.5 `Dinuclmarkov.cpp` — `iterate_patched_pgen` is a clone, not a specialisation

`iterate()` dispatches at `:94–97`; `iterate_patched_pgen` (`:737–966`) duplicates the
VD / DJ / VJ structure of `iterate()` without sharing `iterate_common()`. If the refactor
changes how `iterate()` obtains segments, offsets or the previous nucleotide — and the DJ
reversal convention in particular (`:135–150`) — **`iterate_patched_pgen` will not be updated
by the merge and will silently diverge further.** It already has the DJ direction wrong
(§10.6.3). Strongly consider refactoring it onto a shared `iterate_common`-style helper
*before* merging, so there is one place to fix.

### 13.6 `Singleerrorrate::compute_scenario_error_probability`

The entire pre-existing NT body was re-indented into an `else` branch. A three-way merge will
report the whole function as conflicting even where the NT logic is unchanged. Diff with
whitespace ignored when resolving.

### 13.7 `QuerySequenceContext` ownership

`journaled_query` is a by-value `std::optional<JournaledQuery>` inside a struct that is
otherwise all `const&`, with copy and move deleted. If the refactor branch makes query
contexts copyable, shareable across threads, or cheaply re-bindable per sequence, note that
this member is the only one with a lifetime of its own — and it is potentially large
(three `Int_Str`s of receptor length plus up to 61 alternatives per patch). Sharing it by
`shared_ptr<const JournaledQuery>` is probably the right end state.

### 13.8 `Utils.h` — `Enum_fast_memory_map::set_value` invariant change

`a1ac166` both grows the buffer on demand and **relaxes the layer assertion** so that the
first call (`memory_layer_ptr[key] == -1`) may set any layer. Any refactor that reasons about
memory-layer request/release ordering is reasoning about this invariant. Keep the change —
it fixes a real ASan overflow — but re-validate it against whatever layer discipline the
refactor introduces. `popcountll` in the same header is independent and safe.

### 13.9 Per-object scratch buffers and thread safety

`Deletion::floor_mismatches_vector` and `Gene_choice::no_d_floor_mismatches` are new mutable
members on the event objects, following the existing `mismatches_vector` / `no_d_mismatches`
pattern. Events are shared across the OpenMP parallel region, so these are **not thread-safe**
by the same argument that made `059a27e` necessary. If the refactor moves scenario scratch
state into a per-thread or per-scenario object, these two must move with it — and this is the
right moment to do it.

### 13.10 Smaller items

- `Errorrate::debug_number_scenarios` `int` → `size_t` (`Errorrate.h:172`) — touches anything
  that reads it.
- `GenModel.cpp` gained ~230 lines including the whole of `compute_aa_pgen`; `819ff51` also
  reworked ~48 lines of the inference setup path. Check the inference path specifically, not
  just the new function.
- `tst/igor/Core/CMakeLists.txt` registers two new test targets.
- `EventUtils.h` now includes `Rec_Event.h`; adding the genetic-code utilities to that header
  pulls the genetic code into every TU that includes `EventUtils.h`. If the refactor moves
  `EventUtils`, consider splitting the genetic-code half into its own header.
