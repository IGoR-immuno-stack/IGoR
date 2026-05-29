# AA-Level Pgen Implementation Plan

**Created**: May 2026  
**Status**: Design / Pre-implementation (revised: JournaledQuery design)  
**Relates to**: `ITERATE_METHOD_ANALYSIS.md`, context refactoring (Phases 1–5)

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
