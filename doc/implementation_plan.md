# Staged Plan: Gene_class and Seq_type Separation

This plan refactors the codebase to separate genomic identity (`Gene_class`, which determines genomic/alignment constraints) and sequence-segment identity (`Seq_type`, which specifies constructed sequence segments).

## Architecture Decision Record (ADR)

### Keying Choice in Model Registry (`events_map`)
- **Option A**: Keep current staged approach (no key change until much later, keep double identity on events).
- **Option B (Selected)**: Perform early key migration to `tuple<Event_type, Seq_type, Seq_side>`.
  - *Rationale*: Transitioning directly to this tuple is identifying, avoids a complex intermediate double-identity phase, and is a cleaner path towards tandem-D support. The bridge helpers added in Phase 1e isolate the callsites, making this early migration low-risk.

---

## Tandem-D Seq_type Variant Establishment

Before verifying uniqueness of the tuple key in Phase 2, explicitly establish the `Seq_type` enum variants required for tandem-D models:

1. **Codify new enum values**:
    - `D1_gene_seq` — first D segment in tandem-D arrangement.
    - `D2_gene_seq` — second D segment in tandem-D arrangement.
    - `D1D2_ins_seq` — insertion region between D1 and D2.
2. **Document positional semantics** in model/core headers so these variants are unambiguous relative to standard `V_gene_seq`, `D_gene_seq`, `J_gene_seq`, `VD_ins_seq`, and `DJ_ins_seq`.
3. **Align model I/O support** so internal representations and parsing logic can recognize these variants when tandem-D support is enabled.

This establishment is a prerequisite for the uniqueness proof below.

---

## Uniqueness Proof Checkpoint

We verify below that `tuple<Event_type, Seq_type, Seq_side>` is uniquely identifying for standard and tandem-D models.

### Standard VDJ Model
Standard events and their keys mapping:
*   V Choice: `(GeneChoice_t, V_gene_seq, Undefined_side)`
*   D Choice: `(GeneChoice_t, D_gene_seq, Undefined_side)`
*   J Choice: `(GeneChoice_t, J_gene_seq, Undefined_side)`
*   V 3' Deletion: `(Deletion_t, V_gene_seq, Three_prime)`
*   D 5' Deletion: `(Deletion_t, D_gene_seq, Five_prime)`
*   D 3' Deletion: `(Deletion_t, D_gene_seq, Three_prime)`
*   J 5' Deletion: `(Deletion_t, J_gene_seq, Five_prime)`
*   VD Insertion Length: `(Insertion_t, VD_ins_seq, Undefined_side)`
*   DJ Insertion Length: `(Insertion_t, DJ_ins_seq, Undefined_side)`
*   VD Dinucl Markov: `(Dinuclmarkov_t, VD_ins_seq, Five_prime)` (anchored on V 3' end)
*   DJ Dinucl Markov: `(Dinuclmarkov_t, DJ_ins_seq, Three_prime)` (anchored on J 5' end, traversal 3' to 5')

### Tandem D Model ($V - D_1 - D_2 - J$)
With multiple homologous D segments:
*   $D_1$ Choice: `(GeneChoice_t, D1_gene_seq, Undefined_side)`
*   $D_2$ Choice: `(GeneChoice_t, D2_gene_seq, Undefined_side)`
*   $D_1$ 5' Deletion: `(Deletion_t, D1_gene_seq, Five_prime)`
*   $D_1$ 3' Deletion: `(Deletion_t, D1_gene_seq, Three_prime)`
*   $D_2$ 5' Deletion: `(Deletion_t, D2_gene_seq, Five_prime)`
*   $D_2$ 3' Deletion: `(Deletion_t, D2_gene_seq, Three_prime)`
*   $D_1 D_2$ Insertion: `(Insertion_t, D1D2_ins_seq, Undefined_side)`
*   $D_1 D_2$ Dinucl: `(Dinuclmarkov_t, D1D2_ins_seq, Five_prime)` (anchored on $D_1$ 3' end)

*Conclusion*: `Seq_type` keying prevents all collisions.

---

## Event Key Signature Migration Specification

Lookup APIs in `EventUtils.h` will be updated to take `Seq_type` instead of `Gene_class` without keeping backward-compatibility overloads.

### `try_get_event`
*   **Pre-migration**:
    ```cpp
    bool try_get_event(
        const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Event_type event_type, Gene_class gene_class, Seq_side seq_side, std::shared_ptr<Rec_Event> &event_ptr);
    ```
*   **Post-migration**:
    ```cpp
    bool try_get_event(
        const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Event_type event_type, Seq_type seq_type, Seq_side seq_side, std::shared_ptr<Rec_Event> &event_ptr);
    ```

### `has_insertion_seq_type`
*   **Pre-migration**:
    ```cpp
    bool has_insertion_seq_type(
        const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Seq_type seq_type);
    ```
*   **Post-migration**:
    ```cpp
    bool has_insertion_seq_type(
        const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Seq_type seq_type);
    ```

---

## Phase 1: Semantic Separation in Event and Sampler Internals (Completed)
- Step 1.1: Codify semantic contract and bridge helpers. *(Done)*
- Step 1.2: Refactor fast generation in `FastGenerator.cpp`. *(Done)*
- Step 1.3: Refactor insertion internals in `Insertion.cpp`. *(Done)*
- Step 1.4: Refactor dinucleotide internals in `Dinuclmarkov.cpp`. *(Done)*
- Step 1.5: Refactor deletion internals in `Deletion.cpp`. *(Done)*
- Step 1.6: Refactor ancillary consumers (`Pgencounter.cpp`, `Coverageerrcounter.cpp`, etc.). *(Done)*
- Step 1.7: Add the typed `Seq_type` view in `Model_Parms` and migrate the first production consumer (`Coverage_err_counter`). *(Done)*
- Step 1.8: Add regression coverage for the typed bridge and helper lookups. *(Done)*

---

## Phase 2: Deletion Analysis, Model Key Migration & Dinucleotide Simplification (Active)

### Step 2.1: Deletion Identity Matrix (Completed)
Before executing model structural changes, perform a targeted analysis of every use of `Gene_class` in `Deletion.h`/`Deletion.cpp` to classify as:
- **Genomic Identity** (e.g. overlap safety checks such as `Event_safety::VD_safe`): Keep `Gene_class`.
- **Constructed-segment operations** (e.g. sequence offset retrieval, slicing, mismatch indexing): Use `Seq_type`.
- Write findings to a dedicated `deletion_identity_matrix.md` artifact.

### Step 2.3: Atomic Migration of `events_map` Key and Helper Signatures (In Progress)
This step must be executed as a coordinated atomic operation: `Model_Parms` storage format and `EventUtils` lookup signatures are mutually dependent and must transition together.

- Keep the legacy `events_map()` path available while exposing a typed `Seq_type` view for consumers that are ready to migrate.
- Simultaneously update `EventUtils::try_get_event` and `EventUtils::has_insertion_seq_type` signatures to consume the new key structure (replacing `Gene_class` lookup parameter with `Seq_type`).
- Update unit tests (`tst/igor/Core/test_EventUtils.cpp`, etc.) to verify `Seq_type`-keyed lookups end-to-end.
- Migrate consumers incrementally to the typed view and retire the legacy path only after the remaining call sites are gone.

#### Step 2.3b Consumer Migration Snapshot (Completed for event runtime paths)
- `Deletion`, `Insertion`, `Dinucl_markov`, and `Gene_choice` now expose typed `Seq_type` overloads for `iterate(...)` and `initialize_event(...)`; the legacy `initialize_event(...)` wrappers have been removed, and the `Rec_Event` bridge now centralizes the legacy-to-typed adaptation.
- `Insertion` and `Dinucl_markov` now also expose typed overloads for `initialize_crude_scenario_proba_bound(...)`, again with legacy wrappers.
- Internal lookup code in migrated typed paths now uses `EventUtils` typed helpers (`Seq_type` keys) instead of `Gene_class` tuple lookups.
- Unit validation remains green (`pixi run test_unit`: 31/31).

### Step 2.2: Simplify Dinucleotide Target to 2-Axis (Refinement)
Adjust `Dinuclmarkov` sampler to use a simplified 2-axis targeting:
- `target Seq_type`
- `anchor side`
- Traversal direction is derived implicitly. Remove the 3rd axis (`traversal side`) from the primary state.

---

## Phase 3: Serialization & Legacy Compatibility (Issue #46)
- Perform model file compatibility mapping inside `Model_Parms::read_model_parms`.
- Translate legacy files (using junction classes like `VD_genes`) into the new internal `Seq_type`-keyed events.
- Implement explicit serialization rules (deciding on write-path compatibility) in the scope of Issue #46.

### Phase 3 Acceptance Criteria
- Model files written in legacy format are successfully translated into the internal `Seq_type`-keyed event representation during load.
- Legacy loading is validated by regression and unit loading tests, including round-trip safety where applicable.
- Serialization strategy is explicit and documented: write-path compatibility behavior is defined (legacy vs new format) within Issue #46 scope.

---

## Phase 2 Acceptance Criteria
- Zero compile/link errors.
- 100% of unit tests pass successfully (`pixi run test_unit`).
- All runtime consumers that access `events_map` do so through `EventUtils` helpers (for example `try_get_event` and `has_insertion_seq_type`) and do not perform direct tuple-key lookup calls (`.at(...)`, `.count(...)`). This rule is class-agnostic and applies to future consumers as well.

## Current Migration Snapshot
- Typed `Seq_type` bridge is exposed through `Model_Parms::get_events_map_seq_type()`.
- `ModelContext` now carries a typed `Seq_type`-keyed `events_map` view.
- `GenModel` now constructs context pipelines with the typed view while preserving legacy wrappers for old call paths.
- `GenModel` runtime initialization and downstream probability-bound setup now dispatch through typed `Rec_Event` overloads; the legacy `Model_Parms::get_events_map()` getter is no longer used in the core runtime path.
- `Pgen_counter` and `Coverage_err_counter` context initialization paths now read from typed keys.
- `Error_rate` and hypermutation context bridges no longer rebuild legacy event maps when delegating to legacy compare signatures; those compare implementations are still compatibility-only, but the typed runtime path no longer pays that conversion cost.
- The deprecated `compare_sequences_error_prob(...)` compatibility interface no longer carries an unused legacy `events_map` parameter; both wrappers and remaining implementations now use the narrower signature.
- Hypermutation error-rate classes now route `compute_scenario_error_probability(...)` through internal context-first helpers, with `compare_sequences_error_prob(...)` retained only as a compatibility shim.
- `Rec_Event` legacy iterate wrapper now converts legacy key maps to typed keys before building context objects.
- `Insertion` and `Dinuclmarkov` now rely on inherited `Rec_Event` typed bridges for the remaining pure-forwarding initialization/proba-bound paths, and the redundant subclass `initialize_event(...)` wrappers have been removed.
- `Dinuclmarkov` now also uses a direct typed `initialize_event(...)` path, enabled by a typed `EventUtils::get_insertion_len_max(...)` helper, so that initialization no longer falls back to rebuilding a legacy map on that path.
- The unused legacy `EventUtils::check_gene_choice(...)` helper has been removed; remaining `EventUtils` cleanup is now limited to helpers that still back compatibility paths.
- Validation state: `31/31` unit tests passing.
-- Deferred follow-up:
    The items below are the next follow-up slice, not blockers for the current validated migration batch.
    - Retire legacy overloads/wrappers once all call sites are typed-only (`EventUtils`, event class legacy overloads, and legacy map bridges).
    - Decide and execute end-state for `Model_Parms` registry storage keying (typed-native storage vs legacy+view), then retire the now compatibility-only `get_events_map()` getter.
    - Complete the `Dinuclmarkov` 2-axis simplification refinement.
    - Complete tandem-D enum establishment and integration checkpoints before final key uniqueness closure.

---

## Non-Goals & Boundaries
- **AIRR / Streaming**: [AIRRRearrangementReader.cpp](file:///Users/tkloczko/Development/igor-immuno-stack/igor/src/igor/Streaming/AIRRRearrangementReader.cpp) and [AIRRRearrangementWriter.cpp](file:///Users/tkloczko/Development/igor-immuno-stack/igor/src/igor/Streaming/AIRRRearrangementWriter.cpp) will remain `Gene_class`-oriented unless the external AIRR schema itself changes. They describe V/D/J alignments on the genomic level, not constructed sequence segments.
- **GeneChoice**: Will remain `Gene_class`-oriented because gene choice is alignment-facing.
