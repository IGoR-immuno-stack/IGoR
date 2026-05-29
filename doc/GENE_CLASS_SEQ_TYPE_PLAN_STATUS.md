# Gene_class and Seq_type Separation: Plan Status

Date: 2026-05-29

## Architecture Decision Record (ADR)

### Keying Choice in Model Registry (`events_map`)
- **Option A**: Keep current staged approach (no key change until much later, keep double identity on events).
- **Option B (Selected)**: Perform early key migration to `tuple<Event_type, Seq_type, Seq_side>`.
  - *Rationale*: As pointed out by the IGoR designer, transitioning directly to this tuple is just as identifying, avoids a complex intermediate double-identity phase, and is a cleaner path towards tandem-D support. The bridge helpers added in Phase 1e isolate the callsites, making this early migration low-risk.

---

## Goal of the Plan
The goal is to separate two different concerns that were historically mixed:

- Gene_class: genomic and alignment identity (V, D, J).
- Seq_type: constructed sequence-segment identity (V_gene_seq, D_gene_seq, J_gene_seq, VD_ins_seq, DJ_ins_seq, VJ_ins_seq).

The strategy is staged to reduce risk:

- **Phase 1**: Semantic cleanup in event/runtime internals without changing model identity keys or file formats. *(Completed)*
- **Phase 2**: Focused Deletion analysis + `events_map` key migration to `tuple<Event_type, Seq_type, Seq_side>` (merging Phase 3 keying work).
- **Phase 3**: Serialization compatibility layer (aligned with Issue #46).

---

## Detailed Summary of What Has Been Implemented Until Now
### Overall execution status
- Phase 1a: completed.
- Phase 1b: completed.
- Phase 1c: completed. *(Note: Phase 1c was successfully completed with the initial 3-axis helper logic. The 2-axis targeting is a planned design refinement to simplify the implementation, not a rollback).*
- Phase 1d: completed.
- Phase 1e: completed across the targeted ancillary/runtime consumers.
- Typed `Seq_type` view added in `Model_Parms` and used by the first production consumer.

### Main implementation outcomes
- Added shared bridge helpers in EventUtils to centralize semantic mapping and event lookups:
  - try_gene_class_to_gene_seq_type(...)
  - try_insertion_gene_class_to_seq_type(...)
  - try_insertion_seq_type_to_gene_class(...)
  - has_insertion_seq_type(...)
  - try_get_event(...) (Implements a fallback lookup pattern that will make the events_map key transition transparent to consumers).

- Replaced repeated manual tuple lookup patterns (count/at) with helper-based lookups in ancillary/runtime code paths.

- Reduced duplicated insertion-presence checks by routing through has_insertion_seq_type(...).

- Added Model_Parms::get_events_map_seq_type() as a typed view over the existing registry.

- Switched Coverage_err_counter to the typed registry view and switched GenModel to initialize cloned error rates from the typed path.

- Added a typed Error_rate::initialize(...) overload so typed consumers can forward through legacy error models without blocking the migration.

- Kept events_map keying unchanged (tuple<Event_type, Gene_class, Seq_side>) during all Phase 1 work, as planned.

### Files updated in this implementation stream
- src/igor/Core/EventUtils.h
- src/igor/Core/EventUtils.cpp
- src/igor/Core/Pgencounter.cpp
- src/igor/Core/Coverageerrcounter.cpp
- src/igor/Core/HypermutationfullNmererrorrate.cpp
- src/igor/Core/Hypermutationglobalerrorrate.cpp
- src/igor/Core/Genechoice.cpp
- src/igor/Core/Insertion.cpp
- src/igor/Core/Deletion.cpp
- tst/igor/Core/test_EventUtils.cpp

### Testing and validation executed
- Unit suite run repeatedly with pixi run test_unit.
- Current result: 31/31 unit tests passing.
- Added direct helper regression coverage:
  - EventUtils HasInsertionSeqType
  - EventUtils TryGetEvent
  - Fallback lookup pattern section inside EventUtils TryGetEvent

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

During Phase 2.3, the event lookup helper APIs will be updated. We will *not* keep compatibility overloads in `EventUtils.h` once migration is completed, as Phase 1e centralized all lookups to these helpers.

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

## Deferred Follow-up

The items below are follow-up work for the next slice. They are not blockers for the validated migration batch.

### Phase 2: Deletion Analysis & Model Key Migration (Active)

#### Step 2.1: Deletion Identity Matrix (Completed)
The inventory and classification artifact exists and has been applied to the staged `Deletion` migration.

#### Step 2.2: Typed bridge expansion (Mostly completed for event runtime paths)
- `Deletion`, `Insertion`, `Dinucl_markov`, and `Gene_choice` now have typed overloads (`Seq_type`-keyed maps); the redundant legacy `initialize_event(...)` wrappers have been removed.
- `Insertion` and `Dinucl_markov` also have typed overloads for `initialize_crude_scenario_proba_bound(...)` with compatibility wrappers.

#### Step 2.3: Atomic Migration of `events_map` Key and Helper Signatures
This step must be executed as a coordinated atomic operation: `Model_Parms` storage format and `EventUtils` lookup signatures are mutually dependent and must transition together.

- Keep the legacy getter as a compatibility path until all consumers are migrated. *(Deferred)*
- Core context stack now uses typed map views (`ModelContext`, `GenModel` context construction, `Rec_Event` legacy adapter bridge, `Pgen_counter`, `Coverage_err_counter`, and context-side error rate bridges). *(Completed)*
- `GenModel` no longer relies on `Model_Parms::get_events_map()` in its runtime initialization/exploration path; the remaining legacy getter is now a compatibility surface rather than a core runtime dependency. *(Completed)*
- `Error_rate` and both hypermutation context bridges no longer reconstruct legacy-keyed event maps when delegating to compatibility `compare_sequences_error_prob(...)` implementations. *(Completed)*
- The deprecated `compare_sequences_error_prob(...)` compatibility interface has been narrowed to drop its dead legacy `events_map` parameter across the remaining implementations and call sites. *(Completed)*
- Hypermutation error-rate classes now dispatch `compute_scenario_error_probability(...)` through internal context-first helpers, leaving `compare_sequences_error_prob(...)` as a compatibility shim only. *(Completed)*
- Redundant typed forwarding wrappers have been removed where `Rec_Event` inherited bridges now cover the same dispatch path (`Insertion` probability-bound initialization and `Dinucl_markov` typed init/proba-bound forwarding). *(Completed)*
- `Dinucl_markov` typed initialization no longer depends on reconstructing a legacy-keyed map; a typed `EventUtils::get_insertion_len_max(...)` helper now supports that path directly. *(Completed)*
- The unused legacy `EventUtils::check_gene_choice(...)` helper has been retired; remaining `EventUtils` legacy helpers still correspond to active compatibility bridges. *(Completed)*
- Remove legacy helper signatures in `EventUtils` once no callers remain. *(Deferred)*
- Decide end-state of model registry storage keying in `Model_Parms` and retire `get_events_map()` at the end of migration. *(Deferred)*

#### Step 2.2: Simplify Dinucleotide samplers to 2-axis targeting (Refinement)
Adjust `Dinuclmarkov` to use a 2-axis coordinate:
- `target Seq_type`
- `anchor side`
- Traversal direction is derived from the anchor side.

### Practical Follow-up Worklist (Near-term)
- Remove legacy compatibility overloads and wrappers once all call paths are confirmed typed-only.
- Finalize end-state for `Model_Parms` event registry storage and remove legacy getter path.
- Finalize tandem-D `Seq_type` enum additions (`D1_gene_seq`, `D2_gene_seq`, `D1D2_ins_seq`) and propagate through mapping/IO checkpoints.
- Complete the `Dinuclmarkov` 2-axis targeting refinement.
- Keep serialization compatibility work in Phase 3 / Issue #46 as planned.

---

## Phase 2 Acceptance Criteria
- Zero compile/link errors.
- 100% of unit tests pass successfully (`pixi run test_unit`).
- All runtime consumers that access `events_map` do so through `EventUtils` helpers (for example `try_get_event` and `has_insertion_seq_type`) and do not perform direct tuple-key lookup calls (`.at(...)`, `.count(...)`). This rule is class-agnostic and applies to future consumers as well.

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

## Non-Goals & Boundaries
- **AIRR / Streaming**: [AIRRRearrangementReader.cpp](file:///Users/tkloczko/Development/igor-immuno-stack/igor/src/igor/Streaming/AIRRRearrangementReader.cpp) and [AIRRRearrangementWriter.cpp](file:///Users/tkloczko/Development/igor-immuno-stack/igor/src/igor/Streaming/AIRRRearrangementWriter.cpp) will remain `Gene_class`-oriented unless the external AIRR schema itself changes. They describe V/D/J alignments on the genomic level, not constructed sequence segments.
- **GeneChoice**: Will remain `Gene_class`-oriented because gene choice is alignment-facing.

