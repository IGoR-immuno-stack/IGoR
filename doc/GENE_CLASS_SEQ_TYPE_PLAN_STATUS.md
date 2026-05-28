# Gene_class and Seq_type Separation: Plan Status

Date: 2026-05-28

## Goal of the Plan
The goal is to separate two different concerns that were historically mixed:

- Gene_class: genomic and alignment identity (V, D, J, or mixed legacy classes used in event keys).
- Seq_type: constructed sequence-segment identity (V_gene_seq, D_gene_seq, J_gene_seq, VD_ins_seq, DJ_ins_seq, VJ_ins_seq).

The strategy is staged to reduce risk:

- Phase 1: semantic cleanup in event/runtime internals without changing model identity keys or file formats.
- Phase 2: enrich event representation to hold both identities where needed.
- Phase 3: redesign model identity and serialization only if Phases 1 and 2 are still insufficient for tandem-D uniqueness.

## Detailed Summary of What Has Been Implemented Until Now
### Overall execution status
- Phase 1a: completed.
- Phase 1b: completed.
- Phase 1c: completed.
- Phase 1d: completed.
- Phase 1e: implemented across the targeted ancillary/runtime consumers; ready for sign-off.

### Main implementation outcomes
- Added shared bridge helpers in EventUtils to centralize semantic mapping and event lookups:
  - try_gene_class_to_gene_seq_type(...)
  - try_insertion_gene_class_to_seq_type(...)
  - try_insertion_seq_type_to_gene_class(...)
  - has_insertion_seq_type(...)
  - try_get_event(...)

- Replaced repeated manual tuple lookup patterns (count/at) with helper-based lookups in ancillary/runtime code paths.

- Reduced duplicated insertion-presence checks by routing through has_insertion_seq_type(...).

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
- Current result: 30/30 unit tests passing.
- Added direct helper regression coverage:
  - EventUtils HasInsertionSeqType
  - EventUtils TryGetEvent
  - Fallback lookup pattern section inside EventUtils TryGetEvent

### Notes on non-blocking warnings
- Existing warnings in Deletion.cpp remain known and unchanged in nature.
- No new test failures were introduced by the Phase 1 cleanup slices.

## What Remains To Do
### Immediate closure item
- Mark Phase 1e complete after explicit maintainer sign-off on the latest helper-based cleanup scope.

### Phase 2 (next)
- Add richer per-event identity so events can carry both genomic context and explicit target segment context where required.
- Keep events_map keys unchanged during initial Phase 2 unless demonstrated insufficient.
- Define and implement minimal data-model changes in Rec_Event and relevant subclasses.
- Add dedicated tests covering tandem-D representability under the enriched event identity.

### Phase 3 (conditional)
Execute only if Phase 2 still cannot represent tandem-D identities unambiguously:

- Redesign model event identity and keys.
- Update model serialization and legacy loader compatibility mapping.
- Prune mixed Gene_class values only when compatibility and behavior are proven safe.

## Annexe: Initial Plan and Current Status

### Initial staged structure
- Phase 1: Semantic Separation in Event and Sampler Internals
- Phase 2: Richer Per-Event Identity
- Phase 3: Model Identity and Serialization Redesign (conditional)

### Phase 1: Semantic Separation in Event and Sampler Internals
Goal: disentangle genomic identity and sequence-segment representation inside event execution paths while keeping model maps, serialization, and public API keys unchanged.

#### Step 1.1: Codify semantic contract and bridge helpers
Initial plan content:
- Document Gene_class as genomic/alignment identity and Seq_type as segment-role identity.
- Keep all current Gene_class enum values for now.
- Introduce bridge helpers for mapping semantics.

Status: completed.

#### Step 1.2 (Phase 1a): Refactor fast generation
Initial plan content:
- Use fast generation path as low-risk proving ground.
- Replace direct surrogate usage with bridge-based logic.

Status: completed.

#### Step 1.3 (Phase 1b): Refactor insertion internals
Initial plan content:
- Keep constructor/serialization signatures unchanged.
- Refactor iterate/draw logic to resolve insertion segment identity explicitly.

Status: completed.

#### Step 1.4 (Phase 1c): Refactor dinucleotide internals
Initial plan content:
- Replace overloaded Gene_class traversal meaning with explicit target segment and traversal/anchor side behavior.

Status: completed.

#### Step 1.5 (Phase 1d): Refactor deletion internals
Initial plan content:
- Keep genomic deletion identity in Gene_class.
- Move segment updates to explicit segment-oriented logic.

Status: completed.

#### Step 1.6 (Phase 1e): Refactor ancillary consumers
Initial plan content:
- Update event-map and segment-querying runtime consumers.
- Migrate constructed-segment logic while leaving pure genomic checks intact.
- Target files listed by initial plan: Pgencounter, Coverageerrcounter, HypermutationfullNmererrorrate, Errorrate, EventUtils.

Status: implemented across target runtime consumers and adjacent duplicated callsites; ready for sign-off as complete.

### Phase 2: Richer Per-Event Identity
Goal: support tandem-D readiness by allowing events to hold both genomic and segment identities where relevant.

#### Step 2.1: Evolve Rec_Event and subclasses
Initial plan content:
- Add fields for target sequence segments where both identities are needed.
- Preserve compatibility with existing events_map keys.

Status: not started.

### Phase 3: Model Identity and Serialization Redesign
Goal: redesign event registry keys and model file formats only if Phase 2 is insufficient.

#### Step 3.1: Redesign events_map and file I/O
Initial plan content:
- Redesign keying approach to represent both genomic and segment dimensions robustly.
- Update model readers for legacy compatibility.
- Prune hybrid Gene_class values only after successful migration.

Status: not started (conditional phase).

### Initial verification plan and status
1. Bridge verification tests.
- Status: completed (EventUtils helper tests added and passing).

2. Fast generator verification after Phase 1a.
- Status: completed during Phase 1a progression.

3. Full test suite after each step (1b through 1e).
- Status: completed for implemented slices; current unit suite is green (30/30).

4. Long inference/integration checks after key insertion/dinucleotide changes.
- Status: partially completed in prior progression; rerun as needed at Phase 2 kickoff.

5. Legacy loader regression for Phase 3 compatibility mapping.
- Status: not started (reserved for conditional Phase 3).

6. Tandem-D readiness verification.
- Status: not started (planned for Phase 2/3 evaluation).
