# Model/Legacy

The model and the event classes that inference and generation run on today. "Legacy" marks code
from before the layer refactoring that has not reached its final place yet
([LAYER_REFACTORING_PROPOSAL.md](../../../../doc/LAYER_REFACTORING_PROPOSAL.md) §4). It is not dead
code: every inference and every generation goes through it.

Two documents describe how it works:
- [REC_EVENT_CONTRACT.md](../../../../docs/reference/REC_EVENT_CONTRACT.md): what an event
  implements and exposes, hook by hook, and a checklist for a new event type;
- [ITERATE_INTERNALS.md](../../../../docs/reference/ITERATE_INTERNALS.md): the machinery that calls
  the events during an inference: the segment registry, the layered scenario state, the engine's
  plan, overlap safety and the pruning bounds.

Tests are written as [ITERATE_TEST_GUIDE.md](../../../../docs/reference/ITERATE_TEST_GUIDE.md)
says, in `tst/igor/Model/Legacy/`.

## What is here

| Group | Files | Described in |
|---|---|---|
| The events | `Rec_Event.h / .cpp`, `Genechoice`, `Deletion`, `Insertion`, `Dinuclmarkov`, `EventTypedefs.h` | REC_EVENT_CONTRACT |
| The five contexts of a walk, and the leaf's view | `QuerySequenceContext.h`, `ModelContext.h`, `ScenarioContext.h`, `ExplorationContext.h`, `AccumulationContext.h`, `Scenario.h` | ITERATE_INTERNALS §1, §3 |
| The junction fold's state | `FoldFrontier.h` | ITERATE_INTERNALS §6.3 |
| Error models | `Errorrate`, `Singleerrorrate`, `Hypermutationglobalerrorrate`, `HypermutationfullNmererrorrate` | ITERATE_INTERNALS §7 |
| The counter interface | `Counter.h / .cpp`; the counters themselves are in `Inference/Legacy/` | ITERATE_INTERNALS §7 |
| The model and its parameters | `Model_Parms` (events, edges, segment registry, `finalize()`), `Model_marginals` (the flat marginal array and its layout) | ITERATE_INTERNALS §2.2 |
| Model files | `ModelJson` (the JSON document: structure and parameters), `JsonDetail.h` (what the JSON event constructors share), `LegacyFormat` (the v1 format's gene classes) | *T2* in [the capability plan](../../../../docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md) |
| Generation | `GenerationState.h`, which `construct_realization()` builds into | REC_EVENT_CONTRACT §3.6 |
| Helpers | `EventUtils` | — |
| Forwarding headers | `SafetyMatrix.h`, `JunctionGeometry.h`, `SpanProfile.h`, `UnfilledSegmentLengths.h`, `BoundTightness.h` | they keep the old include paths of headers promoted to `Model/` ([README](../README.md)) |
| Migration code | `gene_to_seqtype_migr.{h,cpp}`, deleted at CAP T6b; `vdj_ordering_migr.{h,cpp}`, deleted once CAP T7b and CNT CT2a have landed | a `_migr` file holds code that only bridges a migration, in namespace `migration`, and is deleted whole at the step it names |

> **Subject to change: LAYER stage 2, MIG L8.** This directory goes. `Rec_Event` keeps the
> definition of an event in `Model/` (MIG Q15), the `iterate()` bodies move to `Inference/Dfs/` as
> `LegacyDfsEngine`'s kernels, `Counter` becomes `Observer`, and `Model_Parms` stays only as the
> importer of v1 files (MIG L5). The steps in between are listed in
> [REC_EVENT_CONTRACT.md](../../../../docs/reference/REC_EVENT_CONTRACT.md) §6.
