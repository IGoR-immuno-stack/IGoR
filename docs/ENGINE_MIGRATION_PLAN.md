# Engine migration plan: from `GenModel` to drivers and engines

**Created**: Oct 16 2026
**Scope**: the legacy runtime: `GenModel` (the EM loop), the dispatch and leaf of `Rec_Event`,
`Model_Parms` and `Model_marginals` as the runtime model, `SequenceGenerator`, the four counters,
and the error model at the leaf. The events' semantics, the error models' internals and the
generation draw have plans of their own (siblings below); this plan moves the runtime around them.
**Starts**: after tandem-D milestone 1 (T5), on its own branch (§4).
**Target**: [ARCHITECTURE_SYNTHESIS.md](ARCHITECTURE_SYNTHESIS.md) §2 and steps 1–3, 6 and 8;
stage 2 of [doc/LAYER_REFACTORING_PROPOSAL.md](../doc/LAYER_REFACTORING_PROPOSAL.md) ("the DFS
becomes an engine").
**Siblings**: [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md)
(CAP: tandem D), [ERROR_MODEL_REWRITE_PLAN.md](ERROR_MODEL_REWRITE_PLAN.md) (ERR),
[GENERATION_REWRITE_PLAN.md](GENERATION_REWRITE_PLAN.md) (GEN),
[ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) (ITER). The design documents
it builds on are the synthesis (SYN), [SEGMENT_DECOMPOSITION_REVIEW.md](SEGMENT_DECOMPOSITION_REVIEW.md)
(REV), [TK_REFACTORING_MERGE_ANALYSIS.md](TK_REFACTORING_MERGE_ANALYSIS.md) (MERGE) and the layer
proposal (LAYER).
**Line numbers** are those of `feature/tandemD` at `8435f98`, where T1 and T3 are merged.
**Prefixes**: L (migration steps), Q (decisions). CT, CX and CR are the counter plan's; E, ER, H,
HR and X are ERR's; T and B are CAP's.
**Ordering across all the plans**: [IMPLEMENTATION_GRAPH.md](IMPLEMENTATION_GRAPH.md).

---

## 0. Summary

The migration moves IGoR's runtime onto the architecture of the synthesis, in this order:
1. **the driver/engine split** (L1): `GenModel::infer_model` becomes a generic EM driver around a
   `LegacyDfsEngine` that runs today's code;
2. **the leaf to the engine, and a scenario record** (L2, L3): the dispatch, the leaf and its
   policies leave `Rec_Event`; counters and error models read a record of the scenario instead of
   event members;
3. **parameters and accumulators in the new model** (L4): `RecombinationModel`'s tensors and the
   handlers replace `Model_marginals` at runtime;
4. **the engine reads `RecombinationModel`** (L5): `Model_Parms` is left as an importer;
5. **generation** (L7), on a lane of its own;
6. **the legacy types go** (L8): stage 2 of the layer proposal.

**Every step is bitwise**, and the claim is checked by a full-precision A/B oracle (L0), because
the goldens print 6 significant digits. **There is no runtime switch**: `GenModel::infer_model`
becomes a shim that each step hollows further, and the CLI changes only at L8.

**The first step needs no other plan.** L1 is code motion over the legacy types: it needs
milestone 1, nothing else (§6). T3, which it builds on, is merged.

**What the migration unblocks rather than includes** (§7): second engines (DFS-B, `ChainDPEngine`,
the DP-bounded DFS, posterior sampling), B10, covariates, D1/D2 parameter tying.

**Real tandem data.** Milestone 1 is a runtime and computation test on a synthetic V-D1-D2-J model
(CAP T5). A fixed VDDJ layout cannot explain the VDJ reads a real repertoire also holds, so
biological use needs milestone 2 (B10). T5's timings decide where B10 goes (Q14): if the legacy DFS
is prohibitive on V-D1-D2-J, as expected, real tandem data waits for the DP engine (SYN step 9),
which comes after the migration, and B10 is placed with it.

---

## 1. Where things stand

### 1.1 The new architecture, as built

Nothing in production uses it: `app/igor` and `app/igor-demo` call `GenModel` and
`SequenceGenerator` only.

| Component | What it does | What it lacks |
|---|---|---|
| `RecombinationModel` ([RecombinationModel.h:70-73](../src/igor/Model/RecombinationModel.h#L70)) | owns a `Topology` and one tensor per event, indexed by uid, parents first and own axis last | no error model, no registry, no fixed flags; `recombination_model_from_files` goes through the bridge, not JSON |
| `Topology` ([Topology.h:72-76](../src/igor/Model/Topology.h#L72)) | holds `legacy::Rec_Event`s, edges sorted by uid, a segment-order vector; `topologicalOrder()` by priority, ties to the smaller uid | no registry, release builds check nothing (asserts only); `topology_from_json` ignores the document's `error_rate` |
| `LegacyBridge` ([LegacyBridge.tpp](../src/igor/Model/LegacyBridge.tpp)) | `Model_Parms` → `Topology` (event copies, edges, order); `Model_marginals` → tensors | `export_to_legacy` (48-63) carries no order, no error rate and calls no `finalize()`, so the exported model cannot even be copied ([Model_Parms.cpp:106](../src/igor/Model/Legacy/Model_Parms.cpp#L106) dereferences the null error rate) and cannot reach `GenModel`; no tensors → marginals direction; lookups by generated name, which collide for D1/D2 until T2 |
| `EventFactory` | builds the four events from a JSON node, keyed by type string | the handler factories key by `EventType` instead |
| `InferenceEngine<T>` ([InferenceEngine.tpp:98-158](../src/igor/Inference/InferenceEngine.tpp#L98)) | a container of handlers: `run(eStep)` = reset, the caller's callback, normalise; `combineAccumulators` by name | no E-step, read loop, threads, likelihood, error model, fixed events or outputs; the merge is unordered; it normalises fixed events too |
| Inference handlers | Categorical and Markov M-steps over the last axis, `double` and `long double` | the zero cut-offs differ from the legacy `sum != 0` ([Rec_Event.cpp:490](../src/igor/Model/Legacy/Rec_Event.cpp#L490)): ≤ 10ε on the 1-D path ([CategoricalInferenceHandler.tpp:83](../src/igor/Inference/CategoricalInferenceHandler.tpp#L83)), 1e-15 for Markov rows ([MarkovInferenceHandler.tpp:86](../src/igor/Inference/MarkovInferenceHandler.tpp#L86)); no `fixed`, no `empty_copy` |
| `SamplingEngine` ([SamplingEngine.tpp:105-146](../src/igor/Generation/SamplingEngine.tpp#L105)) | one index per node, in topological order | no sequence, no `sampleSequence`, no errors, no output; a root Markov draw from row sums is uniform, not anchored; `run()` passes parent indices where `sample()` expects the from-state; `double` only ([SamplingHandlerFactory.cpp:31-36](../src/igor/Generation/SamplingHandlerFactory.cpp#L31)) |
| Tests | `test_ModelInference` runs `GenModel` inside the callback (a hybrid) and discards the inferred error rate; `test_NormalizationCompareLegacy` is the only legacy-vs-new numeric comparison, of the M-step only | no test compares the new engines to a legacy run |

### 1.2 The legacy driver

REV §5.2 splits `GenModel::infer_model` into generic (G) and DFS-specific (D) steps. The split
stands; its line numbers predate the step 1c move. Today
([GenModel.cpp](../src/igor/Inference/Legacy/GenModel.cpp)):

| Step | Lines | G/D |
|---|---|---|
| Option checks; viterbi forces the ratio to 1 | 114-128 | G + D |
| The queue: event order, marginal layout, parent strides | 133 | D (order), G (layout) |
| The traversal plan, built once for every thread, and its admissibility check | 134-141 | D |
| Log and info files, initial dumps | 146-203 | G |
| Counters' master `initialize` (file streams) | 210-223 | G |
| Fresh accumulators; `fast_iter` (best alignments, iteration 0); the chunk plan (#68) | 237-296 | G |
| Per thread: copies of the model, marginals, error rate and counters | 318-332 | G (counters), D (copies exist because events hold state) |
| DFS buffers, index map, crude bounds, queue, `initialize_event`, `next_event_ptr_arr` | 339-422 | D |
| Error model and counters bound to the thread's events | 425-440 | G made D by raw pointers |
| The junction-length fold on one thread, adopted by the others; `update_event_internal_probas` | 449-496 | D |
| Per read: slot, contexts, **`first_event->iterate`** | 505-579 | G, and D at 579 |
| Per-read normalisation, log line, counters, mean-error rejection, chunk drain | 595-637 | G |
| Ordered merge (marginals, error rate); counters merged in thread order | 657-681 | G |
| M-step: `error_rate->update()`, `normalize`, `copy_fixed_events_marginals` | 688-696 | G |
| Per-iteration and final files | 699-709 | G |

The leaf is `Rec_Event::iterate_wrap_up`
([Rec_Event.cpp:328-442](../src/igor/Model/Legacy/Rec_Event.cpp#L328)): the error factor (399),
the threshold and the running maximum (423-425), the `Scenario` snapshot and every counter
(428-432), then `add_to_marginals` on every non-fixed event (435-439).

The event order is `get_model_queue()`: priority first, then the order of the event list, since
`std::list::sort` is stable ([Model_Parms.cpp:482](../src/igor/Model/Legacy/Model_Parms.cpp#L482),
[Rec_Event.h:794-800](../src/igor/Model/Legacy/Rec_Event.h#L794)).

### 1.3 Event state read from outside the event

- The slot accessors, raw pointers into live members:
  `Gene_choice::alignment_offset_slot()` ([Genechoice.h:67](../src/igor/Model/Legacy/Genechoice.h#L67)),
  `Rec_Event::current_realization_index_slot()` ([Rec_Event.h:464](../src/igor/Model/Legacy/Rec_Event.h#L464)),
  `Deletion::deletion_value_slot()` ([Deletion.h:71](../src/igor/Model/Legacy/Deletion.h#L71)).
  `Gene_choice` publishes the first two only for a gene at an end of the ordering
  ([Genechoice.cpp:260](../src/igor/Model/Legacy/Genechoice.cpp#L260), 758).
- `get_current_realizations_index_vec()` ([Rec_Event.h:460](../src/igor/Model/Legacy/Rec_Event.h#L460)),
  the per-event scratch the best-scenarios counter reads.
- The hypermutation models' `friend` access ([Genechoice.h:60-61](../src/igor/Model/Legacy/Genechoice.h#L60),
  [Deletion.h:63-64](../src/igor/Model/Legacy/Deletion.h#L63)); ERR H3 removes it.

### 1.4 What the gate pins

- Five tracks: `cli`, `align`, `inference`, `no_d_align`, `generate`. Only `inference` and
  `no_d_align` produce counter outputs, all on the one human TRB model.
- **The comparison is textual.** Marginal files are written at 6 significant digits, so "bitwise
  on the five tracks" means equal to 6 digits. `inference_info.out` is not compared.
- No VJ model is run by any track.
- No track pins viterbi (`--MLSO`), fixed events (`--fix_err`, fixed marginals), or the
  `fast_iter` data policy on its own; the mean-error read threshold has no CLI option.
- Marginals are independent of the thread count since #68; counter summaries are not (CNT, CX4).

---

## 2. The counters

The counters have a plan of their own, [COUNTER_REWRITE_PLAN.md](COUNTER_REWRITE_PLAN.md) (CNT), which is lane C of §4. What it covers:
- the inventory, and the patterns by which the counters reach into events;
- why `Errors_counter` goes;
- the scenario record and the segment provenance that replace those patterns;
- the requirements a counter checks at initialisation;
- its defects (CX), the tandem-D generalisation, and its steps (CT, CR);
- the output: buffered row sinks, a Parquet backend, the best scenarios as AIRR rearrangements.

Where it meets the migration:

| Counter step | Needs from the migration | The migration needs it for |
|---|---|---|
| CT0, the characterization suite | — | L1: the counters are pinned before the driver split |
| CT1b, `Errors_counter` dropped | L1's `ReadSummary`, softly | — |
| CT2b, the record and the provenance | L3, which adds the record to `ScenarioContext` | — |
| CT4, the lifecycle in the driver | L1, L4a | — |
| CT3, CT4 | — | L8: `Counter` becomes `Observer` |
| CT6, CT7, the sinks and the AIRR output | CT4, so L4a | — |

The segment provenance (CNT §4.2) is written by `Gene_choice` and `Deletion` at their layers, as
part of L3's record. The hypermutation models read it too (ERR H3).

---

## 3. Migration steps

| Step | Needs | Content | Bitwise |
|---|---|---|---|
| L0 | milestone 1 | **The oracle.** (a) An A/B dump at full precision (marginals, likelihoods, error rate, at `max_digits10`) of the base and the candidate on every track: a script, not new goldens. (b) A VJ track (a shipped TRA or IGK model, inference, evaluate, counters). (c) Pins for fixed events, MLSO/viterbi, `fast_iter`, and 1 vs N threads. (d) Bridge parity on every shipped model: each multi-parent tensor's axes against the legacy strides, and the round trip with order and error rate (`[!shouldfail]` until L5). (e) CT0. (f) ERR H0 brought forward, so that the hypermutation models are pinned before L1 touches the leaf | tests only |
| L1 | L0, T7a | **The EM driver and `LegacyDfsEngine`** (SYN step 6, first half). The driver takes §1.2's G steps: iterations, data policy, chunk plan, threads, read loop, `ReadSummary`, per-read normalisation, read policy, ordered merge, the M-step call, outputs, the counter lifecycle. The engine is concrete, with no abstract base while it is the only engine. It implements `check_admissible` (the TraversalPlan and T7a's guards), `prepare_iteration`, `make_workspace` and `process_read`, by calling today's code. An `AccumulatorSet` interface, first implemented over `Model_marginals`, `Error_rate` and the counters. `GenModel::infer_model` becomes a shim. New: driver unit tests against a fake engine (rejection, fixed events, ordered merge, data policy) | yes |
| L2 | L1; ERR E2 | **The leaf and the dispatch to the engine** (REV §10.3, DFS-A). An engine continuation `next(p, bound_facts)` replaces `iterate_wrap_up` and `next_event_ptr_arr`. The leaf — error factor, threshold, running maximum, counter sink, accumulation — is the engine's, and so is the viterbi policy. The unit harness's `RecordingEvent` gives way to a lambda continuation | yes |
| L3 | L2, T2, Q3–Q5 | **The uid and the scenario record** (SYN steps 2–3). uids on the legacy events, in nickname order. `edges`, `Index_map`, `processed_events` and `offset_map` keyed by uid. `ScenarioContext` gains `{uid → index, payload}`, written at each event's layer; the leaf accumulation reads it, with the flat marginal index as payload, so the arithmetic is unchanged. The provenance record and its reader (CNT §4.2). CT2b lands here. The hypermutation `friend`s stay until ERR H3 | yes |
| L4 | L3; ERR E1, E4; Q7 | **Parameters and accumulators in the new model** (SYN step 8, in MERGE §6.2's order). **L4a**: tensor accumulators in the legacy layout, `long double`, chunk slots, the ordered merge; the counters' slots (CT4). **L4b**: the M-step by the handlers' `maximize`, with the `!= 0` cut-off and the legacy summation order; a `fixed` flag per node; the error node (Single's core, an adapter for the hypermutation models); the bound table rebuilt after the M-step; parity tests at exact equality. **L4c**: the DFS reads its parameters from the tensors, and `update_event_internal_probas` from them too: the flat index it computes over one array becomes an offset into the event's own tensor, and since the layouts are the same, the values and the arithmetic are unchanged. `Model_marginals` keeps file I/O only. The Model module's `InferenceEngine<T>` container dissolves — reset into `AccumulatorSet::empty_copy`, combine into its ordered `+=`, update into the driver's `maximize` loop — and the name becomes the strategy interface `LegacyDfsEngine` implements (Q1) | yes |
| L5 | L4; the registry, the ordering and the linear layout owned by `RecombinationModel`, with `finalize()`'s work there (MERGE §4.2); T7b; an error-model JSON reader (ERR X20); Q4 | **The engine reads `RecombinationModel`.** Its plan — queue, layers, adjacency, TraversalPlan — comes from `Topology` and the registry. The model is const, and the engine copies the events into its workspace. `Model_Parms` is left as the v1/v2 importer, and `export_to_legacy` goes | yes |
| L6 | L5 | **Tiers** (REV §13.3): 0a work (copies, `initialize_event`, the plan) once per run, 0b (the fold, the dinucleotide matrix, the bound tables) once per iteration, shared and const. Each move is benchmarked | yes |
| L7 | L0; ERR E1, E4; GEN §6 | **Generation.** A `SamplingDriver`: RNG streams per `(seed, chunk)`, the same output files. `SamplingEngine` gains: the indexed CDF draw; construction by role through `construct_realization`; Dinucl through `sampleSequence(anchor, length)`; errors through the error node. `SequenceGenerator` retires; FastGenerator's fate is decided (T4's second half). It runs in parallel with L1–L5: `SamplingEngine` already reads `RecombinationModel` | the `generate` track moves once (SYN step 1) |
| L8 | L5, L7; ERR H2–H5 or the Q7 adapter; CT3, CT4 | **The legacy types go** (LAYER stage 2). `app/igor` calls the drivers; `GenModel` is deleted; the `iterate()` bodies move to `Inference/Dfs/` as the engine's kernels, a move; `Rec_Event` keeps the definition (Q15); `Counter` becomes `Observer`; the contexts go where Q16 says; the `Seq_type` enum goes with T6b; the bridge is the importer | yes |

**Against the synthesis' order.** SYN step 6 → L1, L2; steps 2 and 3 → L3; step 8 → L4; step 4
→ L5; step 1 → L7. Steps 4b, 5, 7 and 9 come after the migration (§7). L1 goes first because it
needs nothing else and makes the driver testable on its own; step 1 was already deferred by T4's
option (a).

**ERR Part I within it.** E0 → L0 (lane T); ER1 lands before L1, so that no golden moves inside the
migration; E2 before L2; E1 and E4 before L4; E3 before CT3.

---

## 4. Branches and review flow

**One integration branch**, `feature/engine-migration`:
- cut from `feature/tandemD` once milestone 1 is merged, or from `develop` if `feature/tandemD` has
  reached it by then;
- never rebased: its base is merged into it at each checkpoint below;
- its commits carry `Track: ModelRefactoring` (CAP, *Commit tags*).

**One topic branch, one PR, per step** (`mig/L1-driver`, `mig/CT0-counters`, …):
- merged into the integration branch with a merge commit, as layers 1a–1c were (#73);
- rebased onto the integration tip before review;
- mechanical moves and wiring in separate commits; a golden-moving commit alone in its PR;
- the status column of this plan updated in the same PR.

**Lanes by what each step touches.** Steps whose file sets are disjoint run in parallel; steps that
share files run in sequence.

| Lane | Steps | Files | Runs |
|---|---|---|---|
| T, tests | L0, CT0, ERR E0, H0, H1 | `tst/`, `scripts/tests/` | in parallel with everything |
| D, driver and engine | L1 → L2 → L3 → L4 → L5 → L6 → L8 | `GenModel.*`, `Rec_Event.*`, the contexts, the new `Inference/Dfs/` | strictly in sequence, one open PR at a time |
| E, error model | E1, ER1, E2, E3, E4; H2–H5 | `Errorrate*`, `Singleerrorrate*`, `Hypermutation*`, the new core in `src/igor/Model/` | E2 merged before L2 opens; E3 between two D steps; the rest in parallel |
| C, counters ([COUNTER_REWRITE_PLAN.md](COUNTER_REWRITE_PLAN.md)) | CT1, CT1b, CT2a, CR, CT3, CR-P, CT6, CT7 | `Inference/Legacy/*counter*`, `Scenario.h` | in parallel with D, except CT2b (with L3) and CT4 (after L4a) |
| M, Model module | registry ownership, labelled axes (Q3), handler cut-offs, `fixed`, the container's dissolution | `src/igor/Model/`, `src/igor/Inference/` outside `Legacy/` | each as its own PR, ahead of the D step that uses it |
| G, generation | L7, R4 | `src/igor/Generation/`, the sampling handlers; R4 also edits `Dinuclmarkov.cpp`, a lane-D file | in parallel with D; R4 scheduled between two D steps |

```
develop ─────────────────────────────●──────────────●────────────────●───▶
                                      ▲ after L2      ▲ after L5        ▲ after L8
feature/tandemD ──(milestone 1)──┐    │               │                 │
                                 ▼    │               │                 │
feature/engine-migration ────────●────●───────────────●─────────────────●──▶
   lane T  ───(L0, CT0, E0, H0, H1)──▶ merged as they land
   lane D  ───L1──L2──L3──L4──L5──L6──L8   (one open PR at a time)
   lane E  ───E1/ER1──E2──E3──E4──H2…H5
   lane C  ───CT1, CT1b, CT2a, CR──CT3 (best scenarios, Pgen)──CT2b (with L3)──CT3 (coverage)──CT4 (after L4a)──CR-P
   lane M  ───PRs ahead of L4, L5
   lane G  ───────────L7, R4
```

**Before every PR**: `git merge-tree` against the integration tip and against every other open PR
of lanes D, E and C. Two plans' steps once added functions at the same spot of `Rec_Event.cpp` and
conflicted; the check finds that before review.

**The gate of every PR**: the A/B oracle against the integration tip; the five tracks, the VJ track
and the hypermutation track; the unit suites; a benchmark for L4c and L6.

**Checkpoints into `develop`**: after L2 (the driver and the engine exist), after L5 (the engine
reads `RecombinationModel`), after L8 (the legacy types are gone). Every step being bitwise, each
checkpoint is a releasable state.

---

## 5. Running both architectures during the migration

ERR §3.9 states the tensions and the rules: one class with two lifecycles, two homes for the
parameters, no bitwise equality by default, repairs landing twice; logic in engine-agnostic
components, the bridge as the only copy point, parity tests, no new-architecture code including
`GenModel`. For the migration, add:
- **No runtime flag** (Q6). Every step is bitwise, so the shim replaces the old path in place; there
  is never a production choice between two runtimes.
- **The bridge until L5, then the importer.** After L5 nothing exports to `Model_Parms`.
- **The A/B oracle** is the parity test of every lane-D step.
- **`test_ModelInference` becomes a real engine test at L4**, once the driver runs the E-step on the
  tensors; until then it is a hybrid.
- **Parity tests tighten from tolerance to exact equality** when the handlers take over the M-step
  (L4b).

---

## 6. Roadblocks: what other plans must deliver

"Hard" means the step cannot be done, or cannot stay bitwise, without the item; "soft" means it can,
at the cost of redoing work.

| Blocks | Item | Plan | Status (Oct 16 2026) | Hard/soft | Why |
|---|---|---|---|---|---|
| L0 | Milestone 1, on the synthetic V-D1-D2-J model | CAP T5 | pending | hard | the migration starts after it |
| L1 | T3 merged (`SegmentIndex`, `TraversalPlan`, `FlankCheck`) | CAP T3 | merged Oct 20 2026 (#76) | hard | the TraversalPlan is the engine's `check_admissible`; without it L1 splits code that T3 replaces |
| L1 | T7a and its two guards | CAP T7, ERR §4.1 | planned | hard | the guards are admissibility conditions; without T7a pruning is unsound on tandem (X2) |
| L1 | ER1 placed before L1 | ERR §4.2 | after milestone 1 | soft | keeps every golden move outside the migration |
| L1, hypermutation | H0, the hypermutation track | ERR §4.3 | not started | hard for those models | the split of the leaf would be unverified on them |
| L2 | E2: gating and viterbi out of the error model | ERR §4.2 | after milestone 1 | soft | otherwise the engine's leaf inherits the error model's gating (X4, X5) |
| L2 | The unit harness | ITERATE_TEST_GUIDE | — | soft | `RecordingEvent` is a fake successor in `next_event_ptr_arr` |
| L3 | T2: distinct names; the uid | CAP T2, SYN step 2 | T2 pending; uid not started | hard | the record is keyed by uid; the bridge looks up by name |
| L3 | Q3–Q5 | here | open | hard | the uid order must drive neither the layout nor the queue |
| L4 | E1 (the core) and E4 (the node, the bridge carrying r, registrars, the zero-row guard); X20 | ERR §4.2 | after milestone 1 | hard | the error model has no handler otherwise; `maximize` runs on every node |
| L4 | Handler cut-offs at `!= 0`, the legacy summation order, `fixed`, an ordered merge, `long double` | SYN step 8, MERGE §6.2 | not started | hard | bitwise |
| L5 | The registry, ordering and layout owned by `RecombinationModel`, `finalize()` there | MERGE §4.2, SYN §5 | partial: `Topology::seqTypeOrder()` | hard | the engine's plan needs adjacency |
| L5 | T7b: the hypermutation models by name, the pinning gone | CAP T7, ERR Part II (needs H0, H1) | after milestone 1 | soft for L5, hard for L8 | the registry could carry the pinning a while longer, not past L8 |
| L5 | The importer: v1/v2 → `RecombinationModel` with no `Model_Parms` at runtime; the JSON constructors (done via T0); a `topology_to_json` writer (LAYER §6, hand-over); the error-model JSON reader (X20); label-matched parameters | SYN step 4, LAYER | partial | hard (error model), soft (labels) | the engine loads its model through it |
| L7 | The SamplingEngine connection: uid ↔ event, construction by role, Dinucl by `sampleSequence`, the RNG order | GEN §6 | pending | hard | — |
| L7 | The Markov sampler: the anchor nucleotide; parent indices vs the from-state in `run()` | GEN §6 | pending | hard | wrong sequences otherwise |
| L7 | The error node's sampling registrar; `double` or `long double` sampling | ERR §3.5, SYN §5 | pending | hard | `SamplingEngine` throws on an unregistered type |
| L7 | FastGenerator, T4's second half | CAP T4 | pending | soft | — |
| L8 | H2–H5, or the Q7 adapter kept | ERR §4.3 | after milestone 1 | hard | `Error_rate` cannot go before |
| L8 | CT3, CT4 | CNT §8 | — | hard | the counters cannot become observers before |
| L8 | T6b, the enum's removal | CAP T6 | pending | soft | — |

---

## 7. Blocked until the migration

| Item | Plan | Unblocked by | Note |
|---|---|---|---|
| Phase C: structural rules as model validation and engine admissibility | CAP | L1 (`check_admissible`), completed at L5 | `TraversalPlan::check_admissible` is its first piece |
| Phase D: `ChainDPEngine` | CAP, SYN step 9 | L5 + SYN step 7 | also the fallback if T5 measures the DFS as too slow on tandem D |
| The rest of Phase A (footprint reads, `candidates`, domain arity) | CAP | their engine consumers, after L5 | — |
| **B10, milestone 2**: an optional D2, which real repertoires need | CAP B10 | after T5's timings (Q14): L3 at the earliest if the DFS suffices, `ChainDPEngine` otherwise | the guards are model data, adjacency per scenario goes in the record, the layer scopes in the engine |
| B9 step 3 (`VDJ_genes` split) and D1/D2 parameter tying, as a shared handler | CAP, SYN §4 | L4b | — |
| The error node inside the engine; X4 and X5 by construction | ERR | L4b, L5 | — |
| H5's Logistic handler | ERR | L4b | — |
| An error rate per parent realization ([parent, 2]) | ERR §3.4 | L5 | an edge on the error node |
| The SamplingEngine connection; R4; errors on insertions (ERR Q5); FastGenerator as an approximate engine | GEN, ITER, ERR | L7 | — |
| Tier-1a hoists, `prepare(read)`, `ReadCache` | ITER §2.10, REV §13 | L6 + SYN step 7 | — |
| `Events_map` keyed by capability rather than by (type, seq type, side) | ITER §9 | L5 | — |
| Pure semantics, `evaluate` / `EventOutcome`, DFS-B, stateless events | SYN step 7, REV §10.3 | after L5 | — |
| `DpBoundedDfsEngine`, posterior sampling (FFBS), stochastic EM | SYN step 9, REV §5.3 | after `ChainDPEngine` | — |
| Counter classes (CT5); k-best on a DP | REV §5.3 | a second engine | — |
| Covariates (exogenous nodes, LinearModel, Logistic) | SYN, REV §6 | L4, L5 | — |
| `fmt` / `validate` / `diff`, overlays, the accumulators section | REV §14–15 | L4 (`AccumulatorSet`), and SYN 4b | — |
| SYN 4b, the sectioned document with `max_digits10` | SYN | placed after L5 | a golden move kept out of the migration |
| AA Pgen: the patched leaf as the Probe node; the frame as a DP interface variable | AA plan, CAP H8 | L4, L5; `ChainDPEngine` | Q13 |
| Structural events, extra segments, mixtures of histories | VDJ review | the layout DAG and the engines | research |

---

## 8. Golden moves across the plans, and where each sits

| Move | Plan | Placed |
|---|---|---|
| ER1, Single's pooled estimate (X13) | ERR | before L1 |
| HR, the hypermutation repairs | ERR | Part II, its own track |
| The indexed CDF sampler, with the new RNG order | SYN step 1, GEN §6 | L7 |
| The `max_digits10` writer (SYN 4b) | SYN | after L5 |
| CR, CR-P, CT1b | CNT §8 | lane C, each named |
| R16 | ITER, CAP T3 | undecided |

None sits inside L1–L6, which is what lets every lane-D step stay bitwise.

---

## 9. Decisions

| # | Question | Recommendation |
|---|---|---|
| Q1 | What `InferenceEngine` names | **Follow the synthesis** (§5, row `InferenceEngine`): the name becomes the strategy with the four calls; the container's jobs go to `AccumulatorSet` (row 4), the driver's `maximize` (row 2) and the parameterizations (row 4). The swap happens at L4, when the container dissolves; the aliases `LegacyEngine` and `Engine` ([InferenceEngine.h:134, 137](../src/igor/Inference/InferenceEngine.h#L134)) go with it. Until then `LegacyDfsEngine` is concrete |
| Q2 | *Settled.* LAYER §8.3 left open which of two parameter stores the DFS engine reads after stage 2: `Model_marginals` (one flat `long double` array addressed through `Index_map`, what the DFS reads and fills today) or `RecombinationModel`'s tensors. The synthesis already answers it (§5, row 1: one store) | the tensors, from L4c; `Model_marginals` keeps the legacy file I/O |
| Q3 | The order of a tensor's parent axes | labelled and stored per node (REV §15.4), the importer keeping the legacy strides; not derived from the uid order. Today they coincide on TRB by inspection, untested elsewhere (L0d) |
| Q4 | `LegacyDfsEngine`'s event order | the legacy rule — priority, then the event list's order (§1.2) — recorded by the importer as an engine hint; JSON-native models break ties by nickname (SYN §4) |
| Q5 | The uid order | nickname order, as SYN says; safe once Q3 and Q4 hold |
| Q6 | How production switches | no runtime flag (§5) |
| Q7 | Hypermutation models before ERR Part II ends | an adapter node over the legacy object, its `update()` as the handler's M-step, so that they ride through L1–L5; rather than refusing them and keeping `GenModel` alive for them |
| Q8–Q12 | *Moved to the counter plan*: coverage on tandem (CNT QC3), summing counters (QC4), rejected reads (QC5), `Errors_counter` (now dropped, CNT §2), Pgen's modes (QC6) | — |
| Q13 | When AA Pgen lands | not during L2–L5, which rewrite the code it edits; either on the legacy DFS before L2, or onto the engine after L5 |
| Q14 | Where B10 goes | **Re-evaluated once T5's timings are in.** If the legacy DFS is prohibitive on V-D1-D2-J — the expected outcome — real data gains nothing from B10 before the DP engine (SYN step 9, after the migration), and B10 is placed with it. If the DFS suffices, B10 can follow L3, where adjacency per scenario has a home; before L3 it would put the guards and that adjacency into `Model_Parms` and event members, which L3–L5 then move. Either way, B10's representation (CAP B10, its options (b) vs (c)) is decided with its placement |
| Q15 | LAYER tags the events "replace" (copy, rework, delete the original) | promote after extraction instead: the `iterate()` bodies move to kernels and the remaining class is the definition. `Topology` already shares the class with the legacy runtime, so a copy would give it two event types (ERR §3.10) |
| Q16 | Where the contexts go | REV §9.1 — `QuerySequenceContext`, `ModelContext`, `ScenarioContext` and `AccumulationContext` shared by every engine, `ExplorationContext` the DFS's — over LAYER §5, which moves all five into the DFS engine |

---

## 10. Corrections to other documents

Fixed in the commit that adds this plan:
- **GEN §6**: the `MarkovSamplingHandler` row-index defect was fixed by `72a66fc`; what remains open
  is the anchor nucleotide and the from-state passed by `run()`.
- **CAP**: the T0 blocker row is ticked (#70); T5 gains its synthetic fixture; T7b's counter half
  moves to the counter plan's CT2a; the Phase C, D and E rows point here.
- **ERR**: pointers to this plan (header, §3.9, §7, the T7b row).

Stale statements in the untracked design documents, listed for their next revision:

| Document | Statement | Now |
|---|---|---|
| SYN open decision 7 | Single's M-step: keep the mean of per-read frequencies, or pool | decided: pooled (ERR X13, ER1) |
| SYN §5, REV §6 | Single on a Bernoulli handler | a two-cell Categorical, which is the same distribution (ERR §3.5) |
| SYN open decision 3, REV §7.1 | no JSON library in the tree | `nlohmann` is used by `ModelJson` and the JSON constructors |
| SYN §5, "Build" row | the legacy runtime stays in `Core` | it is in `Inference/Legacy/` and `Model/Legacy/` since step 1c |
| REV §6 | Single's M-step must stay as is to keep the goldens bitwise | ER1 moves them once, deliberately |
| REV §7.2 | `Coverage_err_counter` is a `friend` | it reads public slots since step 1c; the hypermutation models are still `friend`s |
| REV §7.2, §9.1 | `pruning_mismatch_floor` is dead code | it becomes the input of the error-model bounds (ERR, decisions) |
| REV §5.2, §11 | line numbers; the Markov defect open | §1.2 here; fixed by `72a66fc` |
| MODEL | a copy of the Model README: own axes first, `read_topology()` | parents first; `read_topology()` deleted on the integration branch |
| CONTEXT_BASED_ITERATE, ITERATE_METHOD_ANALYSIS | the contexts "complete"; the 19-parameter interface | REV §9 |

---

## 11. Hazards

- **The 6-digit gate.** Without L0a, "bitwise" is equality to 6 digits, and a change in the last
  bits of an accumulator passes unseen.
- **Layout coincidence.** The legacy strides follow the queue; the tensors' parent axes follow the
  uid. They agree on TRB by inspection only (L0d, Q3).
- **CX3** (counter plan) may write past an array on long reads with a palindromic J: repair it early, whatever the
  migration's pace.
- **AA Pgen** (`origin/feature/AA_PGEN`) edits the leaf, the error model and `iterate()`: Q13.
- **The synthetic tandem model** exercises the computation only: it says nothing about the aligner
  on real tandem reads, nor about repertoires that mix VDJ and VDDJ.
