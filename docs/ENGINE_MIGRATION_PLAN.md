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
**Line numbers** are those of `feature/tandemD-T7-pinning` (T1 merged in, T3 not). "On T3" marks
`feature/tandemD-T3-plan`.
**Prefixes**: L (migration steps), CT (counter steps), CX (counter defects), CR (counter repairs),
Q (decisions). E, ER, H, HR and X are ERR's; T and B are CAP's.
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
milestone 1 and T3 merged, nothing else (§6).

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
| `LegacyBridge` ([LegacyBridge.tpp](../src/igor/Model/LegacyBridge.tpp)) | `Model_Parms` → `Topology` (event copies, edges, order); `Model_marginals` → tensors | `export_to_legacy` (48-63) carries no order, no error rate and calls no `finalize()`, so the exported model cannot even be copied ([Model_Parms.cpp:105](../src/igor/Model/Legacy/Model_Parms.cpp#L105) dereferences the null error rate) and cannot reach `GenModel`; no tensors → marginals direction; lookups by generated name, which collide for D1/D2 until T2 |
| `EventFactory` | builds the four events from a JSON node, keyed by type string | the handler factories key by `EventType` instead |
| `InferenceEngine<T>` ([InferenceEngine.tpp:98-158](../src/igor/Inference/InferenceEngine.tpp#L98)) | a container of handlers: `run(eStep)` = reset, the caller's callback, normalise; `combineAccumulators` by name | no E-step, read loop, threads, likelihood, error model, fixed events or outputs; the merge is unordered; it normalises fixed events too |
| Inference handlers | Categorical and Markov M-steps over the last axis, `double` and `long double` | the zero cut-offs differ from the legacy `sum != 0` ([Rec_Event.cpp:497](../src/igor/Model/Legacy/Rec_Event.cpp#L497)): ≤ 10ε on the 1-D path ([CategoricalInferenceHandler.tpp:83](../src/igor/Inference/CategoricalInferenceHandler.tpp#L83)), 1e-15 for Markov rows ([MarkovInferenceHandler.tpp:86](../src/igor/Inference/MarkovInferenceHandler.tpp#L86)); no `fixed`, no `empty_copy` |
| `SamplingEngine` ([SamplingEngine.tpp:105-146](../src/igor/Generation/SamplingEngine.tpp#L105)) | one index per node, in topological order | no sequence, no `sampleSequence`, no errors, no output; a root Markov draw from row sums is uniform, not anchored; `run()` passes parent indices where `sample()` expects the from-state; `double` only ([SamplingHandlerFactory.cpp:31-36](../src/igor/Generation/SamplingHandlerFactory.cpp#L31)) |
| Tests | `test_ModelInference` runs `GenModel` inside the callback (a hybrid) and discards the inferred error rate; `test_NormalizationCompareLegacy` is the only legacy-vs-new numeric comparison, of the M-step only | no test compares the new engines to a legacy run |

### 1.2 The legacy driver

REV §5.2 splits `GenModel::infer_model` into generic (G) and DFS-specific (D) steps. The split
stands; its line numbers predate the step 1c move. Today
([GenModel.cpp](../src/igor/Inference/Legacy/GenModel.cpp)):

| Step | Lines | G/D |
|---|---|---|
| Option checks; viterbi forces the ratio to 1 | 103-117 | G + D |
| The queue: event order, marginal layout, parent strides | 122 | D (order), G (layout) |
| Log and info files, initial dumps | 127-184 | G |
| Counters' master `initialize` (file streams) | 191-204 | G |
| Fresh accumulators; `fast_iter` (best alignments, iteration 0); the chunk plan (#68) | 218-277 | G |
| Per thread: copies of the model, marginals, error rate and counters | 299-313 | G (counters), D (copies exist because events hold state) |
| DFS buffers, index map, crude bounds, queue, `initialize_event`, `next_event_ptr_arr` | 320-398 | D |
| Error model and counters bound to the thread's events | 401-416 | G made D by raw pointers |
| The junction-length fold on one thread, adopted by the others; `update_event_internal_probas` | 425-472 | D |
| Per read: slot, contexts, **`first_event->iterate`** | 481-555 | G, and D at 555 |
| Per-read normalisation, log line, counters, mean-error rejection, chunk drain | 571-613 | G |
| Ordered merge (marginals, error rate); counters merged in thread order | 633-657 | G |
| M-step: `error_rate->update()`, `normalize`, `copy_fixed_events_marginals` | 664-672 | G |
| Per-iteration and final files | 675-685 | G |

The leaf is `Rec_Event::iterate_wrap_up`
([Rec_Event.cpp:335-449](../src/igor/Model/Legacy/Rec_Event.cpp#L335)): the error factor (406),
the threshold and the running maximum (430-432), the `Scenario` snapshot and every counter
(435-439), then `add_to_marginals` on every non-fixed event (442-446).

The event order is `get_model_queue()`: priority first, then the order of the event list, since
`std::list::sort` is stable ([Model_Parms.cpp:481](../src/igor/Model/Legacy/Model_Parms.cpp#L481),
[Rec_Event.h:778-784](../src/igor/Model/Legacy/Rec_Event.h#L778)).

### 1.3 Event state read from outside the event

- The slot accessors, raw pointers into live members:
  `Gene_choice::alignment_offset_slot()` ([Genechoice.h:66](../src/igor/Model/Legacy/Genechoice.h#L66)),
  `Rec_Event::current_realization_index_slot()` ([Rec_Event.h:455](../src/igor/Model/Legacy/Rec_Event.h#L455)),
  `Deletion::deletion_value_slot()` ([Deletion.h:70](../src/igor/Model/Legacy/Deletion.h#L70)).
  `Gene_choice` publishes the first two only for a gene at an end of the ordering
  ([Genechoice.cpp:272](../src/igor/Model/Legacy/Genechoice.cpp#L272), 769).
- `get_current_realizations_index_vec()` ([Rec_Event.h:451](../src/igor/Model/Legacy/Rec_Event.h#L451)),
  the per-event scratch the best-scenarios counter reads.
- The hypermutation models' `friend` access ([Genechoice.h:59-60](../src/igor/Model/Legacy/Genechoice.h#L59),
  [Deletion.h:62-63](../src/igor/Model/Legacy/Deletion.h#L62)); ERR H3 removes it.

### 1.4 What the gate pins

- Five tracks: `cli`, `align`, `inference`, `no_d_align`, `generate`. Only `inference` and
  `no_d_align` produce counter outputs, all on the one human TRB model.
- **The comparison is textual.** Marginal files are written at 6 significant digits, so "bitwise
  on the five tracks" means equal to 6 digits. `inference_info.out` is not compared.
- No VJ model is run by any track.
- No track pins viterbi (`--MLSO`), fixed events (`--fix_err`, fixed marginals), or the
  `fast_iter` data policy on its own; the mean-error read threshold has no CLI option.
- Marginals are independent of the thread count since #68; counter summaries are not (§2.4, CX4).

---

## 2. The counters

### 2.1 Inventory

None of the four has a unit test.

| Counter | Reached by | Reads | V/D/J hard-coding (code tokens) | On a tandem model | Output | Tests |
|---|---|---|---|---|---|---|
| `Best_scenarios_counter` | `--scenarios N` ([legacy_main.cpp:913](../app/igor/legacy_main.cpp#L913)) | P(s, read); every event's `get_current_realizations_index_vec()` in queue order ([Bestscenarioscounter.cpp:104](../src/igor/Inference/Legacy/Bestscenarioscounter.cpp#L104)); the V, D, J mismatch lists (193) | 6 `Seq_type` | D1/D2 mismatches silently left out; duplicate column names until T2 | `best_scenarios_counts.csv`, once per read, shared stream | goldens: `inference`, `no_d_align` (sorted, rank dropped); `tie_digits = 12`, `keep_ties` only |
| `Errors_counter` | `-run_demo` and `igor-demo` only; no CLI flag | sizes of the V, D, J sequences and mismatch lists (Errorscounter.cpp:117-140) | 24 `Seq_type` | silently undercounts | `scenarios_background_and_errors.csv`, `sequence_mutation_frequency.csv`, per read, every iteration | none (its sort rules in `test_inference.sh` have no golden) |
| `Pgen_counter` | `--Pgen` (estimator mode only) | `build_full_sequence` (Scenario.h:145-174); P(s) and P(s, read) | 12 `Seq_type`, all without effect on the live path; the snapshot's VDJ/VJ branch | σ assembled as V+J: wrong, silently | `Pgen_counts.csv`, per read | goldens in estimator mode; the per-σ and sequence modes untested |
| `Coverage_err_counter` | `--coverage {V_gene, …, VDJ_genes}` | events by enum, raw slots, the gene choice's realizations map (§2.2) | 24 `Seq_type`, 3 gene class | V and J work while they keep their legacy names; D terminates the process | `{V,D,J}_genes_cov_and_err.csv`, summary per iteration | goldens for V and J; D, `Npoint > 1` and per-read mode untested |

Totals: 89 `Seq_type` and 3 gene-class tokens across the four counters, `Scenario.h` and
`EventUtils.cpp`, as ERR §3.7 counts. The T7a guards that would refuse a tandem model are not in the
code yet.

### 2.2 Event knowledge and direct member access: the patterns, and how to break them

The counters, coverage above all, do not just read the scenario: they know which events exist,
which member holds what, and how the events lay a segment out on the read. Each pattern, what it
is used for, and what replaces it:

| # | Pattern | Where | What it actually needs | Replacement |
|---|---|---|---|---|
| A1 | "The V, D or J gene choice" and "the V 3′, D 5′/3′, J 5′ deletion", found by (event type, `Seq_type` enum, side) | [Coverageerrcounter.cpp:163-265](../src/igor/Inference/Legacy/Coverageerrcounter.cpp#L163); Pgencounter.cpp:93-99 | the event that creates segment s, and the one that modifies each of its ends | `SegmentIndex::content_creators(s)` and `end_modifiers(s, side)` (on T3), resolved once at `initialize` from the segment names the counter is configured with |
| A2 | Raw pointers into live event members, dereferenced at every leaf: `alignment_offset_slot` (→ `&alignment.offset` of the current read's alignment), `current_realization_index_slot`, `deletion_value_slot`. Correct only because each thread owns its events and the DFS has just written them; D publishes neither of the first two | Coverageerrcounter.cpp:164-267; Genechoice.cpp:272-277; Deletion.cpp:325 | which gene, where its template sits on the read, how many nucleotides each end lost or gained | the scenario record (L3): `{uid → realization index, payload}`. A deletion's value is the domain's value at the recorded index, read from the const model, so it needs no payload. A gene choice's payload is the read offset of its template, published for every placement, D's exhaustive scan included |
| A3 | Per-event scratch read through a public accessor, in queue order | Bestscenarioscounter.cpp:104 | every event's realization, Dinucl chains included | the record, in the model's order |
| A4 | A copy of the gene choice's realization map, to size the arrays | Coverageerrcounter.cpp:168-174 | the template lengths | the event's domain in the const model |
| A5 | Geometry by assumption: V is left of the read and J right (296-303, 339-351); template position = read position − alignment offset; P nucleotides told apart by where their mismatches fall (V: `pos − off ≥ len`, J: `pos < off`); end bound `len − max(0, del)` | Coverageerrcounter.cpp:296-351 | for each read position of the segment, the template position it shows, and whether directly or mirrored | the **segment provenance** view, below |
| A6 | What "genomic" means, implied: the V, D and J sequences with their P nucleotides, and their mismatch lists | Errorscounter.cpp:117-140; Bestscenarioscounter.cpp:193 | the segments whose content is compared with the read | the mismatch role (ERR §3.3) and Single's resolved `apply_to` (ERR §3.2): one set, shared with the likelihood, so the counter and the error model agree on what an error is |
| A7 | The snapshot's six fixed slots, filled from ids 0 to 5, and the VDJ/VJ branch of `build_full_sequence` | [Scenario.h:57-90](../src/igor/Model/Legacy/Scenario.h#L57), 145-174 | the constructed segments, in layout order | a snapshot sized from the registry and indexed by `SeqTypeId`; the sequence assembled over `registry.ordering()` |
| A8 | Timing implied by the DFS: a counter runs at the leaf because every layer is written by then; `count_sequence` runs before the error model resets its per-read state | Rec_Event.cpp:435-439; GenModel.cpp:586-605 | a completed scenario; a per-read summary | the engine's sink contract (§2.4) and a `ReadSummary` |
| A9 | The read likelihood taken from the error model | GenModel.cpp:588 | L | `ReadSummary` |

**The segment provenance view.** One engine-agnostic helper in `src/igor/Model/`, under the rules
of ERR §3.9: it depends on neither architecture. From the record, `SegmentIndex` and the const
model, it gives, for one segment of one scenario:
- the template: library and realization index;
- the template interval kept, [t5, t3], after both ends' modifiers;
- the palindromic extension at each end, p5 and p3;
- the read interval the segment occupies;
- for each read position of the segment: (template position, direct or mirrored), or "not
  templated".

One helper then serves every consumer that today re-derives the geometry:
- coverage counts template positions from it;
- errors, best scenarios and the error model select segments by mismatch role and count within them;
- the hypermutation models later read their context through it and the flank segments (ERR H3).

Today's counts are the provenance with the mirrored positions dropped, so moving onto it is
bitwise. The slots and the `friend`s go once their last reader uses it.

### 2.3 Defects and oddities

Categories follow ITER §6.1: *confirmed* defects are pinned at their intended value under
`[!shouldfail]` and fixed once, in a CR commit; *odd* behaviour is pinned as observed.

| # | Finding | Category | Owner |
|---|---|---|---|
| CX1 | D coverage throws a `const char*` ([Coverageerrcounter.cpp:321](../src/igor/Inference/Legacy/Coverageerrcounter.cpp#L321)); `GenModel`'s `catch (exception&)` does not catch it, it leaves the OpenMP region, and the process calls `std::terminate` | confirmed | CR, then CT3 implements D |
| CX2 | `try_get_event`'s result is ignored (163, 198, 242): D coverage on a VJ or tandem model dereferences null | confirmed | CR |
| CX3 | J coverage with a palindromic J (`j5del < 0`): the end bound `max(0, j5del) + (J3′ − J5′ + 1)` (347-351) includes the \|d\| P nucleotides, so it over-counts \|d\| positions at the 3′ end and, when the read covers the whole J, writes past the array. The TRB goldens do not reach it (their reads cover about 13 nt of J) | suspected; confirm in CT0 | CR, early, independently of the migration: it can corrupt memory on long reads |
| CX4 | The coverage summary is merged in thread order (GenModel.cpp:640-657, its FIXME): its last bits depend on the schedule; the 6-digit output hides it | confirmed | CT4 |
| CX5 | `initialize` re-allocates without freeing; a destructor over raw arrays with the implicit copy left in (rule of three); `throw except;` slices; messages copied from the hypermutation model | confirmed | CR |
| CX6 | Pgen's sequence mode writes a header and no rows; its per-σ mode does not print σ (Pgencounter.cpp:181-211) | confirmed | CR, Q12 |
| CX7 | Pgen's V/D/J flags are computed and read only by the dead path | dead | CT1 |
| CX8 | Unguarded divisions: by n_gen (Errorscounter.cpp:176), by L (267-269; best scenarios' `count_sequence`) | confirmed | CR |
| CX9 | `Errors_counter` has no CLI flag and no golden; its top-k breaks ties by enumeration order | gap | CT0, Q11 |
| CX10 | Per-read outputs carry no iteration column: with `last_iter_only = false` (errors, in the demo), rows of every iteration are indistinguishable | odd | CT4 |
| CX11 | Counters run before the mean-error rejection (GenModel.cpp:586-605), so rejected reads are counted | odd | Q10 |
| CX12 | Counters see the scenarios above the *running* maximum, so the set depends on the order the DFS explores | inherent to the DFS | documented; a DP engine differs (REV §5.3, check 1 at threshold 0) |
| CX13 | Viterbi: counters accumulate every successive running maximum, so coverage weights exceed 1 | confirmed | CR, with ERR E2 |
| CX14 | Best scenarios: `initialize` is not idempotent (columns duplicate); P is truncated from `long double` to `double` | confirmed | CR |
| CX15 | The deprecated interface (`initialize_counter`, the 7-argument `count_scenario`) and the four legacy bodies have no caller | dead | CT1 |
| CX16 | A `Scenario` snapshot (three vectors) is allocated at every accepted leaf, counters or not | perf | CT2b |
| CX17 | Thread copies share one stream; rows come in thread-finish order, which the tests sort away | odd | CT4 |
| CX18 | Coverage discards the palindromic nucleotides instead of counting them at their mirrored template position | improvement | CR-P |

### 2.4 What each counter needs from an engine

**Statistic class**, after REV §5.3:

| Counter | Class | An engine other than the DFS would serve it by |
|---|---|---|
| Coverage | posterior expectation over one segment's state (gene, placement, end trims) | forward–backward marginals |
| Errors (a) Σ P·n_gen, (b) Σ P·n_mis | additive | the expectation semiring |
| Errors (c) Σ P·n_mis/n_gen; its top-k | non-linear; k-best | state augmentation; k-best |
| Best scenarios | k-best, with its canonical tie order (40-bit mantissa, then lexicographic on realizations, then mismatches; Bestscenarioscounter.cpp:68-98) | k-best (max, ×) with backtracking |
| Pgen | an arbitrary function of the whole sequence | enumeration or posterior sampling |

**The contract a counter needs**, whatever the engine:
- `initialize` once per worker, against a const model;
- a per-scenario sink receiving P(s), P(s, read), the record, the provenance and the mismatch lists;
- the threshold semantics: it sees exactly the scenarios that make up L;
- a per-read hook with the read index, the iteration, L and whether the read was accepted;
- per-chunk accumulation with an ordered merge, for summaries;
- a summary hook per iteration, honouring `last_iter_only`;
- an output sink that writes rows in read order.

### 2.5 Generalisation for tandem D

Each counter moves from V, D and J to sets resolved from the ordering. VDJ outputs stay unchanged
until CR-P.
- **Best scenarios**: the mismatches of every segment with a mismatch writer, in ordering order (in
  VDJ: V, D, J, as today).
- **Errors**: n_gen and n_mis over Single's `apply_to` set.
- **Pgen**: σ assembled over `registry.ordering()`.
- **Coverage**:
  - keyed by segment, not gene class: `--coverage` accepts segment names, and `D_gene` keeps
    meaning every D segment (Q8);
  - D implemented;
  - geometry from the provenance, for any segment, whichever end anchors it;
  - **palindromic nucleotides counted at their mirrored template position (CX18, CR-P).** A
    template position seen directly and again through its palindrome is counted twice; a mismatch on
    a P nucleotide is an error at the mirrored position. Coverage then reads "expected observations
    per template position", not "probability of being observed", and the output header says so.

### 2.6 Steps

| Step | Needs | Content | Bitwise |
|---|---|---|---|
| CT0 | milestone 1 | Characterization suite at the counter API: `Scenario` and contexts built by hand, no DFS needed. One tiny end-to-end model per counter. VJ and tandem fixtures. Best scenarios: tie settings, viterbi. Pgen: all three modes. Coverage: D (`[!shouldfail]`), `Npoint > 1`, per-read mode, palindromic V and J (CX3). An Errors golden. Thread-count determinism at full precision (`[!shouldfail]` until CT4). Coverage report and mutations (ITER's ITERATE_TEST_GUIDE §7) | tests only |
| CT1 | CT0 | Dead code out (CX7, CX15) | yes |
| CT2a | CT0 | T7b's counter half, on the legacy DFS: the snapshot sized from the registry and indexed by `SeqTypeId` (A7); the events a counter reads resolved once, at `initialize`, through `SegmentIndex` rather than by enum (A1); best scenarios' mismatches in ordering order | yes |
| CT2b | CT2a, L3 | The record and the provenance replace the slots, the scratch reads, the copied realization map and the assumed geometry (A2–A5). The slots go | yes |
| CT3 | CT2a and ERR E3 for best scenarios, errors and Pgen; CT2b for coverage | §2.5, which lifts T7a's guard counter by counter | yes on VDJ; new tandem cases |
| CT4 | L1, L4a | The lifecycle in the driver: the per-read hook takes a `ReadSummary`; an output sink ordered by read index; summing counters in `AccumulatorSet` chunk slots, sparse (Q9); `last_iter_only` applied by the driver | yes, except the coverage summary's last bits (CX4), measured |
| CR | CT0 | One repair per commit: CX1, CX2, CX3, CX5, CX6, CX8, CX13 (with ERR E2), CX14 | the outputs each moves, named |
| CR-P | CT3 | Palindromic positions counted (CX18) | the coverage goldens move once, named |
| CT5 | a second engine | Each counter declares its statistic class and what it reads, for the engine's admissibility check; `Counter` becomes `Observer` (LAYER §5) | §7 |

---

## 3. Migration steps

| Step | Needs | Content | Bitwise |
|---|---|---|---|
| L0 | milestone 1 | **The oracle.** (a) An A/B dump at full precision (marginals, likelihoods, error rate, at `max_digits10`) of the base and the candidate on every track: a script, not new goldens. (b) A VJ track (a shipped TRA or IGK model, inference, evaluate, counters). (c) Pins for fixed events, MLSO/viterbi, `fast_iter`, and 1 vs N threads. (d) Bridge parity on every shipped model: each multi-parent tensor's axes against the legacy strides, and the round trip with order and error rate (`[!shouldfail]` until L5). (e) CT0. (f) ERR H0 brought forward, so that the hypermutation models are pinned before L1 touches the leaf | tests only |
| L1 | L0, T3 merged, T7a | **The EM driver and `LegacyDfsEngine`** (SYN step 6, first half). The driver takes §1.2's G steps: iterations, data policy, chunk plan, threads, read loop, `ReadSummary`, per-read normalisation, read policy, ordered merge, the M-step call, outputs, the counter lifecycle. The engine is concrete, with no abstract base while it is the only engine. It implements `check_admissible` (the TraversalPlan and T7a's guards), `prepare_iteration`, `make_workspace` and `process_read`, by calling today's code. An `AccumulatorSet` interface, first implemented over `Model_marginals`, `Error_rate` and the counters. `GenModel::infer_model` becomes a shim. New: driver unit tests against a fake engine (rejection, fixed events, ordered merge, data policy) | yes |
| L2 | L1; ERR E2 | **The leaf and the dispatch to the engine** (REV §10.3, DFS-A). An engine continuation `next(p, bound_facts)` replaces `iterate_wrap_up` and `next_event_ptr_arr`. The leaf — error factor, threshold, running maximum, counter sink, accumulation — is the engine's, and so is the viterbi policy. The unit harness's `RecordingEvent` gives way to a lambda continuation | yes |
| L3 | L2, T2, Q3–Q5 | **The uid and the scenario record** (SYN steps 2–3). uids on the legacy events, in nickname order. `edges`, `Index_map`, `processed_events` and `offset_map` keyed by uid. `ScenarioContext` gains `{uid → index, payload}`, written at each event's layer; the leaf accumulation reads it, with the flat marginal index as payload, so the arithmetic is unchanged. The provenance helper (§2.2). CT2b lands here. The hypermutation `friend`s stay until ERR H3 | yes |
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
| C, counters | CT1, CT2a, CR, CT3, CR-P | `Inference/Legacy/*counter*`, `Scenario.h` | in parallel with D, except CT2b (with L3) and CT4 (after L4a) |
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
   lane C  ───CT1, CT2a, CR──CT3 (best, errors, Pgen)──CT2b (with L3)──CT3 (coverage)──CT4 (after L4a)──CR-P
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
| L1 | T3 merged (`SegmentIndex`, `TraversalPlan`, `FlankCheck`) | CAP T3 | delivered on `feature/tandemD-T3-plan`, not merged | hard | the TraversalPlan is the engine's `check_admissible`; without it L1 splits code that T3 replaces |
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
| L8 | CT3, CT4 | here | — | hard | the counters cannot become observers before |
| L8 | T6b, the enum's removal | CAP T6 | pending | soft | — |

---

## 7. Blocked until the migration

| Item | Plan | Unblocked by | Note |
|---|---|---|---|
| Phase C: structural rules as model validation and engine admissibility | CAP | L1 (`check_admissible`), completed at L5 | `TraversalPlan::check_admissible` (on T3) is its first piece |
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
| CR and CR-P | here | lane C, each named |
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
| Q8 | Coverage on a tandem model | per segment, one output per segment; `D_gene` means every D segment |
| Q9 | Summing counters under the chunked reduction | chunk slots with sparse per-gene arrays, allocated on first touch: the FIXME's objection is the size of a dense copy per chunk |
| Q10 | Counters and rejected reads (CX11) | per-read outputs for every read, summaries over accepted reads only; to decide |
| Q11 | `Errors_counter` | keep it, with a CLI flag and a golden |
| Q12 | Pgen's two broken modes (CX6) | fix the sequence mode (rows, with σ); drop the per-σ mode that omits σ |
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
  moves to CT2a; the Phase C, D and E rows point here.
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
- **CX3** may write past an array on long reads with a palindromic J: repair it early, whatever the
  migration's pace.
- **T3's line numbers** move when it merges.
- **AA Pgen** (`origin/feature/AA_PGEN`) edits the leaf, the error model and `iterate()`: Q13.
- **The synthetic tandem model** exercises the computation only: it says nothing about the aligner
  on real tandem reads, nor about repertoires that mix VDJ and VDDJ.
