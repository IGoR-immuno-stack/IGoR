# Error model rewrite plan: `Single_error_rate` first, then the hypermutation models

**Created**: Oct 16 2026
**Scope**: the legacy error models, `Error_rate` and its three subclasses in `src/igor/Model/Legacy/`;
their call sites (the leaf in `Rec_Event.cpp`, the bounds written by `Gene_choice` and `Deletion`,
the EM loop in `GenModel`); their file formats. The counters have their own plan ([COUNTER_REWRITE_PLAN.md](COUNTER_REWRITE_PLAN.md)), and
`FastGenerator` is out of scope.
**Parent**: T7 and Phase E in
[REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md).
**Counterpart of**: [ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) and
[GENERATION_REWRITE_PLAN.md](GENERATION_REWRITE_PLAN.md), whose method it follows.
**Toward**: synthesis step 8 ([ARCHITECTURE_SYNTHESIS.md](ARCHITECTURE_SYNTHESIS.md)): the error
model as a Probe event with handlers.
**Scheduled by**: [ENGINE_MIGRATION_PLAN.md](ENGINE_MIGRATION_PLAN.md), which places Part I's steps against the
migration's (its §3) and lists what each needs from the other (its §6). The order across every
plan: [IMPLEMENTATION_GRAPH.md](IMPLEMENTATION_GRAPH.md).
**Line numbers** are those of `feature/tandemD-T7-pinning` (T1 merged in, T3 not). "On T3" marks
`feature/tandemD-T3-plan`; "AA branch" marks `origin/feature/AA_PGEN`, which is not merged.

---

## 0. Summary

The error models are the last part of the inference core still written for one V(D)J layout:
- `Single_error_rate` reads V, D and J at the ids the `Seq_type` enum pins (T7);
- the two hypermutation models find the V, D and J events by name and reach into them through
  `friend` declarations to read the germline template outside the read;
- nothing tests the hypermutation models: no unit case, no regression track.

**Tandem D does not need this plan, apart from T7a.** Milestone 1 runs on the legacy DFS in
`GenModel` (§3.8). From the error models it needs only T7a — `Single_error_rate` charging every
segment that has a mismatch list — and two guards that make the hypermutation models and the
V/D/J-shaped counters refuse a tandem model instead of answering wrongly (§4.1). Everything else
comes after milestone 1, on its own feature branch, in two parts:
- **Part I, `Single_error_rate`.** One engine-agnostic core (§3.5) that a thin legacy adapter and a
  node of the new architecture both call. It ends with Single on the Categorical handler in the new
  `Topology`, and the legacy DFS calling the same arithmetic.
- **Part II, the hypermutation models.** A regression track and a unit suite first, then one
  context model with two parameterizations, flank segments for the context outside the read,
  sound bounds, and the same core-and-adapter shape as Single.

**Decisions taken (Quentin, Oct 15–16 2026):**
- Tandem D first, on the legacy DFS. The ports to the new handlers and engines, and the
  hypermutation repairs, come after milestone 1, on a separate feature branch.
- `Single_error_rate` before the hypermutation models.
- T7 splits: T7a (Single counts over the ordering, plus the two guards) before T5; T7b (the
  hypermutation models by name, unpinning) with Part II, since it needs their tests.
- **The shape of the port is a shared core (§3.5, option (c))**: an engine-agnostic
  `SingleErrorModel` in `src/igor/Model/`, called by a thin legacy `Single_error_rate` adapter and
  by the new-architecture node. Single does not become an event of the legacy queue.
- Single's accumulation is a confirmed defect (X13): the rate is the posterior-weighted sum of
  errors over the sum of positions, computed in `update()`.
- The template context outside the read comes from flank segments in the registry (B3), not from an
  accessor on the events.
- Viterbi and threshold gating belong to the inference engine, not to the error model.
- `pruning_mismatch_floor` is the input of the error-model upper bounds in the tree search.

**The gate.** Every E or H step is bitwise on the five TRB regression tracks and, once H0 exists, on
the hypermutation track. A confirmed defect is pinned by a `[!shouldfail]` case asserting its
intended value (iterate plan §6.1), carried unchanged through the refactor, and fixed in an ER or HR
commit, once, with the outputs it moves named.

---

## 1. What the code does today

### 1.1 The interface

`Error_rate` ([Errorrate.h:111](../src/igor/Model/Legacy/Errorrate.h#L111)) mixes five jobs:

| Job | Methods |
|---|---|
| The leaf factor | `compute_scenario_error_probability(query, model, scenario, exploration)`, called at [Rec_Event.cpp:406](../src/igor/Model/Legacy/Rec_Event.cpp#L406); the deprecated `compare_sequences_error_prob` |
| EM | `update()`; `add_to_norm_counter()`, `clean_seq_counters()` (end of a read); `clear_accumulators()` (end of a chunk); `copy()`, `add_checked()`, `add_to_err_rate()` (ordered merge); `norm_weights_by_seq_likelihood()`; `is_updated()` / `update_value()` (the "fixed" flag) |
| Read summaries | `get_seq_likelihood()`, `get_seq_probability()`, `get_seq_mean_error_number()`, `get_number_non_zero_likelihood_seqs()`, `model_log_likelihood` |
| The pruning bound | `get_err_rate_upper_bound(k, m)` ([Errorrate.h:196](../src/igor/Model/Legacy/Errorrate.h#L196)) over a lazily grown matrix (`upper_bound_entry`, `build_upper_bound_matrix`) |
| Generation and I/O | `generate_errors()`, `write2txt()`, `type()`; `set_viterbi_run()`; `initialize(Events_map)` |

### 1.2 Lifecycle in one EM iteration

[GenModel.cpp](../src/igor/Inference/Legacy/GenModel.cpp):
- 218-221: a master copy of the error model, initialised against the master events; only a merge
  target.
- 299-303, 401: per thread, per iteration, the `Model_Parms` copy carries an error-model copy;
  `initialize(events_map)` binds it to this thread's events.
- 548-605, per read: the leaf calls the error model; 571 normalises the read's marginals by
  `seq_likelihood`; 596-604 keeps the read (`add_to_norm_counter`) or drops it
  (`clean_seq_counters`) on the mean error count.
- 609-638: each chunk drains into a fresh copy; chunks merge in slot order.
- 668: `update()`.

### 1.3 `Single_error_rate`

[Singleerrorrate.cpp](../src/igor/Model/Legacy/Singleerrorrate.cpp):
- 130-148: the sizes of the V, D (if present) and J sequences and mismatch lists, by enumerator.
- 151-154: the factor `P_scenario · (r/3)^k · (1−r)^(n−k)`, n the genomic nucleotides, k the
  mismatches.
- 160: it applies the leaf's threshold itself and returns 0 below it.
- 162-176: per read, Σ P·k/n and the likelihood; in viterbi mode it overwrites instead.
- 316-323, 329-336: the M-step is the mean over reads of each read's posterior k/n.

All five TRB regression tracks use it.

### 1.4 The hypermutation models

Two models over a table p(N-mer context) of mutation probabilities:
- `Hypermutation_global_errorrate`: p from a logistic, `p = s·mu/(1+s·mu)` with
  `s = exp(Σ e_i(nt_i))`; Newton M-step ([Hypermutationglobalerrorrate.cpp:1154](../src/igor/Model/Legacy/Hypermutationglobalerrorrate.cpp#L1154)).
- `Hypermutation_full_Nmer_errorrate`: p stored per context; ratio M-step with a median fallback
  below a count threshold ([HypermutationfullNmererrorrate.cpp:1035](../src/igor/Model/Legacy/HypermutationfullNmererrorrate.cpp#L1035)).

Both assemble the whole sequence, read offsets and mismatch lists by enumerator, find the V, D and
J events by name in `initialize()`, and take N-mer context from the V template before the read and
the J template after it, through `current_realization_index` (friend) and the public
`alignment_offset_slot()` / `deletion_value_slot()`. `learn_on` gates which segments feed the EM
counts; `apply_to` is never read. Details, duplication and bounds: §3.6.

### 1.5 Mismatch lists and the bound slots

- `Gene_choice` points the segment's mismatch list at the alignment's
  ([Genechoice.cpp:289](../src/igor/Model/Legacy/Genechoice.cpp#L289)) and writes
  `surviving_core_bound` into the segment's own key of the downstream bound map (304-307, the
  bound at 680-699).
- `Deletion` trims the list and rewrites the key once the opposite end is settled, 1.0 before
  ([Deletion.cpp:412-424](../src/igor/Model/Legacy/Deletion.cpp#L412)).
- `Insertion` and `Dinucl_markov` write no list; they write the junction's key
  ([Insertion.cpp:225](../src/igor/Model/Legacy/Insertion.cpp#L225),
  [Dinuclmarkov.cpp:286](../src/igor/Model/Legacy/Dinuclmarkov.cpp#L286)).
- The error model reads the lists at the leaf; the counters read them through the `Scenario`
  snapshot.
- `pruning_mismatch_floor` is built ([GenModel.cpp:326](../src/igor/Inference/Legacy/GenModel.cpp#L326))
  and never written or read.

### 1.6 What a mismatch is in AA Pgen

[AA_PGEN_IMPLEMENTATION_PLAN.md](AA_PGEN_IMPLEMENTATION_PLAN.md) defines two tracks per nucleotide
(159-166):
- **floor**: `!comp_nt_int(gene_nt, iupac_union[p])`, a mismatch in every branch of the motif;
- **upper**: `empty_isect[p] || !comp_nt_int(gene_nt, iupac_intersection[p])`, a mismatch in at least
  one branch.

Floor ⊆ upper, and the two coincide for an exact NT query (170-172). The pruning bound is to use
the floor, the leaf the upper track (508-510). In motif mode the leaf counts **patches**: one patch
is one error, and `genomic_nucl` is the number of patches (1006-1008; AA branch
`Singleerrorrate.cpp:131-216`).

Three inconsistencies on the AA branch, which this plan's interface must not inherit:
- its fast path treats an upper-track hit as a confirmed mismatch, which it is not;
- the bound counts nucleotides while the leaf counts patches: for r > 0 the bound falls below the
  leaf value, and pruning is no longer conservative;
- Met and Trp codons make no patch, so a mismatch there is never penalised.

### 1.7 Generation

- The legacy generator calls `generate_errors()` on the assembled sequence (SequenceGenerator.cpp).
- `SamplingEngine` samples one index per node and never calls `sampleSequence`
  ([SamplingEngine.tpp:105-146](../src/igor/Generation/SamplingEngine.tpp#L105)): it builds no
  sequence-valued content at all, Dinucl_markov's included, so it has no error step either.

---

## 2. Defects and oddities

Categories follow the iterate plan's §6.1: *confirmed* defects are pinned at their intended value
under `[!shouldfail]`; *odd* behaviour is pinned as observed, with a comment naming what is
undecided.

### Part I (`Single_error_rate` and the shared machinery)

| # | Finding | Category | Owner |
|---|---|---|---|
| X1 | Single charges ids 0, 2 and 4 only: a tandem model's D1 and D2 mismatches are free | confirmed | T7a |
| X2 | Tandem: D1 and D2 bound slots are below 1 while the leaf charges no D, so the bound product can fall below the leaf value — pruning is unsound | confirmed | T7a; E3 for a segment left out of `apply_to` |
| X4 | Threshold gating inside the error model: Single gates, the hypermutation models do not | confirmed (engine's job) | E2 |
| X5 | Viterbi inside the error model; the hypermutation models keep the last leaf, not the best | confirmed (engine's job) | E2 |
| X13 | Single's M-step averages per-read rates; the EM estimate is Σ posterior errors / Σ posterior positions | confirmed | ER1 |
| X15 | `viterbi_run` uninitialised in the base constructor | confirmed | E2 (viterbi leaves the error model) |
| X20 | JSON exports only Single and reads no error model; the new engine path never infers the rate (test_ModelInference.cpp:567-586 discards the updated copy) | gap | E4 |
| X22 | Dead code: Single's deprecated path (throws on VJ), the base default, `subseq_compare_err_num` | dead | E1 |

### Part II (the hypermutation models)

| # | Finding | Category | Owner |
|---|---|---|---|
| X3 | The bounds are not upper bounds (§3.6) | confirmed | HR |
| X6 | `fix_err` ignored: neither `update()` checks `is_updated()`, so evaluate refits | confirmed | HR |
| X7 | Leaks: `copy()` allocates five 4^N arrays per chunk, `initialize()` the template arrays; `operator+=` by value shares the raw pointers | confirmed | H2 |
| X8 | The J-tail template index is `i−J5'+del`; it should be `i+h−J5'+del` | confirmed | HR |
| X9 | The V head reads a negative template index when the read starts at V position 0 | confirmed | HR |
| X10 | `apply_to` is parsed, written and compared, never applied | odd | H4 |
| X11 | The full N-mer count threshold is not saved; `write2txt` keeps 6 significant digits | odd | HR |
| X12 | `type()` strings differ from the file tags; the full N-mer reader's errors name the global model | cosmetic | H2 |
| X14 | `update()` never invalidates the bound cache (latent: per-iteration copies hide it) | latent | H2 |
| X16 | Iterator comparisons across three different vectors (undefined behaviour) | confirmed | H2 |
| X17 | Tandem: no `D_gene_seq` event, so the sequence is assembled as V+J and read positions misalign | confirmed | refused by T7a's guard; fixed by H4 |
| X18 | Generation draws errors on insertions, the likelihood charges genomic positions only | odd | Q5 |
| X19 | `generate_errors` is const but writes `mutable` scratch: not thread-safe | confirmed | H2 |
| X21 | Junction segments have no mismatch list (the iterate plan's open defect) | known | iterate plan |
| X23 | About 750 duplicated lines between the two classes (§3.6) | duplication | H2 |

---

## 3. Analyses

### 3.1 The interface under the handlers

| Legacy | Becomes |
|---|---|
| `compute_scenario_error_probability` | the event's non-branching evaluation. Storing `scenario_error_w_proba`, the threshold test and the read likelihood go to the engine |
| `update` | `maximize()`; for Single as an event, the normalised marginal row |
| `initialize(Events_map)` | discarded: the scope is resolved from names when the model is built |
| `add_to_norm_counter`, `norm_weights_by_seq_likelihood` | the engine's per-read normalisation of every accumulator |
| `clean_seq_counters` | the engine's read policy (mean-error rejection) |
| `clear_accumulators`, `copy` | `AccumulatorSet::empty_copy()`; parameters are shared, not copied |
| `add_checked`, `add_to_err_rate`, `type()` | the ordered accumulator merge; the type check becomes static |
| `write2txt` | `to_json`, and the `@ErrorRate` lines in the legacy translator |
| `get_err_rate_upper_bound`, `upper_bound_entry`, `build_upper_bound_matrix` | the handler's `bound()`, tabulated after each M-step (tier 0b), taking the mismatch floor |
| `generate_errors` | the event's draw / construction step (generation plan) |
| `seq_likelihood`, `seq_probability`, `seq_mean_error_number`, `model_log_likelihood`, `number_seq` | the engine's per-read summary and run log |
| `set_viterbi_run`, `viterbi_run` | discarded: the engine's semiring |
| `is_updated`, `update_value` | the event's `fixed` flag, as for any event |
| `debug_number_scenarios`, the deprecated `compare_sequences_error_prob`, `subseq_compare_err_num` | discarded |

### 3.2 Scope by configuration

Two separate questions:
- **Scope** — which positions are charged (`apply_to`) and which feed the EM counts (`learn_on`).
  - Both become ordered sets of seq-type names, resolved to ids when the model is built.
  - The v1 translator maps the `Gene_class_legacy` tokens (`V_gene`, `VDJ_genes`, …) through the
    gene class of the `Gene_choice` that creates each segment, so `D_gene` is {D1, D2} in a tandem
    model. A v2 or JSON file names any segments.
  - Defaults: `apply_to` = every segment with a mismatch writer (§3.3); `learn_on` = `apply_to`.
  - Single gets `apply_to` too. No error model looks up V, D or J events any more, and the counters
    iterate the same resolved set.
- **Context** — what each position reads: the Phase A footprint's relative `context_reads()`
  (`{toward, anchor, width}`; width 0 for Single, h = (N−1)/2 for an N-mer model). That is what the
  flank segments, the bound and Phase D need; it is not a segment set.

### 3.3 The mismatch-vector role

**Who writes what today.** In the downstream bound map, a segment's own key holds the error model's
bound on that segment, computed from its mismatch list. A gap's key holds the junction-length
factor. In a VDJ model, with T3's layout:

| Key | Holds | Written by |
|---|---|---|
| `V_gene_seq` | V's error bound | V's gene choice; rewritten by the V 3′ deletion once settled |
| `VD_ins_seq` | the V→D gap's length factor | the insertion and its chain |
| `D_gene_seq` | D's error bound | D's gene choice; the two D deletions |
| `DJ_ins_seq` | the D→J gap's length factor | the insertion and its chain |
| `J_gene_seq` | J's error bound | J's gene choice; the J 5′ deletion |
| `VJ_ins_seq` | nothing (1.0): no event places it in a VDJ model | — |
| 6 | the V→J gap, which encloses D | the events taking part in it (R12) |

A gap reuses its single interior segment's id only when that segment has no error factor of its
own: VD and DJ do, V→J encloses D and gets a fresh key.

**Why `SegmentIndex::has_own_factor` matches by coincidence** (on T3, `SegmentIndex.h:90-112`). It
decides "has an error factor" from the A0 roles: the content is modified, or it is created by the
event that also places the ends. Its own comment states the intent: "its content can still change
against the read". That is a proxy. It is right for the genes (an aligned template can mismatch the
read) and for the junctions only because today's Dinucl_markov copies the read, so its content can
never mismatch.

**Where the proxy fails.** A Dinucl_markov running the forward algorithm, summing over inserted
sequences instead of reading them off the read, generates content that can mismatch. It would write
a mismatch list and an error bound at `VD_ins_seq`'s key — the key the plan gave to the V→D gap,
because by roles the junction still has no factor of its own (its ends are placed by the insertion,
its content created by the chain). The two factors would collide in one slot.

**Proposal: an explicit role.** `get_mismatch_role(id)` returns Creates, Modifies or None:
- `Gene_choice`: Creates; `Deletion`: Modifies; `Insertion`: None;
- `Dinucl_markov`: None today, Creates in a forward version.

It has four consumers:
- `has_own_factor` is defined by it, so a forward chain gets an error slot and its gap a fresh key;
- which segments carry a mismatch list and an error slot;
- the default `apply_to`;
- what the counters iterate.

Admissibility checks, in the engine's plan:
- `apply_to` ⊆ segments with a mismatch writer (else the leaf reads a list nobody wrote);
- `learn_on` ⊆ `apply_to`;
- a segment with a writer stays out of `apply_to` only by explicit opt-out, since its mismatches
  would then be free.

What keeps pruning sound is the bound being segment-aware: 1.0 for a segment outside `apply_to`,
the error model's factor bound inside it.

**Two tracks and a unit** (§1.6). A writer records the floor and the upper list, one aliased vector
in NT mode (hazard H5). The bound takes the floor; the likelihood the upper track. The error model
declares its unit — the nucleotide, or the codon/patch with a frame dependency — and the bound counts
in the same unit, which is what the AA branch gets wrong. The `/3` substitution factor is a
parameter, not a constant.

### 3.4 The error model as an event

**What it shares with the events.**
- Like Dinucl_markov, it is non-branching: given the scenario and the read, its factor is
  determined.
- Each `Rec_Event` method has an equivalent: `iterate` (multiply the factor and hand on),
  `draw_random_realization` (`generate_errors`), `write2txt`.
- With handlers, its parameters and accumulators take the same form as any event's.
- Nothing prevents conditioning on a parent: an error rate or a hypermutation model per V gene is a
  [parent, …, 2] tensor, which also tests whether a hypermutation model captures the sequence
  dependence or a gene effect.

**What is not a reason to treat it differently.**
- The missing `SamplingEngine` error step (§1.7): the engine samples no sequence-valued event, and
  will need a construction stage for Dinucl_markov and the error model alike.
- Per-read normalisation, gating and viterbi: those are the engine's, for every event (X4, X5).

**What is specific.**
- It reads the content of every segment in `apply_to`, so the order must place it after all of them
  are built: an admissibility condition, like the chain's after its insertion.
- Upstream events pre-evaluate its bound segment by segment: the error key of a segment *is* the
  error model's partial evaluation on it. Dinucl_markov is bounded the same way per gap, through
  `affects_proba_of` and `span_proba_factor`
  ([Rec_Event.h:488, 523](../src/igor/Model/Legacy/Rec_Event.h#L488)). The error model's
  counterpart takes (floor mismatches, positions) instead of lengths.

**Not in the code.** The per-segment queries this needs — "is this segment charged", "how wide is
the context", "the factor bound on a window with unknowns" — do not exist; only
`get_err_rate_upper_bound(k, m)` does. They map onto `affects_proba_of`, the footprint's
`context_reads()`, and a per-segment counterpart of `span_proba_factor`.

**The new architecture already requires it.** `Topology` stores `legacy::Rec_Event`s and the
handler factories dispatch on `EventType` ([Topology.h:24, 72](../src/igor/Model/Topology.h#L24),
[InferenceHandlerFactory.tpp:37](../src/igor/Inference/InferenceHandlerFactory.tpp#L37)). An error
model outside the events cannot get a handler.

**A class hierarchy for it?** One option is an abstract `Event` with two abstract daughters,
`RecEvent` and `ErrorEvent`.
- *For:* the types state the intent; `ErrorEvent` no longer answers queries that mean nothing for
  it (offset roles, `length_delta`, `affects_length_of`); `Topology` and the handler factories hold
  the base `Event` and stop requiring a `legacy::Rec_Event`.
- *Against:* the variation is not two-way. `Dinucl_markov` is already non-branching and
  sequence-valued, like the error model; `Insertion` creates no content; a forward-algorithm chain
  would write mismatches. The real axes are capabilities — what an event writes, what it reads,
  whether it branches, its domain — which the A0 roles and the footprint already express as
  queries. A class split freezes one of them as the one that counts.
- *It recreates the old dichotomy* as soon as engines must tell the two daughters apart: an
  `if (ErrorEvent)` in an engine is today's special call at the leaf, moved.
- Splitting the class now would split the one class both architectures share (§3.9).

**Recommendation.** The error model is an event with the Probe role, described by the same
capability queries as every other event, not by a class of its own:
- in the new architecture, one `Event` interface (identity, parents, domain, footprint of writes
  and reads, branching flag, bound participation, apply, I/O) that engines only ever query.
  Convenience bases supplying defaults — a `ProbeEvent` that writes nothing and never branches —
  are fine as long as nothing dispatches on them;
- in the legacy code, no split and no move into the queue: the legacy DFS keeps its leaf call until
  `LegacyDfsEngine` replaces it (§3.5).

### 3.5 The shape of the port: a shared core (decided Oct 16 2026)

Three shapes were weighed for putting Single on the new handlers.

- **(a) An event of the legacy queue.** `Single_error_rate` becomes a `Rec_Event` appended last to
  `Model_Parms`' event list; its rate is a two-cell row [r, 1−r] in `Model_marginals`, filled by
  `add_to_marginals`, so that the machinery every event has (per-read normalisation, the threshold,
  read rejection, viterbi, the ordered merge, the M-step normalisation) does the rest. Rejected: it
  changes the marginals layout and files, rewires `GenModel` and the leaf, is not bitwise (marginals
  and the error model accumulate in different types and orders), and it is a large change to code
  that `LegacyDfsEngine` will replace.
- **(b) A node in the new `Topology` only.** The legacy class is unchanged; `LegacyBridge` appends an
  error node and fills its accumulator from the legacy `Error_rate`. Rejected: two representations of
  one model, and every repair (X4, X5, X13) made twice.
- **(c) A shared core — chosen.** An engine-agnostic `SingleErrorModel` in `src/igor/Model/`, on
  the pattern of `SegmentIndex`, `TraversalPlan` and `FlankCheck`, with small, mostly stateless
  operations:
  - the scope: `apply_to` names resolved to segment ids;
  - the factor (r, k, n), and the per-segment bound (r, floor k, m);
  - one scenario's contribution to the accumulator (P·k, P·n);
  - the M-step: r = Σk / Σn, with a guard for an empty row;
  - errors applied to a sequence, given an RNG.

  Both sides call it:
  - **legacy:** `Single_error_rate` keeps the `Error_rate` interface and becomes a thin adapter. It
    reads the segments from `ScenarioContext` and calls the core: one level of indirection, no
    duplicated arithmetic, bitwise as long as the operations keep their order;
  - **new architecture:** a node holding r in a two-cell tensor. The Categorical handler's M-step,
    which normalises the last axis
    ([CategoricalInferenceHandler.tpp:77-98](../src/igor/Inference/CategoricalInferenceHandler.tpp#L77)),
    is the core's update once X13 is fixed: [E[errors], E[matched positions]] normalises to
    [r, 1−r]. A [parent, 2] shape gives one rate per parent realization.

  With (c), (a) is never needed: the legacy DFS keeps its leaf call until `LegacyDfsEngine`, which
  treats the error node as a node like any other. What remains of (b) is the bridge carrying r
  between the legacy object and the tensor, as it already carries the marginals.

What the node needs:
- an `EventType` value and a `Rec_Event` subclass (`Topology` holds nothing else, §3.4);
- inference registrars for `double` and `long double`;
- a sampling registrar: `SamplingEngine` throws on an unregistered type
  ([SamplingHandlerFactory.tpp:21](../src/igor/Generation/SamplingHandlerFactory.tpp#L21));
- a JSON registrar, and `LegacyBridge` carrying the node (today it copies only
  `get_event_list()`, [LegacyBridge.tpp:25](../src/igor/Model/LegacyBridge.tpp#L25));
- a guard for an all-zero row, which would set 1−r to 0 and zero every scenario;
- a hook after the M-step that rebuilds the bound table.

### 3.6 The hypermutation models (Part II)

**Duplication.** The two classes share, line for line apart from the class name and the inline
sequence assembly (plus two commented-out blocks in the global model):
- the E-step: global 397-1063, full N-mer 431-957;
- `initialize()`: global 1510-1636, full 1085-1201;
- `generate_errors()`: global 1065-1118, full 958-1011;
- `add_to_norm_counter()`, `clean_seq_counters()`, `clear_accumulators()`.

That is about 750 lines, and both likelihoods read one table, `Nmer_mutation_proba[4^N]`.

**Factoring.** One N-mer context model — the context walk, the likelihood, the [4^N, 2]
accumulation, the scope, generation — over the table p(context), with two parameterizations:
- **Logistic**: p from mu and e_i; Newton M-step.
- **Categorical**: p stored; ratio M-step with the median fallback and its count threshold.

These are the synthesis' two handlers.

**Why the bounds are not bounds.** The bound must not be smaller than the factor the leaf will
apply: per position, p(context)/3 for a mismatch and 1 − p(context) for a match.
- The global model uses p̄ = mu/(1+mu), the context in which every e_i is 0
  ([Hypermutationglobalerrorrate.cpp:322, 343](../src/igor/Model/Legacy/Hypermutationglobalerrorrate.cpp#L322)).
- The full N-mer model uses the median p
  ([HypermutationfullNmererrorrate.cpp:316-368](../src/igor/Model/Legacy/HypermutationfullNmererrorrate.cpp#L316)).

A mismatch at a context with p > p̄ costs p/3 > p̄/3; a match at a context with p < p̄ costs
1 − p > 1 − p̄. Either way the true factor exceeds the "bound", and a subtree whose best scenario
is above the threshold can be pruned. Both files carry the sound version commented out.

**What should be done.**
- Minimal: (Pmax/3)^k (1 − Pmin)^m.
- Tighter: the exact factor wherever the whole window is known (template-internal positions), and
  max/min over the completions only where a window reaches an unresolved neighbour.
- The floor count for IUPAC queries (§3.3).

The cost in pruning is measured on the hypermutation track (Q6).

**Context outside the read.** Flank segments (B3): `Gene_choice` writes the template strip before
the read (V) or after it (J) into a flank segment, and the model reads it through the scenario
state. The friend declarations go. Whether flanks have ends the TraversalPlan can check, or are
exempt as content-only segments, is Q7.

### 3.7 Enumerator removal

Occurrences outside comments and string literals, measured Oct 16 2026:

| File | `Seq_type` | gene class | Removed by |
|---|---|---|---|
| `Singleerrorrate.cpp` | 14 | 0 | T7a (the leaf), E1 (the dead path) |
| `Hypermutationglobalerrorrate.cpp` | 39 | 27 | T7b, H2–H4 |
| `HypermutationfullNmererrorrate.cpp` | 45 | 27 | T7b, H2–H4 |
| the four counters, `Scenario.h`, `EventUtils.cpp` | 89 | 3 | CT2a and CT3 of the counter plan; T7b for the enum `build_scenario_sequence()` |

End state: `Gene_class_legacy` only in the v1 reader and writer.

### 3.8 Tandem D on the legacy DFS

Milestone 1 runs on `GenModel`'s DFS, and every tandem step does already:
- T1 (loading) and T3 (the traversal plan, which `GenModel` owns) are written for it;
- T4 is done for the legacy generator, the one T5 generates with;
- T2 (colliding names) can be fixed in the legacy model too: the names key legacy structures
  (`edges`, `Index_map`, `offset_map`), and giving the seq type to the generated name of non-legacy
  segments only leaves every legacy name unchanged.

From this plan, milestone 1 needs T7a and its two guards (§4.1). Where the legacy design does, or
might, stand in the way:

| Concern | On the legacy DFS | Would the new architecture help? |
|---|---|---|
| Correctness | nothing known once T7a and T2 are in | no |
| Cost: scenarios explored, and the fold's build time on V-D1-D2-J | unmeasured; T5 measures the ratio against VDJ. The fold is state by state since R12, so it is polynomial | only if the DFS itself is too slow: a DP engine (synthesis step 9), milestone 2 |
| An absent D (B10) | needs a layout DAG and adjacency per scenario | **yes** — the one real blocker, and it is milestone 2 |
| The V/D/J-shaped counters | degraded, so refused on a tandem model (§4.1) | CT2a and CT3 of the counter plan, on the legacy DFS, for best scenarios and Pgen; coverage waits for the scenario record (the migration's L3); `Errors_counter` is dropped |

### 3.9 Running both architectures at once

The tension is real:
1. **One class, two lifecycles.** `legacy::Rec_Event` is both the DFS event and the `Topology` node.
   Every method the new architecture adds is implemented on the legacy events, and the DFS state
   (memory layers, `iterate`) lives in the node class.
2. **Two homes for parameters.** `Model_marginals` and the `RecombinationModel` tensors, synced by
   `LegacyBridge`. The layouts must stay identical; they are (parents first, own axis last).
3. **Not bitwise by default.** The legacy run merges chunks in a fixed order (#68); `InferenceEngine`
   combines accumulators by name in no fixed order. The file formats differ too (text v1/v2, JSON).
4. **Repairs land twice** unless the two paths share the code: a fix like X13 must be identical on
   both, or comparing them breaks.
5. **The error model is the worst case:** outside the events in the legacy run, a node in the new
   one.

**Rules while both live:**
- Logic goes into engine-agnostic components in `src/igor/Model/`, which depend on neither
  architecture; the legacy code and the new nodes both call them (§3.5, §3.10).
- The bridge is the only place where state is copied between the two.
- A parity test compares the two paths wherever both run (as `test_NormalizationCompareLegacy`
  does for the categorical M-step).
- New-architecture code never includes `GenModel` or a legacy event's DFS internals.

The migration plan adds its own rules for the steps that move the runtime — no runtime flag, a
full-precision A/B oracle — in its §5, and the order in which this plan's steps and its own meet
in its §3 and §6.

### 3.10 The same pattern for the other events

The core-and-adapter shape applies to the four recombination events too, but piece by piece, not as
a parallel class per event:
- **Their parameters already have one home per architecture**, the marginals and the tensors, and
  the handlers' M-steps already match the legacy normalisation. There is no duplication to remove
  there.
- **Their class is already shared** (§3.9, point 1), so a second, non-legacy class per event would
  double the class tree during the transition.
- **What is legacy-specific is the DFS control flow** in `iterate()` — memory layers, recursion,
  pruning — and that belongs to `LegacyDfsEngine`, not to the events.
- **What is engine-agnostic is extracted when a second consumer appears**, as it has been:
  `JunctionGeometry` and `PendingModifierBounds`, the span profiles and the fold, `SegmentIndex`,
  `TraversalPlan` and `FlankCheck` (all in `src/igor/Model/`), and `construct_realization()`, which
  the generation plan split out as the event's `apply`.

Candidates the port will need next:
- **The deletion of a template.** `Deletion` builds a trimmed or palindromic segment twice: in
  `iterate()` on `Int_Str` ([Deletion.cpp:334-391](../src/igor/Model/Legacy/Deletion.cpp#L334)) and in
  `construct_realization()` on `std::string` (536-590). One function, generic over the alphabet,
  serves both and any new engine.
- **The mismatch recomputation** that goes with it (`trim_mismatches`, the palindrome's new
  positions), which belongs with the mismatch role (§3.3).
- **`Gene_choice::surviving_core_bound`**, which is the error model's per-segment bound and moves to
  the core of §3.5.
- **The feasibility tests** of each realization (overlap, read bounds), which every engine that
  enumerates needs.

---

## 4. Steps

### 4.1 On the tandem path (milestone 1)

| Step | Content | Bitwise |
|---|---|---|
| T7a | Single counts every segment that has a mismatch list, over the ordering; its `[tandem_d]` case is written first, `[!shouldfail]` (capability plan) | yes |
| T7a guards | The two hypermutation models and the four V/D/J-shaped counters refuse, at initialisation, a model whose ordering is neither V-D-J nor V-J, with a message naming the segment they cannot place; one case each | yes |

### 4.2 Part I — `Single_error_rate`, after milestone 1, on its own branch

| Step | Content | Bitwise |
|---|---|---|
| E0 | Characterization suite for Single. Harness: `set_error_rate(std::shared_ptr<Error_rate>)` (today `set_error_rate(double)` only, test_utils.h:434-439). Cases: the leaf factor (VDJ, VJ, tandem); gating; viterbi; per-read and cross-read accumulators; the M-step (observed value, and X13's intended value under `[!shouldfail]`); merge order; `copy`; the round trip (v1, v2, JSON export); the bound; seeded `generate_errors`. Coverage report and mutations; ITERATE_TEST_GUIDE.md §7 has the recipe | tests only |
| E1 | The core `SingleErrorModel` in `src/igor/Model/` (§3.5): scope, factor, contribution, M-step (still the observed formula), bound, generation. `Single_error_rate` becomes its adapter; `surviving_core_bound` moves into it; dead code out (X22) | yes |
| ER1 | X13 in the core: Σ errors / Σ positions, normalised per read | every inference golden moves, measured once |
| E2 | Gating and viterbi leave the error model: the leaf evaluates the factor, tests the threshold, then asks the error model to accumulate; viterbi is the engine's (X4, X5, X15). The hypermutation models keep their current behaviour behind a legacy policy flag until Part II | yes |
| E3 | `get_mismatch_role`; `has_own_factor` defined by it; the segment-aware bound taking the floor count (the upper list in NT mode); Single's `apply_to` | yes on shipped models |
| E4 | The node in the new architecture, calling the core: `EventType`, registrars, `LegacyBridge` carrying r, zero-row guard, bound rebuilt after the M-step. test_ModelInference generates reads with errors and infers the rate; a parity case compares the node's M-step with the adapter's | new tests |

### 4.3 Part II — the hypermutation models, after milestone 1

| Step | Content | Bitwise |
|---|---|---|
| H0 | Hypermutation regression track (below) | creates goldens |
| H1 | Unit suite: the leaf (an independent Python oracle for the factors; global(3) = full(3) at the leaf), `update` (Newton, ratio/median), merge order, `copy`, round trip, seeded `generate_errors`, the bound, viterbi, `fix_err`; confirmed defects `[!shouldfail]`; coverage and mutations | tests only |
| T7b | The hypermutation models read by name, a tandem model is refused, the sequence is assembled in ordering order, the pinning goes (capability plan). The counters' half is CT2a of the counter plan | yes, TRB + hypermutation |
| H2 | One context model with two parameterizations (§3.6): `std::vector` storage (X7), `operator+=` by const reference, no `mutable` scratch (X19), tags (X12), cache invalidation (X14), iterators (X16) | yes |
| H3 | Flank segments (B3) replace the friend access | yes |
| H4 | Scope sets (§3.2); `apply_to` applied (X10, a no-op on shipped files, which all say `VDJ_genes`) | yes on shipped files |
| HR | Repairs, one per commit: sound bounds (X3, Q6), `fix_err` (X6), J tail (X8), V head (X9), precision and the saved threshold (X11) | moves the hypermutation goldens |
| H5 | The two parameterizations as handlers: the context model of H2 becomes the core, with a legacy adapter and a node, as Single's in §3.5 | bitwise against the post-HR goldens |

**The hypermutation track (H0).** No BCR reads ship with the repository, but the hypermutation
models generate, so the track seeds itself:
- **Model**: the IGH parms plus the global 5-mer V model
  (`models/human/bcr_heavy/supplementary_models/ind_5mer_V_hyperm_parms.txt`; its 432 event lines
  are the base IGH model's, so the base IGH marginals pair with it), copied into
  `scripts/tests/data/input/` so that submodule bumps do not move the goldens.
- **Reads**: generated once with the legacy generator, a fixed seed and one thread, then frozen with
  their alignments. Trimmed by at least h at both ends (20–40 nt at the 5′ end is realistic): a read
  starting at V position 0 makes the V head read before its template (X9).
- **Runs**: evaluate and 2 iterations on the V model; one iteration each on the D- and J-learning
  models, whose learning loops are separate code; a full N-mer 3-mer fixture, computed to 17
  significant digits from the 3-mer global file's probabilities, evaluated and inferred once.
- **Budget**: measured from 50 to 300 reads. IGH has 97 V and 35 D genes, and the hypermutation
  leaf is O(read length). Stay at 5-mer or below: each `copy()` leaks five 4^N arrays (X7).
- Summaries that depend on thread order are excluded.

---

## 5. Decisions

| # | Question | Recommendation |
|---|---|---|
| Q1 | The shape of the port: (a) an event of the legacy queue, (b) a node of the new `Topology` only, or (c) a shared core | **decided Oct 16 2026: (c)**, §3.5 |
| Q2 | *(Dropped with (c): no step moves Single into the legacy queue, so no summation order changes.)* | — |
| Q3 | The rate when a read set gives an all-zero row | keep the previous value |
| Q4 | The bound's unit in AA mode: patch or nucleotide | the leaf's unit |
| Q5 | Errors on insertions in generation (X18) | to decide with the generation plan |
| Q6 | The sound hypermutation bound, and the pruning it costs | measure on H0 |
| Q7 | Flank segments' ends in the TraversalPlan, or an exemption for content-only segments | with B3 |

---

## 6. Hazards

- **AA Pgen** (`origin/feature/AA_PGEN`, not merged). Its `Single_error_rate` patch mode loops over
  ids 0..5 (AA branch `Singleerrorrate.cpp:137`), another reader of the pinning; it treats the
  upper track as confirmed and counts the bound and the leaf in different units (§1.6). Reconcile
  when AA Pgen resumes, against §3.3.
- **T3 and T1**: E3 needs `SegmentIndex` and the TraversalPlan, which are on T3; this document's
  line numbers move when T3 is merged.
- **Two architectures at once** (§3.9): every step of Parts I and II keeps the legacy DFS as the
  regression oracle, so a port is bitwise against it or names what it moves.

## 7. Out of scope

- The counters: their tests, their move off the enum and their generalisation are
  [COUNTER_REWRITE_PLAN.md](COUNTER_REWRITE_PLAN.md).
- B10, the absent-segment semantics.
- `FastGenerator`, beyond noting X18 and X19.
- `SamplingEngine`'s error step: it waits for the engine's construction stage, which
  Dinucl_markov's sequence needs first ([GENERATION_REWRITE_PLAN.md](GENERATION_REWRITE_PLAN.md)
  §6). The legacy generator applies errors through the core from E1 on.
