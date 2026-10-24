# Scenario exploration internals: registry, scenario state, plan, safety, bounds

**Status**: reference documentation of the code on `feature/tandemD-T5` @ `3c06f52` (Oct 23 2026):
the tandem-D trunk `feature/tandemD` and the branches stacked on it (T2, T7a, T5).
**Audience**: developers who work on the legacy inference engine, `GenModel::infer_model()` and the
`iterate()` bodies of the four events, or who add an event type or a model shape to it.
**Companion**: [REC_EVENT_CONTRACT.md](REC_EVENT_CONTRACT.md) lists what an event implements and
exposes. This document explains the machinery that uses it.
**Not covered**: the Model module's `RecombinationModel`, `Topology`, handlers, `SamplingEngine` and
`InferenceEngine`. That code is merged, but neither inference nor generation calls it yet
(MIG §1.1).

This document describes the code as it is. The plans record how it got there and what comes next.
They are cited by the short names [IMPLEMENTATION_GRAPH.md](../IMPLEMENTATION_GRAPH.md) uses:

| Short | Document |
|---|---|
| CAP | [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](../REC_EVENT_CAPABILITY_REFACTORING_PLAN.md), tandem D (T1–T7, B0–B11) |
| ITER | [ITERATE_GENERIC_REWRITE_PLAN.md](../ITERATE_GENERIC_REWRITE_PLAN.md), the generic `iterate()` bodies and phase R |
| PBM | [PROBA_BOUND_MACHINERY.md](../PROBA_BOUND_MACHINERY.md), the history and measurements of the junction bound |
| MIG | [ENGINE_MIGRATION_PLAN.md](../ENGINE_MIGRATION_PLAN.md), steps L0–L8 |
| ERR | [ERROR_MODEL_REWRITE_PLAN.md](../ERROR_MODEL_REWRITE_PLAN.md), steps E0–E4, H0–H5 |
| CNT | [COUNTER_REWRITE_PLAN.md](../COUNTER_REWRITE_PLAN.md), steps CT0–CT7 |
| GEN | [GENERATION_REWRITE_PLAN.md](../GENERATION_REWRITE_PLAN.md) |
| SYN | [ARCHITECTURE_SYNTHESIS.md](../ARCHITECTURE_SYNTHESIS.md), the target architecture |

**Subject to change.** A part that a planned step changes carries a note naming that step:

> **Subject to change — MIG L3.** What the step changes.

No step is planned for an unmarked part beyond moves and renames. MIG L8 moves the `iterate()` bodies
to `Inference/Dfs/` as the engine's kernels, and retires `GenModel` and the legacy types.

---

## 1. The walk, and what runs when

Inference explores the recombination scenarios of each read depth first. Each level of the tree is
one event, in the order of the model queue (`Model_Parms::get_model_queue()`): at each step, the
highest-priority event whose parents are all queued, ties broken by the order of the event list. A
JSON model sets that list's order with `engine_hints.legacy_dfs.event_order`, or by nickname when it
has no hint. A node is one realization of its event. A leaf is a complete scenario.

The work is split by how often its inputs change:

| When | Where | What |
|---|---|---|
| once per model | `Model_Parms::finalize()` | the registry is frozen, and each event learns its segment id and its neighbours (§2.2) |
| once per run, before any thread starts | `GenModel::infer_model()` | the `TraversalPlan` is built from the queue order and checked (`check_admissible()`, §4). Then `Error_rate::check_admissible()` and each counter's run (§7) |
| each EM iteration, each thread | prologue of the parallel region | the model is copied. The scenario maps are sized from the frozen registry (§3). Every event is given the plan. `initialize_event()` runs in queue order. The next-event array, the error rate and the counters are set up |
| each EM iteration, one thread | `#pragma omp single` | `initialize_Len_proba_bound()` runs in reverse queue order and folds the junction tables (§6.3). The other threads copy the result with `adopt_Len_proba_bound()` |
| each EM iteration, each thread | — | `update_event_internal_probas()`, which builds `Dinucl_markov`'s dinucleotide matrix |
| each read | the read loop | the first event's `iterate()` runs on the read and its alignments |
| each node | `Rec_Event::iterate()` | for each realization: its probability, its writes, the feasibility checks, the downstream bounds, the pruning test, then `iterate_wrap_up()` |
| each leaf | `Rec_Event::iterate_wrap_up()` | the error rate, the threshold, the counters, `add_to_marginals()` |

`iterate_wrap_up()` hands off to `next_event_ptr_arr[event_index]->iterate()`. The pointer is null
after the last event, and the leaf code runs instead.

Every node receives the same five contexts, by reference:

| Context | Holds | Written during the walk |
|---|---|---|
| `QuerySequenceContext` | the read, its integer form, its alignments by gene class (library) | no |
| `ModelContext` | the marginals, the offset map (parent strides), `events_map`, the queue | no |
| `ScenarioContext` | the scenario probability, the constructed sequences, the offsets, the mismatch lists | yes, one layer per event (§3) |
| `ExplorationContext` | the downstream bound map, the best scenario so far and the threshold factor, the index map, the next-event array, the safety matrix, the pruning mismatch floor | yes |
| `AccumulationContext` | the marginals being accumulated, the counters, the error rate | at the leaf |

> **Subject to change — MIG L1, L2, L6.** L1 moves the iteration and read loops into an EM driver
> and a `LegacyDfsEngine`. The engine's `check_admissible()` runs the plan's checks and T7a's. L2
> replaces `iterate_wrap_up()` and the next-event array with an engine continuation, and the leaf
> becomes the engine's. L6 places each set-up step by tier: tier 0a (the copies,
> `initialize_event()`, the plan) once per run, and tier 0b (the fold, the dinucleotide matrix, the
> bound tables) once per iteration, shared and const. ITER §2.10 plans a per-read hook,
> `prepare_for_query()`. It is not built yet.

---

## 2. Segments: the registry and its ids

### 2.1 `SeqTypeRegistry`

A model's segments have names: `V_gene_seq`, `VD_ins_seq`, …, or `D1_gene_seq`, `D1D2_ins_seq` in a
tandem model. The registry ([SeqTypeRegistry.h](../../src/igor/Core/SeqTypeRegistry.h)) gives each
name a dense `SeqTypeId` (`uint16_t`), and holds the **ordering**: the segments from 5′ to 3′.

- **Name and id.** `register_type()`, `id()` and `name()` convert between them. Ids are dense from
  0, so they index arrays directly. The name identifies a segment in files and in messages. The id
  identifies it at run time. The walk never looks a name up.
- **Ordering and neighbours.** `ordering()` lists the ids 5′ to 3′. `left_neighbor()` and
  `right_neighbor()` are array lookups. They return `kNoSeqType` at the ends of the ordering and for
  a type the ordering does not hold.
- **Freezing.** After `freeze()`, registering a new name throws. A map sized from the registry
  requires a frozen one (the `DynamicSequenceMap` constructor checks), so the id space cannot grow
  under a map that is already allocated.
- `contains()` and `index_of()` look at the ordering, not at the id space.

### 2.2 `Model_Parms::finalize()`

It runs once per model, after the events are added:
1. It refuses two events of the same type on the same segment and side. The legacy `Events_map` is
   keyed by `(Event_type, seq_type, Seq_side)` and would silently keep only one of them.
2. It pins the six legacy names to ids 0 to 5 (`register_legacy_seq_types()`).
3. It infers the ordering when the model states none: V-VD-D-DJ-J if a gene choice has class D, and
   V-VJ-J otherwise. A v1 model file states no ordering. A JSON model states it in its `layout`.
4. It freezes the registry.
5. It resolves each event's `seq_type_id` from its seq type name.
6. It tells each event its neighbours in the ordering (`set_adjacent_segments()`). This is the
   only place where the ordering becomes adjacency: events are told who their neighbours are, and
   never look it up.

> **Subject to change — CAP T7b (T7 step 4).** The pinning goes. Ids will follow the ordering, and
> `register_legacy_seq_types()` and `legacy_seq_type_registry()` will be deleted. Until then, every
> model also has ids for the six legacy names, whether its ordering uses them or not. Code that
> visits the segments must walk `ordering()`, not `[0, total_count())`.

> **Subject to change — MIG L5, SYN step 5, CAP B10.** L5 moves the registry, the ordering and the
> layout, with `finalize()`'s work, into `RecombinationModel`. `Model_Parms` then remains only as
> the v1 importer. SYN step 5 makes the layout a DAG of segments with existence guards (milestone 2,
> optional segments). The ordering then becomes one path among several, and a segment's neighbour
> depends on the scenario.

> **Subject to change — ITER §9.** Point 1 above is a limit of the legacy engine, not of the model.
> It goes when the events are queried by capability rather than looked up by key.

### 2.3 Spans and boundaries

[SegmentSpan.h](../../src/igor/Core/SegmentSpan.h):
- `SegmentBoundary{id, side}` is one end of one segment.
- `SegmentSpan{left, right}` is a pair of boundaries. `SegmentSpan::gap(l, r)` is the span between
  two segments, from `l`'s 3′ end to `r`'s 5′ end. No caller builds any other shape.
- `cut_position()` turns an offset into a half-open bound: a 5′ offset stays as it is, a 3′ offset
  gains one. An offset (`Seq_Offset`) is the index of a nucleotide, so a segment is a closed
  interval. Closed intervals do not compose by subtraction, and two past off-by-one bugs came from
  that (ITER §7.1, §7.8). With half-open bounds, a span's length is `end − begin`, and two adjacent
  spans add up.
- `lies_strictly_inside(registry, id, span)` reads the ordering. D lies inside V→J. D1 and D2 both
  lie inside V→J, and D1 lies inside V→D2.

A span has no single length while some of the events that change it are still to run. It has a
**profile**, the best probability for each length it can still reach (§6.3), and the scenario has
committed to one **observed distance** (PBM §2).

> **Subject to change — ITER O11 step 3, CAP Phase D.** A span over a segment's own extent, rather
> than a gap, can be represented but has no consumer yet, so the queries do not handle it. It comes
> with the DP engine (SYN step 9).

---

## 3. Scenario state: layered containers

### 3.1 `LayeredArray`

[LayeredArray.h](../../src/igor/Core/LayeredArray.h) is an integer-keyed array in which every key has
a stack of layers. A layer is a depth of the walk. Writing key `k` at layer `n` leaves the value at
layer `n − 1` intact, so going back up the tree costs nothing.

Each key has two marks:
- **claimed** (`claimed_layer()`): the highest layer granted to an owner, by `request_layer()`
  during initialization. It records ownership.
- **current** (`current_layer()`): the layer the key's data stands at. It moves as the walk writes
  and backtracks. It is −1 for a key never written, and `exists()` is then false.

The rules, enforced at run time:
- **A write lands on a claimed layer.** `set(key, value, layer)` throws if `layer` was never
  claimed, and so does `set_current()` on a key with no claim at all. Layer 0 is no exception: it is
  the layer the first `request_layer()` grants (R3b, ITER decision O10).
- **A read needs a write.** Reading an unwritten key throws, and so does reading a layer above the
  key's current one. A layer that was claimed but not written cannot be read.
- `get()` returns by value: growing the storage reallocates it.
- `restore_layer()` releases the top claim.
- `multiply_all(acc, layers)` multiplies `acc` by one value per key, each read at the layer a
  snapshot names. The pruning bound uses it (§6.1).

The storage is flat, `storage_[key + layer * count]`. The error paths are cold functions outside the
class: building an error message inline was enough to stop GCC inlining the check (CAP D4).

### 3.2 `DynamicSequenceMap` and `SeqOffsetsMap`

[DynamicSequenceMap.h](../../src/igor/Core/DynamicSequenceMap.h) is a `LayeredArray` keyed by
`SeqTypeId` and sized from a frozen registry. The downstream bound map is created larger, with
`TraversalPlan::bound_key_count()` keys (§4.3). For sequences, the map tells three states apart:
- **not processed yet**: never written;
- **absent**: written, as an empty segment;
- **present**: written, and not empty.

`occupied()` is true for a present segment only. `first_occupied_left()` and
`first_occupied_right()` walk the ordering and skip the first two states.

> **Subject to change — CAP B10.** No event calls `first_occupied_*()` today. An insertion and a
> Markov chain read their neighbours from the ordering (§2.2), which is right only while every
> segment is present in every scenario, as milestone 1 requires. Milestone 2 decides how absence is
> written, and whether "neighbour" then means the nearest occupied segment.

[SeqOffsetsMap.h](../../src/igor/Core/SeqOffsetsMap.h) keeps the two ends of every segment in two
`DynamicSequenceMap<SeqOffset>`s, `five_prime` and `three_prime`. Code that knows its side should
address that map directly.

### 3.3 The maps of a walk

| Map | Type | Key | Written by | Read by |
|---|---|---|---|---|
| constructed sequences | `Seq_type_str_p_map` (`DynamicSequenceMap<Int_Str*>`) | segment | gene choice (template), deletion (trimmed copy), Markov chain (the junction it creates) | Markov chain (its anchor), error rate, counters |
| offsets | `Seq_offsets_map` (`SeqOffsetsMap`) | segment, end | gene choice (both ends), insertion (both ends of its junction), deletion (the end it trims) | every geometry check, insertion length, counters |
| mismatch lists | `Mismatch_vectors_map` | segment | gene choice (from the alignment), deletion (trimmed list) | error rate, counters |
| pruning mismatch floor | `Pruning_mismatch_floor_map` | segment | nothing yet (see below) | nothing yet |
| downstream bounds | `Downstream_scenario_proba_bound_map` (`DynamicSequenceMap<double>`) | segment or gap key (§4.3) | every event, for the factors it owns | `compute_upper_bound()` (§6.1) |
| overlap safety | `SafetyMatrix` | pair of segments, by row | gene choice, deletion (§5) | gene choice, deletion |
| index map | `Index_map` (`LayeredArray<size_t>`) | event identifier | `update_parent_tracking()`: a parent moves its children's marginal row | every event's probability lookup |

> **Subject to change — CAP H1–H9 (AA Pgen).** The pruning mismatch floor is allocated, passed to the
> walk in `ExplorationContext` and covered by the layer checks, but no event writes or reads it. It
> holds the mismatches that no reading of an ambiguous (IUPAC or amino-acid) query can remove, the
> count pruning must use for such a query. For a nucleotide query it would equal the mismatch lists.
> It is wired when AA Pgen lands.

### 3.4 The layer contract

Every event follows three rules. The first is enforced by `LayeredArray`; the test harness checks
the second at every hand-off ([ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md) §4.1):
1. **Claim before you write.** `initialize_event()` calls `request_layer()` for every key the event
   will write, and stores the granted layer (`memory_layer_*` members).
2. **Write what you claimed, on every path that hands off.** A reader of `layer − 1` downstream
   would otherwise find nothing. The `no_d_align` path once claimed the safety layers and never
   wrote them (ITER §7.9).
3. **Read what was there before you.** An event reads `claimed − 1` for a key it also writes. A
   Markov chain reads its anchor at the layer the anchor's last writer **before it** claimed
   (`memory_layer_anchor_seq`), not at the current layer: in the tandem model, D1's 5′ deletion
   writes D1 again after the D1D2 chain has run (`01d298a`, ITER §7.12).

Two invariants hold at every leaf. Debug builds assert them in `iterate_wrap_up()`:
- no constructed sequence holds an `int_undefined` nucleotide (`first_unfilled_segment()`);
- both ends of every segment in the ordering are placed (`first_unplaced_segment_end()`).

Since R1 and R3, a Markov chain **creates** the junction it fills, and an insertion only places it.
No segment is ever partly constructed, so the first invariant holds at every hand-off, not only at
the leaf.

> **Subject to change — MIG L3.** `ScenarioContext` gains a realization record (each event's uid,
> index and payload, written at its layer), and the maps keyed by event (`Index_map`, the offset
> map, `processed_events`) are keyed by uid.

> **Subject to change — CAP T7.** `ScenarioContext`'s wrappers still take a `Seq_type`, and events
> cast their id to it. They will take a `SeqTypeId`. Only the types change, with T7 or later.

---

## 4. The engine's plan: `SegmentIndex` and `TraversalPlan`

Which segments an event checks for overlap, and which gaps it bounds, depend on the order the events
run in. That order is the engine's. The engine therefore works these answers out once and tells
each event before `initialize_event()`. Nothing in this section names a gene or a junction (CAP T3).

### 4.1 `SegmentIndex`

[SegmentIndex.h](../../src/igor/Model/SegmentIndex.h) records, for each segment of the ordering:
- which events create or modify each end, from `get_offset_role()`;
- which events create or modify its content, from `get_seq_construction_role()`.

It does not depend on the order, so every engine can share it. `has_own_factor(segment)` is true
when the segment's content can still change against the read: an event modifies it (a deletion), or
the event that creates its content also places it (a templated gene). Content filled into a slot
that another event placed, the Markov chain's nucleotides, has no factor of its own. The gap it
fills carries its bound.

### 4.2 `TraversalPlan` and `EventPlacement`

[TraversalPlan.h](../../src/igor/Model/TraversalPlan.h) walks the engine's order once and records, for
each event, an `EventPlacement`:

| Field | Meaning |
|---|---|
| `target` | the event's segment |
| `placed_left`, `placed_right` | segments whose end facing the target was created before this event runs, nearest first. These are the overlap partners |
| `built_left`, `built_right` | the nearest segment on each side whose content was built before this event runs |
| `target_built` | whether the target's content is built once this event has run |
| `enclosing_gap()` | the gap the target sits in, between `built_left` and `built_right` |
| `gap_toward(side)` | once the target is built, the gap from it to the nearest built content on `side`; until then, the gap it sits in |

*Placed* and *built* coincide for gene choices and deletions in every shipped model. An insertion
places a slot without building it, so the gap around it stays a single gap until the Markov chain
fills the slot.

`InitialState` lists segments that exist before the first event runs. A test fixture uses it when
it writes a segment directly instead of running the event that creates it. A later Phase D cluster
will use it too, with its neighbours' boundaries given.

The plan is built once per run, in `GenModel::infer_model()`, from the queue order. It is never
modified afterwards, so all threads read it concurrently. It files placements by event identifier,
which `copy()` keeps, so one plan serves every thread's copy of the model. An event finds its own
placement with `placement()`.

### 4.3 Keys of the downstream bound map

The bound is a product over the map's keys (§6.1). Each factor needs a key:
- **a segment's own factor** is held at the segment's id;
- **a gap** is held at the id of the single segment inside it, when there is exactly one and it has
  no factor of its own (the slot an insertion places). Otherwise the gap gets a key of its own, past
  the registry's ids, in the order the walk first meets it.

`gap_key(span)` returns the key, and throws for a span that is never an open gap in this order: the
event asking and the plan then disagree about the model. `gaps()` lists the gaps, and
`bound_key_count()` sizes the map. The first rule is a choice of layout, not of meaning. It keeps
the product's order, and so the pruning of every shipped model, bitwise. In VDJ, V→D and D→J stay at
their insertions' ids and V→J gets a key of its own; in VJ, V→J stays at `VJ_ins_seq`.

### 4.4 `check_admissible()`

The plan refuses an order the events cannot run in, with `std::invalid_argument` naming the event
and the reason. It replaces the checks the events used to make on themselves:
- every event's target is a segment of the ordering;
- every end and every content is created by exactly one event;
- an end or a content is modified only after it is created;
- content is filled into a slot only after the slot's ends are created;
- a slot comes after its neighbours' facing ends, and sits inside a gap, with built content on both
  sides;
- every segment is complete when the order ends.

> **Subject to change — MIG L1, L5.** L1 makes `check_admissible()` the engine's. L5 builds the plan
> from `Topology` and the registry owned by `RecombinationModel`, instead of from the `Model_Parms`
> queue. With several engines (SYN §5), the order is each engine's own choice, and a model's
> priority becomes a hint.

> **Subject to change — CAP Phase C (MIG §7).** There is no validation when a model is loaded beyond
> `finalize()` and this check. A malformed model that the order admits can still return wrong
> numbers (CAP decision D2). The model's validation moves into `RecombinationModel` at L5.

---

## 5. Overlap safety

Segments placed on the read must not overlap. When an event places or trims a segment, the ends of
its neighbours may still move: the deletions on them have not all been chosen yet. Each check can
reach one of three verdicts: no outcome can avoid an overlap, every outcome avoids it, or it depends
on events still to come.

### 5.1 How far an end can still move

`JunctionGeometry::PendingModifierBounds`
([JunctionGeometry.h](../../src/igor/Model/JunctionGeometry.h)) is rebuilt by each gene choice and
deletion in `initialize_event()`. For every segment end, it sums `get_offset_delta_bounds()`
over the events not yet processed, and for every segment, it sums `get_length_contribution()`. The
deltas add up, so a model with more than one deletion on an end needs no new code. A range with
`min > max` is refused.

- `reachable(id, side, current)` returns the interval an end can still reach from where it is. It
  is a single point once nothing pending acts on that end.
- `check_overlap(left_3′, right_5′, gap)` returns `Infeasible` when even the best case overlaps,
  `Safe` when even the worst case does not, and `Undetermined` otherwise. `gap` is 0 at every call
  site today.

### 5.2 `SafetyMatrix`

[SafetyMatrix.h](../../src/igor/Model/SafetyMatrix.h) records which pairs of segments are already
known not to overlap in the scenario being explored.
- **Indexed by ordering position**, not by id. `cell(a, b)` returns a `SafetyCell{row, column}`,
  the 5′-most segment first. It is resolved once, in `initialize_event()`.
- **One 32-bit word per row**, layered like a `LayeredArray`. A model can therefore have at most 32
  segments in its ordering, and the constructor refuses more. VDJ has 5, tandem D 7.
- **A safe verdict propagates along the row.** Recording `(A, B)` as safe also marks `(A, C)` for
  every `C` 3′ of `B`. If a complete scenario overlapped `A` and `C`, it would also overlap `B` and
  `C`, and that pair is checked by whichever of the two is placed second. Propagation therefore
  changes where such a scenario is discarded, never whether it is.
- **A write starts from the row as the enclosing depth left it**, so cells nobody touched at this
  depth keep their value.

### 5.3 Who checks what: `FlankCheck`

`resolve_flank_checks()` ([FlankCheck.h](../../src/igor/Model/FlankCheck.h)) builds an event's checks
from its placement, during `initialize_event()`:
- **Nearest placed neighbour only**, on each side requested (ITER §2.3, decision O1). A pair further
  apart follows from the nearer ones, through the propagation above. A gene choice checks both sides
  of the segment it places. A deletion checks the side it trims.
- **One safety-layer claim per check.** The 5′ neighbour's pair is in that neighbour's row and the
  3′ neighbour's pair in the target's row, so the two checks never share a row.
- Each check also records the layer its partner's facing end was claimed at
  (`partner_offset_layer`).

In `iterate()`, an event reads the verdict the enclosing depth left
(`is_overlap_safe(cell, safety_layer − 1)`). If the pair is not already safe, it computes the
verdict from the reachable intervals, discards the realization if it is `Infeasible`, and writes
its verdict at its own layer (`set_overlap_safety()`).

The argument rests on one property: every pair of segments adjacent in the order of placement is
checked by someone. The nearest-only rule was decided in O1 but only enforced in T5 (`199ce6d`).
Before that, an event compared itself with every placed segment on a side. In the tandem joining
order, J's 5′ deletion then compared J with the D1D2 junction, a pair no gene choice checks, and
read a safety layer below the first one.

> **Subject to change — CAP B10.** With optional segments, the nearest placed neighbour depends on
> the scenario, and the propagation argument must then hold on every path of the layout.

> **Subject to change — ITER §8.1, C2.** Insertions and deletions as sequencing errors would widen
> `reachable()` by an indel budget. `reachable()` is a function of the bounds object for that
> reason. They are not scheduled.

---

## 6. Probability bounds and pruning

### 6.1 The pruning test

At every node, after an event has chosen a realization and written its state:

```
upper_bound = scenario_proba × Π_key downstream_proba_map[key]     (compute_upper_bound)
prune       ⇔ upper_bound < seq_max_prob_scenario × proba_threshold_factor
```

The product reads each key at the layer recorded in the event's snapshot,
`current_downstream_proba_memory_layers`. `Rec_Event::initialize_event()` takes that snapshot from
`claimed_layers()` after the event's own claims, so the event sees its own writes and those of the
events above it. `should_prune()` makes the test and reports it to the bound instrumentation
(§6.6). The leaf uses `is_below_threshold()`, the same comparison without the report: there the
value tested is a realized probability, not a bound.

Each key holds one factor of what the scenario can still reach:
- **a segment's own factor**: the error-rate bound for the mismatches that no pending deletion can
  remove (`surviving_core_bound()` for a gene choice, `get_err_rate_upper_bound()` for a deletion
  once the opposite end is settled);
- **a gap's factor**: the junction bound (§6.2).

An event that resolves a factor writes 1.0 at that key and multiplies what it realized into
`scenario_proba` instead. By the leaf every key is 1.0, and the bound equals the realized
probability.

### 6.2 Junction bounds: `JunctionBound` and `SpanProfile`

An event holds three `JunctionBound` slots ([SpanProfile.h](../../src/igor/Model/SpanProfile.h)):
`kLeftJunction`, `kRightJunction` and `kEnclosingJunction`. Three are enough in every layout:
a segment has two flanks and lies within at most one gap. In a tandem model, D1 bounds V→D1 and
D1→D2 and lies within V→J, exactly as D does in VDJ.

`initialize_event()` resolves each slot from the plan, once: its span, its key (`gap_key()`), its
claimed layer and its fold mode.

| Mode | Built during the sweep | Written by `iterate()` |
|---|---|---|
| `Yes` | the profile of the gap (§6.3) | the profile read at the observed distance |
| `No` | nothing | 1.0: the event splits the gap and its two halves carry the bound |
| `Retain` | the decomposition of the enclosing gap (§6.4) | 1.0, as for `No` |

`SpanProfile` maps each distance a span can reach to the best probability for it. It is a dense
`std::vector<double>` indexed by distance, at most about 300 entries on every shipped model.
`best_for(distance)` is one bounds check and one load. It returns nothing when no scenario reaches
that distance, and the caller then discards the realization. A bound of 0.0 is different: it is
pruned on probability.

Which profile an event reads depends on whether its own realization changes the span's length:
- **No** (a gene choice, from which the span is measured): one profile, `profile()`, read at the gap
  as it stands.
- **Yes** (an insertion filling the gap, a deletion moving one of its ends): one profile per
  realization, `profile_for(realization)`, read at the gap as it stood before the event's own
  choice. A deletion reads the same distance for every realization. An insertion reads the gap
  between its neighbours, which it fills without moving them.

### 6.3 The fold

`Rec_Event::initialize_Len_proba_bound()` fills every `Yes` slot of an event, the **reader**, from
the events after it in the queue.

**Participants.** An event after the reader takes part in a span when
`participates_in_span(span) = affects_length_of(span) ∨ affects_proba_of(span)`. The filter runs
once per table, not once per node.
- `affects_length_of()` is one rule for every event, read off its offset roles and the ordering
  (R12): the event changes the span's length if it creates or modifies an end of a segment strictly
  inside the span, or if it modifies one of the span's two ends. Creating one of the two ends does
  not count: that end is where the span is measured from. So a gene inside the span adds its
  template, an insertion inside it adds its length, a deletion adds its trim whether it moves an end
  of the span or an end of a segment inside it, and V and J at the ends of every junction add
  nothing. The rule replaced four tables written per event type. They agreed with it except on D's
  deletions on V→J, which the tables left out: J's V→J table was then not a bound (ITER §7.22, R16).
- `affects_proba_of()` defaults to false. A Markov chain answers true for a span its segment lies
  strictly inside: its nucleotides add a factor to the span without adding length.

**A table never contains its reader.** By the time the reader reads its table it has chosen, and its
probability is already in `scenario_proba`. Counting it again put the bound below the realized
probability and pruned scenarios above the threshold (ITER §7.19, R8). When the reader changes the
span's length, the table is **conditioned** on its realization: one fold per realization, starting
from that realization's length at weight 1 (R13).

**The fold** (`SpanFold::run()`) computes, for every total length, the best probability the
participants can reach together. It does not enumerate paths. It takes the participants one at a
time and keeps a `FoldFrontier`: the best probability for each `FoldState`, which is the part of a
partial path that a later participant can read:
- `length`: the summed `length_delta()` so far, the profile's key in the end;
- `parent_offsets`: where each later participant reads its marginals, given the realizations the
  fold has already chosen for its parents (R6);
- `lengths`: the `UnfilledSegmentLengths` published on this path (below).

Two paths in the same state are interchangeable for everything after them, so only the better one
is kept. The result equals the path-by-path walk bit for bit: every kept value is one path's
left-to-right product, and rounding never reverses the order of two non-negative values. The cost
is the number of states instead of the number of paths. On human IGH, one sweep went from 1.66 s to
0.40 s (R12/R16 stage 1).

**One participant's step** (`fold_step()`):
- An event that changes the span's length enumerates its realizations. For each one, it multiplies
  by `realization_bound()`, adds `length_delta()`, publishes into `lengths`, and moves the
  `parent_offsets` of the later participants it conditions.
- An event that only contributes a factor, the Markov chain, multiplies by `span_proba_factor()`
  once.

**Conditioning** (`SpanConditioning`, R6). An event's marginals hold one row per configuration of
its parents. `realization_bound()` takes the maximum over the rows still open: those of parents the
fold has not chosen. A parent the fold chooses itself (an earlier participant, or the reader of a
conditioned table) is read at its realization. The strides come from the parents'
`memory_and_offsets`, the same strides the walk applies, so no event type is named. This is a
joint maximum over the group (the D gene and its deletions), where the product of the
per-event maxima would be looser.

**Unfilled lengths.** A Markov chain's factor depends on how many nucleotides it will choose, which
an earlier participant decided. `UnfilledSegmentLengths`
([UnfilledSegmentLengths.h](../../src/igor/Model/UnfilledSegmentLengths.h)) carries, per segment, how
many nucleotides its offsets imply that nobody has chosen yet. The publisher is whichever event
creates a segment's offsets but not its content:
- an insertion publishes its length;
- a gene choice publishes 0, since its template fixes its nucleotides;
- a deletion publishes nothing and adds its negative delta to the span.

So a chain states its requirement as a property of the state ("offsets placed, content not yet
created"), never as a link to its insertion. The keying matters: on V→J both junctions' chains take
part at once.

**The Markov chain's factor** is the best chain of `L` steps over the fifteen nucleotide codes a
read position can hold, each pair priced as `iterate()` prices it. `prepare_span_proba_factor()`
computes it once per sweep, before any table the chain takes part in is folded (the sweep runs in
reverse queue order). This is tighter than the former `p^L` wherever the best pair is not a
self-transition (IGK among the shipped models), and remains a bound because ambiguous positions are
priced by an average (R12/R16 stage 3a).

**Once per iteration.** The bound depends on the marginals alone, so one thread folds it
(`#pragma omp single`) and the others copy the profiles (`adopt_Len_proba_bound()`, S4e). Copying
is cheaper than sharing a pointer: `best_for()` runs at every node of the walk.

### 6.4 The retained decomposition

A gene choice between two placed neighbours, which scans positions when the aligner finds nothing
(`no_d_align`, a D gene in every shipped model), needs more than a profile. For each total length,
it needs every way its own realization and the two halves can add up to it. Its enclosing slot is
then `Retain`, and `build_retained_decomposition()` builds a `SpanDecomposition`: placements
`(realization, left distance, right distance, bound)` per total, from the folded left and right
profiles and its own realizations. Each bucket is **sorted by decreasing bound**. The scan stops at
the first pruned placement, which is exact only in that order (ITER §2.6, §6.19).

> **Subject to change — PBM §6 "Planned".** The exhaustive scan visits D's placements in the order
> of their bound, so a tighter table reorders it and moves the `no_d_align` outputs. A sort that
> does not depend on the bound would make every later tightening safe for the order. **Not decided.**

### 6.5 Soundness and what is still loose

- **The bound is a bound.** Since R8 and R16, the instrument finds no node whose bound sits below the
  best leaf under it, on the uniform start or on an inferred model, and an unpruned evaluation agrees
  with the pruned one (PBM §9).
- **It is still a relaxation.** Parents the fold does not choose (above the reader, or the reader of
  a table that is not conditioned on it) are maximized over independently. On the shipped models
  that covers every conditioning edge. The J gene choice carries about 7 of the 15 decades of slack
  at the root (PBM §9).
- **It does not depend on the read**, which is why it is built once per iteration.
- **Tables are built per reader.** Two events reading the same span each fold it (ITER §2.5,
  deferred).

> **Subject to change — MIG L4b, L6.** The bound table is rebuilt after the handlers' M-step (L4b)
> and becomes tier-0b work, shared and const (L6).

> **Subject to change — CAP Q14, SYN step 9.** On the tandem TRB model the fold takes 24.5 s against
> 0.055 s for TRB, and an iteration about 250 times as long (CAP *T5*). Whether real tandem data
> runs on this engine or waits for the `ChainDPEngine` is decision Q14 of MIG, still open.

> **Subject to change — ITER S4d (optional).** A max-convolution of two profiles and Tensor-backed
> containers. Neither has a consumer.

### 6.6 Measuring the bound: `BoundTightness`

[BoundTightness.h](../../src/igor/Model/BoundTightness.h) measures how far the bound sits above what
the walk finds. Its entry points are empty unless `IGOR_BOUND_INSTRUMENTATION` is defined, so the
default build is unchanged.

```
pixi run build_instrumented    # build_instr/, -DENABLE_BOUND_INSTRUMENTATION=ON
build_instr/bin/igor …         # the report goes to stderr at exit
```

It reports, per depth: the bound over the best leaf below each node, the nodes with no leaf below
them and why (`starved`: no child was tested; `pruned`: every child was pruned; `hollow`: the
cause is further down), and the over-estimate split by the event that resolves it. PBM §9 explains
how to read the report, and its caveat: on a uniform model, the numbers describe the shape of the
tree, not the looseness of a converged model.

---

## 7. The leaf, the error model and the readers that refuse

At a leaf, `iterate_wrap_up()` computes the error-weighted probability, tests it against the
threshold, updates the best scenario, calls each counter with a `Scenario` view, and calls
`add_to_marginals()` on every event that is not fixed.

- **`Single_error_rate`** sums the nucleotides and mismatches of every segment of the ordering that
  has a mismatch list, by id: the genes, not the junctions (T7a). The counts are integers, so the
  order of the sum does not matter. Before T7a it read ids 0, 2 and 4 only, and charged no mismatch
  to a tandem model's D segments.
- **`check_admissible(const SeqTypeRegistry &)`** exists on `Error_rate` and on `Counter`, and
  accepts any model by default. The two hypermutation error models and four counters (errors,
  coverage, best scenarios, Pgen) still read V, D and J at their pinned ids. They refuse an ordering
  that is neither V-D-J nor V-J through `migration::require_vdj_or_vj_ordering()`
  ([vdj_ordering_migr.h](../../src/igor/Model/Legacy/vdj_ordering_migr.h)). `GenModel` calls the
  checks before any thread starts, so a refused run writes no file.

> **Subject to change — CAP T7b, CNT CT2a.** The `Scenario` snapshot is sized from the registry and
> indexed by id (CT2a for the counters, T7b for the hypermutation models). Then the guard leaves
> them and `vdj_ordering_migr.{h,cpp}` is deleted.

> **Subject to change — ERR E0–E4, H0–H5.** `Single_error_rate` moves to a core shared by the legacy
> class and a new-architecture error node (Part I). The hypermutation models follow (Part II), and
> SYN step 8 makes the error model a Probe event with handlers. ER1, in E1, moves every inference
> output once.

---

## 8. Where each check lives

| Check | When | Where | Failure |
|---|---|---|---|
| nickname present and unique | event added | `Model_Parms::add_event()`, the readers | exception |
| JSON document and each event node against the schema | load | `ModelSchema`, each JSON event constructor | exception listing violations by JSON pointer |
| one event per (type, segment, side) | `finalize()` | `Model_Parms::finalize()` | exception |
| the order is coherent | before inference | `TraversalPlan::check_admissible()` | `std::invalid_argument` |
| the error model and counters can read the model | before inference | `check_admissible()` (§7) | `std::invalid_argument` |
| the event's own target and set-up | `initialize_event()` | each event: seq type id against the registry, a side for a deletion, a Markov chain for an insertion, a junction for a chain | exception |
| `min ≤ max` for offset deltas and length contributions | `initialize_event()` | `PendingModifierBounds::rebuild()` | `std::logic_error` |
| writes land on claimed layers; reads find a write | run time | `LayeredArray` | `std::out_of_range` |
| at most 32 segments in the ordering | iteration set-up | `SafetyMatrix` | `std::length_error` |
| nothing unfilled, every end placed | each leaf | `iterate_wrap_up()` | assertion, debug builds only |
| claimed layers are written before hand-off; the A0 declarations match the writes | each hand-off | the unit harness, `call_iterate_recording()` | test failure |

---

## 9. Further reading

| Topic | Where |
|---|---|
| why each piece is shaped as it is, step by step | ITER §2 (G1–G10), §6 (delivery notes) |
| the defects found and fixed on the way | ITER §7, phase R in ITER §6.9 |
| the history, costs and measurements of the junction bound | PBM |
| writing tests for an event | [ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md) |
| what each `iterate()` did before the rewrite | [ITERATE_METHOD_ANALYSIS.md](../ITERATE_METHOD_ANALYSIS.md) |
| what comes next, in order | [IMPLEMENTATION_GRAPH.md](../IMPLEMENTATION_GRAPH.md) |
