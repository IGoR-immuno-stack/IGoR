# What a recombination event implements and exposes

**Status**: reference documentation of the code on `feature/tandemD-T5` @ `3c06f52` (Oct 23 2026).
**Audience**: developers who add an event type to IGoR, or change one of the four that exist
(`Gene_choice`, `Deletion`, `Insertion`, `Dinucl_markov`).
**Companion**: [ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) explains the machinery that calls the
event: the registry, the layered scenario state, the engine's plan, overlap safety and the
probability bounds. Read its §1 first.
**Scope**: the contract of `Rec_Event` ([Rec_Event.h](../src/igor/Model/Legacy/Rec_Event.h)) as the
legacy inference engine (`GenModel::infer_model()`) and the legacy generator (`SequenceGenerator`)
use it. The Model module's `Topology`, handlers and engines are merged but not wired, and are not
covered.

The plans are cited by short names (CAP, ITER, PBM, MIG, ERR, CNT, GEN, SYN), listed at the top of
[ITERATE_INTERNALS.md](ITERATE_INTERNALS.md). Parts that a planned step changes carry a note:

> **Subject to change — MIG L3.** What the step changes.

§6 collects these notes in one table.

---

## 0. The contract in one table

| Group | Members | Provided by | Called by |
|---|---|---|---|
| Identity and facts set from outside | nickname, event identifier, seq type and its id, neighbours, traversal plan, priority, side, realizations | the model and the engine (§1) | everything |
| Declarations | `get_offset_role()`, `get_seq_construction_role()`, `get_offset_delta_bounds()`, `get_length_contribution()`, `length_delta()`, `affects_proba_of()`, `span_proba_factor()`, `prepare_span_proba_factor()` | each event type (§2) | `SegmentIndex`, `TraversalPlan`, `PendingModifierBounds`, the fold |
| Set-up | constructors, `copy()`, `initialize_event()` | each event type (§3.1, §3.2) | the readers, `GenModel` |
| Once per iteration | `initialize_Len_proba_bound()`, `adopt_Len_proba_bound()` (base, not virtual); `update_event_internal_probas()`; `set_crude_upper_bound_proba()` | the base, plus the declarations | `GenModel` (§3.3) |
| Exploration | `iterate()` | each event type (§3.4) | the previous event's `iterate_wrap_up()` |
| Leaf | `add_to_marginals()`, `ind_normalize()` | each event type (§3.5) | `iterate_wrap_up()`, the M-step |
| Generation | `draw_realization()`, `construct_realization()`, `propagate_realization()` | the base, plus each event type (§3.6) | `draw_random_realization()` |
| Serialization | JSON constructor, `to_json()`, `write2txt()`, `write2txt_legacy()` | each event type (§3.1) | the readers and writers |

---

## 1. What the model and the engine tell an event

An event never works these out for itself. It is told, once, and reads the answer.

| Fact | Set by | Read through | Notes |
|---|---|---|---|
| **nickname** | the model file; `set_nickname()` | `get_name()`, `get_nickname()` | The event's identity, required and unique in its model (T2). It keys `Model_Parms::edges`, the marginals' index and offset maps and `processed_events`, and names the event in model files and in outputs |
| **event identifier** | `Model_Parms`; `set_event_identifier()` | `get_event_identifier()`, the protected `event_index` | Indexes `Index_map` and the next-event array. `copy()` keeps it, which is how the shared plan finds every thread's copy |
| **seq type** | the constructor | `get_seq_type()` | The name of the segment the event acts on, as the model gives it. The only copy since T1 |
| **seq type id** | `Model_Parms::finalize()` | `get_seq_type_id()` | The run-time handle (ITERATE_INTERNALS §2). `kNoSeqType` until resolved |
| **neighbours** | `Model_Parms::finalize()` | `get_left_adjacent_id()`, `get_right_adjacent_id()` | The segments 5′ and 3′ of the event's own in the ordering. What the event does with them is its own business: an insertion spans both, a Markov chain seeds from the side it runs from |
| **traversal plan** | the engine, before `initialize_event()`; `set_traversal_plan()` | the protected `traversal_plan()`, `placement()` | What is placed and built when the event runs, and the key of each gap (ITERATE_INTERNALS §4). Not carried by `copy()`: the engine sets it on each thread's copy |
| **priority, side** | the model file | `get_priority()`, `get_side()` | A deletion's side is the end it trims; a Markov chain's side is the end of its anchor the chain starts from |
| **realizations** | the model file | `event_realizations`, `realization_at()` | An `Event_realization` holds a name, an `int` value, a string value and its integer form, and an **index**, dense from 0, which is its coordinate in the marginals |
| **gene class** | the model file | `get_class()` | A gene choice's library: which alignments it reads. A deletion's is V, D or J for the three legacy gene segments and `Undefined_gene` otherwise |
| **parent strides** | `Rec_Event::initialize_event()`, from the offset map | `memory_and_offsets` | For each child: its identifier, the index-map layer claimed for it, and the stride by which this event's realization moves the child's marginal row |

> **Subject to change — MIG L3.** The nickname stops keying the legacy maps: they are keyed by an
> integer uid assigned in nickname order. Files and outputs keep the nickname.

> **Subject to change — CAP T2 "Not in T2", T4b.** `event_class` moves off the base into
> `Gene_choice`. It waits for `FastGenerator`, which still finds the segment a deletion trims from
> the deletion's class.

> **Subject to change — CAP T7b.** Until the legacy id pinning goes, ids 0 to 5 belong to the six
> legacy names in every model.

---

## 2. What an event declares

Declarations are answered from the event's own instance data, at initialization, never in the
walk. They let the generic machinery ask "what does this event do to that segment" without naming
the event type. Every consumer reads them through `SegmentIndex`, `TraversalPlan`,
`PendingModifierBounds` or the fold.

### 2.1 Roles and ranges (A0)

All four are pure virtual.

| Query | Answers | Values |
|---|---|---|
| `get_offset_role(id, side)` | what the event does to one end of segment `id` | `None`; `Creates` (sets the initial offset); `Modifies` (shifts an existing one) |
| `get_seq_construction_role(id)` | what it does to the content of segment `id` | `None`; `Creates` (allocates it); `Modifies` (truncates or extends it); `Fills` (fills placeholders in it; no event answers this since R1) |
| `get_offset_delta_bounds(id, side)` | how far the event can still **shift** an end, signed | `OffsetDelta{min, max}`. A 3′ end moves left as nucleotides are removed, a 5′ end right. `{0, 0}` if the event cannot move it |
| `get_length_contribution(id)` | how many nucleotides the event adds to or removes from segment `id`, over all its realizations | `LengthContribution{min, max}`. Positive for content supplied, negative for content removed |

Rules:
- **`min ≤ max`**, always. A provider that flips a sign must reorder the pair.
  `PendingModifierBounds::rebuild()` throws otherwise.
- **The roles describe what `iterate()` writes.** An event that declares `Creates` on an end must
  write that offset, on every path that hands off. The unit harness checks this at every hand-off
  (ITER §2.5, tier 3). `Insertion` used to answer `None` while placing its junction, which hid the
  fact that it recorded no offsets (fixed by R3).
- **Each end and each content of a segment has exactly one creator.** Modifiers run after it.
  `TraversalPlan::check_admissible()` checks both.
- **A segment is complete when the model's order ends**: both ends placed, content built. The two
  can come from different events: an insertion places a junction, and the Markov chain after it
  creates the content (R1, R3, ITER decision O12).
- `SeqConstructionRole::Creates` and `OffsetRole::Creates` coincide for every event today but are
  separate statements. Do not merge them.

### 2.2 Per-realization and span hooks

| Hook | Kind | Answers |
|---|---|---|
| `length_delta(realization)` | pure virtual | The signed length one realization adds to a span it affects: `+` template length for a gene, `+` the value for an insertion, `−` the value for a deletion. Called only on an event whose `affects_length_of()` is true |
| `affects_proba_of(span, registry)` | virtual, default false | Whether the event adds a probability factor to the span without adding length |
| `span_proba_factor(span, unfilled)` | virtual, default 1 | That factor, given the `UnfilledSegmentLengths` published on the fold's path |
| `prepare_span_proba_factor(marginals, index_map)` | virtual, default nothing | Precomputes what `span_proba_factor()` reads, from this iteration's marginals. Runs before any table the event takes part in is folded |

The base class derives the rest. These are not virtual:
- `affects_length_of(span, registry)`: one rule for every event, read off its offset roles and the
  ordering. A realization changes the span's length when the event creates or modifies an end of a
  segment strictly inside the span, or modifies one of the span's two ends (R12).
- `participates_in_span()` is `affects_length_of() or affects_proba_of()`, the filter the fold
  applies.
- `creates_own_offsets()` decides whether the event publishes an unfilled length in the fold. An
  event that creates its segment's offsets but not its content publishes the length it decided; one
  that creates both publishes 0; one that creates neither publishes nothing.

### 2.3 The four events' answers

On their own segment (`id == seq_type_id`); every query answers `None` or `{0, 0}` on any other.

| | `Gene_choice` | `Deletion` | `Insertion` | `Dinucl_markov` |
|---|---|---|---|---|
| `get_offset_role` | `Creates`, both ends | `Modifies`, on its side | `Creates`, both ends | `None` |
| `get_seq_construction_role` | `Creates` | `Modifies` | `None` | `Creates` |
| `get_offset_delta_bounds` | `{0, 0}` | its side: 3′ `{−max, −min}`, 5′ `{min, max}` of its deletion values | `{0, 0}` | `{0, 0}` |
| `get_length_contribution` | `{shortest, longest}` template | `{−max, −min}` | `{fewest, most}` insertions | `{0, 0}` |
| `length_delta(r)` | template length | `−value` | `+value` | 0, never called |
| `affects_proba_of` | false | false | false | its segment lies strictly inside the span |
| `span_proba_factor` | 1 | 1 | 1 | best chain of `L` steps (`L` = the unfilled length of its segment) |
| publishes in the fold | 0 | nothing | its length | nothing |

A negative deletion value is a palindromic extension; the ranges above stay ordered for it.

> **Subject to change — CAP "Out of this plan", SYN §3.** The rest of Phase A lands with Phases C and
> D and with B10, reshaped by the synthesis. What an event **reads** becomes
> `context_reads() → {toward, anchor, width, codon_frame}`, relative to the nearest built content,
> in place of `get_context_dependency()` and `get_context_seq_types()`. `is_branching` becomes an
> engine's `candidates()`, and `is_multi_realization` the arity of the event's domain. Roles,
> deltas and reads together form the event's **footprint**, which belongs to the model and serves
> every engine. No bound-specific capability is planned: the plan derives what bounds need from the
> roles, the order and the R12 rule.

---

## 3. The lifecycle: what each hook must do

### 3.1 Construction, copy and serialization

- **From a JSON node.** A JSON model builds every event through `event_factory::create(nickname,
  node, library)` ([EventFactory.h](../src/igor/Model/EventFactory.h)). An event type provides:
  - a constructor `(const std::string &nickname, const nlohmann::json &node)`, or
    `(nickname, node, library)` if it draws templates from the document's library, which it declares
    by specializing `event_factory::takes_library<T>`;
  - one registration line in its translation unit: `static Registrar<T> r{"TypeName"};`
  - a schema, `src/igor/Model/schemas/0.3/event-<type>.schema.json`, listed in the event union. The
    constructor checks its node against it first (`JsonDetail.h`), so an event built in code is
    checked too.
- **From a v1 model file.** The legacy constructors, through `Model_Parms`' reader. The v1 format
  translates gene classes to seq types and back (`legacy_gene_class_of()`,
  `legacy_gene_class_to_seq_type()`, [LegacyFormat.h](../src/igor/Model/Legacy/LegacyFormat.h)). A
  model the v1 format cannot express (`requires_extended_format()`) is written as JSON.
- **`copy()`** returns an independent event for another thread. It must carry the realizations,
  the seq type and its id, the nickname, the priority, the `fixed` flag and the event identifier.
  It must not carry the traversal plan, which the engine sets on each copy.
- **`to_json()`** is implemented once on the base: type, priority, target, realizations sorted by
  index. Override it only for state that is not a realization, as `Gene_choice` does with its
  library. **`write2txt_legacy()`** writes the v1 event line; **`name_class_token()`** gives the gene
  class token of the v1 name.

> **Subject to change — SYN step 4b, REVIEW §14–15.** The model document gains its `provenance` and
> `accumulators` sections, several targets per event, and layouts other than `Linear`; schema 0.x
> versions replace each other until 1.0.

### 3.2 `initialize_event()`

Called once per thread per EM iteration, in queue order, after the engine has set the plan. On
entry, `processed_events` holds the events before this one in the queue; the event is added to it
at the end. An override does, in this order:

1. **Check its own set-up**, and throw with its nickname if it cannot run. For example: its seq type
   id matches the registry entry for its name (`Insertion`); it names a side (`Deletion`); a Markov
   chain fills its junction (`Insertion::require_dinucl_markov()`); its chain has a neighbour to seed
   from (`Dinucl_markov`). Coherence across events is the plan's, already checked.
2. **Read its placement** (`placement()`), never search `events_map` for neighbours or partners.
3. **Claim every layer it will write**, and store the layer granted:
   - each offset it creates or modifies (`seq_offsets.request_layer(id, side)`);
   - the content it creates or modifies (`constructed_sequences`), and the mismatch list if it has
     one;
   - its own downstream-bound slot, at its id, if it has a factor of its own;
   - one slot per gap it bounds, at `traversal_plan().gap_key(span)`;
   - its safety rows, through `resolve_flank_checks()`, which also records the partners' offset
     layers.
4. **Resolve its junction slots**: `junction_bound(slot).resolve(span, key, layer, fold)`, with the
   gaps from `placement().gap_toward()` or `enclosing_gap()` and a fold mode (`Yes`, `No`,
   `Retain`; ITERATE_INTERNALS §6.2).
5. **Settle the geometry** it needs: `PendingModifierBounds::rebuild()` over the events not yet
   processed, and anchoring read off the ordering (a segment with no neighbour on a side sits at an
   end of the read). Never test a gene class or a segment name.
6. **Precompute what depends on the model only**: a lookup of realizations by value, by index, the
   largest junction a chain can fill.
7. **Call `Rec_Event::initialize_event()` last.** It reserves the current-realization slot, claims
   the index-map layers of its children (`memory_and_offsets`), takes the snapshot of the downstream
   map's claimed layers that `compute_upper_bound()` reads, and adds the event to
   `processed_events`. Called earlier, the snapshot would miss the event's own claims.

> **Subject to change — MIG L6.** This becomes tier-0a work, run once per run on each copy rather
> than once per iteration.

### 3.3 Once per iteration: bounds and internal probabilities

- **`initialize_Len_proba_bound()`** and **`adopt_Len_proba_bound()`** are the base's and not
  virtual. They fold every `Yes` junction slot from the events after this one, or copy another
  thread's result. An event takes part through its declarations alone (§2.2).
- **`update_event_internal_probas()`** is the hook for state derived from the marginals.
  `Dinucl_markov` builds its dinucleotide matrix there. It is also the only set-up hook the legacy
  generator calls.
- **`set_crude_upper_bound_proba()`** stores the event's largest marginal.
  `Dinucl_markov::span_proba_factor()` falls back on it, as `p^L`, beyond the chains it has
  prepared.

> **Subject to change — MIG L4c.** The marginals come from the handlers' tensors, and
> `update_event_internal_probas()` reads them there. The layout and the arithmetic stay the same.

### 3.4 `iterate()`

One node of the walk. Every body has the same shape:

```
base_index = exploration.index_map.get(event_index)
inherited  = scenario.scenario_proba
for each realization the read admits:
    proba = inherited × iterate_common(r.index, base_index, …)   // also moves the children's rows
    write the event's state at its claimed layers                 // offsets, content, mismatches
    for each flank check:                                         // gene choice, deletion
        if not is_overlap_safe(cell, safety_layer − 1):
            verdict = check(reachable intervals)
            Infeasible → next realization
            set_overlap_safety(cell, verdict == Safe, safety_layer)
    set the downstream slots it owns:                             // own factor, junction slots
        a resolved factor → 1.0, its value goes into proba
        a junction → profile.best_for(observed distance); none → next realization
    scenario.scenario_proba = proba
    upper = exploration.compute_upper_bound(proba, current_downstream_proba_memory_layers)
    if not exploration.should_prune(upper):
        iterate_wrap_up(…)                                        // next event, or the leaf
scenario.scenario_proba = inherited
```

What every body must respect:
- **Restart each realization from the inherited probability.** Realizations do not compound.
  `Gene_choice`'s exhaustive position scan once did, and moved the `no_d_align` outputs when fixed
  (ITER §7.16, R7).
- **Write every claimed layer on every path that hands off**, and read the layer below your claim
  for keys you also write (ITERATE_INTERNALS §3.4).
- **Record the realization** in `current_realizations_index_vec`: `add_to_marginals()` and the
  counters read it at the leaf.
- **Use only what `initialize_event()` resolved.** No name lookup, no registry lookup and no span
  question in the walk: slots, keys, layers, cells and partners are handles.
- **Discard what no scenario can complete.** A missing junction bound means no completion reaches
  the observed distance: skip the realization.

> **Subject to change — ITER phase R.** `Insertion::iterate()` throws on a junction length its table
> does not reach, where the other events skip. Changing it is a behaviour change, still open.

> **Subject to change — MIG L2, L8.** L2 replaces `iterate_wrap_up()` with an engine continuation.
> L8 moves the `iterate()` bodies into `Inference/Dfs/` as `LegacyDfsEngine`'s kernels. `Rec_Event`
> keeps the definition of the event (MIG Q15).

### 3.5 At the leaf

- **`add_to_marginals(proba, marginals)`** adds the scenario's error-weighted probability at the
  marginal index of the realization chosen on this path (or sets it, for a Viterbi-like run).
  `iterate_wrap_up()` calls it on every event not fixed.
- **`ind_normalize(marginals, base_index)`** normalizes one row of the marginals, when
  `Model_marginals` normalizes them after the M-step or after reading a file. `Dinucl_markov`
  overrides it: its row holds 4 × 4 transitions, normalized per previous nucleotide.

> **Subject to change — MIG L4a, L4b.** Accumulators become tensors in the legacy layout, and the
> M-step becomes the handlers' `maximize()`.

### 3.6 Generation

`draw_random_realization()` is not virtual. It calls, in order:
- **`draw_realization(marginals, index_map, state, rng)`** draws realization indices and writes
  nothing. The base draws one index from the event's marginal row (`pick_realization()`).
  `Dinucl_markov` overrides it to draw its chain, seeded by its anchor and as long as its junction.
- **`construct_realization(indices, state)`**, pure virtual, builds what the indices stand for into
  the `GenerationState`, keyed by `SeqTypeId`. It uses no RNG, marginals or index map: it is the
  event's *apply*. The sequence is assembled in the order of `registry.ordering()`.
- **`propagate_realization(indices, index_map, offset_map)`** moves each child's row by the drawn
  realization times the child's stride.

A construction refuses what inference never counts: a deletion past its segment's end, a palindrome
longer than its template, a gene at an end of the sequence deleted away (GEN D2, D3).

> **Subject to change — MIG L7, SYN step 1.** A sampling driver and `SamplingEngine` drive
> `construct_realization()` with indices they sampled, and `draw_random_realization()` retires. The
> categorical walk then goes by index rather than in `unordered_map` order, which moves the
> `generate` outputs once. The draw at a zero-mass boundary (GEN D4) is left to the sampling
> handlers.

---

## 4. What an event must not do

- **Search the model** for its neighbours, partners or gaps. The model tells it its neighbours, and
  the engine's plan tells it what is placed and built.
- **Name a gene, a junction or a legacy segment** (`V_gene_seq`, `D_gene`, …). Behaviour that
  depends on position is read off the ordering.
- **Write a layer it did not claim**, or hand off with a claimed layer unwritten.
- **Count its own probability in a bound it reads**, or merge across its own realizations in a
  table it reads (R8, R13). The base's fold already handles both.
- **Hold scenario state across reads** outside the layered maps and the members that
  `initialize_event()` set up.
- **Depend on another event's type.** A Markov chain needs "a segment whose offsets are placed and
  whose content is not created yet", a property of the state, not "my insertion".

---

## 5. Checklist: adding an event type, or a model shape

**A new event type.**
1. The class, deriving from `Rec_Event`: the four A0 roles, `length_delta()`, and the span hooks if
   it adds a factor without length (§2).
2. Its JSON constructor, `Registrar` line and schema (§3.1). A v1 event line only if the v1 format
   can express it.
3. `copy()`, `initialize_event()` as in §3.2, `iterate()` as in §3.4, `add_to_marginals()`.
4. Generation: `construct_realization()`, and `draw_realization()` if its draw is not one
   categorical index.
5. Tests, following [ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md):
   - the ten rows of its §4 matrix, through `call_iterate_recording()`, which checks the layer
     contract and the A0 declarations at every hand-off;
   - its declarations (`test_event_capabilities.cpp`);
   - the plan's admissibility on an order that uses it (`tst/igor/Model/test_traversal_plan.cpp`);
   - a JSON round trip (`test_model_formats.cpp`) and its construction in generation
     (`test_generation_construction.cpp`).
6. The gate: unit and integration suites, the five regression tracks (cli, align, inference,
   no_d_align, generate) bitwise unless the commit says which outputs move and why, and the
   convergence tests.

**A new model shape** with existing events, as tandem D was:
- write it as a JSON document, whose `layout` gives the segment ordering;
- check that `TraversalPlan::check_admissible()` accepts its queue order, and that the ordering has
  at most 32 segments;
- the legacy engine holds one event per (type, segment, side);
- the hypermutation error models and the errors, coverage, best-scenarios and Pgen counters refuse
  any ordering other than V-D-J or V-J, until CNT CT2a and CAP T7b;
- order the events so that every segment a chain or a deletion reads is still there when it runs.
  The tandem TRB model joins its segments D1 to D2 first, then D1D2 to J, then V, each joint's
  deletions right before its insertion (CAP *T5*). If a fully deleted anchor can reach its chain,
  the chain throws, in generation and in inference alike, and inside inference's threads that
  terminates the process, whatever probability the model gives the case (ITER §7.12, GEN D6).

---

## 6. Where the contract is going

| Today | Becomes | Step |
|---|---|---|
| the nickname keys the legacy maps | an integer uid, in nickname order | MIG L3 |
| `iterate_wrap_up()` and the next-event array | an engine continuation; the leaf is the engine's | MIG L2 |
| `GenModel` runs the set-up and the read loop | an EM driver and `LegacyDfsEngine` | MIG L1 |
| the plan is built from `Model_Parms`' queue | built from `Topology` and the registry, owned by `RecombinationModel` | MIG L5 |
| set-up once per iteration per thread | tiers 0a (per run) and 0b (per iteration, shared) | MIG L6 |
| marginals in `Model_marginals`' flat array | handler tensors, same layout | MIG L4 |
| `draw_random_realization()` | `SamplingEngine` driving `construct_realization()` | MIG L7, SYN step 1 |
| the `iterate()` bodies on the event classes | `LegacyDfsEngine`'s kernels in `Inference/Dfs/` | MIG L8 |
| `Error_rate` and the hypermutation models | an error node, then a Probe event with handlers | ERR E0–E4, H0–H5; SYN step 8 |
| `Counter` | `Observer`, reading typed records through sinks | CNT, MIG L8 |
| `event_class` on the base | `Gene_choice`'s library | CAP T2 note, after T4b |
| `ScenarioContext` methods take a `Seq_type` | they take a `SeqTypeId` | CAP T7 |
| ids 0–5 pinned to the legacy names | ids follow the ordering | CAP T7b |
| context dependencies, `is_branching`, `is_multi_realization` | `context_reads()`, the engine's `candidates()`, the domain's arity | CAP Phases C/D, B10 |
| every segment present in every scenario | optional segments in a layout DAG, adjacency per scenario | CAP B10, SYN step 5 |
| one event per (type, segment, side) | events queried by capability | ITER §9 |
| one engine, the DFS | several, among them `ChainDPEngine` | SYN steps 7 and 9; MIG Q14 |
