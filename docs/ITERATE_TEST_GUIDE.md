# Writing `iterate()` tests

**Audience**: anyone adding tests for a `Rec_Event` subclass or an `Error_rate`, and anyone
changing existing ones while refactoring.
**Status**: matches the harness on `feature/tandemD-T5` @ `3c06f52` (Oct 23 2026).
**Reference implementation**: [test_gene_choice_iterate.cpp](../tst/igor/Model/Legacy/test_gene_choice_iterate.cpp).
`Gene_choice` is fully covered, and every pattern below has a worked instance there.
[test_deletion_iterate.cpp](../tst/igor/Model/Legacy/test_deletion_iterate.cpp) shows how sections
stay organised by pattern when one body serves four kinds of deletion.
[test_tandem_junctions.cpp](../tst/igor/Model/Legacy/test_tandem_junctions.cpp) shows a fixture on
a tandem-D layout.
**Harness**: [test_utils.h](../tst/igor/Model/Legacy/test_utils.h) /
[.cpp](../tst/igor/Model/Legacy/test_utils.cpp), in `tst/igor/Model/Legacy/`.
**Related**: [ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) explains the machinery a test exercises
(layers, the plan, safety, bounds), and [REC_EVENT_CONTRACT.md](REC_EVENT_CONTRACT.md) what an event
must do. [ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) records how the
bodies were rewritten and the defects found on the way.

---

## 1. What you are testing

`iterate()` is one node of a depth-first walk over recombination scenarios. Whatever the
event, the body has the same shape:

```
read base index  ->  enumerate realizations  ->  per realization:
      compute probability
      write this event's segment state
      check feasibility against neighbours
      set downstream probability bounds
      prune or hand off to the next event
```

The unit under test is therefore not "the event" but **the hand-off**: for each realization that
survives, what does the next event see? Everything worth asserting is visible at that point.

---

## 2. The harness in one page

```cpp
auto state = create_iterate_state(read);              // owns all five contexts
auto event = make_gene_choice(V_gene, {{"V1", seq}}, /*id=*/0, /*fixed=*/false);
state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, seq.size())});
state.set_marginal(0, 1.0L);

auto rec = call_iterate_recording(event, state);      // runs iterate(), records hand-offs
REQUIRE(rec->call_count() == 1);
CHECK(rec->calls.at(0).five_prime(V_gene_seq) == 0);
```

`IterateTestState` owns one `Storage` struct per context and exposes the contexts over it. It is
neither copyable nor movable: the contexts hold references into its own storage.

`create_iterate_state(read, marginal_size, max_events, registry)` sizes every map from
`registry`. The default is `vdj_seq_type_registry()`; `vj_seq_type_registry()` gives a model with no
D, and a test on another layout builds and freezes its own registry (§6).

| Setter | Effect |
|---|---|
| `set_alignments(gene_class, …)` | the alignments a `Gene_choice` of that library will enumerate |
| `set_marginal(i, p)` | one entry of the model's marginal array |
| `set_base_index(event_id, n)` | where an event reads its marginals from (default 0) |
| `condition_on(child, parent)` | conditions the child's marginals on the parent's realization, as an edge in a model file does |
| `set_scenario_proba(p)` | the probability inherited from upstream (default 1) |
| `set_error_rate(r)` | replaces `Single_error_rate(0.0)` |
| `set_pruning_threshold(t)` | prunes below `t` (default 0, i.e. off) |
| `add_event(ev)` | registers the event in `events_map` only |
| `mark_chosen(ev)` | the event has already run: it is in `processed_events` (so its travel is no longer pending), and it comes before the event under test in the plan's order |
| `preset_segment(id, 5', 3', seq, mismatches)` | writes a segment as its upstream events would have left it, at layer 0, and lists it as built in the plan's `InitialState` |
| `preset_junction(id, 5', length)` | writes a junction as an `Insertion` leaves it: two offsets and **no** sequence. Listed as placed, not built |
| `preset_safety(left, right, value)` | writes an overlap verdict at layer 0. A `true` verdict propagates along the row |
| `add_downstream_event(ev)` | registers the event **and** queues it after the event under test |

`call_iterate()` reproduces `GenModel`'s set-up, in its order:
1. every event in `events_map` gets its neighbours from the registry (standing in for
   `Model_Parms::finalize()`), a base index, a marginal size and a crude bound;
2. the queue: the event under test, then the downstream events. The chosen events are **not**
   initialized (trap 2);
3. a `TraversalPlan` over the order the code under test sees: the chosen events, the event under
   test, the downstream events, then any other event in `events_map`, with the presets as its
   `InitialState`. Every event is given the plan; `state.plan()` returns it afterwards, for
   `gap_key()` lookups in assertions;
4. `initialize_event()` over the queue;
5. the next-event chain: the event under test hands off to the recorder, or to nothing;
6. `initialize_Len_proba_bound()` **in reverse queue order**. This pass is load-bearing (trap 4);
7. the event's `iterate()`.

### `ScenarioSnapshot`

One per surviving realization. Only seq types that were actually written appear, so a missing key
is itself an observation. The maps are keyed by `Seq_type`; a segment of another layout is reached
by casting its id.

```cpp
double scenario_proba;
std::map<Seq_type, std::pair<Seq_Offset, Seq_Offset>> offsets;   // .five_prime(t) / .three_prime(t)
std::map<Seq_type, std::string> sequences;   // ACGT; an ambiguity code as 'N', int_undefined as '.'
std::map<Seq_type, std::vector<std::size_t>> mismatches;
std::map<std::pair<Seq_type, Seq_type>, bool> safety;             // .at({V_gene_seq, D_gene_seq})
std::map<SeqTypeId, double> downstream_bounds;                    // every key, gap keys included
```

`downstream_bounds` holds every key of the downstream map, including a gap's own key past the
registry's ids (ITERATE_INTERNALS §4.3). Look a gap's key up with
`state.plan().gap_key(SegmentSpan::gap(left, right))` rather than assuming it is an insertion's id.

`safety` is the one map where "missing" means something coarser. The overlap verdicts of one left
segment share a word (`SafetyMatrix`), so a pair appears as soon as *any* pair with the same left
segment has been written, with the verdict the enclosing depth left. Absent means nothing touched
that left segment at all.

---

## 3. Choosing an observation instrument

| Instrument | Use for | Cost |
|---|---|---|
| **`RecordingEvent`** (`call_iterate_recording`) | almost everything | one snapshot per surviving realization |
| **Leaf path** (`call_iterate` with no recorder) | `add_to_marginals`, counters | needs a *complete* scenario; see §7 |
| **Direct call** on `compute_scenario_error_probability()` | `Error_rate` | presets only; see §7 |
| **Custom `Counter`** | scenario-level views that mirror production statistics | more setup than the recorder, rarely worth it |

**Default to the recorder.** The leaf path runs the whole leaf: the error rate, the threshold,
`add_to_marginals()` on every event that is not fixed, and, in a Debug build, the two leaf
invariants (every segment of the ordering filled and placed). A single-event fixture rarely
satisfies all of that, so a test through the leaf says little about the event and may abort.

---

## 4. The common test matrix

**Every event needs all ten rows.** Rows 1 to 9 follow the steps of §1; row 10, and the two
checks of §4.1, come with the harness. A missing row is a gap, not a judgement call.

| # | Row | What to assert | `Gene_choice` instance |
|---|---|---|---|
| 1 | **Baseline write** | one clean realization produces the offsets, sequence and mismatch list the next event reads | `baseline writes (G4)` |
| 2 | **Realization enumeration** | N realizations → N hand-offs, in order, each carrying *its own* state (not the last one written) | `realization branching` |
| 3 | **Probability** | `incoming × model_parameters[base_index + realization_index]`; realizations do **not** compound; incoming multiplies through; `base_index` is honoured | `realization branching` |
| 4 | **Feasibility guards** | one section per `continue` in the body, each with a **positive control** differing minimally | `overlap verdicts`, `junction-length bound` |
| 5 | **Neighbour dependence** | behaviour when a neighbour is placed / in the model but not placed yet / absent from the model | `overlap verdicts`, last three sections |
| 6 | **Downstream bounds** | which keys get a bound and what it is | `endogenous-mismatch`, `junction-length bound` |
| 7 | **Pruning** | below threshold → skipped; positive control above it | `pruning` |
| 8 | **Empty vs absent** | a written-but-empty segment has `exists() == true` with empty contents; never-written is absent | `mismatch propagation`, last section |
| 9 | **Memory layering** | the event reads the layer below its claim and writes its own; state written at one layer is still readable at the other | *(see note)* |
| 10 | **Layer contract** | every map the event requests a layer in is *written* at that layer on every path that hands off | **automatic**, see §4.1 |

Notes on the two rows that need care:

**Row 3: "do not compound" is the one a rewrite loses.** Each realization must restart from the
probability the event *inherited*, not from the previous realization's result. Two realizations at
0.5 and 0.25 give 0.5 and 0.25, never 0.5 and 0.125. Assert it explicitly: a body that does
`scenario.scenario_proba *= contribution` on the shared field passes every other row. The
`no_d_align` scan did exactly that until R7 (ITER §7.16).

**Row 9: layering.** `Gene_choice` writes at the layers `initialize_event()` claimed and reads the
alignment, so it has little to assert. `Deletion` and `Dinucl_markov` read what was written below
their claim, and a body that reads the wrong layer sees stale or absent state. Assert that a
segment preset at layer 0 is what the event reads, and that its own write lands where the next
event reads it. A Markov chain reads its anchor at the layer the anchor's last writer before it
claimed, not at the current layer (ITER §7.12, `01d298a`): pin that with a writer of the anchor
initialized after the chain.

### 4.1 Rows checked for you

`call_iterate_recording()` runs three checks at every hand-off and fails the section with a
message naming the map, the key and the layers. **You get them by using the harness; there is
nothing to write.**

**The layer contract (row 10): a claimed layer must be written.** An event requests a layer so it
can write without clobbering the previous value, and so that downstream readers of `layer − 1` see
that value. Request-without-write leaves the next reader on unwritten storage, which `LayeredArray`
refuses. That was ITER §7.9: the `no_d_align` path requested the safety layers and never wrote
them. The check, over `constructed_sequences`, both ends of `seq_offsets`, `mismatches_lists`,
`downstream_proba_map`, `safety_set` and `pruning_mismatch_floor`:

```
for every key this event claimed a layer for:
    current_layer(key) == claimed_layer(key)   at hand-off
```

**Ownership: a written layer must have been claimed.** The converse, over the same maps. Since R3b,
`LayeredArray::set()` also throws on it at run time, so the harness check now mostly localises the
failure. No event is exempt any more.

**Declarations: the event did what its A0 roles say, and nothing else.** Before the event runs,
the recorder reads its `get_seq_construction_role()` and `get_offset_role()`. At each hand-off it
checks that this event wrote every content and every end it declares it `Creates` or `Modifies`
(a `Fills` must also leave no undetermined nucleotide), and that it left untouched every one it
declares `None` (ITER §2.5, tier 3). The check is for consistency, not correctness: an event that
declared `None` everywhere and wrote nothing would pass it. `TraversalPlan::check_admissible()` is
what requires every segment to have a creator.

What the three have in common: **they only see paths a section actually reaches.** They are
dynamic checks on a structural property, so a missing write on a branch no test executes goes
unnoticed. That is what row 4 and coverage are for. `index_map` is excluded: its layering is driven
by parent tracking through the offset map, which single-event fixtures rarely populate.

When one trips, start at the event's `initialize_event()` and find the `request_layer()` whose
write is missing on that path, or the write whose `request_layer()` is.

### Rows that only apply to some events

| Row | Applies when | Note |
|---|---|---|
| **Non-branching** | `Insertion`, `Dinucl_markov` | exactly one hand-off, or zero if the value is out of range. Assert the zero case too |
| **Reverse traversal** | `Dinucl_markov` whose anchor sits 3′ of its junction (`event_side == Five_prime`) | the read window and the filled segment are both handled back to front; assert the *final* orientation |
| **Two-sided state** | `Deletion` on a segment both of whose ends move | whether the opposite end can still move (`opposite_end_still_moves_`) changes the endogenous-mismatch bound |
| **Position-dependent behaviour** | `Gene_choice`, `Deletion` | a segment at an end of the ordering behaves differently from an internal one (overhang clipping, the exhaustive scan, keep-one, the early prune). Cover both positions, and a tandem layout where it matters |

---

## 5. Per-event additions

### `Gene_choice`
- Both realization paths: the alignment path, and the exhaustive position scan used when no
  alignment survives (`no_d_align`), in its two sub-branches. With a placed neighbour on both sides
  the placements come from the retained decomposition, best first, and the scan stops at the first
  prune; with fewer, the window slides against the read.
- The scan runs only for a segment internal to the ordering (`exhaustive_position_fallback_`, ITER
  decision O6): a D in every shipped model, and both Ds of a tandem one. Assert that a segment at an
  end of the ordering never reaches it.
- The template overhang: a segment at the 5′ end of the ordering may overhang the read before it,
  one at the 3′ end after it. B3 will replace the clipping with flank segments.
- Alignments of a template the event does not hold are skipped (`c275c42`): two gene choices on
  one library may hold different templates.

### `Deletion`
Covered in [test_deletion_iterate.cpp](../tst/igor/Model/Legacy/test_deletion_iterate.cpp)
(ITER step 4a, kept through 4b). Read it before adding to it; it covers:

- positive deletion: sequence truncated, offset moved, mismatch list trimmed to a **contiguous
  subrange** of the incoming list (the *prefix* for a 3′ deletion, the *suffix* for a 5′ one);
- negative deletion (palindrome): reverse complement appended (3′) or prepended (5′), new
  mismatches computed against the read, and the 5′ cases re-sort;
- the full-deletion guard: a segment at an end of the ordering keeps at least one nucleotide; an
  internal one may be deleted entirely. Assert the asymmetry: it is a modelling decision, not a
  refactoring one;
- a zero-length segment: written, `exists()` true, offsets in the degenerate convention
  `three_prime == five_prime - 1`;
- the two-stage prune of a segment at an end of the ordering: the first check `break`s (deletions
  are enumerated in decreasing order), the second `continue`s. **The two are not separable by
  observation**: the first bound dominates the second, and both are monotone along the enumeration.
  Pin the stopping, not the mechanism, and say so in the section;
- the flank check: against the nearest placed neighbour **on the side it trims** only.

### `Insertion`
- the insertion count is derived from the neighbours' facing offsets, not from a stored value;
- a count outside the realization range: **zero hand-offs**, no state written;
- it writes **both offsets** of its junction and **no sequence**: the `Dinucl_markov` after it
  creates the segment (R1, R3). Its declarations say so (`OffsetRole::Creates`,
  `SeqConstructionRole::None`), and §4.1 checks them;
- it reads `profile_for(its realization)` at the gap between its neighbours, and **throws** when no
  completion reaches it, where the other events skip the realization. Pin it; changing it is an open
  behaviour change;
- it performs **no** overlap check (ITER §7.7). Pin that too, so that adding one is deliberate.

### `Dinucl_markov`
- it **creates** the junction segment, sized from the offsets the `Insertion` placed: preset those
  with `preset_junction()`, not with a placeholder segment;
- the seed nucleotide comes from the anchor, and which end of it from `event_side`;
- forward and reverse fills produce the same final orientation;
- an empty anchor: the chain **throws** when it has nucleotides to choose (R2), and an empty
  junction never reads its anchor. Both are covered; the skip-to-the-next-occupied-segment
  alternative waits for B10.

---

## 6. Setup recipes

**A neighbour that already ran** (placed, with offsets):

```cpp
auto d_stub = make_gene_choice(D_gene, {{"D1", d_gene}}, /*id=*/1);   // fixed by default
state.add_event(d_stub);
state.mark_chosen(d_stub);
state.preset_segment(D_gene_seq, /*5'=*/12, /*3'=*/15, d_gene);
```

Do **not** put such a neighbour through `initialize_event()` (trap 2).

**Non-zero deletion bounds on a neighbour** (what the feasibility checks read):

```cpp
state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, /*min*/0, /*max*/4, /*id*/3));
```

Registering it downstream also makes it a participant in the junction tables, which is usually
what you want (trap 4).

**A junction an insertion has sized**, for a `Dinucl_markov` under test:

```cpp
state.preset_junction(VD_ins_seq, /*5'=*/8, /*length=*/3);
```

**A layout other than VDJ or VJ.** Build a registry with the legacy names pinned, set the
ordering, freeze it, and move each event onto its segment:

```cpp
SeqTypeRegistry r;
r.register_legacy_seq_types();
r.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq", "D2_gene_seq",
                     "DJ_ins_seq", "J_gene_seq"});
r.freeze();
auto state = create_iterate_state(read, 1000, 32, r);

auto d1 = make_gene_choice(D_gene, {{"D1", seq}}, /*id=*/2);
d1->set_seq_type("D1_gene_seq");
d1->set_seq_type_id(r.id("D1_gene_seq"));
d1->set_nickname(default_nickname(*d1));    // two events cannot share a nickname
```

`test_tandem_junctions.cpp` wraps the last three lines in a helper, `on(event, "D1_gene_seq")`.

**Pruning:**

```cpp
state.set_pruning_threshold(0.5);   // drop anything whose upper bound is below 0.5
```

**A meaningful error bound:** `set_error_rate(0.1)`. The default rate of 0 makes every bound with
a non-zero mismatch count collapse to 0, which asserts nothing.

---

## 7. Testing `Error_rate`

`compute_scenario_error_probability(query, model, scenario, exploration)` runs at the leaf, so the
recorder is the wrong instrument. Call it directly on a state whose segments are preset, as
[test_single_error_rate.cpp](../tst/igor/Model/Legacy/test_single_error_rate.cpp) does.

`Single_error_rate` charges every segment of the ordering that **has a mismatch list**: the genes,
whose gene choices and deletions write one, and not the junctions (T7a). Presence is signalled by
the mismatch list, not by the sequence. `preset_segment()` writes both; `preset_junction()` writes
neither, which is what a junction looks like.

```cpp
state.preset_segment(V_gene_seq, 0, 11, v_seq, {2});
state.preset_segment(D_gene_seq, 13, 16, d_seq, {});
state.preset_segment(J_gene_seq, 18, 25, j_seq, {20});
```

What to assert:

| Row | Assertion |
|---|---|
| **Weighted probability** | `scenario_proba × (r/3)^n_err × (1−r)^(n_genomic − n_err)` |
| **Genomic count** | the sum of the *constructed* gene segment lengths, not the read length |
| **Which segments** | every gene of the ordering, a tandem model's D1 and D2 included; junctions uncharged; an absent D changes both counts, a present-but-empty one contributes 0 |
| **Threshold** | below `seq_max_prob × factor` returns 0 and updates **no** accumulator |
| **Accumulators** | `seq_likelihood`, `seq_weighted_er`, `seq_probability` accumulate across scenarios, and are *assigned* rather than accumulated under `viterbi_run` |
| **Admissibility** | `check_admissible(registry)` accepts the orderings the model can read and refuses the others with a message naming the reader ([test_tandem_guards.cpp](../tst/igor/Model/Legacy/test_tandem_guards.cpp)) |

> **Subject to change — ERR E0, H0, H1.** E0 extends the `Single_error_rate` suite (merge order,
> `copy`, the M-step, the round trip), and `set_error_rate()` gains a form taking any `Error_rate`.
> The hypermutation models get their unit suite (H1) and a regression track (H0). They still read
> V, D and J at the pinned ids until T7b.

---

## 8. Rules that keep tests honest

1. **Derive every expectation by running the code**, never from a comment, a document, or what
   the implementation looks like it should do. That is where the number comes from. What to do when
   the number is *wrong* is a separate question, and the answer is not "pin it": a **confirmed**
   defect is asserted at its intended value under `[!shouldfail]` (§11), and only behaviour whose
   intent is genuinely undecided is pinned as observed, with a comment saying what is undecided.
2. **Every test must discriminate.** Break the thing under test and confirm the section fails. If
   it still passes, the section is testing something else. This is not hypothetical: it happened
   twice in T0 (ITER §7.6).
3. **Every "nothing happened" assertion needs a positive control.** `call_count() == 0` is
   satisfied by a misconfigured fixture. Pair it with a section differing in exactly one value that
   yields `call_count() == 1`.
4. **Assert per hand-off, not in aggregate.** `total_marginal_mass() > 0` passes when the wrong
   realization survives. Check each snapshot's own sequence, offsets and probability.
5. **One knob per section.** If two sections differ in three ways, neither localises a failure.
6. **Name the pattern, not the gene class.** See §10.

---

## 9. Traps in this harness

1. **A preset segment's sequence must be as long as its offsets say.** A segment occupies
   `[off(5'), off(3')]` *inclusively*, so its length is `three_prime - five_prime + 1`. A literal
   string next to a parameterised offset drifts the moment a section passes a different offset,
   leaving the scenario in a state no upstream event could have produced: harmless while the event
   under test ignores the content, and a wrong answer as soon as one does not. Use
   `segment_run(five_prime, three_prime)`, or derive the offsets from the string you want
   (`five + seq.size() - 1`), so the two cannot disagree.
2. **Do not `initialize_event()` a neighbour that already ran.** It requests a layer per map it
   writes, moving that neighbour's claim to layer 1 while `preset_segment()` writes at 0, and the
   event under test then reads an unwritten layer and throws. Register it and mark it chosen; that
   reproduces everything the code under test queries.
3. **Stub events must be `fixed`.** At a leaf, `iterate_wrap_up()` calls `add_to_marginals()` on
   every non-fixed event in `events_map`, and a stub that never iterated writes out of bounds. The
   builders default to `fixed = true`; only the event under test is unfixed.
4. **Without downstream events, every junction table is a single entry.** A table folds the events
   *after* its reader, never the reader itself (R8). With none queued, a gene choice's table is
   `{0: 1.0}`, and a deletion's or an insertion's table for realization `r` holds only its own
   length, `length_delta(r)` (R13). Every geometry whose gap no queued event can absorb is then
   discarded by the junction guard, and a test aimed at another branch quietly becomes a test of
   that guard. Queue the deletions and insertions that fill the gap with `add_downstream_event()`.
5. **`proba_threshold_factor` cannot be changed after construction**: `ExplorationContext` holds it
   by value. The harness pins it at 1 and moves `seq_max_prob` instead, which is held by reference.
   That is what `set_pruning_threshold()` does.
6. **An `OffsetDelta` is signed by the end it moves.** A deletion of `k` moves a 3′ end by `−k` and
   a 5′ end by `+k`, so a 3′ deletion over `[0, 4]` reports `{−4, 0}`, ordered. Getting this
   backwards produces fixtures that look right and test nothing. (The `*_min_del` / `*_max_del`
   scalars that used to carry this went with R10; the bounds come from `PendingModifierBounds`.)
7. **`set_mismatches()` stores a pointer** into `query.gene_alignments`. Anything that outlives or
   aliases that vector is a real bug; assert per-realization lists to catch it.
8. **A conditioned table is read at the gap as it stood before the reader chose.** A deletion's
   table for realization `r` is keyed by everything the fold enumerated, the deletion's own `−r`
   included, and `iterate()` looks it up at the gap *before* the deletion. So realization `r`
   survives the junction guard only if some downstream insertion length `i` gives
   `gap_before = i − r`. Two consequences for a fixture. A boundary realization the downstream
   events cannot absorb is discarded by the junction guard, not by the comparison the section
   aims at, which is how a `<=` versus `<` boundary ends up with no test at all: give the
   downstream insertion a range that covers the gaps the section needs, and check by mutation that
   the comparison is actually reached. And a section aimed at a guard *inside* the body may need no
   junction resolved at all: leave the neighbour on the trimmed side out of the plan (neither
   preset nor chosen), so that the deletion bounds no gap.
9. **The declaration check covers the six legacy ids only.** `call_iterate_recording()` reads the
   event's roles on `V_gene_seq` … `VJ_ins_seq`. On a tandem layout, a segment registered past
   them (`D1_gene_seq`, id 6) is not checked, in either direction: assert its writes in the
   section itself. *Subject to change with CAP T7b, which makes the tests look ids up by name.*

---

## 10. Naming and organisation

**One `TEST_CASE` per generic pattern, one `SECTION` per instance.** Not one per gene class. V/D/J
is the axis the rewrite removed; a file organised on it hides gaps: you cannot see from a
`V_gene` / `D_gene` / `J_gene` layout that the `Infeasible` verdict is covered for V and J but not
for a D's 5′ end.

```
TEST_CASE("Overlap verdicts (G2/G3)")
    SECTION("Infeasible — V 3' vs D 5'")
    SECTION("Safe — min deletions clear")
    SECTION("Undetermined — inside the deletion range")
```

Tags: the event (`[gene_choice]`, `[deletion]`, `[insertion]`, `[dinucl]`), the instrument
(`[iterate]`, `[error_rate]`), and `[tandem_d]` for a case on a tandem layout. A case that runs a
whole inference carries `[integration]`; a slow one that CI must not run starts with `[.]`.

---

## 11. Confirmed defects: `[!shouldfail]`

Do **not** pin a known-wrong value as the expectation under test. Write the assertion the code
*should* satisfy and tag the case `[!shouldfail]`:

```cpp
TEST_CASE("The categorical walk never draws a zero-mass realization", "[generation][!shouldfail]")
```

Catch2 reports an expected failure as a pass, so `ctest` stays green. When the defect is fixed the
case starts passing, which `[!shouldfail]` reports as a **failure**, forcing the tag to be removed
deliberately: no fix lands unnoticed. A test stating the intended behaviour also reads as
documentation, which a pinned wrong number does not. Phase R removed every such tag from the
`iterate()` suites, one per commit, each commit naming the outputs it moved. The two left are
generation's D4 (GEN §4), above.

**Keep these cases section-free.** `[!shouldfail]` is evaluated per test-case *run*, and Catch2
re-runs a case once per leaf section, so a case mixing passing and failing sections reports the
passing ones as unexpected passes.

Reference the plan section in the name, so the fix has an obvious landing site.

**"Confirmed" is a gate.** This applies to a defect that has been diagnosed *and* agreed not to be
the intended behaviour. Behaviour that is merely surprising gets pinned as observed, with a comment
naming what is undecided: an unestablished suspicion tagged `[!shouldfail]` is a guess that reads as
a decision.

**`[!mayfail]` is for a statistical test** whose outcome is not deterministic, such as the
convergence case. It is not a place for a defect.

**A wrong value may still appear as a premise.** `test_event_capabilities.cpp` pins
`get_len_max() == INT16_MIN` on a `Dinucl_markov` to contrast the legacy accessor with the
capability query asserted beside it. That is the defect serving as the case's setup, not as its
expectation, and it is fine.

---

## 12. Checklists

### Adding a new event's tests

- [ ] All ten rows of §4 have at least one section (row 10 and the two other checks of §4.1 are
      automatic)
- [ ] Every `continue` / `break` in the body has a section **and** a positive control
- [ ] Each section mutation-verified: break the code, confirm *this* section fails
- [ ] Confirmed defects assert the **intended** value under `[!shouldfail]`, never pinned wrong
- [ ] Sections live under pattern-named `TEST_CASE`s, not gene-class ones
- [ ] A case on a layout where position matters: an internal segment and one at an end, and a
      tandem layout if the event can sit inside one
- [ ] `pixi run test_unit` green, expected failures reported as expected; and
      `pixi run test_unit_debug`, since the default build defines `NDEBUG` and only Debug runs the
      leaf invariants and the matrix bounds assertions
- [ ] Coverage reported for the event's `iterate()` and `initialize_event()`:
      `pixi run coverage -- -f '[your_event][iterate]'`. Read the uncovered *lines*, not the
      ratio: gcov counts exception edges as branches, so 100 % is not the target. Every uncovered
      branch is either closed with a section or recorded as a deliberate gap

### Changing tests during a refactor

- [ ] A section that now fails is a **behaviour change**: establish it is intended before editing
      the expectation
- [ ] A `[!shouldfail]` case that starts passing means a defect was fixed: remove the tag, and say
      so in the commit
- [ ] Rows do not disappear. If a generic body absorbs three branches, the three sections become
      three instances of one pattern, not one section
- [ ] Update the plan's delivered-coverage table
- [ ] `pixi run build`, then `pixi run test_regression` for anything claiming to be
      behaviour-preserving: the five tracks (cli, align, inference, no_d_align, generate) must stay
      bitwise. Unit tests alone do not establish it
