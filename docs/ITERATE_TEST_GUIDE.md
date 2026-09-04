# Writing `iterate()` tests

**Audience**: anyone adding tests for a `Rec_Event` subclass or an `Error_rate`, and anyone
changing existing ones while refactoring.
**Reference implementation**: [tst/igor/Core/test_gene_choice_iterate.cpp](../tst/igor/Core/test_gene_choice_iterate.cpp)
— `Gene_choice` is fully covered and every pattern below has a worked instance there.
**Harness**: [tst/igor/Core/test_utils.h](../tst/igor/Core/test_utils.h) / [.cpp](../tst/igor/Core/test_utils.cpp)
**Related**: [ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) (what the
rewrite is), [ITERATE_METHOD_ANALYSIS.md](ITERATE_METHOD_ANALYSIS.md) (what each `iterate()`
does today).

---

## 1. What you are testing

`iterate()` is one node of a depth-first walk over recombination scenarios. Whatever the
event, the body is the same shape:

```
read base index  ->  enumerate realizations  ->  per realization:
      compute probability
      write this event's segment state
      check feasibility against neighbours
      set downstream probability bounds
      prune or hand off to the next event
```

The unit under test is therefore not "the event" but **the hand-off**: for each realization
that survives, what does the next event see? Everything worth asserting is visible at that
point.

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

`IterateTestState` owns a `Storage` struct per context and exposes the contexts over it. It
is neither copyable nor movable — the contexts hold references into its own storage.

| Setter | Effect |
|---|---|
| `set_alignments(gene_class, …)` | realizations a `Gene_choice` will enumerate |
| `set_marginal(i, p)` | one entry of the model marginal array |
| `set_base_index(event_id, n)` | where an event reads its marginals from (default 0) |
| `set_scenario_proba(p)` | the probability inherited from upstream (default 1) |
| `set_error_rate(r)` | replaces `Single_error_rate(0.0)` |
| `set_pruning_threshold(t)` | prune below `t` (default 0, i.e. off) |
| `add_event(ev)` | register in `events_map` only |
| `mark_chosen(ev)` | put in `processed_events` — makes `*_chosen` true |
| `preset_segment(seq_type, 5', 3', seq, mismatches)` | write a segment as an upstream event would have left it |
| `add_downstream_event(ev)` | register **and** place after the event under test in the queue |

`call_iterate()` reproduces `GenModel`'s ordering: index map and marginal sizes, crude
bounds, `initialize_event()` over the queue, the next-event chain, then
`initialize_Len_proba_bound()` **in reverse queue order**. That last pass is load-bearing —
see trap 3.

### `ScenarioSnapshot`

One per surviving realization. Only seq_types that were actually written appear, so a
missing key is itself an observation.

```cpp
double scenario_proba;
std::map<Seq_type, std::pair<Seq_Offset, Seq_Offset>> offsets;   // .five_prime(t) / .three_prime(t)
std::map<Seq_type, std::string> sequences;                        // decoded to ACGT, 'N' for placeholders
std::map<Seq_type, std::vector<std::size_t>> mismatches;
std::map<Event_safety, bool> safety;
std::map<Seq_type, double> downstream_bounds;
```

---

## 3. Choosing an observation instrument

| Instrument | Use for | Cost |
|---|---|---|
| **`RecordingEvent`** (`call_iterate_recording`) | almost everything | one snapshot per surviving realization |
| **Leaf path** (`call_iterate` with no recorder) | `Error_rate`, `add_to_marginals`, counters | needs a *complete* scenario; see §7 |
| **Custom `Counter`** | scenario-level views that mirror production statistics | more setup than the recorder, rarely worth it |

**Default to the recorder.** Going through the leaf pulls in `Error_rate`, which
concatenates V, D and J and throws outright if any is missing — so a single-event test
through the leaf says nothing about the event and may not even run.

---

## 4. The common test matrix

**Every event needs all nine rows.** They correspond one-to-one with the steps in §1, and
they are the rows that must stay filled as the generic rewrite lands. A missing row is a
gap, not a judgement call.

| # | Row | What to assert | `Gene_choice` instance |
|---|---|---|---|
| 1 | **Baseline write** | one clean realization produces the offsets, sequence and mismatch list the next event reads | `baseline writes (G4)` |
| 2 | **Realization enumeration** | N realizations → N hand-offs, in order, each carrying *its own* state (not the last one written) | `realization branching` |
| 3 | **Probability** | `incoming × model_parameters[base_index + realization_index]`; realizations do **not** compound; incoming multiplies through; `base_index` is honoured | `realization branching` |
| 4 | **Feasibility guards** | one section per `continue` in the body, each with a **positive control** differing minimally | `overlap verdicts`, `junction-length bound` |
| 5 | **Neighbour dependence** | behaviour when a neighbour is chosen / exists but unchosen / absent from the model | `overlap verdicts`, last three sections |
| 6 | **Downstream bounds** | which seq_types get a bound and what it is | `endogenous-mismatch`, `junction-length bound` |
| 7 | **Pruning** | below threshold → skipped; positive control above it | `pruning` |
| 8 | **Empty vs absent** | a written-but-empty segment has `exists() == true` with empty contents; never-written is absent | `mismatch propagation`, last section |
| 9 | **Memory layering** | the event reads the layer below and writes its own; state written at one layer is still readable at the other | *(see note)* |
| 10 | **Layer contract** | every map the event requests a layer in is *written* at that layer on every path that hands off | **automatic** — see §4.1 |

Notes on the two rows that need care:

**Row 3 — "do not compound" is the one a rewrite loses.** Each realization must restart from
the probability the event *inherited*, not from the previous realization's result. Two
realizations at 0.5 and 0.25 give 0.5 and 0.25, never 0.5 and 0.125. Assert it explicitly;
a body that does `scenario.scenario_proba *= contribution` on the shared field passes every
other row.

**Row 9 — layering.** `Gene_choice` writes at the layer `initialize_event()` requested and
reads the alignment, so it has little to assert. `Deletion`, `Insertion` and `Dinucl_markov`
all read `memory_layer_X - 1` and write `memory_layer_X`, and a rewrite that reads the wrong
layer sees stale or absent state. For those, assert that a segment preset at layer 0 is what
the event reads, and that its own write lands where the next event reads it.

### 4.1 — Row 10 is checked for you

**Requesting a memory layer is a promise to write it before handing off.** An event requests a
layer so it can write without clobbering the previous value, and so downstream readers of
`layer - 1` see that previous value. Request-without-write leaves the next reader on unwritten
storage: `LayeredArray` refuses it and aborts the run, and the pre-B8 containers served it as
uninitialized memory. That is exactly the defect in plan §7.9 — `Gene_choice`'s `no_d_align` path
requested the safety layers and never wrote them — and it survived for as long as it did because
nothing was looking.

`call_iterate_recording()` now checks it on every hand-off, for `constructed_sequences`, both
ends of `seq_offsets`, `mismatches_lists`, `downstream_proba_map`, `safety_set` and
`pruning_mismatch_floor`. **You get it by using the harness; there is nothing to write.**

Two things to know about it, because they decide whether it can see anything:

- **It needs downstream events.** `request_layer()` *sets* a key's current layer to the layer
  requested, so writing at that same layer changes nothing observable on its own. What makes a
  write visible is that events initialized after this one request further layers, pushing the
  current layer above this event's — and only a write pulls it back down. A fixture with no
  `add_downstream_event()` makes the check vacuous. Register the deletions and insertions a real
  model would carry; the junction-length map needs them anyway (trap 3).
- **It only covers paths a section actually reaches.** It is a dynamic check on a structural
  defect, so it finds a missing write only on a branch some test executes. That is what row 4 and
  coverage are for.

`index_map` is deliberately excluded: its layering is driven by parent-realization tracking
through `offset_map`, which single-event fixtures do not populate.

When a section trips it, the message names the map, the key, the layer requested and the layer
actually current — start at the event's `initialize_event()` and find the `request_layer()` whose
matching write is missing on that path.

### Rows that only apply to some events

| Row | Applies when | Note |
|---|---|---|
| **Non-branching** | `Insertion`, `Dinucl_markov` | exactly one hand-off, or zero if the value is out of range. Assert the zero case too. |
| **Reverse traversal** | `Dinucl_markov` with `anchor_side == Five_prime` | the filled sequence and the read window are both reversed; assert the *final* orientation |
| **Two-sided state** | `Deletion` on D | 5' and 3' interact; `d_del_opposite_side_processed` changes the endogenous-mismatch branch |

---

## 5. Per-event additions

### `Gene_choice`
Both realization-enumeration paths: the alignment path, and the `no_d_align` exhaustive
fallback (position-map and sliding sub-branches). The fallback lives inside `case D_gene`,
so V and J never reach it — assert that contrast, because B11 turns it into an explicit
switch. Also the template-overhang clipping (negative offset for V, past-the-read-end for
J), which B3 will replace with flank segments.

### `Deletion`
- positive deletion: sequence truncated, offset moved, mismatch list trimmed to a
  **contiguous subrange** of the incoming list
- negative deletion (palindrome): reverse-complement appended (3') or prepended (5'), new
  mismatches computed against the read, and the 5' case re-sorts
- the full-deletion guards: V forbids it (`>`), D and J allow it (`>=`) — assert the
  asymmetry, it is a modelling decision not a refactoring one
- a zero-length segment: written, `exists()` true, offsets in the degenerate convention
  `three_prime == five_prime - 1`
- the two-stage prune: the first check `break`s (deletions are enumerated in decreasing
  order), the second `continue`s

### `Insertion`
- insertion count derived from the neighbours' facing offsets, not from a stored value
- count outside the realization range → **zero hand-offs**, no state written
- the placeholder sequence: length `n`, all positions unset
- currently performs **no** overlap check; pin that, so adding one is a deliberate change

### `Dinucl_markov`
- the seed nucleotide comes from the neighbour, and which end depends on `anchor_side`
- forward vs reverse fill produce the same final orientation
- an empty anchor segment: today this is undefined behaviour (`.back()` on an empty
  string). Write the test that reaches it as part of B7, not before — it will crash until
  the skip-empty walk lands.

---

## 6. Setup recipes

**A chosen neighbour** (makes `*_chosen` true and gives it offsets):

```cpp
auto d_stub = make_gene_choice(D_gene, {{"D1", d_gene}}, /*id=*/1);   // fixed by default
state.add_event(d_stub);
state.mark_chosen(d_stub);
state.preset_segment(D_gene_seq, /*5'=*/12, /*3'=*/15, d_gene);
```

Do **not** put a chosen neighbour through `initialize_event()` — see trap 1.

**Non-zero deletion bounds on a neighbour** (what the feasibility checks read):

```cpp
state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, /*min*/0, /*max*/4, /*id*/3));
```

Registering it downstream also lets it contribute to the junction-length map, which is
usually what you want — see trap 3.

**Pruning:**

```cpp
state.set_pruning_threshold(0.5);   // drop anything whose upper bound is below 0.5
```

**A meaningful error bound:** `set_error_rate(0.1)`. The default rate of 0 makes every bound
with a non-zero mismatch count collapse to 0, which asserts nothing.

---

## 7. Testing `Error_rate`

`compute_scenario_error_probability(query, model, scenario, exploration)` runs at leaf
nodes, so the recorder is the wrong instrument: use `call_iterate()` **without** a recorder,
or call the method directly on a hand-built `ScenarioContext`.

It needs a **complete** scenario. `Single_error_rate` reads `V_gene_seq` and `J_gene_seq`
unconditionally and `D_gene_seq` only when `mismatches_lists.exists(D_gene_seq)` — so V and
J must always be present, and D's *presence is signalled by its mismatch list*, not by its
sequence. Preset all of them:

```cpp
state.preset_segment(V_gene_seq, 0, 11, v_seq, {2});
state.preset_segment(D_gene_seq, 13, 16, d_seq, {});
state.preset_segment(J_gene_seq, 18, 25, j_seq, {20});
```

What to assert:

| Row | Assertion |
|---|---|
| **Weighted probability** | `scenario_proba × (r/3)^n_err × (1−r)^(n_genomic − n_err)` |
| **Genomic count** | sums the *constructed* segment lengths, not the read length |
| **D optional** | absent D changes both counts; present-but-empty D contributes 0 to both |
| **Threshold** | below `seq_max_prob × factor` returns 0 and updates **no** accumulator |
| **Accumulators** | `seq_likelihood`, `seq_weighted_er`, `seq_probability` accumulate across scenarios, and are *assigned* rather than accumulated under `viterbi_run` |

> `Single_error_rate` hardcodes V/D/J (parent plan, hazard H7). Tests written against that
> shape will need revisiting when it moves to registry-ordered traversal. Prefer asserting
> the *formula* over the enumeration where you can.

---

## 8. Rules that keep tests honest

1. **Characterize, do not idealize.** Assert what the code does, derived by running it — not
   what it ought to do. Several branches are known wrong; pinning the intent there would
   make the suite red from day one.
2. **Every test must discriminate.** Break the thing under test and confirm the section
   fails. If it still passes, the section is testing something else — this is not
   hypothetical, it happened twice in T0 (plan §7.6).
3. **Every "nothing happened" assertion needs a positive control.** `call_count() == 0` is
   satisfied by a misconfigured fixture. Pair it with a section differing in exactly one
   value that yields `call_count() == 1`.
4. **Assert per hand-off, not in aggregate.** `total_marginal_mass() > 0` passes when the
   wrong realization survives. Check each snapshot's own sequence, offsets and probability.
5. **One knob per section.** If two sections differ in three ways, neither localises a
   failure.
6. **Name the pattern, not the gene class.** See §10.

---

## 9. Traps in this harness

1. **Do not `initialize_event()` a chosen neighbour.** It requests a memory layer per
   seq_type it touches, moving that neighbour's current layer to 1 while `preset_segment()`
   writes at 0 — the event under test then reads an unwritten layer and throws. Register and
   mark it; that reproduces everything the code actually queries.
2. **Stub events must be `fixed`.** At a leaf, `iterate_wrap_up()` calls `add_to_marginals()`
   on every non-fixed event in `events_map`, and a stub that never iterated still has
   `new_index == -1`, so it writes out of bounds. `make_gene_choice()` defaults to
   `fixed = true`; only the event under test is unfixed.
3. **Without downstream events the junction-length map collapses to `{0: 1.0}`.** Every
   geometry whose neighbours are not exactly adjacent is then discarded by the junction
   guard, and a test aimed at some other branch quietly becomes a test of that guard. Add
   the deletions and insertions that feed the map via `add_downstream_event()`.
4. **`proba_threshold_factor` cannot be changed after construction** — `ExplorationContext`
   holds it by value. The harness pins it at 1 and moves `seq_max_prob` instead, which is
   held by reference. That is what `set_pruning_threshold()` does.
5. **`*_max_del` and `*_min_del` are negated deletion counts.** They are `Deletion::len_min`
   and `len_max`, so `max_del` is *negative*. Getting this backwards produces fixtures that
   look right and test nothing.
6. **`set_mismatches()` stores a pointer** into `query.gene_alignments`. Anything that
   outlives or aliases that vector is a real bug; assert per-realization lists to catch it.

---

## 10. Naming and organisation

**One `TEST_CASE` per generic pattern, one `SECTION` per instance.** Not one per gene class.
V/D/J is the axis the refactor removes; a file organised on it has to be reorganised after
the collapse, and it hides gaps — you cannot see from a `V_gene` / `D_gene` / `J_gene`
layout that the Infeasible verdict is covered for V and J but not for D-5'.

```
TEST_CASE("Overlap verdicts (G2/G3)")
    SECTION("Infeasible — V 3' vs D 5'")
    SECTION("Safe — min deletions clear")
    SECTION("Undetermined — inside the deletion range")
```

As other events are migrated they add *sections* to these same test cases, and each
pattern's row fills in.

---

## 11. Known defects: `[!shouldfail]`

Do **not** pin a known-wrong value. Write the assertion the code *should* satisfy and tag
the case `[!shouldfail]`:

```cpp
TEST_CASE("DEFECT (plan 7.1): V credits an error-free length larger than its surviving core",
          "[gene_choice][iterate][endogenous][defect][!shouldfail]")
```

Catch2 reports an expected failure as a pass, so `ctest` stays green. When the defect is
fixed the case starts passing, which `[!shouldfail]` reports as a **failure** — forcing the
tag to be removed deliberately, so no fix lands unnoticed. A test stating the intended
behaviour also reads as documentation, which a pinned wrong number does not.

**Keep these cases section-free.** `[!shouldfail]` is evaluated per test-case *run*, and
Catch2 re-runs a case once per leaf section, so a case mixing passing and failing sections
reports the passing ones as unexpected passes.

Reference the plan section in the name, so the fix has an obvious landing site.

---

## 12. Checklists

### Adding a new event's tests

- [ ] All ten rows of §4 have at least one section (row 10 is automatic, but only if the
      fixture registers downstream events — see §4.1)
- [ ] Every `continue` / `break` in the body has a section **and** a positive control
- [ ] Each section mutation-verified: break the code, confirm *this* section fails
- [ ] Known-wrong behaviour is `[!shouldfail]`, not pinned
- [ ] Sections live under pattern-named `TEST_CASE`s, not gene-class ones
- [ ] `pixi run test_unit` green; expected failures reported as expected

### Changing tests during a refactor

- [ ] A section that now fails is a **behaviour change** — establish it is intended before
      editing the expectation
- [ ] A `[!shouldfail]` case that starts passing means a defect was fixed: remove the tag,
      and say so in the commit
- [ ] Rows do not disappear. If a generic body absorbs three branches, the three sections
      become three instances of one pattern — not one section
- [ ] Update the plan's delivered-coverage table
- [ ] `pixi run test_regression` for anything claiming to be behaviour-preserving; unit
      tests alone do not establish that
