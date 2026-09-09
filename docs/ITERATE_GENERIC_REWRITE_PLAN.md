# Generic `iterate()` Rewrite Plan — B5, B6, B7, B11

**Created**: Sep 1 2026
**Scope**: the four `Rec_Event::iterate()` implementations, plus the `initialize_event()` /
`initialize_Len_proba_bound()` machinery they depend on.
**Parent**: [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md)
tasks B5, B6, B7, B11. This document does not replace them; it supplies the design those task
stubs deferred.
**Companion**: [ITERATE_METHOD_ANALYSIS.md](ITERATE_METHOD_ANALYSIS.md) describes what the four
implementations *do*. This document describes what they have *in common*, and how each branch
collapses onto that.

---

## 0. Executive summary

The four `iterate()` bodies contain **13 hardcoded V/D/J or VD/DJ/VJ branches** across ~2 400
lines. They are not four different algorithms with incidental similarities: they are **one
algorithm** — *pick a realization, move or create a segment, check it against its neighbours,
bound the remainder, recurse* — instantiated four times with the topology inlined.

Nine generic patterns account for essentially all of the branching (§2). Three of them are
worth calling out up front because they change the shape of the work:

1. **The safety check is interval-vs-interval, and `Gene_choice` and `Deletion` differ only in
   whether one interval has collapsed to a point** (§2.2). One predicate replaces eight blocks.
2. **Only the *work* collapses, not the storage** (§2.3). Checking the nearest already-chosen
   segment on a side and propagating the verdict along the row is provably equivalent to checking
   all of them — but the *n(n−1)/2* pairs must still be stored, because the pair a flag refers to
   has to be explicit in its key, and event order is model data (`priority`), not a fixed VDJ
   sequence. A row-bitmask over `LayeredArray<uint32_t>` makes propagation O(1) and needs no new
   container.
3. **The junction-length probability machinery generalises to "any ordered pair of gene
   segments"** (§2.5). `VJ_ins_seq` is already used as a key meaning *the whole V→J span* in a
   VDJ model where no `VJ_ins` event exists. Naming that pseudo-junction properly is what makes
   `vj_length_d_position_proba` — the `no_d_align` driver, and the hardest part of B11 —
   generalisable at all.

The recommendation on sequencing (§5): **implement a reduced Phase A first**, because three of
the nine patterns need quantitative bounds (`len_min`/`len_max` per seq_type *and* side) that no
current interface exposes, and because building them ad hoc inside B5 would have to be undone in
Phase A anyway. Before even that, **port the `iterate()` unit-test harness** from
`feature/2_unittests` (§6.1) and **write its missing assertions** — `iterate()` has no unit tests
on this branch, and the sketch is a harness plus a partial coverage map, not a test suite: 31 of
its 69 assertions are in sections that never call `iterate()`, and the 12 sections with no
assertions at all are precisely the safety, `no_d_align` and pruning branches this work rewrites.
Keep the harness, drop most of the tests, and organise what remains by *pattern* rather than by
V/D/J — the sketch's layout reproduces the very `switch` being deleted.

The recommendation on ordering within Phase B (§6): **B11 → B6 → B7 → B5**, not the plan's
B7→B6→B5. B11 is already flagged as the milestone-1 blocker; B5 is the *hardest* and should be
last so it lands on top of finished shared services rather than inventing them.

---

## 1. Method: one event at a time, always runnable

Non-negotiable, per the standing constraint that every step must be independently verifiable:

- **Every step compiles and passes `pixi run test_regression`.** Unit tests alone do not
  discriminate here — they did not catch the VJ id-aliasing trap in B2, and they will not catch a
  pruning-bound change.
- **Bitwise-exact by default.** Each step's definition of done names either "bitwise identical
  Pgen / marginals on the regression corpus" or, if not, *exactly which* numbers change and why.
  §7 lists the four places where a naive generic rewrite would silently change results — those
  are the ones to watch.
- **Shared services land before their first consumer**, each with its own unit tests, in its own
  commit. A service with no caller is acceptable for one commit; a service with three callers
  landing in the same commit as all three is not.
- **Mutation-verify the discriminating tests.** As in B2 `f1d26a4`: weaken the predicate under
  test and confirm exactly the expected tests fail. A safety-check test that passes with the
  check deleted is worthless.
- **No event is migrated while another is half-migrated.** Each event goes from
  fully-switched to fully-generic in one commit (possibly preceded by service commits), so
  `develop` never carries two topology models at once.

### Verification ladder per step

| Level | Command | Gate |
|---|---|---|
| compile | `pixi run build` | always |
| unit | `pixi run test_unit` | always |
| integration | `pixi run test_integration` | always |
| regression | `pixi run test_regression` | always — this is the bitwise gate |
| convergence | `pixi run test_convergence` | **on every step that touches source** — not just pruning bounds. It is excluded from `pixi run test` and `test_unit`, and it is the only gate that caught §7.9 |
| benchmark | `pixi run benchmark` | on B5 and B11 (the two hot paths) |

---

## 2. The generic patterns

Notation used throughout:

- `id` — a `SeqTypeId`; `side` — `Five_prime` / `Three_prime`.
- `off(id, side)` — the current offset of that end, in read coordinates.
- A segment occupies read positions `[off(id,5'), off(id,3')]` inclusive; an empty segment uses
  the degenerate convention `off(id,3') == off(id,5') - 1` (B10).
- "left"/"right" mean 5'-ward / 3'-ward in the registry ordering.

### 2.1 — G1: Pending-modifier bounds (`len_min` / `len_max` per end)

**Where it is today**: 60 lines duplicated verbatim between
[Genechoice.cpp:1233-1299](../src/igor/Core/Genechoice.cpp#L1233-L1299) and
[Deletion.cpp:1537-1598](../src/igor/Core/Deletion.cpp#L1537-L1598) — four
copy-pasted blocks looking up the V-3′, D-5′, D-3′ and J-5′ deletion events in `events_map`,
checking `processed_events`, and caching `get_len_min()` / `get_len_max()` into eight scalars.

**What it actually computes**: for each segment end, *by how much can this offset still move
before the scenario is complete*. An end whose deletion event is already processed cannot move
(`0,0`); an end with no deletion event never moves (`0,0`); otherwise the bound is that event's
realization range.

**The obstacle**: `Rec_Event::len_min` / `len_max` mean three different things.

| Subclass | `len_min` | `len_max` | Units |
|---|---|---|---|
| `Deletion` | `-max_del` | `-min_del` | signed **length delta** |
| `Insertion` | `min_ins` | `max_ins` | non-negative **length** |
| `Gene_choice` | shortest template | longest template | non-negative **length** |
| `Dinucl_markov` | `INT16_MAX` (never set) | `INT16_MIN` (never set) | — |

Three incompatible conventions behind one accessor pair, and a subclass that leaves them at
sentinel values. Every consumer therefore has to know which subclass it is talking to, which is
exactly the coupling B5–B11 are meant to remove. **This is the single strongest argument for
doing Phase A first** (§5).

**Generic replacement**:

```cpp
struct OffsetDelta { int min; int max; };   // signed, applied to an offset
struct LengthRange { int min; int max; };   // non-negative, a segment length

// On Rec_Event (Phase A):
virtual OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const = 0;
virtual LengthRange get_segment_length_bounds(SeqTypeId)         const = 0;
```

`Deletion(X, 3')` returns `{-max_del, -min_del}` for `(X, 3')` and `{0,0}` elsewhere.
`Deletion(X, 5')` returns `{+min_del, +max_del}` for `(X, 5')` — note the **sign flip**, which
is precisely the asymmetry that today lives in the `d_5_min_offset = d_5_offset - d_5_min_del`
vs `v_3_min_offset = v_3_offset + v_3_max_del` idioms and that a reader has to reverse-engineer
each time. `Gene_choice` returns `{0,0}` for offsets and its template-length range for lengths.
`Insertion` returns `{0,0}` / `{min_ins, max_ins}`. `Dinucl_markov` returns `{0,0}` / `{0,0}`.

Built once per event at `initialize_event()` into:

```cpp
class PendingModifierBounds {          // one per Rec_Event instance, rebuilt at initialize_event
    // indexed [id * 2 + side], accumulated over every not-yet-processed event in events_map
    std::vector<OffsetDelta> offset_;
    std::vector<LengthRange> length_;   // indexed [id]
public:
    OffsetDelta offset_delta(SeqTypeId, Seq_side) const;
    LengthRange length(SeqTypeId) const;
};
```

**Deletes**: the eight `*_min_del` / `*_max_del` members from both `Deletion.h` and
`Genechoice.h`, and both 60-line lookup blocks.

#### Delivered (S2) *(Sep 7 2026)*

`JunctionGeometry::PendingModifierBounds` in the new `src/igor/Core/JunctionGeometry.h`, header-only,
no caller yet. Four things settled while writing it that the sketch above left open:

- **It never names a deletion event.** The accumulation loops over *every* unprocessed event in
  `events_map` and asks the A0 capability queries for every `(id, side)`. Only `Deletion` answers
  non-zero, so the result is bitwise the legacy lookup — but a topology with more than one
  modifier per end needs no new code, which is the whole point for tandem D. `Gene_choice` and
  `Insertion` contribute length and never travel; `Dinucl_markov` contributes neither.
- **Contributions sum.** `+=`, not last-write. Under today's topologies exactly one term per end
  is non-zero so the sum is the term; the test that pins the composition uses a synthetic map key,
  since `Model_Parms` keys events by `(type, seq_type, side)` and could not produce two.
- **`{0,0}` deliberately conflates two situations** — an end whose modifier is already processed,
  and an end with no modifier in the model. The legacy scalars conflate them too and no consumer
  distinguishes them. If one ever needs to, that is a separate query, not a third state.
- **The length sum is left unclamped.** A palindromic 5' deletion plus a 3' deletion can drive a D
  segment's lower bound below zero (`{-3, 10}` in the test model). It is a bound, not a reachable
  length; clamping would invent a semantics no consumer has asked for, and the junction-length DP
  applies its own floor. Pinned at the negative value so that adding a clamp later is a visible
  decision rather than a silent one.

The equivalence claim is tested against the legacy arithmetic itself, not against typed-in numbers:
`legacy_offset_delta()` in the test reads `get_len_min()` / `get_len_max()` off the same event
objects and applies the legacy `off + X_max_del` / `off - X_min_del` formulas per side. Literal
expectations sit beside it so that a sign flip on both sides of that comparison cannot pass.

**45 assertions in 7 `TEST_CASE`s. Eight mutations run, all caught**: the `Deletion` 3'/5' sign
flip, ignoring `processed_events`, last-write instead of `+=` (offsets and lengths separately),
`resize` instead of `assign` in `rebuild()`, dropping `side` from the index, and removing either
argument check. Full ladder green including regression and convergence — vacuously, since nothing
includes the header yet.

### 2.2 — G2: Reachable-offset interval, and the overlap predicate

**Where it is today**: **twelve** checks — `Gene_choice` V/D/J (six: V-3′ against D-5′ and J-5′,
D-5′ against V-3′, D-3′ against J-5′, J-5′ against V-3′ and D-3′) and `Deletion` V/D-5′/D-3′/J
(four sites, six checks). *(Counted from the code during S3; this section previously said eight
and five.)* Each computes a `*_min_offset` / `*_max_offset` pair and then applies a three-way
comparison. [Deletion.h:120-146](../src/igor/Core/Deletion.h#L120-L146) holds eight
scalars for this; `Gene_choice` holds the same eight again.

**The interval**:

```cpp
struct OffsetInterval { Seq_Offset lo, hi; };

OffsetInterval reachable(SeqTypeId id, Seq_side side, Seq_Offset current) const {
    const OffsetDelta d = bounds.offset_delta(id, side);
    return { current + std::min(d.min, d.max), current + std::max(d.min, d.max) };
}
```

~~The `min`/`max` is what absorbs the 5′-vs-3′ sign flip.~~ **Wrong, corrected in S3**: the sign
flip is absorbed by the *provider*. `Deletion::get_offset_delta_bounds` returns `{-max_del, -min_del}`
on a 3′ end and `{min_del, max_del}` on a 5′ one, both already ordered, so `reachable()` is a plain
translation and the `min`/`max` was dead code — a mutation removing it changed nothing any test
could see. See the S3 note below. The arithmetic itself checks out at every site:
`d_5_min_offset = d_5_offset - d_5_min_del` with `d_5_min_del = get_len_max() = -min_del`
gives `d_5_offset + min_del` = `current + delta.min` ✓, and symmetrically for the rest.

**The predicate**. Let `L` be the left segment's 3′ interval, `R` the right segment's 5′
interval, and `G` the minimum total length of everything that must sit strictly between them.
The geometric constraint is `L.3' + G < R.5'`.

```cpp
enum class Overlap { Infeasible, Safe, Undetermined };

Overlap check(OffsetInterval L, OffsetInterval R, int G) {
    if (L.lo + G >= R.hi) return Overlap::Infeasible;   // best case still violates
    if (L.hi + G <  R.lo) return Overlap::Safe;         // worst case still satisfies
    return Overlap::Undetermined;
}
```

**The unification.** `Gene_choice` and `Deletion` differ *only* in whether the moving end's own
interval has collapsed:

- In `Gene_choice::iterate`, the segment's own deletion is still pending, so `reachable(own)` is
  a proper interval. Hence `(v_3_off + v_3_max_del) >= d_5_max_offset`
  — an interval-vs-interval comparison.
- In `Deletion::iterate`, the event *is* the pending modifier and is being consumed, so after
  `initialize_event` marks it processed, `bounds.offset_delta(own, own_side) == {0,0}` and
  `reachable(own)` degenerates to the point `[v_3_new_offset, v_3_new_offset]`. Hence
  `v_3_new_offset >= d_5_max_offset` — the *same* predicate, degenerate.

There is at most one deletion event per `(seq_type, side)`, so the collapse is guaranteed, not
incidental. **One predicate replaces all eight blocks.**

`G` is a *lower* bound only; the corresponding upper bound is implicit in the junction-length
map's key set, and whether the two should be unified is recorded as an open item under §2.5.

`G` is `0` in every check performed today (all in-between segments have `len_min == 0`, and
`Insertion` events legitimately allow zero insertions). Keeping `G` explicit costs nothing and
is what lets the predicate stay correct once a tandem-D ordering puts a *gene* segment between
two checked ends. **Bitwise-preservation note**: `G` must be computed from `len_min`, which is
`0` for the current corpus, so this is a no-op today — do not "improve" it to a tighter bound in
the same step (§7.2).

#### Delivered (S3) *(Sep 7 2026)*

`PendingModifierBounds::reachable()` and the free `JunctionGeometry::check_overlap()`, in the same
header, no caller yet. `check` is named `check_overlap` because `JunctionGeometry::check(L, R, 0)`
says nothing at a call site. `gap` is a required argument, not defaulted to `0`: a default is an
invitation to forget it, and §7.2 is precisely about not letting it drift.

**The equivalence is tested by replay, not by restatement.** Writing the legacy comparisons in the
new vocabulary would make the test tautological, so the four idioms the twelve sites are written in
are reproduced in *their own variables* — `own_off + own_max_del >= other_max_offset` and the rest,
with the legacy negated-deletion scalars read off the same event objects — and the generic
predicate must agree at every offset across a sweep that crosses both verdict boundaries. A
separate section counts the verdicts the sweep produces, because two functions that both answer
`Undetermined` everywhere agree perfectly.

**One design change came out of mutation testing.** `reachable()` as sketched wrapped its bounds in
`std::min` / `std::max`; removing that wrapper failed no test, because no provider can return an
unordered `OffsetDelta` — the ordering is established by `Deletion::get_offset_delta_bounds`, not
recovered downstream. Silently re-sorting therefore bought nothing and would have hidden a provider
bug behind a merely-narrow interval. Replaced by:

- `min <= max` stated as an invariant on `OffsetDelta` and `LengthContribution` in `Rec_Event.h`;
- `PendingModifierBounds::rebuild()` throwing `std::logic_error` naming the event, query and
  seq_type when a provider breaks it, once per event per id at initialization;
- `reachable()` reduced to a translation.

Testing the guard needs a deliberately-misbehaving event, since no real subclass can trip it —
`UnorderedDeltaEvent` in the test file — plus a property section asserting that all four real
subclasses honour the invariant on every end of every segment. Without the fake the guard would be
unreachable code that no mutation could reach either, which is how the `min`/`max` got there.

**1441 assertions in 12 `TEST_CASE`s** (S2 and S3 together). Nine further mutations run, all
caught: both comparison operators loosened, each comparison reading the wrong interval end, the
gap dropped, the verdicts swapped, `reachable()` swapping its bounds or ignoring `current`, and the
ordering guard disabled. The eight S2 mutations were re-run against the changed `rebuild()` and are
still caught. Full ladder green, bitwise vacuously — nothing includes the header yet.

### 2.3 — G3: Storage stays *n(n−1)/2*; only the *work* collapses to O(1) per side

*(Revised Sep 1 2026 after review. The first draft of this section proposed collapsing storage to
one flag per left-member `SeqTypeId`. That is wrong — see "Why one slot per row is not enough"
below — and the design recorded here is the reviewed one.)*

**Where it is today**: `Event_safety`'s three values are the three unordered pairs of gene
segments. `V_del` consults both `VD_safe` and `VJ_safe`; `J_del` consults both `VJ_safe` and
`DJ_safe`. The enum enumerates all *n(n−1)/2* pairs, but each event only ever *computes* against
its nearest relevant counterpart — the far entries are filled in by other events.

**The two halves of the generalisation are separate**, and conflating them is what produced the
bad first draft:

| | today | generic |
|---|---|---|
| **storage** | *n(n−1)/2* (`Event_safety`, hardcoded to 3) | *n(n−1)/2*, runtime-sized from the ordering |
| **work per event** | 2 explicit checks against named neighbours | 1 nearest-neighbour resolution per side + a row-suffix write |

**Claim (the near check implies the far ones)**: take `X < Y₁ < Y₂` in the ordering with `Y₁`,
`Y₂` both chosen. Suppose `X.3' + G(X,Y₁) < Y₁.5'` and `Y₁.3' + G(Y₁,Y₂) < Y₂.5'`. Because `Y₁`
is chosen its length is known, and
`G(X,Y₂) ≤ G(X,Y₁) + len_min(Y₁) + G(Y₁,Y₂) ≤ G(X,Y₁) + len(Y₁) + G(Y₁,Y₂)`, so

```
X.3' + G(X,Y₂) ≤ X.3' + G(X,Y₁) + len(Y₁) + G(Y₁,Y₂)
              <  Y₁.5' + len(Y₁) + G(Y₁,Y₂)
              =  Y₁.3' + 1 + G(Y₁,Y₂)
              ≤  Y₂.5'
```

∎ So an event resolves its **nearest already-chosen neighbour** on each side, performs one
`check()` there (G2), and **propagates the verdict along the rest of the row**: establishing
`(V, D)` safe marks `(V, D1D2_ins)`, `(V, D2)`, `(V, J)` … safe in the same write.

**Why propagation never lets an invalid scenario through.** Suppose a completed scenario violates
`C(A,C)`, i.e. `A.3' ≥ C.5'`, on a cell that was filled by propagation from an established
`C(A,B)` with `A < B < C`. `C(A,B)` holds in the completed configuration, so `A.3' < B.5'`.
Offsets satisfy `B.5' ≤ B.3' + 1` — with equality exactly when `B` is empty, under B10's
degenerate convention `three_prime = five_prime − 1`; otherwise `B.5' ≤ B.3'`. Hence
`A.3' < B.5' ≤ B.3' + 1`, so `A.3' ≤ B.3'`, and with the assumption `B.3' ≥ A.3' ≥ C.5'`, i.e.
**`C(B,C)` is violated**. `C(B,C)` is checked by whichever of `B`, `C` is chosen second and
resolved by their deletion events, so the scenario is pruned there. ∎

**Corollary**: propagation can change *where* a scenario is pruned, never *whether*. Marginals
and Pgen are unaffected, which is what makes this eligible for a bitwise-exact step.

Note the empty-segment case is *included*, not an exception — this is one of the few places where
B10's degenerate offset convention earns its keep, and it is worth an explicit test rather than a
comment.

> **Invariant the proof rests on**: every pair adjacent in *chosen order* is checked by someone.
> True today (the second of the two to be chosen performs the check), but it is an assumption
> about model topology and event ordering, not a property of the container. **Phase C should
> validate it**, and B5 should carry a test that constructs an ordering where it would fail and
> asserts the model is rejected.

**Why one slot per row is not enough.** The first draft proposed keying by the left member alone,
on the grounds that "nearest chosen right of V" and "nearest chosen left of J" agree. They agree
*for a fixed chosen set* — but the chosen set grows with recursion depth, and event order is
model data (`priority`), not a fixed VDJ sequence. A single slot would hold "pair {V,J}" when
written at depth 1 and "pair {V,D}" when overwritten at depth 5, with the pair identity implicit
in a chosen set the reader has no access to. The triangular matrix makes the pair identity
explicit in the key, which is the property that makes the flag readable at all.

**Container**. Two shapes were considered; the second is recommended.

- *(a) Triangular-index wrapper over `LayeredArray<bool>`* — `n(n−1)/2` keys, index
  `i·n − i(i+1)/2 + (j−i−1)` for `i<j`. Straightforward, but a row-suffix write touches many keys,
  so an event needs a layer per cell: `memory_layer_safety_1/2` become a vector of
  `(cell, layer)` pairs.
- *(b) Row-bitmask over `LayeredArray<std::uint32_t>`* — **recommended.** Key = ordering position
  `i`; value = a bitmask over columns `j > i`. Storage is `n` words per layer. A row-suffix write
  is one bitwise op, so propagation is **O(1), not O(n)**:

  ```cpp
  // mark every column strictly right of k as safe, in row i
  mask |=  ~((std::uint32_t{1} << (k + 1)) - 1);
  // mark column k unsafe
  mask &= ~(std::uint32_t{1} << k);
  ```

  It needs **no new container** — `LayeredArray` is already tested (`f6f0101`) — and it keeps the
  current per-event shape: one `request_layer` per row touched, at most two rows per event, so
  `memory_layer_safety_1` / `_2` survive as two scalars with unchanged semantics. The cost is a
  documented **limit of 32 seq types** (6 today, 7 for tandem D); assert it at `freeze()`.

**Sizes**: VDJ ordering `[V, VD_ins, D, DJ_ins, J]` → 5 words; tandem D
`[V, VD1_ins, D1, D1D2_ins, D2, D2J_ins, J]` → 7 words. Compare `Event_safety`'s 3 fixed values.

**Amends the parent plan.** B8's *"What B8 does instead"* note says `Safety_bool_map` becomes
`LayeredArray<bool>` with the `Event_safety` enum kept as an opaque dense key of size 3. That
remains the correct B8 state; B5 replaces it with the row-bitmask form above, and `Event_safety`
is deleted then.

### 2.4 — G4: Segment write triple

Every realization loop ends with the same three writes, with the seq_type inlined:

```cpp
scenario.set_sequence_segment(<ID>, &new_str,           memory_layer_cs);
scenario.set_offset          (<ID>, <SIDE>, new_offset, memory_layer_offset_del);
scenario.set_mismatches      (<ID>, &mismatches_vector, memory_layer_mismatches);
```

Collapses trivially once `<ID>` is `this->seq_type_id` and `<SIDE>` is `this->event_side`. This
is the bulk of `Genechoice.cpp`'s 81 hardcoded `*_gene_seq` references and needs no design — but
it is also the part that must not be done piecemeal, because a half-converted branch writing one
of the three to the wrong key fails silently.

### 2.5 — G5: Junction identity and the length-probability bound

**Where it is today**: `vd_length_best_proba_map`, `dj_length_best_proba_map`,
`vj_length_best_proba_map` — three `map<int,double>` members on both `Gene_choice` and
`Deletion`, plus `junction_length_best_proba_map` on `Insertion`, built by
`initialize_Len_proba_bound()` walking the remaining model queue and driven by the
`has_effect_on(Seq_type)` predicate.

**The hidden generalisation, already present**: in a VDJ model there is no `VJ_ins` event, yet
`VJ_ins_seq` is used as a junction key — by `Insertion::has_effect_on`
([Insertion.cpp:430-442](../src/igor/Core/Insertion.cpp#L430-L442)),
`Dinucl_markov::has_effect_on` and `Gene_choice::has_effect_on`, all of which answer `true` for
`VJ_ins_seq` from a VD or DJ event. `VJ_ins_seq` there does not mean "the VJ insertion
sequence"; it means **the whole V→J span**. The code already has the concept of a junction
between two non-adjacent gene segments; it just spells it with an insertion seq_type.

**Generic form**: a junction is an **ordered pair of segments** `(A, B)` with `A < B` in the
ordering. Its length is `off(B,5') - off(A,3') - 1`. An event has effect on `(A,B)` iff it
changes that quantity: it moves `A.3'` or `B.5'`, or it creates/modifies a segment strictly
between them. Both conditions are derivable from the ordering plus Phase A's `OffsetRole` /
`SeqConstructionRole` — **`has_effect_on` becomes a non-virtual base method**, and the three
subclass overrides (all of which encode VDJ) are deleted.

```cpp
// Rec_Event, non-virtual, replaces three virtual overrides:
bool has_effect_on(SeqTypeId left, SeqTypeId right, const SeqTypeRegistry&) const;
```

Storage becomes one map keyed by pair, replacing seven named members:

```cpp
std::map<std::pair<SeqTypeId,SeqTypeId>, std::map<int,double>> junction_len_best_proba_;
```

populated only for the pairs the event's safety checks actually resolve to — which, by G3, is at
most two per event (one per side).

**Bitwise-preservation note**: the walk in
[Rec_Event.cpp:358-392](../src/igor/Core/Rec_Event.cpp#L358-L392) must not change. Only the
*keying* changes.

**Container choice (deferred).** `std::map<int,double>` is a poor fit: the keys are contiguous
small integers (junction lengths from `len_min` to `len_max`), the map is built once and then read
in the hot loop with `.count()` + `.at()`, and every lookup is a red-black-tree descent with a
pointer chase per level. A `std::vector<double>` with a stored offset and a sentinel for "length
not achievable" is a drop-in replacement: the build is a max-accumulation, which is
order-independent, so the result is bitwise identical.

Two constraints on doing it:

- **`vj_length_d_position_proba` is different** — its `map<int, vector<tuple<…>>>` values are
  *sorted by decreasing probability* and the `no_d_align` loop `break`s on the first prune
  ([Genechoice.cpp:622](../src/igor/Core/Genechoice.cpp#L622)). That ordering is load-bearing;
  only the outer `int` key is a candidate for the vector treatment.
- **The Tensor API** expected from a separate refactoring may subsume this. Re-check when it
  lands rather than building a bespoke container now.

Classified as **premature optimisation for this work**: keep `std::map` through steps 1–5, and
measure before changing it. The pair-keyed indirection introduced here is what makes the swap a
one-line change later, which is the actual point.

**The containers stage by dimension** *(Quentin, Sep 10 2026)*. The three objects G5 needs are 1-D,
2-D and 3-D, and each has a different natural home:

| | what | today | target |
|---|---|---|---|
| **1-D** | span profile: length → best proba | `std::map<int,double>` — a red-black descent with a pointer chase per level | `std::vector<double>` + stored offset, or `Matrix<double>(1,N)` ([Utils.h:291](../src/igor/Core/Utils.h#L291)). Build is a max-accumulation, so order-independent and bitwise identical |
| **2-D** | parent-indexed profile: conditioning context × length | does not exist | `Matrix<double>` natively — this is the cross-clique tightening (§6.9 R6) |
| **3-D** | `vj_length_d_position_proba`: total length × `(realization, left_len)` | `map<int, vector<tuple<…>>>` | the **Tensor API**, hence S4d |

Two *different* lookup costs are in play and should not be conflated: the length maps are
`std::map` (tree descent), while the enumeration additionally does an `unordered_map<string,…>`
**hash** lookup per candidate — `event_realizations.at(get<0>(*d_position_iter))` — which is fixed
by carrying a realization index instead of a gene name, independently of any container change.

**S4d is gated, and the gate is unlikely to open first.** The Tensor API is `feature/TensorLinalg`,
which needs `std::mdspan` and therefore C++23, and the parent plan records that the C++23 bump is
**currently blocked** by a sparrow `nullable_variant` visitor interaction (199/201 streaming tests
on C++23 against 201/201 on C++20). So S4d is recorded as a last optional step rather than
scheduled, and the 1-D `Matrix`/vector swap is available meanwhile without waiting for any of it.

#### The frame, the algebra, and what decomposition is retained

*(Established Sep 9 2026 — see §6.10 for the evidence and the incidental findings.)*

**The coordinate frame, which was written down nowhere.** State it in A0's vocabulary rather than
in V/D/J terms, because the VDJ reading bakes in two accidents of the current model set:

> **A span's length is the read-space distance between the *as-created* facing boundaries of the
> two segments that anchor it** — the offsets an `OffsetRole::Creates` event wrote, before any
> `OffsetRole::Modifies` event shifted them.
>
> Equivalently, and this is the same statement: **it is the signed sum of every event's
> per-realization `LengthContribution` over the region** — `Creates` positive, `Modifies` negative.

The two formulations agreeing is exactly why `Δ(r)` is the only thing that differs between the four
`iterate_initialize_Len_proba` bodies (§6.10).

**The anchor criterion is about offsets, not content** — `Insertion` is a content creator too, so
`SeqConstructionRole::Creates` does not pick one out. **But it is not a static capability either,
and two attempts to make it one were wrong.**

*First attempt: `OffsetRole::Creates`.* That only appears to work because `Insertion` currently
reports [`None`](../src/igor/Core/Insertion.cpp#L463), which is the defect R3 fixes (§6.9).

*Second attempt: split `Creates` into `Anchors` (set from an alignment) and `Derives` (computed
from neighbours),* on the argument that an insertion segment must never anchor because
span(V, VD_ins) would be trivially zero-width. **Also wrong** *(Quentin, Sep 10 2026)*, and the
reason matters: whether a particular span is zero-width is a **consumer** question — it asks
whether that consumer should be computing that span at all — not a property of the segment. The
zero-width is an artefact of today's `Insertion`, which requires *both* neighbours to already hold
offsets. An `Insertion` that **enumerated** its realizations could be processed before any D
choice, and span(VD_ins, DJ_ins) would then be non-zero and perfectly meaningful. Reasoning from
the current VDJ model and the current implementation produced a criterion that only holds for them.

**The criterion is dynamic:**

> A segment anchors a span, *for a given consumer*, iff its offsets have already been created at
> that consumer's position in the priority order.

Which is what the code already computes: `v_chosen` / `d_chosen` / `j_chosen` are
`exists && already in processed_events`
([Genechoice.cpp:1111](../src/igor/Core/Genechoice.cpp#L1111)). The reason only gene segments
anchor today is not a capability fact — it is that `Gene_choice` is the only event that writes
offsets at all.

It is still **statically resolvable**, which is what the cache needs: priority order is model data
fixed at load, so the consumer→anchors map is settled in `initialize_event()` and the set of spans
anyone will ever ask for is known before the first read.

So `Insertion` takes plain `OffsetRole::Creates` under R3 and **the enum does not grow**. The
question `Derives` was reaching for — *can this event run before its neighbours?* — is a **context
dependency**, D.9's `get_context_seq_types()`, and belongs there: today's `Insertion` depends on
both neighbours' offsets, an enumerating one would depend on neither, and that difference is
precisely what would let it be scheduled earlier.

**Two capabilities that are redundant today but must stay separate** *(Quentin, Sep 10 2026)*. No
event currently creates a sequence without also creating its offsets, so
`SeqConstructionRole::Creates` and `OffsetRole::Creates` always coincide. They remain different
statements — an event could construct content and leave placement to a downstream event — so keep
both queries rather than collapsing them. What ties them is a leaf invariant:

> **For every seq_type in the registry, both a sequence *and* its offsets must have been created by
> the time a scenario reaches a leaf** — not necessarily by the same event.

Half of that is already asserted: the debug-only leaf check in
[`Rec_Event::iterate_wrap_up`](../src/igor/Core/Rec_Event.cpp#L193) (`1794b5f`) rejects a scenario
carrying `int_undefined` — content allocated but never filled (§7.14). **The offsets half has no
counterpart**, and R3 is what makes it assertable at all: while `Insertion` writes no offsets, an
offsets-complete leaf check would fire on every scenario. Adding it belongs with R3.

**Two enforcement points, catching different things** *(Quentin, Sep 10 2026)*:

| | where | scope | catches |
|---|---|---|---|
| **1 — capability** | `Model_Parms::finalize()` | whole model, static, at load | a registered seq_type *no event declares* it will create — a tandem-D ordering listing `D1D2_ins` with no `Insertion` for it. Sufficient as the **runtime** guard |
| **2 — leaf** | `iterate_wrap_up()`'s debug check | whole scenario, dynamic | any gap at all, but only on paths actually reached, and only in a build carrying asserts |
| **3 — hand-off** | `RecordingEvent::iterate` | **one event**, per realization | the event under test not honouring its own declarations, on every path its unit suite reaches |

Complementary rather than belt-and-braces, and the second covers the class that has actually
bitten: a **path-dependent** gap, where the event writes on most branches and not on one, is
invisible to a static check by construction. That is precisely §7.9 (`no_d_align` recorded no
overlap verdict) and §7.12.

Note the key differs between the two halves. The sequence half is per `SeqTypeId`; the offsets half
must be per **(`SeqTypeId`, `Seq_side`)**, since `get_offset_role` is side-taking and the two ends
of a segment can be created by different events in principle.

**The dynamic check's worth is bounded by branch coverage.** An assert on a branch no test reaches
buys nothing, so it belongs with the *a* steps rather than as an item of its own — and
`Deletion::iterate` at 0 % (§6.4) means it would cover none of that body until 4a lands.

**Neither tier is live today** *(measured Sep 10 2026)*, but the gap is wiring rather than
anything structural:

- The capability check **cannot precede R3**: `Insertion::get_offset_role()` returns `None`, so
  every existing model would fail the offsets half at load.
- The **Debug build links and is green** — a scratch `-DCMAKE_BUILD_TYPE=Debug` tree builds
  clean, `ctest -L unit` passes 220/220 in 1.9 s and `-L integration` 5/5 in 3.2 s, asserts
  enabled throughout. So nothing about the assert is broken.
- But it is **compiled out of every build the §1 ladder runs**: `pixi run configure` passes no
  `CMAKE_BUILD_TYPE`, [CMakeLists.txt:39-42](../CMakeLists.txt#L39) defaults to `RelWithDebInfo`,
  and the configured cache carries `-O2 -g -DNDEBUG`. `test_debug` / `test_unit_debug` exist but
  are not in the ladder.
- And **the unit suite never reaches the leaf branch**, in any build. `RecordingEvent::iterate`
  "record[s] what the real next event would read, then stop[s]. No recursion, no error rate."
  The leaf runs only when the model queue empties, which the harness never lets happen. Only the
  integration tests — the real-inference smoke in particular — get there.

So `1794b5f`'s `int_undefined` check, the precedent the offsets half would follow, is exercised
**only by integration, and only in a build nobody runs.**

Adding `test_integration_debug` to the §1 ladder costs ~3 s per step and activates tier 2 on the
real inference path immediately. But **integration coverage is topology-dependent** — the
regression corpus is one TRB model, so `no_d_align` never fires and whole branches never execute —
so tier 2 cannot be relied on to reach every event. That is what tier 3 is for.

#### Tier 3: check the declaration at the hand-off *(Quentin, Sep 10 2026)*

`call_iterate_recording()` queries the event under test's A0 attributes and passes them to the
`RecordingEvent` constructor; `RecordingEvent::iterate()` — which already receives the scenario
exactly at hand-off — asserts that what the event *declared* it would write, it *did*.

This is strictly better placed than a leaf check for unit testing. A leaf needs every registered
seq_type created, which the single event under test cannot achieve alone; a hand-off check asks
only "did **this** event honour **its own** declarations", which is precisely what a unit test of
that event can establish. It also localises the failure to the event that caused it instead of
reporting that something, somewhere, is incomplete.

The harness is already shaped for it: the layer-contract check runs in the same place, collects
`LayerViolation`s, and `call_iterate_recording()` surfaces them — so *"every test going through
`call_iterate_recording()` is checked automatically, and a new event's sections inherit it without
writing anything"* already holds, and a capability check is a direct sibling.

What it can assert at hand-off:

| declared for a seq_type | assertion |
|---|---|
| `SeqConstructionRole::Creates` | the segment exists, at this event's own claimed layer |
| `SeqConstructionRole::Fills` | the segment exists **and holds no `int_undefined`** |
| `SeqConstructionRole::Modifies` | written at its own layer — *not* that the value changed, since a zero deletion is legal |
| `SeqConstructionRole::None` | this event did **not** advance that key's layer |
| `OffsetRole::Creates` for (id, side) | that end exists, at this event's own layer |
| `OffsetRole::None` for (id, side) | this event did **not** advance that end |

Two rows earn their place beyond the obvious:

- **`Fills` is a better home for the `int_undefined` invariant than the leaf.** At a leaf you learn
  that something, somewhere, is unfilled; at the hand-off you learn that the event which declared
  it would fill *this* segment did not. Same invariant, localised, and no complete fixture needed.
- **The `None` rows catch the converse** — an event touching a segment it never declared. That is
  §7.13's shape exactly: `Dinucl_markov` writing through the pointer `Insertion` stored, at a layer
  it never claimed.

A **static** sibling is available for free at `call_iterate()` time, after `initialize_event()` and
before `iterate()`: the layers an event *requested* must agree with what it *declares*. Requesting
a layer for a seq_type declared `None` is an inconsistency needing no scenario to detect, and the
harness already captures `layer_before_init` and `layer_baseline` at exactly the right moments.

**Scheduling.** Tier 3 is **test-only**, so it is not a behaviour change and escapes O7 — it can
land whenever, and does not wait for R3. Better still, it can land *before* R3 and pass:
`Insertion` currently declares `OffsetRole::None` and writes none, which is self-consistent (wrong,
but consistent). R3 then flips the declaration and the write together, with tier 3 keeping them
honest — which makes R3 self-verifying rather than reliant on its three `[!shouldfail]` cases
alone. Natural home is S4a or 4a, whichever comes first.

**Its honest limit**: tier 3 verifies *consistency*, not *correctness* of the declarations. An
event declaring `None` everywhere would pass trivially. Tier 1 is what makes the declarations
non-vacuous, by requiring some event to declare a creator for every registered seq_type. The three
together say: **someone declares it** (1), **the declarer does it** (3), and **nothing slips
through on an unexercised path** (2).

**Relationship to the existing layer check**: the layer check is *claimed ⇒ written*; tier 3 is
*declared ⇒ claimed ⇒ written* plus *undeclared ⇒ not claimed*. They overlap on declared keys and
diverge on undeclared ones — the layer check still catches an event that claims a layer while
declaring nothing at all. Keep both; neither subsumes the other.

**Consistent with B10, and worth saying so.** B10 already requires an absent segment to be an
*explicit written value* — a null or empty `Int_Str*` at the event's own layer — rather than the
absence of a write, precisely because the claimed/current layer split made the three states
distinguishable. So an absent segment still satisfies "created", and its offsets half is satisfied
by the degenerate convention `3' == 5' − 1`. No segment gets exempted from the invariant on
grounds of being absent.

**Worked instance (VDJ, V–D).** With V's and D's as-created facing boundaries at read positions
`v3` and `d5`, and deletions `dv`, `dd`, the post-deletion gap is `(d5 + dd) − (v3 − dv) − 1`. Under
the current `Insertion`, whose length is *derived* from that gap rather than chosen, this gives
`n_ins = (d5 + dd) − (v3 − dv) − 1`, hence `d5 − v3 − 1 = n_ins − dv − dd` — precisely what the
traversal accumulates (`+n_ins`, `−dv`, `−dd`) and what the consumers look up.

**Why the general form matters, and not only for tidiness** *(Quentin, Sep 9 2026)*. Writing the
frame as "the gap the insertion fills" makes the *shortcut* part of the definition. An `Insertion`
implemented by **enumerating** its realizations instead of deriving the one consistent count would
turn that equation from a definition of `n_ins` into a genuine constraint — `Σ contributions =
as-created gap` — which the derived form silently satisfies and an enumerating form must be
*checked* against. Stated as a signed sum over `LengthContribution`, the frame holds either way, and
the same is true of any future `Modifies` event that is not a `Deletion`.

Note what the frame implies about attribution: **a `Modifies` event belongs to the span, not to the
segment it modifies.** It moves a boundary away from where that boundary was created, so it widens
the span rather than shortening the segment as far as this key is concerned. Any factorisation that
attributes a deletion to its own segment is in a different frame and will not reproduce these keys.

The one visible inconsistency the frame would expose —
`get_deletion_effective_junctions(D_gene_seq, ·)` omits the V→J span while
`Gene_choice(D)::has_effect_on(VJ_ins_seq)` includes it — is harmless only because that map is
never read (§6.10, finding 2).

**The elementary factors** are therefore alternating, not per-seq_type. "Anchor" below is the
dynamic criterion just given. Because the consumer→anchors map is fixed at initialization, the
finest useful partition is nonetheless static: cut the ordering wherever *any* consumer anchors,
and every span anyone asks for is a contiguous fold over those pieces. Today the cuts fall exactly
at the `Gene_choice` segments — but nothing in the algebra depends on that, nor on the anchors
being genomic:

- `T_a[n]` — for each **anchor** segment: the best probability of it contributing as-created
  length `n`;
- `G_i[n]` — for each maximal **inter-anchor gap**: the best probability of that gap having
  read-space width `n`, folding everything strictly inside it — the insertion's contribution, the
  Dinucl `p^L` factor — together with **both** bounding `Modifies` contributions, which belong to
  the gap and not to the anchors they trim (see the frame above).

**The composition.** For anchors `A`, `B` with anchors `C₁…C_k` between them:

> span(A,B) = `G₀ ⊗ T_C₁ ⊗ G₁ ⊗ … ⊗ T_Ck ⊗ G_k`

where `⊗` is **max-product convolution**, `(X ⊗ Y)[n] = maxᵢ₊ⱼ₌ₙ X[i]·Y[j]`. It is associative
over non-negative reals, which is what makes the fold order-free and the intermediate results
cacheable. Checked against the current code: `k = 0` is `vd_length_best_proba_map`; `k = 1` is
`vj_length_d_position_proba`, whose `(gene, vd_len, dj_len)` triple is precisely the index triple
`(T_D, G₀, G₁)`. **The two objects are the same operator at `k = 0` and `k = 1`.**

Two variants, one of which is the forgetful projection of the other:

| | keeps | consumed as |
|---|---|---|
| `⊗ᵐᵃˣ` | the max per total length | the pruning bound |
| `⊗ᵉⁿᵘᵐ` | every composition, sorted by decreasing proba | the enumeration domain (`no_d_align`) |

**The retained decomposition is always three components, for any topology.** An anchor `G`
enumerating exhaustively between anchors `A` and `B` branches on exactly two things — which
realization, and where its 5' end sits — and its 3' remainder is then arithmetic. Everything past
the next anchor is marginalised into `span(G,B)` by `⊗ᵐᵃˣ`, because the events out there do their
own enumeration when they run. So the tuple is

> `( len(span(A,G)), G realization, len(span(G,B)), proba )`

and stays three-plus-proba whether one anchor or four sit between `G` and `B`. The current code
already says so: `get<2>` is never a free dimension — the commented-out block at
[Genechoice.cpp:609-617](../src/igor/Core/Genechoice.cpp#L609) recomputes it as
`j_offset − d_full_3_offset − 1`. It is memoised, not enumerated.

*(An earlier reading of this section claimed a tandem-D D1 would need a five-component tuple. That
was a conflation of the fold that **builds** a span with the decomposition **retained** in the
queried map. Corrected by Quentin, Sep 9 2026.)*

**Cost, restated correctly.** Per enumerating anchor the stored object holds
`|R| × |left range| × |right range|` entries, `|R|` being that anchor's realization count. That is
the same order as today, independent of how many anchors were marginalised into the right span,
because within one length bucket `(realization, left_len)` determines `right_len`. The k-fold product is **build-time work only**: `O(∏ ranges)` time
collapsing into a 1-D array of `O(range)`. Tandem D's exhaustive path is affordable on this
structure.

**What must generalise in the type is not the arity** but (i) the `std::string` gene handle, which
costs a hash lookup per candidate inside the hot enumeration, and (ii) the implicit binding of the
two `int`s to `vd_length_best_proba_map` / `dj_length_best_proba_map` by name — they must index
whichever two spans the enumerating gene actually sits between.

#### Relation to Phase D

G5's structure **is** Phase D's DP with the interface variable collapsed to a single integer and
the semiring relaxed:

| | G5 bound fold | Phase D scenario DP |
|---|---|---|
| semiring | **(max, ×)** | (+, ×) |
| interface variable | **one integer: length** | D.3's tuple — offset, boundary nt, frame phase, window |
| conditional dependencies | **relaxed** — `maxᵢ marginals[base + r.index + i·size()]` | carried exactly |
| sequence content | **discarded** | required (mismatches, error weighting) |
| query dependence | **none** — built once from the marginals | per read |
| separability | exact, everywhere | only where D.3's interface closes a boundary |

The single line that buys all of it is `real_max_proba = maxᵢ marginals[…]`. Maxing over the
conditioning dimension makes a child's contribution independent of its parent's realization, and
that independence is what lets the profile factorise per gap. Discarding sequence content removes
the mismatch coupling; having no query removes the third. Max-product over a superset of
decompositions dominates the true sum-product optimum, so the relaxation is a safe **upper** bound,
which is what pruning requires.

Consequence for how this is built: `fold(⊗, factors)` over a span key means Phase D **widens the
key and swaps the semiring** rather than starting over. That is the concrete reason to reserve
`SegmentSpan` now (§6.10), rather than an aesthetic one.

The one thing that does not carry over: `Dinucl_markov`'s `p^L` is exact under this relaxation only
because it genuinely depends on length alone. Under Phase D it becomes an interface-variable
dependency (D.3's `LeftNt` / `RightNt` rows), so the probability hook should take the accumulator
rather than a bare length — cheap now, and it is the seam Phase D's version attaches to.

#### Open: only the *lower* bound on a junction span is explicit

*(Raised Sep 7 2026 reviewing S3.)*

`check_overlap()` takes `gap` — the minimum number of nucleotides that must sit between two ends —
and encodes nothing about the maximum. That is not an oversight in the port: the legacy code has no
maximum-side check either. The upper bound exists, but it is **implicit in the key set of the
junction-length map**, which holds exactly the achievable gaps and whose largest key is
`max_ins − min_del_left − min_del_right`. A geometrically-possible-but-too-large gap is therefore
rejected at the map lookup, never earlier.

So the two mechanisms are doing overlapping jobs: `check_overlap` duplicates the lower side of the
map's feasibility test as an early-out (§7.6 proves the containment), and the map carries the upper
side alone, plus the probability bound.

**Full unification is not available.** Achievability is *set membership*, not range membership:
a gap is producible iff `L = ins − del_left − del_right` for some legal triple, and if any of the
three realization sets is sparse — nothing in the model format forbids an insertion event with
realizations `{0, 5, 10}` — the achievable set has holes. An interval test can only ever check the
convex hull, so the map stays the authority and any bounds check is strictly weaker than it.

**What *is* worth unifying is the hull.** `[G_min, G_max]` for a junction is derivable from A0 data
already: sum `length()` over the segments strictly between the two ends, and add the two ends'
`offset_delta`. Making that a `span_bounds(left, right)` alongside S4's junction pair key would
give both bounds one derivation and one place to be wrong, instead of one bound in the predicate
and the other buried in a map's keys. It is also the shape a dynamic-programming interface wants:
**an interval to iterate over, plus a per-length oracle** — which is exactly the hull and the map,
named as such.

**Not in this refactor**, for two reasons:

- there is no legacy counterpart to a maximum-side early-out, so it cannot be justified as
  preservation; by §7.6's second consequence these checks are optimizations, and adding one is a
  performance change to be measured rather than assumed;
- moving *when* a scenario is discarded is only bitwise-safe if nothing observable happened in
  between, and in `Gene_choice`'s V branch the `set_overlap_safety` writes sit between the two
  points. Almost certainly harmless — a discarded realization's flag is overwritten by the next
  one — but §7.9 is what happens when a write in one branch is reasoned about instead of checked.

Candidate home: with S4, which is already building the junction identity. Revisit when the Phase D
interface variables are settled.

### 2.6 — G6: Position enumeration for an unanchored segment (`no_d_align`)

**Where it is today**: [Genechoice.cpp:538-846](../src/igor/Core/Genechoice.cpp#L538-L846), two
sub-branches (with and without both flanking genes chosen), driven by `vj_length_d_position_proba`
built at [Genechoice.cpp:1436-1470](../src/igor/Core/Genechoice.cpp#L1436-L1470) as

```
junction_len = d_gene.size() + vd_len + dj_len
```

composed from `vd_length_best_proba_map` × `dj_length_best_proba_map` × the D realization set,
then sorted by decreasing probability so the inner loop can `break` on the first prune.

**What it is generically**: given a segment `X` that has no alignment, and its nearest chosen
neighbours `A` (left) and `B` (right), enumerate the ways `X` can be placed inside the `(A,B)`
span, ordered by decreasing best-achievable probability:

```
span(A,B)  =  len(A→X junction) + len(X) + len(X→B junction)
```

Under G5 this is exactly the composition `junction_len_best_proba_[(A,X)] × length_bounds(X) ×
junction_len_best_proba_[(X,B)]`, indexed by `span`. The identity `junction_len = d.size() +
vd_len + dj_len` — which the parent plan correctly flags as wrong for tandem D — is the `A=V,
X=D, B=J` instance of it.

For tandem D with both D present, `D1`'s neighbours are `(V, D2)` if D2 is chosen first or
`(V, J)` if not, and `D2`'s are `(D1, J)` or `(V, J)`. All four are handled by the same
composition. **The `junction_len` identity is not wrong for tandem D under this form; it was
wrong because `vd_len`/`dj_len` were hardcoded to `VD_ins_seq`/`DJ_ins_seq` rather than derived
from `X`'s own neighbours.**

The second sub-branch (neither neighbour chosen, the sliding-window loop at
[Genechoice.cpp:681-844](../src/igor/Core/Genechoice.cpp#L681-L844)) generalises the same way,
with the span bounded by the read rather than by a chosen neighbour.

**The fallback must become an explicit switch.** Today the exhaustive path exists *only* inside
`case D_gene`, so V and J never reach it — an accident of the switch, not a stated policy.
Generalising G6 would silently extend it to any `Gene_choice` with no alignments, which both
changes results and would be catastrophic for V and J (hundreds of templates × hundreds of
positions).

Add a per-`Gene_choice` boolean:

```cpp
bool exhaustive_position_fallback_ = false;   // Gene_choice member
```

- **Default reproduces legacy exactly**: set at construction — `true` for `Gene_class == D_gene`,
  `false` for V and J. Step 5 is then bitwise by construction.
- **Exposed as a tunable** thereafter, since "should this gene be searched exhaustively when the
  aligner finds nothing" is a per-model, per-locus judgement, not a property of being a D gene.
  Same treatment as the read-visibility switch in G7.
- Resolves open item **O6**: the fixture question dissolves — step 5 is tested by flipping the
  switch on a gene whose alignments are empty, not by defeating the aligner.

### 2.7 — G7: Mismatch-list trimming and palindrome construction

Two mirror-image pairs, four sites, ~180 lines:

| | 3′ end (V-3′, D-3′) | 5′ end (D-5′, J-5′) |
|---|---|---|
| positive deletion | `substr(0, size - k)`; drop mismatches `> new_off` | `substr(k, npos)`; drop mismatches `< new_off` |
| negative deletion | take last `k`, reverse+complement, **append**; new mismatches already ordered | take first `k`, reverse+complement, **prepend**; new mismatches need `sort` |

Collapses to one body parameterised on `event_side`. The `sort` asymmetry is real and must be
kept: appending at 3′ preserves the sorted invariant, prepending at 5′ does not.

#### Two distinct notions are entangled here, and both need switches

"Full deletion" in the parent plan and in the code covers **two independent conditions**, applied
inconsistently across the four branches:

| | (a) template exhausted | (b) nothing of the gene visible on the read |
|---|---|---|
| meaning | the scenario's segment reaches length 0 | the segment's surviving span falls outside the read |
| governed by | the model's allowed deletion range | the sequencing protocol / read coverage |
| V | [Deletion.cpp:308](../src/igor/Core/Deletion.cpp#L308) `size() > value_int` — **forbidden** | [:315](../src/igor/Core/Deletion.cpp#L315) `if (v_3_new_offset < 0) continue;` *"There should be at least one nucleotide of the V in the read"* |
| D 5′ | [:550](../src/igor/Core/Deletion.cpp#L550) `size() >= value_int` — **allowed** | [:558](../src/igor/Core/Deletion.cpp#L558) `if (d_5_new_offset >= int_sequence.size()) continue;` — marked `//THIS IS A TEMPORARY FIX //FIXME` |
| D 3′ | [:781](../src/igor/Core/Deletion.cpp#L781) `size() >= value_int` — **allowed** | *none* |
| J | [:1043](../src/igor/Core/Deletion.cpp#L1043) `size() > value_int` — **forbidden** | *none* |

Neither condition is uniform, and (b) is present for V as settled policy, for D 5′ as an
acknowledged temporary fix, and nowhere else. A generic body cannot pick one behaviour without
changing results, so B5 introduces **two independent, explicitly-named switches**:

```cpp
bool allow_full_template_deletion_;   // (a) per Deletion event, from the model
bool require_visible_nucleotide_;     // (b) per gene / global, from the sequencing setup
```

- **Defaults reproduce the table above verbatim**, including the asymmetries and the two missing
  (b) guards. Step 4 stays bitwise.
- (a) belongs with the model: whether a segment may contribute zero nucleotides is a modelling
  statement, and — per the parent plan's B10 warning — must not be conflated with *absence* of the
  segment.
- (b) belongs with the run: it depends on read length and protocol, exactly as the user notes.
  Making it explicit is what lets the D 5′ `//FIXME` be retired without guessing what it meant.

Filling in the two missing (b) guards, and deciding whether V and J should agree with D on (a),
are **modelling changes** to be made after step 4 with their own evidence — the same handling as
§7.1.

#### Deferred optimisation: precompute palindromes once per (read, gene)

The palindrome construction — take `k` nt from the trimmed end, reverse, complement, and count
mismatches against the read — is re-executed for every scenario reaching that deletion, for every
`k`. It is invariant over the scenario:

- the reverse-complemented string depends only on the gene template and the side. For a 3′
  palindrome after a 5′ deletion (the D case), the last `k` nt are untouched by 5′ trimming as long
  as `k ≤` remaining length, which the guard already enforces — so the template's suffix is the
  right source in every case;
- the mismatch positions depend only on the facing offset, which `Gene_choice` fixes per
  alignment and no deletion on the *other* side moves;
- successive `k` are nested, so a **prefix-mismatch-count array** answers every `k` in O(1).

So one reverse-complement string plus one prefix-count array per `(read, gene realization,
alignment offset, side)`, computed at alignment time or on first use, replaces the inner loop and
the `sort()` at [Deletion.cpp:652](../src/igor/Core/Deletion.cpp#L652) / [:1144](../src/igor/Core/Deletion.cpp#L1144).

**Out of scope here** — it is a pure performance change, it touches the aligner boundary, and it
should be measured on its own. Recorded because the G7 collapse is the moment the four copies
become one, which is the cheapest possible moment to add the precomputation later.

**Bitwise-preservation note**: the V branch guards with `size() > value_int` (full deletion
forbidden) while D 5′ and D 3′ guard with `>=` (allowed). The generic body must carry this via the
switches above, not unify it.

### 2.8 — G8: Endogenous mismatches and the error-rate bound

**What it computes**: the *core* of a segment — the positions that survive every remaining
deletion — is `[max reachable 5′ offset, min reachable 3′ offset]`. Mismatches inside it cannot
be explained away and set a floor on the error probability.

Under G2 this is one line: `core = { reachable(id,5').hi, reachable(id,3').lo }`. All three
current variants are instances:

| Branch | core start | core end |
|---|---|---|
| `Gene_choice` V | `v_5_off` (no pending 5′ del) | `v_3_off + v_3_max_del` ✓ |
| `Gene_choice` D | `d_5_off - d_5_max_del` ✓ | `d_3_off + d_3_max_del` ✓ |
| `Gene_choice` J | `j_5_off - j_5_max_del` ✓ | `j_3_off` (no pending 3′ del) |

The mismatch-counting conditions are all correct. **The credited match length is not** — see
§7.1; this is the one place where the generic form and the current code genuinely disagree.

### 2.9 — G9: Neighbour-derived length, seed and read window (B6 and B7)

`Insertion` and `Dinucl_markov` need no safety machinery at all; they need exactly one thing —
*who is next to me* — and both currently answer it with a hardcoded string comparison.

**Insertion** ([Insertion.cpp:166-210](../src/igor/Core/Insertion.cpp#L166-L210)): three
`std::string` comparisons on `this->seq_type` **inside the per-scenario hot loop**, each
selecting a hardcoded neighbour pair. Generic:

```cpp
const SeqTypeId left  = scenario.constructed_sequences.first_occupied_left (seq_type_id);
const SeqTypeId right = scenario.constructed_sequences.first_occupied_right(seq_type_id);
insertions = scenario.get_offset(right, Five_prime) - scenario.get_offset(left, Three_prime) - 1;
```

This is also a **performance fix**, not only a generality fix: three string compares per scenario
become two array lookups.

**Dinucl_markov** ([Dinuclmarkov.cpp:144-176](../src/igor/Core/Dinuclmarkov.cpp#L144-L176)):
`traversal_specs` already has the right shape — `{target, anchor, anchor_side}` — but the specs
come from a hardcoded `switch` at
[Dinuclmarkov.cpp:35-47](../src/igor/Core/Dinuclmarkov.cpp#L35-L47), and two residual `switch`es
survive inside `iterate()` for the per-junction indices array and memory layer. Generic:

- anchor = `first_occupied_left(target)` when `anchor_side == Three_prime` (forward fill),
  `first_occupied_right(target)` when `Five_prime` (reverse fill);
- seed nt = `anchor_seq.back()` / `.front()` respectively;
- read window from the anchor's facing offset, exactly as today;
- indices array and memory layer become per-instance members rather than a `switch` over three
  fixed arrays.

> **Latent bug this exposes.** `previous_seq.back()` / `.front()` is called with **no
> non-emptiness check**. ~~For VD (anchor V, full deletion forbidden) and DJ this is safe today.~~
> **Corrected in 2a (§7.12): full deletion of V is *not* forbidden — its `Deletion` branch has no
> size guard at all — so this segfaults on the current corpus, not only under tandem D.** It is
> reproduced by a `[.]`-hidden test. The
> `first_occupied_*` walk fixes this by construction, which is precisely the B7 worked example
> already unit-tested in `f1d26a4`. **B7 must include a test that reaches this state**, not only
> the container-level test.

---

### 2.10 — G10: Computation tiers, and the per-read warm-up

The palindrome case (§2.7) is one instance of a pattern that runs through all four bodies:
**computations sitting in the hot path whose inputs do not vary over the hot path**. Naming the
tiers is what turns "we could cache that" into a placement rule the rewrite can follow.

| Tier | Depends on | Runs | Hook today |
|---|---|---|---|
| **0 — model** | model parms, event ordering, `processed_events` | once per model | `initialize_event()`, `initialize_Len_proba_bound()`, `initialize_crude_scenario_proba_bound()` |
| **1 — read × event** | tier 0 + the query sequence and its alignments | once per read | **none — this is the gap** |
| **2 — scenario** | tier 1 + the offsets/sequences written by upstream events | once per scenario branch | `iterate()` |

Tier 0 is well served. Tier 2 is the genuine hot path. **Tier 1 has no home**, so everything in it
is currently computed in tier 2 and repeated once per scenario that reaches the event — which for a
`Gene_choice` or `Deletion` deep in the queue is the branching factor of everything above it.

#### What is currently in tier 2 but belongs in tier 1

| Computation | Site | Actually depends on |
|---|---|---|
| palindrome reverse-complement + its mismatches | [Deletion.cpp:396-410](../src/igor/Core/Deletion.cpp#L396-L410) and 3 mirrors | (gene template, side, k) and the alignment offset — §2.7 |
| endogenous-mismatch count and the credited core length | [Genechoice.cpp:336-348](../src/igor/Core/Genechoice.cpp#L336-L348) and 2 mirrors | (alignment, model deletion bounds) |
| `no_d_mismatches` recomputed at every slid D position | [Genechoice.cpp:627-634](../src/igor/Core/Genechoice.cpp#L627-L634), [:735-742](../src/igor/Core/Genechoice.cpp#L735-L742) | (read, D template, offset) — a sliding-window profile, currently rebuilt in full at each step |
| mismatch-list trimming under a deletion | 4 sites, §2.7 | (alignment, deletion value) — see below |

**The mismatch-list observation.** For a positive deletion the result is always a **contiguous
subrange of the incoming list**: V and D-3′ keep a prefix, D-5′ and J keep a suffix. Composing D's
two sides gives suffix-then-prefix — still contiguous. So for the all-positive case the entire
per-scenario mismatch state is a pair of indices into the alignment's own `std::vector<size_t>`
(which `get_all_mismatches()` already returns **by const reference** — no allocation there), not a
freshly assigned vector per scenario. Negative deletions prepend or append palindrome mismatches,
so the general state is `(base subrange, palindrome length)`; the mixed D case (negative 5′ then
positive 3′) trims a list that already carries a prepended palindrome block and is the one that
does not reduce to a subrange of the alignment list alone.

This is a **substantial** optimisation — it removes a vector assignment from the innermost loop of
the most-executed event — and it is **out of scope here**. It is recorded because the shape of the
G7 collapse determines whether it is a later one-line change or another rewrite.

#### The hook to set up now

A third lifecycle method on `Rec_Event`, sibling to `initialize_event()`:

```cpp
// Called once per query sequence, after alignments are available and before the first
// iterate() of the scenario tree. Default implementation does nothing.
virtual void prepare_for_query(const QuerySequenceContext& query) {}
```

`GenModel` calls it over the model queue once per read, in the same place it currently computes
alignments.

**What this plan commits to**: *placement*, not lifting. Each of steps 1–5 must leave its tier-1
computations (a) expressed as pure functions of their stated inputs, with no reads of
`ScenarioContext`, and (b) grouped so that moving the call site is a cut-and-paste rather than an
untangling. Nothing is actually hoisted into `prepare_for_query()` during steps 1–5 — hoisting
changes evaluation counts and would need its own benchmark and regression evidence, and doing it
concurrently with the topology rewrite would make a regression impossible to attribute.

Add to the definition of done for each of steps 1–5: *every computation in the new body is
annotated with its tier, and no tier-1 computation reads `ScenarioContext`.* That annotation is
the deliverable; the hoist is a follow-up whose cost is then one commit per item.

## 3. Shared services to build

Five new pieces, in dependency order. Each is one commit with its own tests.

| # | Component | Home | Replaces | Consumers |
|---|---|---|---|---|
| S1 | `OffsetDelta` / `LengthRange` + the two Phase A virtuals | `Rec_Event.h` | three conflicting meanings of `len_min`/`len_max` | S2 |
| S2 | `PendingModifierBounds` | new `src/igor/Core/JunctionGeometry.h` | 8 scalars × 2 classes + 120 lines of lookup | S3, B5, B11 |
| S3 | `reachable()` + `Overlap check()` | same | 8 comparison blocks | B5, B11 |
| S4 | Junction key + `has_effect_on(left,right)` non-virtual | `Rec_Event.{h,cpp}` | 3 virtual overrides + 7 named maps | B5, B6, B11 |
| S5 | Safety row-bitmask over `LayeredArray<uint32_t>`, indexed by ordering position (§2.3) | `ExplorationContext.h` | `Event_safety` enum | B5, B11 |

`JunctionGeometry.h` is deliberately a **new header, not an addition to `Utils.h`** — `Utils.h`
went 709 → 582 lines across B8 and should keep shrinking.

S4 is the step to watch for the `Events_map` keying limitation recorded in §9: its `has_effect_on`
should return a *set* of bearing events, not the first match, even though today the set is always
a singleton.

S2/S3 are per-`Rec_Event`-instance state rebuilt at `initialize_event()`, matching the parent
plan's B5 open item: *"plausibly a plain `std::vector<std::pair<Seq_Offset,Seq_Offset>>` indexed
by `SeqTypeId`, with no layers at all, rebuilt per call"*. This design agrees, with one
correction — it is indexed by `(SeqTypeId, Seq_side)`, not `SeqTypeId`, because the two ends of a
segment have independent pending modifiers.

---

## 4. Branch-collapse tables

### B11 — `Gene_choice` (13 branches → 0)

| # | Current branch | Collapses to | Pattern |
|---|---|---|---|
| 1 | `switch(event_class)` V/D/J at [:195](../src/igor/Core/Genechoice.cpp#L195) | — (deleted) | G4 |
| 2–4 | V/D/J "check D choice" / "check J choice" preambles | one loop over the ≤2 nearest chosen neighbours | G1+G3 |
| 5–7 | `vd_check` / `vj_check` / `dj_check` comparison blocks | one `check()` call per side | G2 |
| 8–10 | per-class offset writes (5′ and 3′) | `set_offset(seq_type_id, side, …)` | G4 |
| 11 | per-class junction-bound lookup (`vd_`/`dj_`/`vj_length_best_proba_map`) | `junction_len_best_proba_[(A,B)]` | G5 |
| 12 | per-class endogenous-mismatch window | `core = {reachable(5').hi, reachable(3').lo}` | G8 |
| 13 | `no_d_align` exhaustive path, both sub-branches | neighbour-derived span composition | G6 |
| — | `switch(event_class)` in `initialize_event` at [:1104](../src/igor/Core/Genechoice.cpp#L1104) | one block over `seq_type_id` / `event_side` | G1 |
| — | `switch` in `initialize_Len_proba_bound` at [:1383](../src/igor/Core/Genechoice.cpp#L1383) | pair-keyed build | G5 |
| — | `has_effect_on` switch at [:1311](../src/igor/Core/Genechoice.cpp#L1311) | base-class implementation | G5 |

`event_class` **stays** and remains the alignment-strategy key —
`query.gene_alignments` is keyed by `Gene_class` and D1/D2 correctly share one alignment set.
This is the parent plan's B0 intent and is not a defect to be removed.

### B5 — `Deletion` (7 branches → 0)

| # | Current branch | Collapses to | Pattern |
|---|---|---|---|
| 1 | `switch(target_seq_type)` at [:255](../src/igor/Core/Deletion.cpp#L255) | — | G4 |
| 2 | nested `switch(event_side)` for D at [:515](../src/igor/Core/Deletion.cpp#L515) | `event_side` used directly | G7 |
| 3–5 | three safety preambles + comparison blocks | one `check()` per side over nearest chosen | G2+G3 |
| 6 | positive/negative deletion × 5′/3′ (4 bodies) | one body parameterised on `event_side` | G7 |
| 7 | `get_deletion_effective_junctions()` V/D/J table at [:66](../src/igor/Core/Deletion.cpp#L66) | pair derivation from the ordering | G5 |
| — | `switch` in `initialize_event` at [:1406](../src/igor/Core/Deletion.cpp#L1406) | one block | G1 |
| — | `switch` in `iterate_common` at [:1268](../src/igor/Core/Deletion.cpp#L1268) | one block | G4 |

**Preserve deliberately**: the two-stage prune (first with `proba_contribution == 1` and
`break`, then with it and `continue`) exploits the decreasing-deletion iteration order and is a
real optimisation, not an accident. The `d_del_opposite_side_processed` flag and its two
`//THIS IS A TEMPORARY FIX //FIXME` guards must survive the rewrite unchanged; they are the only
thing preventing a double-counted deletion when both D ends are processed, and untangling them is
separate work.

### B6 — `Insertion` (3 branches → 0)

The whole `if/else` chain at [Insertion.cpp:166-210](../src/igor/Core/Insertion.cpp#L166-L210)
becomes ~12 lines (G9). `initialize_crude_scenario_proba_bound`'s
`switch(ins_seq_type)` at [:393](../src/igor/Core/Insertion.cpp#L393), which locates the paired
`Dinucl_markov` event, becomes an `events_map` lookup on `seq_type_id`.
`has_effect_on` at [:430](../src/igor/Core/Insertion.cpp#L430) goes to the base class (G5).
`update_event_name`'s `switch` at [:524](../src/igor/Core/Insertion.cpp#L524) — which still emits
the *legacy* names `"VD_genes"`/`"DJ_gene"`/`"VJ_gene"` — is B9's business, not B6's; leave it.

**Smallest of the four, and the natural first migration.**

### B7 — `Dinucl_markov` (4 residual switches → 0)

`get_dinucl_traversal_specs` becomes a registry walk; the `indices_array_for_target` and
`memory_layer_for_target` lambdas at [:127](../src/igor/Core/Dinuclmarkov.cpp#L127) and
[:140](../src/igor/Core/Dinuclmarkov.cpp#L140) become per-spec members; the `switch` caching
`vd_seq_size`/`dj_seq_size`/`vj_seq_size` at [:148](../src/igor/Core/Dinuclmarkov.cpp#L148)
becomes a per-spec field; `has_effect_on` and `iterate_initialize_Len_proba`'s three-way string
cascade at [:557](../src/igor/Core/Dinuclmarkov.cpp#L557) go to the base class (G5).

The three fixed `int[max_*_ins]` arrays allocated in `initialize_event` become one per spec —
note this is the *same* class of fixed-size-buffer bug that B2 removed from `Rec_Event`
(`int current_downstream_proba_memory_layers[6]`), and it should be closed the same way.

---

## 5. The Phase A question

The user's intuition is right, and it is stronger than "would benefit": **three of the nine
patterns cannot be expressed cleanly without it.**

G1 needs per-`(seq_type, side)` offset-delta bounds and per-`seq_type` length bounds. Today the
only interface is `get_len_min()`/`get_len_max()`, which carries three incompatible meanings and
is unset on `Dinucl_markov`. G5's `has_effect_on(left, right)` needs to know whether an event
*modifies an offset* or *creates a segment* — that is literally `OffsetRole` and
`SeqConstructionRole`. G6 needs both.

Building these ad hoc inside B5 means building the same accessors twice — and the existing
accessors are not merely overloaded but fragile: the accumulation that produces them is
order-dependent by construction (§7.4), latent today only because `unordered_map`'s hash order
happens to cooperate.

**Recommendation: implement a reduced Phase A first** — the two quantitative virtuals (S1) plus
`get_seq_construction_roles()` and `get_offset_roles()`. Defer `is_branching()`,
`is_multi_realization()`, `get_context_dependency()` and `get_context_seq_types()`: they gate
Phases C/D/E, not this work, and `SeqContextDependency` still needs the `CodonFrame` value that
hazard H8 identified.

This is a **change to the parent plan's D2 decision**, which reads *"Phase A is not a
prerequisite for tandem D"*. That remains true of Phase A *as a whole* — tandem D does not need
the capability matrix. It is not true of the two quantitative virtuals, which B5 and B11 both
require. Recording the narrowed dependency:

> **D2 amendment (proposed)**: Phase A remains off the tandem-D critical path except for
> `get_offset_delta_bounds()` and `get_segment_length_bounds()`, which B5/B11 require and which
> should land as a standalone step (A0) before S2.

Cost estimate: A0 is four small overrides plus the base declarations — smaller than either of
the 60-line lookup blocks it deletes.

---

## 6. Migration sequence

Ordering rationale: **services before consumers; smallest consumer first; hardest last.** This
inverts the parent plan's B7→B6→B5 in one respect — B11 moves to the front, because it is
already flagged as the milestone-1 blocker and because `Gene_choice` is the only event that
*creates* offsets, so every other event's neighbour queries read what it wrote.

| Step | Content | Gate | Bitwise? |
|---|---|---|---|
| **T0** | ✅ **done** — harness ported from `feature/2_unittests`, 11 non-`iterate()` cases dropped, 39 pattern-keyed sections + 4 `[!shouldfail]` defect cases (§6.1) | unit + mutation | n/a — tests only |
| **A0** | ✅ **done** — `OffsetDelta` / `LengthContribution` + four capability virtuals on all four subclasses | unit | yes — no caller yet |
| **S2** | ✅ **done** — `JunctionGeometry::PendingModifierBounds` in the new `JunctionGeometry.h`, unit-tested against a VDJ and a VJ model | unit + mutation | yes — no caller yet |
| **S3** | ✅ **done** — `reachable()` + `check_overlap()` in `JunctionGeometry.h`, replayed against all twelve current comparison sites | unit + mutation | yes — no caller yet |
| **1a** | ✅ **done** — `Insertion` characterization, 8 `TEST_CASE`s against the **unmodified** event | unit + mutation | n/a — tests only |
| **1b** | ✅ **done** — **B6**, `Insertion::iterate` generic (G9). Smallest, one hot-loop win. | full ladder | **yes** |
| **2a** | ✅ **done** — `Dinucl_markov` characterization, 12 `TEST_CASE`s; the empty-anchor case is `[.]`-hidden because it segfaults (§7.12) | unit + mutation | n/a — tests only |
| **2b** | ✅ **done** — **B7**, specs from the registry (G9) and per-spec buffers. Skip-empty walk **deferred to phase R**: it is §7.12's fix, not a refactor (§7.11) | full ladder | **yes** |
| **S4a** | `SegmentSpan`; `affects_length_of(SegmentSpan)` replacing `has_effect_on`; the queue-level filter restored and the per-body self-filter removed (§6.10). Also lands the **tier-3 hand-off capability check** in the harness (§2.5) — test-only, so it carries no bitwise risk | full ladder | **yes** |
| **S4b** | `length_delta` + `span_proba_factor` as **group** hooks; the four `iterate_initialize_Len_proba` bodies → one non-virtual traversal; `SpanAccumulator` replaces the `constructed_sequences` side channel (§6.10) | full ladder | **yes** |
| **S4c** | Span-keyed structure owned by the model; `⊗ᵐᵃˣ`; six members → one; each span built once instead of five times; `initialize_Len_proba_bound` de-virtualised; findings 2–3's dead code deleted. **Removes the tandem-D enum ceiling** — on the milestone-1 critical path | full ladder + init-time measurement | **yes** |
| **S4d** | Tensor-backed containers for the 3-D `no_d_align` structure — **gated on the Tensor API**, itself blocked on the C++23 bump (§2.5). Optional, performance only | full ladder + benchmark | **yes** |
| **3** | **B11a** — `Gene_choice` alignment path generic (G4, G2, G8, and G5 via S4a-c). Characterization already delivered by T0. **First production consumer of S3** | full ladder + benchmark | **yes**, except §7.1 |
| **4a** | `Deletion` characterization sections, including the zero-length junction T0 deferred. **Moved ahead of S5** (§6.8, F4) | unit + mutation | n/a — tests only |
| **S5** | Safety row-bitmask; row-suffix propagation; `Event_safety` deleted | full ladder + the empty-segment transitivity test + 4a's sections unchanged | **yes** (§2.3 corollary) |
| **4b** | **B5** — `Deletion::iterate` generic (all patterns). **First production consumer of S2** | full ladder + benchmark + convergence | **yes** |
| **5a** | `no_d_align` characterization beyond T0's G6 sections, on a fixture that *forces* the path. **Must raise `Gene_choice::iterate` block coverage** — see §6.4. Also lands the `bound / realized_proba` instrumentation (§6.10) | unit + mutation | n/a — tests only |
| **5b** | **B11b** — `no_d_align` exhaustive path generic (G6), including `⊗ᵉⁿᵘᵐ` — the retained decomposition, always three components (§2.5) | full ladder + a fixture that *forces* the path | **yes** |
| **R1–R4** | **Repair phase** (§6.9) — the decided behaviour changes, held here so everything above is idempotent end to end: §7.13, §7.12, `Insertion`'s three `[!shouldfail]` defects, `dinuc_proba_matrix` → `initialize_event()` | full ladder, per commit | **no** — golden data may move; each commit names which outputs and why |
| **R5** | §7.1 and §7.8 off-by-one corrections, per decision O4; the four `[!shouldfail]` tags come off | full ladder + the corrected-core unit tests | **no** — same |
| **R6** | Within-clique joint max in the span fold (§6.9); optional cross-clique parent indexing | full ladder, **convergence weighted heavily** | **no** — a tighter bound prunes more |

**The split is an ordering requirement, not bookkeeping.** The *a* commit lands before the *b*
commit and is mutation-verified there, where mutation-verification means something: a
characterization suite written after its collapse pins the new behaviour and agrees with the
refactor because it was derived from it. §7.9 is the concrete argument — the regression corpus has
a TRB topology that never fires `no_d_align`, so "collapse, then assert non-regression, then test"
has a hole exactly where the risk is. Per-event unit tests are the only cover for branches the
corpus does not reach, and they are only evidence if they predate the change.

**Each *a* commit reports coverage of the two functions it is characterizing**, via
`pixi run coverage` (see §6.4). The suite is only evidence for the branches it reaches, so an
uncovered branch is a gap to close deliberately or to record deliberately — not something to
discover after the collapse. Writing 1a's report immediately exposed one: the discard path was
covered on the VD arm only, so a rewrite could have lost two of the three guards undetected.

A *b* commit's definition of done is "its *a* sections pass unchanged" — or, where the collapse
legitimately reorganises a section, each adapted assertion is named in the commit message. Some
adaptation is expected: the *a* sections are written against the legacy switch structure, and the
pattern-keyed layout (§6.1) minimises but does not eliminate the churn.

B11 is split: the alignment path (3) is what milestone 1 needs and is straightforward; the
exhaustive path (5) is the hardest single piece in the whole set and depends on S4 being settled.
The parent plan's guidance to build the milestone-1 fixture so `no_d_align` never fires stands —
step 5 is exercised by handing the D `Gene_choice` an empty alignment list, which needs no aligner
manipulation. The G6 fallback switch is still built — it is what preserves the V/J-vs-D asymmetry
once the `case D_gene` that carried it is gone (§2.6, O6).

**Every behaviour change is held to phase R, after 5b.** The rows above are bitwise-exact by
construction, so `pixi run test_regression` means the same thing at every one of them and the
golden data never moves inside the refactoring block. Interleaving a fix — even a cheap one landing
on an event that is freshly migrated and fully covered — would require regenerating golden outputs
partway through, and from that point on a "regression" is no longer a single unambiguous signal:
every later step's verdict has to be read against which baseline it was taken on. The cost of
deferring is re-establishing context on `Insertion` and `Dinucl_markov` later; the cost of not
deferring is the gate itself, which §6.2 already shows is the fragile part. R is where golden data
is allowed to move, once, deliberately, with each commit naming the outputs it changes.

### 6.1 — T0: build a focused `iterate()` characterization suite

`iterate()` has **no unit tests on this branch**, and the migration needs branch-level
non-regression that the corpus run cannot give: a Pgen delta says a step broke something, not which
branch.

#### What `feature/2_unittests` actually provides

A harness sketch and a partial coverage map — **not** a test suite. Measured on
`test_gene_choice_iterate.cpp` (1005 lines, 12 `TEST_CASE`s, 33 `SECTION`s, plus 4 `TEST_CASE`s
that are names with no body):

| | count |
|---|---|
| live `REQUIRE` | 69 |
| … of which are in sections that **do not call `iterate()` at all** | **31** |
| commented-out `REQUIRE` | 24 |
| `TODO` markers | 129 |
| sections with zero live assertions | 12 of 33 |

Nearly half the live assertions test constructors, the mock helpers, or the harness itself. That
is the "trash" fraction, and porting it would both dilute the effort and leave a suite whose
apparent coverage is mostly self-referential.

#### Port list

**Take the harness** — `IterateTestState`, the five `*Storage` structs, `call_iterate()`, and the
accessors (`get_seq_offset`, `get_constructed_sequence`, `get_mismatches`, `is_safe`,
`create_stub_gene_choice`, `create_events_map`, `create_perfect_alignment`). This is the reusable
design and the only part worth taking wholesale. Adaptation needed, since the branch forks at
`f393974`, before B2/B8:

- the `*Storage` structs need a `SeqTypeRegistry` — `Seq_type_str_p_map`, `Mismatch_vectors_map`,
  `Downstream_scenario_proba_bound_map` and `Seq_offsets_map` are all registry-constructed since B8;
- `ExplorationStorage` gains `Pruning_mismatch_floor_map`;
- `Events_map` is keyed by the seq_type **name string** since B4, not `Gene_class`;
- `create_stub_gene_choice` / `create_events_map` must set `seq_type` and `seq_type_id`.

**Drop outright** — 11 `TEST_CASE`s worth, none of which exercises `iterate()`:

| Dropped | Why |
|---|---|
| `Gene_choice event construction and realizations` (2 sec, 14 asserts) | tests `add_realization`; not this code path |
| `Gene_choice draw_random_realization` (1 sec, 2 asserts) | generation path, not inference |
| `Alignment_data construction` (2 sec, 9 asserts) | tests the mock helper — tests the test |
| `IterateTestState creation` (2 sec, 6 asserts) | tests the harness — tests the test |
| `Single V gene with single alignment …` (4 asserts) | duplicate of `Basic V alignment - no other genes` |
| `V gene with multiple alignments`, `Multiple alignments with different positions`, `Multiple J alignments` | three sections for one thing: loop iteration over the alignment list. Keep **one**, folded into the baseline |
| `VD safety` (0 asserts, 1 TODO) | an empty wrapper around its two child sections; keep the children |
| the 4 body-less `TEST_CASE`s (`empty alignment list`, `sequence boundary conditions`, `zero-length junction`, `all three genes with complex junctions`) | names with no content. Two of the *behaviours* are worth having and are written fresh below; do not port empty shells |

#### Structure the suite on the axis of the rewrite, not on V/D/J

The sketch is organised as `V_gene iterate` / `D_gene iterate` / `J_gene iterate` — i.e. it
reproduces exactly the `switch` this work deletes. Ported as-is it would have to be reorganised
after step 3, and it hides gaps: it is not visible from that layout that the Infeasible verdict is
covered for V and J but not for D-5′.

**Organise each `TEST_CASE` by pattern, each `SECTION` by instance.** The file then survives the
collapse unchanged, and a missing instance is a visible hole in a row.

| `TEST_CASE` | Pattern | Sections (instances) | From |
|---|---|---|---|
| **Overlap verdicts** (G2/G3) | the three-way `check()` outcome | `Infeasible — V 3' vs D 5'`; `Infeasible — D 3' vs J 5'`; `Infeasible — J 5' vs V 3'`; `Safe — min deletions clear`; `Undetermined — inside the deletion range`; `pair not adjacent in the ordering (V vs J across D)`; `counterpart not yet chosen ⇒ no check` | `VD safety check: overlap detected`, `DJ safety overlap detection`, `Overlap causes skip`, `safe with min deletions` (4 asserts), `unsafe in deletion range` (4 asserts), `VJ safety checks (no D gene)`, `VJ safety check (no D)`, `DJ safety check`, `Alignment-based with VD safety check` |
| **Junction-length bound** (G5) | length lookup miss ⇒ `continue`; hit ⇒ bound set | `achievable length sets the bound`; `unachievable length discards the scenario` | `Junction length probability check` |
| **Exhaustive position fallback** (G6) | the `no_d_align` path and its switch | `fallback off ⇒ no realization enumerated` *(new — pins the G6 default for V/J)*; `both neighbours chosen — position map path`; `neither chosen — sliding path`; `mismatch profile per position` | `Basic D exhaustive search…`, `Exhaustive search with V and J chosen`, `Exhaustive search sliding D position`, `Mismatch computation in exhaustive search` |
| **Endogenous mismatch bound** (G8, §7.1) | the surviving-core window and its credited length | `V — the surviving core`, `J — the same slip`, `D — both ends truncated` (the §7.1 defect asserted at its intended value, `[!shouldfail]`) | `Alignment with mismatches` (1 assert), `Both ends truncated (max deletions)` (1 assert) |
| **Pruning** | `should_prune` at the realization boundary | `below threshold ⇒ realization skipped` | `Probability threshold filtering` |
| **Baseline writes** (G4) | the offset/sequence/mismatch triple | `V`; `D`; `J`; `negative alignment offset (pre-alignment strip)`; `zero-length junction between adjacent segments` *(new)* | `Basic V/D/J alignment - no other genes` (5/3/3 asserts), `V gene with negative offset` |

**≈22 sections in 6 `TEST_CASE`s**, against 33 sections in 12 `TEST_CASE`s. Two are written fresh
(the G6 switch default, and the zero-length junction — the latter because B10's degenerate offset
convention is what makes §2.3's propagation proof hold at the equality case, and nothing currently
pins it).

#### Assertions come from the running code — but a confirmed defect is asserted at its intended value

Characterization here means the expectation is derived by **running the unmodified
implementation**, never from the sketch's comments and never from what the code looks like it
should do. That settles where the number comes from. It does not settle what to assert once the
number turns out to be wrong, and the two questions are separate:

| The behaviour is | Assert | Marked |
|---|---|---|
| correct | the observed value | — |
| odd, but intent not yet established | the observed value | a comment naming precisely what is undecided |
| a **confirmed** defect — diagnosed, and agreed not to be intended | the value the code *should* produce | `[!shouldfail]`, section-free |

**Never pin a value that is known to be wrong as the expectation under test.** Two reasons:

- **It misreads in the diff.** A pinned wrong value has to be *edited* when the defect is fixed,
  and an edited assertion inside a refactoring diff is indistinguishable from a test accommodated
  to a behaviour change — the one signal this suite exists to keep trustworthy. Removing a
  `[!shouldfail]` tag is a one-line, self-describing change that cannot be misread as anything
  else.
- **It inverts what the test teaches.** A case asserting the intended value is a specification
  that happens to fail today; a case asserting the wrong value hands the next reader the wrong
  invariant, stated with the authority of a green suite.

The objection this rule replaces — that asserting intent would make the suite red before step 1 —
does not hold. Catch2 reports an expected failure as a pass, so `ctest` stays green; the case turns
red only when the defect is fixed and the tag is not removed, which is the point.

Two boundaries on it. First, **confirmed** is a real gate: an unestablished suspicion asserted
under `[!shouldfail]` is a guess that will be read as a decision, so behaviour that is merely
strange gets the middle row, not the last. Second, a wrong value may still appear as a **premise**
inside a case whose subject asserts intent — `test_event_capabilities.cpp` pins
`get_len_max() == INT16_MIN` precisely to contrast the legacy accessor with the correct capability
query beside it (§7.4). What the rule forbids is the wrong value standing as the expectation.

The 21 live assertions in the sections being kept are a starting point, not a baseline — re-derive
each against the current implementation, since the branch predates B2/B8.

#### Delivered

| `TEST_CASE` | Sections | Notes |
|---|---:|---|
| Baseline writes (G4) | 3 | V, D, J |
| Alignment mismatch propagation | 6 | V, D, J verbatim; outside-core kept; per-realization; empty-not-absent |
| Realization branching and probability | 6 | one call per alignment; no compounding; incoming proba; `base_index + realization`; one gene at two placements; zero-probability realization |
| Template overhangs (G4 → B3) | 3 | V negative offset, J past read end, and the asymmetry itself |
| Overlap verdicts (G2/G3) | 6 | Infeasible + control, Safe, Undetermined, not-chosen, no-counterpart |
| Junction-length bound (G5) | 2 | achievable / unachievable |
| Endogenous-mismatch bound (G8, §7.1) | 4 | V slip, V none, J slip, D off-by-one |
| Pruning | 2 | below threshold + control |
| Exhaustive fallback (G6) | 9 | V off, J off, sliding baseline, sliding with 5' / 3' deletions and both clamps, position map, per-position mismatches |

**39 sections in 9 `TEST_CASE`s plus 4 `[!shouldfail]` defect cases**, 160 assertions. Nine
mutations run: the safe-verdict force, the junction-guard removal, the §7.1 sign correction,
compounding the scenario probability across realizations, dropping `base_index` from the marginal
read, trimming the mismatch list to the surviving core, and removing either sliding-window clamp
are all caught; the overlap `continue` removal is not, which is §7.6.

Known defects are **not** pinned at their wrong values. Each gets a section-free
`[!shouldfail]` case asserting what the code should do — §7.1 for V, J and D, and §7.8 for the
position-path offset. They read as specifications, they fail today, and fixing the defect makes
them pass, which `[!shouldfail]` reports as a failure until the tag is removed.

Deferred to step 4 rather than written here: the zero-length-junction section. `Gene_choice`
never writes an empty segment — a genomic template is never empty — so B10's degenerate offset
convention can only be exercised from `Deletion`.

> **Writing the sections is documented separately.** [ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md)
> carries the harness reference, the nine-row matrix every event's tests must fill, the setup
> recipes, the traps, and the `[!shouldfail]` convention. Steps 1, 2 and 4 follow it rather than
> re-deriving the shape.

#### Scope discipline

T0 covers `Gene_choice` only, matching the sketch's reach. Sections for `Insertion`,
`Dinucl_markov` and `Deletion` arrive as the *a* commit of steps 1, 2 and 4 — each step contributes
the instances for the branches it is about to collapse, into the pattern-keyed `TEST_CASE`s above.
The suite grows with the migration rather than becoming a project of its own, and the row for each
pattern fills in as the events are migrated.

Staging them per event is about **harness churn, not about deferring the testing**: T0 rewrote the
layer-contract probe three times before the `LayeredArray` split made it work (§7.10), and had four
event suites existed by then, all four would have carried the broken probe. `Dinucl_markov` in
particular needs fixture support the harness does not have yet. What is *not* deferred is the
ordering — see the note under the §6 table: each event's sections land **before** its collapse, in
their own commit, mutation-verified against the unmodified event.

### 6.3 — Delivered (1a): the `Insertion` characterization *(Sep 7 2026)*

`tst/igor/Core/test_insertion_iterate.cpp`, **80 assertions in 11 `TEST_CASE`s** (3 of them
`[!shouldfail]`), written against the unmodified event. All ten matrix rows filled; two change shape because the event never
branches — row 2 becomes "one hand-off, or none", and row 3's *do not compound* has nothing to
say. Row 5 has no chosen/unchosen distinction either: `Insertion` never asks whether a neighbour
was chosen, it reads the offsets and assumes, so the section pinning what happens when they are
absent stands in for it (it throws, which is `LayeredArray` turning the assumption into an error
rather than a read of uninitialized storage).

Four things pinned that B6 has to reproduce, none of them obvious from the body:

- **It creates a segment and assigns it neither offsets nor a mismatch list — confirmed defects,
  not behaviour to reproduce.** The span being *derived* from the neighbours is the shortcut G9 generalises, and it
  holds only while the error model forbids indels; but a derived offset is still an offset, and
  requiring a consumer to know that insertion segments have no offsets is exactly the coupling this
  refactor exists to remove. Two `[!shouldfail]` cases state the intended behaviour: the junction
  occupies the positions strictly between its neighbours, and an empty one uses B10's degenerate
  `off(3') == off(5') - 1` convention rather than absence — which is the three-state problem B10
  has to solve, appearing here concretely. Nothing in Core reads or writes `seq_offsets` for an
  insertion seq_type today (checked), so **B6 can fix this without touching any consumer**, and the
  fix is expected to stay bitwise.
- **The two `new_index` derivations agree.** VD and DJ compute it as
  `base_index + event_realizations.at(to_string(n)).index` — a string conversion and a hash lookup
  in the hot loop, carrying its own `FIXME` — while VJ uses the `realization_index` that
  `iterate_common()` has already resolved. Pinned equal through `add_to_marginals()`, the only
  reader of `new_index`, so B6 can keep the second and delete the first.
- **A zero-probability length is discarded exactly like an unreachable one.** `proba_contribution
  != 0` is the single gate for both. A generic body that separates "not a realization" from
  "probability zero" changes which scenarios reach the next event.
- **No layer promise is owed on a discard path.** The event claims its `downstream_proba_map`
  layer once at `initialize_event()` and leaves it unwritten when the scenario is dropped, which
  is sound only because there is no hand-off. Claiming per call instead of once would turn this
  into the §7.9 defect.

**The missing mismatch list is the same defect, and is confirmed too** *(Sep 7 2026)*. A constructed
sequence implies a comparison against the read, so a constructed segment should carry a list —
empty at this point, since the junction holds placeholders and nothing is decided yet. Absence
forces every consumer to know that insertion segments are exempt, which is the coupling being
removed. It is load-bearing rather than tidy: under amino-acid Pgen a placeholder position scores
differently under ceiling and floor mismatch semantics, and either choice needs a list to write
into — an absent list cannot express *no mismatches yet* as distinct from *not compared*, the same
three-state problem B10 has for sequences. Third `[!shouldfail]` case; B6 fixes all three together.

§7.7's shortcut is pinned as observed, not as a defect: a negative junction length is discarded by
the *realization lookup*, not by geometry, and the section says so. The guard B6 replaces is a
set-membership test, which is the same distinction §2.5's open item turns on.

Harness additions: `make_insertion()` and `make_dinucl_markov()`. The latter is not optional —
`initialize_crude_scenario_proba_bound()` looks the Dinucl_markov up in `events_map` and throws
without it, so an `Insertion` cannot be initialized alone.

**Nine mutations run, all caught**: the junction length off by one, the neighbour pair swapped,
out-of-range lengths no longer discarded, the scenario probability not updated, `base_index`
dropped, the placeholder count wrong, pruning disabled, the downstream bound replaced by a
constant, and the VJ arm reading the wrong span.

### 6.4 — Coverage of the functions being migrated *(measured Sep 7 2026)*

`pixi run coverage` builds the instrumented tree and reports line, branch and block coverage per
function for the four `Rec_Event` subclasses' `iterate()` and `initialize_event()`. Scope it with a
Catch2 spec: `python3 scripts/tests/coverage_report.py -f '[insertion][iterate]'`.

**The tooling was broken and is fixed here.** `ENABLE_COVERAGE` added `-fprofile-generate`, which is
PGO instrumentation: it emits `.gcda` arc counts but no `.gcno` notes, so `gcov` had nothing to map
them onto. The `coverage` task then called `llvm-profdata` / `llvm-cov`, which are neither installed
nor matched to the GCC toolchain. Every report this project could have produced was empty. Now:
`--coverage` under GNU, the instrumented-Clang flags under Clang — branching on the *compiler*, not
the OS, which a Linux Clang build also got wrong.

Baseline under the unit suite (`~[integration]~[slow]`), before any event is collapsed:

| Function | Lines | Branch | Blocks | Calls |
|---|---:|---:|---:|---:|
| `Gene_choice::iterate` | 71.4% | 46.5% | 59.2% | 42 |
| `Gene_choice::initialize_event` | 86.8% | 51.8% | 75.3% | 42 |
| `Insertion::iterate` | 97.6% | 65.4% | 71.9% | 25 |
| `Insertion::initialize_event` | 100.0% | 62.5% | 77.8% | 26 |
| **`Deletion::iterate`** | **0.0%** | **0.0%** | **0.0%** | **0** |
| `Deletion::initialize_event` | 93.0% | 57.4% | 79.8% | 36 |
| **`Dinucl_markov::iterate`** | **0.0%** | **0.0%** | **0.0%** | **0** |
| `Dinucl_markov::initialize_event` | 94.7% | 43.3% | 54.1% | 25 |

Two readings matter more than the numbers themselves:

- **`Deletion::iterate` and `Dinucl_markov::iterate` are at zero.** No unit test executes them at
  all — steps 2a and 4a start from nothing, and until they land the *only* thing standing between a
  change in those bodies and a wrong answer is the regression corpus. Their `initialize_event`
  figures are non-zero only because other events' fixtures initialize them as neighbours.
- **Branch percentages read low and should not be chased to 100%.** gcov counts an exception edge
  as a branch, so every `.at()`, every string temporary and every destructor contributes an
  untakeable arc to the denominator. Block coverage is the more honest single number, and the
  useful artefact is the *list* of uncovered lines, not the ratio.

**Superseded for two rows since (Sep 8 2026).** 1b took `Insertion::iterate` to 100 % lines,
branches and blocks and `Insertion::initialize_event` back to 100 % lines; 2a/2b took
`Dinucl_markov::iterate` from 0 % to 100 % lines and blocks. The table above is kept as the
*pre-migration baseline* — it is what the delivered figures are measured against — but the two open
rows are now only `Gene_choice::iterate` (59.2 % blocks) and `Deletion::iterate` (**still 0 %**).
4a is the step that closes the second, and it is now scheduled before S5 (§6.8, F4).

`Gene_choice::iterate` at 59.2% blocks after T0 is the figure to watch, and **raising it is part of
5a's definition of done**. The uncovered remainder is mostly the `no_d_align` exhaustive path, which
is the hardest single piece in the set *and* the one the regression corpus never reaches — §7.9 is
what that combination already cost once. T0 covered the path's entry conditions; 5a has to reach the
position enumeration itself, on a fixture that forces it. Until then that body has no gate at all
below the convergence tests, which are unseeded and can only catch a crash.

### 6.5 — Delivered (1b): B6 *(Sep 7 2026)*

The three-arm `if/else` chain is one body. **All of 1a's sections pass unchanged**; one *fixture*
was adapted — the VJ sections now build their state on a VJ-ordered registry, because a junction
resolves its neighbours from the ordering and `VJ_ins_seq` has none in a VDJ one. No assertion
changed. Full ladder green, regression bitwise.

**Coverage: `Insertion::iterate` went to 100% lines, branches and blocks** (from 97.6 / 65.4 / 71.9).
The collapse deleted exactly the arms nothing could reach — the three-way seq_type chain and the
`throw invalid_argument` backstop that `initialize_event()` already made unreachable. That is the
clearest statement of what a branch collapse buys: the untestable code is gone rather than covered.

Three things beyond the chain:

- **The neighbour ids are resolved once at `initialize_event()`**, not per scenario — G10's tier 1.
  Three `std::string` comparisons per scenario become two member reads. It adds one branch B6 owes:
  an insertion at the end of the ordering has nothing to bound its junction, rejected at
  initialization rather than surfacing as a `kNoSeqType` subscript in the hot loop. Covered by its
  own section, which brought `initialize_event` back to 100% lines.
- **`initialize_crude_scenario_proba_bound`'s `switch` is one `events_map` lookup.** The switch
  converted `ins_seq_type` back into the very string the map is keyed by.
- **`initialize_event` addresses `seq_type_id` directly**, instead of converting the name back to a
  `Seq_type` through a three-way string chain. The scenario maps are `SeqTypeId`-keyed and the id is
  the only identity a non-legacy seq_type can have, so the conversion was both redundant and a hard
  ceiling: a tandem-D `D1D2_ins` could not have initialized at all. The chain that remains feeds
  genuinely `Seq_type`-keyed APIs — the generation path's `unordered_map<Seq_type, string>` (B9) and
  the Len_proba machinery (G5/S4) — and now goes through a helper that *returns* the value instead
  of an out-parameter seeded with `VD_ins_seq`, an idiom that reads as "defaults to the VD junction"
  when the value is in fact dead on every path that does not throw.
- **The name and the id are checked to agree**, once, at `initialize_event()`. This replaces a
  weaker guard: rejecting an *unrecognised* name caught less and did it by enumerating the three
  legacy junctions. A disagreement between the two identities does not fail — everything downstream
  is keyed by the id, so it aliases the event's segment, offsets and bounds onto another seq_type's
  keys and returns a wrong answer. That is the shape of the trap B2 hit with VJ, and the constructor
  fix above is the same class of bug caught one layer down.
- **The `Insertion(Seq_type)` constructor now sets `seq_type`.** It previously left the string empty
  while `ins_seq_type` held the enum — two identities of the same fact, disagreeing until a caller
  happened to set one. Harmless while every lookup went through the enum; a latent trap the moment
  one goes through the name, which is how `test_EventUtils.cpp`'s bridge test caught it.

The three `[!shouldfail]` defects are **not** fixed here: adding offsets and a mismatch list is a
behaviour change, and §1 keeps fixes out of refactoring commits. They land next, in their own
commit, with their own regression run.

**Harness change worth knowing about**: `IterateTestState` now takes a `SeqTypeRegistry`, defaulting
to a VDJ ordering, with `vj_seq_type_registry()` for a model with no D. The harness previously built
every map from `legacy_seq_type_registry()`, which registers the six seq_types but sets **no
ordering** — so `left_neighbor()` answered `kNoSeqType` for everything and *no generic body could be
tested at all*. Steps 2a and 4a would have hit this too.

### 6.6 — Delivered (2a): the `Dinucl_markov` characterization *(Sep 7 2026)*

`tst/igor/Core/test_dinucl_markov_iterate.cpp`, **75 assertions in 12 `TEST_CASE`s**, against the
unmodified event. Coverage of `Dinucl_markov::iterate` went from **0% to 98% lines / 69.8% blocks** —
it had no unit test of any kind before this.

What the sections pin, beyond the matrix rows:

- **The seed is the anchor's own nucleotide; every later pair is read-to-read.** The first term of
  the product is `(anchor's last nt, first read nt)` and the rest are `(read[i-1], read[i])` — the
  anchor sequence is used *once*, and nowhere else. A body that seeded from the read instead passes
  every other assertion, so it gets a dedicated section that changes only the anchor's last base.
- **The reverse traversal's double reversal.** DJ reverses the read window, fills, and reverses the
  result back; asserting only the filled content would pass with one of the two reversals dropped,
  so the section asserts the final orientation *and* that it is not the singly-reversed one.
- **The two placeholder guards are separate.** The first position has its own copy, because its
  `previous` nucleotide comes from the anchor. A mutation deleting that copy survived the whole
  suite until a section pre-filled position 0 specifically — the second guard's section had been
  overwriting position 1.
- **The junction is mutated in place, through the pointer the `Insertion` stored** — a confirmed
  defect, §7.13, pinned by a `[!shouldfail]` case rather than asserted.
- **An ambiguous read position is averaged, not indexed** — `dinuc_proba_matrix` instead of a
  marginal lookup, and `-1` in the realization indices. The matrix is built by
  `update_event_internal_probas()`, which `GenModel` calls and `initialize_event()` does not, so the
  fixture builds it explicitly. **Follow-up (agreed, not scheduled here):** that construction belongs
  in `initialize_event()`. It is a tier-0 computation in G10's terms — it depends only on the
  marginals — and leaving it to an out-of-band call means the matrix is whatever the last caller
  left, or zero. Any consumer that forgets it gets probability 0 for every ambiguous position, with
  no diagnostic. Moving it is a behaviour-preserving change for `GenModel` and removes a way to hold
  the event wrong.

**Ten mutations run, all caught** after the first-position gap was closed: both window offsets, the
final reversal, the seed source, both placeholder guards, the conditional index losing its
`previous` nucleotide, the scenario probability, pruning, and the downstream bound.

The one line left uncovered is `iterate()`'s `throw invalid_argument` for an unknown seq_type —
unreachable, because `initialize_event()` performs the same check first. Exactly the situation
`Insertion` was in before B6, and B7 deletes it the same way.

Harness addition: `preset_placeholders()`, which writes a junction as an `Insertion` leaves it —
placeholder nucleotides and **no offsets**. The missing offsets are deliberate: they are the defect
§6.3 pins, and a fixture that supplied them would test against a scenario the production code cannot
produce.

### 6.7 — Delivered (2b): B7 *(Sep 8 2026)*

The hardcoded traversal-spec table, the two residual `switch`es in `iterate()`, the one in
`add_to_marginals()`, three raw `new int[]` buffers and nine named members are gone.
**All of 2a's sections pass unchanged.** `Dinucl_markov::iterate` reaches **100% lines and blocks**.

The load-bearing decision: **the registry supplies *which* segment is adjacent, the event supplies
*which side* seeds the chain.** They are different kinds of fact — the ordering is topology, the
direction is a property of the Markov chain — and the old table conflated them. `event_side`
already carried the direction; `Model_Parms` derives it from the gene class for legacy files and
reads it from the file for v2 ones, so B7 changes no model format. A test builds a VD junction
seeded from its *right* neighbour, which IGoR ships no model for, to show the anchor is not being
recovered from the name.

Three supporting changes:

- **`Model_Parms::finalize()` tells every event what sits next to it** — `Rec_Event` gains two ids
  and a non-virtual `set_adjacent_segments()`. Topology is the model's fact: the ordering is read in
  one place, at one moment, and events are *told* rather than handed a registry to ask.
  `Dinucl_markov` derives its junction from those two ids plus its own `event_side` and stores
  nothing, so there is no resolution step to run twice, no state to keep in sync and nothing to
  make idempotent. `Insertion` uses the same two ids instead of resolving its own.

  *(This replaces a first attempt that put a `resolve_topology()` virtual on `Rec_Event`, called
  from `finalize()` and again from `initialize_event()`. Quentin's objection was right on three
  counts: it made a model-level fact the event's responsibility, the second call was compensation
  for un-finalized models rather than design, and it sat among inference-only lifecycle methods
  while also serving generation.)*
- **Finalization is no longer optional.** `Model_Parms` tracks whether an event has been added since
  the last `finalize()` and asserts on it where the model is handed out, and the copy constructor
  re-finalizes — a deep copy rebuilds events through `copy()`, which does not carry resolved state,
  and inference makes one copy per thread. Fixtures that build an `events_map` without a
  `Model_Parms` now do the adjacency pass themselves, which is the honest place for it: the harness
  is what stands in for `finalize()`.
- **Per-spec scratch state.** The index buffer and memory layer live in the spec. A model with
  several junctions per event needs no new members, and the buffer is freed with the event rather
  than by three hand-written `delete[]`s. It is `clear()`ed and `push_back`-filled per scenario with
  its capacity reserved once, rather than written at fixed offsets alongside a separate filled
  count — which also removes a latent fault the fixed-width array carried: a position that arrives
  already filled contributes no probability term, but its *previous* scenario's index stayed in the
  array and was still credited by `add_to_marginals()`. Unreachable today, because `Insertion`
  re-allocates the junction as placeholders on every call; reachable the moment §7.13's shared
  buffer gains a sibling.
- **`iterate()`'s empty-specs `throw` is deleted** — unreachable, since `initialize_event()` refuses
  to leave them empty. The same dead backstop B6 removed from `Insertion`.

**Deferred, deliberately**: G9's `first_occupied_*` walk. It differs from the ordering neighbour
exactly when the anchor is empty, which is §7.12's crash, so swapping it in *is* that fix rather
than a refactor (§7.11). **The behaviour is decided (Sep 8 2026): throw.** Silently discarding a
scenario whose anchor was fully deleted would bias the model — those scenarios are geometrically
legitimate and their absence would not be visible anywhere — whereas a scenario that cannot be
scored is a modelling error the user must see. It still changes behaviour relative to today's
segfault, so it lands in its own commit with its own regression run, not in a refactoring step —
scheduled as R2 in §6.9.

A tandem-D `D1D2_ins_seq` junction, which no `Seq_type` enum names, resolves correctly today; the
test checks it through `resolve_topology()` rather than `iterate()`, because the harness cannot yet
build events keyed by `SeqTypeId` alone.

**Two fixture changes, no assertion changes.** `make_dinucl_markov()` now sets the chain direction
the model file carries — the fixture had been under-specifying the event and relying on the
hardcoded table to supply it. And `IterateTestState::add_event()` now keys events exactly as
`Model_Parms::get_events_map()` does, with `Undefined_side` for `Dinucl_markov`; the harness had
been keying by the event's own side, which only worked while that side was always `Undefined_side`.

### 6.8 — Re-assessment after 2b *(Sep 9 2026)*

Two of four events are collapsed. By mass that is not half:

| Event | `iterate()` lines | Status |
|---|---:|---|
| `Insertion` | ~110 | ✅ 1b |
| `Dinucl_markov` | ~180 | ✅ 2b |
| `Gene_choice` | **845** ([Genechoice.cpp:185-1029](../src/igor/Core/Genechoice.cpp#L185)) | steps 3 and 5b |
| `Deletion` | **1004** ([Deletion.cpp:239-1242](../src/igor/Core/Deletion.cpp#L239)) | step 4b |

**~77 % of the mass is ahead, in two of the remaining rows.** Smallest-first was the right order and
the two delivered steps validate the method; they do not yet validate the services.

**F1 — A0, S2 and S3 have no production consumer.** `JunctionGeometry.h` is referenced by its own
test and by one doc comment in [Rec_Event.h:113](../src/igor/Core/Rec_Event.h#L113); the four A0
capability virtuals are called from exactly one place, `PendingModifierBounds::rebuild()`, itself
uncalled. §1 allows a service to sit without a caller for one commit; this is three steps.

The risk is specific. B6 and B7 consumed *none* of it — G9 came out instead as
`Rec_Event::set_adjacent_segments()`, a model-tells-event mechanism invented during 2b and not in
this plan. So G1/G2/G3 remain the only §2 patterns no production code exercises, and S2/S3 were
mutation-verified against a specification rather than against a consumer. If B11a finds
`check_overlap()`'s shape wrong, S3's tests get rewritten *inside* a refactoring diff — the one
signal §6.1 exists to protect. The response is to get consumers onto them, not to build more
services: hence steps 3 and 4b are labelled with which service they are the first consumer of, and
S5 no longer lands before the characterization that would exercise it.

**F3 — S4 was under-scoped, and it is a tandem-D blocker.** As originally written S4 deleted three
`has_effect_on` overrides, which is cheap and unblocks nothing. G5's actual weight is elsewhere:
`vd_`/`vj_`/`dj_length_best_proba_map` are duplicated as six members across
[Genechoice.h:195-197](../src/igor/Core/Genechoice.h#L195) and
[Deletion.h:214-216](../src/igor/Core/Deletion.h#L214), plus `vj_length_d_position_proba`, and the
whole machinery hangs off `iterate_initialize_Len_proba(Seq_type considered_junction, …)` — a pure
virtual on `Rec_Event`, keyed by the enum.

That last point is the finding that matters for the parent plan. `Insertion` reaches it through
`insertion_seq_type_or_throw()`, which throws on any name the `Seq_type` enum does not know, so **a
tandem-D `D1D2_ins` junction cannot pass `initialize_Len_proba_bound()` at all** — it throws before
inference starts. The parent plan's milestone-1 path (`B2 → B8 → B11 → B7 → B6 → B5`) does not
mention the Len_proba machinery, and is incomplete without it. S4 is therefore the pair-keyed
junction structure, on the critical path rather than beside it, and it is also C2's stated landing
point: one accessor, so the exact-match lookup becomes a range scan in one place instead of twenty.
Per §9 its query returns a *set* of bearing events, not the first match. **§6.10 works the
analysis through and proposes an S4a/S4b/S4c split.**

**F4 — 4a moves ahead of S5.** `Deletion::iterate` is at 0 % coverage and is the largest single
untackled job in the plan, bigger than 1a and 2a combined. Scheduling it after S5 would mean
characterizing a body whose safety mechanism S5 had just replaced. Running it first also gives S5 a
consumer instead of making it a third uncalled service, and its definition of done becomes "4a's
sections pass unchanged", the same shape every *b* commit already has.

**F5 — O6's switch is required; only its *fixture* rationale was overstated** *(corrected
Sep 10 2026)*. An earlier version of this finding claimed O6 was "dissolved" because the exhaustive
path is reachable from production code without a new switch. That conflated two things, and §2.6
already had it right.

*True*: `no_d_align` is not "the aligner returned nothing" — it is set `false` only when an
alignment **survives pruning** ([Genechoice.cpp:530](../src/igor/Core/Genechoice.cpp#L530)). The
path is therefore reachable two ways, an empty alignment list or a pruning threshold above every
alignment's bound, and T0's G6 sections already use the first. Step 5's fixture does not need to
defeat the aligner, and does not need the switch in order to *reach* the path.

*False*: that the switch is therefore unnecessary. **It is required by the collapse itself.** The
fallback lives inside `case D_gene`, so V and J never reach it — an accident of the switch
statement rather than a stated policy. Once the body is generic there is no `case D_gene` left to
confine it to, and the fallback then fires either for every gene class — a behaviour change, and
catastrophic for V and J at hundreds of templates × hundreds of positions — or for none, also a
behaviour change, losing D's fallback. The per-`Gene_choice` boolean is what carries the asymmetry
after the statement that used to carry it is gone. §2.6 specifies it; T0's sections *"V with no
alignments enumerates nothing"* and *"J with no alignments enumerates nothing"* already pin the
behaviour it must preserve, and that test file's own comment says exactly this.

Same migration pattern as B7's `event_side`: the hardcoded table held information that had to be
**relocated into data**, not deleted. Legacy models derive the flag from the gene class, so no
model-format change; whether v2 should carry it explicitly is open.

What genuinely remains of O6 as a *policy* question is whether V and J should ever fall back. The
default answers it conservatively either way.

### 6.9 — Phase R: the repair queue

Decided behaviour changes, none of which had a row in §6 before this re-assessment. They land
**after 5b**, in the order below, each in its own commit with its own regression run.

| # | Fix | Decided | Touches | Expected regression effect |
|---|---|---|---|---|
| R1 | §7.13 — give `Dinucl_markov` its own constructed-sequence layer | Sep 7 | `Dinucl_markov` | none expected; the write lands at a different layer, same content |
| R2 | §7.12 — `first_occupied_*` walk, **throw** on an empty anchor | Sep 8 | `Dinucl_markov` | none on the corpus (no model produces an empty anchor); removes the `[.]` tag from the reproducer |
| R3 | `Insertion` writes offsets and a mismatch list (three `[!shouldfail]`), **its `get_offset_role` stops reporting `None`**, and the leaf invariant's offsets half becomes assertable (see below) | Sep 7 | `Insertion` | none expected — neither `Insertion::iterate` nor `Dinucl_markov::iterate` touches `seq_offsets` at all |
| R4 | `dinuc_proba_matrix` construction moves into `initialize_event()` | Sep 7 | `Dinucl_markov` | none — `GenModel` already calls the out-of-band builder |
| R5 | §7.1 and §7.8 off-by-one corrections (decision O4) | Sep 1 | `Gene_choice` | **golden data moves**; the credited core length changes |
| R6 | Within-clique **joint** max in the span fold, using S4b's group hook; optionally cross-clique parent indexing after it | Sep 10 | the span fold (all events) | **golden data may move** — a tighter bound prunes more, so fewer scenarios are summed. Needs the **convergence** gate, not just regression |

R1–R4 are each expected to be bitwise-neutral despite being behaviour changes — they close paths
the corpus does not reach. That expectation is the thing to *test*, not to assume: a surprise here
is a finding about the corpus, not a reason to accept the diff. R5 and R6 are the two expected to
move numbers, and they are last for exactly that reason. They are independent of each other.

**R6 in more detail** *(Quentin, Sep 10 2026)*. §6.10 shows the span fold accumulates
`∏ₑ maxᵢ Pₑ(rₑ|i)`, a product of per-event maxima, and that taking the max **jointly** over a
conditioned clique — the D block, in every model IGoR ships — removes the dominant slack term at no
storage cost. It is nonetheless a **behaviour change**: a tighter upper bound prunes more, so
scenarios that were explored before are now discarded, and the summed marginals move by whatever
those scenarios contributed. With the probability-ratio threshold set low enough the shift should be
negligible, but "should be" is what the convergence gate exists to check — this is the one repair
whose effect is on *which scenarios are summed* rather than on a single value, so regression
bitwise-equality is the wrong question to ask of it alone.

Optional second stage, same commit or the next: **cross-clique parent indexing**, for `V → D` and
`J → D`, where the conditioning parent is already chosen at consumption time but unknown when the
profile is built. Estimated ≈ 270 kB per span per thread at human BCR-H dimensions, against the
memory S4c's deduplication frees (five builds → one). Strictly optional, and only worth doing if
R6's first stage shows the remaining slack still matters.

**Naming**: R1–R6 are *repairs*; F1–F5 in §6.8 are the *findings* of the re-assessment. Different
sequences, deliberately different letters.

#### R3 also has to correct A0, not just `iterate()` *(Quentin, Sep 10 2026)*

`Insertion::get_offset_role()` returns `OffsetRole::None`, and the comment beside
`get_offset_delta_bounds` justifies it: *"An insertion writes no offsets at all: its span is
derived from where its neighbours already sit, which is exactly what makes the generic B6 rule
possible."* Two claims are conflated there. **The span is derived from the neighbours** is true and
is B6's whole point. **Therefore no offsets are recorded** is the defect.

This matters more than the missing write on its own. A capability query is what a consumer is
supposed to ask *instead of* knowing which event produced a segment — that is A0's entire purpose —
so a query reporting the defect propagates exactly the coupling it exists to remove. A0 documented
the current behaviour as if it were the design.

R3 therefore has three parts, not two: write the offsets, write the mismatch list, and make
`get_offset_role` report the truth — plain `OffsetRole::Creates`, with **no enum change**; §2.5
records why the `Anchors`/`Derives` split first proposed here was the wrong answer.

A fourth part is available only once R3 lands: the **offsets half of the leaf invariant** (§2.5).
`iterate_wrap_up`'s debug check already rejects a leaf carrying unfilled *content*; the matching
check for unplaced *offsets* cannot be added while `Insertion` writes none, because it would fire
on every scenario.

Each fix removes a `[!shouldfail]` tag or a `[.]` tag as part of its definition of done. A fix that
lands without its tag coming off leaves a case that now passes under `[!shouldfail]`, which Catch2
reports as a failure — the mechanism §6.1 chose deliberately, and the reason the queue can be
deferred this far without being forgotten.

### 6.10 — G5 analysis: what the Len_proba machinery is, and what S4 must do *(Sep 9 2026)*

Written to settle S4's scope (finding F3, §6.8). Read: all four `iterate_initialize_Len_proba` /
`initialize_Len_proba_bound` pairs, the base traversal, and every consumer.

#### The four bodies are one body modulo one scalar

`Gene_choice`, `Deletion` and `Insertion` are structurally identical:

```cpp
if (has_effect_on(J)) {
    base_index_map.set_current_layer(event_index, 0);
    base_index = base_index_map.get(event_index);
    for (auto& r : event_realizations) {
        real_max_proba = maxᵢ marginals[base_index + r.index + i*size()];
        wrap_up(..., scenario_proba * real_max_proba, seq_len + Δ(r));
    }
} else { wrap_up(..., scenario_proba, seq_len); }
```

| Event | `Δ(r)` | Enumerates? | Extra |
|---|---|---|---|
| `Gene_choice` | `+ r.value_str.length()` | yes | — |
| `Deletion` | `− r.value_int` | yes | — |
| `Insertion` | `+ r.value_int` | yes | writes `constructed_sequences` as a **side channel** |
| `Dinucl_markov` | **0** | **no** | multiplies `p^L`, reading `L` back out of that side channel |

`Dinucl_markov` differs in kind: it contributes probability as a function of an already-accumulated
length, and learns that length only because `Insertion` stashed a dummy string in
`constructed_sequences`. That is the `//TODO constructed sequences should not be used but it is
useful to compute the dinucl contribution` on
[Rec_Event.cpp:374](../src/igor/Core/Rec_Event.cpp#L374).

#### What `has_effect_on` means

Every implementation answers exactly one question: *does a realization of this event change the
accumulated length of the span named by this argument?* **Length only** — not offsets, not content,
not probability. `Dinucl_markov::has_effect_on` returning `true` while contributing zero length, to
gate a *probability* factor, is the one place the name actively misleads.

And the argument is a **span**, not a seq_type: `VJ_ins_seq` means the VJ insertion segment in a VJ
model and the whole V→J span in a VDJ one. That is §2.5's "hidden generalisation", confirmed at
every call site.

#### Five incidental findings

1. **The queue-level filter is commented out.** [Rec_Event.cpp:387-397](../src/igor/Core/Rec_Event.cpp#L387)
   carries `//if(next_event_p->has_effect_on(considered_junction)){` with
   `//TODO fix this and find a way not to loop over all events`. The traversal therefore visits
   every event in the model and each self-filters at the top of its own override — the predicate
   exists twice over, at the wrong level.
2. **`Gene_choice::has_effect_on` never returns `true` into a value anyone reads.** It is `true`
   only for `D_gene` on `VJ_ins_seq`. Every consumer is guarded `if (d_chosen) {adjacent}
   else if (other_chosen) {vj}` — [Genechoice.cpp:321](../src/igor/Core/Genechoice.cpp#L321),
   [:975](../src/igor/Core/Genechoice.cpp#L975),
   [Deletion.cpp:453](../src/igor/Core/Deletion.cpp#L453) — so the VJ map is read only when there
   is no D, and when there is no D there is no D `Gene_choice` to fire. Dead in both topologies.
3. **`Deletion` builds `vj_length_best_proba_map` in every VDJ model and nothing reads it.**
   `get_deletion_effective_junctions(V_gene_seq, ·)` returns `{VD, VJ}` unconditionally. Wasted
   initialization; also what hides the asymmetry noted in §2.5.
4. **The VD span profile is built five times per model, per thread** — `Gene_choice(V)`,
   `Gene_choice(D)`, `Deletion(V,3')`, `Deletion(D,5')`, and `Insertion(VD)`'s own
   `junction_length_best_proba_map`. Five identical traversals, five stored copies. Same for DJ.
   Six named members hold three logical maps, across two classes.
5. **`Insertion::initialize_Len_proba_bound` runs the whole traversal `|R|` times where once would
   do.** [Insertion.cpp:549-557](../src/igor/Core/Insertion.cpp#L549) loops over its own
   realizations *outside* the traversal purely to set `inserted_str` for the Dinucl side channel,
   but `Insertion::iterate_initialize_Len_proba` re-enumerates the same realizations *inside* and
   overwrites it; `wrap_up` takes `model_queue` by value so the queue survives each pass. With ~40
   realizations that is a 40× init cost producing an identical map. **High confidence, confirm by
   measurement in S4b** before removing.

`chosen` is a model-level fact, not a per-scenario one — `EventUtils::check_gene_choice(…,
processed_events)` at [Genechoice.cpp:1111](../src/igor/Core/Genechoice.cpp#L1111) resolves it once
during priority-ordered `initialize_event`. Each gene's anchor pair is therefore fixed at init,
which is what makes a cache keyable at all.

#### API naming — reserve the family now

The ambiguity in `has_effect_on` is that three orthogonal axes get one name. A0 already occupies
four cells; naming the axes shows which are empty.

| | addressed by a **segment** | addressed by a **span** |
|---|---|---|
| offset | `get_offset_role` ✅ A0 · `get_offset_delta_bounds` ✅ A0 | — |
| length (bounds) | `get_length_contribution` ✅ A0 | **`affects_length_of(SegmentSpan)`** ← replaces `has_effect_on` |
| length (per realization) | **`length_delta(const Event_realization&)`** ← missing, and it is the hook | — |
| sequence content | `get_seq_construction_role` ✅ A0 | — |
| probability | — | **`span_proba_factor(SegmentSpan, const SpanAccumulator&)`** ← `Dinucl_markov`'s cell |

```cpp
/// An ordered, anchor-exclusive range of the registry ordering. The addressing unit for
/// junction-length bounds (G5), for effect queries, and for Phase D's cluster boundaries
/// (D.3) -- the same object in all three, named once.
struct SegmentSpan { SeqTypeId left, right; };
```

`affects_length_of` is deliberately narrow rather than a general `affects()`: a general predicate
would need a "what" argument and would immediately be a worse `has_effect_on`. The narrow name
leaves `affects_offsets_of` / `affects_content_of` free, and D.9's `get_context_seq_types()`,
`get_context_dependency()`, `is_branching()`, `is_multi_realization()` stay uncollided.

#### Would more Phase A capabilities help?

Yes, and specifically two — both already implied by the table above:

- **A per-realization length delta.** A0 gives `LengthContribution{min,max}`, the bound over all
  realizations, which is what `PendingModifierBounds` needs. The traversal needs the value for
  *one* realization. This single accessor is what collapses four bodies into one.
- **A probability hook for non-enumerating events.** Phase D already names the shape —
  D.9's `is_multi_realization()`, "contributes a sum of paths to a single subscenario (e.g.
  `Dinucl_markov`)". The concept is reserved; G5 is where it first has a caller.

**A modifier-type enum would not help.** The four `Δ` implementations differ by sign and by which
field of `Event_realization` they read (`value_str.size()` vs `value_int`); a per-realization
accessor states that directly. A taxonomy tagging events insertion-like / deletion-like would
re-create `Gene_class` one level up — the mistake B1 has just finished undoing. The *qualitative*
taxonomy already exists as `SeqConstructionRole{None, Creates, Modifies, Fills}`.

#### Bound tightness: where the relaxation costs pruning

`maxᵢ marginals[base + r.index + i·size()]` is applied **per event** and the results multiplied, so
the accumulated bound is `∏ₑ maxᵢ Pₑ(rₑ|i)` against a truth of `max_decomp ∏ₑ Pₑ(rₑ|actual)`.
Product-of-maxes ≥ max-of-products, and the slack **compounds with the number of marginalised
events**.

That matters because of what gets marginalised. Every VDJ model IGoR ships makes the D block the
densest conditioning cluster in the graph:

```
%GeneChoice_V_gene…;GeneChoice_D_gene…            (human TRB, human BCR-H)
%GeneChoice_J_gene…;GeneChoice_D_gene…            (human TRB, human BCR-H, mouse TRB)
%GeneChoice_D_gene…;Deletion_D_gene_Three_prime…
%GeneChoice_D_gene…;Deletion_D_gene_Five_prime…
%Deletion_D_gene_Five_prime…;Deletion_D_gene_Three_prime…   (mouse TRB)
```

— and that block is exactly what `span(D1,J)` marginalises when D1 enumerates.

**Within-clique conditioning has a free fix.** For the D block the conditioning parent is *inside*
the marginalised span, so a **joint** max over `(gene, del5, del3)` is available at build time: no
extra storage, ~68k evaluations once per model at human BCR-H dimensions (35 × 44 × 44). Today's
structure cannot express it, because `real_max_proba` is computed inside each event's own body.
**This is why the S4b hook must fold over contribution *groups* rather than events** — a group
being one event today and a conditioned clique tomorrow. Reserving that seam is not insurance; it
is the mitigation.

**Cross-clique conditioning is live but priced.** `V → D` and `J → D` differ: at *consumption* time
V and J are already chosen, but at *build* time — `initialize_Len_proba_bound` runs once before any
read — they are not, so the bound discards information the scenario has. Indexing the span profile
by conditioning context costs `|V| × |J|` copies: ≈ 97 × 7 × 50 lengths × 8 B ≈ **270 kB per span
per thread** at human BCR-H dimensions *(estimate, worth checking)*. Not obviously prohibitive —
and S4c's deduplication (finding 4) frees roughly what it would spend, so the two belong on the
same page rather than being decided apart.

**Measure before milestone 2 commits.** `vj_length_d_position_proba` already marginalises one gene,
two gaps and two Dinucl — the same mechanism at roughly half the event count, on the existing
corpus, today. Instrumenting `bound / realized_proba` at leaves gives a measured baseline instead of
an argument. **Scheduled into 5a**, which already owns raising coverage on that path and already
needs a fixture that forces it.

#### Proposed S4 split *(not yet approved)*

**Approved Sep 10 2026**, with the joint-max tightening moved out to phase R.

| | Content | Justified by |
|---|---|---|
| **S4a** | `SegmentSpan`; `affects_length_of(SegmentSpan)` replacing `has_effect_on`; the filter moved to the queue level (closing finding 1) and removed from the four bodies. Also the tier-3 harness check (§2.5) | naming; no behaviour change |
| **S4b** | `length_delta` + `span_proba_factor`, both as **group** hooks; the four bodies → one non-virtual traversal; `SpanAccumulator` (carrying per-segment lengths) replaces the `constructed_sequences` side channel and its `Seq_type_str_p_map` parameter | the collapse; findings 1 and 5 |
| **S4c** | Span-keyed structure owned by the model; `⊗ᵐᵃˣ`; six members → one; each span built once; `initialize_Len_proba_bound` de-virtualised; **the dead code of findings 2 and 3 deleted here** | **removes the tandem-D enum ceiling**; findings 2, 3, 4 |
| **S4d** | Tensor-backed containers, **gated on the Tensor API landing** — see the container note in §2.5 | performance only; strictly optional |
| **→ 5b** | `⊗ᵉⁿᵘᵐ` — bucketing `(realization, left_len)` pairs by total, sorted | `no_d_align` only |
| **→ R6** | The within-clique **joint** max, and optionally cross-clique parent indexing | bound tightening — **changes results**, §6.9 |

**S4b keeps the group-shaped hook but not the joint max.** A group of one is exactly today's
behaviour — each event's own `maxᵢ`, multiplied in the same order — so S4b stays bitwise. Taking
the max *jointly* over a conditioned clique tightens the bound and therefore changes which
scenarios survive pruning, which is a behaviour change and belongs in phase R (R6). Deferring the
tightening is not the same as deferring the seam: retrofitting a per-event hook into a per-group
one touches every implementation, so the shape lands in S4b and R6 merely uses it.

**Finding 2 and 3's dead code is deleted in S4c, not deferred to phase R** *(decided Sep 10 2026)*.
Provably-unread code is not a behaviour change in any observable sense, and porting it into the new
structure only to remove it later is worse than deleting it at the point of the move. The claim
"provably unread" is what the S4c regression run has to bear out: if the analysis in findings 2–3 is
wrong, the run is no longer bitwise, which is exactly the signal wanted.

`⊗ᵐᵃˣ` belongs in S4c, not 5b: every gene needs it to build `span(G,B)` past marginalised genes.
Only the *retention* of the decomposition is 5b's.

`⊗` must accept a **scalar weight per operand**, not just profiles — see the B10 note in the parent
plan. Cheap now, and it is what lets an explicitly-weighted absence branch tighten the bound.


### 6.2 — The regression gate has a flaky output

*(Observed Sep 2 2026 during A0.)*

`pixi run test_regression` is designated the bitwise gate for every step of this plan. One of
its outputs is **not deterministic**: `best_scenarios_counts.csv` mismatched once on an A0 build
that adds only uncalled virtuals, then matched on the next three runs (two `test_inference`, one
full `test_regression`), while the pre-A0 baseline matched on its single run.

**Mechanism (diagnosed properly since; the row-ordering explanation this section used to give was
wrong — the comparator sorts by sequence index before comparing).** The marginals themselves are
non-deterministic: `#pragma omp for schedule(dynamic) nowait` over query sequences, per-thread
`single_thread_marginals +=`, merged under `#pragma omp critical`. Floating-point addition is not
associative, so the sums differ by a few ULP between runs — invisible in the six-digit
`iteration_N.txt` dumps, which is why only this one output moves.

Those ULP differences decide a strict `>` in `Best_scenarios_counter` between realizations that are
*exactly* degenerate: TRBV3-1\*01, TRBV3-2\*01 and TRBV3-2\*02 have identical
`P(V=g) · P(v_3_del=d | V=g)` for d ≤ 7 under a uniform initialization with indistinguishable
alignments. It only surfaces when that degenerate group straddles the `output.scenarios` cutoff of
10; with 15 scenarios all three are kept and the comparison passes.

**Discriminating a flake from a regression**: re-run with `OMP_NUM_THREADS=1`. One thread is
deterministic, so a single-threaded pass plus a multi-threaded failure confined to
`best_scenarios_counts.csv` is the flake. This is what 2b's first regression run turned out to be.

**This matters more than its size suggests.** A gate that fails intermittently trains its readers
to re-run rather than investigate, which is exactly how a real regression gets waved through
during a multi-step refactor. Every "bitwise" claim in §6 depends on this gate meaning what it
says.

Fix shape: either make the reduction deterministic (fixed chunks accumulated in index order) or
give the counter a canonical tie-break (relative epsilon, then realization-vector ordering).
Pinning the thread count alone is not enough with `schedule(dynamic)`. Sorting rows in the
comparator — this section's original suggestion — would **not** help, since the rows differ in
content, not order.

## 7. Where a generic rewrite would silently change results

These are the traps. Each must be preserved bit-for-bit in the step that touches it, then fixed
in a separate, explicitly-labelled commit — never folded into a refactor.

### 7.1 — The credited match length in `Gene_choice` V and J (a real bound bug)

`Gene_choice` passes `n_error_free` to `get_err_rate_upper_bound()` as:

| Branch | code | correct core length |
|---|---|---|
| V [:346](../src/igor/Core/Genechoice.cpp#L346) | `gene_seq.size() - v_3_max_del - endo` | `gene_seq.size() + v_3_max_del - endo` |
| J [:979](../src/igor/Core/Genechoice.cpp#L979) | `gene_seq.size() - j_5_max_del - endo` | `gene_seq.size() + j_5_max_del - endo` |
| D [:513](../src/igor/Core/Genechoice.cpp#L513) | `(d_3_off + d_3_max_del) - (d_5_off - d_5_max_del) - endo` | same, minus 1 |

`*_max_del` is `get_len_min()`, i.e. **negative** (`-max_del_value`). The D branch therefore has
the sign right (off by one on the inclusive count); V and J have it **inverted**, crediting
`size + max_del` error-free positions where at most `size − max_del` can survive.

**Measured** (T0, `[endogenous]` sections), with rate 0.1 and one endogenous mismatch:

| Branch | credited exponent | correct | direction |
|---|---:|---:|---|
| V, 12 nt template, max 4 deletions | 15 | 7 | bound too **small** → over-prunes |
| J, 8 nt template, max 4 deletions | 11 | 3 | bound too **small** → over-prunes |
| D, 8 nt template, max 2 per end | 2 | 3 | bound too **large** → under-prunes, harmless |

So it is two defects, not one: a sign inversion in V and J, and an inclusive-count off-by-one in
D. Only the first can discard scenarios that should have been kept.

Since `get_err_rate_upper_bound(i,j) = (r/3)^i · (1−r)^j`
([Singleerrorrate.cpp:91](../src/igor/Core/Singleerrorrate.cpp#L91)) is strictly decreasing in
`j`, over-crediting yields a bound that is **too small** — i.e. *more aggressive pruning than the
model justifies*. At `r = 0.01` and `max_del = 16` the V and J downstream bounds are low by a
factor `0.99^32 ≈ 0.73`.

This is not a refactoring artefact; it is present on `develop` today. The generic G8 form
computes the core correctly and would therefore **change Pgen values** on the regression corpus.

**Handling (decided, O4)**: reproduce first, fix at the very end. Step 3 must carry the current
arithmetic verbatim — **parameterise the sign, do not derive it** — and every step through 5 keeps
it, so the whole migration stays regression-testable against the existing corpus.

**How it is tested (revised after review)**: the T0 sections assert the **correct** value and
carry Catch2's `[!shouldfail]` tag, rather than pinning the buggy value. Catch2 reports an
expected failure as a pass, so the suite stays green, and the moment the defect is fixed the case
starts passing — which `[!shouldfail]` turns into a failure, forcing the tag to be removed
deliberately. A test that states the intended behaviour also reads as documentation, which a
pinned wrong number does not. Each such case is section-free, because `[!shouldfail]` is evaluated
per test-case run and Catch2 re-runs a case once per leaf section. Expected
direction of change: slightly *less* pruning, slightly slower, marginally more complete scenario
sets. It is a modelling-visible change, so it ships with its own convergence evidence and its own
commit message — never folded into a refactor.

### 7.2 — The in-between gap `G`

Every current check uses `G = 0` implicitly. G2 makes `G` explicit; computing it from `len_min`
gives `0` for the current corpus (insertions allow 0, D allows full deletion), so it is a no-op —
**but only if computed from `len_min`, not from a tighter estimate.** Do not "improve" it in the
same step.

### 7.3 — `d_del_opposite_side_processed`

Two `//THIS IS A TEMPORARY FIX //FIXME` guards
([Deletion.cpp:596](../src/igor/Core/Deletion.cpp#L596),
[:818](../src/igor/Core/Deletion.cpp#L818)) plus the two endogenous-mismatch branches that key
off it. They interact with the order of the two D deletion events, and the `else` branch has a
commented-out body plus a `//TODO finish this part`. A generic rewrite that "cleans this up"
changes results. Carry it forward literally.

### 7.4 — `len_min` / `len_max` accumulation is order-dependent (latent)

Both file-loading constructors — `Deletion(seq_type, side, realizations)`
([Deletion.cpp:175-189](../src/igor/Core/Deletion.cpp#L175-L189)) and
`Insertion(seq_type, realizations)`
([Insertion.cpp:97-110](../src/igor/Core/Insertion.cpp#L97-L110)) — accumulate the two bounds
with an `if / else if` over an **unordered** map:

```cpp
if (v > (-len_min))      len_min = -v;
else if (v < (-len_max)) len_max = -v;      // only reached when v is not a new maximum
```

Because the second test sits behind `else`, a strictly ascending iteration order never reaches
it and leaves `len_max` at its `INT16_MIN` sentinel. The range constructors
(`Deletion(…, pair<int,int>)`, `Insertion(…, pair<int,int>)`) pre-seed both bounds and are
unaffected; the constructors the model reader actually uses do not.

**Measured, not assumed**: replaying the exact accumulation over `std::unordered_map<std::string,…>`
for the realistic key sets (deletions −4..16 and 0..16, insertions 0..30, and the 2- and
3-realization minimal-model cases) gives the correct bounds in every case on this libstdc++ —
the hash order places a near-maximal key early enough. **The bug is latent, not active.**

It is worth recording anyway, because the failure mode is silent and severe: `get_len_max()` is
read as `*_min_del`, so a sentinel value would push `v_3_max_offset` far negative and mark every
overlap check "safe", skipping real constraints. A0 removes the whole class of hazard by
computing the bounds once, correctly, in one place. Add a unit test asserting the bounds against
`std::minmax_element` over the realization set as part of A0.

### 7.6 — In `Gene_choice`, the overlap early-out is subsumed by the junction-length guard

*(Found Sep 1 2026 while writing T0's overlap sections; the first draft of them passed for the
wrong reason.)*

With D already chosen, `Gene_choice`'s V branch can drop a realization in two places:

```cpp
// (a) early, Genechoice.cpp:265
if ((v_3_off + v_3_max_del) >= d_5_max_offset) continue;
// (b) later, Genechoice.cpp:322
if (vd_length_best_proba_map.count(d_offset - v_3_off - 1) <= 0) continue;
```

**(a)** reads *"even with V deleted as far as it can go and D deleted as far as it can go, V's 3'
end still reaches D's 5' end."* **(b)** reads *"this V–D gap is not a gap the model can produce."*

They are not independent. Write `L = d_5_off − v_3_off − 1` for the gap, and substitute
`v_3_max_del = −max_del_V`, `d_5_max_del = −max_del_D5`:

```
(a) fires  ⟺  v_3_off − max_del_V ≥ d_5_off + max_del_D5
           ⟺  L ≤ −max_del_V − max_del_D5 − 1
```

The map in (b) is keyed by achievable gaps, and a gap is achievable when
`L = ins − del_V − del_D5` for some legal triple — so its **smallest key is
`min_ins − max_del_V − max_del_D5`**. (a) therefore fires only on gaps *strictly below the
smallest producible gap*, and such a gap is by definition not in the map. **(a) ⊆ (b).**

Verified by mutation on the T0 overlap sections: deleting (a), or deleting (b), each leaves every
section green; only removing both changes the outcome.

Three consequences:

1. **It is not dead code and must be kept.** (a) fires before `iterate_common()`, the mismatch
   scan and the error-rate bound, so it skips real work on the most-executed event.
2. **It is an optimization, not a constraint** — which is *reassuring* for B11. If the generic
   `check()` shifts the Infeasible boundary slightly, results cannot change, because (b) still
   catches whatever (a) missed. The correctness burden sits entirely on the junction guard.
3. **`Deletion` must not inherit this conclusion.** At `Deletion::iterate` time the moving end's
   interval has already collapsed to a point (§2.2), so the two conditions no longer line up the
   same way. Step 4 re-runs this experiment rather than assuming.

### 7.7 — `Insertion`'s missing safety check

`Insertion` performs **no** overlap check at all; it computes a length from two offsets and
discards the scenario if that length is outside the realization range
([Insertion.cpp:241-247](../src/igor/Core/Insertion.cpp#L241-L247)). The user is right that this
is a shortcut. Under the generic form it is *tempting* to add the check — do not, in B6. A
negative computed insertion count is currently caught by the range lookup returning 0; adding an
explicit geometric check would change which scenarios are enumerated even if it changes no final
probability. If it is worth adding, it is a separate step after B5 with its own regression run.

---

### 7.8 — The two D realization paths disagree on what a junction length means

*(Found Sep 1 2026 during T0.)*

```
alignment path, Genechoice.cpp:322 :  L = d_5_off − v_3_off − 1   ⟹  d_5_off = v_3_off + L + 1
position  path, Genechoice.cpp:556 :  d_5_off = v_3_off + L
```

`L` in the position path comes from `vj_length_d_position_proba`, which is composed from the
**same** `vd_length_best_proba_map` the alignment path's guard consults — so the two are using one
quantity under two different conventions, and the position path is one short. At `L = 0`, meaning
"no VD insertions", it places D's 5' end *on* V's 3' end rather than immediately after it,
overlapping V's last nucleotide.

Reachable only through `no_d_align`, i.e. only for a D gene the aligner found nothing for. Pinned
by a `[!shouldfail]` case in the T0 suite.

**This is for the maintainer to adjudicate, not for the refactor to decide.** It changes inference
results on any model where the exhaustive path fires, so it is a modelling-visible fix and belongs
with §7.1 at the end, not inside B11. B11 must reproduce it verbatim and generalise it unchanged.

### 7.9 — The `no_d_align` path wrote no overlap verdict, and the B8 port turned that into a crash

*(Found Sep 2 2026 while verifying A0. Fixed in the same commit; recorded because the shape
recurs.)*

`Gene_choice`'s D branch writes `VD_safe` / `DJ_safe` in two places: the preamble, when the
neighbour has **not** been chosen, and inside the alignment loop under `vd_check` / `dj_check`,
when it has. The `no_d_align` block runs neither — there was not one `set_overlap_safety()` call
between [Genechoice.cpp:538](../src/igor/Core/Genechoice.cpp#L538) and
[:846](../src/igor/Core/Genechoice.cpp#L846) — so with a chosen neighbour the flag was never
written at this event's layer.

A downstream `Deletion` reads `memory_layer_safety - 1`, i.e. exactly that layer. Confirmed by
backtrace on the `[convergence]` inference test:

```
V_choice (Genechoice.cpp:362) → J_choice (:995)
  → D_choice, no_d_align position map (:678)
    → Deletion::iterate V case (Deletion.cpp:267)   ← throws
```

**Why it surfaced only now.** Before B8's `a58808b`, that read was
`Enum_fast_memory_map::at(key, layer)`, which accepted `layer <= current + 1` — one layer *above*
current — set it current, and returned the slot. The storage was `new bool[]`, i.e.
uninitialized, so the overlap verdict for every exhaustive-path scenario was whatever happened to
be in memory. `LayeredArray::get()` refuses that read. The port converted silent undefined
behaviour into an abort.

Note what "current layer" means here: `request_layer()` advances it, but `set()` pulls it back to
the layer written, so after `iterate()` it tracks the last **write**. `V_choice` writing at layer
0 is what leaves `D_choice`'s layer 1 unreadable.

**Fix**: write the conservative verdict `false` ("not established safe") at the top of the
`no_d_align` block, guarded by `v_chosen` / `j_chosen` so it fills exactly the gap the preamble
leaves. Every downstream deletion then performs its own check rather than skipping it — the same
value the alignment loop writes whenever the verdict is undetermined. There is no prior behaviour
to preserve: the old value was uninitialized memory.

**A harness-level guard now covers this class of defect.** `call_iterate_recording()` checks the
layer contract — *requesting a layer is a promise to write it before handing off* — on every
hand-off, across all six seq_type-keyed maps, for every section in the suite. Verified against
this very bug: reverting the fix makes the exhaustive-path sections report
`safety_set key 0 / key 1: requested layer 0, current layer 2`. See
[ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md) §4.1. Its one limit is that it only sees branches
a test executes.

**And the container no longer permits the ambiguity at all** (§7.10): `LayeredArray` now tracks
the claimed layer and the written layer separately, so reading a requested-but-unwritten layer
always throws instead of serving value-initialized storage.

Worth noting for B5/B11: `Gene_choice::initialize_event` requests the safety layers with its
`if (d_chosen)` / `if (j_chosen)` guards **commented out** at
[Genechoice.cpp:1135-1225](../src/igor/Core/Genechoice.cpp#L1135-L1225) — six unconditional
requests against conditional writes. Restoring those guards is the structural fix; the
conservative write above is the stopgap.

**B11 owes the real verdict here.** `d_5_off` and `d_full_3_offset` are known per position, so the
exhaustive path can compute the same three-way outcome the alignment path does. The conservative
write is a stopgap that makes the path defined, not the right long-term answer.

**Process consequences**, both larger than the bug:

1. **`[convergence]` must join the verification ladder** for any step touching source. It is
   excluded from both `pixi run test` (`-LE convergence`) and `test_unit`, which is how
   `a58808b` shipped: that commit recorded "218/218 unit and integration tests, and all four
   regression suites" — and none of those run it. §1's ladder is updated accordingly.
2. **`[!mayfail]` does not contain a crash.** The convergence case carries it, but a SIGABRT
   takes the process down regardless, so an abort there costs the whole suite, not one test.
3. `a58808b`'s claim that "all six `is_overlap_safe()` reads are the read-at-layer-1-then-write
   shape already established" was reasoning by analogy with `Index_map`, where the property was
   actually proved. It was not checked for this map. When porting a container whose accessor
   tightens a precondition, each call site needs the argument made, not inherited.

### 7.10 — `LayeredArray` conflated "claimed" with "written", making detection order-dependent

*(Raised in review Sep 2 2026, implemented the same day.)*

`request_layer()` advanced `layer_of_`, the same counter `get()` validated against and `exists()`
reported. Three consequences, none of them visible until §7.9 forced a look:

1. **A requested-but-unwritten layer was readable**, returning value-initialized storage. §7.9
   threw only because another event's write had pulled `layer_of_` back below the read; under a
   different interleaving the identical missing write returns a default and nothing notices.
   **Whether the container caught a missing write depended on the order of unrelated writes.**
2. **`exists()` meant "requested or written".** Every caller — the `Scenario` view,
   `Single_error_rate`, `Errorscounter`, `DynamicSequenceMap::occupied()` — guards a dereference
   with it, and for a pointer-valued map the value-initialized default is `nullptr`.
3. **B10's three-state distinction was not implementable.** "Not yet processed (layer −1,
   `exists()` false) / actively absent / present" cannot hold when a request makes a key read
   back as written-with-a-default, which for `Int_Str*` is indistinguishable from actively absent.

**The split**: the ownership mark is raised by `request_layer()` and reported
by `claimed_layer()` / `claimed_layers()`; the data mark is moved by `set()` and reported by
`exists()` and `current_layer()`. Invariant `current_layer() <= claimed_layer()`; writing at a
layer claims it.

All 74 ownership call sites across the four events are in `initialize_event` and mean
*"which layer do I own"* (the single hit inside `iterate` is a commented-out `cout`), so they take
the ownership mark unchanged. `current_layers()` likewise: `Rec_Event` snapshots it at init for
`multiply_all`, and it has to name the layers the event owns rather than what happened to be
written when the snapshot was taken.

**This is what makes row 10 of the test guide unconditional.** The harness check no longer needs
downstream events in the fixture to create an observable gap; it compares written against claimed
directly.

### 7.11 — The occupancy-skipping neighbour walk is not the ordering neighbour

*(Found Sep 7 2026 while implementing B6.)*

G9 sketches the generic junction length as

```cpp
const SeqTypeId left = scenario.constructed_sequences.first_occupied_left(seq_type_id);
```

`first_occupied_left()` skips segments that are **written but empty**, which is not what the
hardcoded pairs did. The two agree only while no gene segment can be empty — and one can:
`Deletion` guards its 5' branch with `if (value_int > previous_str.size()) continue`, a strict `>`,
so deleting exactly the whole segment is legal and `substr(size, npos)` writes an empty `Int_Str`.
The V 3' branch has no guard at all.

With an empty D in a VDJ model the two disagree concretely. `DJ_ins_seq`'s ordering neighbour to the
left is `D_gene_seq`, whose offsets are still correct — degenerate, but correct. The occupancy walk
skips D, reaches `VD_ins_seq`, and asks for *its* 3' offset — which no insertion writes (§6.3), so it
throws. Skip that too and it reaches `V_gene_seq`, giving the merged V→J span and double-counting
against the VD junction.

**B6 therefore uses `registry.left_neighbor()` / `right_neighbor()`**, which reproduce the hardcoded
pairs exactly for every topology where all segments are present — including tandem D, where
`D1D2_ins`'s ordering neighbours are `D1` and `D2`. Occupancy skipping is the right answer *once
absence has a defined meaning*, which is B10's decision and explicitly out of scope here (§9).

Two consequences:

- **B7 must make the same choice deliberately.** G9's `Dinucl_markov` anchor lookup has the same
  shape, and its note about a fully-deleted `D1` anchor is the *same* hazard read from the other
  side: the walk is what makes that case safe, and also what changes today's answers.
- **When B10 switches to the occupancy walk, the insertion offsets fixed alongside B6 are a
  prerequisite** — the walk lands on a junction segment as soon as a gene segment is skipped, and a
  junction with no offsets cannot answer.

### 7.12 — `Dinucl_markov` reads a seed nucleotide from an empty anchor, and segfaults

*(Found and reproduced Sep 7 2026 during 2a.)*

`Dinucl_markov::iterate` takes its seed from `previous_seq.back()` (forward) or `.front()`
(reverse), with **no non-emptiness check**. `Int_Str` is a `std::vector<int>`, so on an empty
anchor `data()` is null and `back()` dereferences `nullptr - 1`. Verified: `SIGSEGV`, not a wrong
answer.

§2.9 already flagged this, but concluded it was safe on the current corpus and only reachable under
tandem D. **That is wrong.** The VD junction anchors on V, and an empty V is producible today:
`Deletion`'s V-3' branch computes `previous_str.size() - value_int` with **no size guard at all**,
and the D-5' branch guards with a strict `>` so deleting exactly the whole segment is a legal
realization. Nothing between there and here rejects the resulting scenario — the junction length
stays in range, since the degenerate offsets are consistent.

This is the second reachable crash this work has turned up in a branch no test executed, after
§7.9. Both share a shape: a value that is only *usually* present, read without asking. Together they
are a plausible source of the non-reproducible segfaults reported against the legacy code — the
scenario has to be enumerated in the right order, on the right read, for either to fire.

**Handling.** The reproducer is in the suite as a `[.]`-hidden case (`[empty_anchor]`), because a
segfault aborts the run rather than failing a test, which `[!shouldfail]` cannot express and CI
cannot survive. It asserts what the fix owes: an anchor carrying no nucleotide must be rejected,
never read. **Decided (Sep 8 2026): throw, do not discard.** A scenario whose anchor was fully
deleted is geometrically legitimate, so dropping it silently removes probability mass the model
should account for and leaves no trace that it happened; a scenario that cannot be scored is a
modelling error, and the user has to see it. G9's `first_occupied_*` walk is the fix — but
per §7.11 the walk is *also* a behaviour change, so B7 must land the two together and say which
scenarios move. When it does, the `[.]` comes off.

### 7.13 — `Dinucl_markov` writes through the `Insertion`'s pointer, claiming no layer of its own

*(Raised Sep 8 2026 reviewing 2a.)*

Every other event claims a memory layer for what it writes, so a sibling scenario restores the
previous value on backtracking. `Dinucl_markov` does not: it takes the `Int_Str *` the `Insertion`
left in the constructed-sequence map and fills that buffer in place, claiming nothing and never
writing the map. The filled junction therefore stands at the *Insertion's* layer.

**This is sound only because the two events behave as one.** Neither branches — an insertion's
realization is determined by its neighbours and a Dinucl_markov's by the read — so there is never a
sibling scenario to corrupt, and the `Insertion` re-assigns the buffer on its next call regardless.
That is a property of the current pair, not of the layer contract.

It breaks as soon as either side gains a branch, and both are plausible:

- an `Insertion` that enumerates lengths rather than deriving one — which is what an indel-aware
  error model needs, and which §8.1 already lists as a coming requirement;
- a `Dinucl_markov` that branched over the nucleotides an ambiguous read position stands for.

In either case two sibling scenarios share one buffer, and the second reads the first's nucleotides
where it expects placeholders. Silently: the fill is guarded by `ins_seq.at(i) == -1`, so an
already-written position reads as *someone has filled this*, which is exactly what the second
sibling must not conclude.

**Pinned as a `[!shouldfail]` case** stating the intended behaviour — the filled junction at this
event's own layer, the placeholders still readable at the layer below. Not fixed in B7: it changes
which layer a downstream reader finds the junction at, so it needs its own commit (R1, §6.9) and its own
regression run. It is also a prerequisite for either branching change, not a tidy-up.

### 7.14 — "Not filled yet" was a bare `-1`, on the same axis as `int_N`

*(Raised and addressed Sep 8 2026 reviewing 2a.)*

`Int_Str` carried three states on one axis, only two of them named:

| State | Was | Meaning |
|---|---|---|
| position absent | not in the string | the segment is shorter, or not there |
| present, undetermined | `-1` | allocated by `Insertion`, not yet filled |
| present, ambiguous | `int_N = 14` | filled; the read does not say which base |

The middle one existed **only as a private protocol between `Insertion` and `Dinucl_markov`** —
written as a bare literal in five places, produced by one event and consumed by one other, with no
name, no type support and no consumer contract. It is the same three-state problem as B10's
absent-segment semantics, one level down.

**Addressed by naming it, at 15 rather than −1** — `int_undefined`, immediately past the real codes,
with `kIntNtCount` defined *from* it so the two cannot drift. The position is what makes it more
than a rename:

- `dinuc_proba_matrix` is `kIntNtCount` square, so `matrix(int_undefined, j)` trips its bounds
  assertion, where `matrix(-1, j)` read one row *before* the array — silently, since
  `Matrix::operator()` asserts only the upper bound and asserts are compiled out under `NDEBUG`;
- the guard `if ((first_nt_index < 4) & (sec_nt_index < 4))` now means what it reads as. With `-1`
  an undefined nucleotide **satisfied** `x < 4` and went down the *unambiguous* path, indexing the
  marginal array at `base_index + (-1) * 4 + sec` — four entries before the intended block, no
  crash, no diagnostic. With 15 it takes the ambiguous path and lands on a bounds-checkable index.

Unreachable today, because every anchor is a gene segment. It becomes reachable the moment §7.11's
occupancy walk lands: skip an empty D and `DJ_ins_seq`'s anchor is `VD_ins_seq`, whose last position
is undefined until its own `Dinucl_markov` has run. §7.11, §7.12 and this entry are one cluster.

**Enforced at the one boundary where it has to hold.** `Rec_Event::iterate_wrap_up`'s leaf branch
asserts that no segment reaching the error rate still carries a placeholder, via
`first_unfilled_segment()`. Under `#ifndef NDEBUG` only: the walk is linear in the scenario's
length, and the default build is RelWithDebInfo, so release pays nothing — verified by the message
string being absent from `libigorCore.so`. The predicate is exposed rather than buried in the assert
so it can be unit-tested without a debug build, and the check itself was verified by forcing the
guard on with the fill disabled, which reports the offending seq_type and aborts.

Tests render an undefined position as `.` and an ambiguity code as `N`, so a failure message cannot
blur the two. That change alone corrected four assertions that described a freshly-allocated
junction as `"NNN"` — it is `"..."`.

**Still open**: whether to *drop* the state rather than name it. `Insertion` knows both neighbour
offsets, so it could write the read window itself and leave `Dinucl_markov` computing only
probability — no undetermined state anywhere. That forecloses a branching `Dinucl_markov`, and it
depends on §7.13, since sharing a buffer is what makes writing earlier equivalent. Decide it with
B10: "allocated but undetermined" and "processed but absent" are the same question asked of a
nucleotide and of a segment, and deciding them apart risks deciding them in opposite directions.

## 8. Decisions taken

O1–O6 from the Sep 1 2026 review; O7–O9 from the Sep 9 2026 re-assessment (§6.8).

| # | Question | Decision |
|---|---|---|
| O1 | Safety-flag key and layering | **Storage stays *n(n−1)/2*.** Triangular matrix over the ordering, filled by nearest-neighbour check + row-suffix propagation; row-bitmask container so propagation is O(1) and the per-event layer count stays at two. §2.3 rewritten. |
| O2 | Is A0 an accepted amendment to D2? | **Yes.** Record the amendment in the parent plan when A0 lands. |
| O3 | Where does `PendingModifierBounds` live? | **Per event instance.** Not only to match current behaviour: its content depends on `processed_events`, i.e. on where the event sits in the recursion, so it is not shareable via `ModelContext`. |
| O4 | Does §7.1 get fixed here? | **Reproduce first, fix last.** The sign slip is preserved verbatim through steps 1–5 for regression-testing purposes, then fixed as the final commit with unit tests pinning the corrected core length. |
| O5 | Milestone-2 absence semantics (a)/(b)/(c) | **Defer; update the parent plan once step 5 is carried.** G5's pair-keyed junctions make (b) cheaper than the parent plan's estimate — re-score it then, not now. |
| O6 | Second `no_d_align` fixture | **The switch is required; only the fixture rationale changes** *(Sep 1 2026; a Sep 9 amendment claiming it was dissolved was wrong and is withdrawn — §6.8 F5)*. The per-`Gene_choice` boolean is not optional: per §2.6 the generic body has no `case D_gene` left to confine the fallback to, so the flag is what carries the V/J-vs-D asymmetry — default `true` for `D_gene`, `false` for V and J, bitwise by construction. What the amendment got right is narrower: the *fixture* need not defeat the aligner, since an empty alignment list (or a pruning threshold above every bound) reaches the path through production code, as T0's G6 sections already do. Whether V and J *should* fall back stays open policy. |
| O7 | Where do the decided behaviour fixes land? | **All of them after 5b, as phase R** (§6.9). Landing a fix mid-sequence would move the golden data partway through, after which "bitwise" no longer means one thing across the remaining steps and every verdict has to be read against which baseline it was taken on. The refactoring block stays idempotent end to end; F is the one place golden data may move, once, with each commit naming the outputs it changes. Cost accepted: re-establishing context on `Insertion` and `Dinucl_markov` later. |
| O8 | Is S4 just the `has_effect_on` overrides? | **No — S4 is the pair-keyed junction structure** (§6.8 F3). The six `*_length_best_proba_map` members plus `vj_length_d_position_proba` collapse to one map keyed by an ordered `(SeqTypeId, SeqTypeId)` behind a single accessor. This is not scope creep: `iterate_initialize_Len_proba` is enum-keyed, so a tandem-D junction throws before inference starts, which puts S4 **on the milestone-1 critical path**. It is also C2's stated landing point. Per §9 the query returns a set, not the first match. **Split into S4a/S4b/S4c by the §6.10 analysis, approved Sep 10 2026**, with an optional S4d for Tensor-backed containers and the joint-max bound tightening moved out to R6 as a behaviour change. The composition operator divides across the split: `⊗ᵐᵃˣ` in S4c because every gene needs it, `⊗ᵉⁿᵘᵐ` in 5b because only `no_d_align` retains the decomposition. |
| O9 | 4a before or after S5? | **Before.** S5 replaces the safety mechanism `Deletion::iterate` reads; characterizing against a body S5 has already moved is the wrong order. It also gives S5 a consumer rather than making it a third service with none (§6.8 F1), and S5's definition of done becomes "4a's sections pass unchanged". |

### 8.1 — Two standing design constraints

Neither lands in this work; both change what "done" looks like for it.

**C1 — the rewrite must make `iterate()` unit-testable.** It currently has no unit tests, and the
`feature/2_unittests` sketch is a harness, not a suite (§6.1). This is not a by-product to hope
for: T0 lands the harness *and its missing assertions* before any production change, and every
subsequent step extends it with the branches it collapses. Concretely, the generic form helps three ways the switch
form cannot — a single body means one test per *pattern* rather than one per V/D/J branch;
`PendingModifierBounds` and `reachable()`/`check()` are free functions over plain values, testable
with no model loaded at all (the same property that made `LayeredArray` testable in B2); and the
G6/G7 switches turn three behaviours that are currently reachable only by choosing the right gene
class into direct inputs.

**C2 — in/del error models are coming, and the hook points are already marked.** IGoR does not
support insertions/deletions as *errors* today. Twenty comments say so, and they are all on the
same construct — the junction-length feasibility guard:

> `continue; //This means no scenario can lead to a correct solution, would need to be changed for
> Error models with in/dels`

13 in [Genechoice.cpp](../src/igor/Core/Genechoice.cpp), 6 in
[Deletion.cpp](../src/igor/Core/Deletion.cpp), 1 in
[Dinuclmarkov.cpp](../src/igor/Core/Dinuclmarkov.cpp). Support must **not** land here, but the two
places it will touch are exactly the two this plan consolidates, so leaving them in a shape that
accepts it costs nothing now and a rewrite later:

| Affected | Today | With in/dels | Consequence for this plan |
|---|---|---|---|
| junction-length lookup (G5) | exact-length lookup; miss ⇒ `continue` | the read-space span and the constructed length differ by the net indel count, so a span `L` is compatible with lengths in `[L − max_ins, L + max_del]`, indel-weighted ⇒ a **range query** | keep the lookup behind **one accessor** on the pair-keyed structure, so the exact-match becomes a range scan in one place instead of 20 |
| reachable-offset interval (G2) | width comes from pending deletions only | widens by the indel budget between the two ends | keep `reachable()` a **function of the bounds object**, not of `Deletion` internals, so an indel term is an added contribution rather than a new call site |

Worth knowing: the **alignment layer already carries the data**. `Alignment_data` exposes
`get_all_insertions()` / `get_all_deletions()` alongside `get_all_mismatches()`, and the aligner
populates them ([Aligner.cpp:989-995](../src/igor/Core/Aligner.cpp#L989-L995),
[:1091-1095](../src/igor/Core/Aligner.cpp#L1091-L1095)). It is the scenario layer that cannot
consume them — which is why the 20 comments all sit on the same guard.

Neither shaping decision costs anything today; both are already implied by S2/S3/S4 as specified.
Recording them so the eventual indel work has a stated landing point rather than 20 comments.

## 9. What this plan does not cover

- **B3 (flank seq types)** — independent; the `reachable()`/`check()` services are indifferent to
  flanks because flanks carry no offsets.
- **B10 (absent-segment semantics)** — milestone 2. G9's skip-empty walk is the mechanism B10
  needs, but the *modelling* decision (parent plan, "Milestone 2 design decision") is untouched
  here.
- **The `iterate()` docstrings** — the parent plan's *Documentation debt* section already makes
  rewriting each one part of the definition of done for B5/B6/B7/B11. Still applies; each step
  above inherits it.
- **Phase D decomposition** — G5's junction pair key is plausibly the interface variable Phase D
  needs (`D.3 — Interface Variables`), but that is not established here.
- **Re-keying `Events_map`** — see below. S2 exposed the limitation; fixing it is a `Model_Parms`
  change with its own blast radius, and nothing in B5–B11 needs it.

### The `Events_map` key cannot express more than one modifier per end

`Model_Parms::get_events_map()` keys by `(Event_type, seq_type, Seq_side)` and inserts with
`emplace`, so a second event of the same type bearing on the same end is **silently dropped**, not
added. No topology IGoR builds today produces one, but nothing in the model formalism forbids it:
a junction with two offset modifiers is a perfectly reasonable thing for a future model to
declare, and tandem D is the direction that makes it likely.

S2 is unaffected — `PendingModifierBounds` accumulates by addition and never queries the map by
key, so it composes correctly the moment the map can hold both. The test that pins that composition
has to insert under a synthetic key precisely because the real keying cannot produce the input
(`test_junction_geometry.cpp`, "two modifiers on one end sum").

**The eventual fix is to query by capability rather than by key** — the map becomes a flat list of
events and consumers ask `get_offset_delta_bounds` / `get_offset_role` / `has_effect_on` which ones
bear on the end in question, which is what A0 exists to make possible. Every consumer S4 touches
should be written so that it does not care how many events answer. Recorded as future work, not
scheduled: it is a `Model_Parms` change, and every remaining `try_get_event` caller is a site that
would have to stop assuming a unique answer.
