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

### 2.2 — G2: Reachable-offset interval, and the overlap predicate

**Where it is today**: eight blocks — `Gene_choice` V/D/J (three), `Deletion` V/D-5′/D-3′/J (four
sites, five checks), each computing a `*_min_offset` / `*_max_offset` pair and then applying a
three-way comparison. [Deletion.h:120-146](../src/igor/Core/Deletion.h#L120-L146) holds eight
scalars for this; `Gene_choice` holds the same eight again.

**The interval**:

```cpp
struct OffsetInterval { Seq_Offset lo, hi; };

OffsetInterval reachable(SeqTypeId id, Seq_side side, Seq_Offset current) const {
    const OffsetDelta d = bounds.offset_delta(id, side);
    return { current + std::min(d.min, d.max), current + std::max(d.min, d.max) };
}
```

The `min`/`max` is what absorbs the 5′-vs-3′ sign flip. Verified against all eight sites:
`d_5_min_offset = d_5_offset - d_5_min_del` with `d_5_min_del = get_len_max() = -min_del`
gives `d_5_offset + min_del` = `current + delta.min` ✓, and symmetrically for the other seven.

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

`G` is `0` in every check performed today (all in-between segments have `len_min == 0`, and
`Insertion` events legitimately allow zero insertions). Keeping `G` explicit costs nothing and
is what lets the predicate stay correct once a tandem-D ordering puts a *gene* segment between
two checked ends. **Bitwise-preservation note**: `G` must be computed from `len_min`, which is
`0` for the current corpus, so this is a no-op today — do not "improve" it to a tighter bound in
the same step (§7.2).

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
> non-emptiness check**. For VD (anchor V, full deletion forbidden) and DJ (anchor J, `.front()`
> on a J that has only been 5′-deleted) this is safe today. For a tandem-D `D1D2_ins` whose
> anchor `D1` has been fully deleted — legal, `>=` guard — it is undefined behaviour. The
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
| **T0** | Port the harness from `feature/2_unittests`; drop the 11 non-`iterate()` test cases; write ~22 pattern-keyed sections against observed behaviour (§6.1) | unit | n/a — tests only |
| **A0** | ✅ **done** — `OffsetDelta` / `LengthContribution` + four capability virtuals on all four subclasses | unit | yes — no caller yet |
| **S2** | `PendingModifierBounds` in `JunctionGeometry.h`; unit-tested against a VDJ and a VJ model | unit + mutation | yes — no caller yet |
| **S3** | `reachable()` + `Overlap check()`; unit-tested on all eight current comparison shapes | unit + mutation | yes — no caller yet |
| **1** | **B6** — `Insertion::iterate` generic (G9). Smallest, one hot-loop win. | full ladder | **yes** |
| **2** | **B7** — `Dinucl_markov` specs from the registry (G9); skip-empty walk; per-spec buffers | full ladder + the empty-anchor test | **yes** |
| **S4** | Junction pair key; `has_effect_on(left,right)` in the base; the three overrides deleted | full ladder | **yes** |
| **3** | **B11a** — `Gene_choice` alignment path generic (G4, G2, G8, and G5 via S4) | full ladder + benchmark | **yes**, except §7.1 |
| **S5** | Safety row-bitmask; row-suffix propagation; `Event_safety` deleted | full ladder + the empty-segment transitivity test | **yes** (§2.3 corollary) |
| **4** | **B5** — `Deletion::iterate` generic (all patterns) | full ladder + benchmark + convergence | **yes** |
| **5** | **B11b** — `no_d_align` exhaustive path generic (G6) | full ladder + a fixture that *forces* the path | **yes** |

B11 is split: the alignment path (3) is what milestone 1 needs and is straightforward; the
exhaustive path (5) is the hardest single piece in the whole set and depends on S4 being settled.
The parent plan's guidance to build the milestone-1 fixture so `no_d_align` never fires stands —
step 5 is exercised instead by flipping the G6 fallback switch on a gene with no alignments, which
needs no aligner manipulation.

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
| **Endogenous mismatch bound** (G8, §7.1) | the surviving-core window and its credited length | `V — pins the match-length sign slip (§7.1)`; `J — pins the same slip`; `D — both ends truncated` | `Alignment with mismatches` (1 assert), `Both ends truncated (max deletions)` (1 assert) |
| **Pruning** | `should_prune` at the realization boundary | `below threshold ⇒ realization skipped` | `Probability threshold filtering` |
| **Baseline writes** (G4) | the offset/sequence/mismatch triple | `V`; `D`; `J`; `negative alignment offset (pre-alignment strip)`; `zero-length junction between adjacent segments` *(new)* | `Basic V/D/J alignment - no other genes` (5/3/3 asserts), `V gene with negative offset` |

**≈22 sections in 6 `TEST_CASE`s**, against 33 sections in 12 `TEST_CASE`s. Two are written fresh
(the G6 switch default, and the zero-length junction — the latter because B10's degenerate offset
convention is what makes §2.3's propagation proof hold at the equality case, and nothing currently
pins it).

#### Assertions come from observed behaviour, not from the sketch

Characterization in the strict sense: run the **unmodified** implementation, capture what it does,
assert that. Not from the sketch's comments, and not from what the code *ought* to do — several of
these branches are exactly where §7 says it does the wrong thing, and pinning intended behaviour
would make T0 fail before step 1 starts.

Where a pinned value is known-wrong, name it so the eventual fix has an obvious landing site:

```cpp
SECTION("V — pins the match-length sign slip (plan §7.1)") { … }
```

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

T0 covers `Gene_choice` only, matching the sketch's reach. Sections for `Deletion`, `Insertion` and
`Dinucl_markov` are added as part of steps 4, 1 and 2 — each step contributes the instances for the
branches it collapses, into the pattern-keyed `TEST_CASE`s above. The suite grows with the
migration rather than becoming a project of its own, and the row for each pattern fills in as the
events are migrated.

### 6.2 — The regression gate has a flaky output

*(Observed Sep 2 2026 during A0.)*

`pixi run test_regression` is designated the bitwise gate for every step of this plan. One of
its outputs is **not deterministic**: `best_scenarios_counts.csv` mismatched once on an A0 build
that adds only uncalled virtuals, then matched on the next three runs (two `test_inference`, one
full `test_regression`), while the pre-A0 baseline matched on its single run.

Mechanism: `Best_scenarios_counter::dump_sequence_data()` is called from
[GenModel.cpp:512-514](../src/igor/Core/GenModel.cpp#L512-L514) inside
`#pragma omp critical(dump_counters)`, within the parallel loop over query sequences. `critical`
serialises but does **not** order, so the rows land in whatever order threads finish sequences.
Row *content* is deterministic — the tie-break at
[Bestscenarioscounter.cpp:133](../src/igor/Core/Bestscenarioscounter.cpp#L133) is a strict `>`,
so ties keep enumeration order, and enumeration within one sequence is single-threaded — but row
*order* is not.

**This matters more than its size suggests.** A gate that fails intermittently trains its readers
to re-run rather than investigate, which is exactly how a real regression gets waved through
during a multi-step refactor. Every "bitwise" claim in §6 depends on this gate meaning what it
says.

Fix shape (not done here — out of A0's scope): either sort rows by `seq_index` before comparing,
or make the dump ordered. Confirm first whether the comparator is a plain `diff`; if it is,
sorting in the comparator is the smaller change and does not touch inference code. Worth doing
**before step 1**, since steps 1–5 all lean on this gate.

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

## 8. Decisions taken (Sep 1 2026 review)

| # | Question | Decision |
|---|---|---|
| O1 | Safety-flag key and layering | **Storage stays *n(n−1)/2*.** Triangular matrix over the ordering, filled by nearest-neighbour check + row-suffix propagation; row-bitmask container so propagation is O(1) and the per-event layer count stays at two. §2.3 rewritten. |
| O2 | Is A0 an accepted amendment to D2? | **Yes.** Record the amendment in the parent plan when A0 lands. |
| O3 | Where does `PendingModifierBounds` live? | **Per event instance.** Not only to match current behaviour: its content depends on `processed_events`, i.e. on where the event sits in the recursion, so it is not shareable via `ModelContext`. |
| O4 | Does §7.1 get fixed here? | **Reproduce first, fix last.** The sign slip is preserved verbatim through steps 1–5 for regression-testing purposes, then fixed as the final commit with unit tests pinning the corrected core length. |
| O5 | Milestone-2 absence semantics (a)/(b)/(c) | **Defer; update the parent plan once step 5 is carried.** G5's pair-keyed junctions make (b) cheaper than the parent plan's estimate — re-score it then, not now. |
| O6 | Second `no_d_align` fixture | **Dissolved into the G6 switch.** The fallback becomes an explicit per-`Gene_choice` boolean, defaulted to reproduce legacy (D on, V/J off); step 5 is tested by flipping it, not by defeating the aligner. |

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
