# The junction-length probability bound: what it was, what it is becoming

Orientation for reviewing the S4 series. The migration steps live in
[ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) §6; this document explains the
machinery they change and states, per step, **what each one actually buys**.

Four labels are used throughout, and every delivered item is tagged with all that apply:

| tag | meaning |
|---|---|
| **REFACTOR** | duplication removed, behaviour bitwise-identical |
| **GENERALIZE** | drops a legacy VDJ `Seq_type` assumption — admits richer topologies, or breaks a direct event-to-event linkage |
| **PERF** | faster initialization, or faster bound access in `iterate()`'s hot loop |
| **CORRECT** | fixes a defect that produced wrong results, crashed, or refused a legal model |

---

## 1. What the machinery is for

IGoR explores a tree of recombination scenarios per read. The tree is astronomically large — 10⁸ to
10¹⁰ nodes for human IGH — so it is pruned: at each node, an **upper bound** on everything still
downstream is multiplied into the probability accumulated so far, and if the product falls below a
threshold the whole subtree is discarded.

`ExplorationContext::compute_upper_bound()` forms that product, and one of its factors is the
subject here: **given the geometry the scenario has already committed to, what is the best
probability the model can assign to the events that must fill the gap?**

A worked instance. A scenario has chosen a V gene and a D gene, so it knows where V's 3′ end and
D's 5′ end sit on the read. Between them is a fixed number of positions, and that number has to be
produced by *some* combination of a V deletion, a D deletion, and a VD insertion. The bound asks:
over every such combination, what is the highest probability the model gives? If even the best one
is too improbable, no scenario through this node can survive, and the subtree goes.

So the machinery has two halves, with very different cost profiles:

| half | when | cost |
|---|---|---|
| **build** | once per thread per EM iteration, before the per-sequence loop | **model-only**, so it amortises over a batch — but is paid in full by a one-sequence query. Measured below |
| **consume** | at every scenario node | the hottest loop in IGoR |

Both halves are addressed, but they are not equally worth optimising, and the plan is explicit about
which improvements land where.

### How expensive is the build, really?

*(An earlier draft of this document called the build "effectively free". That is true of a large
batch and false of a small one — a CLI user querying a sequence at a time pays the whole thing per
invocation, and it is model-only so no batch size makes it cheaper per model. Measured Sep 11 2026;
`tst/igor/Core/test_proba_bound_benchmark.cpp` keeps it honest.)*

Single-threaded, the sweep alone, uniform marginals — `igor_tests "[.benchmark][proba_bound]"`:

| model | topology | before S4c | after S4c |
|---|---|---:|---:|
| human TCR-α | VJ | 4.26 ms | 4.06 ms |
| TRB regression corpus | VDJ | 203 ms | **143 ms** |
| human TCR-β | VDJ | 205 ms | **146 ms** |
| human BCR-heavy | VDJ | 10.46 s | **9.74 s** |

> **The first version of this table was wrong, by a factor of three.** `d56d1b1`'s fixture ran the
> bound sweep without the `initialize_event()` pass that precedes it in `GenModel::infer_model`, so
> `Gene_choice`'s `v_chosen` / `d_chosen` / `j_chosen` flags were read uninitialised and **no gene
> choice folded anything at all** — the numbers were `Deletion` and `Insertion` only. It went
> unnoticed because before S4c the sweep picked its junctions from an enum switch and so ran
> *something* regardless; S4c made the sweep depend on init, which turned a wrong measurement into
> an obviously empty one. The fixture now runs the init pass, and both columns above are measured
> with it.

Two later changes moved it again, and neither is comparable to the table above: the machine reads
~25 % slower on Sep 16 than it did on Sep 11, so both are measured **interleaved against each
other** in one session (20 samples, three repetitions, rebuilding between each) rather than against
the Sep 11 column.

| model | dense profile (`18cf2e4`) | **+ O14**, the flat participant array | |
|---|---:|---:|---:|
| human TCR-α | 4.26 ms | **1.40 ms** | 3.0× |
| TRB regression corpus | 152 ms | **46.1 ms** | 3.3× |
| human TCR-β | 160 ms | **49.8 ms** | 3.2× |
| human BCR-heavy | 8.70 s | **2.59 s** | 3.4× |

Those are the *coldest* pair of each — the first repetition after a rebuild. Across the three
repetitions the old code drifts upward (TRB 152 → 187 → 209 ms) and the new one does not (46.1 →
44.5 → 44.8 ms), so **3× is the conservative reading of this change and 4× the warm one**. The
asymmetry is itself informative: the old fold's cost is dominated by allocation and refcount
traffic, which is what degrades as the machine warms.

**In a real inference run it is worse than the table suggests.** Timed in situ on the TRB corpus
with 22 threads: **706 ms mean, 1214 ms max — per thread**, against the 203 ms the same work takes
alone. The ~3.5× is contention, and every thread is computing **the identical result**: the bound is
model-only and thread-invariant, so N threads redundantly build N copies of one answer. That is a
second redundancy stacked on the per-consumer one, and a larger one. *(An earlier draft put this
ratio at 12×, by comparing the in-situ figure against the broken benchmark's 58.5 ms. S4e's case
rests on the N-copies redundancy, which is untouched by the correction; only the contention
multiplier moves.)*

*(That in-situ figure predates S4c. The sweep is ~4× cheaper single-threaded since, so the absolute
per-thread cost is now a fraction of 706 ms — but S4e's case never rested on the absolute: N threads
still compute the identical model-only answer N times, and that is unchanged.)*

Where the VDJ time goes, per event (in situ, TRB):

| event | mean | live? |
|---|---:|---|
| `GeneChoice_J_gene` | 501 ms | **live** — J is priority 7 and the D gene 6, so `d_chosen` is false at J and it legitimately anchors on V, needing the full V→J span |
| `Deletion_V_gene_Three_prime` | 176 ms | **dead** — `Deletion(V,3')` is priority 5, *below* the D gene, so `d_chosen` is **true** and its `iterate` reads the VD map. The V→J map it spends this building is never read |
| everything else | ~30 ms | — |

So **§6.10 finding 3's dead map is ≈25 % of VDJ initialization**, and deleting it in S4c is a
measurable saving rather than tidiness. **Borne out** *(S4c, Sep 11 2026)*: the benchmark's two TRB
models each lost 29 % of the sweep. BCR-heavy lost only 7 %, because there the cost is dominated by
`Gene_choice(D)`'s retained decomposition — `|D| × |VD| × |DJ|` — which S4c does not touch. The leaf counts corroborate the split: `GC(J)`'s V→J fold is
≈1.39 M leaves against `Deletion(V)`'s ≈464 k, and 501/176 ≈ 2.8 matches the 3×.

A tandem-D row is still missing from the benchmark, but no longer because the bound machinery
refuses it: after S4c an event holds junction *handles* and nothing enumerates VD / DJ / VJ. What a
D1D2 row now waits on is a model that declares the segments — the registry and event wiring, not
this structure.

---

## 2. Three things that "the length of a span" runs together

The single most useful correction to come out of this work, because the legacy naming hid it.

| | what it is | scalar? |
|---|---|---|
| the **span** | a region — a pair of boundaries | never |
| its **profile** at a consumer | `{length → best probability}` over contributors the scenario has **not yet** committed to | only when nothing is left unresolved |
| the **observed distance** | one integer, read off the boundary slots as they currently stand | always |

The bound is `profile[observed]` — the third used as a key into the second. The profile is what the
model still allows; the observed distance is what the **scenario** has already committed to — the scenario being the confrontation of model and read, not
the read itself, which is why the observed distance differs between consumers while the read plainly
does not.

Two consequences that the legacy code embodies without stating:

- **A boundary is a slot, not a constant.** `Gene_choice` creates an offset; `Deletion` moves it. So
  the same span has a different profile depending on where in the priority ordering it is read.
- **The profiles are a suffix family.** `GenModel` pops the initialization stack in *reverse*
  priority order and hands each event the queue *after* itself, so every event's profile folds
  itself plus its suffix. Upstream contributors are excluded because the scenario has already
  committed to them — which is exactly what keeps the observed distance and the profile in the same
  frame, with no bookkeeping.

That second point means the six per-event maps are **not** six copies of three things. They are
different objects, keyed by different quantities. See §6.10 finding 4.

---

## 3. The legacy machinery

### Build

```
GenModel: for each event, in reverse priority order
    event->initialize_Len_proba_bound(queue_after_this_event, marginals, index_map)
        -> Rec_Event::iterate_initialize_Len_proba(Seq_type junction, map&, ...)
            -> virtual iterate_initialize_Len_proba   (4 overrides)
                 for each realization: multiply maxᵢ P(r|i), add Δ(r), recurse
            -> iterate_initialize_Len_proba_wrap_up   (pops the next event, recurses)
            -> at the leaf: map[accumulated_length] = max(existing, accumulated_proba)
```

Storage, across two classes:

| holder | members |
|---|---|
| `Gene_choice` | `vd_`, `dj_`, `vj_length_best_proba_map`, plus `vj_length_d_position_proba` |
| `Deletion` | `vd_`, `dj_`, `vj_length_best_proba_map` |
| `Insertion` | `junction_length_best_proba_map` |

### Consume

At every scenario node, roughly twenty sites of the shape

```cpp
if (vd_length_best_proba_map.count(d_offset - v_3_off - 1) <= 0) { /* discard */ }
downstream_proba_map.set(VD_ins_seq, vd_length_best_proba_map.at(d_offset - v_3_off - 1), layer);
```

### The supporting cast

- `has_effect_on(Seq_type)` — a predicate each of the four overrides called on itself.
- A `Seq_type_str_p_map` threaded through the whole fold as a **side channel**: `Insertion` stashed a
  dummy `Int_Str` of the right length in it so that `Dinucl_markov`, running later in the same
  traversal, could read `->size()` back out and raise its per-nucleotide probability to that power.
  The `//TODO constructed sequences should not be used but it is useful to compute the dinucl
  contribution` on `Rec_Event.cpp` is exactly this.

---

## 4. What was wrong with it

| # | problem | kind |
|---|---|---|
| 1 | **`iterate_initialize_Len_proba` is `Seq_type`-keyed.** A tandem-D junction has no enum value, so a tandem-D model throws before inference starts | **blocker** |
| 2 | `has_effect_on` answered *one* question — does this change the span's **length** — under a name suggesting three, and its argument was a **span**, not a seq_type: `VJ_ins_seq` means a segment in a VJ model and the whole V→J span in a VDJ one | naming / generality |
| 3 | The queue-level filter was **commented out** behind `//TODO fix this and find a way not to loop over all events`, so the traversal visited every event in the model and each self-filtered — the predicate existed twice, at the wrong level | duplication / init cost |
| 4 | Four near-identical overrides differing only in one scalar | duplication |
| 5 | A whole sequence map carried through the fold to move **one integer** between two events | coupling |
| 6 | `Insertion::initialize_Len_proba_bound` ran the entire traversal `\|R\|` times where once would do | init cost |
| 7 | Every consumption site pays **two** red-black descents (`count` then `at`) on the same key, per scenario node | hot-loop cost |
| 8 | `Insertion.cpp:209` calls `.at()` with no `count` guard, unlike the other ~18 sites, so it throws where they discard | latent defect |
| 9 | `Deletion` builds `vj_length_best_proba_map` in every VDJ model and nothing reads it; `Gene_choice::has_effect_on` is `true` only in a case no consumer reaches | dead code |
| 10 | Six named members hold what the enum naming presents as three logical maps, across two classes | duplication / generality |

---

## 5. What the new machinery is

### Addressing

```cpp
struct SegmentBoundary { SeqTypeId id; Seq_side side; };   // a half-open range bound
struct SegmentSpan     { SegmentBoundary left, right; };   // with a gap(l, r) factory
```

`cut_position()` states the coordinate conversion once: `Seq_Offset` is a *closed*-interval
nucleotide index, and closed intervals do not compose by subtraction. As half-open bounds
(`5' → offset`, `3' → offset + 1`) a span is a plain `end - begin` and adjacent spans compose
because `[a,b)` and `[b,c)` make `[a,c)`.

### Queries

| | |
|---|---|
| `affects_length_of(SegmentSpan)` | does a realization change the span's accumulated length? |
| `affects_proba_of(SegmentSpan)` | does this event contribute a probability factor to it? |
| `participates_in_span(...)` | the disjunction the traversal filters on |

The split was forced, not cosmetic: `Dinucl_markov` contributes **zero** length and answered `true`
purely to gate its `p^L` factor. A length-only predicate at the queue would have silently dropped
that factor from every bound.

### The fold

One non-virtual body on `Rec_Event`, behind two hooks:

| | |
|---|---|
| `length_delta(const Event_realization&)` | the scalar the four bodies differed by |
| `span_proba_factor(SegmentSpan, const SpanAccumulator&)` | the factor, defaulted to 1 |

`SpanAccumulator` carries per-segment lengths, replacing the side channel. Only a segment's
**creator** publishes (`SeqConstructionRole::Creates`), so each key has exactly one writer per path.

---

## 6. What each step buys

### Delivered

| step | REFACTOR | GENERALIZE | PERF | CORRECT |
|---|---|---|---|---|
| **S4a** predicate + queue filter (`b432247`) | four per-body self-filters deleted — problem 3, 4 | `SegmentSpan` names the addressing unit; the `Seq_type`-as-span overload stops being implicit — problem 2 | traversal depth now proportional to *contributing* events, not model size. **Init only** | — (bitwise) |
| **S4a** tier-3 harness check (`43b3efd`) | — | — | — | test-only, but **found a real defect**: see below |
| **S4a** layer-ownership rule (`f568bd4`) | — | — | — | test-only; same finding |
| **O11** boundary-addressed spans (`9af1367`) | — | boundaries can name a **segment's own extent**, not only the gap between two segments — which a per-`Seq_type` Phase D decomposition needs and the segment-pair form could not express at all | — | preventive: `cut_position()` states the ±1 convention once, in the area §7.1 and §7.8 are both off-by-one bugs in |
| **S4b** the collapse (`ce3e4b0`) | four bodies → one; the side channel and its `Seq_type_str_p_map` parameter gone — problems 4, 5 | the fold takes a `SegmentSpan` and keys the accumulator by `SeqTypeId`; `Dinucl_markov` no longer reads a length out of a map `Insertion` wrote, so that **event-to-event linkage is broken** | problem 6 removed — but **worth ~0 ms**, see below | — (bitwise) |
| **S4c** the re-keying | six enum-named members plus three `memory_layer_proba_map_junction*` scalars → `std::array<JunctionBound,3>` on `Rec_Event`; four `initialize_Len_proba_bound` overrides → one non-virtual driver plus a `finalize` hook only `Gene_choice(D)` uses — problem 10 | **the tandem-D enum ceiling is gone — problem 1**: nothing in the bound machinery enumerates VD / DJ / VJ, and an event holds at most a left, a right and an enclosing junction whatever the topology | value-or-absent accessor replaces `count`+`at` at all 19 consumption sites (problem 7); `record()` is one descent where the fold took up to three; problem 9's dead fold deleted — **29 % of the sweep on both TRB models**, 7 % on BCR-heavy | deletes problem 9's dead code; problem 8's guarded / unguarded asymmetry is gone as a class — (bitwise) |
| **O14** the flat participant array | the fold's three signatures stop naming a `std::queue` at all; `SpanParticipants` is a plain `vector<const Rec_Event *>` | the participant list is **immutable and shareable**, which the by-value queue was not — the precondition S4e needs to hand one fold's work to every thread | the queue copy (a deque allocation and a `shared_ptr` refcount pair, per event, per node) and the per-node participation filter both collapse to once per junction: **3–4× on the whole sweep**, every model | — (bitwise, plus a 17-digit profile-equality probe) |
| **S4c** the dense profile | — | — | `SpanProfile` becomes a `std::vector<double>` indexed by distance, replacing the `std::map` S4c kept: **`best_for` falls from 8 % to 2 % of `iterate()`**, inference 20.2 → 18.3 s on the N=1000 pipeline, and the init sweep 9.74 → 8.33 s on BCR-heavy. The key space is the sumset of the contributors' realization ranges, so it is contiguous and at most 2 kB — problem 7, finished | — (bitwise) |

**Nothing delivered so far is a correctness fix**, and that is by design — every step above is
bitwise on the regression corpus. What the S4a test infrastructure *did* do is **surface** a defect:
`Insertion` writes the segment it creates without ever requesting a layer for it, the only violation
of "a written layer must have been requested" among 49 measured writes. It is repaired in **R3**, and
**R3b** then hardens `LayeredArray::set()` so the rule holds at runtime rather than only under test.

**S4c is the first step that touches the hot loop**, and the two PERF entries before it are worth
less than they look. Both are initialization, and the per-event breakdown above shows where
initialization actually goes: `Insertion`'s own sweep costs **0.0 ms** on the TRB corpus, so
finding 6's `\|R\|`-fold removal — which this document's first draft sold as a performance win —
deleted 41 repetitions of something that was not measurably costing anything. It was redundant work
and deserved to go; it was not the bottleneck. **The measured costs are `GeneChoice_J_gene` and
`Deletion_V_gene_Three_prime`**; S4c removes the second of the two.

> **Provenance of the S4b measurement.** Finding 6's `\|R\|`-fold redundancy was confirmed by
> instrumenting the **integration** inference path, which loads the human TCR-α model — a **VJ**
> model with 41 insertion realizations. 438 invocations, the map final after the first pass in every
> one. The reasoning is structural and applies to a VDJ model equally, but the measurement covers the
> VJ case only. (`ce3e4b0`'s commit message still says "the regression model", which is wrong; §6.10
> finding 5 and the comment in `Insertion.cpp` were corrected afterwards.)

### Planned

| step | REFACTOR | GENERALIZE | PERF | CORRECT |
|---|---|---|---|---|
| **S4e** | — | — | hoists the sweep **out of the OpenMP region**: model-only and thread-invariant, so 22 threads currently build 22 copies of one answer. S4c left the profiles owned **by the event**, not by the model, so what S4e owes is making them shareable across `Rec_Event::copy()` — a `shared_ptr<const SpanProfile>` in `JunctionBound`, or a model-side arena. The crude-bound pass stays per-thread because its `updated_proba_bounds_list` points into per-event mutable state | — |
| **S4d** (deferred here) | — | `⊗ᵐᵃˣ`, the max-convolution of two profiles. **Dropped from S4c**: it has no production consumer — every current query is gap-bounded and `Gene_choice(D)` keeps the *retained* decomposition, not a max-folded one — and §2.5's own recommendation is not to add query semantics before a consumer exists. Lands with `⊗ᵉⁿᵘᵐ` in 5b | — | — |
| **S4d** | — | — | Tensor-backed containers for the 3-D structure. **Gated** on the Tensor API, itself blocked on the C++23 bump | — |
| **5b** | — | `⊗ᵉⁿᵘᵐ` — the retained decomposition for the `no_d_align` enumeration, always three components whatever the topology | — | — |
| **R1** | — | — | — | **`Dinucl_markov` creates the insertion segment** (O12 (a′)) — no partially-constructed segment, `int_undefined` leaves constructed sequences, §7.13 dissolves, `SpanAccumulator` deleted |
| **R3** | — | — | — | `Insertion` writes its offsets and **requests the layer it writes**; its mismatch-list defect is *dissolved* by R1's rescope rather than fixed |
| **R3b** | — | — | — | `LayeredArray::set()` requires a prior claim instead of raising it silently |
| **R6** | — | — | — | tighter bound via a within-clique **joint** max. **Changes which scenarios survive pruning**, so it is gated on convergence rather than bitwise regression |

---

## 7. Decided: `Dinucl_markov` will create the segment

*(Quentin, Sep 11 2026 — decision O12. Analysis kept because the *why* outlives the decision.)*

`SpanAccumulator` narrows the `Insertion` ↔ `Dinucl_markov` handshake from a shared pointer into a
sequence map down to a declared publish/read of one integer. But they remain **the only two users**,
so it is still a channel with exactly two ends — and it carries an **unstated assumption**:
`Dinucl_markov` raises its per-nucleotide probability to the published *length*, which is only right
because that length is **entirely undetermined**. Nothing says so.

Two ways out. They are not complementary — the second makes the first vacuous.

### (b) Publish the undetermined count, not the length

The creator publishes *how many positions are still to be chosen*: 0 for a gene template, `n` for an
insertion. `Dinucl_markov`'s factor becomes `p^(undetermined)`, which is the correct general form;
`p^(length)` is today's coincidence.

| | |
|---|---|
| **for** | small and local — one publisher rule, one exponent. No change to `iterate()`, to layers, or to model semantics |
| | bitwise today, since undetermined ≡ length for insertions |
| | moves the invariant out of an assumption and into the published quantity, which is exactly the defect named |
| **against** | keeps the placeholder state, so §7.13's write-through-a-borrowed-pointer and R1 both stand |
| | keeps two events writing one segment; the coupling narrows, it does not go |
| | introduces a distinction that is **identically zero everywhere today** — unexercised generality, of the kind O11 step 3 was deferred to avoid |

### (a′) Move creation to `Dinucl_markov`

`Insertion` decides length and offsets only (declaring `SeqConstructionRole::None` for the
sequence); `Dinucl_markov` **creates** the segment rather than filling it. The partially-constructed
state never exists.

| | |
|---|---|
| **for** | **§7.13 dissolves.** "Dinucl writes through the pointer Insertion stored" stops being a defect when Dinucl owns the segment — R1 becomes unnecessary rather than pending |
| | **one of R3's three defects dissolves too**: Insertion creates no sequence, so it needs no mismatch list for one |
| | `int_undefined` leaves constructed sequences entirely, so *"no undetermined nucleotide at any hand-off"* becomes a **global** invariant rather than a `Fills`-conditional one — and `1794b5f`'s leaf assert strengthens into a hand-off assert |
| | tier 3's `Fills` row, and §2.7's ceiling/floor question for placeholder positions, both stop applying to insertions |
| **against** | changes `iterate()` for two events — the **hot path**, not init — and moves layer claims, so it is not obviously bitwise |
| | the segment is **transiently absent** between the two events. Today nothing runs in that window, because every model orders a Dinucl immediately after its Insertion — but that is a model-ordering fact, not an invariant, and B10 wants absence to be an explicit written value |
| | it **breaks the accumulator's publisher rule**: with `Insertion` declaring `None` for sequence construction, "only creators publish" means it publishes nothing and Dinucl has nothing to read. The rule has to move to *whoever decides the length* — the offsets creator, after R3. Coherent, but a second change |

### Decision: (a′) — and `SpanAccumulator` stays, reinterpreted

**(a′) subsumes (b).** If the creator is also the filler, "undetermined" and "length" are the same
number by construction and (b)'s distinction has nothing left to express. The strongest argument for
(b) would have been a future case of genuine partial determination, and **amino-acid Pgen is not
one** — its ambiguity lives in the *query* (`iupac_union` / `iupac_intersection` / `patches` on
`JournaledQuery`), not in constructed sequences, so gene segments stay fully determined by their
creators.

So (a′) is carried, and it rescopes two repairs rather than adding one: **R1** becomes *"Dinucl
creates the segment"* instead of *"give Dinucl its own layer"*, and **R3** loses its mismatch-list
item — under (a′) `Insertion` creates no sequence, so that `[!shouldfail]` case is **deleted rather
than made to pass**.

**`SpanAccumulator` stays, reinterpreted** *(Quentin, Sep 11 2026, reversing an earlier call to
delete it)*. Its content generalises to

> **the number of nucleotides implied by a segment's offsets for which no constructed sequence
> exists yet** — equivalently, how many nucleotides of that segment are *still to be chosen*.

The two readings coincide on every case that exists: before (a′) an insertion's segment is present
but all placeholder, after (a′) it is absent entirely, and either way the count is `n`. A gene
segment's count is **0** — its nucleotides are fixed by the template.

That phrasing is what earns the object its keep, and it is worth being precise about why. It lets
`Dinucl_markov` state its requirement as *"a segment whose offsets are already placed and whose
sequence is not yet created"* — **a property of the state, not a reference to `Insertion`**. An
explicit lookup from a `Dinucl_markov` to its `Insertion` would be the coupling this refactor exists
to remove; a declared precondition over the state says nothing about which event established it, or
how.

### Why it cannot simply be deleted

Two reasons, and the first is the one that settles it.

**The need does not go away with (a′).** In the fold there are no offsets — only accumulated
lengths — so `Dinucl_markov` still has to learn from somewhere how many nucleotides it will choose.
(a′) removes the *reason* for the link, not the link.

**The keying is load-bearing, not decoration.** On a V→J span **both** insertion/dinucl pairs
participate at once — `affects_proba_of` answers true for `VD_ins_seq` and for `DJ_ins_seq` alike
when the junction is `VJ_ins_seq` — so a single unkeyed scalar would alias the two.

Grouping each `Insertion` with its `Dinucl_markov` in the fold *would* remove the channel, since the
length would be in scope within one group iteration. That remains available, and the pair is a
genuine functional clique, so R6's group mechanism would get a second consumer. But it is an
optimisation of the mechanism, not a reason to drop the concept.

### What changes, and when

| | today | after R1 / R3 |
|---|---|---|
| **publisher** | the segment's *sequence* creator, publishing its length | an event that **creates the offsets but not the sequence**, publishing the length they imply |
| **derivable from** | `get_seq_construction_role` | `get_offset_role` **and** `get_seq_construction_role` — both already exist, so **no new capability is needed** |
| **`Gene_choice` publishes** | its template length — semantically wrong under the new reading, harmless because nothing reads it | **0** |
| **name** | `SpanAccumulator` — names the context, not the content | something content-shaped: `UnfilledSegmentLengths` reads well at the call site (`unfilled.length_of(id)`) |

One thing this exposes that A0 has no vocabulary for: `Dinucl_markov`'s *"offsets placed, sequence
absent"* is a **precondition**, and every A0 query so far states what an event **provides**. The
reserved slot for requirements is D.9's `get_context_dependency()` / `get_context_seq_types()`, so
this is the first concrete consumer for it — worth knowing before Phase D designs that API rather
than after.

---

## 8. What has *not* changed, and is not scheduled to

Worth stating so the scope is not over-read:

- **The bound is still a relaxation.** It maximises each event's conditional probability
  independently, which is looser than maximising jointly over a conditioned clique. R6 tightens it;
  until then the slack is real and §6.10 quantifies where.
- **The bound is still query-independent.** It is built from the model alone and knows nothing about
  the read, which is what lets it be built once per thread rather than per sequence.
- **The maps are still built per consumer**, and redundantly. The cache-and-invalidate scheme that
  would remove that is recorded in the iterate plan §2.5 as a deferred optimisation, deliberately not
  scheduled — it is the *smaller* of the two redundancies. The larger one, every thread rebuilding
  the identical answer, **is** scheduled, as S4e.
- **`initialize_Len_proba_bound` is still virtual** and still selects its map by enum. S4c is what
  changes that; **S4b did not remove the enum ceiling**, and a tandem-D model still throws today.

---

## 9. Where to read further

| | |
|---|---|
| the composition algebra, the frame, boundary-addressed spans | §2.5 of the iterate plan |
| what the legacy machinery is, finding by finding | §6.10 |
| delivered notes per step | §6.11 (S4a), §6.12 (S4b) |
| the repair queue and its gates | §6.9 |
| decisions O8 (S4 scope), O10 (layer ownership), O11 (boundary spans) | §8 |
