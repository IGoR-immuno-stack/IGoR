# Counter rewrite plan: tests, a scenario record, and tandem D

**Created**: Oct 19 2026, split from [ENGINE_MIGRATION_PLAN.md](ENGINE_MIGRATION_PLAN.md) (MIG), whose
lane C it is.
**Scope**: the `Counter` base class ([Counter.h](../src/igor/Model/Legacy/Counter.h)) and its four
subclasses in `src/igor/Inference/Legacy/`; the `Scenario` snapshot they read
([Scenario.h](../src/igor/Model/Legacy/Scenario.h)); their CLI options and output files.
**Siblings**: [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md)
(CAP: T7a's guards, T7b), [ERROR_MODEL_REWRITE_PLAN.md](ERROR_MODEL_REWRITE_PLAN.md) (ERR: the
mismatch role, the error model's scope, the hypermutation models' context), MIG (the scenario
record, the driver). The order across all plans is in
[IMPLEMENTATION_GRAPH.md](IMPLEMENTATION_GRAPH.md).
**Line numbers** are those of `feature/tandemD-T7-pinning` (T1 merged in, T3 not). "On T3" marks
`feature/tandemD-T3-plan`.
**Prefixes**: CT (steps), CX (defects), CR (repairs), A (access patterns), QC (decisions).

---

## 0. Summary

- **None of the four counters has a unit test.** Three are pinned by sorted 6-digit goldens on two
  regression tracks, coverage for V and J only; `Errors_counter` is pinned by nothing.
- **On a tandem-D model** they return wrong values without any error, and D coverage terminates
  the process. T7a's guards will refuse such a model until this plan lands.
- **`Errors_counter` should go** (§2). It is not contained in the coverage counter, but nearly all
  of it is already computed elsewhere: by the error model per read, and by the best-scenarios
  counter per scenario. What remains moves into the per-read summary.
- **The coverage counter's coupling is to a segment's template, not to `Gene_choice`** (§3). Its
  direct reads of event members give way to two records: the realization index of every event,
  and a provenance record per segment, both written at the writer's layer as offsets are today,
  with no strings in the hot loop and only when something reads them.
- **The duplicated `Deletion` semantics can be removed** (§4). The events that build a segment write
  what they did in segment terms, as a provenance record; coverage and the hypermutation models
  read that instead of re-deriving trims and palindromes from deletion values. What assumptions
  remain are declared by each counter as requirements, checked when it is initialised.
- **The outputs move to buffered sinks** (§9): the text files stay byte for byte, a Parquet backend
  writes one row group per chunk, and the best scenarios can also be written as AIRR rearrangements.

---

## 1. Inventory

| Counter | Reached by | Reads | V/D/J hard-coding (code tokens) | On a tandem model | Output | Tests |
|---|---|---|---|---|---|---|
| `Best_scenarios_counter` | `--scenarios N` ([legacy_main.cpp:913](../app/igor/legacy_main.cpp#L913)) | P(s, read); every event's `get_current_realizations_index_vec()` in queue order ([Bestscenarioscounter.cpp:104](../src/igor/Inference/Legacy/Bestscenarioscounter.cpp#L104)); the V, D, J mismatch lists (193) | 6 `Seq_type` | D1/D2 mismatches silently left out; duplicate column names until T2 | `best_scenarios_counts.csv`, once per read, shared stream | goldens: `inference`, `no_d_align` (sorted, rank dropped); `tie_digits = 12`, `keep_ties` only |
| `Errors_counter` | `-run_demo` and `igor-demo` only; no CLI flag | sizes of the V, D, J sequences and mismatch lists ([Errorscounter.cpp:117-140](../src/igor/Inference/Legacy/Errorscounter.cpp#L117)) | 24 `Seq_type` | silently undercounts | `scenarios_background_and_errors.csv`, `sequence_mutation_frequency.csv`, per read, every iteration | none (its sort rules in `test_inference.sh` have no golden) |
| `Pgen_counter` | `--Pgen` (estimator mode only) | `build_full_sequence` (Scenario.h:145-174); P(s) and P(s, read) | 12 `Seq_type`, all without effect on the live path; the snapshot's VDJ/VJ branch | σ assembled as V+J: wrong, silently | `Pgen_counts.csv`, per read | goldens in estimator mode; the per-σ and sequence modes untested |
| `Coverage_err_counter` | `--coverage {V_gene, …, VDJ_genes}` | events by enum, raw slots, the gene choice's realizations map (§3.1) | 24 `Seq_type`, 3 gene class | V and J work while they keep their legacy names; D terminates the process | `{V,D,J}_genes_cov_and_err.csv`, summary per iteration | goldens for V and J; D, `Npoint > 1` and per-read mode untested |

Totals: 89 `Seq_type` and 3 gene-class tokens across the four counters, `Scenario.h` and
`EventUtils.cpp`, as ERR §3.7 counts.

**Lifecycle, as `GenModel` drives it** ([GenModel.cpp](../src/igor/Inference/Legacy/GenModel.cpp)):
- the master copies are initialised once, opening the files (191-204);
- per iteration and per thread, they are copied and initialised against the thread's events
  (306-313, 404-416);
- at every accepted leaf, `count_scenario` (Rec_Event.cpp:435-439);
- per read, `count_sequence` and `dump_sequence_data`, before the mean-error rejection (586-605);
- after the iteration, a merge in thread order and `dump_data_summary` (640-662).

---

## 2. `Errors_counter` and `Coverage_err_counter`: is the first included in the second?

### 2.1 What each computes

| | `Errors_counter` | `Coverage_err_counter` |
|---|---|---|
| Unit | the scenario: n_gen = V + D + J sequence lengths, n_mis = their mismatch counts | the template position of one gene, or N-tuples of positions |
| P nucleotides | counted in n_gen, and their mismatches in n_mis (they are part of the constructed sequence) | discarded (Coverageerrcounter.cpp:298-316, 342-364) |
| Segments | V, D, J | V or J; D terminates (CX1) |
| Per read | E[n_gen], E[n_mis], E[n_mis/n_gen]; the top N scenarios with (P, n_gen, n_mis) | per-gene arrays, in a per-read mode the CLI cannot reach |
| Over the reads | — | Σ over reads of P(s \| read)·[position covered] and ·[position mismatched], per iteration |

### 2.2 Where each part of `Errors_counter` already exists

- **E[n_mis] per read** is already a column of `inference_logs.txt`: `get_seq_mean_error_number()`
  ([GenModel.cpp:579-580](../src/igor/Inference/Legacy/GenModel.cpp#L579), Errorrate.cpp:80-87). It
  uses the same definition of genomic positions (V, D and J sizes and mismatch lists,
  [Singleerrorrate.cpp:128-148](../src/igor/Model/Legacy/Singleerrorrate.cpp#L128)), and the
  hypermutation models fill the same field.
- **E[n_mis/n_gen] per read** is computed by `Single_error_rate` as `seq_weighted_er`
  (Singleerrorrate.cpp:167-174), which guards n_gen = 0 where the counter does not (CX8). It is
  the quantity Single's legacy M-step averages (ERR X13), and it is not printed.
- **E[n_gen] per read** is computed by neither. Summing the coverage over positions gives it for V
  and J only, without P nucleotides, and only in the per-read mode.
- **The top N scenarios with their (P, n_mis)** are what the best-scenarios counter lists, with the
  mismatch positions themselves and a canonical tie order, where this counter breaks ties by
  enumeration order (CX9). n_gen is not printed there.

So `Errors_counter` is **not** contained in the coverage counter: coverage works per template
position, sums over reads, covers V and J only and drops P nucleotides. It is, however, almost
entirely contained in what the error model and the best-scenarios counter already compute.

### 2.3 Recommendation: drop it

- **Delete `Errors_counter`**, its two demo instances (legacy_main.cpp:1773, igor-demo
  main.cpp:271) and the dead sort rules in `test_inference.sh`.
- **The per-read summary gains E[n_gen] and E[n_mis/n_gen]**, next to the E[n_mis] it already has.
  - They come from the error model, which counts k and n for every scenario anyway: in the core of
    ERR §3.5 for Single, and in the context model of ERR §3.6 for the hypermutation models.
  - So there is one definition of a genomic position, the error model's scope (ERR §3.2), instead
    of the counter's V/D/J literal.
  - After L1 they are fields of the driver's `ReadSummary`. Where they are written is QC1.
- **The top-N list with n_gen is dropped.** If n_gen per scenario is wanted later, it is one more
  column of the best-scenarios output, computed from the record (§3).
- **This settles three defects without code:** CX8's division by n_gen, CX9's missing flag and
  golden, and its enumeration-order ties.
- **The analyses it fed** (REV §5.3 rows (a)–(c), LAYER's `ErrorStatsObserver`) then refer to the
  per-read summary.

---

## 3. The coverage counter's coupling, and a generic record

### 3.1 How it reaches into the events

The coverage counter does not just read the scenario. It knows which events exist, which member
holds what, and how the events lay a segment out on the read. The other counters share some of
these patterns:

| # | Pattern | Where | What it actually needs |
|---|---|---|---|
| A1 | "The V, D or J gene choice" and "the V 3′, D 5′/3′, J 5′ deletion", found by (event type, `Seq_type` enum, side) | [Coverageerrcounter.cpp:163-265](../src/igor/Inference/Legacy/Coverageerrcounter.cpp#L163); Pgencounter.cpp:93-99 | the event that creates segment s, and the one that modifies each of its ends |
| A2 | Raw pointers into live event members, dereferenced at every leaf: `alignment_offset_slot()` (→ `&alignment.offset` of the current read's alignment), `current_realization_index_slot()`, `deletion_value_slot()`. Correct only because each thread owns its events and the DFS has just written them; `Gene_choice` publishes the first two only for a gene at an end of the ordering ([Genechoice.cpp:272-277](../src/igor/Model/Legacy/Genechoice.cpp#L272)) | Coverageerrcounter.cpp:164-267 | which gene, where its template sits on the read, how many nucleotides each end lost or gained |
| A3 | Per-event scratch read through a public accessor, in queue order | Bestscenarioscounter.cpp:104 | every event's realization, Dinucl chains included |
| A4 | A copy of the gene choice's realization map, to size the arrays | Coverageerrcounter.cpp:168-174 | the template lengths |
| A5 | `Deletion`'s semantics re-derived (§4.1) | Coverageerrcounter.cpp:282-366 | for each read position of the segment, the template position it shows, and whether directly or mirrored |
| A6 | What "genomic" means, implied: the V, D and J sequences with their P nucleotides | Errorscounter.cpp:117-140; Bestscenarioscounter.cpp:193 | the segments compared with the read |
| A7 | The snapshot's six fixed slots, filled from ids 0 to 5, and the VDJ/VJ branch of `build_full_sequence` | [Scenario.h:57-90](../src/igor/Model/Legacy/Scenario.h#L57), 145-174 | the constructed segments, in layout order |
| A8 | Timing implied by the DFS: a counter runs at the leaf because every layer is written by then; `count_sequence` runs before the error model resets its per-read state | Rec_Event.cpp:435-439; GenModel.cpp:586-605 | a completed scenario; a per-read summary |
| A9 | The read likelihood taken from the error model | GenModel.cpp:588 | L |

### 3.2 Is the coverage counter coupled to `Gene_choice` by design?

The claim: coverage measures something specific to the gene choice (which template, and where), so
coupling to that event is by design. The arguments against it:
1. **It is already coupled to two event types, not one.** `Gene_choice` decides the template and
   where it sits. `Deletion` decides how much of it is kept and whether a palindrome extends it.
   The counter reads both (A2, A5). What it measures belongs to the *segment's construction*, which
   several events share.
2. **The plans already foresee other templated content:**
   - the flank segments (CAP B3, ERR H3), template outside the read;
   - a D inversion, modelled as a content modifier (VDJ review);
   - the P nucleotides themselves, templated content a modifier creates, which CR-P counts.

   A counter keyed to `Gene_choice`'s members follows none of them.
3. **The same mapping has other readers.** The hypermutation models map read positions to template
   positions for their context, and the error model's per-segment bound counts within the same
   segments. Each reader re-deriving it is how the copies of §4.1 arose, and the defects sit in
   the copies.
4. **Another engine could not serve it.** A DP engine computes coverage exactly from posterior
   marginals over (gene, placement, end trims) (REV §5.3), but only if the counter states what it
   needs. Member reads tie it to the DFS that wrote the members.
5. **The new architecture forbids it.** The synthesis has engines read a const model (SYN §5, row
   2) and wants events without per-read state (SYN §4, row 5.6). A slot into a per-thread event copy
   cannot survive that.
6. **It cannot be tested alone.** Without a full model and a DFS to write the members (A8), there is
   nothing to read.

What survives: coverage *is* about genomic templates. A segment with no template has no coverage,
so a requirement that the segment have a templated creator is legitimate (§4.3). That is a
coupling to a capability, not to a class or its members.

### 3.3 What a generic record can hold, and what stays typed

The idea: each event dumps its state for the scenario into something counters read generically.
The requirements are:
- written as the scenario is built;
- read at the leaf;
- no allocation per realization: the DFS visits about 1,200 complete scenarios per TRB read and
  many more partial ones (REV §13.5), each with about ten events;
- introspectable, so a dump can print it.

| Option | Cost in the hot loop | Who interprets it | Verdict |
|---|---|---|---|
| (a) A JSON string per event per realization | an allocation and formatting per realization, parsing per read | the reader, by key | rejected for the transport; JSON is the right format for a *dump* |
| (b) `std::any` per event | an allocation for anything larger than the small buffer | the reader, by casting to a type it must know | rejected |
| (c) A typed struct per event type | none | the reader, which names the event type | the honest form for a reader that really is about one event type |
| (d) Named integer fields that each event declares (`template_offset`, …) | one or two stores per realization | the reader, by field name | rejected: see below |
| (e) **One universal field per event, and a typed record per segment** | one store per realization; a few for the events that build a segment | the event writes, in a vocabulary fixed by the segment | **recommended** (QC2) |

**Why (d) adds nothing over (c).** A named field is a contract: the reader asking the creator of a
segment for `template_offset` must know that this creator has a field of that name, and what it
means. If two events that can create a segment name or define it differently, the reader must
know each event's names, which is (c) with strings for keys. (d) is generic only if the names form
a vocabulary every event shares, with one meaning each. Such a vocabulary exists, but it belongs to
the segment, not to the events, and that is (e).

**(e), what a reader actually needs.** None of today's readers needs event-specific state beyond
two things:
- **The realization record.** One integer per event, the realization index: universal, since every
  event has one.
  - It is MIG L3's scenario record, and the engine's own accumulation reads it at the leaf in place
    of `add_to_marginals` reading event members.
  - A sequence-valued event (`Dinucl_markov`) needs no separate value: its chain of dinucleotide
    indices follows from the content of the segment it fills and the anchor nucleotide.
  - A deletion's number of nucleotides, or a gene's name, is not stored: the reader takes the index
    from the record and the value from the event's domain in the const model.
  - It generalises `current_realizations_index_vec`, which the best-scenarios counter already reads
    for every event without knowing its type (A3).
- **The segment provenance** (§4.2). A fixed struct per segment: template, kept interval,
  palindromic extensions, read interval. Its fields mean the same whichever event wrote them,
  because they describe the segment.

Coverage needs the provenance. Best scenarios needs the realization record. Pgen needs the
constructed sequence. The hypermutation models need the provenance and the flank segments.

**A dump** — the best scenarios with gene names and deletion values, a debug trace — is rendered
off the hot path. Each event type turns a realization index into a description from its domain
(`describe(index) → json`). No event state is copied for it.

**When (c) is right.** If a future counter really is about one event type's internals, it uses a
typed struct that event writes, and it declares the event type as a requirement (§4.3). The
coupling is then explicit and checked, instead of a slot read in passing.

### 3.4 Emitted on demand

Neither record allocates. Their slots are claimed when the model is initialised, in layered
arrays like the offsets, and a realization writes integers at its own layer. What a record costs
is the stores. And the `Scenario` snapshot built at every accepted leaf costs three vectors,
counters or not (CX16).

Emitting only what someone reads removes even that:
- **Declared reads.** At initialisation, every consumer declares what it reads:
  - the counters, through `check_admissible` (§4.3) and later their declared classes (CT5);
  - the error model, its mismatch lists;
  - the engine, the realization record.
- **Resolution.** The model resolves, per event and per segment, whether anything reads each
  record. An event writes a record only if something does, behind a flag set once at
  initialisation, as `publishes_alignment_state_` is today (Genechoice.cpp:769). An unread record
  claims no layer.
- **The realization record is always on.** The engine's accumulation reads it, so it replaces reads
  of event members rather than adding work.
- **The provenance is on only when read:** by coverage, by the hypermutation models (ERR H3), by
  the AIRR output (§9.3) or by a debug dump.
- **The snapshot** is built only when a counter is attached, as a view over the registry-sized maps
  (CT2a), not as three vectors per leaf.
- **Measured.** The cost is measured as CAP D4 measures containers: the five tracks, with and
  without a provenance reader.

---

## 4. The duplicated `Deletion` semantics, and the requirements that remain

### 4.1 The copies

`Deletion` builds a segment's end
([Deletion.cpp:334-391](../src/igor/Model/Legacy/Deletion.cpp#L334)):
- a positive value trims the template;
- a negative value is a palindrome: the |k| nucleotides nearest the trimmed end, reversed,
  complemented and put back on the far side, with their mismatches against the read appended.

A 3′ palindrome nucleotide i therefore mirrors template position t3 − i, and a 5′ one mirrors
t5 + (k − 1 − i).

The same semantics is written five times:

| Copy | Where | Known defects in the copy |
|---|---|---|
| `Deletion::iterate` | Deletion.cpp:334-391 | — (the owner) |
| `Deletion::construct_realization` | Deletion.cpp:536-590 | — (ERR §3.10 already plans to share it with the first) |
| The coverage counter: V's covered range `max(0, −offset) .. len − max(0, del)`, and P-nucleotide mismatches told apart by position; J's mirror image | Coverageerrcounter.cpp:282-366 | CX3: J's end bound counts the palindrome, and can write past the array |
| `Hypermutation_global_errorrate`, the V head and J tail from `v_3_del` and `j_5_del` | Hypermutationglobalerrorrate.cpp:539-945 | ERR X8 (J tail index), X9 (V head index) |
| `Hypermutation_full_Nmer_errorrate`, the same | HypermutationfullNmererrorrate.cpp:591-881 | the same |

The copy also builds in a second assumption: that a segment end has **one** modifier, and that it
is a `Deletion`. The `Events_map` key cannot express more than one modifier per end anyway (ITER
§9).

### 4.2 Removing the duplication: the event that changes a segment says what it did

The semantics stays with its one owner: **each event that builds or changes a segment writes, at
its layer, what it did in segment terms**. Readers consume that and never see a deletion value.
- **A segment provenance record**, one per segment, layered like the offsets:
  - the template: library and realization index;
  - the kept template interval [t5, t3];
  - the palindromic extension at each end, p5 and p3;
  - the read interval [r5, r3].
- **Who writes what:**
  - `Gene_choice` writes the template, its placement and the full interval, already clipped at the
    read where the template overhangs it (`clip_template_before_read_`, `_after_read_`);
  - `Deletion` writes t3 (or t5) for a trim, and p3 (or p5) for a palindrome;
  - a future modifier writes the same fields by its own semantics.
- **One reader-side helper, in `src/igor/Model/`** (engine-agnostic, ERR §3.9's rules), maps each
  read position of the segment to (template position, direct or mirrored), or to "not templated".
  It serves coverage, the hypermutation models' context (ERR H3, with the flank segments) and the
  error model's per-segment bound.
- **The shared construction function** that ERR §3.10 plans for `Deletion` returns the content and
  the provenance update together, so `iterate()` and `construct_realization()` cannot disagree.
- **Bitwise.** The values are the same. Today's coverage is the provenance with the mirrored
  positions dropped, and CR-P then adds them, deliberately.
- **Cost.** A few integer writes per realization of the events that build a segment, and none when
  nothing reads the provenance (§3.4).

With it, the five copies become one owner and one reader, CX3, X8 and X9 disappear with their
copies, and several modifiers on one end compose, each updating the record in turn.

### 4.3 What a counter still assumes, and how to check it

Even with the provenance, a counter's numbers mean something only on models with the right shape.
Coverage needs:
- a segment with a template: one content creator with a template library;
- every event that changes that segment writing the provenance;
- offsets created for it.

Today, before the provenance exists, it assumes much more: exactly one `Gene_choice`, at most one
`Deletion` per end, V at the left end of the read and J at the right. Those assumptions are
unchecked: D ends in `std::terminate` (CX1), and a missing gene choice in a null dereference
(CX2).

**Requirements, declared by the counter and checked at initialisation:**
- **Interface.** `Status check_admissible(const ModelView &) const` on the counter. `ModelView` gives
  the registry, `SegmentIndex` (on T3) and the events' capability answers. The shape is that of
  `TraversalPlan::check_admissible` (on T3) and of the engine's own check (MIG L1). A failure
  refuses the run, with a message naming the counter, the segment and the requirement, which is
  what T7a's guards do by hand today.
- **A small vocabulary.**

  | Requirement | Checked through |
  |---|---|
  | `templated(s)`: exactly one content creator of s, and it has a template library | `SegmentIndex::content_creators(s)`, the creator's domain |
  | `end_modifiers(s, side) ⊆ K`, with at most n of them | `SegmentIndex::end_modifiers(s, side)`, their event types |
  | `writes_provenance(s)`: every creator and modifier of s writes the provenance | a capability query on each event |
  | `anchored(s, side)`: s is at that end of the ordering, so that end is the read's | `registry.ordering()` |
  | `has_mismatch_writer(s)` | the mismatch role (ERR §3.3) |

- **Coverage before the provenance (from CT2a):** `templated(s)`; `end_modifiers(s, side) ⊆
  {Deletion}`, at most one; `anchored(s, Five_prime)` for V-style geometry or
  `anchored(s, Three_prime)` for J-style. A tandem model's D1 is then refused cleanly instead of
  crashing.
- **Coverage after the provenance (from CT2b):** `templated(s)` and `writes_provenance(s)`. No event
  type is named, and nothing requires the segment to be anchored.
- **Check levels.**
  - C1 at initialisation, as above.
  - C2 at the leaf, in debug builds: the provenance's length equals the length of the segment's
    constructed sequence, which ties the provenance to the content it describes.
- **The other users.** The same mechanism carries T7a's guards for the other counters and for the
  hypermutation models (ERR X17). It is also how a counter will declare its needs to an engine
  other than the DFS (CT5).

---

## 5. Defects and oddities

Categories follow ITER §6.1: *confirmed* defects are pinned at their intended value under
`[!shouldfail]` and fixed once, in a CR commit; *odd* behaviour is pinned as observed.

| # | Finding | Category | Owner |
|---|---|---|---|
| CX1 | D coverage throws a `const char*` ([Coverageerrcounter.cpp:321](../src/igor/Inference/Legacy/Coverageerrcounter.cpp#L321)); `GenModel`'s `catch (exception&)` does not catch it, it leaves the OpenMP region, and the process calls `std::terminate` | confirmed | CR; then CT3 implements D |
| CX2 | `try_get_event`'s result is ignored (163, 198, 242): D coverage on a VJ or tandem model dereferences null | confirmed | CT2a's requirements refuse it |
| CX3 | J coverage with a palindromic J (`j5del < 0`): the end bound `max(0, j5del) + (J3′ − J5′ + 1)` (347-351) includes the \|d\| P nucleotides, so it over-counts \|d\| positions at the 3′ end and, when the read covers the whole J, writes past the array. The TRB goldens do not reach it (their reads cover about 13 nt of J) | suspected; confirm in CT0 | CR, early: it can corrupt memory on long reads |
| CX4 | The coverage summary is merged in thread order (GenModel.cpp:640-657, its FIXME): its last bits depend on the schedule; the 6-digit output hides it | confirmed | CT4 |
| CX5 | `initialize` re-allocates without freeing; a destructor over raw arrays with the implicit copy left in (rule of three); `throw except;` slices; messages copied from the hypermutation model | confirmed | CR |
| CX6 | Pgen's sequence mode writes a header and no rows; its per-σ mode does not print σ (Pgencounter.cpp:181-211) | confirmed | CR, QC6 |
| CX7 | Pgen's V/D/J flags are computed and read only by the dead path | dead | CT1 |
| CX8 | Unguarded divisions: by n_gen (Errorscounter.cpp:176), by L (267-269; best scenarios' `count_sequence`) | confirmed | the first goes with `Errors_counter` (CT1b); CR for best scenarios |
| CX9 | `Errors_counter` has no CLI flag and no golden; its top-k breaks ties by enumeration order | gap | CT1b |
| CX10 | Per-read outputs carry no iteration column: with `last_iter_only = false`, rows of every iteration are indistinguishable | odd | CT4 |
| CX11 | Counters run before the mean-error rejection (GenModel.cpp:586-605), so rejected reads are counted | odd | QC5 |
| CX12 | Counters see the scenarios above the *running* maximum, so the set depends on the order the DFS explores | inherent to the DFS | documented; a DP engine differs (REV §5.3, check 1 at threshold 0) |
| CX13 | Viterbi: counters accumulate every successive running maximum, so coverage weights exceed 1 | confirmed | CR, with ERR E2 |
| CX14 | Best scenarios: `initialize` is not idempotent (columns duplicate); P is truncated from `long double` to `double` | confirmed | CR |
| CX15 | The deprecated interface (`initialize_counter`, the 7-argument `count_scenario`) and the four legacy bodies have no caller | dead | CT1 |
| CX16 | A `Scenario` snapshot (three vectors) is allocated at every accepted leaf, counters or not | perf | CT2b |
| CX17 | Thread copies share one stream; rows come in thread-finish order, which the tests sort away | odd | CT4 |
| CX18 | Coverage discards the palindromic nucleotides instead of counting them at their mirrored template position | improvement | CR-P |

---

## 6. What each counter needs from an engine

| Counter | Statistic class (REV §5.3) | An engine other than the DFS would serve it by |
|---|---|---|
| Coverage | posterior expectation over one segment's state (gene, placement, end trims) | forward–backward marginals |
| Best scenarios | k-best, with its canonical tie order (40-bit mantissa, then lexicographic on realizations, then mismatches; Bestscenarioscounter.cpp:68-98) | k-best (max, ×) with backtracking |
| Pgen | an arbitrary function of the whole sequence | enumeration or posterior sampling |
| The per-read summary, once `Errors_counter` is gone: E[n_gen], E[n_mis] (additive), E[n_mis/n_gen] (non-linear) | — | the expectation semiring; state augmentation for the ratio |

**The contract a counter needs**, whatever the engine:
- `initialize` once per worker, against a const model, with `check_admissible` (§4.3);
- a per-scenario sink receiving P(s), P(s, read), the record, the provenance and the mismatch lists;
- the threshold semantics: it sees exactly the scenarios that make up L;
- a per-read hook with the read index, the iteration, L and whether the read was accepted;
- per-chunk accumulation with an ordered merge, for summaries;
- a summary hook per iteration, honouring `last_iter_only`;
- an output sink that writes rows in read order.

---

## 7. Generalisation for tandem D

Each counter moves from V, D and J to sets resolved from the ordering. VDJ outputs stay unchanged
until CR-P.
- **Best scenarios**: the mismatches of every segment with a mismatch writer, in ordering order (in
  VDJ: V, D, J, as today).
- **The per-read summary**: n_gen and n_mis over the error model's scope.
- **Pgen**: σ assembled over `registry.ordering()`.
- **Coverage**:
  - keyed by segment, not gene class: `--coverage` accepts segment names, and `D_gene` keeps
    meaning every D segment (QC3);
  - D implemented;
  - geometry from the provenance, for any segment, whichever end anchors it;
  - **palindromic nucleotides counted at their mirrored template position (CX18, CR-P).** A
    template position seen directly and again through its palindrome is counted twice; a mismatch on
    a P nucleotide is an error at the mirrored position. Coverage then reads "expected observations
    per template position", not "probability of being observed", and the output header says so.

---

## 8. Steps

| Step | Needs | Content | Bitwise |
|---|---|---|---|
| CT0 | milestone 1 | Characterization suite at the counter API: `Scenario` and contexts built by hand, no DFS needed. One tiny end-to-end model per counter. VJ and tandem fixtures. Best scenarios: tie settings, viterbi. Pgen: all three modes. Coverage: D (`[!shouldfail]`), `Npoint > 1`, per-read mode, palindromic V and J (CX3). `Errors_counter`'s per-read values on the tiny model, as the oracle for CT1b. Thread-count determinism at full precision (`[!shouldfail]` until CT4). Coverage report and mutations (ITERATE_TEST_GUIDE §7) | tests only |
| CT1 | CT0 | Dead code out (CX7, CX15) | yes |
| CT1b | CT0; MIG L1 for `ReadSummary`, softly | `Errors_counter` deleted (§2.3); the per-read summary gains E[n_gen] and E[n_mis/n_gen], checked against CT0's oracle | the per-read output moves once (QC1) |
| CT2a | CT0, T3 | T7b's counter half, on the legacy DFS: the snapshot sized from the registry and indexed by `SeqTypeId` (A7); the events a counter reads resolved through `SegmentIndex` rather than by enum (A1); best scenarios' mismatches in ordering order; `check_admissible` with coverage's interim requirements (§4.3), which take over T7a's guards for the counters | yes |
| CT2b | CT2a, MIG L3 | The record (§3.3) and the provenance (§4.2): `Gene_choice` and `Deletion` write the provenance; the counters read the record and the provenance instead of the slots, the scratch reads, the copied realization map and the re-derived geometry (A2–A5); coverage's requirements take their capability form | yes |
| CT3 | CT2a and ERR E3 for best scenarios and Pgen; CT2b for coverage | §7, which lifts the remaining guards counter by counter | yes on VDJ; new tandem cases |
| CT4 | MIG L1, L4a | The lifecycle in the driver: the per-read hook takes a `ReadSummary`; an output sink ordered by read index; summing counters in `AccumulatorSet` chunk slots, sparse (QC4); `last_iter_only` applied by the driver | yes, except the coverage summary's last bits (CX4), measured |
| CR | CT0 | One repair per commit: CX1, CX3, CX5, CX6, the best-scenarios half of CX8, CX13 (with ERR E2), CX14 | the outputs each moves, named |
| CR-P | CT2b, CT3 | Palindromic positions counted (CX18) | the coverage goldens move once, named |
| CT6 | CT4 | The row sinks (§9.2): schemas per output, per-chunk buffers flushed in chunk order, the text backend byte for byte, the Parquet backend in Streaming, `output.format` | yes for text; Parquet new tests |
| CT7 | CT2b, CT3, CT6 | The best scenarios as AIRR rearrangements (§9.3), optional, in TSV or Parquet; the Streaming AIRR writer gains the missing columns | new output, new tests |
| CT5 | a second engine | Each counter declares its statistic class and what it reads, for the engine's admissibility check; `Counter` becomes `Observer` (LAYER §5) | MIG §7 |

What the migration needs from these steps, and the reverse: CT0 before MIG L1; CT2b with L3; CT4
after L4a; CT3 and CT4 before L8. CT6 and CT7 follow CT4 and need nothing else from the migration.

---

## 9. Output: buffered sinks, Parquet, and AIRR

### 9.1 Today

- **The files.** Each counter opens its own `ofstream`, shared by every thread's copy, and writes
  semicolon-separated text in a critical section, once per read.
- **What follows from that:**
  - rows come in the order threads finish, which the tests sort away (CX17);
  - numbers carry the default 6 significant digits;
  - per-read files have no iteration column (CX10);
  - lists are packed into strings such as `(c0,c1,…)` and `(m1,m2,…)`;
  - one mode writes a header and no rows (CX6).
- **What exists to build on.**
  - Arrow, Parquet and sparrow are already required dependencies
    ([CMakeLists.txt:125-127](../CMakeLists.txt#L125)).
  - `ParquetWriter::write_batch` writes one record batch to a whole file through
    `parquet::arrow::FileWriter` ([ParquetWriter.cpp:79-121](../src/igor/Streaming/ParquetWriter.cpp#L79)).
  - The AIRR writer serialises the aligner's alignments, not scenarios, with 18 of the schema's
    columns and nothing for a second D (AIRRRearrangementWriter.h, AIRRCommon.cpp).

### 9.2 Sinks

**A row sink per output.**
- At initialisation, a counter declares the schema of each file it writes: column names and types
  (integer, double, string, list of integers or doubles).
- Rows are appended into a per-chunk buffer, and the driver flushes the buffers in chunk order.
  This is the ordered output sink of CT4.
  - Rows come in read order, whatever the schedule (CX17).
  - There is no critical section per read.
  - Every per-read row carries the iteration (CX10).

**Two backends behind the sink.**
- **Text, the default**: byte for byte today's files, so no golden moves. The header-only mode is
  fixed by CR (CX6).
- **Parquet**: one file per output, one row group per chunk, written in chunk order, so the file
  is deterministic.
  - Doubles are stored in full precision.
  - Lists are native, for example `list<int32>` for mismatch positions and for each event's
    realization.
  - Columns compress well, and polars, pandas and arrow read them directly.
  - It needs an incremental writer: a `FileWriter` kept open, with one `WriteRecordBatch` per row
    group, beside the whole-file one.

**Where the code lives** (QC8).
- The sink interface and the text backend live with the counters, in Inference.
- The Parquet backend lives in Streaming, and the application chooses it. Inference links only
  Model today ([src/igor/Inference/CMakeLists.txt:79-82](../src/igor/Inference/CMakeLists.txt#L79));
  injecting the backend avoids an Inference → Streaming dependency, and lets the unit tests use an
  in-memory sink.

**The CLI.** One `output.format` key (`text`, `parquet`) beside the existing `output.*` keys
([app/igor/main.cpp:258-264](../app/igor/main.cpp#L258)).

**Schemas.**

| Output | Columns |
|---|---|
| Best scenarios | `seq_index`, `iteration`, `rank`, `proba_cond_seq` (double), one `list<int32>` column per event named by its nickname, `mismatches` (`list<int32>`) |
| Pgen | `seq_index`, `iteration`, `pgen_estimate`; in the sequence mode `sequence`, `pgen`, `p_joint` (QC6) |
| Coverage | `iteration`, `segment`, `gene_index`, `gene_name`, `coverage` and `errors` as `list<double>` over positions (QC10) |
| The per-read summary (CT1b) | `seq_index`, `iteration`, the likelihood, E[n_gen], E[n_mis], E[n_mis/n_gen] |

### 9.3 The best scenarios as AIRR rearrangements

A second, optional output of the best-scenarios counter: one AIRR Rearrangement row per recorded
scenario, in TSV, or in Parquet with the same column names. The difficulty is that every column has
to be read off the model, not off the counter's own state.

**Which segment is which.** The AIRR roles come from the ordering:
- the templated segment at the left end is V, the one at the right end J;
- the internal templated segments, in order, fill `d_call` and `d2_call`;
- the junctions, in order, fill `np1`, `np2` and `np3`.

In a VJ model `d_call` is null and `np1` runs from V to J. The schema has exactly two D slots, so a
model with a third templated internal segment is refused (§4.3).

**From the record and the provenance (CT2b):**

| AIRR fields | Source |
|---|---|
| `v_call`, `d_call`, `d2_call`, `j_call` | the gene name: the realization index in the record, the name from the event's domain |
| `*_sequence_start`, `*_sequence_end` (1-based, closed) | the provenance's read interval, without the palindromic extensions |
| `*_germline_start`, `*_germline_end` | the provenance's kept template interval |
| `p3v_length`, `p5d_length`, `p3d_length`, `p5d2_length`, `p3d2_length`, `p5j_length` | the provenance's palindromic extensions |
| `np1`, `np2`, `np3` and their lengths | the junction's inserted sequence plus the flanking P nucleotides: AIRR counts them in the combined N/P region, where IGoR's segments carry them inside the gene segment, so they are moved across |
| `n1_length`, `n2_length`, `n3_length` | the junction's inserted length alone |
| `sequence_alignment`, `germline_alignment` | the read from V's start to J's end; the germline assembled from the kept templates and the palindromes, with `N` at untemplated positions. A scenario places each template at one offset and the events read no alignment indels, so the two have equal length |
| `v_cigar`, `d_cigar`, `d2_cigar`, `j_cigar` | ungapped, from the provenance, following the schema's CIGAR conventions (the scenario's alignment, not the aligner's) |
| `junction`, `junction_aa`, `productive` | the CDR3 anchors (`CDR3SeqData::v_anchor`, `j_anchor`, loaded with `-set_CDR3_anchors`) mapped through V's and J's provenance; null when no anchors are loaded |
| `sequence_id`, `sequence`, `rev_comp` | the read; `rev_comp` false |

**Several scenarios per read.**
- The schema has no rank or probability field, and identifies a rearrangement by its
  `sequence_id`.
- Each recorded scenario is a row with the read's `sequence_id`, plus two extra columns:
  `igor_scenario_rank` and `igor_scenario_proba` (P(s | read)). The schema allows extra columns
  whose names collide with none of its fields.
- Tools that expect one row per sequence take rank 1, and a `top 1` option writes only those rows
  (QC9).

**What it needs.**
- The provenance, so the output is reachable only from CT2b on.
- The tandem generalisation (CT3), for `d2_call` and `np3`.
- The sinks (§9.2).
- The Streaming AIRR writer extended with the missing columns.

**Tests.** Read back with the repository's `AIRRRearrangementReader`, plus a check of the
required columns. Then a VDJ, a VJ and a tandem fixture, each pinned against hand-computed rows.

---

## 10. Decisions

| # | Question | Recommendation |
|---|---|---|
| QC1 | Where the per-read summary's new values go | appended to `inference_logs.txt`, after the existing columns, with `test_inference.sh` extended to keep them: one per-read file instead of two, and one named golden move. The alternative is an opt-in per-read file, which moves nothing |
| QC2 | How events expose their state per scenario | §3.3 (e): the realization index of every event, and the provenance per segment, emitted only when read (§3.4); a typed struct (c) for a reader genuinely about one event type, with its requirement declared |
| QC3 | Coverage on a tandem model | per segment, one output per segment; `D_gene` means every D segment |
| QC4 | Summing counters under the chunked reduction | chunk slots with sparse per-gene arrays, allocated on first touch: the FIXME's objection is the size of a dense copy per chunk |
| QC5 | Counters and rejected reads (CX11) | per-read outputs for every read, summaries over accepted reads only; to decide |
| QC6 | Pgen's two broken modes (CX6) | fix the sequence mode (rows, with σ); drop the per-σ mode that omits σ |
| QC7 | The default output format | text, so that no golden moves; Parquet by `output.format` |
| QC8 | Where the Parquet backend lives | in Streaming, injected by the application; the sink interface and the text backend in Inference |
| QC9 | Several scenarios per read in the AIRR output | one row per scenario, with `igor_scenario_rank` and `igor_scenario_proba`; a `top 1` option for tools that expect one row per sequence |
| QC10 | Coverage in Parquet | one row per (iteration, segment, gene), positions as list columns; a long format (one row per position) is the alternative, at about 300 times the rows for V |
