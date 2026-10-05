# Legacy generator rewrite plan: draw, construct, and no V/D/J names

**Created**: Oct 4 2026
**Scope**: the legacy generator: `GenModel::generate_unique_sequence()` and the four
`Rec_Event::draw_random_realization()` overrides. `FastGenerator` (`generate.fast`) is out of scope,
and so is T6b, which waits on it.
**Parent**: T4 in [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md)
(option (a): the legacy generator comes off the `Seq_type` enum, bitwise).
**Counterpart of**: [ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md), whose method
it follows. That plan made `iterate()` generic; this one does the same for generation.
**Toward**: the synthesis' sampling path ([ARCHITECTURE_SYNTHESIS.md](ARCHITECTURE_SYNTHESIS.md) §5,
"Data flow: generation"), in which an engine samples realization indices and each event's `apply`
builds the sequence from them.

---

## 0. Summary

Each `draw_random_realization()` override does four things in one loop:

1. walks the CDF of its marginal row, in `unordered_map` order, until `prob_count >= rand`;
2. writes its piece of the sequence;
3. records the realization's index;
4. moves its children's rows in `index_map` by `index × stride`.

Steps 1, 3 and 4 are the same code in `Gene_choice`, `Deletion` and `Insertion`. **Every V/D/J
switch is in step 2**, the construction. The plan therefore:

- **splits** the step into a *draw* (indices from marginals and the RNG) and a *construction*
  (sequence from indices, with no RNG, marginals or `index_map`). The construction is the
  synthesis' `Event::apply`, and the half a `SamplingEngine` will later drive with
  `SampledEvent::indices`;
- moves the generator onto a **`SeqTypeId`-indexed container** from the same family as inference's
  `constructed_sequences`;
- makes each **construction generic**, one event type per commit, as the iterate plan did for
  `iterate()`.

Every step is bitwise on all five regression tracks. A defect found in construction is pinned by a
`[!shouldfail]` test that asserts the intended value, carried unchanged through the refactor, and
fixed in the R phase (§5), once, with the outputs it moves named.

## 1. What the code does today

| Event | Construction (step 2) | Draw (step 1) |
|---|---|---|
| `Gene_choice` (Genechoice.cpp:713) | `switch (event_class)`, three arms, each writing the template at V, D or J | categorical |
| `Deletion` (Deletion.cpp:545) | `switch (target_seq_type)`, four arms. The V and J arms ignore `event_side`. Scratch strings are `mutable` members of a `const` method | categorical |
| `Insertion` (Insertion.cpp:303) | `insertion_seq_type_str_to_enum()`; a name the enum does not know writes nothing | categorical |
| `Dinucl_markov` (Dinuclmarkov.cpp:311) | throws unless `get_junction().legacy_enums_valid`; fills the insertion's `'I'` placeholders | a chain: one uniform per placeholder, from `dinuc_proba_matrix` (tier 0b), seeded by the anchor's end nucleotide |

- **Assembly**: `V + VJ + VD + D + DJ + J`, by name (GenModel.cpp:860).
- **Dinucl_markov's draw reads constructed content**: the anchor nucleotide, and how many
  placeholders there are. It is the only draw that does. A pure "sample every index, then build"
  pass cannot draw it before construction. tk's `SamplingHandler::sampleSequence(rng, first_state,
  n_steps)` already has the right shape for it, with `first_state` the anchor nucleotide and
  `n_steps` the insertion length.
- **The queue order is a construction order.** Deletions must follow their gene choice, and a
  chain must follow its insertion and the deletions of its anchor. Legacy priorities give that
  order, and nothing checks it.
- **`txt2marginals()` reads values by position, not by header.** It reads the values in sequence
  into the flat array, in model-queue order, ignoring each block's `@nickname`. A marginals file
  written in another order loads silently onto the wrong events. Found while writing G0's fixtures;
  outside this plan's scope, recorded for the model-document work (synthesis §15).

## 2. The split

On `Rec_Event`:

| Member | Kind | Content |
|---|---|---|
| `draw_random_realization(...)` | **non-virtual** | `draw_realization` → `propagate_realization` → `construct_realization`; returns the indices, as today |
| `draw_realization(marginals, index_map, const segments&, rng)` | virtual | The base implements the categorical walk once. The walk is a pure function of the uniform, which makes D4 testable. `Dinucl_markov` overrides it to draw its chain |
| `propagate_realization(indices, index_map, offset_map)` | virtual | The default is the loop that was copied in three places. `Dinucl_markov` overrides it with nothing: a chain is sequence-valued and has no single index to condition a child on, and it never moved one |
| `construct_realization(indices, segments)` | pure virtual | No RNG, marginals or `index_map`. The later `apply` |

The indices are an index list shaped like `SampledEvent::indices`. Construction finds a realization
by index with `realization_at()`, which scans the realizations. They are keyed by name, and a cache
by index would have to stay coherent through `copy()` and through the two constructors that insert
realizations directly. That is not worth it for a scan this size. An index-ordered domain is the
minimal Event's (REVIEW §7.1).

A chain's draw returns one entry per placeholder, with `kNoRealization` where the walk chose
nothing, so that construction puts each nucleotide back at the position it was drawn for. The
realizations file still records only the nucleotides drawn, as before.

## 3. The container

Generation uses the container family behind inference's `constructed_sequences`, rather than a new
type, so that the later shared `ScenarioContext::Layer` (REVIEW §9.2) does not have to absorb a
second "segments by `SeqTypeId`" type.

| | Inference `constructed_sequences` | Generation |
|---|---|---|
| Type | `DynamicSequenceMap<Int_Str*>` | `DynamicSequenceMap<std::string*>` |
| Layers | one per claiming event, for DFS backtracking | one: no backtracking |
| Values | pointers into event-owned `Int_Str` buffers | pointers into one `std::string` buffer per seq type, owned by the generator loop. A creator points the map at its buffer; a modifier edits the buffer in place |
| Lifetime | per thread, reused across reads | a fresh map per generated sequence; the buffers can be reused |
| "Not drawn" vs "empty" | `exists()` vs `occupied()` | every key's layer is requested when the map is built, so a segment nobody created does not `exist()`, and reading it throws, as `.at()` does today |

The map and the buffers it points into travel together in `GenerationState`
(Model/Legacy/GenerationState.h), the synthesis' "generation state" an event's `apply` writes into.
It is a thin owner, not a second container: the segments are the `DynamicSequenceMap`. Its surface
says what an event does to a segment: `create` (a constructing event), `modify` (in place;
throws when nothing has created it), `read`, `built`, and `assemble` in `registry.ordering()`.

Pointers keep the stored values small and trivially copyable (LayeredArray.h:158), and construction
copies nothing. Copies would cost little today, because the legacy draw dominates: it re-sums the
row and does a string hash on `index_map.at(get_name())` for every realization it visits. They
would matter once an indexed draw drives the construction (§6).

**B10 is not a prerequisite.** Generation never walks past an empty segment:
- a chain reads its anchor at its fixed ordering neighbour (`get_junction().anchor_id`);
- an empty anchor throws (D6), as inference does after R2.

## 4. Defects and oddities

Found in G0, against the unmodified code. "Inference" is `Deletion::iterate` on the same segment;
G0 compares the two directly.

| # | Finding | Inference | Generation | Status |
|---|---|---|---|---|
| D1 | V 5′ / J 3′ deletion | trims the side the model names | the V arm always trims 3′, the J arm always 5′ | **fixed** (R): generation trims the side the model names. No model has either deletion, so nothing moved |
| D2 | deletion past the segment's end; palindrome longer than its template | drops the realization | a 5′ trim clamps; a 3′ trim throws `std::out_of_range` | **fixed** (R): refused on either side, as inference drops it. No trained model reaches it (§5.8) |
| D3 | V or J deleted to empty | keeps one nucleotide: drops the realization | an empty segment | **fixed** with D2: a segment on an end of the sequence keeps one nucleotide, read off the ordering as inference reads it |
| D4 | CDF walk boundaries | — | `>=` picks a zero-mass realization at u = 0. When the row's mass stays below u, nothing is written, `()` is recorded, and the next reader throws. G0's fixture hit the second case for real: a misordered marginals file gave an insertion an all-zero row, and its chain then threw on the missing segment | **left** to the sampling handlers: a `FIXME` in the code, the two `[!shouldfail]` cases stay (§5.8) |
| D5 | an insertion with no Markov chain | — | its `'I'` placeholders reach the output | **kept**: sound as is; the parent plan's Phase C may refuse such a model at load |
| D6 | an empty chain anchor (a fully deleted D1 anchoring D1D2) | throws (§7.12, R2) | throws `out_of_range` from `substr` | consistent; waits on B10, which decides for both |

None of these can move the regression corpus. The `generate` track's model is inferred, so it puts
no mass on scenarios inference drops, and no shipped model has a V 5′ or a J 3′ deletion.

## 5. Steps

Gate for every step: `pixi run build`, `test_unit`, `test_integration`, `test_regression`,
`test_generate`, and `test_convergence` on every step that touches source. The convergence round
trip generates through this path (test_inference.cpp:513, test_ModelInference.cpp:484).

| Step | Content | Bitwise? |
|---|---|---|
| **G0** | This document; characterization suites against the unmodified code (§5.1) | n/a: tests only |
| **G1** | The split of §2. The test stubs that override `draw_random_realization` move to `construct_realization`. Container unchanged | **yes** |
| **G2** | `GenerationState`: a `DynamicSequenceMap<std::string*>` with the buffers it points into, and assembly in `registry.ordering()`. No caller | **yes** |
| **G3** | The generator adopts the container. The remaining switches address it through the legacy ids pinned to `Seq_type`. Assembly by ordering, which gives the same string for VDJ and VJ, since no segment outside the ordering is ever written. An empty ordering throws | **yes** |
| **G4a** | `Gene_choice` writes at `seq_type_id`; then `Insertion`. `insertion_seq_type_str_to_enum` goes if nothing else calls it | **yes** |
| **G4b** | `Deletion`: one body on `event_side`. D1 is carried generically: a segment with no left neighbour is trimmed at its 3′ end and one with no right neighbour at its 5′ end, whatever the side. This is B11a's anchored-segment boolean. D2's clamp/throw asymmetry falls out of the shared body unchanged | **yes** |
| **G4c** | `Dinucl_markov` reads at `get_junction()`'s ids. `legacy_enums_valid`, the spec's `Seq_type` fields and `kLegacySeqTypeCount` retire. The tandem end-to-end case turns green | **yes** |
| **R** | ✅ D1 (trim by side). ✅ D2 and D3 (refused, as inference drops them). D4 left to the sampling handlers, D5 kept, D6 to B10 (§5.8). One commit each, each removing its tag | each names what it moves; bitwise expected on the corpus |

Each G4 commit removes the `[!shouldfail]` tags of the tandem cases it makes pass.

### 5.1 G0, delivered

- **`tst/igor/Model/Legacy/test_generation_construction.cpp`**: 15 `TEST_CASE`s organized by
  pattern, not by gene.
  - Every draw puts all of an event's mass on one realization, so nothing depends on the
    `unordered_map` walk order.
  - RNG consumption is pinned by comparing the generator against a copy advanced by the expected
    number of uniforms.
  - Segments are addressed by name through one fixture class, the only code that knows the
    container, so that G3 changes the fixture and not the cases.
  - Patterns:
    - the categorical draw;
    - parent propagation;
    - template write;
    - trim, including D2 and D3;
    - palindrome;
    - generation against inference over every deletion value inference keeps, on all four shipped
      arms, plus the premise of D1;
    - D1 (`[!shouldfail]` ×2);
    - insertion placeholders;
    - the Markov fill: 3′ and 5′ anchors, VJ, empty insertion, placeholder mask, D6;
    - tandem V-D1-D2-J (`[!shouldfail]` ×5: gene choice, insertion, deletion, chain, end to end).
- **`tst/igor/Inference/Legacy/test_generation.cpp`**: 5 `TEST_CASE`s through
  `GenModel::generate_sequences`:
  - point-mass VDJ, with J conditioned on V, so the row a parent selects is pinned end to end;
  - point-mass VJ;
  - D5;
  - seeded pins of the shipped TRA and TRB models (20 sequences each, FNV-1a over both output
    files). **Nothing covered VJ generation before**: the `generate` track is VDJ only.
- Results on the unmodified code: 13 cases pass, and the 7 tagged ones fail as expected.
- Mutations: see §5.2.

### 5.2 G0 mutation results

Eight mutations of the unmodified generator, each built and run against both suites:

| Mutation | Caught by |
|---|---|
| a D 5′ deletion trims the 3′ end | trim and inference-agreement cases; three whole-sequence cases |
| a palindrome not reversed | palindrome and inference-agreement cases |
| a palindrome not complemented | palindrome and inference-agreement cases |
| a chain seeded from the wrong end of its anchor | the Markov fill; three whole-sequence cases |
| a 5′-anchored chain not written reversed | the Markov fill; two whole-sequence cases |
| a chain overwriting non-placeholder characters | the placeholder-mask section only, the one written for it |
| a gene choice that moves no child row | parent propagation; four whole-sequence cases |
| assembly with D before VD | three whole-sequence cases |

All eight caught; the sources were restored and rebuilt after the run.

### 5.3 G1, delivered

- `draw_random_realization()` is non-virtual on `Rec_Event`: draw, then construct, then
  propagate. The four overrides became `construct_realization()`. `Dinucl_markov` also overrides
  the draw (its chain) and the propagation (none).
- The categorical walk is written once, as `pick_realization(marginals, base, u)`. It keeps the
  legacy arithmetic: a `double` running sum of `long double` entries, `>=`, and
  `event_realizations` order. The uniform is drawn before the row is looked up, as each override
  did.
- `Deletion`'s `mutable` scratch strings are gone; the palindrome is a local.
- Tests: a case on the split itself (the draw writes nothing; the construction consumes no RNG,
  for each event type), and D4 as two `[!shouldfail]` cases through `pick_realization()`. Every
  G0 case passes unchanged.

### 5.4 G2, delivered

`GenerationState`, header-only, with 5 cases in tst/igor/Model/Legacy/test_generation_state.cpp:
- not built: reading or modifying throws;
- built and empty is distinct from not built;
- a modifier edits in place;
- assembly in the ordering: VDJ built out of order, a segment nobody built, a registered segment
  outside the ordering, tandem D;
- a registry that is not frozen is refused.

Three mutations, all caught: assembly over every id instead of the ordering; every segment built
up front; `create` appending.

### 5.5 G3, delivered

- `GenModel::generate_unique_sequence()` builds each sequence in a `GenerationState` sized from
  the model's registry, and assembles it in `registry.ordering()`. A model with no ordering is
  refused. Both programmatic callers of `generate_sequences` (igor-demo, the legacy demo) load
  from a file, so `finalize()` has set one.
- The construction switches stayed. They address the state through
  `Rec_Event::legacy_segment()`, which looks a legacy name up in the generator's registry rather
  than casting the enum, so a v2 file is addressed correctly whatever its ordering.
- The test fixture now holds a real `GenerationState`; no case changed.
- Gates: unit 480/480, integration 10/10, regression bitwise on all five tracks, convergence 2/2;
  the seeded TRA and TRB pins unchanged.

### 5.6 G4, delivered

Three commits, each bitwise on all five tracks, with the seeded TRA and TRB pins unchanged:

| Commit | What went | Tandem cases turned green |
|---|---|---|
| **G4a** | `Gene_choice`'s three-arm `switch (event_class)`; `Insertion`'s `insertion_seq_type_str_to_enum()`, which had no caller left. Each creates the segment at its own `seq_type_id` | gene choice, insertion |
| **G4b** | `Deletion`'s four-arm `switch (target_seq_type)`: one body on the segment at its `seq_type_id` and the end its side names, with D1 carried by position in the ordering (§5) | deletion |
| **G4c** | `Dinucl_markov`'s enum handles: it reads and writes at `get_junction()`'s ids. `DinuclTraversalSpec::target_seq`, `anchor_seq` and `legacy_enums_valid`, `kLegacySeqTypeCount` and `Rec_Event::legacy_segment()` are deleted | the D1-D2 chain, end to end |

- **One departure, in G4b.** A deletion of an internal segment that names neither side used to
  write nothing; it now throws, with the message inference's `initialize_event()` gives the same
  model.
- **Two test adaptations.** G4b's trim helper places its deletion in a VDJ ordering, as
  `finalize()` places every event, because anchoring is now read from the ordering. G4c deletes
  `test_dinucl_markov_iterate.cpp`'s assertion on `legacy_enums_valid`, together with the field.
- **The D1 cases are a live pin.** A mutation that drops the carried rule makes them pass, which
  Catch reports as two failures.

The generation path now names no gene. `draw_realization`, `pick_realization`, the four
`construct_realization` overrides, `draw_chain` and `generate_unique_sequence` address segments
only by `SeqTypeId`, and the milestone-1 layout generates end to end at unit level. Generating
from a tandem model *file* still waits for T1, which makes such a file loadable.

Left in this plan: the R phase (§5). `Deletion::target_seq_type`, `Insertion::ins_seq_type` and
`Dinucl_markov::ins_seq_type` are no longer read by generation. Their readers left are the
constructors, `copy()` and the generated names, which are T1's and T2's ground, as the parent
plan anticipated ("T3 and T4 shrink T1").

### 5.7 R, D1 delivered

**Decision (Quentin, Oct 5 2026):** a deletion trims the end its event names, with no special case
for the segments at the ends of the ordering. No model to date has a V 5′ or a J 3′ deletion, so
carrying the legacy V and J behaviour kept, by hand, a rule no output could show. `Deletion`'s
body reads `event_side` alone, and refuses a deletion that names neither end. The two D1 cases
lost their tags.

### 5.8 R, D2 and D3 delivered; D4, D5 and D6 decided

**Decisions (Quentin, Oct 5 2026):**
- **D2: refuse what inference drops.** A deletion past its segment's end, or a palindrome longer
  than the template it mirrors, now throws `std::out_of_range` on either side, naming the event
  and the lengths. Before, a 5′ trim clamped and a 3′ trim threw from `erase()` or `substr()`.
  The throw is what is coherent with the event's inference side. The fuller answer is a
  neighbour-type question, left open: whether a deletion may run on into the next junction
  (deleting insertions there), and whether a palindrome may use nucleotides past the gene it is
  cut from.
- **D3: solved by D2.** With D1 fixed, D3 is the keep-one case of the same refusal. A segment on
  an end of the sequence (no neighbour on one side of the ordering, as `finalize()` sets it)
  keeps at least one nucleotide, whichever side is trimmed. This is how inference resolves
  `keep_one_nucleotide_`. An internal segment may still be deleted away, leaving it built and
  empty.
- **D4: left to the sampling handlers.** Fixing the walk changes which realization a given
  uniform draws, a small behaviour change, and the walk itself goes away when the sampling
  handlers take over the draw. `pick_realization()` and `Dinucl_markov::draw_chain()`, which has
  the same `>=`, carry a `FIXME`, and the two `[!shouldfail]` cases stay.
- **D5: kept.** The behaviour is sound. A load-time rule, the parent plan's Phase C, may refuse
  such a model.
- **D6: waits on B10.** Generation and inference behave the same, and B10 decides for both.

**Tests.**
- The trim and palindrome cases assert the refusal and its message, on both sides, for internal
  and anchored segments.
- The cross-check against inference now runs every value from a palindrome one longer than the
  template to a deletion one past it, on the four shipped arms. Every value inference keeps comes
  out of generation identical, and every value it drops is refused.
- Three mutations, all caught: drop the keep-one clause, refuse an internal segment deleted to
  empty, and drop the palindrome bound.

**Does the refusal bite?** No. The legacy generator (`generate.fast false`) produced 2M sequences,
seed 20261005, from each of nine models:
- human IGH, IGK, IGL, TRA and TRB;
- the two supplementary TRB models, `naive_1` and `naive_4`;
- mouse TRB;
- the `generate` track's inferred model.

None threw. Every D model conditions `d_3_del` on `d_5_del`, so an inferred model puts no mass on
a pair of D trims longer than its gene.

## 6. Left for the SamplingEngine connection

Recorded here, not done in this plan:
- the uid ↔ `Rec_Event` mapping (T2, `LegacyBridge`);
- a construction order by role (Creates before Modifies, then the sequence-valued readers),
  instead of the queue;
- `Dinucl_markov` drawn after construction, through
  `sampleSequence(first_state = anchor nucleotide, n_steps = length)`;
- the RNG order changes there, so the goldens move once (synthesis step 1);
- the `MarkovSamplingHandler` row-index defect (REVIEW §11).

**`FastGenerator` gets its speed from the draw**:
- CDFs precomputed per condition, sampled by binary search or the alias method (FastSampling.h);
- OpenMP chunks;
- batched I/O.

Its construction is the same kind of code as the legacy one, over the same
`unordered_map<Seq_type, string>`. Once the draw is its own function, an indexed draw plus this
construction covers what `FastGenerator` does. Whether it is then retired or kept as an approximate
engine (synthesis §4) is outside this plan.
