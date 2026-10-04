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
| `propagate_realization(indices, index_map, offset_map)` | non-virtual | The loop now copied in three places. `Dinucl_markov` has none, as today |
| `construct_realization(indices, segments)` | pure virtual | No RNG, marginals or `index_map`. The later `apply` |

The indices are an index list shaped like `SampledEvent::indices`. Construction finds a realization
by index through a lookup on `Rec_Event`, rebuilt when realizations change.

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
| D1 | V 5′ / J 3′ deletion | trims the side the model names | the V arm always trims 3′, the J arm always 5′ | **confirmed defect**: two `[!shouldfail]` cases; fixed in R |
| D2 | deletion past the segment's end; palindrome longer than its template | drops the realization | a 5′ trim clamps; a 3′ trim throws `std::out_of_range` | observed, intent undecided |
| D3 | V or J deleted to empty | keeps one nucleotide: drops the realization | an empty segment | observed, intent undecided |
| D4 | CDF walk boundaries | — | `>=` picks a zero-mass realization at u = 0. When the row's mass stays below u, nothing is written, `()` is recorded, and the next reader throws. G0's fixture hit the second case for real: a misordered marginals file gave an insertion an all-zero row, and its chain then threw on the missing segment | to pin once the walk is a function of u (G1) |
| D5 | an insertion with no Markov chain | — | its `'I'` placeholders reach the output | observed, intent undecided |
| D6 | an empty chain anchor (a fully deleted D1 anchoring D1D2) | throws (§7.12, R2) | throws `out_of_range` from `substr` | consistent; B10 decides for both |

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
| **G2** | `Generated_seq_p_map` = `DynamicSequenceMap<std::string*>`, and assembly in `registry.ordering()` next to `build_scenario_sequence(registry, …)`. No caller | **yes** |
| **G3** | The generator adopts the container. The remaining switches address it through the legacy ids pinned to `Seq_type`. Assembly by ordering, which gives the same string for VDJ and VJ, since no segment outside the ordering is ever written. An empty ordering throws | **yes** |
| **G4a** | `Gene_choice` writes at `seq_type_id`; then `Insertion`. `insertion_seq_type_str_to_enum` goes if nothing else calls it | **yes** |
| **G4b** | `Deletion`: one body on `event_side`. D1 is carried generically: a segment with no left neighbour is trimmed at its 3′ end and one with no right neighbour at its 5′ end, whatever the side. This is B11a's anchored-segment boolean. D2's clamp/throw asymmetry falls out of the shared body unchanged | **yes** |
| **G4c** | `Dinucl_markov` reads at `get_junction()`'s ids. `legacy_enums_valid`, the spec's `Seq_type` fields and `kLegacySeqTypeCount` retire. The tandem end-to-end case turns green | **yes** |
| **R** | D1 (trim by side); D4 (walk by `u < cdf`, with the tail mapped to the last realization with mass); then D2, D3 and D5 once their intent is decided. One commit each, each removing its tag | each names what it moves; bitwise expected on the corpus |

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
