# Segment decomposition: review, and the minimal Event

**Status:** response to `SEGMENT_DECOMPOSITION_PROPOSAL.md`, for discussion.
**Date:** 2026-10-01.
**Source:** Q. Marcou's remarks on the proposal, plus an analysis of the code on `feature/tandemD`.
**Related:**
- `SEGMENT_DECOMPOSITION_PROPOSAL.md`, cited as "the proposal", with its section numbers.
- `MODEL_TOPOLOGY_REFACTORING.md`, a copy of tk's `src/igor/Model/README.md`.
- `TK_REFACTORING_MERGE_ANALYSIS.md`, cited as MERGE_ANALYSIS.
- `REC_EVENT_CAPABILITY_REFACTORING_PLAN.md`.
- `ITERATE_GENERIC_REWRITE_PLAN.md`, cited as ITERATE_GENERIC.
- `PROBA_BOUND_MACHINERY.md`.
- `CONTEXT_BASED_ITERATE.md`.
- `ITERATE_TEST_GUIDE.md`.
- `ARCHITECTURE_SYNTHESIS.md`, the summary that merges this review with the proposal and the Model module design.

**Code references:**
- Unmarked references are on `feature/tandemD` under `src/igor/Core/`.
- References marked **TK** are on the `origin-feature-tk_refactoring` worktree under `src/igor/Model/`.

**Literature:** the biology references were checked against PubMed on 2026-10-01 and are listed with their DOIs in §16. The algorithm references are outside PubMed and still need checking.

## 0. Summary

The proposal is a good description of one inference engine. It is not a good definition of the model.

- **What is kept:**
  - the registry/layout split;
  - the conditioning graph as model data;
  - removing the generated name;
  - the indexed sampler;
  - the exhaustive-run validation.
- **What is changed:**
  - Events keep their explicit attributes: targets, side, and library (today's `gene_class`).
  - The segment kinds and cardinality rules become admissibility conditions of one engine.
  - The layout becomes a DAG.
  - Handlers are declared, not inferred from tensor rank.
  - Error models become ordinary nodes with their own handlers.
  - Per-segment tables are one engine (`ChainDPEngine`) among several, next to the legacy DFS, which is kept.

**Minimal Event.** The event holds structure and semantics: identity, domain, targets, footprint, and the effect of a realization. Handlers hold the numbers. Engines own order, pruning, indexing, layers, caching, accumulation and threads.

Computation is placed by **tier** (§13), according to what it depends on:
- **0a** structure;
- **0b** parameters;
- **1a** read;
- **1b** read × parameters;
- **2** scenario.

Only tier 1a can be cached across EM iterations, and that cache is the proposal's per-segment table.

**Model file** (§14, §15):
- Events are keyed by nickname. `type` selects the constructor, and the order of entries is not content.
- A realization is identified by its label. Its index exists only where parameter values exist, and parameters are matched by label.
- An authoring form and a canonical form of one schema keep the file editable by hand and reproducible bit for bit.
- Parameters live in their own section of one sectioned document. Each section is inline or referenced and fingerprinted, and keyed by nickname, with labelled axes and round-trip precision.
- Any subset of sections and nodes can be exported, then overlaid onto another model.

Most of what `Rec_Event` holds today can move to an engine. That includes the whole probability-bound machinery, `iterate_wrap_up`, the layers, and the per-scenario state. The five existing Contexts are kept: what one realization produces is framed as *one layer of each Context it writes*. That choice lets the same event code serve both a zero-copy DFS and an orchestrator-driven DFS.

---

## 1. Agreements

- **The registry and the layout are separate objects.**
  - `SeqTypeRegistry` becomes a pure name ↔ `SeqTypeId` translator (MERGE_ANALYSIS §0.1).
  - The ordering, today built into the registry, becomes the layout (§4).
- **The conditioning graph is model data** (proposal §2.4): parents, a topological order, and an acyclicity check at load time.
- **The generated name goes**, replaced by a UID assigned by the model.
  - TK already has `uid()` / `setUid()`, and MERGE_ANALYSIS §5 proposes `uid == event_identifier`.
  - The nickname stays as the human handle in files.
  - The generated name is the key of `Model_Parms::edges`, `get_index_map`, `processed_events` and `offset_map` (Rec_Event.cpp:104-108, Model_Parms.cpp:148-191, Model_marginals.cpp:868).
  - It **already collides on tandem D**: `d1_choice` and `d2_choice` both become `GeneChoice_D_gene_Undefined_side_prio6_size<N>`.
  - `operator==` (Rec_Event.cpp:84-102) ignores `seq_type`.
- **An indexed CDF sampler for generation** (proposal §4). It is independent of everything else and can land first, as one deliberate commit that moves the goldens.
- **Validation check 1** (proposal §7) runs the new computation against the exhaustive legacy run at threshold 0. It is the acceptance test for *any* new engine.
- **Handlers are declared per node in the model file** (`"parameterization": {...}`), with a default per event type. Rank becomes a consistency check, not a selector. Details in §6.
- **Open points kept as stated:** parameter tying (D1/D2) and per-scenario outputs (proposal §9).

## 2. Event attributes: stated on the event, or implied by where it is filed

### 2.1 What the proposal does

"An event no longer carries a gene class, a side or a generated name" (proposal §2.1). Yet the JSON of proposal §3 still carries side and target, only as positions in a schema:

- `"modify": {"3'": ...}` gives the side;
- membership in a segment gives the target;
- `gene_class` disappears entirely.

The information has not gone away: it has become implicit. Part of it is lost outright.

**`gene_class` is not `seq_type`, and tandem D shows why.**
- `gene_class` names the **template library and its alignment set**: `query.gene_alignments.at(this->event_class)` (Genechoice.cpp:195).
- `seq_type` names the **slot in the layout**.
- D1 and D2 are two slots that share one library. In the proposal, segment "D1" would have to learn, implicitly, that it uses the D alignments.

**Side carries most of Deletion's semantics:**
- `trims_three_prime_`, the sign of `OffsetDelta`, the palindrome side, the direction in which the mismatch list is trimmed, and the junction slot it fills (Deletion.cpp:331-361, 477-518, 863-866, 994).
- For Dinucl_markov, side gives the direction of the chain and the neighbour that seeds it (Dinuclmarkov.cpp:94-101).

The current refactoring goes the other way. Model facts are *stated* on the event, and the subclass *answers* capability queries (MERGE_ANALYSIS §0.2). Examples are the A0 queries `get_offset_role`, `get_seq_construction_role`, `get_offset_delta_bounds` and `get_length_contribution` (Rec_Event.h:237-246). The proposal derives an event's meaning from where it is filed instead.

### 2.2 Pros and cons

| | **Stated** (explicit attributes plus answered capabilities, the current direction) | **Filed** (meaning from schema position, the proposal) |
|---|---|---|
| **Pros** | Each event describes itself. One event can have several targets, and one end can have several modifiers. A new event type needs no schema change. An event can be tested on its own (today's `make_deletion(...)` fixtures). The conditioning graph and the handlers refer to every event the same way. | Well-formed by construction: the schema *is* Phase C. One source of truth. Compact files for standard models. Engines can rely on the structure without checking it. |
| **Cons** | Redundancy allows inconsistency (an event targets D1, but nothing creates D1), so Phase C validation is required. The per-segment view has to be derived. It is easier to write a model that makes no sense. **The modeller implicitly has to know event internals** to write a consistent model: that a Deletion moves an offset, that an Insertion needs a Dinucl, which ends an event reads. | A closed world: each new arrangement changes the schema, the reader and every engine. Semantics become implicit (side from a key, library from a segment name). An event cannot be described or tested outside its slot. Relations across segments cannot be expressed. The modeller has to know the schema instead, though it is documented in one place. |

### 2.3 Middle ground

- **Canonical form:** the attributes stated on the event.
- **Derived index:** the model builds a `SegmentIndex`, listing for each segment and each end its creators and its modifiers. Engines and validation read it.
- **Shorthand:** the file format **also accepts the proposal's `segments` block, as a shorthand** that the reader expands into stated attributes. Standard models stay compact, and the general form stays available.
- **Reducing the "internals" cost.** What a modeller needs to know is exactly the answered capabilities (§7.1, ingredient 4): roles, deltas, and what an event reads.
  - **Each event type publishes its contract.** Because the capabilities are declared data, they can be printed. An `igor describe-event Deletion` lists what the type creates, modifies, reads and requires, so the modeller learns a contract, not code.
  - **Phase C diagnostics speak the modeller's language:** "segment D2 has no creator", "junction D1D2 has a length event but no content event". They never refer to how `iterate()` behaves.
  - **The shorthand covers the standard patterns** (a gene with trimmed ends; a junction with a length and a content). A modeller only meets the contracts when leaving those patterns.

## 3. Segment kinds and cardinality rules: what they exclude

Two rules in the proposal describe the models that exist today:
- "A gene segment has one Create event and at most one Modify event per end."
- "Junction content depends on its two neighbours and on nothing else."

Each case below is excluded by them, and each is real biology:

- **V(H) replacement: two V-templated contents in one read.**
  - An upstream VH recombines into an already rearranged VHDJH. It uses a cryptic RSS heptamer (TACTGTG) near the 3′ end of the rearranged VH.
  - It leaves a **footprint** inside CDR3: a few nucleotides of the old VH 3′ end, often with N additions.
  - The embedded heptamer is conserved in over 70 % of known VH genes (Kleinfield et al. 1986). Replacement footprints are found in normal human IgH sequences (Zhang et al. 2003).
  - Refs: Reth et al. 1986, *Nature* 322:840; Kleinfield et al. 1986, *Nature* 322:843; Zhang et al. 2003, *Immunity* 19:21; Zhang, Burrows & Cooper 2004, *Immunol Rev* 197:231.
  - In IGoR terms, this is a second V GeneChoice that creates a footprint segment, or a path through an extra segment (§4).
- **Receptor editing (Igκ) and successive TRA rearrangements.**
  - Secondary Vκ–Jκ joins use an upstream V and a downstream J, and delete the primary join. Refs: Tiegs, Russell & Nemazee 1993, *J Exp Med* 177:1009; Gay et al. 1993, *J Exp Med* 177:999; Radic et al. 1993, *J Exp Med* 177:1165; review in Nemazee 2017, *Nat Rev Immunol* 17:281.
  - TRA goes through several successive Vα–Jα rounds on one allele. Ref: Petrie et al. 1993, *J Exp Med* 178:615.
  - The read still shows one V and one J, so this is **not** two creators. It is a **mixture of histories** (primary or secondary), each with its own (V, J) statistics and possibly its own trimming statistics. That is a structural choice (§4) plus conditioning, not a cardinality.
- **D–D fusions and D inversions (IGH).**
  - D–D fusions occur as direct or as **inverted** rearrangements. Ref: Meek, Hasemann & Capra 1989, *J Exp Med* 170:39, "Inversions and fusions add to IgH somatic diversity".
  - Frequency of V(DD)J: about 1 in 800 circulating human B cells (Briney et al. 2012, *Immunology* 137:56). Tandem CDR3s are present in every immunosequencing dataset analysed by Safonova & Pevzner 2019, *Front Immunol* 10:987.
  - An inversion is a Modify that changes the **orientation of content**, not its extent.
- **TRD with up to three D genes** (VDDJ, VDDDJ). Ref: Lefranc & Lefranc, *The T Cell Receptor FactsBook*, 2001.
- **P nucleotides.**
  - They come from asymmetric hairpin opening, which happens before exonuclease trimming. Ref: Lafaille et al. 1989, *Cell* 59:859.
  - IGoR encodes them as negative deletions (Marcou, Mora & Walczak 2018, *Nat Commun* 9:561; Murugan et al. 2012, *PNAS* 109:16161). That is one modelling choice.
  - Hairpin opening and trimming as two events on one end is another choice, and the cardinality rule forbids it.
- **Context that crosses a segment boundary.** N-mer hypermutation windows span boundaries; proposal §5.5 already concedes this.

**Conclusion.** These rules become the **admissibility condition of an engine** (§5), not the definition of the model. Phase C validation stays at model level. The rules are then per-engine requirements, or warnings.

## 4. Layouts with several paths

### 4.1 Uses beyond tandem D

A layout with several paths is not one feature for one event. Uses:

1. **Number of D segments:** IGH V(DD)J; TRD with 0 to 3 D.
2. **D-less joins:** "D absent" kept distinct from "D fully trimmed" (IGH, TRB, TRD).
3. **Locus mixtures in one dataset:** unsorted IGK and IGL light chains; TRA/TRD with shared TRAV/DV genes; TRG/TRD.
4. **V replacement footprint:** `V → Vfp → VD → D → DJ → J`.
5. **D inversion** as an alternative segment `D_rc`, if it is not modelled as a content Modify.
6. **Optional technical or extra segments:** the constant region or isotype for reads that cover C; leader/UTR; UMI, primer or adapter.
7. **Inter-locus trans-rearrangements,** e.g. the TRG–TRD chimeras Vγ-(Dδ)-Jδ and Vδ-(Dδ)-Jγ found in normal lymphoid tissue (Tycko, Palmer & Sklar 1989, *Science* 245:1242).
8. **Primary versus secondary rearrangement histories** (§3).

### 4.2 Not tied to a single "path event"

The general form is a **layout DAG whose segments carry existence guards** on **structural events**.

- A guard is a predicate on realizations, for example `present_if: {n_d: [2]}`.
- Several independent structural events can each control a branch point: locus × number of D × D orientation.
- A path is the set of segments whose guards hold.
- **Validity rule:** every joint assignment of the structural events yields exactly one source → sink chain.
- Structural events are ordinary events. They have a handler, take part in EM, and can be conditioned, e.g. P(two D | V).
- Events attached to a guarded segment exist only when the guard holds. This *existence dependency* is distinct from the conditioning graph.
- **Adjacency becomes a per-scenario query** (`left_neighbour(seg)` on the scenario state). It is no longer stated once at `finalize()` (Rec_Event.h:260).
- **Variant:** a segment can be made optional by an "absent" realization of its creator. That is more compact. But rewiring the junctions (VD1 + D1D2 + D2J collapsing to a VJ-like junction) still needs the alternative junction segments in the DAG, so guards stay the general form.

### 4.3 A dedicated data structure, or hard-coded in an experimental engine

| | **A lightweight `SegmentLayout` DAG in the model** | **Hard-coded in an experimental engine** |
|---|---|---|
| **Pros** | The model stays declarative and serialisable. Every engine and sampler shares one definition, so inference and generation agree by construction. Structural events go through handlers and EM, and can be conditioned. The file format is decided once. The structure is small: an adjacency list, guards, `paths()`, and acyclicity and validity checks. | The fastest route to a working tandem D. We learn what is actually needed (sharing junctions, cost) before committing to a format. No risk to the legacy code. |
| **Cons** | A design commitment made before the experiments, with a risk of over-generalising: the mirror image of the v1 error. Every engine must accept a general DAG or explicitly reject it. | The engine carries part of the model's meaning, so the model file no longer describes itself. The path probability is a parameter outside handlers and EM, and cannot be conditioned. The sampler has to repeat the same hard-coding, which gives two sources of truth. The temporary form tends to become the format (the v1 history). |

**Middle ground (recommended):**
1. **Now:** the model holds the small DAG, its guards and its validity check, and the file format carries it. This costs little. Structural choices are ordinary events.
2. **During the experiment:** the experimental engine **reads** that DAG, but declares a narrow admissibility condition. For example, "only the two-path tandem-D layout", enumerated as a mixture of linear chains. That is where the hard-coding lives.
3. **The deciding test:** a sampler and two inference engines must agree on one model while sharing nothing but the model.

## 5. Engines

### 5.1 The family

The per-segment table and forward-backward pass of proposal §5 make a good **`ChainDPEngine`**. It is admissible for models where every dependency between segments goes through the boundary interface, and where every per-scenario output requested can be written as a semiring.

It should not replace the legacy enumeration:
- On small combinatorics (TRA, IGK), the legacy DFS computes exact probabilities along a topological order very efficiently.
- It is the only engine today that supports arbitrary per-scenario `Counter`s (`count_scenario` at every leaf, Rec_Event.cpp:268-277).
- Its goldens are bitwise. Keeping it intact keeps the regression gate.

| Inference engine | Use | Status |
|---|---|---|
| `LegacyDfsEngine` | TRA/IGK, any counter, bitwise goldens | Today's `iterate()` and bound fold |
| `ChainDPEngine` | IGH, TRB, tandem D | Proposal §5 |
| `DpBoundedDfsEngine` (semi-DP, §5.3) | Exact likelihood together with arbitrary counters | New |
| Posterior sampling / stochastic EM | Very large models | FFBS on the DP messages (§5.3) |
| Viterbi / k-best | Best scenarios | A flag of the DFS today; (max,×) on the DP |

| Sampling engine | Use |
|---|---|
| Exact topological sampler | The reference, with an indexed CDF |
| Two-loop sampler (proposal §4) | Admissible only when no scalar event depends on constructed content |
| `FastGenerator`-like approximate sampler | Speed |

TK `InferenceEngine` has no E-step (InferenceEngine.tpp:144-151): `run(eStep)` only brackets a callback with `resetAccumulators()` and `updateParameters()`. Its test hands the E-step to the legacy `GenModel`. The E-step is precisely the part that should vary from one engine to the next.

### 5.2 The common workflow, read from `GenModel::infer_model`

`GenModel::infer_model` (GenModel.cpp:87-680) is today's only E-step implementation. Each step is marked **G** (generic, shared by every engine) or **D** (specific to the DFS).

| # | Step (lines) | G/D | Destination |
|---|---|---|---|
| 0.1 | Option checks; Viterbi forces the ratio to 1 (95-110) | G + D | The driver validates; the ratio is a DFS pruning parameter |
| 0.2 | Model plan: `get_model_queue`, `get_index_map`, `get_inverse_offset_map` (115-118) | D (order), G (parameter layout) | Order goes to the engine; parameter layout goes to the handlers |
| 0.3 | Logs, initial dumps, lists of fixed and inferred events (112-177) | G | Driver / output layer |
| 0.4 | Counters `initialize` once, for the file streams (184-197) | G | Driver (counter lifecycle) |
| 1.1 | Fresh accumulators: `new_marginals`, `error_rate_copy` plus `initialize(events_map)` (211-215) | G, except the raw-pointer `initialize` | `AccumulatorSet::empty_copy()` |
| 1.2 | `fast_iter`: best alignments only on iteration 0 (224-238) | G | A driver **data policy**, independent of the engine |
| 1.3 | Deterministic chunk plan under a memory budget (240-270, #68) | G | Driver reduction plan |
| 1.4 | Per-thread copies of `Model_Parms`, marginals, error rate and counters, filtered by `last_iter_only` (292-306) | G (counters), D (the copies exist only because events hold state) | The driver gives each thread a const model and an accumulator; the engine builds a workspace |
| 1.5 | Scenario state: SafetyMatrix, sequence, mismatch and offset maps, downstream map (313-327) | D | DFS workspace (its Scenario and Exploration contexts) |
| 1.6 | Index map, event sizes, crude bound (329-348) | D | DFS engine / handlers |
| 1.7 | `initialize_event` in queue order (layers, partners), Viterbi flag (359-369) | D | DFS engine plan |
| 1.8 | `next_event_ptr_arr` (377-391) | D | DFS engine dispatch |
| 1.9 | Error-rate `initialize(events_map)` (raw pointers), counters `initialize` (394-409) | G, made D by the raw pointers | Driver, once the scenario record exists (§7.2) |
| 1.10 | Bound fold on one thread, then adopt (418-452) | D | `prepare_iteration()` of the DFS engine (per iteration, shared, read-only) |
| 1.11 | `update_event_internal_probas` (459-464) | D/handler | Handler materialisation in `prepare_iteration` |
| 2.1 | Per read: fresh accumulator, thresholds, `nt2int`, contexts (488-544) | G (accumulator, encoding), D (thresholds) | Driver: encoding and per-read accumulator. Engine: its contexts |
| 2.2 | **E-step kernel**: `first_event->iterate(...)` (547) | **D** | `engine.process_read(...)`, the only call that differs between engines |
| 2.3 | Logging of the exception context (551-560) | G | Driver |
| 2.4 | Normalisation by the read likelihood, `norm_weights_by_seq_likelihood`, held by the error rate (563) | G | Accumulator. The likelihood becomes an output of a `ReadSummary`, not state of the error rate |
| 2.5 | Per-read log line (566-577) | G | Driver, from `ReadSummary {likelihood, mean_errors, n_scenarios, best}` |
| 2.6 | `count_sequence`, `dump_sequence_data` (578-586) | G | Driver (per-read counter hooks) |
| 2.7 | Read rejection on mean error count (588-597) | G | Driver read policy |
| 2.8 | Fold into the chunk slot; drain the error rate at chunk end (593, 599-605) | G | Driver reduction |
| 3.1 | Merge chunks in order; counters merged in thread order (the FIXME at 632-649) | G | Driver. Counters get chunk slots once accumulators are uniform |
| 3.2 | M-step: `error_rate->update()`, `normalize`, `copy_fixed_events_marginals` (660-664) | G | `maximize()` on every node, error model included; `fixed` becomes a handler flag |
| 3.3 | Likelihood file, per-iteration dumps, counter summaries (651-670) | G | Driver / output layer |

**Of about 600 lines, the DFS-specific part is small:** step 0.2 (order), steps 1.5 to 1.8, steps 1.10 and 1.11, and step 2.2. The rest is a generic **EM driver**:

```
driver.run(model, reads, engine, counters, policies):
  engine.check_admissible(model, counters)
  for it in iterations:
    acc  = AccumulatorSet::empty_copy(model)                 // handlers + error model + counters
    plan = engine.prepare_iteration(model)                   // DFS: bound fold; DP: refresh handler rows
    parallel over fixed chunks (#68):
      ws = engine.make_workspace(plan)                       // DFS: layered contexts; DP: message buffers
      for read in chunk:
        r_acc   = acc.empty_copy()
        summary = engine.process_read(policy(read, it), ws, r_acc, counters.scenario_sink())
        r_acc.normalize(summary.likelihood); counters.count_sequence(summary, r_acc)
        if policy.accept(summary): chunk_acc += r_acc
    acc = ordered_merge(chunk_accs); model.maximize(acc); output(it, acc, summaries)
```

- **Every engine implements four calls:** `check_admissible`, `prepare_iteration`, `make_workspace` and `process_read`.
- **An optional fifth call, `prepare_read(read) -> ReadCache`.** It produces a per-read cache that survives from one iteration to the next: the `prepare(read)` hook of §7.1, which is also the proposal's per-segment table, reused by every EM iteration. The driver holds the cache, under a memory policy.
- **Generation mirrors this with a `SamplingDriver`:**
  - per-chunk RNG streams seeded by `(seed, chunk)`, for determinism;
  - an error-application step;
  - the output writers;
  - all around `engine.sample(rng) -> ScenarioRecord`, followed by the events' `apply`.

### 5.3 Dynamic programming and arbitrary Counters

Below, P(s) is the probability of scenario s and L = Σ_s P(s). Every counter today falls into one of four classes:

| Statistic class | Counters today | DP method |
|---|---|---|
| Posterior expectation of a function of one segment state | `Coverage_err_counter` (allele, alignment and deletion of one segment); the error-rate sufficient statistics; the marginals; `Errors_counter` (a) Σ P·n_gen and (b) Σ P·n_mis | (+,×) forward-backward posterior marginals, or the expectation semiring (Eisner 2002; Li & Eisner 2009). Exact. |
| Non-linear function of additive quantities | `Errors_counter` (c) Σ P·n_mis/n_gen; histograms of error counts | State augmentation (carry the value's range). Cost × range. |
| Top-k | `Best_scenarios_counter`, with ties and `keep_ties` | k-best (max,×) with backtracking (Huang & Chiang 2005). The tie order and the unbounded tie group need care. |
| Arbitrary function of the whole scenario | `Pgen_counter` (keyed by the full constructed sequence, then a log); any user counter; per-scenario dumps | Not decomposable: the key grows with the mismatches. Needs enumeration or sampling. |

**Three "semi-DP" designs keep every counter:**

1. **DFS bounded by DP.** Compute the per-read tables and the DP backward messages first. Then enumerate depth-first along the chain, using the messages as the downstream bound.
   - **Max-product messages are *exact* upper bounds.** No subtree that cannot produce a leaf above the threshold is ever entered.
   - **Sum-product messages give each subtree's exact mass.** The pruned mass is therefore known, so the error is quantified, and L is exact.
   - **The counter callback is unchanged** and sees every scenario visited.
   - **This replaces the whole bound fold** (PROBA_BOUND_MACHINERY) with DP messages.
   - **Cost:** the DP, plus an enumeration whose cost is proportional to its output.
2. **Posterior sampling by FFBS (forward-filtering, backward-sampling).**
   - **Algorithm.** Take the chain of boundary states `x_0 … x_K` (proposal §5.3).
     - The forward pass `α_k(x) = Σ_{x'} α_{k-1}(x') ψ_k(x', x)` is exactly the (+,×) pass the likelihood already needs, with `Z = Σ_x α_K(x) = L`.
     - Then sample backwards: `x_K ~ α_K(·)/Z`, and for `k = K-1 … 0`, `x_k ~ α_k(x_k) ψ_{k+1}(x_k, x_{k+1})`. Each draw is a categorical over one table row, so one sample costs `O(K × |states|)`, with no recomputation of α.
     - Then fill each segment's interior given the sampled boundary states. For example, D's pair `(del5, del3)` given D's two offsets, drawn in proportion to its local factor.
     - On the outer loop over gene tuples: draw the tuple with probability ∝ `P(tuple) · Z_tuple`, then run FFBS inside it.
     - The result is an **exact** draw from `P(s | read)`. The same scheme extends to a tree of segments: draw the root, then its children.
     - Refs: Carter & Kohn 1994, *Biometrika* 81:541; Frühwirth-Schnatter 1994, *J Time Ser Anal* 15:183; Scott 2002, *JASA* 97:337.
   - **Uses:**
     - (a) **Arbitrary counters by Monte Carlo.** `E[f(s) | read] ≈ (1/M) Σ f(s_m)` is unbiased, with an error of `O(1/√M)` however many scenarios exist. For `Pgen_counter`, `E[log Pgen(σ)]` is estimated directly; the per-σ table only approximately.
     - (b) **Stochastic, or Monte Carlo, EM.** Accumulate counts from M samples per read. This is cheaper when exact expectations would need augmented state. SAEM handles convergence (Delyon, Lavielle & Moulines 1999, *Ann Stat* 27:94). This is the "stochastic EM" engine.
     - (c) **Posterior predictive checks.**
     - (d) **Importance sampling for models the DP does not admit.**
       - Run FFBS on an admissible approximation `q`, for example one that ignores an N-mer error window across a boundary.
       - Reweight each sample by `p(s, read) / q(s, read)`, evaluated through the generic `evaluate` path (§10).
       - The approximation is then corrected exactly in expectation, so the DP stays useful beyond the models it admits.
   - **Limits:**
     - The samples are exact only if the DP is exact (no pruning of the tables).
     - In (d), the variance grows when `q` is far from `p`.
     - Determinism needs one RNG stream per read, derived from `(seed, read index)`.
3. **Counters declare their class:** `posterior_marginals(events | segments)`, `additive(local f_e)`, `k_best(k)`, or `scenarios(threshold | samples)`. The engine serves each class with the cheapest exact method, and uses designs 1 or 2 only for the last class. The declared classes are also an input to each engine's admissibility check.

**A pure DP with no enumeration loses:**
- per-scenario callbacks;
- the exact Pgen of each sequence;
- best scenarios ordered by the legacy tie rules;
- every counter outside the first three rows of the table above.

Today's counters only see scenarios above the threshold, so exact DP counters will differ from them. That is consistent with validation check 1 at threshold 0.

**Prerequisite for every engine:** counters stop reading event members through `friend` raw pointers, and read a scenario record instead (§7.2).

## 6. Parameters: handlers, accumulators and error models

**The error models are parametric handlers in disguise.**

- **`Hypermutation_full_Nmer_errorrate`** is a categorical over {mutated, not mutated}, conditioned on the N-mer context.
  - Its tensor is `[4^N, 2]`.
  - Its accumulators `Nmer_N_SHM` and `Nmer_N_bg` are exactly the count tensor.
  - Its M-step is a regularised ratio: `N_SHM/N_bg` above a count threshold, the median otherwise (HypermutationfullNmererrorrate.cpp:1029-1077).
- **`Hypermutation_global_errorrate`** is a **logistic regression** on one-hot position × nucleotide features: `logit p_k = log mu + Σ_i e_i(nt_i)`.
  - It has the **same sufficient statistics** as the full N-mer model.
  - Its M-step is Newton plus an Armijo line search on the binomial log-likelihood (Hypermutationglobalerrorrate.cpp:1148-1425, 1697-1710, 1830-1839).
- **`Single_error_rate`** is a Bernoulli.
  - Its M-step is the mean of per-sequence frequencies, not the pooled MLE (Singleerrorrate.cpp:311-318). It must be kept as is for the goldens to stay bitwise.

So the three error models are **one Probe event with three handlers**. The event is sequence-valued, with a context window. The handlers are `Categorical`, `Logistic` (features = the one-hot N-mer) and `Bernoulli`. The accumulator has the same shape in all three (context × domain); only the M-step differs. That is the case for breaking the distinction between how event parameters and error-model parameters are built.

**Covariates** (patient age, clonal abundance) are exogenous observed nodes in the conditioning graph. A `LinearModel` or `Logistic` handler, for instance P(insertion length | age, abundance), cannot use a dense tensor over continuous covariates as its accumulator. The accumulator has to be **owned by the handler**, e.g. Σ x xᵀ w, or per-read posterior rows stored with x.

**Handler contract.** This generalises TK's `InferenceHandler` and `SamplingHandler`:

```
distribution(context)        -> row over the event's domain   // context = parents' realizations + covariates
accumulate(context, posterior row | (realization, weight))
merge(const Accumulator&)    // ordered, #68
maximize()
sample(context, rng)
bound(partial context)       -> max over the unresolved context // what DFS pruning needs
```

- `bound` replaces `realization_bound()` (Rec_Event.cpp:447-460), and Dinucl_markov's use of `event_upper_bound_proba` in `span_proba_factor`.
- Accumulators are objects separate from the handlers, with `empty_copy()` and an ordered `+=` (MERGE_ANALYSIS §6.2).
- Tensor-based engines get dense per-read rows by calling `distribution()`.
- The flat-array indexing of the legacy engine (`Index_map`, `iterate_common`) stays valid only for tensor handlers. That limit belongs to the engine, not to the model.

**Handler selection.**
- TK selects the handler by `Event_type` (InferenceHandlerFactory.cpp:20-33).
- The proposal selects it by rank (proposal §2.5).
- This review: neither. The handler is declared per node in the model file, with a default per event type. Rank is a consistency check.

## 7. The minimal Event

**Principle:**
- The **event** is structure and semantics: what it is, what it touches, and what one of its realizations does.
- The **handler** holds the numbers.
- The **engine** owns order, pruning, indexing, layers, caching, accumulation and threads.

### 7.1 What stays on the event

| # | Ingredient | Today | Form kept |
|---|---|---|---|
| 1 | Identity | `name`, `nickname`, `type`, `event_index` | `uid` (assigned by the model), `nickname`, `type`. No generated name. |
| 2 | Domain | `event_realizations` (`unordered_map<string, Event_realization>`) | Realizations indexed `0..n-1`, each with a payload (gene name and sequence, integer, nucleotide). Plus an arity: scalar, or sequence-valued over the positions of a segment (Dinucl_markov, the Probe). The domain shape (`inherent_shape()`) is distinct from the parameter shape. |
| 3 | Targets (stated) | `seq_type` / `seq_type_id`, `event_side`, `event_class` | A list of `(seq_type, side?)` targets, which generalises "one seq_type, one side" and allows events on several segments. A `library` (today's `gene_class`) for templated events. Optionally `structural: true` (§4). Explicit and serialised. |
| 4 | Footprint (answered) | A0: `get_offset_role`, `get_seq_construction_role`, `get_offset_delta_bounds`, `get_length_contribution`; `affects_length_of/proba_of`, `length_delta` | Kept, and extended with the planned `SeqContextDependency`, i.e. what the event *reads* (REC_EVENT_CAPABILITY plan §A). Made independent of the registry: no R12 enum tables, no `kGeneSegments` (Genechoice.cpp:768), no `legacy_junction_of`. |
| 5 | Semantics | Buried in `iterate()` and `draw_random_realization()` | `candidates(contexts)`: the domain, in order, with a flag saying whether the bound is monotone along that order. `evaluate(r, const Query&, const Model&, const Scenario&, writer) -> Status`, which writes one `EventOutcome` (§9.2). `apply(r, generation state)` for generation. `prepare(read) -> ReadCache`, the tier-1a hook (§13; ITERATE_GENERIC §2.10). |
| 6 | Construction | Typed constructors per subclass; EventFactory default-constructs, then calls setters | **One generic constructor from a JSON description** (`EventClass(const json&)` / `from_json`), a matching `to_json`, and registration by type string in `EventFactory`. The `type` field dispatches; the nickname identifies (§14). |

The generic constructor (ingredient 6) replaces three paths:
- **The typed constructors** `Deletion(Seq_type, Seq_side, …)`, `Insertion(Seq_type, …)` and `Dinucl_markov(Seq_type)`. Their `str2SeqType` throws on `D1_gene_seq` (Model_Parms.cpp:904, 923, 953).
- **EventFactory's "default-construct, then setters"** path. It cannot reach subclass members, so `copy()` rebuilds `d_3_del` as a V deletion (MERGE_ANALYSIS §4.1).
- **`write2txt*`.**

The v1 and v2 readers become importers that emit JSON. Tests, `LegacyBridge` and the factory then share one construction path.

There is no JSON library in the C++ tree today: no `nlohmann`, and no `model_parms_to_json` was found. Adding one is a dependency decision.

### 7.2 What leaves `Rec_Event`

| Destination | Members and methods today |
|---|---|
| **Engine: dispatch and leaf** | `iterate_wrap_up` (Rec_Event.cpp:173-287: calls the next event through `next_event_ptr_arr`; at the leaf, the error rate, threshold, counters, and `add_to_marginals` over all of `events_map`). Also `event_index`; queue order, `Event_comparator`, `processed_events`; `initialize_event*` (tier-0 initialisation that depends on queue position, decision O3). |
| **Engine: layers** | Every `memory_layer_*`, `safety_layer`, `partner_offset_layer`; the layer claims; the claimed-layer snapshot (Rec_Event.cpp:322-324). |
| **Bound service** | `junction_bounds_`, `initialize_Len_proba_bound`, `adopt_Len_proba_bound`, `build_retained_decomposition`, `SpanFold`, `fold_step`, `realization_bound`, `creates_own_offsets`, `event_upper_bound_proba`, `set_crude_upper_bound_proba`, `scenario_upper_bound_proba`; per type `write_junction_bounds`, `surviving_core_bound`, `pending_` (`PendingModifierBounds`). Its only inputs are the footprint queries (§7.1, ingredient 4) and the handlers' `bound()`. |
| **Geometry and safety service** | `flank_checks_`; the overlap loops (Genechoice.cpp:217-240, nearly identical to Deletion.cpp:297-317); `SafetyMatrix`. |
| **Indexing (handler / engine)** | `memory_and_offsets`, `iterate_common`, `update_parent_tracking`, `Index_map`, `event_marginal_size`, `ind_normalize`, `update_event_internal_probas`. |
| **Scenario record** | `new_index`, `realization_indices`, `current_realizations_index_vec`, `current_realization_index`, `deletion_value`, `alignment_offset_p`; the scratch buffers. Replaced by `{uid → realization}` (TK `SampledEvent`) in `ScenarioContext`. **Prerequisite:** the `friend` access from `Coverage_err_counter` and both hypermutation models (Genechoice.h:54-56, Deletion.h:57-59; Coverageerrcounter.cpp:157-261; Hypermutationglobalerrorrate.cpp:1504-1630) moves first. |
| **Accumulator, sampler, serialiser, validation** | `add_to_marginals`, `viterbi_run`, `fixed`; `draw_random_realization`; `write2txt*`, `update_event_name`, `get_v2_name`, `get_legacy_name`; `require_dinucl_markov` and the side and registry checks in `initialize_event`. |
| **Delete (dead code)** | `compare_sequences`, `Gene_choice::add_realization(int)` and `Dinucl_markov::add_realization(int)` (declared, never defined); `initialized`; Insertion's `new_scenario_proba`, `previous_index`, `insertions_str`; Dinucl's `compute_nt_freq`, `previous_seq_size`, `total_nucl_count`; `ExplorationContext::pruning_mismatch_floor` (never read); `ModelContext::offset_map` inside `iterate` (unused); `len_min`/`len_max` once `EventUtils::get_insertion_len_max` uses `get_length_contribution`. |

The recent commits already go in this direction:
- R14 removed the crude-bound runtime.
- R15 turned the fold state into one object.
- 5b left no virtual call in the sweep.

The fold driver is no longer virtual, so moving it into an engine-side table keyed by uid is mechanical.

## 8. Anatomy of the refactored `iterate()`

### 8.1 The shared skeleton

All four bodies follow the same skeleton. ITERATE_GENERIC.md:20-22 summarises it as "pick a realization, move or create a segment, check it against its neighbours, bound the remainder, recurse". In the code:

```
once:  base_index ← index_map; base_p ← scenario_proba; read upstream state at layer-1
for r in candidates:                                   R  domain and order: event-specific
    admissible(r)?  else continue                      S  guards (specific) + overlap loop (generic)
    effect ← compute(r)                                K  event-specific
    p ← base_p · P(r|parents)  (+ index_map write)     P  generic (iterate_common)
    write seq / offsets / mismatches at own layers     W  generic (segment write triple G4)
    write bound slots (junction best_for, error bound) S+P+W  generic shape, specific inputs
    prune(p · Π slots)?  else continue                 P  generic
    scenario_proba ← p; iterate_wrap_up                W+R generic
```

Legend:
- **S**: safety and geometric admissibility.
- **K**: what the event knows about the scenario and computes from it.
- **P**: probability, bounds and pruning.
- **W**: writes to the scenario state.
- **R**: recursion and loop.

### 8.2 Where each event departs from it

- **Gene_choice** (Genechoice.cpp:159-557):
  - The candidates come from the read (its alignments), with two fallbacks that do not depend on the read: the retained decomposition and the sliding window.
  - It writes its sequence before the safety check.
  - It prunes a second time after scoring.
  - It publishes alignment pointers into its own members.
  - It copies `gene_seq` again for each alignment of each scenario (Genechoice.cpp:196-209), although the copy depends only on the read.
- **Deletion** (Deletion.cpp:214-467):
  - It inherits the safety verdict from layer-1, then re-asserts it.
  - Its two-stage prune `break`s out of the loop. That is valid only because the candidates are sorted by decreasing deletion (Deletion.cpp:422-456; `del_numb_compare`, 956-959).
  - It builds and guards palindromes.
- **Insertion** (Insertion.cpp:168-244):
  - It has a single candidate, fixed by the geometry.
  - Its admissibility test is "the realization exists".
  - It **throws** on an unreachable bound where the others `continue`.
  - It multiplies in place (`*=`) without a snapshot.
- **Dinucl_markov** (Dinuclmarkov.cpp:146-268):
  - It has one composite realization.
  - It copies the anchor segment.
  - It stores the junction pointer before filling the junction.
  - It sets its bound slot to 1.0 and multiplies in place.

### 8.3 What is specific to an event, and where each phase goes

| Phase | Specific to the event: stays in `evaluate` | Generic: moves to the orchestrator and its services |
|---|---|---|
| **S** safety | Semantic validity: keep one nucleotide, offset inside the read, palindrome guards, "the realization exists" | Overlap and flank checks between segments, which are a relation of the **layout**; the `SafetyMatrix` layers; inherited verdicts |
| **K** knowledge | Template clipping, mismatch trimming, palindrome strings, junction length, Dinucl's read window | — but **most of K is tier 1a** (§13: read × realization, independent of the scenario and of the parameters): `gene_seq` per alignment, palindromes, endogenous mismatch counts, D's sliding windows (ITERATE_GENERIC §2.10, L1466-1474). It moves to `prepare(read)`: the event computes it, the engine caches it. **This cache is exactly the proposal's per-segment table (proposal §5.2).** |
| **P** probability | The structure of the factor (scalar lookup or chain product); the *bound facts* (junction distances fixed, surviving-core counts); the candidate order and whether the bound is monotone along it | Marginal lookup (handler row), parent tracking, `compute_upper_bound` / `should_prune`, turning facts into factors (`best_for`, the error bound), Deletion's `break` (driven by the declared monotonicity) |
| **W** writes | *What* is written (the output of K) | *Where and how*: layers, claims, the segment write triple |
| **R** recursion | The candidate set | The loop driver, dispatch, the leaf |

## 9. Contexts: keep all five, finish the Context refactoring

The five contexts already split the state along the right lines:
- `QuerySequenceContext`: the read;
- `ModelContext`: the model;
- `ScenarioContext`: the current path;
- `ExplorationContext`: the exploration policy;
- `AccumulationContext`: the results.

The proposal here is **not** to replace them. It is to finish what CONTEXT_BASED_ITERATE.md already sets out (L536-541): mutable state goes into the contexts, and `ModelContext` is const. A large part of the scenario still lives in event members today (§7.2). The changes are about *access discipline* and *completion*, not new types.

### 9.1 What each context becomes

| Context | Shared by all engines? | Event access (`evaluate` / kernel) | Changes |
|---|---|---|---|
| `QuerySequenceContext` | yes | `const` | Gains the per-read caches produced by `prepare(read)` (tier 1a, §13). `gene_alignments` is selected by **library (`gene_class`)**, more evidence for §2. `journaled_query` is not read by events. |
| `ModelContext` | yes | `const` | `model_parameters` gives way to the handlers: the event receives **the row for its current parent context**, resolved by the engine, and never sees `index_map`. `events_map` and `model_queue` stay for counters and engines. `offset_map` is unused by events (Genechoice.cpp:662-670, Deletion.cpp:523-530). |
| `ScenarioContext` | yes | **Read: `const`, restricted to the event's footprint. Write: through a layer writer of the `ScenarioContext`** (§9.2) | Gains the **realization record** `{uid → realization}`, which today lives on events as `new_index`, `realization_indices`, `current_realizations_index_vec`, `deletion_value` and `alignment_offset_p`. Also gains the remaining reach of neighbouring ends (today's `pending_`, computed by the engine from what is still undecided) and the per-scenario adjacency of multi-path layouts (§4). Counters and error models read the record instead of `friend` pointers. |
| `ExplorationContext` | **no: it is the DFS engine's own context** | Not passed to `evaluate`. A DFS kernel reaches it only through the continuation (§10.3) | Contents unchanged (`downstream_proba_map`, `seq_max_prob_scenario`, threshold, `index_map`, `next_event_ptr_arr`, `safety_set`). Formally, it becomes the exploration state of `LegacyDfsEngine`. A DP engine has its own (message buffers); a sampler needs none. `pruning_mismatch_floor` is dead. |
| `AccumulationContext` | yes | Not passed to `evaluate` | Unchanged, apart from the handler-owned accumulators (§6). `error_rate->get_err_rate_upper_bound` is a model query and moves to the bound service. It is the only reason events touch this context today (Genechoice.cpp:281, 588; Deletion.cpp:416). |

### 9.2 What a realization produces, in Context terms

One realization of one event produces exactly **one layer of each context it can write**, plus a status:

```
struct EventOutcome {                                  // the "Delta", or "subscenario"
    Status                         status;              // ok | infeasible(reason) -- today a silent `continue`
    ScenarioContext::Layer         scenario;            // offsets, sequence handle and mismatch handle of the
                                                        //   footprint segments, the probability factor,
                                                        //   and the realization record entry {uid -> r}
    ExplorationContext::BoundFacts bound_facts;         // junction distances fixed, surviving-core counts:
                                                        //   what today's code turns into downstream_proba_map slots
};
```

- `ScenarioContext::Layer` is what an event writes today at its claimed layer: the segment write triple G4 and its `scenario_proba` factor.
- `BoundFacts` are the *inputs* of what an event writes today into `downstream_proba_map`. The bound service turns them into factors.
- The realization record entry is what `add_to_marginals` and the counters read at the leaf today.

**One event code path, two engines.** The event writes its outcome through a writer that has the shape of `ScenarioContext::Layer`, and the engine chooses which writer:
- **DFS-A** passes a writer bound to the layered `ScenarioContext`, at the event's claimed layer. Writes happen in place: **zero-copy, as today**.
- **DFS-B, the DP engines and the samplers** pass a writer into a detached `ScenarioContext::Layer`. The engine splices it in at a layer it owns, or stores it as a table row.

"Zero-copy versus subscenario" thus becomes the engine's choice of writer, not a property of the event (§10.3).

`CONTEXT_BASED_ITERATE.md` is out of date and should be refreshed as part of this work:
- It puts `index_map`, `downstream_proba_map` and `safety_set` in `ScenarioContext` (L116-118).
- It lists `seq_likelihood` and `seq_probability`, which do not exist, in `AccumulationContext` (L186-187).
- It uses the old `allowed_realizations` and memory-layer API (L215, L446-458).

## 10. Zero-copy legacy DFS versus orchestrator-driven DFS

### 10.1 How zero-copy works today

- **Layers.** `LayeredArray` holds one slot per (key, layer). Layers are claimed in queue order at initialisation (`request_layer`), and every write requires a claim (R3b, LayeredArray.h:216-266). An event's layer for a key is the number of earlier events in the queue that claimed that key.
- **Nothing is ever undone.**
  - `restore_layer` is never called during the search.
  - A sibling realization overwrites the same layer.
  - On return, deeper layers simply go stale.
- **Why stale values stay invisible:**
  - reads that name `layer-1` explicitly;
  - `set()` resets `current_layer_`;
  - bounds read through a snapshot of the layers claimed at initialisation (Rec_Event.cpp:322-324).
- **Maps hold pointers, not copies.** The sequence and mismatch maps store **pointers into buffers owned by events**, or into the alignment data, and those buffers are rewritten in place.
- **`scenario_proba` is a shared reference.** Ancestors that loop restore it; Insertion and Dinucl_markov rely on that.
- **The call stack *is* the DFS stack,** through `iterate_wrap_up` and `next_event_ptr_arr`.

**Invariants this imposes:**
- The layer index is fixed by queue position.
- Every claimed layer is written before handing off.
- Reads of the current layer assume no later writer.
- `pending_`, `partner_chosen` and the suffix of the junction profiles depend on queue position.
- A pointer lives for one realization.
- Event members are part of the scenario state.

### 10.2 What the test harness reveals

The harness is `RecordingEvent` (tst/igor/Core/test_utils.h:550-618, test_utils.cpp:196-448; ITERATE_TEST_GUIDE.md).

- **How it works.** A fake successor in `next_event_ptr_arr` (id 31, within `max_events` = 32) takes a snapshot of the state by deep copy, because `iterate` returns `void`.
- **A unit test has to assemble a mini-model:**
  - a registry, and `set_adjacent_segments` called by hand;
  - a queue, with `initialize_event` in order;
  - the bound fold over downstream events (an Insertion cannot even initialise without its Dinucl_markov);
  - layer-0 claims for the preset neighbours.
- **What the harness cannot observe:**
  - *why* a realization was dropped: overlap, guard, prune and range all look the same, and Deletion's `break` cannot be told from its `continue` (ITERATE_TEST_GUIDE.md:208-211);
  - `scenario_upper_bound_proba`;
  - `index_map`;
  - the state in event members, except through `add_to_marginals` after a single hand-off;
  - anything at the leaf.
- **The leaf pulls in the whole model,** including `Single_error_rate`, which dereferences V and J unconditionally.

### 10.3 Two DFS designs on one semantic core

**DFS-A: push, legacy, zero-copy.** The event keeps its loop, but `iterate_wrap_up` gives way to a **continuation provided by the engine**:

```
iterate(const QuerySequenceContext&, const ModelContext&, ScenarioContext& /* writer at own layer */, Next& next)
```

- The engine binds the `ScenarioContext` writer to this event's claimed layers. It is the same writer interface as `evaluate` uses (§9.2), bound in place.
- `next(p, bound_facts)` runs the engine's prune check, then descends or handles the leaf.
- `next.reject(reason)` reports a dropped realization.

| Keeps | Gains | Costs |
|---|---|---|
| Zero-copy (layers and in-place buffers). The call stack as the DFS stack. The candidate order and Deletion's `break`. Bitwise arithmetic. | The engine owns dispatch, leaf, counters, accumulation and the pruning policy. Tests replace `RecordingEvent` with a lambda continuation: no fake event, no id limit, no leaf, and the reasons for drops become observable. | Layer *indices* still come from queue order (now owned by the engine, still fixed per engine plan). Events still hold scratch buffers, so the per-thread `copy()` remains. |

**DFS-B: pull, driven by the orchestrator.** The orchestrator calls `candidates(...)`, then `evaluate(r, …)` for each candidate into a detached `EventOutcome` (§9.2). It splices the outcome into the layered contexts it owns, asks the bound service, and recurses.

| Gains | Costs and requirements |
|---|---|
| Events are const and hold no state: shared across threads, no per-thread `copy()` or `finalize()`. They are tested as pure functions (contexts in, outcome out): no queue, no layers, no `RecordingEvent`. The engine may **reorder events per read** (most constrained first) and **memoise** `evaluate` and `prepare`. The **same** `evaluate` builds the DP tables and drives posterior sampling, so one semantic implementation serves every engine. Reasons for drops arrive in `status`. | The outcome must not point into mutable event scratch. Dinucl_markov mutates its buffer after publishing the pointer, and Gene_choice reuses one `placement_mismatches` pointer. Content must therefore be views into immutable per-read data: a template `string_view` plus a cut; mismatch lists as **subranges of the alignment's sorted list** (ITERATE_GENERIC L1501-1513); a palindrome as (template, extension), computed in `prepare`. Otherwise it needs a per-depth arena. The per-candidate overhead (one virtual call, plus building and splicing the outcome, for about 20 Deletion candidates per call) **must be measured**. Order-dependent pruning has to be declared (a sorted order plus a monotonicity flag). Gene_choice's three enumeration modes produce heterogeneous candidates (placement × realization). Bitwise equality is possible only if the order of multiplications is reproduced, which is fragile. |

**Recommendation:**
1. Extract each event type's semantics into pure functions: the event-specific part of S, then K, the structure of the factor, and the bound facts.
2. DFS-A's `iterate` becomes a thin kernel of `LegacyDfsEngine` that calls those functions and writes in place, keeping bitwise equality. This makes concrete the idea that "engines own per-event-type kernels": kernels are registered per (engine, event type), instead of adding one virtual per engine to `Rec_Event`.
3. DFS-B and the DP engine call the same functions through `evaluate`, with a detached-layer writer. The semantic code is shared; only the writer differs.
4. Measure DFS-B against DFS-A on TRA and TRB before choosing a default.

## 11. Defects found during the analysis

| Defect | Location |
|---|---|
| `MarkovSamplingHandler` indexes rows differently in `precomputeCDF` (`ps * n_from + from`) and in `sample` (`from_state * m_parent_slice_count + ps`). The two agree only without parents. | TK MarkovSamplingHandler.tpp:85 and :151 |
| In Viterbi mode, `add_to_marginals` overwrites the current entry without clearing earlier bests, so the result is not one-hot after normalisation. | Genechoice.cpp:910-916, Deletion.cpp:892-899, Insertion.cpp:404-411 |
| `Pgen_counter` in sequence-output mode writes a header but no rows. | Pgencounter.cpp:181-188 |
| `Errors_counter` divides by `n_gen` without a guard. | Errorscounter.cpp:169-170 |
| `Coverage_err_counter` throws on D. | Coverageerrcounter.cpp:314-316 |
| `Hypermutation_global_errorrate` allocates GSL vectors at each Newton iteration and never frees them. | Hypermutationglobalerrorrate.cpp:1342-1351 |
| The pruning bound `mu/(3(1+mu))` of `Hypermutation_global_errorrate` is not a bound when some `S_k > 1`. | Hypermutationglobalerrorrate.cpp:297-327 |
| Comments and README in the TK handlers describe the old axis order (own axes first). The code puts parents first. | TK handlers, `Model/README.md` |

## 12. Suggested order of work

1. **Indexed CDF sampler.** One deliberate commit that moves the generation goldens.
2. **The UID replaces the generated name,** which fixes the tandem-D collision.
3. **Scenario record in `ScenarioContext`;** the `friend` raw pointers are removed.
4. **Generic JSON constructor for events, and `to_json`.** The v1 and v2 readers become importers.
5. **Registry/layout split;** layout DAG with guards; adjacency per scenario.
6. **Split `GenModel::infer_model` into the generic EM driver and `LegacyDfsEngine`** (§5.2), with `LegacyDfsEngine` as the only engine.
   - DFS-A: a continuation replaces `iterate_wrap_up`, and the bound fold and the layers move into the engine.
   - The goldens stay bitwise.
7. **Pure semantics per type and `prepare(read)`;** `evaluate` and `EventOutcome`; a DFS-B prototype, then the measurement.
8. **Handler contract with external accumulators.**
   - The error models become a Probe event with `Categorical`, `Logistic` or `Bernoulli` handlers.
   - The code is moved, not rewritten, so the goldens stay bitwise.
9. **`ChainDPEngine` prototype on the `prepare` caches** (proposal §8, step 1), validated by check 1. Then `DpBoundedDfsEngine` for the counters.

Steps 4 and 7 take the file layout of §14 and §15 and the tier annotations of §13 as their definition of done. The parameter-store writer of §15 (round-trip precision) is its own golden-moving commit, placed next to step 1.

## 13. Computation tiers

### 13.1 What exists

ITERATE_GENERIC §2.10 (G10) defines three **computation tiers** by what a computation depends on:

| G10 tier | Depends on | Runs | Hook today |
|---|---|---|---|
| 0 — model | model parameters, event ordering, `processed_events` | once per model | `initialize_event()`, `initialize_Len_proba_bound()` |
| 1 — read × event | tier 0 + the read and its alignments | once per read | **none** |
| 2 — scenario | tier 1 + what upstream events wrote | once per scenario branch | `iterate()` |

G10's placement rule:
- Every computation is annotated with its tier.
- No tier-1 computation reads `ScenarioContext`.
- The hook `prepare_for_query(const QuerySequenceContext&)` is proposed but not built. Nothing is hoisted yet, because hoisting changes evaluation counts and needs its own benchmark.

**Name collision.** ITERATE_GENERIC §2.5 (L625-671) and the Phase A status of the REC_EVENT_CAPABILITY plan also say "tier 1/2/3" for the three *checks* of the segment-completeness invariant:
- capability at `finalize()`;
- the leaf assert;
- the hand-off check in `RecordingEvent`.

Those are check levels, not computation tiers. This document proposes calling them **check levels C1–C3** from now on, and keeping "tier" for computation.

### 13.2 Refining the tiers

G10's tiers 0 and 1 each mix two cadences. A computation depends on some subset of:
- **S**, the structure: events, layout, conditioning, engine order;
- **θ**, the parameter values;
- **R**, the read: its sequence, alignments and covariates;
- **P**, the path: what upstream events decided on the current branch.

Five tiers follow, each compatible with G10's numbering:

| Tier | Depends on | Cadence | Shared by | Driver / engine hook | Event-side hook | Examples today |
|---|---|---|---|---|---|---|
| **0a structure** | S | once per model load | all threads, all iterations | model construction; `check_admissible` and the engine plan | footprint queries; the JSON constructor | adjacency and layout, `SegmentIndex`, queue order, layer plan and claims, `pending_` reach (depends on order), `partner_chosen`, flank-check wiring, junction-slot resolution, tensor shapes |
| **0b parameters** | S + θ | once per EM iteration (once per run for sampling) | all threads | `prepare_iteration` | none: the handlers materialise | bound fold (junction profiles), `dinuc_proba_matrix` (`update_event_internal_probas`), sampling CDFs, handler `bound()` tables, N-mer mutation probabilities (`update_Nmers_proba`) |
| **1a read** | S + R | once per read, **valid across iterations** | immutable, per read | `prepare_read` → `ReadCache`, which the driver may keep | `prepare(read)` | alignments; clipped templates per alignment; palindrome strings and their mismatches; endogenous mismatch counts; D's sliding mismatch profile; mismatch subranges per deletion; SEGMENT §5.2's segment tables (offsets, end nucleotides, match and mismatch counts) |
| **1b read × parameters** | S + θ + R | once per read per iteration | one thread's workspace | prologue of `process_read` | none (engine and handler) | per-read handler rows with covariates; local likelihoods from stored counts (`(1-e)^m · e^k`); the surviving-core error bound; DP forward and backward messages; the bounds of the DP-bounded DFS; initial `seq_max_prob` |
| **2 scenario** | S + θ + R + P | once per branch | the path (stack) | the engine's inner loop | `evaluate` / kernel | offset and sequence writes, the factor, bound facts, overlap verdicts, Dinucl_markov's chain over junction content |

Mapping: G10 tier 0 = 0a + 0b; G10 tier 1 = 1a + 1b; G10 tier 2 = 2.

### 13.3 Why the split matters

- **0a vs 0b: only 0b changes between EM iterations.**
  - Today `infer_model` redoes all of 0a inside the iteration loop, on every thread: the copies of `Model_Parms`, `initialize_event`, the queue and `next_event_ptr_arr` (GenModel.cpp:292-409).
  - 0a could run once per run.
- **1a vs 1b: only 1a can be cached across iterations.**
  - SEGMENT §5.2's "computed once per read and reused by every EM iteration" is exactly 1a.
  - Its "recomputed from the stored counts when `e` changes" is exactly 1b.
- **Placement rule, extending G10.**
  - A tier-k computation reads only inputs of tier ≤ k.
  - No tier-0 or tier-1 computation reads `ScenarioContext`.
  - No 0a or 1a computation reads parameter values.
  - The rule is testable: a 1a function takes only `(const Event&, const QuerySequenceContext&)`.
- **Determinism.** Moving a pure function to a lower tier changes how often it runs, not its result. The goldens stay bitwise if the arithmetic is unchanged. G10's caution still holds: each hoist needs its own benchmark.
- **Memory.** A 1a cache for every read may not fit. The driver chooses a policy: keep all, keep within a budget, or recompute. §13.5 estimates the sizes and recommends recomputing by default.
- **Threads:**
  - 0a and 0b are shared and const;
  - 1a is immutable per read, so it can be shared;
  - 1b lives in one thread's workspace;
  - 2 lives on the path's stack.
- **Generation** has 0a, 0b (the CDFs; `FastGenerator::initialize` is already a model 0b hook, ITERATE_GENERIC L1495) and a per-sample tier 2 (`sample` then `apply`). It has no tier 1, because there is no read.

### 13.4 How the tiers line up with the rest

| | 0a | 0b | 1a | 1b | 2 |
|---|---|---|---|---|---|
| Context | `ModelContext` (structure) | `ModelContext` (handler rows) | `QuerySequenceContext` + `ReadCache` | engine workspace (DFS: the per-read part of `ExplorationContext`) | `ScenarioContext` (+ the DFS's layered `ExplorationContext`) |
| Engine call | `check_admissible`, plan | `prepare_iteration` | `prepare_read` | `process_read` prologue | `process_read` inner loop |
| Event method | footprint, constructor | — | `prepare(read)` | — | `candidates`, `evaluate` / kernel, `apply` |
| Handler method | shape, family | `maximize`, materialisation, CDFs | — | `distribution(context)` with covariates, `bound` | `distribution` lookup |

An event has no tier-1b method. Tier 1b is where parameters meet the read, and it belongs to the engine and the handlers. This is one more reason the probability bounds leave the event (§7.2).

### 13.5 What can be cached, and how much storage it takes

**Calibration** (from the tree; everything else is an estimate):
- **Reads and alignments.** The regression run `default_evaluate` uses the human TRB test model, with 300 reads of 60 nt. The aligner keeps, per read:
  - **22 V alignments** on average (max 64);
  - **12 J** (max 19);
  - **39 D** (max 90).

  It explores **1,226 complete scenarios** per read on average (max 9,349); see `inference_logs.txt` and `aligns/*.csv` in `regression_reports/`.
- **Domains.** Human TRB: V 89, D 3, J 15, deletions 21/21/21/23, insertions 31. Human IGH (`models/human/bcr_heavy`): V 97, D 35, J 7, deletions 27 (V), 44 and 44 (D), 27 (J), insertions 61.
- **IGH reads.** Estimated at about 400 nt, with 40 V, 100 D and 7 J alignments. There are no IGH alignment data in the tree.
- **`Alignment_data` (Aligner.h:83).** About 230 B fixed (a string, six vectors, offsets, score), plus the heap: the IMGT header as `gene_name` (about 100 characters) and the mismatch vector (about 8 entries of 8 B in TRB). That is roughly **400 B per alignment**, so about **30 kB per TRB read**. That memory is held today for every read for the whole run, in the `sequences` vector of `infer_model`.

| # | Value | Tier | Engines | Representation | Per read, TRB | Per read, IGH (est.) |
|---|---|---|---|---|---|---|
| 1 | Alignments (exist today) | input of 1a | all | `Alignment_data` | ~30 kB | ~90 kB (more mismatches) |
| 2 | Mismatch prefix profile over the range deletions can reach (plus palindrome positions) | 1a | all | `uint16` per position per alignment: V 21+4, J 23+4, D template + 4 per end | ~3.6 kB | ~11 kB |
| 3 | Clipped template extents (start, end; a view into the shared library template) | 1a | all | 12 B per alignment | ~0.9 kB | ~1.8 kB |
| 4 | Palindrome extensions (reverse complement of ≤ 4 nt per end) | 1a | all | derivable from the template; ≤ 8 B per alignment if stored | ≤ 0.6 kB | ≤ 1.2 kB |
| 5 | N-mer context index per covered position per alignment (hypermutation models only) | 1a | all | `uint16` (4^5 = 1024 contexts) | — | ~31 kB |
| | **Core 1a, total** (items 2–5) | | | | **~5 kB** | **~13 kB** (single error rate), **~45 kB** (N-mer) |
| 6 | Dense segment table (SEGMENT §5.2): per (alignment, deletion[s]), two offsets, end nucleotide, match and mismatch counts | 1a | ChainDP | ~12 B per entry. TRB: V 462, J 276, **D 39 × 21 × 21 = 17,199** entries. IGH: **D 100 × 44 × 44 = 193,600** | ~215 kB | ~2.3 MB |
| 7 | Mismatch subranges per (alignment, single deletion); D composed from (d5) then (d3) | 1a | LegacyDfs | 2 × `uint16` per entry | ~9.5 kB | ~40 kB |
| 8 | Powers `e^k`, `(1-e)^m` (single error rate) | **0b** | all | 2 × read length doubles, per process | ~1 kB per process | ~6 kB per process |
| 9 | Junction chain: prefix log-products of the Markov kernel along the read, per junction type (the seed nucleotide only changes the first factor) | 1b | ChainDP, DP-bounded DFS, FFBS | 8 B per position per junction type | ~1 kB | ~6.4 kB |
| 10 | Local likelihood per dense table entry | 1b | ChainDP | 0 if computed on the fly from items 2 and 8; 8 B per entry if materialised | 0 – 143 kB | 0 – 1.5 MB |
| 11 | Hypermutation log-likelihood, prefix along positions per alignment | 1b | all | 8 B per (alignment, position) | — | ~125 kB |
| 12 | DP messages for the gene tuple in flight: (read offset × end nucleotide) states per boundary, forward and backward | 1b | ChainDP, DP-bounded DFS, FFBS | 4 boundaries × 4·L states × 2 × 8 B | ~15 kB | ~100 kB |
| 13 | Handler rows for covariate-dependent nodes | 1b | all | one row per such node | < 1 kB | < 1 kB |
| 14 | Surviving-core error bound per alignment end; initial `seq_max_prob` | 1b | LegacyDfs | 8 B per alignment end | ~1 kB | ~2 kB |

**What follows.**
- **Tier 1b is never cached beyond the read.** It changes every iteration, and each read needs it once. It lives in the thread's workspace:
  - TRB: ≤ ~250 kB per thread, even with the dense table (items 6, 10, 12);
  - IGH: about 1–4 MB per thread;
  - so at 64 threads, under about 256 MB in total.
- **Tier 1a pays off inside the read, not across iterations.**
  - **Within a read**, it hoists work out of the scenario loop. Today, a Deletion deep in the queue redoes its trimming once per upstream branch that reaches it. That repetition factor is in the hundreds, judging by the 1,226 complete scenarios per read.
  - **Across iterations**, persisting 1a only saves recomputing it once per EM iteration (about 10). Computing core 1a costs about Σ alignments × reachable positions ≈ 73 × 25 ≈ 2·10³ elementary operations per TRB read. That is small next to an E-step: more than 10⁴ event evaluations in the DFS (1,226 complete scenarios × ~10 events, before counting pruned branches), and about 2.5·10⁶ operations in the DP (SEGMENT §5.3).
  - Persistence buys a few per cent at most, and costs, for 10⁵ reads:
    - **0.5 GB** for core 1a;
    - **21 GB** for the dense TRB table (**230 GB** for IGH).
- **Recommendation:**
  - Compute 1a at the start of `process_read`, into the thread's workspace.
  - Persist core 1a across iterations only under a memory budget (as `REDUCTION_MEMORY_BUDGET_BYTES` does for the reduction), and only if a measurement shows that it pays.
  - Never persist engine-specific dense tables or tier 1b.
- **The dense DP table is a speed choice, not a necessity.** With item 2, every entry of item 6 is computable in O(1). A DP engine can therefore materialise only the rows it is about to use.
- **The real memory lever is item 1.** The alignments already cost about 30 kB per read today, for every read, for the whole run. Compacting them would bring that to about 100 B per alignment, i.e. about 7 kB per TRB read, which saves more than any 1a cache costs:
  - a library index instead of the IMGT header string;
  - `uint16` instead of `size_t`;
  - no vectors for the usually empty extended fields.
- **Is caching specific to an engine? Partly.**
  - **Core 1a (items 2–5) is event semantics.** It is produced by the events' `prepare(read)` and shared by every engine.
  - **Derived tables (items 6, 7) belong to each engine.** They are built on top of core 1a in that engine's `prepare_read`.
  - **Tier 1b belongs to the engine and the handlers.**
  - In short: the event defines *what is known* about a read; the engine decides *what it materialises* and *for how long*.
- **Tandem D.** D1 and D2 draw from the same library, so they share the D alignments and the core 1a computed from them: once per library, not once per event. This is another place where library and `seq_type` differ (§2.1).
- **Caveats.** The TRB figures come from 60-nt test reads; Adaptive TRB reads are 87 nt, so about 45 % more positions. The alignment counts depend on the aligner's thresholds. The IGH figures are estimates. Once `prepare` exists, a per-read size counter will replace this table.

## 14. Model file: how entries are keyed, and how the constructor is chosen

### 14.1 The question

In SEGMENT §3, `"events"` is an object keyed by nickname, with `type` as a property:

```json
"events": { "v_choice": { "type": "GeneChoice", "realizations": [ ... ] } }
```

Is it a problem to call the right constructor afterwards? Is the type just a property?

The type is a **discriminator**: a property whose value selects the decoder. This is the standard "internally tagged union":
- JSON Schema `oneOf` with `discriminator` in OpenAPI 3;
- serde's `#[serde(tag = "type")]`;
- an nlohmann `from_json` that dispatches on `j.at("type")`.

Every C++ JSON library in common use parses into a DOM first, so the whole object is in memory before decoding, and reading `type` first costs nothing:

```cpp
auto event = EventFactory::create(j.at("type").get<std::string>(), j);   // unknown type -> error listing known types
```

**Dispatch is therefore not the problem.** The type cannot be the primary key, because several events share a type (three `GeneChoice`). Grouping by type would be a hierarchy, not an identity.

### 14.2 Order: where it is content, and where it is not

The order in which events are listed carries no meaning, and it should not:
- **Priority and order belong to the engines.** Every order an engine uses (topological order, DFS queue, DP chain) is the engine's decision (§5).
- **Links between objects use the nickname.** These are conditioning, layout guards, tying and parameters. Converting nicknames to dense indices is internal.

Today the opposite holds. `Model_Parms::get_model_queue` sorts with a stable `std::list::sort` (Model_Parms.cpp:474), so events of **equal priority keep their order in the file**. Swapping two equal-priority events in `model_parms.txt` changes the DFS queue, and with it the layer plan, the bound fold and the last bits of the results. The new rule removes that dependency:
- **No engine reads the order of the file.**
- **Ties are broken by nickname.**

Order is content in exactly four places:

| Where | Why it is content |
|---|---|
| Path segments in the layout (`segments`, `edges`) | They define the geometry |
| `values` of a parameter tensor | The layout of the numbers, row-major over the declared `axes` |
| The `axes` list of a tensor | Which axis is which |
| Realization coordinates | They fix the summation order inside a tensor (normalisation, marginalisation). Numerically this is visible in the last bits, so for exact reproducibility the coordinate must be fixed by the file's **content**, not by how a parser iterates it (§14.4) |

With that settled, **an object keyed by nickname is the better layout for events**, as SEGMENT §3 proposes:
- the nickname is the key, so uniqueness is part of the structure;
- references read naturally;
- it is what a modeller edits.

What remains is one loader requirement and one internal rule:
- **Duplicate keys are rejected.** RFC 8259 only says names *should* be unique, and `nlohmann::json` silently keeps the last one. The loader detects duplicates through the parser callback, which sees every key, and fails with both locations.
- **Member order is never read.** `nlohmann::json` stores objects in a `std::map` (alphabetical). `nlohmann::ordered_json` keeps the order of the file, but other tools that rewrite the file (Python, `jq`, editors) need not. Either way the loader must not depend on it.
- **Internal uid.** Dense indices assigned in **byte-wise nickname order**: a pure function of content, identical whatever the parser and whatever the layout of the file.

### 14.3 Options for the events container, revisited

| | Layout | Verdict |
|---|---|---|
| **A** | Object keyed by nickname, `type` inside (SEGMENT §3) | **Recommended**, with duplicate rejection and the uid rule above |
| **B** | Grouped by type: `{"GeneChoice": {"v_choice": {...}}, ...}` | Rejected: the kind becomes a position (the "filed" problem of §2); nickname uniqueness must be checked across groups; every new kind adds a group |
| **C** | Array of objects with `nickname` and `type` | Equivalent once order carries no meaning, but it gives up structural uniqueness of the nickname, and its order would *look* meaningful to a reader |
| **D** | Externally tagged array | Rejected: verbose, awkward to edit |

The other rules stand:
- **Every polymorphic object is tagged by `type`:** events, parameterizations, guards, layout, Probe, covariates.
- **Every object that can be referenced is keyed by its nickname,** in one namespace shared by events and covariates.
- **Two registries,** `EventFactory` and `ParameterizationFactory`. Each registration provides the constructor from JSON, `to_json`, a JSON Schema fragment, and the contract description of §2.3.
- **A two-phase load:**
  1. Parse, rejecting duplicates; validate against the schema assembled from the registries.
  2. Construct the events, assigning uid by nickname order; resolve references; run Phase C validation.

### 14.4 Realizations: identity, index, and events without parameters

A realization has two different things that are easy to conflate:
- **Its identity: the label.** This is the gene name, the deletion or insertion value, or the nucleotide. It is unique within its event, and always required.
- **Its index: the coordinate in the parameter tensor.** It is only meaningful when parameter values exist.

There are three ways to link parameter values to realizations:

| | Link | Robust to hand edits? | Cost |
|---|---|---|---|
| (i) | Position in the realization list | No: reordering two lines silently permutes the parameters | — |
| (ii) | Explicit `index` property on the realization (Q. Marcou's proposal; the legacy files' `;57`) | Yes to reordering. Adding a realization means choosing an index; removing one leaves a hole | Redundant with position; needs validation (dense `0..n-1`, unique) |
| (iii) | **Parameters keyed by label**: each tensor axis lists the labels of its coordinates | Yes to reordering, adding and removing. A mismatch is detected and named | Labels repeated in the parameter section |

**Recommendation: combine (ii) and (iii) under one "index authority" rule.**

1. **The label is the identity** everywhere a realization is referenced, including in parameter sections (§15.4).
2. **The internal index is decided by the structure**, in this order of authority:
   1. an explicit `index` on the realization, if present;
   2. otherwise the next free index, in canonical label order: numeric for integer-valued realizations, byte-wise for names.

   An event may mix the two. Realizations given an `index` keep it; those without one are appended. This is what makes hand-adding easy: add a line with no index, and `igor model fmt` assigns it.
3. **Indices must be dense and unique** after assignment. A hole left by a removal is an error, and `fmt --renumber` closes it, rewriting the parameter section accordingly.
4. **Parameters carry labels per axis.**
   - On load, values are permuted from the file's label order into the internal index order, matched **by label**.
   - A parameter file therefore never depends on the structure's indices, and the reverse.
   - The canonical writer emits labels in internal index order, so a round trip is the identity.
5. **Realizations without parameters need no index.**
   - This covers skeleton models that are not yet fitted, and fixed domains (Dinucl's `A,C,G,T`).
   - Ranges are allowed as shorthand: `{"range": [-4, 16]}` expands in ascending order.
   - This is Q. Marcou's "no index without tensor values", generalised: **an index exists exactly where values exist.**
6. **A realization present in the structure but absent from the parameters** (the hand-added case) is handled by an explicit load policy:
   - `on_missing = error` (the default);
   - `pseudocount:x` or `uniform`, with renormalisation of the affected rows.

   The policy and the rows it touched are recorded in `provenance`. This keeps the normalisation caveat visible instead of silent.
7. **A realization present in the parameters but absent from the structure** (a removal) is an error by default. With `on_extra = drop`, its mass is dropped and the rows are renormalised, and that is recorded too.

**Legacy files carry explicit indices,** and the importer keeps them. The internal tensor layout is therefore the legacy layout, which is MERGE_ANALYSIS's invariant, and the legacy engine stays bitwise.

### 14.5 Editable by hand and exactly reproducible

These two goals conflict only if one form has to serve both. The answer is **two forms of one schema**, with tools that convert between them:

| | **Authoring form** (written by the modeller) | **Canonical form** (written by IGoR) |
|---|---|---|
| Shorthand | Allowed: `segments` block (§2.3), deletion `range`, linear `layout` | Expanded |
| Defaults | May be omitted | **Written explicitly** (with `schema_version`), so a model keeps its meaning if a default changes in a later release |
| Indices | Optional (§14.4) | Present wherever parameters exist |
| Member order | Free | Keys sorted (the RFC 8785 ordering, pretty-printed) |
| Array order | Free for realizations (indices decide) | Realizations by index; labels in index order |
| Comments | JSONC accepted on input (nlohmann's `ignore_comments`) | Comments are dropped, so the schema provides a `description` field on every object, which survives |
| Numbers | Any valid literal | Integers exact. Floats as the shortest decimal that round-trips for the stored type (`std::to_chars`): about 17 digits for `double`, 21 for x86 `long double` |

**Contracts, each a test:**
- `fmt(fmt(x)) == fmt(x)`: canonicalising is idempotent.
- `load(fmt(x))` is **bitwise** the same in-memory model as `load(x)`.
- For a canonical `c`, `save(load(c)) == c`, **byte for byte**.
- The **structure fingerprint** is computed on the canonical form, after removing `description` fields. Cosmetic edits (whitespace, key or line order, comments) leave it unchanged; semantic edits change it.

**The floating-point type is a decision, not a detail.** The legacy runtime holds marginals as `long double`. The default `nlohmann::json` parses numbers into `double`, which loses the extra bits of an 80-bit value. There are three options:
- configure `basic_json` with `long double` as its float type;
- store values as strings at 21 digits;
- decide that `double` is the storage type, and move the goldens once.

Whichever is chosen, the reader and the writer must use the same type, or the second contract fails.

**Tools:**
- `igor model fmt`: canonicalise, assign missing indices, `--renumber`.
- `igor model validate`: schema plus Phase C, with messages by nickname and label.
- `igor model diff`: a semantic diff, by nickname and label, independent of order.
- `igor describe-event <type>`.
- A published JSON Schema, generated from the registries, so editors validate and autocomplete.

**Authoring form, as a modeller would write it:**

```jsonc
{
  "schema_version": 3,
  "structure": {
    "events": {
      // order of entries is irrelevant; the nickname is the key
      "v_choice":  { "type": "GeneChoice", "library": "TRBV", "targets": [ { "segment": "V" } ],
                     "realizations": [ { "label": "TRBV5-8*01" }, { "label": "TRBV7-2*01" } ] },
      "v_3_del":   { "type": "Deletion", "targets": [ { "segment": "V", "side": "3'" } ],
                     "realizations": { "range": [-4, 16] } },
      "vd_dinucl": { "type": "DinucMarkov", "targets": [ { "segment": "VD", "side": "5'" } ] }
    },
    "layout": { "type": "Linear", "segments": ["V", "VD", "D", "DJ", "J"] },
    "conditioning": { "v_3_del": ["v_choice"] }
  },
  "library": { "$ref": "imgt_human_TRB_2026-07.json", "sha256": "…" }
}
```

**Canonical form of one event after a fit** (excerpt):

```json
"v_choice": {
  "description": "",
  "library": "TRBV",
  "parameterization": { "type": "Categorical" },
  "realizations": [ { "index": 0, "label": "TRBV5-8*01" }, { "index": 1, "label": "TRBV7-2*01" } ],
  "targets": [ { "segment": "V" } ],
  "type": "GeneChoice"
}
```

- Template sequences live in the `library` section, keyed by label, so the structure stays short enough to edit.
- `layout` is tagged: `Linear` is shorthand for a one-path DAG; `Dag` carries `segments`, `edges` and `guards` (§4).
- A missing `parameterization` means the default for the event's type.

## 15. Parameter store: where parameters are written, and partial exports

### 15.1 Today

- `model_parms.txt` holds:
  - the structure;
  - the germline templates;
  - **and the error model's parameter values**: `@ErrorRate` / `#SingleErrorRate` followed by the rate.
- `model_marginals.txt` holds the event parameters:
  - keyed by nickname (`@j_choice`), with a dimension line (`$Dim[69,7]`) and rows labelled by parent coordinates (`#[v_choice,k]`).
  - It is written at the stream's **default precision, 6 significant digits** (Model_marginals.cpp:1119-1160, no `setprecision`), from `long double` values in memory, and renormalised at load (Model_marginals.cpp:1235).
  - So a model reloaded from its own output is not the model that wrote it.
- Each EM iteration writes `iteration_N.txt` and `iteration_N_parms.txt`. The second repeats the whole structure, templates included.
- Size is not an issue: human IGL is about 23 kB of parameters-and-structure plus 18 kB of marginals.

### 15.2 What a model contains, by rate of change

| Section | Content | Changes when | Shared by |
|---|---|---|---|
| `structure` | events (domains, targets, library references), layout, conditioning, covariate declarations, parameterization families and their hyperparameters (N-mer size, count threshold) | the modeller edits | every fit of that model |
| `library` | germline templates (names, sequences) | a reference release (IMGT) | every model of that locus and species |
| `parameters` | values per node: tensors, logistic coefficients, error rate | every EM iteration, every dataset | — |
| `provenance` | software version, date, input-data hash, iterations, likelihood trace, parent model | each fit | — |
| `accumulators` (optional) | sufficient statistics per node | each E-step | — |

### 15.3 Options

| | Placement | Pros | Cons |
|---|---|---|---|
| **A** | Everything in one JSON (SEGMENT §3, "same document") | Self-contained, one file, atomic | Structure and numbers churn together (diffs, version control); each iteration repeats the structure; the library is duplicated in every model |
| **B** | Two files, structure and parameters (today's split, modernised) | One structure, many parameter sets (iterations, bootstraps, individuals); clean diffs | Pairing errors, unless a fingerprint ties them |
| **C** | **One sectioned document; each section inline or referenced:** `{"$ref": "file", "sha256": "..."}` | One schema and one loader for every split; self-contained by default; split on demand | Reference resolution and hash checks to implement |
| **D** | Binary container (HDF5, or a zip of `.npy` = `.npz`) | Exact numbers, large tensors, direct numpy/h5py access | A dependency (HDF5 is heavy); not editable by hand |

**Recommendation: C.** Numeric payloads are inline JSON by default. A `.npy` sidecar is allowed for large tensors, later, if a model ever needs it.

### 15.4 The `parameters` section

```json
"parameters": {
  "structure_fingerprint": "sha256:...",
  "nodes": {
    "j_choice": {
      "axes": [ { "event": "v_choice", "labels": ["TRBV1*01", "TRBV10-1*01", "..."] },
                { "event": "j_choice", "labels": ["TRBJ1-1*01", "TRBJ1-2*01", "..."] } ],
      "encoding": "decimal",
      "values": [ [0, 0.0672312, 0.0812928, "..."],
                  [0.2043, 0.152668, 0.1652, "..."] ] },
    "errors": { "type": "Logistic", "values": { "log_mu": -5.1, "e": [[ "..." ]] } }
  }
}
```

- **Keyed by nickname** (§14). The uid is runtime only.
- **`axes` lists the axes in order,** `[parents…, own]`, each with the **labels** of its coordinates.
  - The section describes its own layout and is matched to the structure **by label** (§14.4).
  - Reordering, adding or removing realizations in the structure is detected and handled by the explicit load policy, never silently.
  - An edit to `conditioning` that changes a node's parents is caught at load, because the axes no longer match.
  - This keeps the one strength of the legacy file, its `#[v_choice,k]` labels.
- **`values` are nested arrays,** one level per axis: a row per parent combination. They are readable, and the shape is implied and checked.
- **Compact variant for large tensors:** an axis may give `{"event": "v_choice", "fingerprint": "sha256:..."}` instead of its labels. It then refers to the structure's index order, and the fingerprint covers the (index, label) pairs. That axis then loses the protection of label matching, so the canonical writer uses labels unless asked otherwise.
- **`structure_fingerprint`** is a hash of the canonical structure section. It produces a warning in partial imports, where only the per-node checks are binding, and an error in strict mode.
- **Encoding.**
  - `"decimal"` with the round-trip precision of the stored type: `max_digits10`, i.e. 17 digits for `double` and 21 for the x86 80-bit `long double`. Use `std::to_chars` where it is available.
  - A file written and read back is then **bitwise** the same model.
  - `"hexfloat"` (`%La`) is an exact alternative, but the 80-bit format is not portable.
  - Changing the writer moves the goldens that compare marginal files, once, in a deliberate commit (§12).
- **Non-tensor parameterizations** write their own payload through their registered `to_json` (§14.3).
- **Accumulators use the same layout,** with `"kind": "accumulator"`. That enables:
  - a distributed E-step: shards export their accumulators, which are merged in a fixed order before the M-step;
  - checkpointing and restarting a fit.

### 15.5 Partial exports and imports

An export is a choice of sections × a choice of nodes × inline or referenced:

| Use | Export | Import |
|---|---|---|
| Share a finished model | all sections, inline | load |
| EM run output | `structure` + `library` once; per iteration, `parameters` (referencing the structure) + `provenance` | load any iteration with its structure |
| Model skeleton | `structure` (+ `library` reference), no parameters | the loader initialises (uniform, or by rule); start inference |
| Transfer one component (an error model fitted elsewhere; one person's V usage) | `parameters` for selected nodes only: a **fragment** | **overlay** onto a base model, with the per-node domain check; overlaid nodes optionally marked fixed (today's "fixed events") |
| Germline update | `library` alone | swap, then remap the affected nodes by realization name |
| Distribution or checkpoint | `accumulators` per shard | ordered merge, then `maximize` |
| Analysis | derived summaries (marginal P(V), P(insertion length)) | **never re-imported**: a report, not a model |

```cpp
ModelDocument      export_model(const RecombinationModel&, const ExportSpec& spec);  // sections, nodes, inline|ref, encoding
RecombinationModel load_model(const ModelDocument& base, const std::vector<ModelDocument>& overlays = {});
```

The CLI mirrors this:
- `igor model export --sections structure,parameters --nodes v_choice,j_choice -o fragment.json`
- `igor model overlay base.json fragment.json -o combined.json`

**Legacy mapping** (the importer of §7.1):
- `model_parms.txt` becomes `structure` plus `library`. Its `@ErrorRate` value moves to `parameters`.
- `model_marginals.txt` becomes `parameters`. Its nickname keys and `#[parent,k]` labels map directly onto `axes`.
- The importer computes the fingerprints.
- The `iteration_N` files become parameter documents that reference one structure.

## 16. References

All biology references below were checked against PubMed on 2026-10-01; the DOIs are the ones PubMed returns. The algorithm references are outside PubMed's scope and are **still to verify**.

**One correction from that check.** Tycko, Palmer & Sklar 1989 report **TRG–TRD** chimeras: Vγ-(Dδ)-Jδ and Vδ-(Dδ)-Jγ. It is not a TRGV–TRBJ example; §4.1 is corrected accordingly.

| Reference | DOI | Cited for |
|---|---|---|
| Reth, Gehrmann, Petrac & Wiese 1986, *Nature* 322:840 | [10.1038/322840a0](https://doi.org/10.1038/322840a0) | VH replacement in pre-B cells |
| Kleinfield et al. 1986, *Nature* 322:843 | [10.1038/322843a0](https://doi.org/10.1038/322843a0) | VH replacement; embedded heptamer TACTGTG, conserved in > 70 % of VH genes |
| Zhang et al. 2003, *Immunity* 19:21 | [10.1016/s1074-7613(03)00170-5](https://doi.org/10.1016/s1074-7613(03)00170-5) | VH replacement footprints in human IgH |
| Zhang, Burrows & Cooper 2004, *Immunol Rev* 197:231 | [10.1111/j.0105-2896.2004.0107.x](https://doi.org/10.1111/j.0105-2896.2004.0107.x) | Review of VH replacement; footprint in CDR3 |
| Tiegs, Russell & Nemazee 1993, *J Exp Med* 177:1009 | [10.1084/jem.177.4.1009](https://doi.org/10.1084/jem.177.4.1009) | Receptor editing |
| Gay, Saunders, Camper & Weigert 1993, *J Exp Med* 177:999 | [10.1084/jem.177.4.999](https://doi.org/10.1084/jem.177.4.999) | Receptor editing |
| Radic, Erikson, Litwin & Weigert 1993, *J Exp Med* 177:1165 | [10.1084/jem.177.4.1165](https://doi.org/10.1084/jem.177.4.1165) | Skewed Vκ/Jκ usage (Jκ5) after editing |
| Nemazee 2017, *Nat Rev Immunol* 17:281 | [10.1038/nri.2017.19](https://doi.org/10.1038/nri.2017.19) | Review: central tolerance, secondary recombination |
| Petrie et al. 1993, *J Exp Med* 178:615 | [10.1084/jem.178.2.615](https://doi.org/10.1084/jem.178.2.615) | Successive TRA rearrangements |
| Meek, Hasemann & Capra 1989, *J Exp Med* 170:39 | [10.1084/jem.170.1.39](https://doi.org/10.1084/jem.170.1.39) | D–D fusions, direct or inverted |
| Briney et al. 2012, *Immunology* 137:56 | [10.1111/j.1365-2567.2012.03605.x](https://doi.org/10.1111/j.1365-2567.2012.03605.x) | V(DD)J in about 1 in 800 circulating B cells |
| Safonova & Pevzner 2019, *Front Immunol* 10:987 | [10.3389/fimmu.2019.00987](https://doi.org/10.3389/fimmu.2019.00987) | Tandem D in CDR3s |
| Lafaille et al. 1989, *Cell* 59:859 | [10.1016/0092-8674(89)90609-0](https://doi.org/10.1016/0092-8674(89)90609-0) | P nucleotides |
| Marcou, Mora & Walczak 2018, *Nat Commun* 9:561 | [10.1038/s41467-018-02832-w](https://doi.org/10.1038/s41467-018-02832-w) | IGoR |
| Murugan, Mora, Walczak & Callan 2012, *PNAS* 109:16161 | [10.1073/pnas.1212755109](https://doi.org/10.1073/pnas.1212755109) | Generative model of TRB recombination |
| Tycko, Palmer & Sklar 1989, *Science* 245:1242 | [10.1126/science.2551037](https://doi.org/10.1126/science.2551037) | Inter-locus (TRG–TRD) chimeric rearrangements |
| Lefranc & Lefranc 2001, *The T Cell Receptor FactsBook* (book) | — (not in PubMed) | Up to three D in TRD |

**Not in PubMed; still to verify:**
- Eisner 2002 (expectation semiring);
- Li & Eisner 2009 (first- and second-order expectation semirings);
- Huang & Chiang 2005 (k-best);
- Carter & Kohn 1994, *Biometrika* 81:541;
- Frühwirth-Schnatter 1994, *J Time Ser Anal* 15:183;
- Scott 2002, *JASA* 97:337;
- Delyon, Lavielle & Moulines 1999, *Ann Stat* 27:94.
