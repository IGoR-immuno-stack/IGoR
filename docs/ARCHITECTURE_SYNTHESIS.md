# Architecture synthesis: model, engines and drivers

Status: summary, for discussion.
Date: 2026-10-01.

This document merges three designs into one target architecture:

| Short name | Document | What it is |
|---|---|---|
| **SEGMENT** | `SEGMENT_DECOMPOSITION_PROPOSAL.md` | T. Kloczko and Q. Marcou's discussion of 2026-10-01: events, segments, layout, conditioning graph; an exact per-segment DP for inference |
| **MODEL** | `MODEL_TOPOLOGY_REFACTORING.md` | tk's Model module, as built: `RecombinationModel`, `Topology`, handlers, engines, factories, `LegacyBridge`. A copy of `src/igor/Model/README.md` on `feature/tk_refactoring` |
| **REVIEW** | `SEGMENT_DECOMPOSITION_REVIEW.md` | Q. Marcou's remarks on SEGMENT, plus an analysis of the code. The reasoning for every change below is there; the "review §" column points to it |

`TK_REFACTORING_MERGE_ANALYSIS.md` (MERGE_ANALYSIS) is the background: the two graphs, the ownership points to agree on, and determinism.

Sections 4 and 5 list, item by item, what changes in SEGMENT and in MODEL. Changes are marked with a **Status**:

| Status | Meaning |
|---|---|
| *Kept* | Unchanged |
| ***Changed*** | Modified; the changed cell is in bold |
| ***New*** | Not in that document |
| ***Dropped*** | Removed |

## 1. Principles

1. **The model describes; engines compute.** Everything needed to define a probability distribution over scenarios is model data. That covers events, the layout, the conditioning graph and the parameterizations. Order, pruning, caching, indexing and the computation scheme belong to an engine.
2. **Several engines, one model.** An engine declares which models and counters it admits. A restriction that helps one algorithm becomes that engine's admissibility condition, not a rule of the model.
3. **The event holds structure and semantics; the handler holds the numbers.** Error models and covariate models are ordinary nodes with handlers.
4. **Stated attributes, derived views.** Events carry their targets, side and library explicitly. Per-segment views are derived and checked. The compact per-segment notation is accepted as shorthand.
5. **Place every computation by tier.** A computation's tier is fixed by what it depends on: 0a structure, 0b parameters, 1a read, 1b read × parameters, 2 scenario. It reads nothing from a higher tier, and it is cached for exactly its tier's lifetime (review §13).
6. **One file grammar.** Every polymorphic object is tagged by `type`, and every object that can be referenced has a `nickname`. A model is one sectioned document (`structure`, `library`, `parameters`, `provenance`, `accumulators`); each section is inline or referenced by hash (review §14, §15).
7. **Keep what works.** The legacy DFS stays as one engine, with its bitwise goldens. The five Contexts stay, and the Context refactoring is completed.

## 2. Target architecture

```
┌───────────────────────────────── Model (RecombinationModel) ─────────────────────────────────┐
│                                                                                              │
│  SeqTypeRegistry        SegmentLayout                 ConditioningGraph (today: Topology)    │
│  name <-> SeqTypeId     DAG of segments, guards on    events + exogenous covariate nodes;    │
│                         structural events; paths()    parents, topological order, acyclic    │
│                                                                                              │
│  Events (minimal): uid, nickname, type, domain, targets (seq_type, side), library,           │
│                    footprint (roles, deltas, reads), semantics (candidates/evaluate/apply/   │
│                    prepare); built from JSON                                                 │
│  Probe event:      the error model (sequence-valued, context window)                         │
│                                                                                              │
│  Parameterization per node (declared): Categorical | Markov | Bernoulli | Logistic |         │
│                    LinearModel | ...  -> distribution / maximize / sample / bound            │
│  Parameter store:  one entry per node (a tensor for tensor handlers, a vector otherwise)     │
│  Derived:          SegmentIndex (creators and modifiers per segment and end), validation     │
└───────────────────────────────────────────────┬──────────────────────────────────────────────┘
                                                │ const
                   ┌────────────────────────────┴───────────────────────────┐
                   ▼                                                        ▼
┌──────────────── EM driver ─────────────────────┐  ┌────────────── Sampling driver ──────────────┐
│ iterations, data and read policies,            │  │ RNG streams per chunk, error application,   │
│ chunked deterministic reduction (#68),         │  │ output writers                              │
│ per-read normalisation, counter lifecycle,     │  │                                             │
│ M-step (maximize on every node), output        │  │                                             │
└───────────────┬────────────────────────────────┘  └───────────────┬─────────────────────────────┘
                │ check_admissible / prepare_iteration /            │ sample(rng) -> ScenarioRecord
                │ make_workspace / process_read [/ prepare_read]    │ + Event::apply
                ▼                                                   ▼
┌──────────────── Inference engines ─────────────┐  ┌────────────── Sampling engines ─────────────┐
│ LegacyDfsEngine     (DFS-A, bitwise)           │  │ Exact topological sampler (indexed CDF)     │
│ ChainDPEngine       (SEGMENT §5)               │  │ Two-loop sampler          (SEGMENT §4)      │
│ DpBoundedDfsEngine  (semi-DP, all counters)    │  │ Fast approximate sampler  (FastGenerator)   │
│ Posterior sampling / stochastic EM (FFBS)      │  │                                             │
│ (DFS-B: orchestrator-driven, to be measured)   │  │                                             │
└───────────────┬────────────────────────────────┘  └─────────────────────────────────────────────┘
                │ reads and writes
                ▼
┌──────────────── Contexts (all five kept) ─────────────────────────────────────────────────────┐
│ QuerySequenceContext (const, + per-read caches)  ModelContext (const, handler rows)           │
│ ScenarioContext (+ realization record, per-scenario adjacency; written one layer at a time)   │
│ ExplorationContext (per engine; today the DFS's)  AccumulationContext (+ AccumulatorSet)      │
└───────────────────────────────────────────────────────────────────────────────────────────────┘
```

How the diagram maps onto the computation tiers (review §13):

| Tier | Where it lives in the diagram |
|---|---|
| 0a | model structure and the engine plan (`check_admissible`) |
| 0b | handlers, materialised in `prepare_iteration` |
| 1a | `QuerySequenceContext` + `ReadCache` (`prepare_read`), held by the EM driver |
| 1b | the engine workspace, in the prologue of `process_read` |
| 2 | `ScenarioContext`, in the inner loop of `process_read` |

The model document of review §15 serialises the Model box. `structure` and `library` cover 0a; `parameters` covers 0b's inputs; `accumulators` covers what the EM driver merges.

## 3. Concepts

| Term | Definition | Review § |
|---|---|---|
| Event | A random variable of the model. It holds identity, domain, targets, footprint and semantics. It has no parents, no parameters and no scenario state. | 7.1 |
| Domain | The indexed realizations `0..n-1` with their payloads. Scalar, or sequence-valued over the positions of a segment. Distinct from the parameter shape. | 7.1 |
| Target | A `(seq_type, side?)` pair that the event acts on. An event may have several. | 2, 7.1 |
| Library | The template set and alignment set that a templated event draws from (today's `gene_class`). D1 and D2 share one. | 2.1 |
| Footprint | What an event creates, modifies and reads, and by how much. These are the answered capabilities (A0, `affects_*`, `SeqContextDependency`). | 7.1 |
| Segment | One piece of the receptor: a slot in the layout, identified by `SeqTypeId`. | 2 |
| Layout | A DAG of segments with existence guards. A path is one source → sink chain. A linear layout is a special case. | 4 |
| Structural event | An ordinary event whose realizations drive guards, for example the number of D segments or the locus. | 4.2 |
| Conditioning graph | A DAG over events and exogenous covariates. It holds the statistical dependencies. | 1 |
| Existence dependency | An event on a guarded segment exists only when the guard holds. Distinct from conditioning. | 4.2 |
| Parameterization / handler | The declared family of P(event \| context): `distribution`, `maximize`, `sample`, `bound`. | 6 |
| Accumulator | Sufficient statistics owned by the handler, with `empty_copy` and an ordered `+=`. | 6 |
| Engine | One computation scheme. It declares what it admits and implements the driver calls. | 5 |
| Driver | The generic workflow around engines: EM, or sampling. | 5.2 |
| Computation tier | What a computation depends on, which fixes its cadence and cache lifetime: **0a** structure (once per model), **0b** parameters (once per EM iteration), **1a** read (once per read, valid across iterations), **1b** read × parameters (once per read per iteration), **2** scenario (per branch). Refines G10's tiers 0/1/2 | 13 |
| Check level | C1 capability at load, C2 leaf assert, C3 hand-off in tests: the checks of the segment-completeness invariant. Called "tiers" in older docs; renamed here to avoid the collision | 13.1 |
| `prepare(read)` / ReadCache | Tier-1a work, i.e. per read and independent of scenario and parameters. This is also SEGMENT's per-segment table. | 8.3, 13 |
| Discriminator | The `type` field of a polymorphic object. It selects the constructor from a registry and is not an identity | 14 |
| Nickname / uid | The nickname is the persistent identity and the only link between objects. The uid is a runtime dense index assigned in byte-wise nickname order, so it never depends on the file's layout. Engines break ties by nickname | 14.2 |
| Realization label / index | The label (gene name, value, nucleotide) is the identity. The index is the tensor coordinate: it exists only where parameter values exist, and is decided by the structure (explicit `index`, otherwise the next free index in label order). Parameters are matched to realizations by label | 14.4 |
| Authoring / canonical form | Two forms of one schema: lenient for the modeller (shorthand, optional indices, comments), canonical as written by IGoR (sorted keys, explicit defaults, round-trip numbers). `fmt` converts between them; reloads are bitwise identical | 14.5 |
| Model document | One sectioned JSON: `structure`, `library`, `parameters`, `provenance`, `accumulators`. Each section is inline or `{"$ref", "sha256"}` | 15 |
| Fragment / overlay | A partial export (some sections, some nodes). An overlay applies it onto a base model, after a per-node domain check | 15.5 |
| EventOutcome | What one realization produces: one `ScenarioContext` layer, the bound facts, and a status. | 9.2 |
| Scenario record | `{uid → realization}` in `ScenarioContext`. Read by counters and error models. | 7.2, 9.1 |

## 4. Changes to SEGMENT (`SEGMENT_DECOMPOSITION_PROPOSAL.md`)

| SEGMENT § | SEGMENT says | Synthesis | Status | Review § |
|---|---|---|---|---|
| 1 | Inference *is* a sum over the layout with per-segment tables; the interleaved exploration is removed | **This is one engine (`ChainDPEngine`) among several; the legacy DFS stays** | *Changed* | 5.1 |
| 2.1 | Roles Create / Modify / Probe | Kept, as **answered capabilities per target; an event may have several. New: structural events** | *Changed* | 4.2, 7.1 |
| 2.1 | No gene class, no side, no generated name; identity = position in a segment | **Targets, side and library kept, stated on the event**; generated name dropped; identity = UID | *Changed* | 1, 2 |
| 2.1 | Parameterization = own axes (`inherent_shape()`) | Domain shape **distinct from the parameter shape**; the parameterization is declared | *Changed* | 6, 7.1 |
| 2.2 | Segment = an ordered composite of events | **Derived `SegmentIndex`**: a view, not storage | *Changed* | 2.3 |
| 2.2 | Gene/junction kinds; one Create and at most one Modify per end; junction content depends only on its neighbours | **Admissibility conditions of an engine, not rules of the model** | *Changed* | 3 |
| 2.2 | Boundary interface: offset and end nucleotides | Kept, as the DP engine's interface; its width is a parameter | *Kept* | 5.1 |
| 2.3 | Layout = ordered list; adjacency = position ± 1 | **A DAG with guards; adjacency per scenario** | *Changed* | 4 |
| 2.3 | — | **Lightweight DAG in the model now; an experimental engine may admit only the two-path tandem-D layout** | *New* | 4.3 |
| 2.4 | Conditioning graph = model data; parents, topological order, acyclicity | Kept | *Kept* | 1 |
| 2.4 | — | **Exogenous covariate nodes; existence dependency** | *New* | 4.2, 6 |
| 2.5 | Two handler kinds, chosen by the rank of the own axes | **Declared per node, with a default per event type; rank is a check** | *Changed* | 6 |
| 2.5 | Error model outside events (until §8 step 7) | **Probe event with Categorical, Logistic or Bernoulli handlers; one accumulator shape** | *Changed* | 6 |
| 3 | `events` is an object keyed by nickname, with `type` as a property | Kept: `type` is a discriminator that selects the constructor; the nickname is the key. **The order of events is not content: no engine reads it, the uid follows nickname order, and ties are broken by nickname (legacy ties follow file order, Model_Parms.cpp:474). Duplicate keys are rejected at parse** | *Changed* | 14.1–14.3 |
| 3 | Realizations carry an explicit `index` | **Label = identity; index = tensor coordinate, present only where parameters exist; may be omitted by hand and assigned by `fmt`; parameters matched by label; explicit policy for missing or extra realizations (renormalisation recorded in provenance)** | *Changed* | 14.4 |
| 3 | — | **Authoring and canonical forms, `fmt` / `validate` / `diff`, round-trip contracts, a fingerprint over the canonical form** | *New* | 14.5 |
| 3 | Marginals in the same document or a sibling file, as one array per event | **A `parameters` section of a sectioned document, inline or referenced by hash. Keyed by nickname, with labelled `axes` and per-axis domain fingerprints; round-trip precision. The error-model value moves there from the structure** | *Changed* | 15 |
| 3 | — | **Partial exports (sections × nodes) and overlays; `library` and `accumulators` as their own sections** | *New* | 15.5 |
| 3 | JSON v3 with `layout`, `segments`, `conditioning`, `events`, `error_model` | JSON v3 kept. **Events carry targets, side and library; `layout` is a DAG; a parameterization per node; the error model is an event. `segments` becomes accepted shorthand** | *Changed* | 2.3, 4, 7.1 |
| 3 | v1 and v2 read by an importer | Kept. **The importer emits JSON; one generic JSON constructor per event type** | *Changed* | 7.1 |
| 4 | Generation = topological outer loop + inner layout loop | **One sampling engine, admissible when no scalar event depends on constructed content** | *Changed* | 5.1 |
| 4 | Indexed CDF sampler; goldens regenerated once | Kept; lands first | *Kept* | 1, 12 |
| 4 | `FastGenerator` and `draw_random_realization()` retired | `draw_random_realization` retired. **`FastGenerator` kept as an approximate sampling engine** | *Changed* | 5.1 |
| 5.1–5.3 | Factorisation, segment tables, forward-backward assembly | Kept, as `ChainDPEngine`. **Segment tables = the `prepare(read)` caches of events** | *Changed* | 8.3 |
| 5.3 | Three semirings; counters through extended state | Kept, **plus counter classes, DP-bounded DFS and FFBS, so that every counter is supported** | *Changed* | 5.3 |
| 5.3 | Pruning heuristics disappear | They disappear in the DP engines. **They stay in `LegacyDfsEngine`; in `DpBoundedDfsEngine` the DP messages become the exact bounds** | *Changed* | 5.3 |
| 5.4 | EM from posterior marginals; `maximizeLikelihood()`; #68 reduction | Kept, **inside a generic EM driver extracted from `GenModel::infer_model`** | *Changed* | 5.2 |
| 5.5 | Boundary windows for N-mer error models, later | Kept. **Importance sampling over an admissible approximation is another option** | *Changed* | 5.3 |
| 5.6 | Events hold no per-read state; no per-thread copies | Kept as the goal. **Reached by completing the Contexts (scenario record) and stateless events** | *Changed* | 9, 10.3 |
| 6 | `iterate()` bodies, exploration context, spans, offsets, `SafetyMatrix`, `LayeredArray`, `Index_map` removed | **Not removed: they become `LegacyDfsEngine`'s kernels and workspace.** Only the orchestration leaves `Rec_Event` | *Changed* | 7.2, 10.3 |
| 6 | `Topology` shrinks to a parents table | Kept, plus covariates. Edge inversion and ancestors are not needed at runtime | *Kept* | — |
| 7 | Bitwise gate dropped; three checks replace it | **Bitwise gate kept for `LegacyDfsEngine`.** The three checks apply to the new engines | *Changed* | 1, 5.1 |
| 8 | Order: DP prototype first | **Replaced by the merged order (§8 below)**; the DP prototype comes after the engine split | *Changed* | 12 |
| 9 | Priority: maybe kept as data | **An optional engine hint** | *Changed* | 7.1 |
| 9 | Parameter tying for D1/D2 | **A shared handler (one parameter entry for two nodes); the M-step sums their accumulators** | *Changed* | 6 |
| 9 | Per-scenario outputs: to decide | **Counter classes** decide which engine serves which output | *Changed* | 5.3 |
| 9 | Negative deletions: extended content in the table | Kept, **computed in `prepare(read)`** | *Changed* | 8.3 |

## 5. Changes to MODEL (`MODEL_TOPOLOGY_REFACTORING.md`)

| MODEL item | MODEL today | Synthesis | Status | Review § |
|---|---|---|---|---|
| 1. Single source of truth | `RecombinationModel` owns `Topology` and one `Tensor<T>` per event; handlers borrow | Kept. **The parameter store holds one entry per node and per parameterization (a tensor, or a vector for Logistic/LinearModel). `RecombinationModel` also owns the registry, the `SegmentLayout` and the error-model node** | *Changed* | 1, 4, 6 |
| 2. Ownership | `InferenceEngine` holds a **mutable** model; `SamplingEngine` a const one | **All engines read a const model. The EM driver applies `maximize` between iterations** | *Changed* | 5.2 |
| 3. Borrowed references | Handlers borrow their tensors, never copy | Kept | *Kept* | — |
| 4. Two handler hierarchies | `InferenceHandler` (weights + accumulator, `maximizeLikelihood`) and `SamplingHandler` (weights + CDFs, `sample`), each Categorical or Markov | **Split by role, not by engine: a parameterization (`distribution`, `maximize`, `sample`, `bound`), separate accumulators (`empty_copy`, ordered `+=`), and sampler caches (CDFs). New families: Bernoulli, Logistic, LinearModel** | *Changed* | 6 |
| 4. Accumulator | Owned by the inference handler | **Separate object, owned by an `AccumulatorSet`, chunked for #68** | *Changed* | 6 |
| 5. Factories | Self-registering `Registrar` keyed by `Event_type` | Self-registration kept. **Keyed by the declared parameterization; a default per event type; rank is a check** | *Changed* | 6 |
| File loading | `read_topology()` (v1 only) + `read_parameters()`; `recombination_model_from_files(parms, marginals)` | **`load_model(document, overlays)` over the sectioned JSON. Two-phase load: schema validation from the registries, then construction with uids in nickname order and resolution of references by nickname; parameters matched by label. `export_model(model, spec)` for full or partial export. v1/v2 become importers** | *Changed* | 14, 15 |
| Parameter precision | Legacy marginals are written with 6 significant digits from `long double` and renormalised at load | **Written with the stored type's `max_digits10`, so a reload is bitwise identical (goldens over marginal files move once)** | *Changed* | 15.4 |
| — | No tier discipline | **Each engine call and each event or handler method is tied to a tier. 0a work leaves the EM iteration loop; the 1a cache lives with the driver** | *New* | 13.4 |
| 5. Tensor shape | README: `[own_dims, parent_dims…]`; normalise axis 0 | **`[parents…, own…]`, as the code already does (commit `ca97599`); the README and handler comments are stale** | *Changed* | 11 |
| 5. `EventFactory` | Default-constructs, then calls base setters (cannot reach subclass members) | **One generic JSON constructor and `to_json` per event type** | *Changed* | 7.1 |
| 6. Navigator | Index-based views: parents, children, topological order | Kept. **Topological order is one engine's order, not *the* order** | *Kept* | 5 |
| 7. LegacyBridge | `Model_Parms ⇄ Topology`, `Model_marginals → RecombinationModel` | **Becomes the v1/v2 → JSON importer**; export to `Model_Parms` only while the legacy runtime needs it | *Changed* | 7.1 |
| 8. Scenario value objects | `SampledEvent{uid, indices}`, `SampledScenario`, for sampling only | **Promoted to the scenario record in `ScenarioContext`, used by inference, counters and error models** (replacing the `friend` raw pointers) | *Changed* | 7.2, 9.1 |
| `InferenceEngine` | A container: `run(eStep)` = reset → callback → update; no E-step of its own | **An engine strategy: `check_admissible`, `prepare_iteration`, `make_workspace`, `process_read` [, `prepare_read`], under the generic EM driver. Several implementations** | *Changed* | 5.1, 5.2 |
| `SamplingEngine::run` | Topological order; one index per parent; `sampleSequence` never called | **Calls `sampleSequence` for sequence-valued events and the event's `apply()`; parents may contribute several indices** | *Changed* | 5.1 |
| `MarkovSamplingHandler` | Row index differs between `precomputeCDF` and `sample` | **Fixed** | *Changed* | 11 |
| Sampling `long double` | Not registered (double only) | **To register, or document as unsupported** | *Changed* | — |
| `Topology` | Keyed by nickname; edge inversion and ancestors | **Renamed `ConditioningGraph`; keyed by UID; covariate nodes**; inversion and ancestors only off the runtime path | *Changed* | 1 |
| — | No layout | **`SegmentLayout` (a DAG with guards) and `SeqTypeRegistry` in `RecombinationModel`** | *New* | 4 |
| — | No error model | **The error model as a Probe node with a handler** | *New* | 6 |
| — | No driver | **EM driver and sampling driver** | *New* | 5.2 |
| Data flow: inference | `read_topology` + `read_parameters` → `RecombinationModel` → `InferenceEngine.run(E-step)` | **JSON → `RecombinationModel` (const) → EM driver → engine.`process_read` per read → `AccumulatorSet` → `maximize`** | *Changed* | 5.2 |
| Data flow: generation | `RecombinationModel` → `SamplingEngine.run(rng)` → `SampledScenario` | **JSON → `RecombinationModel` → sampling driver → engine.`sample` → `ScenarioRecord` → `Event::apply` → sequence + errors** | *Changed* | 5.2 |
| Build | `igor::Model` links `igor::Core` and `igor::Math` | Kept. **Events move toward `Model`; the legacy runtime stays in `Core` as `LegacyDfsEngine`** | *Changed* | 10.3 |

## 6. Kept from both

- The conditioning graph as model data, with a topological order and an acyclicity check.
- A UID per event; the nickname as a handle.
- The indexed CDF sampler.
- Borrowed parameter references and a single store.
- Self-registering factories.
- The `[parents…, own…]` tensor layout.
- The deterministic chunked reduction (#68).
- The JSON v3 file, with v1 and v2 as import formats.
- Validation check 1: any new engine against the exhaustive legacy run at threshold 0.
- The boundary interface of segments for the DP engine.
- The (+,×), (max,×) and extended-state semirings.

## 7. Open decisions

| # | Decision | Options | Review § |
|---|---|---|---|
| 1 | When the layout DAG comes in | Recommended: the DAG and guards in the model now, with a narrow-admission experimental engine. Alternative: hard-code it in the experiment first | 4.3 |
| 2 | Which DFS design | DFS-A (continuation, zero-copy, bitwise) is certain. DFS-B (orchestrator-driven `evaluate`) depends on the measurement | 10.3 |
| 3 | JSON library | No JSON in the C++ tree today; choose a dependency (e.g. header-only) | 7.1 |
| 4 | Which counters to keep | Pgen per sequence and best scenarios with the legacy tie rules need enumeration or sampling | 5.3 |
| 5 | Priority | Dropped from the model, or kept as an optional engine hint | 7.1 |
| 6 | Viterbi mode | Fix the overwrite without clearing, or retire the flag in favour of k-best | 11 |
| 8 | Numeric encoding and storage type of parameters | Decimal at `max_digits10` (recommended); hexfloat (exact, but 80-bit is not portable); `.npy` sidecar for large tensors. And the storage type: keep `long double` (needs a `long double` JSON float type, or string values) or move to `double` (goldens move once) | 14.5, 15.4 |
| 9 | Persisting the 1a cache | Recommended: recompute per read in the thread workspace (core 1a ≈ 5 kB per TRB read, a few per cent of an E-step); persist only under a memory budget, if measured worthwhile; never persist dense DP tables (≈ 215 kB per TRB read, ≈ 2.3 MB per IGH read) or tier 1b | 13.5 |
| 11 | Compacting `Alignment_data` | ≈ 30 kB per TRB read today, held for every read; ≈ 7 kB after compaction (library index, `uint16`, no empty vectors). The largest memory saving available | 13.5 |
| 10 | Strictness of the structure fingerprint | A warning (only per-node checks bind) or an error, for parameter documents and overlays | 15.4 |
| 7 | `Single_error_rate` M-step | Keep the mean of per-sequence frequencies (bitwise), or move to the pooled MLE (goldens move) | 6 |

## 8. Merged order of work

Each step is tagged with the design it changes: **S** = SEGMENT, **M** = MODEL, **C** = Core/legacy runtime.

| # | Step | Tags | Goldens |
|---|---|---|---|
| 1 | Indexed CDF sampler | S, M | generation goldens move once |
| 2 | UID replaces the generated name; fixes the tandem-D collision | S, M, C | unchanged |
| 3 | Scenario record in `ScenarioContext`; `friend` raw pointers removed | M, C | unchanged |
| 4 | Generic JSON event constructor and `to_json`; type dispatch, label-matched parameters and the index rule (§14); v1/v2 importers; `EventFactory` and `LegacyBridge` rewired | S, M | unchanged |
| 4b | Sectioned model document: `parameters` with axes and fingerprints, `library`, `provenance`; partial export and overlay; round-trip precision (§15) | S, M | marginal-file goldens move once |
| 5 | Registry/layout split; `SegmentLayout` DAG with guards; adjacency per scenario | S, M, C | unchanged |
| 6 | EM driver extracted from `GenModel::infer_model`; `LegacyDfsEngine` (DFS-A: continuation, bound fold and layers inside the engine) | M, C | **bitwise** |
| 7 | Pure semantics per type, annotated by tier (§13); `prepare(read)` (1a); 0a work moved out of the iteration loop; `evaluate` and `EventOutcome`; DFS-B prototype and measurement | C | bitwise for DFS-A |
| 8 | Handler contract: parameterizations, `AccumulatorSet`, declared selection; error models as a Probe node (code moved, not rewritten) | S, M | **bitwise** |
| 9 | `ChainDPEngine` on the `prepare` caches, validated by check 1; then `DpBoundedDfsEngine` and FFBS | S, M | tolerance |
