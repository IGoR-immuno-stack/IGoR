# Layer refactoring proposal

Status: proposal, 2026-10-02, revised the same day with the Legacy policy (§4) and the
three-step stage 1 (§6). Steps 1a and 1b are merged; the first pass of 1c is recorded in §6.
Based on `feature/model-on-tandemD` @ `5ac1972`.
Related: `ARCHITECTURE_SYNTHESIS.md` (Q. Marcou, 2026-10-01, not yet in the repo),
`docs/SEGMENT_DECOMPOSITION_PROPOSAL.md` (superseded), `doc/CORE_TO_MODEL_MIGRATION_PLAN.md`
(obsolete: it assumed Model would replace Core).

## 1. Why the current layers no longer fit

`Model` and `Math` were cut so that the tensor-based model could be built next to the legacy
one without touching it. That constraint is gone: there will be one IGoR, and the synthesis
settles the target (events keep their attributes, the layout is a DAG, handlers are declared
per node, the legacy DFS is one inference engine among several).

What is left is a layering problem, not a coexistence problem:

| Layer | Lines | Depends on | Content today |
|---|---:|---|---|
| Core | 29 550 | GSL, nlohmann_json, OpenMP | everything legacy: enums, containers, aligner, events, error rates, EM driver, counters, generator, model I/O, DFS contexts and geometry |
| Math | 2 310 | nothing | Tensor, Linalg, HybridBuffer, mdspan polyfill |
| Model | 2 838 | Core, Math | Topology, RecombinationModel, inference and sampling handlers and engines, EventFactory, LegacyBridge |
| Streaming | 3 734 | Core | Parquet and AIRR readers and writers |

Core is not a "common concepts" layer: it is the whole legacy program. The symptoms:

- Streaming links the whole of Core for two things: `Alignment_data` and the `Gene_class` enum.
- `igor-convert` links Core although it never aligns, infers or generates.
- Model's `LegacyBridge` is the only reason the tensor engines know `Model_Parms` and
  `Model_marginals`, and the name says "two worlds" when the target is one.
- `Utils.h` holds `Matrix<T>`, the hashes, the DFS scratch-map typedefs, string helpers and
  `gethostname()`. Quentin's iterate plan already says it must keep shrinking.
- Almost every ex-Core class sits in the global namespace, in `Snake_case` (`Rec_Event`,
  `Model_Parms`, `Best_scenarios_counter`), while Model is `igor::model` and PascalCase.

## 2. Target layers

The rule stays the one you always use: one `Core` with the shared vocabulary, thematic layers
around it, a DAG of dependencies, no layer depends on a layer above it.

```
                 app/igor   app/igor-convert   app/igor-demo
                     │             │                │
   ┌─────────────────┼─────────────┼────────────────┤
   ▼                 ▼             ▼                ▼
Inference        Generation     Streaming       Alignment
   │  │              │  │           │                │
   │  └──────┬───────┘  │           │                │
   ▼         ▼          ▼           │                │
 Model ──► Math         │           │                │
   │                    │           │                │
   └────────────────────┴───────────┴────────────────┴──► Core
```

| Layer | Namespace | Role | May depend on | External deps |
|---|---|---|---|---|
| **Core** | `igor::core` | the vocabulary: event and segment enums, `SeqTypeRegistry`, `SegmentSpan`, `Int_Str`, genetic code, `Alignment_data`, the layered containers, string helpers, export and config headers | nothing | none |
| **Math** | `igor::math` | numerics with no IGoR semantics: `Tensor`, `Linalg`, `HybridBuffer`, mdspan compat, `Matrix<T>`, categorical samplers (`FastSampling`) | nothing | mdspan |
| **Alignment** | `igor::alignment` | algorithms on sequences: Smith-Waterman aligner, CDR3 and feature extraction, `JournaledQuery` and AA motifs | Core, Math | OpenMP |
| **Model** | `igor::model` | the probabilistic model and its persistence: `Rec_Event` and its four kinds, error-rate models, graph (`Model_Parms`/`Topology`), parameters (`Model_marginals`/`RecombinationModel<T>`), `EventFactory`, txt v1/v2 and JSON readers and writers | Core, Math | nlohmann_json, GSL |
| **Inference** | `igor::inference` | computing with a model on observed sequences: EM driver, the inference engines (legacy DFS, tensor handlers, later the per-segment DP), the observers (counters: Pgen, best scenarios, coverage, errors) | Model, Math, Core | OpenMP |
| **Generation** | `igor::generation` | sampling from a model: `FastGenerator`, `SamplingEngine` and its handlers, error injection | Model, Math, Core | none |
| **Streaming** | `igor::streaming` | sequence and alignment data I/O: Parquet, AIRR | Core | Arrow, Parquet, sparrow |

Two decisions embedded in the table:

- **Inference and Generation are two layers, not one `Engines` layer.** They share Model and
  the scenario value types, not code. `igor-convert` and a future Pgen service link only what
  they need. The synthesis describes generation as "the same structure run forward", which is a
  symmetry of design, not a shared implementation.
- **Model I/O stays in Model, not in Streaming.** Streaming is about sequence data; the model
  file formats (txt v1/v2, JSON) are the model's own persistence. Keeping them next to the
  events avoids a Streaming → Model dependency and keeps the JSON work where it is.

## 3. Conventions

They are the ones Model already follows; the other layers align on them.

| Topic | Rule | Example |
|---|---|---|
| Directory | `src/igor/<Layer>/` for the layer, `src/igor/<Layer>/Legacy/` for legacy files | `src/igor/Model/Legacy/Rec_Event.h` |
| Namespace | `igor::<layer>` for the layer, `igor::<layer>::legacy` for `Legacy/` | `igor::model::legacy::Rec_Event` |
| Sub-namespaces | lower_snake, for factories and detail only | `igor::model::event_factory`, `igor::math::detail` |
| Include path | full path from `src`, angle brackets, including `Legacy/` | `#include <igor/Model/Legacy/Rec_Event.h>` |
| Types | PascalCase, one type per file, file named after the type | `RecombinationModel.h`, `SamplingEngine.tpp` |
| Methods, free functions | camelCase | `topologicalOrder()`, `addEdge()` |
| Members | `m_` prefix | `m_weights` |
| Template parameter | `T` for the scalar type | `InferenceHandler<T>` |
| Export | one macro per layer, from `igor/<Layer>/Export.h`; `Legacy/` uses the same macro | `MODEL_EXPORT` |
| CMake | one target per layer; `Legacy/` is a subdirectory of that target, never a target of its own | `igor::Model` |
| Tests | `tst/igor/<Layer>/` and `tst/igor/<Layer>/Legacy/` mirror the sources | `tst/igor/Model/Legacy/test_model_json.cpp` |

`Legacy/` is not a CMake target because legacy and new code of one layer reference each other
in both directions during the transition (`Topology` holds `shared_ptr<Rec_Event>`,
`Model_Parms::finalize()` will one day call `Topology`). A separate target would force a
dependency direction that does not exist yet.

The include path carries `Legacy/` on purpose: the remaining legacy debt of any consumer is
one `grep -r "Legacy/"`, and the end state is checkable in CI: no file outside a `Legacy/`
directory includes a `Legacy/` header.

## 4. Legacy policy

Every ex-Core file lands first in the `Legacy/` directory of its layer, inside
`igor::<layer>::legacy`, with no exception and no judgement at that point (step 1a). Files
already written to the conventions (Math, the tk Model files, Streaming) go straight to the
layer directory. From there, each legacy file has one of three fates, decided file by file in
§5 and executed case by case in step 1c:

| Tag | Meaning | Mechanism | Ends when |
|---|---|---|---|
| **keep** | the concept is obsolete; the file lives only for its remaining consumers | nothing moves | the last consumer outside `Legacy/` is gone; file and its tests are deleted |
| **promote** | the concept lives on and the design is not changing; only name, namespace and style are wrong | the file moves out of `Legacy/` into the layer directory and namespace; a rename is done in place, with a `using` alias left in `legacy` for old consumers | step 1c, one file or cluster per PR |
| **replace** | the concept lives on but gets a new design | the file is **copied** into the layer directory under its new name, reworked there; the legacy original stays untouched for its consumers | the copy has an owner and a target date; the original is deleted when its last consumer switches |

Three rules that follow from the tags:

1. **Vocabulary is promoted, never replaced.** Enums, typedefs, `SeqTypeRegistry`,
   `SegmentSpan`, `Int_Str`, `Alignment_data` are shared by legacy and new code. A copy would
   give two incompatible types and conversions at every boundary. Renames of these are done in
   place with a `legacy` alias:
   ```cpp
   namespace igor::core { class IntStr : public std::vector<int> { /* ... */ }; }
   namespace igor::core::legacy { using Int_Str = igor::core::IntStr; }
   ```
2. **A replace copy is a ticket, not a parking spot.** Two copies of `iterate()` or of the
   hypermutation models alive for months is the coexistence we are leaving. A copy is made
   when someone starts the rework, and the original is deleted when that rework lands. No
   copy of a file Quentin is currently rewriting (the four events, the contexts) before his
   R-queue and B5/B11 are closed: those are tagged replace but **deferred to stage 2**.
3. **Quentin's 2026 headers are the first promotions.** `LayeredArray`,
   `DynamicSequenceMap`, `SeqOffsetsMap`, `SeqTypeRegistry`, `SegmentSpan`, `SpanProfile`,
   `JunctionGeometry`, `SafetyMatrix`, `UnfilledSegmentLengths`, `BoundTightness`,
   `JournaledQuery` are the new design, written this year, mostly already PascalCase. They
   go through `Legacy/` like everything else in step 1a, but they head the step 1c queue,
   so that their author does not read `Legacy/` on his own work for longer than one PR.

## 5. Where each file goes

Lines are the current `.h` + `.cpp` + `.tpp` counts. "From" is the current directory, the
layer is where the file lands in step 1a (under `Legacy/` for every ex-Core file). The tag is
the step 1c decision; "Target name" is the name after promote or replace, blank means
unchanged.

### Core (≈ 3 200 lines)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Core | `CoreEnums.h`, `StdTypedefs.h`, `Typedef.h` | promote | merged into `Core/Types.h`; enum values renamed in place later |
| Core | `IntStr.*` | promote | `IntStr`, alias `legacy::Int_Str` |
| Core | `GeneticCode.*` | promote | |
| Core | `SeqTypeRegistry.h`, `SegmentSpan.h` | promote | |
| Core | `LayeredArray.h`, `DynamicSequenceMap.h`, `SeqOffsetsMap.h` | promote | `SeqOffsetsMap`, alias `legacy::Seq_offsets_map` |
| Core | `Alignment_data` (extracted from `Aligner.h` in step 1c) | promote | `Core/AlignmentData.h`, `AlignmentData`, alias `legacy::Alignment_data`; until then it lives in `Alignment/Legacy/Aligner.h` |
| Core | `Utils.*` | keep | `Core/Legacy/Utils.h`; pieces leave on demand (`Matrix` → Math now, string helpers → `Core/Strings.h` when a new consumer needs them) |
| Core | `Config.h.in`, export header | promote | `Core/Export.h` |

### Math (≈ 3 300 lines)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Math | `Tensor.*`, `TensorCreation.*`, `Linalg.*`, `HybridBuffer.*`, `MdspanCompat.h` | new | unchanged |
| Core | `Matrix<T>` (from `Utils.h`) | keep | `Math/Legacy/Matrix.h`; replaced by `Tensor` in its five consumers over time |
| Core | `FastSampling.*` | promote | `CategoricalSampler`, `ConditionalSampler` into `igor::math`; Math stops being header-only or the file becomes a `.tpp` |

### Alignment (≈ 4 300 lines)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Core | `Aligner.*`, `AlignerInternal.h` | promote | 3 000 lines of Smith-Waterman nobody plans to rewrite; namespace and method names only. `AlignerInternal.h` and its `swalign` testing export stay: they exist for the three whitebox DP tests of `test_aligner.cpp` (decided 2026-10-03) |
| Core | `ExtractFeatures.*`, `CDR3SeqData.*` | promote | |
| Core | `JournaledQuery.*` | promote | Nicolas's AA_PGEN work, already in style |

### Model (≈ 18 500 lines in stage 1, ≈ 8 000 after stage 2)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Core | `Rec_Event.*`, `Genechoice.*`, `Deletion.*`, `Insertion.*`, `Dinuclmarkov.*` | replace, **stage 2** | `Event`, `GeneChoice`, `Deletion`, `Insertion`, `DinuclMarkov`: definition only, `iterate()` goes to Inference |
| Core | `EventUtils.*`, `gene_to_seqtype_migr.*` | keep | migration shims; die with the legacy formats |
| Core | `JsonDetail.h`, `ModelJson.*` | promote | ours, recent; keyed on `Model_Parms` until `Topology` replaces it |
| Core | `Errorrate.*`, `Singleerrorrate.*`, `Hypermutationglobalerrorrate.*`, `HypermutationfullNmererrorrate.*` | replace | `ErrorModel<T>`, `BernoulliErrorModel<T>`, `ContextErrorModel<T>`; not before stage 2 (they take the contexts) |
| Core | `Model_Parms.*` | keep | replaced by `Topology`, which exists; deleted when the DFS engine reads `Topology` |
| Core | `Model_marginals.*` | keep | replaced by `RecombinationModel<T>`, which exists; deleted when the DFS engine reads tensors (the open marginal-storage contract) |
| Model | `EventFactory.*`, `Topology.*`, `RecombinationModel.*`, `Scenario.h` | new | unchanged |
| Model | `LegacyBridge.*` | new, rename | `MarginalsBridge`: it converts between two parameter layouts of one model |
| Core | `Counter.h/.cpp` (base class) | replace, stage 2 | `Observer<T>`; in `Model/Legacy/` until then because `iterate_wrap_up()` calls it |
| Core | `QuerySequenceContext.h`, `ModelContext.h`, `ScenarioContext.h`, `ExplorationContext.h`, `AccumulationContext.h`, `Scenario.h` | replace, stage 2 | DFS engine state; `Model/Legacy/` only because `Rec_Event::iterate()` takes them |
| Core | `JunctionGeometry.h`, `SpanProfile.h`, `SafetyMatrix.h`, `UnfilledSegmentLengths.h`, `BoundTightness.h` | promote | head of the 1c queue, into `Model/` and `igor::model`; relocated to `Inference/` with `iterate()` in stage 2 |

### Inference (≈ 5 000 lines in stage 1, ≈ 15 000 after stage 2)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Core | `GenModel.*` | replace | `EmDriver`; `load_genmodel()`/`readtxt()` go to Model's persistence |
| Core | `Pgencounter.*`, `Bestscenarioscounter.*`, `Coverageerrcounter.*`, `Errorscounter.*` | replace | `PgenObserver`, `BestScenariosObserver`, `CoverageObserver`, `ErrorStatsObserver` |
| Model | `InferenceEngine.*`, `InferenceHandler.*`, `CategoricalInferenceHandler.*`, `MarkovInferenceHandler.*`, `InferenceHandlerFactory.*`, `Navigator.h` | new | namespace becomes `igor::inference` |
| stage 2 | `iterate()` bodies, the five contexts, geometry, `SafetyMatrix`, `SpanProfile`, `BoundTightness`, `Counter` base | | the legacy DFS as one engine under `Inference/Dfs/`, see §7 |

### Generation (≈ 1 800 lines)

| From | File | Tag | Target name / note |
|---|---|---|---|
| Core | `FastGenerator.*` | keep | replaced by `SamplingEngine`, which exists; deleted when `app/igor` switches |
| Core | `GenModel::generate_sequences()`, `generate_sequences_fast()` | keep | free functions in `Generation/Legacy/`, split off in step 1c; until then `GenModel` stays whole in `Inference/Legacy/` |
| Model | `SamplingEngine.*`, `SamplingHandler.*`, `CategoricalSamplingHandler.*`, `MarkovSamplingHandler.*`, `SamplingHandlerFactory.*` | new | namespace becomes `igor::generation` |

### Streaming (unchanged, ≈ 3 700 lines)

Already in style; step 1b renames `igor` → `igor::streaming` and `igor::airr` →
`igor::streaming::airr`. Until `AlignmentData.h` exists (step 1c) it includes
`igor/Alignment/Legacy/Aligner.h` and links Alignment; afterwards it links Core only.

### Tests and apps

- `tst/igor/<Layer>/` and `tst/igor/<Layer>/Legacy/` mirror the sources; deleting a legacy
  file deletes its tests with it. The CTest labels added in `5ac1972` follow the layers.
- `app/igor` links Inference, Generation, Streaming, Alignment. It is the main consumer of
  `Legacy/` headers and the thing that keeps most **keep** files alive.
- `app/igor-convert` links Streaming only.
- `app/igor-demo` links Inference and Alignment.

## 6. Stage 1: the move, in three steps

Right after the `tk_refactoring` + `tandemD` merge lands. Steps 1a and 1b are two PRs, each
mechanical, each produced by a committed, re-runnable script (`scripts/migrate_layers/`), so
that Quentin applies the same scripts to `feature/TensorLinalg` and
`feature/MarginalRefactoring` instead of rebasing a 200-file rename by hand. Step 1c is a
queue of small PRs.

### 1a. Everyone into Legacy

`git mv` only, no content change except the include lines the move breaks.

- Every ex-Core file goes to `<Layer>/Legacy/` of the layer §5 assigns it to. No judgement at
  this step: `SafetyMatrix.h` and `Utils.h` travel the same way.
- Math, the tk Model files and Streaming stay where they are; the tk engine files move from
  `Model/` to `Inference/` and `Generation/` (layer directory, not `Legacy/`: they are new
  code).
- Old `Core` → `Model/Legacy` is a directory rename, and the other files move out of it,
  rather than the reverse: the Rec_Event cluster, which is where Quentin's branches live,
  keeps its history under one rename.
- `tst/igor/<Layer>/Legacy/` mirrors the sources.
- Each layer gets its CMake target, alias `igor::<Layer>`, and `target_link_libraries` per the
  "May depend on" column of §2. Three dependencies are transitional and listed as such in the
  CMake files, each removed by a named 1c item:

  | Transitional link | Reason | Removed by |
  |---|---|---|
  | Streaming → Alignment | `Alignment_data` is still in `Aligner.h` | `AlignmentData.h` extraction |
  | Inference → Generation | `GenModel.cpp` calls `FastGenerator` for `generate_sequences*` | the `GenModel` split |
  | Model → Alignment | `Rec_Event.h` and `QuerySequenceContext.h` include `Aligner.h` | `AlignmentData.h` extraction, then stage 2 |

Check: the build is identical, the Linux bitwise regression passes, `ctest -LE convergence`
passes per layer.

### 1b. Conventions: namespaces, includes, CMake

One PR, scripted, no file moves.

- Every file under `<Layer>/Legacy/` is wrapped in `namespace igor::<layer>::legacy { }`;
  the files already in a layer directory take `igor::<layer>` (Streaming: `igor` →
  `igor::streaming`, `igor::airr` → `igor::streaming::airr`; the tk engine files: `igor::model`
  → `igor::inference` / `igor::generation`). `NamespaceIndentation: None` in `.clang-format`,
  so the diff per file is two lines plus the includes, and Quentin's branches still rebase.
- Consumers (`app/`, `tst/`, the other layers) get `using namespace` or qualified names,
  whichever is less noise per file; the global namespace is empty of IGoR symbols at the end.
- Every include becomes `<igor/<Layer>/File.h>` or `<igor/<Layer>/Legacy/File.h>`. No
  `"File.h"`, no `<igor/Core/...>` for something that left Core.
- CMake made consistent with the two points above and with Model's current layout:
  - each layer's `CMakeLists.txt` is `project(<Layer>)`, one `add_library`, one
    `FILE_SET HEADERS` with `BASE_DIRS ${CMAKE_SOURCE_DIR}/src` listing the `Legacy/` headers
    next to the others, install rules identical across layers;
  - one export header per layer generated as `igor/<Layer>/Export.h` and included by that
    path, with macro `<LAYER>_EXPORT`. Core today generates `igorCoreExport.h` and includes
    it unqualified; it aligns on Model's scheme;
  - `src/igor/CMakeLists.txt` adds the subdirectories in dependency order: Core, Math,
    Alignment, Model, Inference, Generation, Streaming;
  - `tst/igor/<Layer>/CMakeLists.txt` labels its tests with the layer name, as `5ac1972` does
    for Model, Math and Streaming.

Check: same as 1a, plus `grep -rn "^#include \"" src` is empty and
`grep -rn "Legacy/" src app tst | wc -l` is recorded in the PR as the starting debt.

What executing 1b added to the rules above (2026-10-03, `scripts/migrate_layers/step1b_namespaces.py`):

- **Visibility between legacy namespaces** goes through using-directives placed right after each
  legacy namespace opening (`namespace igor::model::legacy { using namespace igor::core::legacy;
  using namespace igor::alignment::legacy; ... }`), following the DAG. They are transitional and
  leave with the 1c promotions. Consumers (tests, apps) nominate the same namespaces after their
  includes; new code qualifies instead: `legacy::X` from `igor::model`, `model::legacy::X` and
  `core::legacy::X` from the engines, fully qualified in Streaming, where `alignment::` would
  name `igor::streaming::airr::alignment`.
- **Core names no Model type any more.** `Next_event_ptr`, `Events_map` and
  `inverse_offset_comparator` moved from `Utils.h` to `Model/Legacy/EventTypedefs.h`.
- **`EventUtils` was three namespaces** (GeneticCode.h in Core, JournaledQuery.h in Alignment,
  EventUtils.h in Model); nominated together they are ambiguous. Core's is now `genetic_code`,
  Alignment's `journaled_query`, Model keeps `EventUtils`. The only rename of 1b.
- **Friendship across layers**: `Deletion` and `Gene_choice` befriend `Coverage_err_counter`,
  which is in Inference. The friend is named in full behind a forward declaration, an upward
  reference from Model to Inference that the observer rework of 1c removes.
- **`Model/Forward.h`** forward-declares the Model types the engines name (`RecombinationModel`,
  `Navigator`, `Topology`, `SampledScenario`, `SampledEvent`); the engines bring them in with
  using-declarations.
- **Two exceptions stay in the global namespace**: the `portable_getpid` / `portable_gethostid`
  shims of `Utils.h`, which sit in the `#if defined(_WIN32)` block that includes `winsock2.h`
  (a system header cannot be included inside a namespace), and `Math/HybridBuffer.h:78`, the one
  quoted include left, in a layer 1b does not touch. Both are 1c items.
- `using std::to_string;` inside `igor::core::legacy`: the legacy `to_string(Gene_class)`
  overloads would otherwise hide `std::to_string` from every namespace that nominates Core.

### 1c. Case by case

A queue of small changes, one commit each, applying the §5 tags. First pass executed on
2026-10-03 (branch `feature/layers-1c`), in this order:

| # | Item | Result |
|---|---|---|
| 1 | `Alignment_data` → `igor::core::AlignmentData` (`Core/AlignmentData.{h,cpp}`), with its CIGAR writer (`parseCigar`, `toCoreCigar`, `toExtendedCigar`) | Streaming no longer links Alignment |
| 2 | Quentin's 2026 headers promoted: `LayeredArray`, `SeqTypeRegistry`, `DynamicSequenceMap`, `SeqOffsetsMap`, `SegmentSpan` (Core); `SpanProfile`, `SafetyMatrix`, `UnfilledSegmentLengths`, `BoundTightness`, `JunctionGeometry` (Model) | ten headers out of `Legacy/`, ten stubs |
| 3 | `GenModel` split: generation becomes `igor::generation::legacy::SequenceGenerator` | Inference no longer links Generation; generated sequences bitwise identical |
| 4 | Core vocabulary: `Core/Types.h` (`EventType`, `SeqSide`, `SeqType`, `SeqTypeString`, `EventName`, `MarginalArrayPtr`, `SeqOffset`, `CodonTable`, `index_type`), `Core/IntStr`, `Core/GeneticCode`, `Core/Platform.h` | Core's only real legacy file left is `Utils`; nothing of IGoR at global scope |
| 5 | `test_RecombinationModel` back in the Model tests; `Matrix<T>` → `Math/Legacy/Matrix.h`; last quoted include fixed | Core's `Utils.h` defines no container |
| 6 | The events no longer befriend `Coverage_err_counter`: three public slot accessors | no upward reference left from Model to Inference |

How a promotion works (`scripts/migrate_layers/promote.py`): the header moves to the layer
directory and namespace; a **stub** stays at the old `Legacy/` path, which includes the promoted
header and re-declares its names in `igor::<layer>::legacy`, so legacy consumers keep their
include line and their unqualified names. A stub goes when its last consumer switches.
`check_directives.py` removes the using-directives that a promotion leaves without a target.

Decisions taken while executing:

- **`LegacyBridge` keeps its name.** The rename to `MarginalsBridge` was proposed before the
  Legacy policy existed. Now that `igor::model::legacy` is a namespace, a bridge between
  `legacy::Model_Parms` / `legacy::Model_marginals` and `Topology` / `RecombinationModel` is
  exactly what the name says.
- **`FastSampling` stays in `Generation/Legacy`.** The file mixes the samplers with buffered file
  writers and thread pools; moving it whole would put I/O in Math, and splitting it is a
  replace, not a promote.
- **Alignment links Math**, for the substitution `Matrix`. The DAG gains Alignment → Math.
- **Method and free-function names are unchanged** in the promoted headers (still snake_case).
  Renaming them would touch the `iterate()` bodies; it goes with stage 2.
- The promoted `AlignerInternal.h` question: kept, see §5.

What remains transitional: Model → Alignment, for `JournaledQuery` held by
`QuerySequenceContext` and the nucleotide helpers the events call. It goes with stage 2.

What remains in `Legacy/` after this pass, by tag: **keep** (`Utils`, `Matrix`, `EventUtils`,
`gene_to_seqtype_migr`, `Model_Parms`, `Model_marginals`, `FastGenerator`, `FastSampling`);
**promote, not done yet** (`Aligner`, `ExtractFeatures`, `CDR3SeqData`, `JournaledQuery`,
`JsonDetail`, `ModelJson`, and the `Gene_class` enums of `Utils.h`); **replace, stage 2** (the
events, the contexts, the error rates, `Counter` and the four counters, `GenModel`,
`SequenceGenerator`).

## 7. Stage 2: the DFS becomes an engine

This is not a layering task. It is the handler extraction the synthesis calls for ("handlers
are declared per node", "the legacy DFS is one inference engine among several"). The layering
only gives it a destination: `Rec_Event::iterate()` and the state it threads (the five
contexts, geometry, safety matrix, span profiles, bound instrumentation, the `Counter`
interface) leave Model for `Inference/Dfs/`. `Rec_Event` is then a definition: realizations,
attributes, capabilities, JSON. That is when the **replace** copies of the events, contexts
and error models are made. Model shrinks to about 8 000 lines and no longer sees
`AlignmentData`.

Who and when: Quentin's zone, after his repair queue (R-steps) and B5/B11 are stable, since it
touches the same four `iterate()` bodies. Not before.

## 8. Open points

1. `Alignment` or `Sequence` for the aligner layer. `Alignment` is proposed: the layer is about
   the algorithms, the sequence types themselves are vocabulary and live in Core.
2. Whether `Math` stays header-only. Moving `FastSampling` in makes it a compiled library unless
   it is templated. Either is fine; decide at move time.
3. `Model_marginals` vs `RecombinationModel<T>`: both survive stage 1 as two parameter layouts
   of one model. Which one the DFS engine reads after stage 2 is the marginal-storage contract
   still open between the branches (`feature/MarginalRefactoring`); the layering does not
   decide it and does not need to.
