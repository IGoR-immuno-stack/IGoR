# Model Module

The Model module implements the modern probabilistic-graph layer of IGoR.
It owns the **graph topology**, the **parameter tensors**, and the **engines**
that operate on them (inference and generation).

---

## Architecture Overview

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                Model module                                  │
│                                                                              │
│                      ┌─────────────────────────────────┐                     │
│                      │    RecombinationModel<T>        │                     │
│                      │                                 │                     │
│                      │  unique_ptr<const Topology>     │                     │
│                      │  ┌────────────┐                 │                     │
│                      │  │  Topology  │  DAG of         │                     │
│                      │  │            │  Rec_Event      │                     │
│                      │  └────────────┘                 │                     │
│                      │                                 │                     │
│                      │  vector<Tensor<T>>  m_weights   │                     │
│                      │  ┌────┐ ┌────┐ ┌────┐  ...      │                     │
│                      │  │ T₀ │ │ T₁ │ │ T₂ │           │                     │
│                      │  └──┬─┘ └──┬─┘ └──┬─┘           │                     │
│                      │     │      │      │             │                     │
│                      │     └──────┼──────┘             │                     │
│                      └────────────┼────────────────────┘                     │
│                                   │                                          │
│          ┌─── borrow ─────────────┴───────────── borrow ─────┐               │
│          │   Tensor<T>&                     const Tensor<T>& │               │
│          │   (mutable)                      (immutable)      │               │
│          ▼                                                   ▼               │
│  ┌─────────────────────────────┐  ┌─────────────────────────────┐            │
│  │     InferenceEngine<T>      │  │     SamplingEngine<T>       │            │
│  │                             │  │                             │            │
│  │  shared_ptr<                │  │  shared_ptr<const           │            │
│  │    RecombinationModel<T>>   │  │    RecombinationModel<T>>   │            │
│  │                             │  │                             │            │
│  │  owns N handlers:           │  │  owns N handlers:           │            │
│  │  ┌────────────────────────┐ │  │  ┌────────────────────────┐ │            │
│  │  │ InferenceHandler<T>    │ │  │  │ SamplingHandler<T>     │ │            │
│  │  │ (abstract)             │ │  │  │ (abstract)             │ │            │
│  │  └───┬───────────┬────────┘ │  │  └───┬───────────┬────────┘ │            │
│  │      │           │          │  │      │           │          │            │
│  │      ▼           ▼          │  │      ▼           ▼          │            │
│  │  ┌─────────┐ ┌──────────┐   │  │  ┌─────────┐ ┌──────────┐   │            │
│  │  │Categori-│ │ Markov   │   │  │  │Categori-│ │ Markov   │   │            │
│  │  │cal      │ │ Inference│   │  │  │cal      │ │ Sampling │   │            │
│  │  │Inference│ │ Handler  │   │  │  │Sampling │ │ Handler  │   │            │
│  │  │Handler  │ │ <T>      │   │  │  │Handler  │ │ <T>      │   │            │
│  │  └─────────┘ └──────────┘   │  │  └─────────┘ └──────────┘   │            │
│  └─────────────────────────────┘  └─────────────────────────────┘            │
│                                                                              │
│  ┌──────────────┐  ┌───────────────────┐  ┌───────────────────┐              │
│  │  Navigator   │  │ Handler Factories │  │   Event Factory   │              │
│  │  <NodeType>  │  │ (Sampling +       │  │                   │              │
│  │              │  │  Inference)       │  │                   │              │
│  └──────────────┘  └───────────────────┘  └───────────────────┘              │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────┐         │
│  │                          LegacyBridge                           │         │
│  │  Model_Parms  ⇄  Topology                                       │         │
│  │  Model_marginals → RecombinationModel                           │         │
│  └─────────────────────────────────────────────────────────────────┘         │
│                                                                              │
│           Depends on:  igor::Core  (Rec_Event, Model_Parms, Model_marginals) │
│                        igor::Math  (Tensor<T>)                               │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Design Choices

### 1. Single Source of Truth — `RecombinationModel<T>`

All probability tensors live in one place: `RecombinationModel<T>`.
Engines and handlers never copy these tensors — they **borrow references**
to them. This guarantees that when the inference M-step writes normalised
values, every component sees the update immediately.

```
RecombinationModel<T>
├── unique_ptr<const Topology>        (exclusive ownership)
├── vector<Tensor<T>> m_weights       (one per event, indexed by UID)
│       ▲ mutable ref          ▲ const ref
│       │                      │
│  InferenceHandler<T>    SamplingHandler<T>
│  (writes back during    (reads only, builds
│   M-step)                CDF tables)
```

**Example — constructing a model from files:**

```cpp
#include <igor/Model/RecombinationModel.h>

auto model = igor::model::recombination_model_from_files<double>(
    "models/mouse_tcr_beta/model_parms.txt",
    "models/mouse_tcr_beta/model_marginals.txt");

// Access the probability tensor for a specific event
auto& v_weights = model.weight("v_choice");
std::cout << "V-gene tensor shape: " << v_weights.shape() << "\n";

// Iterate over all weights in topological order
for (const auto& tensor : model.orderedWeights()) {
    std::cout << tensor.size() << " elements\n";
}
```

### 2. Strict Ownership Hierarchy

Ownership is modelled with standard smart pointers and follows a strict tree:

| Owner | Owns | Via |
|---|---|---|
| `RecombinationModel<T>` | `Topology` | `unique_ptr<const Topology>` |
| `RecombinationModel<T>` | Weight tensors | `vector<Tensor<T>>` |
| `SamplingEngine<T>` | `RecombinationModel` | `shared_ptr<const RecombinationModel<T>>` |
| `SamplingEngine<T>` | Sampling handlers | `vector<unique_ptr<SamplingHandler<T>>>` |
| `InferenceEngine<T>` | `RecombinationModel` | `shared_ptr<RecombinationModel<T>>` (mutable) |
| `InferenceEngine<T>` | Inference handlers | `vector<unique_ptr<InferenceHandler<T>>>` |

The key asymmetry: `SamplingEngine` holds a `const` shared pointer (read-only
access to the model), while `InferenceEngine` holds a **non-const** shared
pointer (so the M-step can write normalised values back into the model's
tensors through the borrowed mutable references).

**Example — sharing a model between engines:**

```cpp
auto model = std::make_shared<RecombinationModel<double>>(
    recombination_model_from_files("model_parms.txt", "model_marginals.txt"));

// Inference engine borrows mutable references
InferenceEngine<double> inference(model);

// Sampling engine borrows const references — accepts shared_ptr<RecombinationModel<T>>
// directly, the const conversion is handled internally by the constructor.
SamplingEngine<double> sampling(model);

// After an M-step in the inference engine, the sampling engine's handlers
// still point to the same tensors — they see the updated values.
inference.updateParameters();
```

### 3. Borrowed-Reference Pattern

Handlers never copy the probability tensors. Instead, each handler stores a
**reference** to its tensor in the `RecombinationModel`. The handler's
constructor receives this reference, and the handler's lifetime is tied to
its owning engine, which in turn holds a `shared_ptr` to the model — so the
reference is always valid.

For inference handlers, the reference is mutable (`Tensor<T>&`):

```cpp
class CategoricalInferenceHandler : public InferenceHandler<T> {
    CategoricalInferenceHandler(std::string name, index_type uid,
                                 math::Tensor<T>& weights);  // mutable ref
private:
    math::Tensor<T>& m_weights;      // borrowed from RecombinationModel
    math::Tensor<T>  m_accumulator;   // owned — only this is the handler's own state
};
```

For sampling handlers, the reference is const (`const Tensor<T>&`):

```cpp
class CategoricalSamplingHandler : public SamplingHandler<T> {
    CategoricalSamplingHandler(std::string name, index_type uid,
                                const math::Tensor<T>& weights);  // const ref
private:
    const math::Tensor<T>& m_weights;  // borrowed, read-only
    math::Tensor<T>        m_cdfs;     // owned — precomputed CDF tables
};
```

### 4. Two Parallel Handler Hierarchies

The module provides **two independent abstract hierarchies** with
intentionally different APIs:

| | Inference side | Sampling side |
|---|---|---|
| **Base class** | `InferenceHandler<T>` | `SamplingHandler<T>` |
| **Categorical** | `CategoricalInferenceHandler<T>` | `CategoricalSamplingHandler<T>` |
| **Markov** | `MarkovInferenceHandler<T>` | `MarkovSamplingHandler<T>` |

Why two hierarchies rather than one?

* **Different responsibilities.**
  The inference handler holds a mutable reference to the weight tensor *and*
  owns an accumulator tensor, exposes `resetAccumulator()` /
  `maximizeLikelihood()` for the EM loop. The sampling handler holds a const
  reference to the weight tensor plus precomputed CDF tables and exposes
  `sample()` / `sampleSequence()`.

* **Different constness.**
  A `SamplingHandler` is conceptually frozen after `precomputeCDF()` is called;
  the generation loop never writes to the model. An `InferenceHandler` is
  mutated on every E-step and writes back during the M-step.

**Example — inference handler EM cycle:**

```cpp
InferenceEngine<double> engine(model);

// E-step preparation
engine.resetAccumulators();

// … accumulate counts into each handler's accumulator tensor …

// M-step: normalise accumulators → write into model weights
engine.updateParameters();
```

**Example — sampling handler generation cycle:**

```cpp
SamplingEngine<double> engine(model);

std::mt19937_64 rng(42);
auto scenario = engine.run(rng);

// scenario.events[uid].indices holds the sampled values
std::size_t v_gene_index = scenario.index_of(v_choice_uid);
```

### 5. Topology-Driven Construction via Factories

No handler is created manually. Both engines delegate construction to factory
functions that inspect the `Topology` and each event's type and shape.

The factory inspects each node in the topology, computes the correct tensor
shape (`[own_dims, parent1_dims, parent2_dims, ...]`), looks up the
registered creator for that event type, and constructs the appropriate handler:

```
RecombinationModel<T>
        │
        ├─► inference_handler_factory::build(model)
        │       for each topology node:
        │           shape = event.inherent_shape() + parent_shapes
        │           handler = create<T>(event_type, event, weight_ref)
        │       → vector<unique_ptr<InferenceHandler<T>>>
        │
        └─► sampling_handler_factory::build(model)
                for each topology node:
                    shape = event.inherent_shape() + parent_shapes
                    handler = create<T>(event_type, event, const_weight_ref)
                → vector<unique_ptr<SamplingHandler<T>>>
```

Each concrete handler type registers itself at static-initialisation time via
a `Registrar` object — no central switch statement:

```cpp
// InferenceHandlerFactory.cpp
static Registrar<double, CategoricalInferenceHandler<double>>
    categorical_registrar{ GeneChoice_t, Deletion_t, Insertion_t };

static Registrar<double, MarkovInferenceHandler<double>>
    markov_registrar{ Dinuclmarkov_t };
```

Adding a new event type requires only a new handler class and a static
`Registrar` — no existing code is modified (open/closed principle).

`EventFactory` follows the same self-registering pattern for `Rec_Event`
subclasses (`Gene_choice`, `Deletion`, `Insertion`, `Dinucl_markov`), used
by `read_topology()` when parsing model parameter files.

### 6. Navigator — Generic Index-Based View

`Navigator<NodeType, PtrType>` is a lightweight, non-owning **view** over a
subset of nodes in a parallel vector. It is constructed from a reference to
the full node vector and an index list (parent IDs, child IDs, or topological
order). It provides a full random-access iterator and supports range-based
`for` loops.

This single template serves four different node types:

| Instantiation | Used by |
|---|---|
| `Navigator<Rec_Event>` | `Topology::parents()`, `Topology::children()` |
| `Navigator<SamplingHandler<T>, unique_ptr<…>>` | `SamplingEngine::orderedHandlers()`, `parents()`, `children()` |
| `Navigator<InferenceHandler<T>, unique_ptr<…>>` | `InferenceEngine::orderedHandlers()`, `parents()`, `children()` |
| `Navigator<Tensor<T>, Tensor<T>>` | `RecombinationModel::orderedWeights()` |

**Example — iterating handlers in topological order:**

```cpp
InferenceEngine<double> engine(model);

// orderedHandlers() returns a Navigator — parents always appear before children
for (const auto& handler_ptr : engine.orderedHandlers()) {
    std::cout << handler_ptr->name() << "  uid=" << handler_ptr->uid() << "\n";
}
```

**Example — resolving parent realizations during sampling:**

```cpp
// Inside SamplingEngine::run() — simplified
for (const auto& handler_ptr : orderedHandlers()) {
    auto uid = handler_ptr->uid();

    // parents(uid) returns a Navigator over only the parent handlers
    std::vector<std::size_t> parent_reals;
    for (const auto& parent_ptr : parents(uid)) {
        parent_reals.push_back(scenario.index_of(parent_ptr->uid()));
    }

    std::size_t val = handler_ptr->sample(rng, parent_reals);
    scenario.events[uid] = { uid, {val} };
}
```

**Example — walking the Topology graph:**

```cpp
auto topology = read_topology("model_parms.txt");

// Children of event 0
for (const auto& child : topology->children(0)) {
    std::cout << child->get_nickname() << "\n";
}

// Random access into parents
auto p = topology->parents(3);
if (!p.empty()) {
    std::cout << "first parent: " << p[0]->get_nickname() << "\n";
}
```

### 7. LegacyBridge — Adapter to Core Types

`LegacyBridge` is a set of **free functions** that convert between the legacy
`Core` types and the modern `Model` types. It is the **only place** where
`Model` code includes `Model_Parms.h` and `Model_marginals.h`.

| Direction | Function | Purpose |
|---|---|---|
| Legacy → Modern | `import_from_legacy(const Model_Parms&)` | Build a `Topology` from legacy event graph |
| Modern → Legacy | `export_to_legacy(const Topology&)` | Export a `Topology` back to `Model_Parms` |
| Legacy → Model  | `import_from_legacy(RecombinationModel<T>&, const Model_marginals&)` | Fill tensors from legacy flat array |

**Example — round-trip through the bridge:**

```cpp
// Legacy → Modern
Model_Parms legacy_parms;
// … populate legacy_parms …
auto topology = import_from_legacy(legacy_parms);

// Modern → Legacy
auto parms_copy = export_to_legacy(*topology);

// Load marginals into a RecombinationModel
RecombinationModel<double> model(std::move(topology));
import_from_legacy(model, legacy_marginals);
```

### 8. Scenario — Value Objects

`SampledEvent` and `SampledScenario` are plain value types that decouple the
sampling result from the engine that produced it.

```cpp
struct SampledEvent {
    index_type               event_id;  // topology node UID
    std::vector<std::size_t> indices;   // sampled realization indices
    // 1 element for categorical events, N for Markov chains
};

struct SampledScenario {
    std::vector<SampledEvent> events;   // indexed by topology UID
    std::size_t index_of(index_type i) const;  // shorthand for events[i].indices[0]
};
```

**Example:**

```cpp
auto scenario = sampling_engine.run(rng);

// Categorical event: single index
std::size_t v_gene = scenario.index_of(v_choice_uid);

// Markov event: full nucleotide chain
auto& chain = scenario.events[dinucl_uid].indices;
// chain = {first_nucleotide, nt2, nt3, ...}
```

---

## Concrete Handler Details

### Categorical Handlers

Used for `Gene_choice`, `Deletion`, and `Insertion` events.

**Tensor shape:** `[n_realizations, parent1_dim, parent2_dim, ...]`

| Aspect | Inference | Sampling |
|---|---|---|
| **Weight ref** | mutable `Tensor<T>&` | const `Tensor<T>&` |
| **Owned state** | `m_accumulator` (same shape) | `m_cdfs` (precomputed CDF table) |
| **M-step** | Normalise axis 0: each column of own realizations sums to 1 per parent combination | N/A |
| **Sampling** | N/A | Binary search on CDF row for the selected parent slice |
| **Key method** | `maximizeLikelihood()` | `sample(rng, parent_indices)` |
| **Accessor** | `realizationCount()` | — |

### Markov Handlers

Used for `Dinucl_markov` events (dinucleotide transition matrices).

**Tensor shape:** `[n_states, n_states, parent1_dim, ...]` — `(from, to, parents)`

| Aspect | Inference | Sampling |
|---|---|---|
| **Weight ref** | mutable `Tensor<T>&` | const `Tensor<T>&` |
| **Owned state** | `m_accumulator` (same shape) | `m_row_cdfs` + `m_first_cdf` (stationary marginal) |
| **M-step** | Normalise axis 1: each row of "to" states sums to 1 per "from" state and parent combination | N/A |
| **Sampling** | N/A | Two modes: (1) empty parents → sample first nucleotide from marginal; (2) `parent_indices[0]` = from_state → sample next state from row CDF |
| **Key methods** | `maximizeLikelihood()` | `sample()`, `sampleSequence(rng, first_state, n_steps, ...)` |
| **Accessor** | `stateCount()` | — |

---

## File Inventory

| File | Role |
|---|---|
| `Topology.h / .cpp` | DAG of `Rec_Event` nodes; Kahn's topological sort; edge operations; `read_topology()` file parser |
| `Navigator.h` | Generic index-based view with full random-access iterator |
| `RecombinationModel.h / .tpp` | Pairs a `Topology` (unique ownership) with one `Tensor<T>` per node; `orderedWeights()` Navigator; `read_parameters()` and `recombination_model_from_files()` |
| `InferenceHandler.h` | Abstract base for inference handlers (mutable weight ref + accumulator) |
| `CategoricalInferenceHandler.h / .tpp` | Categorical distribution handler (inference): axis-0 normalisation |
| `MarkovInferenceHandler.h / .tpp` | Markov transition matrix handler (inference): axis-1 normalisation |
| `InferenceEngine.h / .tpp` | Owns `InferenceHandler<T>` vector; drives EM loop; Navigator-based iteration |
| `InferenceHandlerFactory.h / .tpp / .cpp` | Self-registering abstract factory for inference handlers |
| `SamplingHandler.h / .tpp` | Abstract base for generation handlers (const weight ref + CDF tables) |
| `CategoricalSamplingHandler.h / .tpp` | CDF-based categorical sampler; precomputed per parent slice |
| `MarkovSamplingHandler.h / .tpp` | Row-CDF Markov chain sampler; `sampleSequence()` for chains |
| `SamplingEngine.h / .tpp` | Owns `SamplingHandler<T>` vector; `run()` generates a `SampledScenario`; Navigator-based iteration |
| `SamplingHandlerFactory.h / .tpp / .cpp` | Self-registering abstract factory for sampling handlers |
| `EventFactory.h / .cpp` | Self-registering abstract factory for `Rec_Event` subclasses |
| `Scenario.h` | `SampledEvent` and `SampledScenario` value types |
| `LegacyBridge.h / .tpp` | Adapter: `Topology ⇄ Model_Parms`, `Model_marginals → RecombinationModel` |

---

## Data Flow

### Inference Path

```
model_parms.txt ──► read_topology() ──► Topology
                                            │
                                     RecombinationModel<T>
                                            │
model_marginals.txt ──► read_parameters() ──┘
                                            │
                                    InferenceEngine<T>
                                     (handlers borrow mutable refs)
                                           │
                                    ┌──────┴──────┐
                                    │  EM Loop    │
                                    │             │
                              resetAccumulators() │
                                    │             │
                              [ E-step: accumulate│counts ]
                                    │             │
                              updateParameters()  │
                                    │  (normalise │accumulators
                                    │   → write   │into model)
                                    └──────┬──────┘
                                           │
                                    Model weights updated in place.
                                    All handlers see the new values.
```

### Generation Path

```
model_parms.txt ──► read_topology() ──► Topology
                                           │
                                    RecombinationModel<T>
                                           │
model_marginals.txt ──► read_parameters() ──┘
                                           │
                                    SamplingEngine<T>
                                     (handlers borrow const refs,
                                      precompute CDFs)
                                           │
                                    engine.run(rng)
                                           │
                                    ┌──────┴──────┐
                                    │  per event  │
                                    │  (topo      │
                                    │   order):   │
                                    │             │
                                    │  resolve    │
                                    │  parent     │
                                    │  realizations│
                                    │       │     │
                                    │  sample()   │  ◄─ polymorphic
                                    │  or sample- │     (categorical
                                    │  Sequence() │      or Markov)
                                    └──────┬──────┘
                                           │
                                           ▼
                                    SampledScenario
```

### One-Step Model Loading

```cpp
// Convenience function that combines read_topology + read_parameters
auto model = recombination_model_from_files<double>(
    "model_parms.txt", "model_marginals.txt");

// Both engines can be constructed directly from the model
auto shared = std::make_shared<RecombinationModel<double>>(std::move(model));
InferenceEngine<double> inference(shared);
```

---

## Build

The module produces a shared library `igorModel` (CMake target `igor::Model`)
that links publicly against `igor::Core` and `igor::Math`.

```
igor::Core ◄── igor::Model ──► igor::Math
```

Compiled sources: `EventFactory.cpp`, `InferenceHandlerFactory.cpp`,
`SamplingHandlerFactory.cpp`, `Topology.cpp`. All other code is
header-only (templates in `.h` / `.tpp` files).
