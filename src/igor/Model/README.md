# Model Module

The Model module implements the modern probabilistic-graph layer of IGoR.
It owns the **graph topology**, the **parameter tensors**, the **engines** that
operate on them (inference, generation), and the **bridge** that connects
everything to the legacy `Core` types (`Model_Parms`, `Model_marginals`).

---

## Architecture Overview

```
┌──────────────────────────────────────────────────────────────────────────┐
│                           Model module                                   │
│                                                                          │
│                         ┌────────────┐                                   │
│                         │  Topology  │                                   │
│                         │  (DAG of   │                                   │
│                         │ Rec_Event) │                                   │
│                         └─────┬──────┘                                   │
│                               │  shared_ptr<const Topology>              │
│                    ┌──────────┴──────────┐                               │
│                    │                     │                               │
│                    ▼                     ▼                               │
│       ┌────────────────────┐   ┌────────────────────┐                    │
│       │ InferenceEngine<T> │   │ SamplingEngine<T>  │                    │
│       │                    │   │                    │                    │
│       └───────────┬────────┘   └────────┬───────────┘                    │
│                   │                     │                                │
│          owns N handlers       owns N handlers                           │
│                   │                     │                                │
│                   ▼                     ▼                                │
│       ┌────────────────────┐   ┌────────────────────┐                    │
│       │ MarginalHandler<T> │   │ SamplingHandler<T> │                    │
│       │  (abstract)        │   │  (abstract)        │                    │
│       └───┬───────────┬────┘   └─────┬─────────┬────┘                    │
│           │           │              │         │                         │
│           ▼           ▼              ▼         ▼                         │
│ ┌──────────────┐ ┌──────────┐ ┌──────────────┐ ┌──────────────┐          │
│ │ Categorical  │ │ Markov   │ │ Categorical  │ │   Markov     │          │
│ │ Handler<T>   │ │Handler<T>│ │ Sampling     │ │  Sampling    │          │
│ │              │ │          │ │ Handler<T>   │ │  Handler<T>  │          │
│ └──────────────┘ └──────────┘ └──────────────┘ └──────────────┘          │
│                                                                          │
│  ┌──────────────┐  ┌───────────────────┐  ┌──────────────┐               │
│  │  Navigator   │  │ SamplingHandler   │  │   Event      │               │
│  │  <NodeType>  │  │ Factory           │  │   Factory    │               │
│  └──────────────┘  └───────────────────┘  └──────────────┘               │
│                                                                          │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │                    LegacyBridge                          │            │
│  │                                                          │            │
│  │  Model_Parms  ⇄  Topology                                │            │
│  │  Model_marginals → InferenceEngine                       │            │
│  │  Model_marginals → SamplingEngine                        │            │
│  └──────────────────────────────────────────────────────────┘            │
└──────────────────────────────────────────────────────────────────────────┘

       Depends on:  igor::Core  (Rec_Event, Model_Parms, Model_marginals)
                    igor::Math  (Tensor<T>)
```

---

## Design Philosophy

### 1. Separation of Structure and Computation

The **model graph** (which events exist and how they condition each other) is
strictly separated from the **numerical engines** (inference, generation) that
operate on it:

| Concern | Class | Owns |
|---|---|---|
| Graph structure | `Topology` | `Rec_Event` nodes, adjacency lists |
| Inference weights + EM | `InferenceEngine<T>` | `MarginalHandler<T>` instances |
| Generation / sampling | `SamplingEngine<T>` | `SamplingHandler<T>` instances |

Both engines receive a `shared_ptr<const Topology>` at construction and never
mutate it. This makes the topology safely shareable between engines and across
threads.

### 2. Two Parallel Handler Hierarchies

The module deliberately provides **two independent abstract hierarchies** with
intentionally different APIs:

| | Inference side | Sampling side |
|---|---|---|
| **Base class** | `MarginalHandler<T>` | `SamplingHandler<T>` |
| **Categorical** | `CategoricalHandler<T>` | `CategoricalSamplingHandler<T>` |
| **Markov** | `MarkovHandler<T>` | `MarkovSamplingHandler<T>` |

Why two hierarchies rather than one?

* **Different responsibilities.**
  The inference handler owns a `parameters` tensor *and* an `accumulator`
  tensor, exposes `reset_accumulator()` / `maximize_likelihood()` for the E-M
  loop, and serialises to the legacy text format. The sampling handler owns
  only the probability tensor plus precomputed CDF tables and exposes
  `sample()` / `sampleSequence()`.

* **Different scalar types.**
  Inference typically runs in `long double` (to match the legacy
  `Marginal_array_p` format), while sampling runs in `double`.

* **Immutable after setup.**
  A `SamplingHandler` is conceptually frozen after `precomputeCDF()` is called;
  the generation loop never writes to it. An inference handler is mutated on
  every E-step.

### 3. Topology-Driven Construction via Factories

No handler is created manually. Both engines delegate construction to factory
functions that inspect the `Topology` and each event's type and shape:

```
Topology  ──►  SamplingHandlerFactory::build<T>(topology)
                   │
                   ├─ for each node:
                   │    shape = event.inherent_shape()
                   │            + [parent₁.inherent_shape()…]
                   │
                   │    handler = factory.create<T>(event.type, event, shape)
                   │    handler.setUid(uid)
                   │
                   └─ returns vector<unique_ptr<SamplingHandler<T>>>
```

Both `EventFactory` (for `Rec_Event` construction) and
`SamplingHandlerFactory` (for `SamplingHandler` construction) use the
**self-registering factory** pattern described below.

---

## Design Patterns

### Strategy — Handler Polymorphism

`MarginalHandler<T>` and `SamplingHandler<T>` are abstract strategy
interfaces. The concrete strategies (`Categorical*`, `Markov*`) encapsulate
the distribution-specific logic (normalisation axis, CDF layout, transition
matrices) behind a uniform API.

The engines iterate over handlers polymorphically and never need to know which
distribution family an event uses:

```cpp
// SamplingEngine::run()
for (auto& handler : orderedHandlers()) {
    auto parent_reals = resolve_parents(handler->uid());
    std::size_t val = handler->sample(rng, parent_reals);   // strategy call
}
```

### Abstract Factory — Self-Registering Creators

`SamplingHandlerFactory` and `EventFactory` implement the **Abstract Factory**
pattern with **static self-registration**.

Each concrete type registers a creator lambda at static-initialisation time:

```cpp
// SamplingHandlerFactory.cpp
static Registrar<double, CategoricalSamplingHandler<double>>
    categorical_registrar{ GeneChoice_t, Deletion_t, Insertion_t };

static Registrar<double, MarkovSamplingHandler<double>>
    markov_registrar{ Dinuclmarkov_t };
```

This makes the factory **open for extension**: adding a new event type
requires only a new handler class and a static `Registrar`; no existing code
is modified.

**Example — `SamplingHandlerFactory::build()` constructs all handlers from a
`Topology`:**

```cpp
template <typename T>
std::vector<HandlerPtr<T>> build(const Topology& topology)
{
    std::vector<HandlerPtr<T>> handlers(topology.size());

    for (index_type uid = 0; uid < static_cast<index_type>(topology.size()); ++uid) {

        auto event = topology.event(uid);

        // Shape = event's own dimensions + each parent's dimensions
        std::vector<std::size_t> shape = event->inherent_shape();
        for (index_type parent_uid : topology.parentsIds(uid)) {
            auto parent_shape = topology.event(parent_uid)->inherent_shape();
            shape.insert(shape.end(), parent_shape.begin(), parent_shape.end());
        }

        // Dispatches to the registered creator for event->get_type()
        // (CategoricalSamplingHandler or MarkovSamplingHandler)
        auto handler = create<T>(event->get_type(), event, shape);
        handler->setUid(uid);

        handlers[uid] = std::move(handler);
    }
    return handlers;
}
```

**Example — `EventFactory` self-registration and use in `read_topology()`:**

```cpp
// EventFactory.cpp  — registration at static-init time
static Registrar<GeneChoice_t,  Gene_choice>   gene_choice_registrar;
static Registrar<Deletion_t,    Deletion>       deletion_registrar;
static Registrar<Insertion_t,   Insertion>       insertion_registrar;
static Registrar<Dinuclmarkov_t, Dinucl_markov> dinucl_markov_registrar;

// Topology.cpp  — called while parsing model_parms.txt
auto event = event_factory::create(GeneChoice_t);   // returns shared_ptr<Rec_Event>
event->set_nickname("v_choice");
// … populate realizations …
index_type uid = topology->addEvent(event);
```

### Concept-Constrained Generic View — Navigator

`Navigator<NodeType, PtrType>` is a lightweight, non-owning **view** over a
subset of nodes in a parallel vector. It is parameterised by a C++20 concept
`HasUid`:

```cpp
template <typename T>
concept HasUid = requires(T& t, const T& ct, igor::index_type id) {
    { ct.uid()     } -> std::same_as<igor::index_type>;
    { t.setUid(id) } -> std::same_as<void>;
};
```

This single template serves three different node types:

| Instantiation | Used by |
|---|---|
| `Navigator<Rec_Event>` | `Topology::parents()`, `Topology::children()` |
| `Navigator<MarginalHandler<T>>` | (reserved for future `InferenceEngine` graph traversal) |
| `Navigator<SamplingHandler<T>, unique_ptr<…>>` | `SamplingEngine::parents()`, `orderedHandlers()` |

A `Navigator` is constructed from a reference to the full node vector and an
index list (parent IDs, child IDs, or topological order). It satisfies the
random-access-range concept and supports range-based `for` loops.

**Example — resolving parent realizations in `SamplingEngine::run()`:**

```cpp
SampledScenario SamplingEngine<T>::run(std::mt19937_64& rng) const
{
    SampledScenario result(m_topology->size());

    // orderedHandlers() returns a Navigator over the topological order
    for (const auto& handler_ptr : this->orderedHandlers()) {

        const auto& current = *handler_ptr;
        index_type  uid     = current.uid();

        // parents(uid) returns a Navigator<SamplingHandler<T>, unique_ptr<…>>
        // over only the parent handlers of this node
        auto parents_nav = this->parents(uid);

        std::vector<std::size_t> parent_reals;
        parent_reals.reserve(parents_nav.size());

        for (const auto& parent_ptr : parents_nav) {
            // Each element is a const unique_ptr<SamplingHandler<T>>&
            parent_reals.push_back(result.index_of(parent_ptr->uid()));
        }

        // Polymorphic sample: categorical or Markov, decided at runtime
        std::size_t val = current.sample(rng, parent_reals);
        result.events[uid] = { uid, {val} };
    }
    return result;
}
```

**Example — walking the `Topology` graph with `Navigator<Rec_Event>`:**

```cpp
auto topology = read_topology("model_parms.txt");

// Iterate over children of event 0
for (const auto& child : topology->children(0)) {
    std::cout << child->get_nickname()                    // range-for works
              << "  uid=" << child->uid() << "\n";
}

// Random access
auto parent_nav = topology->parents(3);
if (!parent_nav.empty()) {
    auto& first_parent = parent_nav[0];   // subscript operator
    std::cout << "first parent: " << first_parent->get_nickname() << "\n";
}
```

### Adapter — LegacyBridge

`LegacyBridge` is a **two-way adapter** between the legacy `Core` types and
the modern `Model` types:

| Direction | Function |
|---|---|
| Legacy → Modern | `import_from_legacy(Model_Parms) → Topology` |
| Modern → Legacy | `export_to_legacy(Topology) → Model_Parms` |
| Legacy → Inference | `import_from_legacy(InferenceEngine, Model_marginals, Model_Parms)` |
| Inference → Legacy | `export_to_legacy(InferenceEngine, Model_marginals, Model_Parms)` |
| Legacy → Sampling | `import_from_legacy(SamplingEngine, Model_marginals, Topology)` |
| File → Sampling | `read_parameters(filename, SamplingEngine)` — standalone parser |

The bridge is the **only place** where `Model` code includes `Core` headers
such as `Model_Parms.h` and `Model_marginals.h`. Engines and handlers are
agnostic to the legacy flat-array representation.

### Value Object — Scenario

`SampledScenario` and `SampledEvent` are plain value types with no behaviour
beyond storage and the `index_of()` convenience accessor. They decouple the
sampling result from the engine that produced it.

---

## File Inventory

| File | Role |
|---|---|
| `Topology.h / .cpp` | DAG of `Rec_Event` nodes; Kahn's topological sort; edge operations; `read_topology()` file parser |
| `Navigator.h` | Generic index-based view with full random-access iterator; `HasUid` concept |
| `MarginalHandler.h` | Abstract base for inference handlers; `EventDescriptor` struct |
| `CategoricalHandler.h / .tpp` | Categorical distribution handler (inference): axis-0 normalisation |
| `MarkovHandler.h / .tpp` | Markov transition matrix handler (inference): axis-1 normalisation |
| `InferenceEngine.h / .tpp` | Owns `MarginalHandler<T>` vector; drives E-M loop operations |
| `SamplingHandler.h / .tpp` | Abstract base for generation handlers; `rawData()` / `rawDataSize()` API |
| `CategoricalSamplingHandler.h / .tpp` | CDF-based categorical sampler; precomputed per parent slice |
| `MarkovSamplingHandler.h / .tpp` | Row-CDF Markov chain sampler; `sampleSequence()` for chains |
| `SamplingEngine.h / .tpp` | Owns `SamplingHandler<T>` vector; `run()` generates a `SampledScenario`; `read_parameters()` file parser |
| `SamplingHandlerFactory.h / .tpp / .cpp` | Self-registering abstract factory for sampling handlers; `build<T>(Topology)` |
| `EventFactory.h / .cpp` | Self-registering abstract factory for `Rec_Event` subclasses |
| `Scenario.h` | `SampledEvent` and `SampledScenario` value types |
| `LegacyBridge.h / .tpp` | Bidirectional adapter: `Topology ⇄ Model_Parms`, `Engine ⇄ Model_marginals` |
| `CMakeLists.txt` | Shared library target `igor::Model`, depends on `igor::Core` and `igor::Math` |

---

## Data Flow

### Inference Path (Legacy-Backed)

```
model_parms.txt ──► Model_Parms ──► extract_event_descriptors()
                                          │
                                          ▼
                                    InferenceEngine<long double>
                                          │
model_marginals.txt ──► Model_marginals ──┘  import_from_legacy()
                                          │
                                    [ E-step: reset / accumulate ]
                                    [ M-step: maximize_likelihood ]
                                          │
                                    export_to_legacy()
                                          │
                                          ▼
                                    Model_marginals ──► model_marginals.txt
```

### Generation Path (Modern)

```
model_parms.txt ──► read_topology() ──► Topology (shared)
                                             │
                                   SamplingEngine<double>(topology)
                                             │
                                   ┌─────────┴──────────┐
                                   │  Two loading paths  │
                                   │                     │
                          read_parameters()     import_from_legacy()
                          (direct file parse)   (via Model_marginals)
                                   │                     │
                                   └─────────┬───────────┘
                                             │
                                      precomputeCDF()
                                             │
                                      engine.run(rng)
                                             │
                                             ▼
                                      SampledScenario
```

---

## HasUid Protocol

The `HasUid` concept is the **structural contract** that unifies nodes across
all layers. Any type exposing `uid() → index_type` and `setUid(index_type)`
can be stored in a `Navigator` and indexed by a `Topology`-assigned UID.

Current satisfying types:

| Type | Where UID is assigned |
|---|---|
| `Rec_Event` | `Topology::addEvent()` |
| `MarginalHandler<T>` | `InferenceEngine::register_handler()` |
| `SamplingHandler<T>` | `SamplingHandlerFactory::build()` |

This protocol enables the same `Navigator` template to traverse parents and
children regardless of whether the underlying nodes are raw events, inference
handlers, or sampling handlers.

---

## Build

The module produces a shared library `igorModel` (CMake target `igor::Model`)
that links publicly against `igor::Core` and `igor::Math`.

```
igor::Core ◄── igor::Model ──► igor::Math
```
