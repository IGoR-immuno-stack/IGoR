# Model_marginals Refactoring Design

## Executive Summary

Replace the monolithic `Model_marginals` class (a single flattened `long double[]` array with manual offset arithmetic) with an **event-centric handler architecture** where each recombination event owns its parameter and accumulator `Tensor` objects.

The new class, `InferenceEngine`, lives in `igor::model` namespace and operates alongside the legacy `Model_marginals` during migration.

**Prototype branch**: `feature/TensorLinalg` (based on `develop`).
Tandem D support (`feature/refactoring`) deferred to post-prototype integration.

---

## Current Problems

### 1. Single Flattened Array

```cpp
// Current: one giant array, manual offset management
std::unique_ptr<long double[]> marginal_array_smart_p;

// Which event does index 137 belong to? Impossible to know without calculation.
marginal_array[137] += scenario_probability;
```

**Issues**: Fragile indexing, no bounds checking, poor cache locality, error-prone offset arithmetic.

### 2. Forced Dimensional Uniformity

All events are forced into the same flat 1D layout, even when their natural structure is different:

| Event Type | Natural Shape | Current Layout |
|:---|:---|:---|
| Gene Choice (V/D/J) | `[n_genes]` | ✅ 1D — fits |
| Deletion | `[n_deletions]` | ✅ 1D — fits |
| Insertion length | `[max_insertions]` | ✅ 1D — fits |
| Dinucleotide Markov | `[4, 4]` transition matrix | ❌ Flattened to `[16]` |
| Error rate | `[1]` parameter, `[2]` statistics | ❌ Forced into same-size slot |

### 3. Scattered Event Logic

Event-specific logic (accumulation, normalization, update) is spread across `Model_marginals`, `Rec_Event::iterate()`, and `Rec_Event::add_to_marginals()`. Adding or modifying an event type requires changes in 5+ places.

### 4. Type Erasure

```cpp
// These are completely different but have the same type:
typedef std::unique_ptr<long double[]> Marginal_array_p;

Marginal_array_p v_gene_probs;     // 1D: [n_v_genes]
Marginal_array_p markov_matrix;    // 2D: [4,4] flattened to [16]
```

No compile-time safety. No dimensional information.

---

## Proposed Architecture: Handler-Based

### Design Principles

1. **Each event is a self-contained handler** that owns its tensors and implements its EM logic.
2. **No `dynamic_cast`** — each handler knows its own data type statically.
3. **No Data-Process-View separation** — intentionally avoided for simplicity. Data and logic co-located in handlers.
4. **Backward-compatible I/O** — can read/write legacy `.txt` marginal files.

### Event Inventory (`develop` baseline)

On `develop`, the `events_map` key is `tuple<Event_type, Gene_class, Seq_side>`. The concrete `Rec_Event` subclasses are:

| Legacy Class | Event Type Enum | Handler Type | Tensor Shape |
|:---|:---|:---|:---|
| `Gene_choice` | `GeneChoice_t` | `CategoricalHandler` | `[n_genes]` |
| `Deletion` | `Deletion_t` | `CategoricalHandler` | `[n_deletions]` |
| `Insertion` | `Insertion_t` | `CategoricalHandler` | `[max_insertions]` |
| `Dinucl_markov` | `Dinuclmarkov_t` | `MarkovHandler` | `[4, 4]` |
| Error models | N/A (`Error_rate` hierarchy) | Out of scope | event-specific |

**Key observation**: Most events (Gene_choice, Deletion, Insertion) are **categorical distributions** and share identical EM logic (normalize to sum to 1). Only `Dinucl_markov` needs a specialized handler.

### Class Hierarchy

```cpp
namespace igor::model {

// ─── Bridge to Legacy Code ─────────────────────────────────────────────

/// Describes an event extracted from Model_Parms + Rec_Event
struct EventDescriptor {
    Rec_Event_name name;        // e.g. "v_choice"
    Event_type type;            // GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t
    Gene_class gene_class;      // V_gene, D_gene, J_gene, VD_genes, VJ_genes, DJ_genes
    Seq_side side;              // Five_prime, Three_prime, Undefined_side
    std::vector<size_t> shape;  // Tensor dimensions derived from Rec_Event::size()
};

/// Build descriptors from existing Model_Parms
std::vector<EventDescriptor> extract_event_descriptors(const Model_Parms& parms);

// ─── Handler Interface ─────────────────────────────────────────────────

class MarginalHandler {
public:
    virtual ~MarginalHandler() = default;

    // Identity
    const Rec_Event_name& name() const { return name_; }

    // Tensor access
    virtual const math::Tensor<double>& parameters() const = 0;
    virtual math::Tensor<double>& parameters() = 0;
    virtual const math::Tensor<double>& accumulator() const = 0;
    virtual math::Tensor<double>& accumulator() = 0;

    // EM operations
    virtual void reset_accumulator() = 0;
    virtual void maximize_likelihood() = 0;

    // I/O (backward compatible with legacy txt format)
    virtual void write_parameters(std::ostream& out) const = 0;
    virtual void read_parameters(std::istream& in) = 0;

protected:
    explicit MarginalHandler(Rec_Event_name name) : name_(std::move(name)) {}
    Rec_Event_name name_;
};

// ─── Concrete Handlers ─────────────────────────────────────────────────

/// Handles Gene_choice, Deletion, Insertion
/// (any event whose M-step is: normalize accumulator to sum to 1)
class CategoricalHandler : public MarginalHandler {
public:
    CategoricalHandler(Rec_Event_name name, size_t n_realizations)
        : MarginalHandler(std::move(name))
        , parameters_({n_realizations})
        , accumulator_({n_realizations})
    {
        double uniform = 1.0 / static_cast<double>(n_realizations);
        std::fill(parameters_.begin(), parameters_.end(), uniform);
        std::fill(accumulator_.begin(), accumulator_.end(), 0.0);
    }

    const math::Tensor<double>& parameters() const override { return parameters_; }
    math::Tensor<double>& parameters() override { return parameters_; }
    const math::Tensor<double>& accumulator() const override { return accumulator_; }
    math::Tensor<double>& accumulator() override { return accumulator_; }

    void reset_accumulator() override {
        std::fill(accumulator_.begin(), accumulator_.end(), 0.0);
    }

    void maximize_likelihood() override {
        double total = math::linalg::sum(accumulator_);
        if (total > 1e-15) {
            for (size_t i = 0; i < parameters_.size(); ++i) {
                parameters_.data()[i] = accumulator_.data()[i] / total;
            }
        }
    }

    // Accumulate a single realization
    void accumulate(size_t realization_index, double probability) {
        accumulator_(realization_index) += probability;
    }

    void write_parameters(std::ostream& out) const override { /* ... */ }
    void read_parameters(std::istream& in) override { /* ... */ }

private:
    math::Tensor<double> parameters_;   // [n_realizations]
    math::Tensor<double> accumulator_;  // [n_realizations]
};

/// Handles Dinucl_markov (insertion nucleotide identity)
/// M-step: row-wise normalization of transition matrix
class MarkovHandler : public MarginalHandler {
public:
    MarkovHandler(Rec_Event_name name, size_t n_states)
        : MarginalHandler(std::move(name))
        , n_states_(n_states)
        , parameters_({n_states, n_states})
        , accumulator_({n_states, n_states})
    {
        double uniform = 1.0 / static_cast<double>(n_states);
        std::fill(parameters_.begin(), parameters_.end(), uniform);
        std::fill(accumulator_.begin(), accumulator_.end(), 0.0);
    }

    const math::Tensor<double>& parameters() const override { return parameters_; }
    math::Tensor<double>& parameters() override { return parameters_; }
    const math::Tensor<double>& accumulator() const override { return accumulator_; }
    math::Tensor<double>& accumulator() override { return accumulator_; }

    void reset_accumulator() override {
        std::fill(accumulator_.begin(), accumulator_.end(), 0.0);
    }

    void maximize_likelihood() override {
        // Row-wise normalization: P(to|from) = count[from,to] / sum_to(count[from,to])
        for (size_t from = 0; from < n_states_; ++from) {
            double row_sum = 0.0;
            for (size_t to = 0; to < n_states_; ++to) {
                row_sum += accumulator_(from, to);
            }
            if (row_sum > 1e-15) {
                for (size_t to = 0; to < n_states_; ++to) {
                    parameters_(from, to) = accumulator_(from, to) / row_sum;
                }
            }
        }
    }

    // Accumulate a transition observation
    void accumulate(size_t from_state, size_t to_state, double probability) {
        accumulator_(from_state, to_state) += probability;
    }

    void write_parameters(std::ostream& out) const override { /* ... */ }
    void read_parameters(std::istream& in) override { /* ... */ }

private:
    size_t n_states_;
    math::Tensor<double> parameters_;   // [n_states, n_states]
    math::Tensor<double> accumulator_;  // [n_states, n_states]
};

// ─── Engine ────────────────────────────────────────────────────────────

class InferenceEngine {
public:
    /// Construct from existing model parameters
    explicit InferenceEngine(const Model_Parms& parms);

    /// Register a handler for a named event
    void register_handler(Rec_Event_name name,
                         std::unique_ptr<MarginalHandler> handler);

    /// Access handler by event name
    MarginalHandler& handler(const Rec_Event_name& name);
    const MarginalHandler& handler(const Rec_Event_name& name) const;

    /// EM: Reset all accumulators
    void reset_accumulators();

    /// EM: Update all parameters (M-step)
    void update_parameters();

    /// I/O: Read/write marginals in legacy text format
    void write_marginals(const std::string& filename) const;
    void read_marginals(const std::string& filename);

    /// Convert to/from legacy flat array for validation
    void export_to_legacy(Model_marginals& legacy, const Model_Parms& parms) const;
    void import_from_legacy(const Model_marginals& legacy, const Model_Parms& parms);

    /// Ordered iteration (follows Model_Parms topological order)
    const std::vector<Rec_Event_name>& event_names() const { return event_order_; }

private:
    std::unordered_map<Rec_Event_name, std::unique_ptr<MarginalHandler>> handlers_;
    std::vector<Rec_Event_name> event_order_;
};

} // namespace igor::model
```

### Integration with Legacy Code

#### Building from `Model_Parms`

```cpp
std::vector<EventDescriptor> extract_event_descriptors(const Model_Parms& parms) {
    std::vector<EventDescriptor> descriptors;

    for (const auto& event_ptr : parms.get_event_list()) {
        EventDescriptor desc;
        desc.name = event_ptr->get_name();
        desc.type = event_ptr->get_type();
        desc.gene_class = event_ptr->get_class();
        desc.side = event_ptr->get_side();

        // Derive shape from event type:
        switch (desc.type) {
            case GeneChoice_t:
            case Deletion_t:
            case Insertion_t:
                desc.shape = {static_cast<size_t>(event_ptr->size())};
                break;
            case Dinuclmarkov_t:
                desc.shape = {4, 4};  // Always 4 nucleotides
                break;
            default:
                desc.shape = {static_cast<size_t>(event_ptr->size())};
                break;
        }

        descriptors.push_back(std::move(desc));
    }

    return descriptors;
}
```

#### Constructing `InferenceEngine` from `Model_Parms`

```cpp
InferenceEngine::InferenceEngine(const Model_Parms& parms) {
    auto descriptors = extract_event_descriptors(parms);

    for (const auto& desc : descriptors) {
        std::unique_ptr<MarginalHandler> handler;

        switch (desc.type) {
            case GeneChoice_t:
            case Deletion_t:
            case Insertion_t:
                handler = std::make_unique<CategoricalHandler>(
                    desc.name, desc.shape[0]);
                break;
            case Dinuclmarkov_t:
                handler = std::make_unique<MarkovHandler>(
                    desc.name, desc.shape[0]);
                break;
            default:
                handler = std::make_unique<CategoricalHandler>(
                    desc.name, desc.shape[0]);
                break;
        }

        register_handler(desc.name, std::move(handler));
    }
}
```

#### Validation Against Legacy

```cpp
/// Export new tensors back to legacy flat array for numerical comparison
void InferenceEngine::export_to_legacy(Model_marginals& legacy,
                                        const Model_Parms& parms) const {
    auto index_map = legacy.get_index_map(parms);

    for (const auto& event_name : event_order_) {
        const auto& h = *handlers_.at(event_name);
        int offset = index_map.at(event_name);
        const auto& params = h.parameters();

        for (size_t i = 0; i < params.size(); ++i) {
            legacy.marginal_array_smart_p[offset + i] =
                static_cast<long double>(params.data()[i]);
        }
    }
}

/// Import from legacy flat array (e.g., after loading from file)
void InferenceEngine::import_from_legacy(const Model_marginals& legacy,
                                          const Model_Parms& parms) {
    auto index_map = legacy.get_index_map(parms);

    for (const auto& event_name : event_order_) {
        auto& h = *handlers_.at(event_name);
        int offset = index_map.at(event_name);
        auto& params = h.parameters();

        for (size_t i = 0; i < params.size(); ++i) {
            params.data()[i] =
                static_cast<double>(legacy.marginal_array_smart_p[offset + i]);
        }
    }
}
```

---

## Migration Plan

### Prototype Phase (on `feature/TensorLinalg`, based on `develop`)

Objective: Prove that the handler architecture works end-to-end with the standard VDJ model.

**Step 1: Handlers** (target: compile + unit tests pass)
1. Add `src/igor/Model/MarginalHandler.h` — base class + `EventDescriptor`.
2. Add `src/igor/Model/CategoricalHandler.h` — categorical distribution handler.
3. Add `src/igor/Model/MarkovHandler.h` — transition matrix handler.
4. Add `tst/igor/Model/test_Handlers.cpp` — unit tests for normalize, accumulate, reset.

**Step 2: Engine** (target: integration test passes)
1. Add `src/igor/Model/InferenceEngine.h` + `.cpp`.
2. Add `tst/igor/Model/test_InferenceEngine.cpp` — construct from `Model_Parms`, round-trip test.

**Step 3: Legacy Bridge** (target: numerical validation passes)
1. Implement `import_from_legacy()` / `export_to_legacy()`.
2. Load a real IGoR model (from `IGoR-models` repo), import into `InferenceEngine`, export back, compare.

### Post-Prototype

Once the prototype validates on `develop`:

1. **Merge into `feature/refactoring`**: Adapt `EventDescriptor` for tandem D (`int` gene class, `SequenceTypeRegistry`).
2. **Phase 2: Side-by-Side EM**: Run parallel EM iterations comparing both implementations.
3. **Phase 3: Integration**: Refactor `Rec_Event::add_to_marginals()` to accept handler references.
4. **Phase 4: Cleanup**: Remove `Model_marginals`, finalize `InferenceEngine` as canonical.

---

## Testing Strategy

### Unit Tests (per Handler)

```cpp
TEST_CASE("CategoricalHandler normalize", "[model][handler]") {
    CategoricalHandler handler("v_gene", 3);

    handler.accumulate(0, 10.0);
    handler.accumulate(1, 20.0);
    handler.accumulate(2, 30.0);

    handler.maximize_likelihood();

    REQUIRE_THAT(handler.parameters()(0), WithinRel(1.0/6.0, 1e-10));
    REQUIRE_THAT(handler.parameters()(1), WithinRel(2.0/6.0, 1e-10));
    REQUIRE_THAT(handler.parameters()(2), WithinRel(3.0/6.0, 1e-10));
}

TEST_CASE("MarkovHandler row normalize", "[model][handler]") {
    MarkovHandler handler("vd_dinucl", 4);

    handler.accumulate(0, 0, 10.0);  // A → A
    handler.accumulate(0, 1, 30.0);  // A → C
    handler.maximize_likelihood();

    // Row 0 should sum to 1
    REQUIRE_THAT(handler.parameters()(0, 0), WithinRel(0.25, 1e-10));
    REQUIRE_THAT(handler.parameters()(0, 1), WithinRel(0.75, 1e-10));
}
```

### Integration Tests

```cpp
TEST_CASE("InferenceEngine round-trip with legacy", "[model][integration]") {
    // Load legacy model
    Model_Parms parms;
    parms.read_model_parms("test_model_parms.txt");
    Model_marginals legacy(parms);
    legacy.txt2marginals("test_marginals.txt", parms);

    // Import into engine
    InferenceEngine engine(parms);
    engine.import_from_legacy(legacy, parms);

    // Export back
    Model_marginals roundtrip(parms);
    engine.export_to_legacy(roundtrip, parms);

    // Compare
    for (size_t i = 0; i < legacy.get_length(); ++i) {
        REQUIRE_THAT(static_cast<double>(roundtrip.marginal_array_smart_p[i]),
                     WithinRel(static_cast<double>(legacy.marginal_array_smart_p[i]),
                               1e-10));
    }
}
```

### Property Tests

- All categorical parameters sum to 1.0 after `maximize_likelihood()`.
- All Markov rows sum to 1.0 after `maximize_likelihood()`.
- Zero accumulator → parameters unchanged.
- Import/export round-trip preserves all values.
- Event count matches `Model_Parms::get_event_list().size()`.

---

## Open Questions

1. **Error models**: `Singleerrorrate`, `Hypermutation_global_errorrate`, `Hypermutation_full_Nmer_errorrate` are managed by `Model_Parms::error_rate`. Should they become handlers or remain as-is?
2. **Conditional marginals**: The Bayesian network encodes parent-child dependencies between events. Currently `get_offsets_map` and `get_inverse_offset_map` handle this. How should handlers represent conditional dependencies?
3. **Thread safety**: Per-thread accumulators or atomic operations for parallel inference. Deferred to post-prototype.

---

*Last Updated: 10 February 2026*
*Baseline: `develop` @ `83aa113` (via `feature/TensorLinalg`)*
*Status: Draft — Prototype Phase*
