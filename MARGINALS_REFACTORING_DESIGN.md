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
//
// Template parameter T is the scalar type for storage:
//   - long double : backward compatible with legacy Marginal_array_p
//   - double      : sufficient precision, 2x less memory, better vectorization
//

template<typename T = double>
class MarginalHandler {
public:
    using scalar_type = T;

    virtual ~MarginalHandler() = default;

    // Identity
    const Rec_Event_name& name() const { return name_; }

    // Tensor access
    virtual const math::Tensor<T>& parameters() const = 0;
    virtual math::Tensor<T>& parameters() = 0;
    virtual const math::Tensor<T>& accumulator() const = 0;
    virtual math::Tensor<T>& accumulator() = 0;

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
template<typename T = double>
class CategoricalHandler : public MarginalHandler<T> {
public:
    CategoricalHandler(Rec_Event_name name, size_t n_realizations)
        : MarginalHandler<T>(std::move(name))
        , parameters_({n_realizations})
        , accumulator_({n_realizations})
    {
        T uniform = T(1) / static_cast<T>(n_realizations);
        std::fill(parameters_.begin(), parameters_.end(), uniform);
        std::fill(accumulator_.begin(), accumulator_.end(), T(0));
    }

    const math::Tensor<T>& parameters() const override { return parameters_; }
    math::Tensor<T>& parameters() override { return parameters_; }
    const math::Tensor<T>& accumulator() const override { return accumulator_; }
    math::Tensor<T>& accumulator() override { return accumulator_; }

    void reset_accumulator() override {
        std::fill(accumulator_.begin(), accumulator_.end(), T(0));
    }

    void maximize_likelihood() override {
        T total = math::linalg::sum(accumulator_);
        if (total > T(1e-15)) {
            for (size_t i = 0; i < parameters_.size(); ++i) {
                parameters_.data()[i] = accumulator_.data()[i] / total;
            }
        }
    }

    // Accumulate a single realization
    void accumulate(size_t realization_index, T probability) {
        accumulator_(realization_index) += probability;
    }

    void write_parameters(std::ostream& out) const override { /* ... */ }
    void read_parameters(std::istream& in) override { /* ... */ }

private:
    math::Tensor<T> parameters_;   // [n_realizations]
    math::Tensor<T> accumulator_;  // [n_realizations]
};

/// Handles Dinucl_markov (insertion nucleotide identity)
/// M-step: row-wise normalization of transition matrix
template<typename T = double>
class MarkovHandler : public MarginalHandler<T> {
public:
    MarkovHandler(Rec_Event_name name, size_t n_states)
        : MarginalHandler<T>(std::move(name))
        , n_states_(n_states)
        , parameters_({n_states, n_states})
        , accumulator_({n_states, n_states})
    {
        T uniform = T(1) / static_cast<T>(n_states);
        std::fill(parameters_.begin(), parameters_.end(), uniform);
        std::fill(accumulator_.begin(), accumulator_.end(), T(0));
    }

    const math::Tensor<T>& parameters() const override { return parameters_; }
    math::Tensor<T>& parameters() override { return parameters_; }
    const math::Tensor<T>& accumulator() const override { return accumulator_; }
    math::Tensor<T>& accumulator() override { return accumulator_; }

    void reset_accumulator() override {
        std::fill(accumulator_.begin(), accumulator_.end(), T(0));
    }

    void maximize_likelihood() override {
        // Row-wise normalization: P(to|from) = count[from,to] / sum_to(count[from,to])
        for (size_t from = 0; from < n_states_; ++from) {
            T row_sum = T(0);
            for (size_t to = 0; to < n_states_; ++to) {
                row_sum += accumulator_(from, to);
            }
            if (row_sum > T(1e-15)) {
                for (size_t to = 0; to < n_states_; ++to) {
                    parameters_(from, to) = accumulator_(from, to) / row_sum;
                }
            }
        }
    }

    // Accumulate a transition observation
    void accumulate(size_t from_state, size_t to_state, T probability) {
        accumulator_(from_state, to_state) += probability;
    }

    void write_parameters(std::ostream& out) const override { /* ... */ }
    void read_parameters(std::istream& in) override { /* ... */ }

private:
    size_t n_states_;
    math::Tensor<T> parameters_;   // [n_states, n_states]
    math::Tensor<T> accumulator_;  // [n_states, n_states]
};

// ─── Engine ────────────────────────────────────────────────────────────

template<typename T = double>
class InferenceEngine {
public:
    using scalar_type = T;

    /// Construct from existing model parameters
    explicit InferenceEngine(const Model_Parms& parms);

    /// Register a handler for a named event
    void register_handler(Rec_Event_name name,
                         std::unique_ptr<MarginalHandler<T>> handler);

    /// Access handler by event name
    MarginalHandler<T>& handler(const Rec_Event_name& name);
    const MarginalHandler<T>& handler(const Rec_Event_name& name) const;

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
    std::unordered_map<Rec_Event_name, std::unique_ptr<MarginalHandler<T>>> handlers_;
    std::vector<Rec_Event_name> event_order_;
};

// ─── Type Aliases ──────────────────────────────────────────────────────

/// Legacy-compatible engine using long double (matches Marginal_array_p)
using LegacyEngine = InferenceEngine<long double>;

/// Default engine using double (better performance, sufficient precision)
using Engine = InferenceEngine<double>;

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

#### Constructing `InferenceEngine<T>` from `Model_Parms`

```cpp
template<typename T>
InferenceEngine<T>::InferenceEngine(const Model_Parms& parms) {
    auto descriptors = extract_event_descriptors(parms);

    for (const auto& desc : descriptors) {
        std::unique_ptr<MarginalHandler<T>> handler;

        switch (desc.type) {
            case GeneChoice_t:
            case Deletion_t:
            case Insertion_t:
                handler = std::make_unique<CategoricalHandler<T>>(
                    desc.name, desc.shape[0]);
                break;
            case Dinuclmarkov_t:
                handler = std::make_unique<MarkovHandler<T>>(
                    desc.name, desc.shape[0]);
                break;
            default:
                handler = std::make_unique<CategoricalHandler<T>>(
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
template<typename T>
void InferenceEngine<T>::export_to_legacy(Model_marginals& legacy,
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
template<typename T>
void InferenceEngine<T>::import_from_legacy(const Model_marginals& legacy,
                                             const Model_Parms& parms) {
    auto index_map = legacy.get_index_map(parms);

    for (const auto& event_name : event_order_) {
        auto& h = *handlers_.at(event_name);
        int offset = index_map.at(event_name);
        auto& params = h.parameters();

        for (size_t i = 0; i < params.size(); ++i) {
            params.data()[i] =
                static_cast<T>(legacy.marginal_array_smart_p[offset + i]);
        }
    }
}
```

---

## Migration Plan

### Prototype Phase (on `feature/TensorLinalg`, based on `develop`)

Objective: Prove that the handler architecture works end-to-end with the standard VDJ model.

**Step 1: Tensor API** ✅
1. ~~Compound assignment operators (`+=`, `-=`, `*=`, `/=`) for tensor and scalar.~~
2. ~~Exact equality (`==`, `!=`) and approximate comparison (`allclose`).~~
3. ~~Deep copy semantics for `Tensor` (always owning after copy).~~
4. ~~Unit tests: `tst/igor/Math/test_TensorOps.cpp` (820 assertions).~~

**Step 2: Handlers** ✅
1. ~~`src/igor/Model/MarginalHandler.h` — base class + `EventDescriptor`.~~
2. ~~`src/igor/Model/CategoricalHandler.h/.tpp` — categorical distribution handler.~~
3. ~~`src/igor/Model/MarkovHandler.h/.tpp` — transition matrix handler.~~
4. ~~`tst/igor/Model/test_Handlers.cpp` — unit tests for normalize, accumulate, reset, merge.~~

*Design note*: `MarginalHandler` has no `combine_accumulator` method — accumulator
merging is done directly via `Tensor::operator+=` at the call site (`InferenceEngine`
or user code). This keeps the handler interface minimal.

**Step 3: Engine** ✅
1. ~~`src/igor/Model/InferenceEngine.h/.tpp` — orchestrator.~~
2. ~~`tst/igor/Model/test_InferenceEngine.cpp` — construction, EM cycle, I/O, handler access.~~
3. ~~`combine_accumulators()` uses `Tensor::operator+=` directly on handler accumulator tensors.~~

**Step 4: Legacy Bridge** ✅
1. ~~`src/igor/Core/LegacyBridge.h/.tpp` — `extract_event_descriptors`, `import_from_legacy`, `export_to_legacy`.~~
2. ~~`tst/igor/Core/test_LegacyBridge.cpp` — basic tests.~~
3. ~~Load a real IGoR model (Mouse TCR beta), import into `InferenceEngine`, export back, compare.~~

### Post-Prototype

Once the prototype validates on `develop`:

**Phase 2: Side-by-Side EM Validation**

Run parallel EM iterations comparing legacy `Model_marginals` against new `InferenceEngine`
to ensure numerical equivalence before replacing the inference code.

**Validation setup**:
```cpp
// Load same model in both systems
Model_Parms parms;
parms.read_model_parms("model_parms.txt");

Model_marginals legacy(parms);
legacy.txt2marginals("initial_marginals.txt", parms);

InferenceEngine<long double> engine(parms);
import_from_legacy(engine, legacy, parms);

// Verify initial state matches
Model_marginals exported(parms);
export_to_legacy(engine, exported, parms);
assert_arrays_equal(legacy.marginal_array_smart_p.get(),
                   exported.marginal_array_smart_p.get(),
                   legacy.get_length());
```

**EM iteration comparison**:
```cpp
// Process same sequence dataset through both implementations
vector<string> sequences = load_sequences("test_data.txt");

for (int iteration = 0; iteration < max_iterations; ++iteration) {
    // LEGACY PATH
    legacy.initialize();  // Zero accumulators
    double legacy_likelihood = 0.0;

    for (const auto& seq : sequences) {
        auto scenarios = generate_scenarios(seq, parms);
        for (const auto& scenario : scenarios) {
            double prob = scenario.probability;
            legacy_likelihood += log(prob);

            // Accumulate into legacy flat array
            for (const auto& event : scenario.events) {
                event->add_to_marginals(legacy, parms, prob);
            }
        }
    }
    legacy.normalize(parms);  // M-step

    // NEW PATH
    engine.reset_accumulators();
    double engine_likelihood = 0.0;

    for (const auto& seq : sequences) {
        auto scenarios = generate_scenarios(seq, parms);
        for (const auto& scenario : scenarios) {
            double prob = scenario.probability;
            engine_likelihood += log(prob);

            // Accumulate into handlers
            for (const auto& event : scenario.events) {
                accumulate_to_engine(engine, event, prob);
            }
        }
    }
    engine.update_parameters();  // M-step

    // COMPARE
    REQUIRE_THAT(engine_likelihood, WithinRel(legacy_likelihood, 1e-12));

    // Export and compare all parameters
    export_to_legacy(engine, exported, parms);
    for (size_t i = 0; i < legacy.get_length(); ++i) {
        REQUIRE_THAT(static_cast<double>(exported.marginal_array_smart_p[i]),
                     WithinRel(static_cast<double>(legacy.marginal_array_smart_p[i]),
                               1e-12));
    }

    // Check convergence (both should converge identically)
    if (fabs(engine_likelihood - prev_likelihood) < convergence_threshold) {
        break;
    }
}
```

**Adapter for event accumulation** (temporary bridge code):
```cpp
void accumulate_to_engine(InferenceEngine<long double>& engine,
                         const Rec_Event* event,
                         double probability) {
    const auto& name = event->get_nickname();

    switch (event->get_type()) {
        case GeneChoice_t:
        case Deletion_t:
        case Insertion_t: {
            auto& handler = engine.get_handler<CategoricalHandler<long double>>(name);
            int realization = event->get_current_realization_index();
            handler.accumulator().data()[realization] += probability;
            break;
        }
        case Dinuclmarkov_t: {
            auto& handler = engine.get_handler<MarkovHandler<long double>>(name);
            // Get flat indices from legacy realization queue
            const auto& indices = event->get_realization_indices();
            for (int flat_idx : indices) {
                if (flat_idx >= 0) {
                    handler.accumulator().data()[flat_idx] += probability;
                }
            }
            break;
        }
    }
}
```

**Success criteria for Phase 2**:
- ✅ Log-likelihood trajectories identical within `1e-12` relative tolerance
- ✅ Parameter values match at each iteration within `1e-12`
- ✅ Both implementations converge in same number of iterations
- ✅ Final parameters produce identical sequence probabilities
- ✅ Validation passes on multiple model types (TRA, TRB, IGH, IGL)

**Phase 3: Integration — Refactor `Rec_Event::add_to_marginals()`**:

   **Key design point — Markov accumulation adapter**:

   In the legacy code, `Dinucl_markov::add_to_marginals()` does **not** use 2D `(from, to)`
   indexing. Instead, during `iterate_common()`, each inserted nucleotide position produces a
   **flat realization index** (`base_index + from * 4 + to`) stored in `vd_realizations_indices[i]`.
   Then `add_to_marginals()` iterates over all nucleotide positions in the insertion and
   accumulates `scenario_proba` into each corresponding flat index:

   ```cpp
   // Legacy: one scenario contributes MULTIPLE transitions (one per inserted nucleotide)
   for (size_t i = 0; i != vd_seq_size; ++i) {
       if (vd_realizations_indices[i] >= 0) {
           updated_marginals[vd_realizations_indices[i]] += scenario_proba;
       }
   }
   ```

   This is compatible with the `MarkovHandler` tensor layout because both use **row-major
   contiguous** storage: `accumulator_.data()[from * 4 + to]` accesses the same cell as
   `accumulator_(from, to)`. The adapter for the new system will be:

   ```cpp
   // New: same loop, but writing into the handler's tensor buffer directly
   auto* acc = handler.accumulator().data();
   for (size_t i = 0; i < insertion_length; ++i) {
       int from_nt = previous_nucleotide(i);
       int to_nt   = current_nucleotide(i);
       if (from_nt < 4 && to_nt < 4) {
           acc[from_nt * 4 + to_nt] += scenario_proba;
       }
   }
   ```

   The accumulation can only happen once the full scenario is explored (all nucleotide
   positions assigned). This is identical to the current behavior — `add_to_marginals()`
   is called after `iterate()` completes the scenario.

**Phase 4: Performance & Optimization**:

   After numerical validation in Phase 2, measure and optimize:

   - **Cache performance**: Handlers group related parameters → better locality than flat array
   - **Vectorization**: Contiguous `Tensor` storage enables SIMD operations
   - **Parallel accumulation**: Per-thread `InferenceEngine` with `combine_accumulators()`
   - **Memory footprint**: Compare `double` vs `long double` (2x reduction with sufficient precision)

   Expected improvements:
   - 10-50% faster EM iterations (better cache hits, vectorization)
   - 50% memory reduction if using `double` instead of `long double`
   - Linear scaling with thread count for sequence processing

**Phase 5: Cleanup & Deprecation**:

   Once Phase 2-4 complete and production validation passes:

   1. Mark `Model_marginals` as `[[deprecated]]`
   2. Route all inference through `InferenceEngine`
   3. Remove legacy flat array code and offset arithmetic
   4. Update documentation and examples
   5. Archive migration design documents

---

## Testing Strategy

### Unit Tests (per Handler, parameterized on scalar type)

```cpp
TEMPLATE_TEST_CASE("CategoricalHandler normalize", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("v_gene", 3);

    handler.accumulator()(0) = TestType(10);
    handler.accumulator()(1) = TestType(20);
    handler.accumulator()(2) = TestType(30);

    handler.maximize_likelihood();

    REQUIRE_THAT(double(handler.parameters()(0)), WithinRel(1.0/6.0, 1e-10));
    REQUIRE_THAT(double(handler.parameters()(1)), WithinRel(2.0/6.0, 1e-10));
    REQUIRE_THAT(double(handler.parameters()(2)), WithinRel(3.0/6.0, 1e-10));
}

TEMPLATE_TEST_CASE("MarkovHandler row normalize", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> handler("vd_dinucl", 4);

    handler.accumulator()(0, 0) = TestType(10);  // A → A
    handler.accumulator()(0, 1) = TestType(30);  // A → C
    handler.maximize_likelihood();

    // Row 0 should sum to 1
    REQUIRE_THAT(double(handler.parameters()(0, 0)), WithinRel(0.25, 1e-10));
    REQUIRE_THAT(double(handler.parameters()(0, 1)), WithinRel(0.75, 1e-10));
}

TEMPLATE_TEST_CASE("MarkovHandler flat accumulation (legacy-compatible)", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> handler("vd_dinucl", 4);

    // Simulate legacy-style flat accumulation: A→T = index 3 = (0*4 + 3)
    handler.accumulator().data()[0 * 4 + 3] += TestType(1.0);
    // Same as handler.accumulator()(0, 3) += TestType(1.0)

    REQUIRE(handler.accumulator()(0, 3) == TestType(1.0));
}
```

### Integration Tests

```cpp
TEST_CASE("LegacyEngine round-trip with Model_marginals", "[model][integration]") {
    // Load legacy model
    Model_Parms parms;
    parms.read_model_parms("test_model_parms.txt");
    Model_marginals legacy(parms);
    legacy.txt2marginals("test_marginals.txt", parms);

    // Import into engine (long double for exact match)
    InferenceEngine<long double> engine(extract_event_descriptors(parms));
    import_from_legacy(engine, legacy, parms);

    // Export back
    Model_marginals roundtrip(parms);
    export_to_legacy(engine, roundtrip, parms);

    // Compare — should be exact match with long double
    for (size_t i = 0; i < legacy.get_length(); ++i) {
        REQUIRE(roundtrip.marginal_array_smart_p[i] ==
                legacy.marginal_array_smart_p[i]);
    }
}

TEST_CASE("Engine (double) round-trip within tolerance", "[model][integration]") {
    Model_Parms parms;
    parms.read_model_parms("test_model_parms.txt");
    Model_marginals legacy(parms);
    legacy.txt2marginals("test_marginals.txt", parms);

    // Import into engine (double — may lose some precision)
    InferenceEngine<double> engine(extract_event_descriptors(parms));
    import_from_legacy(engine, legacy, parms);

    // Export back
    Model_marginals roundtrip(parms);
    export_to_legacy(engine, roundtrip, parms);

    // Compare — within double precision tolerance
    for (size_t i = 0; i < legacy.get_length(); ++i) {
        REQUIRE_THAT(static_cast<double>(roundtrip.marginal_array_smart_p[i]),
                     WithinRel(static_cast<double>(legacy.marginal_array_smart_p[i]),
                               1e-15));
    }
}
```

### Property Tests

- All categorical parameters sum to 1.0 after `maximize_likelihood()`.
- All Markov rows sum to 1.0 after `maximize_likelihood()`.
- Zero accumulator → parameters unchanged.
- Import/export round-trip preserves all values.
- Event count matches `Model_Parms::get_event_list().size()`.
- `InferenceEngine<long double>` round-trip is **exact** (no precision loss).
- `InferenceEngine<double>` round-trip is within `1e-15` relative tolerance.
- Flat accumulation `data()[from*4+to]` is identical to 2D accumulation `(from, to)`.

---

## Open Questions

1. **Error models**: `Singleerrorrate`, `Hypermutation_global_errorrate`, `Hypermutation_full_Nmer_errorrate` are managed by `Model_Parms::error_rate`. Should they become handlers or remain as-is?
2. **Conditional marginals**: The Bayesian network encodes parent-child dependencies between events. Currently `get_offsets_map` and `get_inverse_offset_map` handle this. How should handlers represent conditional dependencies?
3. **Thread safety**: Per-thread accumulators or atomic operations for parallel inference. Deferred to post-prototype.

---

*Last Updated: 11 February 2026*
*Baseline: `develop` @ `83aa113` (via `feature/TensorLinalg`)*
*Status: Prototype Phase — COMPLETE (All core components implemented & validated)*

