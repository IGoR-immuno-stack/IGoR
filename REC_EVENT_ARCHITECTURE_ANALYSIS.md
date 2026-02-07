# Rec_Event Architecture Analysis

**Date:** 6 February 2026
**Status:** Design Review
**Related:** EventFactory implementation, Core module refactoring

---

## Executive Summary

The `Rec_Event` class hierarchy exhibits a **three-dimensional architecture** where each class simultaneously handles:
1. **Event Data** - domain model (what the event represents)
2. **Inference Process** - computing probabilities from observed sequences
3. **Generation Process** - sampling sequences from probability distributions

This analysis explores the architectural implications of separating these concerns and proposes three strategic options for refactoring.

---

## Current Architecture: Three Dimensions in One Hierarchy

### Dimension 1: Event Data (What the event IS)

**Responsibility:** Describe what the recombination event is and what values it can take.

```cpp
// Identity and realizations
Event_type type;                                    // GeneChoice_t, Deletion_t, etc.
Gene_class event_class;                            // V_gene, D_gene, J_gene
Seq_side event_side;                               // Five_prime, Three_prime
std::unordered_map<std::string, Event_realization> event_realizations;
std::string nickname;
int priority;

// Data access methods
Event_type get_type() const;
const unordered_map<string, Event_realization> get_realizations_map() const;
std::shared_ptr<Rec_Event> copy();
int size() const;
```

**Clean concern:** Pure domain model with no algorithmic dependencies.

---

### Dimension 2: Inference Process (Compute probabilities FROM sequences)

**Responsibility:** Given an observed sequence, compute the probability distribution over event realizations.

```cpp
// Core inference methods
void iterate(double& scenario_proba, ...);          // Explore all scenarios
void add_to_marginals(long double, Marginal_array_p&) const;  // Accumulate posteriors
void initialize_event(...);                         // Setup for inference
bool has_effect_on(Seq_type) const;                // Which sequences affected
void iterate_initialize_Len_proba(...);             // Length probability bounds

// Inference state (problematic mutable members)
mutable int base_index;                            // Where to read marginals
double new_scenario_proba;                         // Probability computation
int memory_layer_cs;                               // Memory addressing
Safety_bool_map& safety_set;                       // Offset checks
```

**Complex concern:** Heavy algorithmic logic with external dependencies (marginals, alignments, error rates).

---

### Dimension 3: Generation Process (Generate sequences FROM model)

**Responsibility:** Sample from the probability distribution to generate synthetic sequences.

```cpp
// Core generation method
std::queue<int> draw_random_realization(
    const Marginal_array_p&,
    unordered_map<Rec_Event_name, int>&,
    ...) const;

void write2txt(ofstream&);                         // Serialize for generation

// Uses:
// - Event realizations (data dimension)
// - Model marginals (probability distribution)
```

**Focused concern:** Sampling algorithm with statistical dependencies.

---

## Problem: Single Hierarchy for Three Dimensions

### Current Implementation

```
                    Rec_Event (abstract)
                         |
        +----------------+----------------+----------------+
        |                |                |                |
   Gene_choice       Deletion        Insertion      Dinucl_markov

   Each class implements ALL THREE dimensions:
   - iterate() differently          (Inference)
   - draw_random_realization()      (Generation)
   - add_to_marginals()             (Inference)
   - Stores data + algorithm state  (Mixed)
```

### Issues

1. **Single Responsibility Violation** - each class has 3 distinct jobs
2. **Tight Coupling** - changes in one dimension affect others
3. **Testing Complexity** - cannot test dimensions independently
4. **Code Duplication** - similar patterns across inference/generation
5. **God Classes** - 40+ members in Gene_choice, mixing data and algorithm state

---

## Strategic Options for Refactoring

### Option 1: Three Separate Hierarchies (Recommended)

**Core Principle:** Independent dimensions get independent polymorphism.

```cpp
// DIMENSION 1: Event Data (domain model)
// Minimal or NO hierarchy - data is just data
struct EventData {
    Event_type type;
    Gene_class gene_class;
    Seq_side side;
    std::unordered_map<std::string, EventRealization> realizations;
    std::string nickname;
};

// No inheritance needed! Data differences are in:
// - Type enum (GeneChoice_t vs Deletion_t)
// - Realization content (strings vs ints)
// - Not in behavior


// DIMENSION 2: Inference Strategy Hierarchy
class InferenceStrategy {
public:
    virtual void iterate(
        const EventData& event,
        const Sequence& seq,
        InferenceContext& ctx) = 0;
    virtual void accumulate_posterior(
        const EventData& event,
        double prob,
        Marginals& m) = 0;
    virtual ~InferenceStrategy() = default;
};

class GeneChoiceInference : public InferenceStrategy {
    void iterate(...) override {
        // Gene choice specific logic with alignments
    }
};

class DeletionInference : public InferenceStrategy {
    void iterate(...) override {
        // Deletion specific logic with offset calculations
    }
};

class InsertionInference : public InferenceStrategy {
    void iterate(...) override {
        // Insertion specific logic
    }
};

class DinuclMarkovInference : public InferenceStrategy {
    void iterate(...) override {
        // Markov chain specific logic
    }
};


// DIMENSION 3: Generation Strategy Hierarchy
class GenerationStrategy {
public:
    virtual int draw_realization(
        const EventData& event,
        const Marginals& marginals,
        RNG& rng) const = 0;
    virtual ~GenerationStrategy() = default;
};

class GeneChoiceGeneration : public GenerationStrategy {
    int draw_realization(...) const override {
        // Sample from gene choice distribution
    }
};

class DeletionGeneration : public GenerationStrategy {
    int draw_realization(...) const override {
        // Sample deletion length
    }
};

// ... etc.


// Coordinator (no longer needs to be polymorphic!)
class RecombinationEvent {
    EventData data_;
    std::unique_ptr<InferenceStrategy> inference_;
    std::unique_ptr<GenerationStrategy> generation_;

public:
    // Factory creates with correct strategies
    static std::unique_ptr<RecombinationEvent> create(Event_type type) {
        EventData data{type, ...};

        auto inference = create_inference_strategy(type);
        auto generation = create_generation_strategy(type);

        return std::make_unique<RecombinationEvent>(
            std::move(data),
            std::move(inference),
            std::move(generation)
        );
    }

    // Polymorphism delegated to strategies
    void infer(const Sequence& seq, InferenceContext& ctx) {
        inference_->iterate(data_, seq, ctx);
    }

    int generate(const Marginals& m, RNG& rng) {
        return generation_->draw_realization(data_, m, rng);
    }
};
```

#### Impact Summary
- **1 hierarchy → 3 hierarchies** (one per dimension)
- **Event classes eliminated** - replaced by strategy combinations
- **Polymorphism is explicit** - clear what varies and why
- **Lower coupling** - strategies don't know about each other

#### Benefits
- **Independent testing** - test data, inference, generation separately
- **Clear dependencies** - EventData has no algorithm dependencies
- **Reusability** - mix and match strategies
- **Easier maintenance** - changes isolated to relevant strategy

---

### Option 2: Hybrid (Keep Event Classes, Add Strategies)

**Conservative approach:** Maintain existing classes but delegate to strategies.

```cpp
// Base class becomes thinner
class Rec_Event {
protected:
    EventData data_;
    std::unique_ptr<InferenceStrategy> inference_strategy_;
    std::unique_ptr<GenerationStrategy> generation_strategy_;

public:
    // Template method pattern
    void iterate(InferenceContext& ctx) final {
        inference_strategy_->iterate(data_, ctx);
    }

    int draw_random_realization(GenerationContext& ctx) final {
        return generation_strategy_->draw_realization(data_, ctx);
    }

    // Subclasses configure strategies
    virtual void configure_strategies() = 0;
};

class Gene_choice : public Rec_Event {
    void configure_strategies() override {
        inference_strategy_ = std::make_unique<GeneChoiceInference>();
        generation_strategy_ = std::make_unique<GeneChoiceGeneration>();
    }
};

class Deletion : public Rec_Event {
    void configure_strategies() override {
        inference_strategy_ = std::make_unique<DeletionInference>();
        generation_strategy_ = std::make_unique<DeletionGeneration>();
    }
};
```

#### Impact Summary
- **1 hierarchy + 2 new hierarchies** (backward compatible)
- **Existing classes preserved** - easier migration
- **Strategies encapsulate algorithms** - can swap implementations
- **Still some coupling** - event classes know about strategies

#### Benefits
- **Backward compatible** - existing code continues to work
- **Gradual migration** - can refactor incrementally
- **Lower risk** - proven pattern (Strategy + Template Method)

---

### Option 3: Data-Driven (No Hierarchies)

**Most radical:** Eliminate polymorphism entirely, use data + function dispatch.

```cpp
// Pure data structure (no virtuals)
struct EventData {
    Event_type type;
    Gene_class gene_class;
    Seq_side side;
    std::unordered_map<std::string, EventRealization> realizations;
    std::string nickname;
};

// Strategy registries (similar to your EventFactory pattern)
class InferenceEngine {
    std::unordered_map<Event_type,
        std::function<void(const EventData&, InferenceContext&)>> strategies_;

public:
    void register_strategy(Event_type type, auto strategy) {
        strategies_[type] = strategy;
    }

    void infer(const EventData& event, InferenceContext& ctx) {
        strategies_.at(event.type)(event, ctx);
    }
};

class GenerationEngine {
    std::unordered_map<Event_type,
        std::function<int(const EventData&, const Marginals&, RNG&)>> strategies_;

public:
    int generate(const EventData& event, const Marginals& m, RNG& rng) {
        return strategies_.at(event.type)(event, m, rng);
    }
};

// Usage (similar to your factory registration)
namespace {
    void register_gene_choice_strategies() {
        inference_engine.register_strategy(
            GeneChoice_t,
            [](const EventData& e, InferenceContext& ctx) {
                // Gene choice inference logic
            }
        );

        generation_engine.register_strategy(
            GeneChoice_t,
            [](const EventData& e, const Marginals& m, RNG& rng) {
                // Gene choice generation logic
                return sample_gene(e, m, rng);
            }
        );
    }
}
```

#### Impact Summary
- **No hierarchies at all** - pure data + function dispatch
- **Maximum flexibility** - strategies are just functions
- **Similar to EventFactory pattern** - consistent with existing approach
- **Potential runtime overhead** - std::function vs virtual (usually negligible)

#### Benefits
- **Simplest design** - no inheritance complexity
- **Runtime flexibility** - can change strategies dynamically
- **Consistent with EventFactory** - same registration pattern
- **Easy to extend** - just register new functions

---

## Comparison Matrix

| Aspect | Current | Option 1 (3 Hierarchies) | Option 2 (Hybrid) | Option 3 (Data-Driven) |
|--------|---------|--------------------------|-------------------|------------------------|
| **Number of hierarchies** | 1 | 3 | 3 | 0 |
| **Polymorphism** | Event classes | Strategy classes | Both | Function dispatch |
| **Backward compatibility** | N/A | Breaking | High | Breaking |
| **Code migration effort** | N/A | High | Low-Medium | High |
| **Testing complexity** | High | Low | Medium | Low |
| **Runtime flexibility** | Low | Medium | Medium | High |
| **Type safety** | Compile-time | Compile-time | Compile-time | Runtime checks |
| **Coupling** | High | Low | Medium | Low |
| **Cognitive load** | High | Medium | Medium | Low |

---

## Testing Impact

### Current: Monolithic Testing

```cpp
// Must test all dimensions together
TEST_CASE("Gene_choice works correctly") {
    Gene_choice gc(V_gene);

    // Setup massive context for inference
    InferenceState state;
    Marginals marginals;
    ErrorRate error_rate;
    // ... 14 parameters

    gc.iterate(scenario_proba, downstream_map, seq, int_seq,
               index_map, offset_map, next_event, marginals,
               model_params, alignments, constructed_seqs,
               seq_offsets, error_rate, counters, events_map,
               safety_set, mismatches, max_prob, threshold);

    // What are we actually testing?
}
```

### With Separation: Focused Testing

```cpp
// Test data independently
TEST_CASE("EventData stores realizations correctly") {
    EventData gc_data(GeneChoice_t, V_gene);
    gc_data.add_realization("TRBV1", "ATGC");
    REQUIRE(gc_data.size() == 1);
}

// Test inference without generation
TEST_CASE("GeneChoiceInference computes correct probabilities") {
    EventData data = create_test_gene_choice();
    GeneChoiceInference inference;
    InferenceContext ctx = create_minimal_context();  // Much simpler!

    inference.iterate(data, "ATGC", ctx);

    REQUIRE(ctx.probability == expected_prob);
}

// Test generation without inference
TEST_CASE("GeneChoiceGeneration samples correctly") {
    EventData data = create_test_gene_choice();
    GeneChoiceGeneration gen;

    auto realization = gen.draw_realization(data, test_marginals, rng);

    REQUIRE(is_valid_realization(realization));
}
```

---

## Key Insights from Code Analysis

### What Actually Varies Between Event Types?

```cpp
// Gene_choice::iterate() - needs alignments, complex offset logic
// Deletion::iterate() - simple offset math
// Insertion::iterate() - length sampling
// Dinucl_markov::iterate() - Markov chain transitions

// Gene_choice::draw_random_realization() - sample from genomic templates
// Deletion::draw_random_realization() - sample deletion length
// Insertion::draw_random_realization() - sample insertion length
// Dinucl_markov::draw_random_realization() - sample dinucleotide sequence
```

**Critical Observation:** The algorithms are **completely different**, not variations of a common algorithm.

This suggests:
- ❌ **Weak justification for shared base class** - no code reuse
- ✅ **Strong justification for separate strategies** - fundamentally different approaches

---

## Recommended Migration Path

Given that the codebase already has:
- ✅ **EventFactory pattern** - data-driven creation
- ✅ **Event_type enum** - runtime type identification
- ✅ **Registrar template** - automatic registration

**Recommendation: Option 3 (Data-Driven) with gradual migration**

### Phase 1: Add Strategy Registries (Non-Breaking)

```cpp
// New inference engine (co-exists with current code)
namespace igor::inference {
    class Engine {
        using InferFn = std::function<void(const EventData&, Context&)>;
        std::unordered_map<Event_type, InferFn> strategies_;
    public:
        void register_strategy(Event_type type, InferFn fn);
        void infer(const EventData& event, Context& ctx);
    };

    // Use same Registrar pattern!
    template<Event_type Type>
    struct StrategyRegistrar {
        StrategyRegistrar(Engine::InferFn fn) {
            inference_engine.register_strategy(Type, fn);
        }
    };
}

// In Genechoice.cpp - familiar pattern!
namespace {
    static igor::inference::StrategyRegistrar<GeneChoice_t>
        gene_choice_inference(gene_choice_inference_impl);
}
```

### Phase 2: Extract EventData (Non-Breaking)

```cpp
// Add accessor to Rec_Event
class Rec_Event {
public:
    EventData get_data() const {
        return {type, event_class, event_side, event_realizations, nickname};
    }

    // Old virtual methods still work
    virtual void iterate(...) = 0;
};
```

### Phase 3: Dual Implementation (Validation)

```cpp
// Call both old and new, compare results
void iterate(...) override {
    // Old implementation (for verification)
    auto old_result = iterate_legacy(...);

    // New implementation
    EventData data = get_data();
    InferenceContext ctx = ...;
    inference_engine.infer(data, ctx);

    assert(results_match(old_result, ctx));
}
```

### Phase 4: Remove Hierarchy (Breaking, v2.0)

- Replace `std::shared_ptr<Rec_Event>` with `EventData`
- Remove virtual methods
- Clients call engines directly

---

## Impact Summary

### Immediate Changes

**Option 1 or 3:**
- **4 concrete classes → 0 concrete classes** (or just EventData struct)
- **1 abstract base → 0-2 abstract bases** (strategy interfaces)
- **Virtual dispatch → Strategy dispatch** (function or virtual)

**Option 2:**
- **4 concrete classes → 4 concrete classes** (preserved)
- **1 abstract base → 1 abstract base + 2 strategy interfaces**
- **Virtual dispatch → Virtual + Strategy dispatch**

### Migration Effort

- **Low** if using Option 2 (Hybrid) - backward compatible
- **Medium** if using Option 3 (Data-Driven) - gradual migration possible
- **High** if using Option 1 (Clean Slate) - rewrite from scratch

### Long-term Benefits

✅ **Independent evolution** of data, inference, and generation
✅ **Easier testing** - each dimension tested separately
✅ **Better performance** - no virtual calls if using function dispatch
✅ **More flexible** - mix and match strategies at runtime
✅ **Lower coupling** - dimensions don't depend on each other
✅ **Clearer code** - each component has single responsibility

---

## Architectural Principles Validated

This separation aligns with:

- ✅ **Single Responsibility Principle** - each dimension has one job
- ✅ **Strategy Pattern** - algorithms are pluggable
- ✅ **Open/Closed Principle** - extend behavior without modifying data
- ✅ **Dependency Inversion** - algorithms depend on data abstraction
- ✅ **Interface Segregation** - clients depend only on what they need

---

## Next Steps

1. **Review and discuss** this analysis with team
2. **Choose strategic option** based on constraints and priorities
3. **Create prototype** for chosen approach
4. **Measure impact** on build time, runtime, test coverage
5. **Plan incremental migration** if Option 2 or 3 chosen
6. **Update architecture documentation** with chosen design

---

## References

- **Related Files:**
  - `src/igor/Core/Rec_Event.h` - Current base class
  - `src/igor/Core/Genechoice.h` - Largest derived class (40+ members)
  - `src/igor/Core/EventFactory.h` - Factory pattern already in use
  - `EVENT_FACTORY_PLAN.md` - Factory implementation details

- **Design Patterns:**
  - Strategy Pattern - algorithm encapsulation
  - Template Method Pattern - workflow control
  - Registry Pattern - function dispatch (already used in EventFactory)

- **Related Work:**
  - EventFactory implementation (completed 5 Feb 2026)
  - ModelIO Phase 4 integration (uses factory for event creation)

---

**Document Status:** Draft for review
**Last Updated:** 6 February 2026
