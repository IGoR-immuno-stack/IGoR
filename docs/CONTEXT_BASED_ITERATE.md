# Context-Based Iterate Pattern

**Last Updated**: April 30, 2026  
**Status**: ✅ Complete - All components migrated to context-based architecture

---

## Overview

The `iterate()` method has been refactored from a 19-parameter interface to use 5 semantic context objects. This refactoring improves:

- **Readability**: Related parameters are logically grouped into contexts
- **Maintainability**: Adding new features requires minimal signature changes
- **Testability**: Contexts can be mocked or partially initialized for isolated testing
- **Performance**: 3% faster than baseline (better memory locality)

### Before (Legacy - 19 parameters)
```cpp
void iterate(
    double& scenario_proba,
    Downstream_scenario_proba_bound_map& downstream_proba_map,
    const std::string& sequence,
    const Int_Str& int_sequence,
    Index_map& index_map,
    const std::unordered_map<Rec_Event_name, 
        std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>& offset_map,
    std::shared_ptr<Next_event_ptr>& next_event_ptr_arr,
    Marginal_array_p& updated_marginals,
    const Marginal_array_p& model_parameters_point,
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>& allowed_realizations,
    Seq_type_str_p_map& constructed_sequences,
    Seq_offsets_map& seq_offsets,
    std::shared_ptr<Error_rate>& error_rate_p,
    std::map<size_t, std::shared_ptr<Counter>>& counters_list,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, 
        std::shared_ptr<Rec_Event>>& events_map,
    Safety_bool_map& safety_set,
    Mismatch_vectors_map& mismatches_lists,
    double& seq_max_prob_scenario,
    const double& proba_threshold_factor
);
```

### After (Context-based - 5 contexts)
```cpp
void iterate(
    QuerySequenceContext& query,
    const ModelContext& model,
    ScenarioContext& scenario,
    ExplorationContext& exploration,
    AccumulationContext& accumulation
);
```

---

## Context Objects

### 1. QuerySequenceContext

**Purpose**: Provides read-only access to the input sequence being analyzed.

**Contains**:
- `sequence`: Original DNA sequence string
- `int_sequence`: Integer-encoded sequence
- `allowed_realizations`: Per-gene alignment data from alignment step

**Usage**:
```cpp
// Read the original sequence
std::string original = query.sequence;

// Access alignments for V gene
const auto& v_alignments = query.allowed_realizations.at(V_gene);
```

**Lifetime**: Created once per sequence, passed by reference through recursion.

---

### 2. ModelContext

**Purpose**: Provides **const** (read-only) access to model structure and parameters.

**Contains**:
- `model_parameters`: Marginal probability array (read-only)
- `offset_map`: Event dependency relationships
- `events_map`: All recombination events (V/D/J choice, deletions, insertions)
- `model_queue`: Event processing order

**Usage**:
```cpp
// Read event probability
double event_prob = model.model_parameters[param_index];

// Lookup event by type
auto v_choice = model.events_map.at(std::make_tuple(GeneChoice_t, V_gene, Undefined_side));

// Check event dependencies
const auto& dependencies = model.offset_map.at(this->event_name);
```

**Key Design**: **ModelContext is const** to enforce read-only model access during iteration. Mutable state belongs in other contexts.

---

### 3. ScenarioContext

**Purpose**: Tracks mutable state as scenarios are recursively explored.

**Contains**:
- `scenario_proba`: Running scenario probability
- `constructed_sequences`: Sequence fragments built by events
- `seq_offsets`: Position of each fragment in final sequence
- `mismatches_lists`: Error positions for each gene segment
- `index_map`: Parent event realization tracking
- `downstream_proba_map`: Probability pruning bounds
- `safety_set`: Gene overlap detection flags
- `scenario_error_w_proba`: Error-weighted probability (computed at leaf nodes)

**Usage**:
```cpp
// Update probability
scenario.scenario_proba *= event_probability;

// Record sequence position
scenario.seq_offsets.set_value(V_gene_seq, Five_prime, v_5p_del, memory_layer);

// Access constructed sequences
Int_Str& v_seq = *scenario.constructed_sequences[V_gene_seq];

// Check gene overlap safety
if (scenario.safety_set.value(D_gene, J_gene)) {
    // D and J genes are compatible
}
```

**Key Design**: ScenarioContext provides **helper methods** (added refactoring) for querying scenario state without directly accessing memory layer maps.

---

### 4. ExplorationContext

**Purpose**: Controls scenario space exploration (pruning, recursion, safety).

**Contains**:
- `next_event_ptr_arr`: Recursion control (which event to call next)
- `seq_max_prob_scenario`: Best scenario probability found so far (for pruning)
- `proba_threshold_factor`: Relative pruning threshold
- `index_map`: Parent realization tracking (shared with ScenarioContext)
- `safety_set`: Gene overlap safety (shared with ScenarioContext)

**Usage**:
```cpp
// Check if scenario should be pruned
double upper_bound = scenario.scenario_proba * max_possible_downstream_prob;
if (upper_bound < exploration.seq_max_prob_scenario * exploration.proba_threshold_factor) {
    continue;  // Prune this branch
}

// Update best scenario
if (scenario.scenario_error_w_proba > exploration.seq_max_prob_scenario) {
    exploration.seq_max_prob_scenario = scenario.scenario_error_w_proba;
}

// Recurse to next event
if (exploration.next_event_ptr_arr.get()[this->event_index]) {
    exploration.next_event_ptr_arr.get()[this->event_index]->iterate(
        query, model, scenario, exploration, accumulation
    );
}
```

**Key Design**: Separates exploration control (pruning, recursion) from scenario data (ScenarioContext).

---

### 5. AccumulationContext

**Purpose**: Accumulates results across scenarios and sequences.

**Contains**:
- `updated_marginals`: EM algorithm marginal updates (accumulated)
- `counters`: Statistics collectors (best scenarios, coverage, pgen, errors)
- `error_rate`: Error model (mutable - accumulates per-sequence statistics)
- `seq_likelihood`: Per-sequence likelihood accumulator
- `seq_probability`: Per-sequence probability accumulator

**Usage**:
```cpp
// Update marginals for EM algorithm
accumulation.updated_marginals[param_index] += scenario.scenario_error_w_proba;

// Call counters to record scenario
for (auto& [id, counter] : accumulation.counters) {
    counter->count_scenario(Scenario(scenario), query, model);
}

// Compute error-weighted probability
double error_w_proba = accumulation.error_rate->compute_scenario_error_probability(
    query, model, scenario, exploration
);
```

**Key Design**: **Error_rate lives here (not ModelContext)** because it maintains mutable per-sequence state. See architectural discussion in [ITERATE_REFACTORING_PLAN.md](ITERATE_REFACTORING_PLAN.md#50-architectural-analysis---error_rate-position).

---

## Context Field Mapping (Legacy → Context)

| Legacy Parameter | Context | Field |
|-----------------|---------|-------|
| `sequence` | `query` | `.sequence` |
| `int_sequence` | `query` | `.int_sequence` |
| `allowed_realizations` | `query` | `.allowed_realizations` |
| `model_parameters_point` | `model` | `.model_parameters` |
| `offset_map` | `model` | `.offset_map` |
| `events_map` | `model` | `.events_map` |
| `scenario_proba` | `scenario` | `.scenario_proba` |
| `constructed_sequences` | `scenario` | `.constructed_sequences` |
| `seq_offsets` | `scenario` | `.seq_offsets` |
| `mismatches_lists` | `scenario` | `.mismatches_lists` |
| `downstream_proba_map` | `scenario` | `.downstream_proba_map` |
| `index_map` | `exploration` | `.index_map` |
| `next_event_ptr_arr` | `exploration` | `.next_event_ptr_arr` |
| `safety_set` | `exploration` | `.safety_set` |
| `seq_max_prob_scenario` | `exploration` | `.seq_max_prob_scenario` |
| `proba_threshold_factor` | `exploration` | `.proba_threshold_factor` |
| `updated_marginals` | `accumulation` | `.updated_marginals` |
| `counters_list` | `accumulation` | `.counters` |
| `error_rate_p` | `accumulation` | `.error_rate` |

---

## Usage Examples

### Example 1: Implementing a New Event

```cpp
class MyCustomEvent : public Rec_Event {
public:
    void iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation
    ) override {
        // Get my probability parameters
        size_t param_start = this->marginal_array_index;
        
        // Loop over realizations
        for (size_t i = 0; i < num_realizations; ++i) {
            // Compute probability from model
            double event_prob = model.model_parameters[param_start + i];
            
            // Save current state
            double saved_proba = scenario.scenario_proba;
            
            // Update scenario
            scenario.scenario_proba *= event_prob;
            scenario.index_map.set_value(this->event_index, i, 0);
            
            // Update constructed sequences (example)
            scenario.constructed_sequences[V_gene_seq]->append("ACGT");
            scenario.seq_offsets.set_value(V_gene_seq, Five_prime, 0, 0);
            
            // Check pruning
            if (scenario.scenario_proba < exploration.seq_max_prob_scenario * 
                exploration.proba_threshold_factor) {
                scenario.scenario_proba = saved_proba;  // Restore
                continue;  // Prune this branch
            }
            
            // Recurse to next event
            this->iterate_wrap_up(query, model, scenario, exploration, accumulation);
            
            // Restore state for next iteration
            scenario.scenario_proba = saved_proba;
        }
    }
};
```

### Example 2: Testing an Event in Isolation

```cpp
TEST_CASE("Test Deletion event with contexts") {
    // Create minimal contexts for testing
    QuerySequenceContext query("ACGTACGT", Int_Str("01230123"), {});
    
    // Model context with dummy parameters
    Marginal_array_p params = std::make_unique<long double[]>(100);
    params[0] = 0.5;  // Deletion probability
    std::unordered_map<Rec_Event_name, 
        std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map;
    std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, 
        std::shared_ptr<Rec_Event>> events_map;
    std::vector<std::shared_ptr<Rec_Event>> model_queue;
    
    ModelContext model(params, offset_map, events_map, model_queue);
    
    // Scenario context with initial state
    double proba = 1.0;
    Downstream_scenario_proba_bound_map proba_map(6);
    Index_map idx_map(10);
    Seq_type_str_p_map sequences;
    Seq_offsets_map offsets;
    Safety_bool_map safety;
    Mismatch_vectors_map mismatches;
    
    ScenarioContext scenario(proba, proba_map, idx_map, sequences, 
                            offsets, safety, mismatches);
    
    // Exploration context
    std::shared_ptr<Next_event_ptr> next_ptr;
    ExplorationContext exploration(proba_map, 1.0, 0.01, idx_map, 
                                   next_ptr, safety);
    
    // Accumulation context
    Marginal_array_p marginals = std::make_unique<long double[]>(100);
    std::map<size_t, std::shared_ptr<Counter>> counters;
    std::shared_ptr<Error_rate> error_rate;
    AccumulationContext accumulation(marginals, counters, error_rate);
    
    // Test the event
    auto deletion = std::make_shared<Deletion>(/* ... */);
    deletion->iterate(query, model, scenario, exploration, accumulation);
    
    // Verify results
    REQUIRE(scenario.scenario_proba < 1.0);
    REQUIRE(scenario.index_map.get_value(deletion->event_index, 0) >= 0);
}
```

### Example 3: Implementing a Counter

Counters use a flattened **Scenario view** instead of raw ScenarioContext to simplify access:

```cpp
class MyCustomCounter : public Counter {
public:
    void initialize(const ModelContext& model) override {
        // Store event references needed for counting
        v_gene_event = model.events_map.at(
            std::make_tuple(GeneChoice_t, V_gene, Undefined_side)
        );
    }
    
    void count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model
    ) override {
        // Get V gene realization
        int v_real = scenario.get_realization_index(GeneChoice_t, V_gene);
        
        // Get V gene offsets
        auto [v_5p, v_3p] = scenario.get_offsets(V_gene_seq);
        
        // Get mismatches
        const std::vector<int>& v_mismatches = scenario.get_mismatches(V_gene_seq);
        
        // Record statistics
        stats[v_real] += scenario.scenario_error_w_proba;
    }
    
private:
    std::shared_ptr<Rec_Event> v_gene_event;
    std::map<int, double> stats;
};
```

---

## Migration Guide

### For Event Developers

If you're implementing a new event or modifying an existing one:

1. **Implement**: `iterate(QuerySequenceContext&, const ModelContext&, ScenarioContext&, ExplorationContext&, AccumulationContext&)`

2. **Access Parameters**: Use context fields instead of direct parameters
   - `scenario_proba` → `scenario.scenario_proba`
   - `seq_offsets` → `scenario.seq_offsets`
   - `sequence` → `query.sequence`
   - `model_parameters` → `model.model_parameters`

3. **Recursion**: Call `this->iterate_wrap_up(query, model, scenario, exploration, accumulation)` to recurse to next event

4. **State Management**: Save/restore `scenario.scenario_proba` when looping over realizations

### For Counter Developers

If you're implementing a new counter:

1. **Implement**: 
   - `initialize(const ModelContext& model)` - Store event references
   - `count_scenario(const Scenario&, const QuerySequenceContext&, const ModelContext&)` - Record statistics

2. **Use Scenario Accessors**:
   - `scenario.get_realization_index(event_type, gene)` - Get realization
   - `scenario.get_offsets(seq_type)` - Get sequence positions
   - `scenario.get_mismatches(seq_type)` - Get error positions
   - `scenario.get_sequence_segment(seq_type)` - Get sequence fragment

3. **Don't Access Memory Layers**: The Scenario view provides flattened access without memory layer details

---

## Architecture Highlights

### Zero-Copy Design

Contexts hold **references** to existing data structures, not copies:

```cpp
struct ScenarioContext {
    double& scenario_proba;                         // Reference
    Seq_type_str_p_map& constructed_sequences;      // Reference
    Seq_offsets_map& seq_offsets;                   // Reference
    // ...
};
```

**Benefits**:
- No allocation overhead
- Changes to context fields modify original data
- Minimal memory footprint

### Const-Correctness

**ModelContext is const** to enforce read-only model access:

```cpp
void iterate(
    QuerySequenceContext& query,              // Input (mutable for caching)
    const ModelContext& model,                // Model (read-only)
    ScenarioContext& scenario,                // Scenario state (mutable)
    ExplorationContext& exploration,          // Exploration (mutable)
    AccumulationContext& accumulation         // Results (mutable)
);
```

### Memory Layer Abstraction

**ScenarioContext** exposes memory layers explicitly:
```cpp
// Direct memory layer access
Int_Str& v_seq = *scenario.constructed_sequences[V_gene_seq];
```

**Scenario view** (for Counters) hides memory layers:
```cpp
// Flattened access (memory layer 0 implicit)
const Int_Str* v_seq = scenario.get_sequence_segment(V_gene_seq);
```

### Error_rate Placement

Error_rate is in **AccumulationContext** (not ModelContext) because:
- It maintains **mutable per-sequence state** (`seq_likelihood`, `seq_weighted_er`)
- ModelContext must be **const** (read-only model structure)
- Error_rate is both a model (has parameters) and an accumulator (collects statistics)

See detailed discussion in [ITERATE_REFACTORING_PLAN.md](ITERATE_REFACTORING_PLAN.md#50-architectural-analysis---error_rate-position).

---

## Performance Notes

### Benchmark Results (April 29, 2026)

Context-based architecture is **3% faster** than legacy implementation:

| Configuration | Legacy Baseline | Context-Based | Improvement |
|--------------|----------------|---------------|-------------|
| N=500, T=1, 2 iterations | 18.8254s | 18.2459s | **3.1% faster** |
| N=1000, T=1, 2 iterations | - | 35.8809s | - |
| N=500, T=4, 2 iterations | - | 5.2480s | - |

**Why faster?**
- Better memory locality (grouped parameters)
- Compiler optimization opportunities (const contexts)
- Reduced parameter passing overhead

---

## Legacy Compatibility

### Current State (April 2026)

- ✅ All `Rec_Event` subclasses use context-based interface
- ✅ All `Counter` subclasses use context-based interface
- ✅ All `Error_rate` subclasses use context-based interface
- ⚠️ Legacy interfaces preserved in base classes (deprecated, not yet removed)

### Hypermutation Models

Hypermutation error models currently use a **bridge pattern**:
- Context-based interface unpacks contexts
- Delegates to legacy implementation
- Preserves complex N-mer mutation logic

**Rationale**: Hypermutation models are complex (4^N parameters, sliding windows, boundary nucleotides). Full rewrite deferred until regression tests available.

---

## Further Reading

- **Refactoring Plan**: [ITERATE_REFACTORING_PLAN.md](ITERATE_REFACTORING_PLAN.md) - Complete refactoring history and decisions
- **Method Analysis**: [ITERATE_METHOD_ANALYSIS.md](ITERATE_METHOD_ANALYSIS.md) - Technical analysis of iterate() execution
- **Testing Plan**: [testing-refactoring-plan.md](testing-refactoring-plan.md) - Testing strategy and coverage

---

## Questions & Support

### How do I add a new parameter to iteration?

1. Identify which context it belongs to (query, model, scenario, exploration, accumulation)
2. Add field to appropriate context struct
3. Update context constructor
4. No need to change any iterate() signatures!

### Can I mix legacy and context code?

Yes, during migration. Legacy adapters exist in:
- `Rec_Event::iterate()` - NEW→OLD adapter (already removed in subclasses)
- `Counter::count_scenario()` - NEW→OLD adapter (already removed in subclasses)
- `Error_rate::compute_scenario_error_probability()` - NEW→OLD adapter (in hypermutation bridge implementations)

### Where should I put mutable state?

- **Per-sequence results**: `AccumulationContext` (marginals, counters, error_rate)
- **Per-scenario state**: `ScenarioContext` (proba, sequences, offsets)
- **Exploration control**: `ExplorationContext` (pruning, recursion)
- **Model structure**: Never! ModelContext is const

### How do I test contexts?

See Example 2 above. Create minimal context instances with only the fields your test needs. No need to fully initialize every field.

---

**Status**: ✅ Refactoring complete - All components use context-based architecture  
**Date**: April 30, 2026
