# Rec_Event::iterate Method Analysis

**Last Updated**: April 28, 2026 (Phase 1-5 Context Refactoring Complete)

## Overview

The `iterate` method is the core of IGoR's scenario exploration algorithm. It recursively explores all possible realizations of recombination events for a given sequence, computing scenario probabilities and performing tree pruning to eliminate improbable scenarios.

**Major Update (April 2026)**: The iterate() interface has been refactored from 19 individual parameters to 5 semantic context objects. This document reflects the current context-based architecture while preserving information about the legacy interface.

## Context Refactoring (Phase 1-5, Complete April 2026)

**Major Update**: The iterate() interface has been refactored to use semantic context objects instead of 19 individual parameters. This improves maintainability, type safety, and code clarity.

### Context Architecture

**5 Context Structs** organize parameters by semantic role:

1. **QuerySequenceContext** (input query) - read-only
   - `sequence`, `int_sequence`, `gene_alignments`

2. **ModelContext** (model structure) - const, read-only
   - `model_parameters`, `offset_map`, `events_map`, `model_queue`

3. **ScenarioContext** (scenario state) - mutable
   - `scenario_proba`, `constructed_sequences`, `seq_offsets`, `mismatches_lists`

4. **ExplorationContext** (exploration policy) - mutable
   - `downstream_proba_map`, `seq_max_prob_scenario`, `proba_threshold_factor`, `index_map`, `next_event_ptr_arr`, `safety_set`

5. **AccumulationContext** (result accumulation) - mutable
   - `updated_marginals`, `counters`, `error_rate`

**Benefits**:
- ✅ Clear semantic grouping (input vs model vs state vs policy vs output)
- ✅ Const-correctness enforced (ModelContext is const)
- ✅ Easier testing (can mock individual contexts)
- ✅ Reduced parameter coupling
- ✅ Future-proof (easier to extend contexts than add parameters)

### Migration Status

- ✅ **Phase 1**: Context definitions and legacy adapter
- ✅ **Phase 2**: All Rec_Event subclasses migrated
- ✅ **Phase 3**: Context helper methods (Scenario view, accessors)
- ✅ **Phase 4**: Counter interface refactored (uses Scenario view)
- ✅ **Phase 5**: Error_rate interface refactored (uses contexts)

**Current State**: Both interfaces coexist:
- **New (context-based)**: `iterate(query, model, scenario, exploration, accumulation)`
- **Old (legacy)**: `iterate(19 parameters...)` - delegates to new interface

---

## Method Signatures

### New Context-Based Signature (Phase 2, Complete)

```cpp
virtual void iterate(
    QuerySequenceContext& query,         // Input sequence and alignments
    const ModelContext& model,           // Model structure (const)
    ScenarioContext& scenario,           // Scenario state (mutable)
    ExplorationContext& exploration,     // Exploration policy (mutable)
    AccumulationContext& accumulation    // Result accumulation (mutable)
) = 0;
```

### Legacy Signature (Deprecated)

```cpp
virtual void iterate(
    double& scenario_proba,                          // [in,out] Current scenario probability
    Downstream_scenario_proba_bound_map& downstream_proba_map,  // [in,out] Probability bounds for downstream events
    const std::string& sequence,                     // [in] Nucleotide sequence
    const Int_Str& int_sequence,                     // [in] Integer-encoded sequence
    Index_map& base_index_map,                       // [in,out] Maps event indices to marginal array positions
    const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& offset_map,  // [in] Parent-child offset relationships
    std::shared_ptr<Next_event_ptr>& next_event_ptr_arr,  // [in] Pointer to next event in chain
    Marginal_array_p& updated_marginals_point,       // [in,out] Accumulates scenario posteriors
    const Marginal_array_p& model_parameters_point,  // [in] Current probability distribution
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>& allowed_realizations,  // [in] Pre-computed alignments
    Seq_type_str_p_map& constructed_sequences,       // [in,out] Partially constructed sequences
    Seq_offsets_map& seq_offsets,                    // [in,out] 5' and 3' offsets for each sequence segment
    std::shared_ptr<Error_rate>& error_rate_p,       // [in] Error model
    std::map<size_t,std::shared_ptr<Counter>>& counters_list,  // [in] Counters for statistics
    const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& events_map,  // [in] All model events
    Safety_bool_map& safety_set,                     // [in,out] Tracks overlap safety checks
    Mismatch_vectors_map& mismatches_lists,          // [in,out] Mismatch positions
    double& seq_max_prob_scenario,                   // [in] Best scenario probability (for pruning)
    double& proba_threshold_factor                   // [in] Pruning threshold multiplier
) = 0;
```

---

## Context Field Mapping Reference

Quick reference for migrating from legacy to context-based code:

| Legacy Parameter | Context | Field |
|-----------------|---------|-------|
| `sequence` | `query` | `.sequence` |
| `int_sequence` | `query` | `.int_sequence` |
| `allowed_realizations` | `query` | `.gene_alignments` |
| `model_parameters_point` | `model` | `.model_parameters` |
| `offset_map` | `model` | `.offset_map` |
| `events_map` | `model` | `.events_map` |
| `scenario_proba` | `scenario` | `.scenario_proba` |
| `constructed_sequences` | `scenario` | `.constructed_sequences` |
| `seq_offsets` | `scenario` | `.seq_offsets` |
| `mismatches_lists` | `scenario` | `.mismatches_lists` |
| `downstream_proba_map` | `exploration` | `.downstream_proba_map` |
| `seq_max_prob_scenario` | `exploration` | `.seq_max_prob_scenario` |
| `proba_threshold_factor` | `exploration` | `.proba_threshold_factor` |
| `base_index_map` | `exploration` | `.index_map` |
| `next_event_ptr_arr` | `exploration` | `.next_event_ptr_arr` |
| `safety_set` | `exploration` | `.safety_set` |
| `updated_marginals_point` | `accumulation` | `.updated_marginals` |
| `counters_list` | `accumulation` | `.counters` |
| `error_rate_p` | `accumulation` | `.error_rate` |

---

## Common Logic Across All Implementations

### 1. **Index Retrieval**
All implementations start by retrieving the base index for reading from marginal arrays:
```cpp
base_index = base_index_map.at(this->event_index);
```

### 2. **Realization Loop**
Each event iterates over its possible realizations (alignments for Gene_choice, deletion counts for Deletion, etc.)

### 3. **Probability Computation**
- Calls `iterate_common()` helper to:
  - Compute realization index in marginal array
  - Fetch probability from `model.model_parameters` (context-based) or `model_parameters_point` (legacy)
  - Update `exploration.index_map` (context) or `base_index_map` (legacy) for child events
- Multiplies `scenario.scenario_proba` (context) or `scenario_proba` (legacy) by `proba_contribution`

### 4. **Downstream Probability Bound**
Computes upper bound on downstream scenario probability:
```cpp
scenario_upper_bound_proba = new_scenario_proba;
downstream_proba_map.multiply_all(scenario_upper_bound_proba, current_downstream_proba_memory_layers);
```

### 5. **Tree Pruning**
Checks if scenario is probable enough to continue:

**Context-based**:
```cpp
if(!exploration.should_prune(scenario_upper_bound_proba)) {
    iterate_wrap_up(query, model, scenario, exploration, accumulation);
}
```

**Legacy**:
```cpp
if(scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
    iterate_wrap_up(...);  // Recursively call next event
}
```

### 6. **Recursive Call**
`iterate_wrap_up()` calls the next event's `iterate()` method, or if this is the last event (leaf node):

**Context-based (Phase 4-5)**:
- Calls `accumulation.error_rate->compute_scenario_error_probability(query, model, scenario, exploration)`
- Creates `Scenario` view: `Scenario scenario_view(scenario)`
- Calls counters: `counter->count_scenario(scenario_view, query, model)`
- Accumulates marginals: `event->add_to_marginals(scenario.scenario_error_w_proba, accumulation.updated_marginals)`

**Legacy**:
- Records complete scenario in `updated_marginals_point`
- Calls `error_rate_p->compare_sequences_error_prob()` with 8 parameters
- Calls `counter->count_scenario()` with 7 parameters

---

## Parameter Usage Details

### **scenario_proba** (double&)
- **Input**: Probability of the incomplete scenario up to this event
- **Modified**: Multiplied by this event's realization probability
- **Output**: Probability including this event's realization
- **Usage**: Core probability accumulator through recursive chain

### **downstream_proba_map** (Downstream_scenario_proba_bound_map&)
- **Purpose**: Maintains upper bounds on probabilities for downstream events
- **Key Operations**:
  - `set_value(seq_type, proba, memory_layer)`: Sets bound for a sequence segment
  - `multiply_all(result, memory_layers)`: Multiplies all bounds to compute total downstream bound
- **Example**: After Gene_choice, sets error rate upper bound for that gene

### **sequence** (const std::string&) & **int_sequence** (const Int_Str&)
- **sequence**: Nucleotide sequence ("ACGTACGT...")
- **int_sequence**: Integer-encoded sequence (A=0, C=1, G=2, T=3)
- **Usage**: 
  - Gene_choice: Aligns against templates
  - Deletion: Checks overlap constraints
  - Dinucl_markov: Extracts inserted sequence for probability calculation

### **base_index_map** (Index_map&)
- **Purpose**: Maps each event to its current index in the marginal probability arrays
- **Key Operations**:
  - `at(event_index)`: Get current base index
  - Modified to include current realization index for child events
- **Example**: If V_choice has 50 genes and gene #7 is chosen, child events read from `base_index + 7`

### **offset_map** (const unordered_map&)
- **Purpose**: Defines parent-child relationships and how parent realizations offset child indices
- **Structure**: `event_name -> [(parent_event_ptr, offset_value)]`
- **Usage**: `iterate_common()` uses this to update `base_index_map` correctly
- **Helper Function**: See `EventUtils::initialize_offset_memory()` for memory layer initialization pattern

### **next_event_ptr_arr** (shared_ptr<Next_event_ptr>&)
- **Purpose**: Points to the next event to call in the recursion chain
- **Usage**: `iterate_wrap_up()` uses this to determine which event's `iterate()` to call next

### **updated_marginals_point** (Marginal_array_p&)
- **Purpose**: Accumulates posterior probabilities for complete scenarios
- **Usage**: When recursion reaches the end, scenario probability is added to marginals
- **Type**: Array pointer to marginal distributions being updated during inference

### **model_parameters_point** (const Marginal_array_p&)
- **Purpose**: Current recombination probability distributions (prior or current estimate)
- **Usage**: `iterate_common()` reads from here using computed indices to get P(realization)

### **allowed_realizations** (const unordered_map&)
- **Purpose**: Pre-computed sequence alignments for Gene_choice events
- **Structure**: `gene_class -> [Alignment_data{gene_name, offset, mismatches}]`
- **Usage**: Only used by Gene_choice to iterate over plausible gene alignments

### **constructed_sequences** (Seq_type_str_p_map&)
- **Purpose**: Stores pointers to sequence segments constructed so far
- **Key Types**: V_gene_seq, D_gene_seq, J_gene_seq, VD_ins_seq, DJ_ins_seq, VJ_ins_seq
- **Usage**:
  - Gene_choice: Sets gene sequence pointer
  - Deletion: Creates truncated sequence
  - Insertion: Creates placeholder for inserted nucleotides
  - Dinucl_markov: Fills in inserted nucleotides

### **seq_offsets** (Seq_offsets_map&)
- **Purpose**: Tracks 5' and 3' alignment offsets for each sequence segment
- **Operations**: `at(seq_type, side, memory_layer)`, `set_value(...)`
- **Usage**: 
  - Gene_choice: Sets gene alignment offsets
  - Deletion: Modifies 3' offset after deletion
  - Insertion: Computes insertion length from offset gaps

### **error_rate_p** (shared_ptr<Error_rate>&) → **accumulation.error_rate**
- **Purpose**: Error model for computing mismatch probabilities
- **Context Location**: AccumulationContext (mutable - maintains per-sequence accumulators)
- **New Interface (Phase 5)**:
  ```cpp
  double compute_scenario_error_probability(
      const QuerySequenceContext& query,
      const ModelContext& model,
      ScenarioContext& scenario,
      ExplorationContext& exploration
  );
  ```
- **Implementation Status**:
  - ✅ Single_error_rate: Direct context-based implementation
  - ✅ Hypermutation_global_errorrate: Bridge to legacy (preserves N-mer logic)
  - ✅ Hypermutation_full_Nmer_errorrate: Bridge to legacy (preserves N-mer logic)
- **Usage**: 
  - Called at leaf nodes by `iterate_wrap_up()`
  - Computes error-weighted probability: `P(scenario) * P(errors | model)`
  - Gene_choice: Computes upper bound on error probability given mismatches
  - Used in downstream probability bounds
- **Design Note**: Error_rate stays in AccumulationContext (not ModelContext) because it maintains mutable per-sequence accumulators (seq_likelihood, seq_weighted_er)

### **counters_list** (map<size_t, shared_ptr<Counter>>&) → **accumulation.counters**
- **Purpose**: Collects statistics during iteration (e.g., edge activation counters, Pgen computation)
- **Context Location**: AccumulationContext
- **New Interface (Phase 4)**:
  ```cpp
  void count_scenario(
      const Scenario& scenario,           // Scenario view (flattened)
      const QuerySequenceContext& query,
      const ModelContext& model
  );
  ```
- **Scenario View**: Lightweight read-only accessor over ScenarioContext
  - Provides `get_realization_index()`, `get_offsets()`, `get_sequence_segment()`, etc.
  - Hides memory layer complexity from Counters
- **Usage**: 
  - Called at leaf nodes by `iterate_wrap_up()` for each counter
  - Passed through recursion, updated when complete scenarios are found
  - Examples: Best_scenarios_counter (Pgen), Coverage_errorrate_counter

### **events_map** (const unordered_map&)
- **Purpose**: Global map of all events in the model
- **Key**: `tuple<Event_type, Gene_class, Seq_side>`
- **Usage**: Events can query other events (e.g., check if D gene exists)
- **Helper Functions**: See `EventUtils::check_gene_choice()` and `EventUtils::get_insertion_len_max()` for common query patterns

### **safety_set** (Safety_bool_map&)
- **Purpose**: Tracks whether overlap safety checks have been resolved
- **Values**: VD_safe, VJ_safe, DJ_safe
- **Usage**: 
  - Gene_choice: Sets safety flags based on offset constraints
  - Deletion: Checks/updates safety based on deletion amount

### **mismatches_lists** (Mismatch_vectors_map&)
- **Purpose**: Stores positions of mismatches between constructed sequence and observed sequence
- **Usage**:
  - Gene_choice: Initializes from alignment mismatches
  - Deletion: Updates by removing mismatches in deleted regions
  - Used for error rate calculations

### **seq_max_prob_scenario** (double&)
- **Purpose**: Probability of the most likely complete scenario found so far
- **Usage**: Threshold for tree pruning - scenarios with `upper_bound < max * threshold` are discarded

### **proba_threshold_factor** (double&)
- **Purpose**: Multiplier for pruning threshold (typically 10^-5 to 10^-10)
- **Usage**: Controls trade-off between exhaustive search vs computational speed

---

## EventUtils Helper Functions

**Namespace**: `EventUtils` (defined in `src/igor/Core/EventUtils.h/cpp`)

Several common patterns have been extracted into reusable helper functions. These utilities reduce code duplication and make event implementations more readable.

### **check_gene_choice()**

```cpp
struct GeneChoiceStatus {
  bool exists;   // Does this gene class have a choice event?
  bool chosen;   // Has the choice been made (processed)?
  std::shared_ptr<const Rec_Event> event_ptr;  // Pointer to the event
};

GeneChoiceStatus check_gene_choice(
    Gene_class gene,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>>& events_map,
    const std::unordered_set<Rec_Event_name>& processed_events
);
```

**Purpose**: Check if a gene choice event (V, D, or J) exists in the model and whether it has been processed.

**Returns**: `GeneChoiceStatus` struct with three fields:
- `exists`: true if gene choice event is in events_map
- `chosen`: true if event has been processed (in processed_events set)
- `event_ptr`: pointer to the event (nullptr if doesn't exist)

**Usage Example** (from Gene_choice safety checks):
```cpp
auto d_status = EventUtils::check_gene_choice(D_gene, events_map, processed_events);
if (d_status.exists && !d_status.chosen) {
    // D gene event exists but not yet chosen - need to check safety
}
```

**Refactoring Benefit**: Replaces scattered `events_map.count()` and `processed_events.count()` checks with a single semantic function.

---

### **build_scenario_sequence()**

```cpp
Int_Str build_scenario_sequence(
    Seq_type_str_p_map& constructed_sequences,
    bool has_v, bool has_d, bool has_j,
    bool has_vd_ins, bool has_dj_ins, bool has_vj_ins
);
```

**Purpose**: Concatenate constructed sequence segments in the correct order to form the complete recombined sequence.

**Parameters**:
- `constructed_sequences`: Map of constructed sequence segments (V_gene_seq, VD_ins_seq, D_gene_seq, etc.)
- Boolean flags indicating which segments exist

**Returns**: Complete `Int_Str` sequence in proper order:
- **VDJ model**: V → VD_ins → D → DJ_ins → J
- **VJ model**: V → VJ_ins → J

**Implementation** (from EventUtils.cpp):
```cpp
Int_Str scenario_resulting_sequence;
if (has_v) {
    scenario_resulting_sequence += (*constructed_sequences[V_gene_seq]);
}
if (has_d) {
    if (has_vd_ins) {
        scenario_resulting_sequence += (*constructed_sequences[VD_ins_seq]);
    }
    scenario_resulting_sequence += (*constructed_sequences[D_gene_seq]);
    if (has_dj_ins) {
        scenario_resulting_sequence += (*constructed_sequences[DJ_ins_seq]);
    }
} else {
    if (has_vj_ins) {
        scenario_resulting_sequence += (*constructed_sequences[VJ_ins_seq]);
    }
}
if (has_j) {
    scenario_resulting_sequence += (*constructed_sequences[J_gene_seq]);
}
return scenario_resulting_sequence;
```

**Usage**: Called at leaf nodes (complete scenarios) to build the final sequence for comparison with observed sequence.

**Refactoring Benefit**: Centralizes sequence concatenation logic that was previously duplicated in iterate_wrap_up implementations.

---

### **get_insertion_len_max()**

```cpp
int get_insertion_len_max(
    Gene_class gene_pair,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>>& events_map
);
```

**Purpose**: Query the maximum insertion length for a junction (VD, DJ, or VJ).

**Parameters**:
- `gene_pair`: Gene class (VD_genes, DJ_genes, or VJ_genes)
- `events_map`: Global event map

**Returns**: Maximum insertion count for that junction, or 0 if insertion event doesn't exist.

**Implementation**:
```cpp
auto key = std::make_tuple(Insertion_t, gene_pair, Undefined_side);
if (events_map.count(key) != 0) {
    return events_map.at(key)->get_len_max();
}
return 0;
```

**Usage Example** (from Gene_choice computing junction length bounds):
```cpp
int vd_ins_len_max = EventUtils::get_insertion_len_max(VD_genes, events_map);
// Use to compute feasible offset ranges
```

**Refactoring Benefit**: Replaces direct events_map queries with semantic accessor.

---

### **initialize_offset_memory()**

```cpp
void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>& offset_vector,
    Index_map& index_map,
    std::forward_list<std::tuple<int, int, int>>& memory_and_offsets
);
```

**Purpose**: Initialize memory layers in Index_map for tracking parent event realizations, and populate memory_and_offsets list.

**Parameters**:
- `offset_vector`: Vector from offset_map (parent events and their contributions)
- `index_map`: Index map to request memory layers from
- `memory_and_offsets`: Output list of (event_id, memory_layer, offset_value) tuples

**Implementation**:
```cpp
for (auto iter = offset_vector.begin(); iter != offset_vector.end(); ++iter) {
    int event_identifier = (*iter).first->get_event_identifier();
    index_map.request_memory_layer(event_identifier);
    memory_and_offsets.emplace_front(
        event_identifier,
        index_map.get_current_memory_layer(event_identifier),
        (*iter).second
    );
}
```

**Usage**: Called during iterate_common() to set up memory layers for marginal index computation.

**Refactoring Benefit**: Encapsulates memory layer initialization pattern used by all events.

---

### Future EventUtils Extensions (Post-Context Refactoring)

Once contexts are introduced (see Phase 3), additional EventUtils functions could be added:
- `compute_marginal_index()` - Consolidate marginal navigation logic
- `check_safety_constraints()` - Centralized safety checking
- `update_downstream_bounds()` - Unified downstream probability updates

These would complement the context abstraction layer (Phase 3.5) by providing standalone helper functions that operate across multiple contexts.

---

## Event-Specific Implementations

### 1. **Gene_choice** (Genechoice.cpp)

**Purpose**: Models gene choice (V, D, or J) based on pre-computed sequence alignments. For D genes, if there's no alignment available, all putative sequence offsets are tried. 

**Specific Logic**:

1. **Switch on Gene Class** (V_gene, D_gene, J_gene):
   - Different constraints for each gene class

2. **Safety Checks**:
   - V_gene: Checks if V-D and V-J spacing allows valid deletions
   - D_gene: Checks if D fits between V and J
   - J_gene: Similar checks for J positioning

3. **Alignment Iteration**:
   ```cpp
   for(auto& iter : allowed_realizations.at(V_gene)) {
       gene_seq = event_realizations.at(iter.gene_name).value_str_int;
       // Handle negative offsets (partial gene visibility)
       if(iter.offset >= 0) {
           v_5_off = iter.offset;
       } else {
           gene_seq = gene_seq.substr(-iter.offset);
           v_5_off = 0;
       }
   ```

4. **Overlap Detection**:
   - Uses min/max deletion bounds to check if gene placements would overlap
   - Example: V with max deletions shouldn't reach D with min deletions
   ```cpp
   if((v_3_off + v_3_max_del) >= (d_5_max_offset)) {
       continue;  // Bad alignment, skip
   }
   ```

5. **Mismatch Handling**:
   - Copies mismatches from alignment data
   - Computes **endogenous mismatches**: These are mismatches between the gene template and the observed sequence that occur in the portion of the gene that will remain visible even after the maximum possible deletion. Since these mismatches cannot be "explained away" by deletions, they must be accounted for by the error model, which provides a lower bound on the error probability for this gene choice.
     ```cpp
     // Example: V gene aligned at position 0-100, with max 3' deletion = 15
     // Mismatches at positions 0-85 are endogenous (will remain after max deletion)
     // Mismatches at positions 86-100 might be removed by deletions
     endogeneous_mismatches = 0;
     mism_iter = iter->mismatches.begin();
     while((mism_iter != iter->mismatches.end()) && ((*mism_iter) <= v_3_off + v_3_max_del)) {
         ++endogeneous_mismatches;
         ++mism_iter;
     }
     ```
   - Updates downstream probability bound based on error rate using the count of endogenous mismatches and the length of the conserved region

6. **Junction Length Probability**:
   - Pre-computes best probability for each possible junction length
   - Sets downstream bound: `downstream_proba_map.set_value(VD_ins_seq, vd_length_best_proba_map.at(junction_length))`

**Key Private Members**:
- `vd_check`, `vj_check`, `dj_check`: Boolean flags for safety checks
- Offset variables: `d_5_min_offset`, `j_5_max_offset`, etc.
- `iterate_common()`: Reads gene choice probability from marginals

**D Gene Exhaustive Search**:

D genes require special handling because they can be heavily deleted on both ends during VDJ recombination, making them difficult or impossible to recognize via standard sequence alignment:

1. **Problem**: 
   - D genes are typically short (10-20bp)
   - Can have 0-10+ deletions on each end
   - After deletions, may only have 1-5bp visible in the read
   - Standard aligners (Smith-Waterman, BLAST) often fail to find such short matches

2. **Fallback Strategy**:
   - When `allowed_realizations.at(D_gene)` is empty (no alignments found), iterate over all possible positions
   - For each D gene template in the model:
     - Try every offset position between V_3' and J_5' offsets
     - Apply overlap safety constraints
     - Compute mismatches for each position/gene combination
   ```cpp
   if(no_d_align) {  // No alignments available
       // Exhaustive search: try all D genes at all plausible positions
       for(int d_pos = v_3_min_offset; d_pos <= j_5_max_offset; ++d_pos) {
           for(each D gene template) {
               // Try this gene at this position
               // Compute mismatches, check constraints
           }
       }
   }
   ```

3. **Computational Implications**:
   - **With alignments**: O(#alignments) typically 1-5 alignments per gene class
   - **Without alignments**: O(#D_genes × #possible_positions) typically 20-50 genes × 10-30 positions = 200-1500 scenarios
   - This makes D_choice the most expensive event in the recursion tree
   - Tree pruning is critical to keep this tractable
   - Pre-computed `vd_length_best_proba_map` and `dj_length_best_proba_map` help prune early

4. **Safety Constraints**:
   - D position must be far enough from V to allow minimum VD insertions
   - D position must be far enough from J to allow minimum DJ insertions
   - D gene length (minus max deletions) must fit in available space
   - These constraints dramatically reduce the search space

5. **Testing Consideration**:
   - Tests should cover both alignment-based and exhaustive search modes
   - Exhaustive search tests should verify all positions are tried (with single-realization models)
   - Performance tests should confirm pruning is effective

---

### 2. **Deletion** (Deletion.cpp)

**Purpose**: Models deletions (or palindromic insertions) at gene ends.

**Specific Logic**:

1. **Switch on Gene Class and Side**:
   - V_gene (3' deletions), D_gene (5' and 3'), J_gene (5' deletions)

2. **Feasibility Check**:
   ```cpp
   if((int)previous_str.size() > iter.value_int) {  // Don't delete entire gene
   ```

3. **Offset Update**:
   - Positive deletions: `v_3_new_offset = v_3_offset - deletion_value`
   - Checks if offset stays within read boundaries

4. **Safety Check Updates**:
   - If deletion amount moves offsets into "safe" or "unsafe" ranges, updates `safety_set`
   ```cpp
   if(v_3_new_offset < d_5_min_offset) {
       safety_set.set_value(Event_safety::VD_safe, true);
   } else if(v_3_new_offset >= d_5_max_offset) {
       continue;  // Overlap even with max D deletions
   }
   ```

5. **Sequence Modification**:
   - **Positive deletions**: Truncate sequence
     ```cpp
     new_str = previous_str.substr(0, previous_str.size() - deletion_value);
     ```
   - **Negative deletions (palindromes)**: Append reverse-complement
     ```cpp
     tmp_str = previous_str.substr(previous_str.size() + deletion_value);
     reverse(tmp_str.begin(), tmp_str.end());
     make_transversions(tmp_str);  // Complement
     new_str = previous_str + tmp_str;
     ```

6. **Mismatch List Update**:
   - Positive deletions: Remove mismatches beyond new offset
   - Negative deletions: Add mismatches for palindromic nucleotides

7. **Endogenous Mismatches**:
   - Counts mismatches that can't be removed by further deletions
   - Updates error rate bound in `downstream_proba_map`

**Key Private Members**:
- `int_value_and_index`: Forward list of deletion realizations (integer values)
- `iterate_common()`: Reads deletion probability, updates index map
- Safety check flags: `vd_check`, `vj_check`

**Loop Order**:
- Iterates deletions in **decreasing order** (max deletions first)
- This ordering helps with safety checks and pruning

---

### 3. **Insertion** (Insertion.cpp)

**Purpose**: Models the **number** of insertions at junctions (VD, DJ, VJ).

**Specific Logic**:

1. **Compute Insertion Count from Offsets**:
   ```cpp
   // VD junction
   insertions = seq_offsets.at(D_gene_seq, Five_prime) 
              - seq_offsets.at(V_gene_seq, Three_prime) - 1;
   
   // DJ junction  
   insertions = seq_offsets.at(J_gene_seq, Five_prime)
              - seq_offsets.at(D_gene_seq, Three_prime) - 1;
   
   // VJ junction (no D)
   insertions = seq_offsets.at(J_gene_seq, Five_prime)
              - seq_offsets.at(V_gene_seq, Three_prime) - 1;
   ```

2. **Probability Lookup**:
   - `iterate_common()` checks if insertion count is in allowed range
   - Returns 0 if out of range (discards scenario)
   ```cpp
   if(ordered_realization_map.count(insertions) > 0) {
       realization_index = ordered_realization_map.at(insertions).index;
       proba_contribution = model_parameters_point[base_index + realization_index];
   } else {
       proba_contribution = 0;  // Discard scenario
   }
   ```

3. **Placeholder Creation**:
   - Creates an integer string filled with -1 (unknown nucleotides)
   ```cpp
   inserted_str.assign(insertions, -1);
   constructed_sequences[VD_ins_seq] = &inserted_str;
   ```

4. **Downstream Dinucleotide Bound Update**:
   - Sets upper bound for dinucleotide model probability
   ```cpp
   (*dinuc_updated_bound) = upper_bound_per_ins.at(insertions);
   ```

5. **Junction Length Probability Map**:
   - Uses pre-computed `junction_length_best_proba_map`
   - Sets downstream bound: `downstream_proba_map.set_value(VD_ins_seq, junction_length_best_proba_map.at(insertions))`

**Key Private Members**:
- `ordered_realization_map`: Map from insertion count to realization
- `junction_length_best_proba_map`: Pre-computed best probabilities for each length
- `inserted_str`: Mutable Int_Str for inserted sequence placeholder
- `upper_bound_per_ins`: Upper bounds on dinucleotide probability per insertion count

**No Looping**:
- Insertion count is **deterministic** given the offsets
- Only one realization per scenario (unlike Gene_choice or Deletion)

---

### 4. **Dinucl_markov** (Dinuclmarkov.cpp)

**Purpose**: Models the **identity** of inserted nucleotides using a dinucleotide Markov chain.

**Specific Logic**:

1. **Switch on Gene Class** (VD_genes, DJ_genes, VJ_genes, VDJ_genes):
   - Can handle multiple junctions if D gene exists

2. **Extract Relevant Sequences**:
   - Gets previous gene sequence (to start Markov chain)
   - Gets insertion placeholder (created by Insertion event)
   - Extracts corresponding read subsequence
   ```cpp
   // VD junction
   previous_seq = (*constructed_sequences.at(V_gene_seq));
   vd_seq = (*constructed_sequences.at(VD_ins_seq));
   vd_seq_size = vd_seq.size();
   previous_nt_str = previous_seq.back();  // Last nucleotide of V
   data_seq_substr = int_sequence.substr(seq_offsets.at(V_gene_seq, Five_prime) 
                                       + previous_seq_size, vd_seq.size());
   ```

3. **Markov Chain Heuristic**:
   - Instead of full forward algorithm, uses heuristic:
     - Assumes low error rate
     - Picks most likely nucleotide at each position given previous nucleotide
   ```cpp
   iterate_common(vd_realizations_indices, previous_nt_str, vd_seq, model_parameters_point);
   ```
   - `iterate_common()` fills in `vd_seq` with nucleotide identities and computes probability

4. **DJ Junction Reverse Processing**:
   - DJ insertions are processed from J backwards (3' to 5')
   ```cpp
   // Reverse the data sequence
   reverse(data_seq_substr.begin(), data_seq_substr.end());
   iterate_common(dj_realizations_indices, previous_nt_str, dj_seq, model_parameters_point);
   // Reverse the result back
   reverse(dj_seq.begin(), dj_seq.end());
   ```

5. **Realization Indices**:
   - Stores indices for each inserted nucleotide
   - Used to read from dinucleotide transition probability matrix
   - `current_realizations_index_vec` holds all nucleotide indices

6. **Downstream Probability**:
   - Sets downstream bound to 1.0 (no further uncertainty)
   ```cpp
   downstream_proba_map.set_value(VD_ins_seq, 1.0, memory_layer_proba_map_junction);
   ```

**Key Private Members**:
- `dinuc_proba_matrix`: Matrix of dinucleotide transition probabilities
- `vd_realizations_indices`, `dj_realizations_indices`, `vj_realizations_indices`: Arrays storing nucleotide indices
- `iterate_common()`: Fills in nucleotide sequence using Markov heuristic

**No Traditional Looping**:
- Uses heuristic to pick single most likely nucleotide sequence
- Doesn't explore all possible insertion sequences (would be 4^n exponential)

---

## Execution Flow Example

Consider a VDJ recombination with the event order: **V_choice → V_del → D_choice → D_del_5 → D_del_3 → VD_ins → VD_dinuc → J_choice → J_del → DJ_ins → DJ_dinuc**

**Context-Based Flow (Phase 2+)**:

1. **V_choice::iterate(query, model, scenario, exploration, accumulation)** is called first:
   - Loops over V gene alignments from `query.gene_alignments`
   - For each alignment:
     - Computes `scenario.scenario_proba *= P(V_gene_i)` (from `model.model_parameters`)
     - Sets `scenario.constructed_sequences[V_gene_seq]`
     - Sets `scenario.seq_offsets[V_gene_seq]`
     - Calls `iterate_wrap_up(query, model, scenario, exploration, accumulation)` → recursively calls **V_del::iterate()**

2. **V_del::iterate(...)** (V 3' deletion):
   - Loops over deletion counts (e.g., 0 to 15)
   - For each deletion:
     - Computes `scenario.scenario_proba *= P(del=k | V_gene_i)`
     - Modifies `scenario.constructed_sequences[V_gene_seq]` (truncates)
     - Updates `scenario.seq_offsets[V_gene_seq][Three_prime]`
     - Checks safety (does V still fit before D?) using `exploration.safety_set`
     - Calls `iterate_wrap_up(...)` → recursively calls **D_choice::iterate()**

3. **D_choice::iterate(...)**:
   - Loops over D gene alignments from `query.gene_alignments`
   - Similar to V_choice but with D-specific constraints
   - Calls **D_del_5::iterate(...)**

4. **VD_ins::iterate(...)** (VD insertion count):
   - **No loop** - insertion count is deterministic from offsets
   - Computes `insertions = d_5_offset - v_3_offset - 1`
   - Looks up `P(n_ins = insertions)` from `model.model_parameters`
   - Creates placeholder sequence of length `insertions`
   - Calls **VD_dinuc::iterate(...)**

5. **VD_dinuc::iterate(...)** (VD insertion identities):
   - **No traditional loop** - uses Markov heuristic
   - Fills in nucleotide identities in placeholder
   - Computes `scenario.scenario_proba *= P(sequence | dinucleotide model)`
   - Calls **J_choice::iterate(...)**

6. **Final Event (Leaf Node)** - `iterate_wrap_up()` detects end of chain:
   - **Phase 5**: Calls `accumulation.error_rate->compute_scenario_error_probability(query, model, scenario, exploration)`
     - Computes error-weighted probability
     - Stores result in `scenario.scenario_error_w_proba`
   - **Phase 4**: Checks pruning threshold via `exploration.should_prune()`
   - Creates Scenario view: `Scenario scenario_view(scenario)` 
   - **Phase 4**: Calls counters: `accumulation.counters[i]->count_scenario(scenario_view, query, model)`
   - Accumulates marginals: `event->add_to_marginals(scenario.scenario_error_w_proba, accumulation.updated_marginals)`
   - Returns back up the recursion stack

7. **Tree Pruning**:
   - At each step: `if(exploration.should_prune(scenario_upper_bound_proba))`
     - Skip recursive call
     - Return immediately (prune this branch)
   - Helper method: `exploration.update_max_prob(scenario_error_w_proba)` updates best scenario probability

**Legacy Flow**: Same logic but with 19 individual parameters instead of 5 context objects. Legacy interface delegates to context-based implementation via adapter.

---

## Memory Layers

Many maps use "memory layers" to allow multiple concurrent explorations without copying:
- `constructed_sequences.at(V_gene_seq, memory_layer)`
- `seq_offsets.at(V_gene_seq, Five_prime, memory_layer)`
- `safety_set.at(Event_safety::VD_safe, memory_layer)`

Each recursive level uses a different memory layer to avoid conflicts.

---

## Testing Strategy for `iterate` Method

**Context Refactoring Impact**: Testing is now easier with semantic context objects. Can mock individual contexts rather than initializing 19 parameters.

Given this analysis, unit tests should focus on:

### 1. **Gene_choice Tests**:
- Single alignment, verify offset setting
- Multiple alignments, verify probability sum
- Negative offsets (partial gene visibility)
- Overlap detection (V-D, V-J constraints)
- Mismatch counting
- Downstream bound setting
- **Context-based**: Mock `query.gene_alignments` with test alignments

### 2. **Deletion Tests**:
- Positive deletions: sequence truncation, offset update
- Negative deletions (palindromes): sequence extension, complementation
- Safety check updates (VD_safe, VJ_safe transitions)
- Mismatch list updates
- Boundary cases (0 deletions, maximum deletions)
- **Context-based**: Verify `scenario.seq_offsets` and `exploration.safety_set` modifications

### 3. **Insertion Tests**:
- Offset-based insertion count computation
- Out-of-range insertion counts (should discard)
- Placeholder sequence creation
- **Context-based**: Verify `scenario.constructed_sequences` placeholder creation
- Downstream bound updates for dinucleotide model

### 4. **Dinucl_markov Tests**:
- VD, DJ, VJ junction handling
- Reverse processing for DJ junction
- Nucleotide sequence filling
- Probability computation via Markov chain

### 5. **Integration Tests**:
- Multi-event chains (V → V_del → D → ...)
- Tree pruning effectiveness
- Complete scenario probability accumulation
- Error rate interaction

---

## Testing Approaches: Detailed Analysis

### Approach 1: Direct iterate() Calls (Unit Testing)

#### Objects Mutated During iterate()

**Context Refactoring Note**: With Phase 1-5 complete, these parameters are now organized into semantic context objects, making testing easier to set up and reason about.

When `iterate()` is called, it explores multiple scenarios recursively. The following objects are mutated and track the exploration state:

**Per-Scenario State (ScenarioContext + ExplorationContext):**
- `scenario.scenario_proba` (double&): Accumulates probability as realizations are chosen
- `exploration.index_map` (Index_map&): Updated to include current realization indices for child events
- `scenario.constructed_sequences` (Seq_type_str_p_map&): Sequence segments progressively built (V_gene_seq, D_gene_seq, etc.)
- `scenario.seq_offsets` (Seq_offsets_map&): 5' and 3' offsets for each segment
- `exploration.safety_set` (Safety_bool_map&): Boolean flags for overlap safety (VD_safe, VJ_safe, DJ_safe)
- `scenario.mismatches_lists` (Mismatch_vectors_map&): Mismatch positions per sequence type
- `exploration.downstream_proba_map` (Downstream_scenario_proba_bound_map&): Upper bounds on downstream probabilities

**Accumulated Results (AccumulationContext + ExplorationContext):**
- `accumulation.updated_marginals` (Marginal_array_p&): Accumulates posterior probabilities for complete scenarios
- `exploration.seq_max_prob_scenario` (double&): Tracks the highest scenario probability found (updated by iterate_wrap_up at leaf nodes)

**Memory Layer Architecture:**
- Most maps use memory layers: `map.at(key, memory_layer)` or `map.at(key1, key2, memory_layer)`
- Each recursion level uses a different layer to avoid conflicts
- Allows concurrent exploration without copying entire maps
- Example: `constructed_sequences.at(V_gene_seq, layer=0)` for first event, `layer=1` for second event, etc.

#### Extracting Subscenario Information

**Challenge: iterate_wrap_up is not easily mockable**
- `iterate_wrap_up()` is a protected/private method called internally by each event
- It handles the recursive call to the next event: `exploration.next_event_ptr_arr->iterate(...)`
- Cannot easily intercept without modifying source code or using inheritance

**Strategy 1: Inspect State After iterate() Returns (Context-Based)**
```cpp
// After calling iterate(), inspect mutated context objects:
REQUIRE(scenario.constructed_sequences.at(V_gene_seq, 0) != nullptr);
REQUIRE(scenario.seq_offsets.at(V_gene_seq, Five_prime, 0) == expected_offset);

// For single event test (next_event_ptr = nullptr):
// iterate() explores all realizations but doesn't recurse
// Can examine which realizations were tried by checking accumulation.updated_marginals
```

**Strategy 2: Pass nullptr for next_event_ptr**
- Stops recursion at the first event
- iterate() loops over realizations but doesn't call iterate_wrap_up recursively
- **Problem**: Many implementations check `if(exploration.next_event_ptr_arr != nullptr)` and behave differently
- **Result**: Tests single event logic but not the full recursion

**Strategy 3: Create a Mock Event as "Next Event"**
```cpp
class MockNextEvent : public Rec_Event {
public:
    int call_count = 0;
    std::vector<double> observed_probabilities;
    std::vector<std::string> observed_sequences;
    
    // Context-based interface (Phase 2)
    void iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation
    ) override {
        call_count++;
        observed_probabilities.push_back(scenario.scenario_proba);
        // Capture scenario state
        // Don't recurse further
    }
};

// In test:
auto mock_next = std::make_shared<MockNextEvent>();
exploration.next_event_ptr_arr.get()[v_choice_ptr->get_event_identifier()] = mock_next;
v_choice_ptr->iterate(query, model, scenario, exploration, accumulation);

// Verify:
REQUIRE(mock_next->call_count == expected_subscenarios);
REQUIRE(mock_next->observed_probabilities[0] == Approx(0.5));
```
**Limitation**: Requires careful setup of events_map and next_event_ptr relationships

**Strategy 4: Use Custom Counter to Capture State (Recommended)**

Counters are designed for this purpose! The `Counter` class has been refactored (Phase 4):

**New Context-Based Interface**:
```cpp
virtual void count_scenario(
    const Scenario& scenario,           // Phase 3 - Scenario view  
    const QuerySequenceContext& query,
    const ModelContext& model
);
```

**Legacy Interface (Deprecated)**:
```cpp
virtual void count_scenario(
    long double scenario_proba,
    double scenario_seq_probability,
    const std::string& original_sequence,
    Seq_type_str_p_map& constructed_sequences,
    const Seq_offsets_map& seq_offsets,
    const std::unordered_map<...>& events_map,
    Mismatch_vectors_map& mismatches
);
```

This is called by `iterate_wrap_up()` when reaching a leaf node (complete scenario).

**Custom Counter Implementation:**
```cpp
class TestObserverCounter : public Counter {
public:
    struct ScenarioSnapshot {
        double proba;
        std::string v_sequence;
        int v_5_offset;
        int v_3_offset;
        std::vector<int> v_mismatches;
    };
    
    std::vector<ScenarioSnapshot> observed_scenarios;
    
    void count_scenario(
        long double scenario_proba,
        double scenario_seq_probability,
        const std::string& original_sequence,
        Seq_type_str_p_map& constructed_sequences,
        const Seq_offsets_map& seq_offsets,
        const std::unordered_map<...>& events_map,
        Mismatch_vectors_map& mismatches
    ) override {
**Custom Counter Implementation (Phase 4 - Context-Based)**:
```cpp
class TestObserverCounter : public Counter {
public:
    struct ScenarioSnapshot {
        double proba;
        std::string v_sequence;
        int v_5_offset;
        int v_3_offset;
        std::vector<int> v_mismatches;
    };
    
    std::vector<ScenarioSnapshot> observed_scenarios;
    
    // Phase 4: New context-based interface
    void count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model
    ) override {
        ScenarioSnapshot snapshot;
        snapshot.proba = scenario.scenario_error_w_proba;
        
        // Use Scenario view accessors (Phase 3)
        const Int_Str* v_seq = scenario.get_sequence_segment(V_gene_seq);
        if (v_seq) {
            snapshot.v_sequence = *v_seq;
        }
        
        auto [v_5_off, v_3_off] = scenario.get_offsets(V_gene_seq);
        snapshot.v_5_offset = v_5_off;
        snapshot.v_3_offset = v_3_off;
        
        const std::vector<int>* v_mismatches = scenario.get_mismatches(V_gene_seq);
        if (v_mismatches) {
            snapshot.v_mismatches = *v_mismatches;
        }
        
        observed_scenarios.push_back(snapshot);
    }
};

// In test (context-based):
auto observer = std::make_shared<TestObserverCounter>();
accumulation.counters[0] = observer;

v_choice_ptr->iterate(query, model, scenario, exploration, accumulation);

// Assertions:
REQUIRE(observer->observed_scenarios.size() == expected_count);
REQUIRE(observer->observed_scenarios[0].proba == Approx(0.5));
REQUIRE(observer->observed_scenarios[0].v_5_offset == 0);
```

**Advantages of Counter Approach:**
- ✅ No source code modification needed
- ✅ Captures complete scenario state at leaf nodes via Scenario view (Phase 3)
- ✅ Can test multi-event chains (V→V_del→D→...)
- ✅ Matches how IGoR actually collects statistics
- ✅ Can verify probability accumulation, sequence construction, offset computation
- ✅ Clean interface with Scenario accessor methods (hides memory layer complexity)

**Limitations of Counter Approach:**
- ❌ Only captures **complete** scenarios (at leaf nodes), not intermediate events
- ❌ Doesn't capture subscenarios that get pruned before completion
- ❌ Can't test single event in isolation (needs full chain to leaf)

#### Recommendation for Direct iterate() Testing

**Best Hybrid Approach:**

1. **For single-event logic** (unit tests):
   - Inspect context state immediately after iterate() returns
   - Check `scenario.constructed_sequences`, `scenario.seq_offsets`, `exploration.index_map`
   - Verify `accumulation.updated_marginals` has non-zero values at expected indices
   - Use short chains (V_choice only, or V_choice→V_del)

2. **For multi-event chains** (integration tests):
   - Use custom Counter subclass to capture complete scenarios (context-based interface)
   - Verify scenario count, probabilities, sequence construction
   - Test conditional dependencies (e.g., VD safety affects D placement)
   - Use Scenario view accessors for clean, testable code

3. **For subscenario exploration** (advanced):
   - Create mock events for "next" pointer (with context-based interface)
   - Or instrument source code with logging (for debugging only)

---

### Approach 2: GenModel-Based Testing (Integration Testing)

#### How GenModel Orchestrates iterate()

```cpp
// In GenModel::compute_Pgen_given_seq():
1. Read model from files (Model_Parms::initialize(...))
2. Build event queue: model_queue = model_parms.get_model_queue()
3. Initialize maps: index_map, offset_map, etc.
4. Compute alignments: aligner.compute_alignments(sequence)
5. Initialize events: first_event->initialize_event(...)
6. Call first event: first_event->iterate(...) // Legacy adapter wraps to context-based
7. Collect results from updated_marginals
```

The key advantage: **All setup is automatic and realistic**.

#### Partial Model Feasibility

**Question: Can we create models with only V_choice + V_del (no D, no J)?**

Let's analyze the constraints:

**Event Dependencies (from source code analysis):**

1. **Gene_choice constraints:**
   - V_choice checks for D and J events via `events_map` to compute safety bounds
   - Code: `if(events_map.count(std::make_tuple(GeneChoice_t, D_gene, Undefined_side)))`
   - **Refactored**: Now uses `EventUtils::check_gene_choice(D_gene, events_map, processed_events)`
   - **Benefit**: Returns struct with `exists`, `chosen`, `event_ptr` fields (clearer than count checks)
   - **Workaround**: Code uses `.count()` to check existence, so should handle missing events gracefully

2. **Deletion constraints:**
   - V_del references D_gene offsets for safety checks: `d_5_min_offset`, `d_5_max_offset`
   - Code: `seq_offsets.at(D_gene_seq, Five_prime, ...)`
   - **Problem**: If D doesn't exist, these offsets are undefined
   - **Workaround**: Code might use `.count()` or try/catch, but likely will crash

3. **Junction length probability maps:**
   - `iterate_initialize_Len_proba()` computes best probabilities for each junction length
   - Requires downstream events (VD_ins, D_choice, etc.)
   - **Problem**: If chain is incomplete, these maps aren't computed correctly

4. **Error rate calculations:**
   - Downstream probability bounds assume complete sequences
   - Missing events might leave `downstream_proba_map` partially uninitialized

**Potential Solutions for Partial Models:**

**Option A: Stub Events**
- Create minimal D_choice, J_choice events with single realization
- Set probabilities to 1.0
- Allows V_del to access offsets without crashing
- Tests focus on V logic but model is still "complete"

**Option B: Modify Safety Checks**
- Would require source code changes to make safety checks optional
- Check if D/J exist before accessing their offsets
- Not practical for testing without code modification

**Option C: VJ Model (skip D entirely)**
- Many TCR beta chains don't use D genes (VJ recombination)
- IGoR supports VJ models: `events_map` checks `if (D_gene exists)`
- Should work if model files are structured correctly
- Test: V_choice → V_del → VJ_ins → VJ_dinuc → J_choice → J_del

**Option D: End-to-End Test with Single Realizations**
- Keep complete VDJ model but with 1 realization per event
- Deterministic outcomes, easy to verify
- Still tests full chain including dependencies

**Limitations of Current Code for Partial Models:**

Based on source analysis, likely issues:
1. ❌ V_del alone (without J) will fail: needs J_5_offset for safety
2. ❌ V_choice alone might work but won't compute junction probabilities correctly
3. ✅ VJ model (no D) should work: code has conditional checks for D existence
4. ✅ Complete VDJ with single realizations works: deterministic and testable

**Testing via Model Files:**

You mentioned creating model files for VDJ/VJ. Key considerations:

```
models/
  ├── minimal_vj/          # No D gene
  │   ├── v_gene_CDR3_anchors.csv
  │   ├── j_gene_CDR3_anchors.csv
  │   ├── v_3_del.csv
  │   ├── j_5_del.csv
  │   ├── vj_ins.csv
  │   ├── vj_dinuc.csv
  │   └── marginals.txt
  │
  ├── minimal_vdj/         # Complete but 1-2 realizations each
  │   ├── v_gene_CDR3_anchors.csv  (1 gene)
  │   ├── d_gene_CDR3_anchors.csv  (1 gene)
  │   ├── j_gene_CDR3_anchors.csv  (1 gene)
  │   ├── v_3_del.csv              (2-3 deletion values)
  │   ├── d_5_del.csv              (2-3 values)
  │   └── ...
  │
  └── test_v_only/         # ❌ Likely won't work
      ├── v_gene_CDR3_anchors.csv
      └── v_3_del.csv
```

**Recommendation**: Use **minimal_vdj** with single/few realizations for deterministic GenModel tests.

---

### Comparison: Direct iterate() vs GenModel

| Aspect | Direct iterate() | GenModel |
|--------|------------------|----------|
| **Setup Complexity** | ⚠️ High - manual map initialization | ✅ Low - automatic from files |
| **Execution Speed** | ✅ Fast - no file I/O | ⚠️ Slower - reads files, computes alignments |
| **Isolation** | ✅ Can test single events | ❌ Tests full chain |
| **Edge Cases** | ✅ Easy to create artificial scenarios | ⚠️ Limited by realistic models |
| **Debugging** | ✅ Fine-grained control, easier to isolate failures | ⚠️ Harder to pinpoint issues in long chains |
| **Realism** | ⚠️ Might miss real-world interactions | ✅ Tests actual usage patterns |
| **Conditional Dependencies** | ⚠️ Hard to test (need mock chains) | ✅ Natural - dependencies built into model |
| **Subscenario Observation** | ⚠️ Requires custom Counters or state inspection | ✅ Counter-based observation works well |
| **Probability Verification** | ⚠️ Manual probability calculation needed | ✅ Can compare against known Pgen values |
| **Maintenance** | ⚠️ Breaks if iterate() signature changes | ✅ Stable - relies on public API |

---

### Recommended Testing Strategy

**Tier 1: Direct iterate() - Unit Tests (30% of tests)**
- **Purpose**: Verify individual event logic, with single realizations based tests for the event of interest.
- **Scope**: Single event
- **Tools**: `call_iterate()` wrapper, state inspection
- **Examples**:
  - Gene_choice: offset setting, mismatch recording, safety checks
  - Deletion: sequence truncation, offset updates, palindrome handling
  - Insertion: deterministic insertion count computation
  - Dinucl_markov: nucleotide filling, reverse processing

**Tier 2: Direct iterate() with Counters - Chain Tests (40% of tests)**
- **Purpose**: Verify individual event logic with several realizations, Verify simple multi-event interactions
- **Scope**: V→V_del→D→D_del chains (3-5 events)
- **Tools**: Custom Counter subclass to capture scenarios
- **Examples**:
  - V_choice → V_del: verify deletions modify V offsets correctly
  - V→V_del→VD_ins: verify insertion count computed from offsets
  - D exhaustive search with pruning
  - Conditional dependencies: VD safety affects D placement

**Tier 3: GenModel - Integration Tests (30% of tests)**
- **Purpose**: Verify complete workflows and realistic scenarios
- **Scope**: Full VDJ or VJ recombination chains
- **Tools**: Model files with 1-2 realizations per event
- **Examples**:
  - Compute Pgen for known sequences
  - Verify scenario exploration count (with pruning vs without)
  - Test alignment-based vs exhaustive D search
  - Validate marginal accumulation correctness

---

### Example Test Cases

**Direct iterate() - Single Event:**
```cpp
TEST_CASE("V_choice sets offsets correctly") {
    Gene_choice v_choice(V_gene);
    v_choice.add_realization("TRBV1", "ACGTACGTACGT");
    
    auto state = create_iterate_state("ACGTACGTACGTNNNNNN");
    state.model_marginals[0] = 1.0;
    
    // Create single alignment at offset 0
    Alignment_data align = create_mock_alignment_data("TRBV1", 0, 0, 11, {});
    std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed;
    allowed[V_gene] = {align};
    
    // Call iterate
    call_iterate(v_choice_ptr, state, allowed, events_map);
    
    // Verify offsets were set
    REQUIRE(state.seq_offsets.at(V_gene_seq, Five_prime, 0) == 0);
    REQUIRE(state.seq_offsets.at(V_gene_seq, Three_prime, 0) == 11);
    REQUIRE(state.constructed_sequences.at(V_gene_seq, 0) != nullptr);
}
```

**Direct iterate() with Counter - Multi-Event Chain:**
```cpp
TEST_CASE("V_choice followed by V_del modifies offsets") {
    auto observer = std::make_shared<TestObserverCounter>();
    state.counters_list[0] = observer;
    
    // Set up V_choice → V_del chain
    // ... (create events, set marginals, alignments)
    
    call_iterate(v_choice_ptr, state, allowed, events_map);
    
    // Verify multiple scenarios were explored
    REQUIRE(observer->observed_scenarios.size() == 3); // 3 deletion values
    
    // Verify deletion modified offsets
    REQUIRE(observer->observed_scenarios[0].v_3_offset == 11); // del=0
    REQUIRE(observer->observed_scenarios[1].v_3_offset == 10); // del=1
    REQUIRE(observer->observed_scenarios[2].v_3_offset == 9);  // del=2
    
    // Verify probabilities sum correctly
    double total_proba = std::accumulate(
        observer->observed_scenarios.begin(),
        observer->observed_scenarios.end(),
        0.0,
        [](double sum, const auto& s) { return sum + s.proba; }
    );
    REQUIRE(total_proba == Approx(1.0));
}
```

**GenModel - Integration Test:**
```cpp
TEST_CASE("Compute Pgen for deterministic VDJ scenario") {
    // Load minimal model: 1 V, 1 D, 1 J, 2 deletion values each
    GenModel model;
    model.load_model("models/minimal_vdj/");
    
    // Construct sequence: V(del=0) + VD_ins(3bp) + D(del=0,0) + DJ_ins(2bp) + J(del=0)
    std::string sequence = "ACGTACGTACGT" + "NNN" + "TGCA" + "NN" + "GGCCGGCC";
    
    double pgen = model.compute_Pgen_given_seq(sequence);
    
    // Expected: P(V) * P(V_del=0) * P(VD_ins=3) * P(VD_dinuc) * 
    //           P(D) * P(D_5_del=0) * P(D_3_del=0) * P(DJ_ins=2) * P(DJ_dinuc) *
    //           P(J) * P(J_del=0) * P(error_rate)
    double expected_pgen = 1.0 * 0.5 * 0.1 * (dinuc^3) * 1.0 * 0.5 * 0.5 * 0.1 * (dinuc^2) * 1.0 * 0.5 * 0.99;
    
    REQUIRE(pgen == Approx(expected_pgen).epsilon(0.01));
}
```

---

## Key Insights for Test Design

1. **Isolation is Challenging**: `iterate` requires significant setup:
   - `allowed_realizations` (for Gene_choice)
   - `constructed_sequences`, `seq_offsets` (from previous events)
   - `base_index_map` (properly initialized)
   - `downstream_proba_map`, `safety_set`, `mismatches_lists`

2. **Mock Setup Utility**: Need comprehensive test utilities to:
   - Create minimal valid state for testing single events
   - Mock alignment data
   - Initialize maps with correct structure

3. **Focus on Event Logic, Not I/O**:
   - Don't test model file parsing (Model_Parms, GenModel)
   - Test the core logic: probability computation, sequence construction, pruning

4. **Single Realization Models**:
   - Start with models where each event has 1-2 realizations
   - This makes expected outcomes deterministic and testable

5. **Incremental Complexity**:
   - Test isolated events first (with proper mocks)
   - Then test 2-event chains
   - Finally test complete VDJ chains

6. **Use Counters for Scenario Observation**:
   - Custom Counter subclasses are the best tool for capturing complete scenario state
   - Can verify probability accumulation, sequence construction, conditional dependencies
   - Only captures leaf nodes, not intermediate events

7. **Partial Models Have Limitations**:
   - V_choice + V_del alone likely won't work (needs J offsets)
   - VJ models (no D) should work with proper model files
   - For deterministic tests, use complete VDJ with 1-2 realizations per event

---

## Summary Table (Context-Based)

| Event Class      | Loops Over                  | Modifies (Context.field)        | Key Checks              | Downstream Updates    |
|------------------|-----------------------------|--------------------------------|-------------------------|-----------------------|
| **Gene_choice**  | Alignments (`query.gene_alignments`) | `scenario.constructed_sequences` | Overlap safety (VD, VJ) | Error rate bounds     |
|                  |                             | `scenario.seq_offsets`         | Offset feasibility      | Junction length probs |
|                  |                             | `scenario.mismatches_lists`    |                         |                       |
| **Deletion**     | Deletion counts (integer)   | `scenario.constructed_sequences` | Overlap safety updates  | Error rate bounds     |
|                  |                             | `scenario.seq_offsets`         | Palindrome feasibility  | (adjusted for dels)   |
|                  |                             | `scenario.mismatches_lists`    |                         |                       |
| **Insertion**    | (No loop - deterministic)   | `scenario.constructed_sequences` | Insertion count in range| Dinucleotide bounds   |
|                  |                             | (placeholder)                  |                         | Junction length probs |
| **Dinucl_markov**| (Heuristic - single path)   | `scenario.constructed_sequences` | Gene class validity     | Sets probs to 1.0     |
|                  |                             | (fills nucleotides)            |                         | (no more uncertainty) |

**Note**: All events also update `scenario.scenario_proba` and `exploration.index_map` via `iterate_common()`.

---

## Conclusion

The `iterate` method implements a **depth-first recursive exploration** of the recombination scenario tree with **probabilistic pruning**. Each event type has specific logic for:
1. **Enumerating realizations** (alignments, deletion counts, etc.)
2. **Computing probabilities** via marginal array indexing (from `model.model_parameters`)
3. **Constructing sequences** incrementally (in `scenario.constructed_sequences`)
4. **Checking constraints** (overlaps, safety via `exploration.safety_set`)
5. **Updating downstream bounds** for pruning (in `exploration.downstream_proba_map`)
6. **Recursively calling** the next event (via `exploration.next_event_ptr_arr`)

### Context Refactoring Benefits (Phase 1-5, Complete April 2026)

The migration to context-based architecture provides significant improvements:

- **✅ Maintainability**: 5 semantic contexts replace 19 individual parameters
- **✅ Type Safety**: Const-correctness enforced (ModelContext is immutable)
- **✅ Testability**: Easier to mock individual contexts, cleaner test setup
- **✅ Clarity**: Parameter roles explicit (input vs model vs state vs policy vs output)
- **✅ Extensibility**: Can add fields to contexts without breaking all event signatures
- **✅ Consistency**: Error_rate and Counter interfaces match iterate() pattern

### Current State (April 2026)

- **All Rec_Event subclasses**: Migrated to context-based iterate()
- **All Counter subclasses**: Migrated to context-based count_scenario() with Scenario view
- **All Error_rate subclasses**: Migrated to context-based compute_scenario_error_probability()
- **Legacy interfaces**: Preserved for backward compatibility, delegate to context-based implementation
- **Testing**: Should focus on **context-based** interfaces for new tests

Testing should focus on **isolated event logic** with context mocking infrastructure, rather than full I/O or end-to-end inference workflows. The context architecture makes testing significantly easier by providing semantic grouping and clear interfaces.
