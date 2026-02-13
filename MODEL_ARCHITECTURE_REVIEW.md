# Model Architecture Review & Refactoring Plan

## 1. Current Architecture Analysis

The current model architecture relies on three tightly coupled components:

*   **`Model_Parms` (Structure)**: Defines the graph topology (DAG).
*   **`Model_marginals` (Data)**: Stores all probability tables (Condition Probability Distributions - CPDs) in a single contiguous `double` array ("Arena Allocation").
*   **`Error_rate` (Observation)**: Handles sequence comparison and mismatch penalization as a special-case "leaf callback".

### Identified Issues

1.  **The "Giant Array" Pattern (`Model_marginals`)**:
    *   **Fragility**: Access requires manual pointer arithmetic (`base_index + i * distinct_offset`). A single calculation error leads to silent data corruption (reading neighbor event's data).
    *   **Coupling**: `Rec_Event` must ask `Model_marginals` for offsets to navigate its own logical data structure.
    *   **Opacity**: Debugging is difficult; the array is a soup of numbers without boundaries.

2.  **Implicit Schema**:
    *   The "shape" of the data (dimensions of the CPDs) is implicit in the execution order of `Rec_Event::iterate`. There is no explicit object describing the tensor shape for a given event.

3.  **Observation as a Special Case (`Error_rate`)**:
    *   `Error_rate` logic is invoked explicitly at the end of recursion.
    *   **Optimization Barrier**: Mismatches cannot be evaluated until the full sequence is constructed, preventing "Early Pruning" (discarding bad branches as soon as they diverge from observation).
    *   **Leaky Abstraction**: The core recursion loop knows about "mismatches" and "offsets", concepts that belong to the observation model, not the graph traversal.

## 2. Proposed Architecture

### A. Data Layout: Split & Typed Containers
Instead of one monolithic array, `Model_marginals` will act as a registry of **Event Tables**.

*   **Structure**: `std::map<std::string, std::vector<double>> event_tables;`
*   **Logic**: Each event owns its own contiguous block of memory. Bounds checking becomes possible.

### B. Access Pattern: Multi-Dimensional Views (`std::mdspan`)
We will replace raw pointer arithmetic with modern C++ views.

*   **Concept**: Use `std::mdspan` (C++23) or a lightweight backport/implementation to provide a multi-dimensional view over the flat `std::vector`.
*   **Usage**:
    ```cpp
    // Old
    double p = array[base + parent_val * offset + child_val];
    // New
    double p = event_view(parent_val, child_val);
    ```

### C. Unified Topology: Observation as a Node
The "Error Rate" logic will be encapsulated into a standard Graph Node (e.g., `ObservationEvent`).

*   **Topology**: The graph will formally include an `ObservationEvent` as a child of the final sequence construction events.
*   **Recursion**: The `iterate` loop naturally flows into this node. It computes the mismatch probability and returns it, just like any other probabilistic event.
*   **Optimization**: This enables "Intermediate Observation Nodes" (e.g., checking V-gene matches immediately after V-gene choice) for massive performance gains via early connection pruning.

### E. Algebraic Capabilities
`Model_marginals` is not just storage; it is a mathematical object used in the EM algorithm. The refactored design must support:

*   **Thread-Safe Accumulation**: Efficient `operator+=` to merge thread-local arrays into global counters.
*   **Element-wise Operations**: Normalization (scalar division) and scaling.
*   **Tensor Arithmetic**: Support for the marginalization and product operations required by Bayesian Inversion (e.g., `invert_edge`).

## 3. Refactoring Plan


This refactoring is a prerequisite for the broader `Rec_Event` modernization.

### Phase 1: Safe Storage (Internal Refactor)
**Goal**: Break the "Giant Array" without changing external APIs significantly.

1.  **Modify `Model_marginals`**:
    *   Replace `unique_ptr<double[]> marginal_array_smart_p` with `unordered_map<string, vector<double>> storage`.
    *   Implement `get_event_data_pointer(name)` to return `vector.data()`.
2.  **Update `Rec_Event` Initialization**:
    *   Pass `local_ptr` (0-based) instead of `global_ptr + offset`.
    *   Update initialization logic to request data by name.
3.  **Validation**:
    *   Ensure exact numerical output match with current version.

### Phase 2: Structural Descriptors & Views
**Goal**: Explicitly define the shape of probabilities.

1.  **Define `EventShape`**:
    *   Create a struct describing dimensions: `shape = {parent1_size, parent2_size, self_size}`.
2.  **Integrate `mdspan` (or equivalent)**:
    *   Update `Rec_Event::iterate` to accept a View object instead of raw pointers.
    *   Refactor internal math to use `view(i, j, k)`.

### Phase 3: The "Observation Node" (Topology Change)
**Goal**: Unify `Error_rate` into the graph.

1.  **Create `class ObservationEvent : public Rec_Event`**:
    *   Move `compare_sequences_error_prob` logic into `ObservationEvent::iterate`.
2.  **Update `Model_Parms`**:
    *   Automatically append this node to the graph during model loading.
3.  **Simplify Recursion**:
    *   Remove `if (leaf) { error_rate->... }` from the generic traversal code.

## 4. Immediate Next Step
Begin **Phase 1**: implementing the split storage in `Model_marginals`. This provides immediate safety benefits and paves the way for `mdspan` without requiring a "Big Bang" rewrite of the algorithm.
