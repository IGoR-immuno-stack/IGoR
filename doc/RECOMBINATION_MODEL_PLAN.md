# Recombination Model Refactoring Plan

**Target**: `src/igor/Model/` namespace `igor::model`

## 1. Objective

Replace the monolithic `Model_Parms` with a modular architecture that separates:
1.  **Topology**: The graph structure (nodes, edges, order) managed by `ModelTopology`.
2.  **Parameters**: The probability distributions managed by `RecombinationModel` (via `ParameterHandler`).
3.  **Identity**: Efficient Integer IDs (UIDs) for identifying events, enabling O(1) lookups and cache-friendly vector storage.

## 2. Architecture Principles

### A. The UID Contract
*   Every `Rec_Event` added to the model is assigned a unique `size_t uid`.
*   This `uid` corresponds exactly to the event's index in the internal dense vectors of `ModelTopology` and `RecombinationModel`.
*   **Constraint**: An event instance belongs to one model topology.

### B. Class Separation

#### 1. `ModelTopology`
*   **Responsibility**: Manages the Directed Acyclic Graph (DAG).
*   **Storage**:
    *   `std::vector<std::shared_ptr<Rec_Event>> events_`: Dense storage indexed by UID.
    *   `std::vector<std::vector<size_t>> adjacency_`: Children UIDs for every event.
    *   `std::unordered_map<string, size_t> name_map_`: Name-to-UID resolution.
*   **Algorithms**: Topological sort, cycle detection, ancestry checks.
*   **Performance**: All internal graph traversals use `size_t` indices, not string lookups or pointer chasing.

#### 2. `RecombinationModel`
*   **Responsibility**: The central container for the model state.
*   **Structure**:
    *   Owns `ModelTopology`.
    *   Owns `std::vector<std::unique_ptr<MarginalHandler>> parameters_`: Parallel vector to topology events.
    *   Owns `ErrorRate` model.
*   **Access**:
    *   `get_event(uid)` -> O(1)
    *   `get_handler(uid)` -> O(1)

## 3. Data Data Structures

### `Rec_Event` (Modifications)
Ensure the existing `event_identifier` (or a new strict `uid` field) is used consistently.

```cpp
class Rec_Event {
    size_t uid = -1; // Assigned by ModelTopology
    // ...
public:
    void set_uid(size_t id) { uid = id; }
    size_t get_uid() const { return uid; }
};
```

### `ModelTopology` Interface
```cpp
class ModelTopology {
public:
    using EventId = std::size_t;
    
    // Core Graph Construction
    EventId add_event(std::shared_ptr<Rec_Event> event);
    void add_edge(EventId parent, EventId child);
    
    // Efficient Access
    const Rec_Event& get_event(EventId id) const;
    const std::vector<EventId>& get_children(EventId id) const;
    const std::vector<EventId>& get_parents(EventId id) const; // Maintained for backward traversal
    
    // Algorithms
    const std::vector<EventId>& get_topological_order() const; // Cached
    
private:
    std::vector<std::shared_ptr<Rec_Event>> events_;
    std::vector<std::vector<EventId>> children_;
    std::vector<std::vector<EventId>> parents_; // Optional, useful for views
    std::unordered_map<std::string, EventId> name_to_id_;
};
```

### `RecombinationModel` Interface
```cpp
class RecombinationModel {
public:
    // Lifecycle
    ModelTopology& topology();
    const ModelTopology& topology() const;

    // Parameters (Synced with Topology UIDs)
    void set_handler(EventId id, std::unique_ptr<MarginalHandler> handler);
    MarginalHandler* get_handler(EventId id);

    // Operations
    void normalize_parameters(); // Iterates dense vector
    
private:
    ModelTopology topology_;
    std::vector<std::unique_ptr<MarginalHandler>> handlers_; // handlers_[i] corresponds to topology_.events_[i]
    std::unique_ptr<Single_error_rate> error_rate_;
};
```

## 4. Implementation Steps

### Step 1: `ModelTopology` Implementation
*   Create `src/igor/Model/ModelTopology.h/cpp`.
*   Implement dense vector storage and integer-based adjacency lists.
*   Implement 'Add node' and 'Add edge' logic with cycle detection checks.
*   Implement topological sort (Kahn's or DFS) returning a `vector<EventId>`.

### Step 2: `RecombinationModel` Skeleton
*   Create `src/igor/Model/RecombinationModel.h/cpp`.
*   Implement aggregation of `ModelTopology`.
*   Implement `handlers_` vector management (resize on event addition).

### Step 3: Migration / Interop
*   Crucially, we must allow `Model_Parms` to co-exist or import from this structure initially.
*   Create a factory/converter: `RecombinationModel::from_legacy(const Model_Parms&)` to allow testing the new structure with existing data.

## 5. Benefits vs Current `Model_Parms`

| Feature | `Model_Parms` | New `RecombinationModel` |
| :--- | :--- | :--- |
| **Storage** | `std::list` (linked list) | `std::vector` (contiguous memory) |
| **Lookup** | `O(N)` linear search | `O(1)` direct indexing |
| **Edges** | Separate `Adjacency_list` map | Internal `vector<vector<size_t>>` |
| **Traversal** | Pointer chasing | Cache-friendly integer iteration |
| **Parameters** | Mixed in `Model_marginals` | Decoupled `ParameterHandler` |

