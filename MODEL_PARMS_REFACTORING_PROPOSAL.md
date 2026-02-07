# Model_Parms Refactoring Proposal

This document outlines specific architectural improvements for the `Model_Parms` class, focusing on structural definitions and graph optimization.

## 1. Structural Definition: Typed DAG

The current graph implementation in `Model_Parms` relies heavily on string hashing (`Rec_Event_name`) and ad-hoc traversal logic, which impacts performance and type safety.

### Proposed Architecture

*   **Integer IDs**: Map string names to `EventID` (uint32_t) at load time.
    *   Internal storage should use flat vectors indexed by ID (e.g., `std::vector<Event> nodes`).
    *   This eliminates repetitive string hashing during graph traversal and inference.
*   **Explicit Graph Structure**: Replace `unordered_map<string, AdjacencyList>` with a dedicated DAG structure.
    *   Optimizes parent/child lookups to simple vector indexing.
    *   Guarantees topological validity by construction.
*   **Immutability**: Once the model is loaded, the structure should ideally become read-only.
    *   Allows concurrent readers (e.g., parallel inference threads) without complex locking mechanisms.

## 2. Refactoring Plan

### Phase 4: Graph Optimization (`Model_Parms`)

**Goal**: Speed up topology queries and improve type safety.

1.  **Indexer System**:
    *   Introduce an internal `EventID` system to replace raw pointers and strings in hot paths.
    *   Maintain a bidirectional `String <-> ID` registry for IO and debugging.

2.  **Vectorization**:
    *   Replace internal `std::unordered_map` storage with `std::vector` indexed by `EventID`.
    *   This improves cache locality and lookup speed.

3.  **Robust Custom Implementation**:
    *   Implement standard graph operations (cycle detection, topological sort) directly on the vector structure, avoiding external dependencies.
