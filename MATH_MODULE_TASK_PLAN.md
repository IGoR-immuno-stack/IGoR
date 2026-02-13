# Task Plan: Math Module Implementation

This plan outlines the steps to implement the standalone `Math` module in IGoR, containing the `Tensor` and `Linalg` classes as specified in `TENSOR_LINALG_PROPOSAL_V2.md`.

## 1. Module Structure Setup
*   [x] **Create Directory Structure**
    *   Ensure `src/igor/Math/` exists.
    *   Create `tst/igor/Math/` for unit tests.
*   [x] **Build System Configuration**
    *   Create `src/igor/Math/CMakeLists.txt` defining a header-only library `igor_math`.
    *   Update `src/igor/CMakeLists.txt` to include the Math subdirectory.
    *   Update `tst/CMakeLists.txt` to discover and build Math tests.

## 2. Implementation: Storage Layer
*   [x] **HybridBuffer (Phase 1)**
    *   File: `src/igor/Math/HybridBuffer.h`
    *   Implement `HybridBuffer<T>` with owning/borrowing logic.
    *   Implement move semantics (constructor/assignment).
    *   **Test**: `tst/igor/Math/test_HybridBuffer.cpp` (Verify alloc vs wrap, lifecycle).
*   [x] **Tensor Class (Phase 1)**
    *   File: `src/igor/Math/Tensor.h`
    *   Implement `Tensor<T>` with dynamic runtime rank (`std::vector` shape).
    *   Implement `view<Rank>()` returning `std::mdspan` for optimization.
    *   **Test**: `tst/igor/Math/test_Tensor.cpp` (Verify dynamic shapes, accessors).

## 3. Implementation: Algorithm Layer
*   [ ] **Linear Algebra Operations (Phase 2)**
    *   File: `src/igor/Math/Linalg.h`
    *   Implement Element-wise ops: `element_wise_add`, `multiply`, `divide`, `scale`.
    *   Implement `sum` (reduction).
    *   **Test**: `tst/igor/Math/test_Linalg_Elementwise.cpp`
*   [ ] **Probabilistic Operations (Phase 3)**
    *   Extend `Linalg.h`.
    *   Implement `normalize` (dim-specific).
    *   Implement `marginalize` (tensor contraction).
    *   **Test**: `tst/igor/Math/test_Linalg_Prob.cpp` (validate sums to 1, correct dimension reduction).

## 4. Integration Verification
*   [ ] **Benchmarks**
    *   Create a simple benchmark comparing `Tensor` loops vs `std::vector` raw pointers to ensure <5% overhead (optional).

## 5. Branch Strategy
*   All work will be performed in the current branch.
*   Commits should be atomic per component (Buffer, Tensor, Linalg).
