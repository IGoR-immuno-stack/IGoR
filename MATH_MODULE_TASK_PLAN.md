# Task Plan: Math Module Implementation

This plan outlines the steps to implement the standalone `Math` module in IGoR, containing the `Tensor` and `Linalg` classes as specified in `TENSOR_LINALG_PROPOSAL_V2.md`.

## 1. Module Structure Setup
*   [x] **Create Directory Structure**
    *   Ensure `src/igor/Math/` exists.
    *   Create `tst/igor/Math/` for unit tests.
*   [x] **Build System Configuration**
    *   Create `src/igor/Math/CMakeLists.txt` defining a header-only library `igor::Math`.
    *   Update `src/igor/CMakeLists.txt` to include the Math subdirectory.
    *   Update `tst/CMakeLists.txt` to discover and build Math tests.

## 2. Implementation: Storage Layer (Complete)
*   [x] **HybridBuffer**
    *   File: `src/igor/Math/HybridBuffer.h`, `src/igor/Math/HybridBuffer.tpp`
    *   Implement `HybridBuffer<T>` with owning/borrowing logic.
    *   Implement move semantics.
    *   **Test**: `tst/igor/Math/test_HybridBuffer.cpp`.
*   [x] **Tensor Class**
    *   File: `src/igor/Math/Tensor.h`, `src/igor/Math/Tensor.tpp`
    *   Implement `Tensor<T>` with dynamic runtime rank.
    *   Implement `view<Rank>()` returning `std::mdspan`.
    *   Implement `ndim()` (was rank).
    *   Implement Accessors: Variadic, Span, Initializer List.
    *   Implement Iterators and Visitor Pattern (`apply`).
    *   **Test**: `tst/igor/Math/test_Tensor.cpp` (Full API coverage).

## 3. Implementation: Algorithm Layer (Linalg)
*   [x] **Basic Arithmetic & Reductions (Phase 2a)**
    *   File: `src/igor/Math/Linalg.h`
    *   Implement Element-wise: `add`, `subtract`, `multiply`, `divide` (on mdspan views).
    *   Implement `scale` (scalar multiplication).
    *   Implement `sum` (full reduction).
    *   Implement `min`, `max`, `argmax`.
    *   **Test**: `tst/igor/Math/test_Linalg_Basic.cpp`
*   [x] **Probabilistic Operations (Phase 2b)**
    *   Implement `normalize` (sum to 1).
    *   Implement `dot` (vector-vector).
    *   Implement `matmul` (matrix-matrix).
    *   **Test**: `tst/igor/Math/test_Linalg_Probabilistic.cpp`.
*   [x] **Log-Space & Stability (Phase 2c)**
    *   Implement `log_add_exp`.
    *   Implement `center` (for numerical stability).
    *   Implement `log_sum_exp` (stable reduction).
    *   Implement `log_normalize` (stable normalization).
    *   **Test**: `tst/igor/Math/test_Linalg_Probabilistic.cpp`.
*   [x] **Broadcasting Support (Phase 2d)**
    *   **Approach**: Minimal design using free functions only (no Tensor modifications).
    *   Implement `broadcast_to` (creates strided view with 0-strides for broadcasting).
    *   Implement `broadcast_multiply` (convenience wrapper combining broadcast + multiply).
    *   **Test**: `tst/igor/Math/test_Linalg_Broadcast.cpp` (98 assertions).
    *   **Features**: Zero-copy broadcasting, cache-efficient, works with existing operations.
## 4. Documentation & Integration
*   [x] **Documentation: README.md**
    *   Created `src/igor/Math/README.md` with features and usage examples.
*   [ ] **Benchmarks**
    *   Create a simple benchmark comparing `Tensor` loops vs `std::vector` raw pointers.
*   [ ] **Integration with Core** (Future)
    *   Refactor `Model_Parms` to use `Tensor`.

## 5. Branch Strategy
*   All work will be performed in the current branch.
*   Commits should be atomic per component.
