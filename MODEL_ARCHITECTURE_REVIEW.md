# Math Module & Model Marginals Modernization Review

## 1. Math Module Summary

The **Math Module** is a high-performance, C++23-compliant library designed as the new mathematical foundation for IGoR. It provides N-dimensional tensor operations, linear algebra primitives, and zero-copy data handling, serving as a direct replacement for the legacy "giant array" and raw pointer arithmetic currently used in `Model_Marginals`.

The module is built on three core pillars:

### Core Pillars

1.  **`Tensor<T>` & `HybridBuffer`**
    *   A unified container that can either **own** its memory (new allocation) or **borrow** existing memory (zero-copy view).
    *   This allows it to wrap legacy data structures without duplication, facilitating a smooth transition.

2.  **`std::mdspan` Architecture**
    *   All operations work on structural views (`std::mdspan`) rather than specific containers.
    *   This decouples algorithms from data storage, following modern C++26 design principles (P1673).

3.  **`Linalg` Suite**
    *   A collection of stateless, optimized algorithms for:
        *   **Element-wise ops**: Addition, subtraction, multiplication, division, scaling.
        *   **Reductions**: Full tensor summation, min/max, argmax.
        *   **Probabilistic ops**: Dot product, matrix multiplication, full tensor normalization ($\sum P = 1$).
        *   **Log-space stability**: `log_add_exp`, `log_sum_exp`, `log_normalize` for numerically stable inference.
        *   **Broadcasting**: Zero-copy dimension expansion via stride manipulation (NumPy-style).

---

## 2. The Link: Modernizing `Model_Marginals`

The `Model_Marginals` class currently suffers from high complexity due to manual memory management and "unsafe" pointer arithmetic. The Math module explicitly solves these issues, providing a path to a safer, more maintainable, and parallelizable architecture.

### Comparison: Legacy vs. Modern Approach

| Feature | Current `Model_Marginals` Approach | Modernized `Math` Approach |
| :--- | :--- | :--- |
| **Storage** | `std::unique_ptr<long double[]> marginal_array_smart_p`<br>*(A single monolithic 1D array)* | **`Tensor<long double>`**<br>*(Structured N-dimensional container with runtime rank)* |
| **Access** | Manual index calculation:<br>`array[base + offset_i + j]` | **Direct Multi-dim Indexing**:<br>`tensor(event_idx, i, j)` or `tensor.view()` |
| **Safety** | Raw pointer arithmetic (prone to segfaults & off-by-one errors) | **`std::mdspan` Views**<br>*(Type-safe, bounds-checked in debug, zero-overhead in release)* |
| **Logic** | Nested `for` loops interspersed with logic | **`linalg::normalize` / `linalg::broadcast_multiply`**<br>*(Expressive, optimized standard algorithms)* |
| **Stability** | `long double` (80-bit) dependence | **Log-Space Support**<br>*(Generic math supporting `double` with `log_add_exp`, `log_normalize`)* |
| **Parallelization** | Single-threaded giant array | **Thread-Local Tensors**<br>*(Lock-free parallel accumulation via `operator+=`)* |

---

## 3. Implementation Pathway

### Current Status: Math Module Complete (268 tests passing)

**Implemented:**
- ✅ Tensor<T> with runtime rank (1-5 dimensions)
- ✅ HybridBuffer with owning/borrowing semantics
- ✅ Element-wise operations (add, subtract, multiply, divide, scale)
- ✅ Reductions (sum, min, max, argmax)
- ✅ Probabilistic ops (dot, matmul, normalize - full tensor only)
- ✅ Log-space stability (log_add_exp, log_sum_exp, log_normalize)
- ✅ Broadcasting (broadcast_to, broadcast_multiply)

**Critical Gap for Model_Marginals Migration:**
- ❌ **Axis-specific operations** (normalize along dimension, sum along axis)
- ❌ **Slicing operations** (extract sub-tensors along specific dimensions)

### Migration Strategy (5 Phases)

The transition requires careful phasing due to Model_Marginals complexity (956 lines, recursive event structure):

#### **Phase 0: Fill Math Module Gaps** (2-3 weeks)
*Pre-requisite before any migration*

1.  **Implement axis-specific reductions**
    *   `sum_axis(tensor, axis)` - sum along specific dimension
    *   `normalize_axis(tensor, out, axis)` - normalize along dimension
    *   Critical for per-event normalization in Model_Marginals

2.  **Implement slicing/indexing**
    *   `slice(tensor, dim, index)` - extract sub-tensor
    *   Needed for navigating complex event hierarchies

#### **Phase 1: Borrowing Wrapper** (1 week)
*Zero risk - wraps existing array*

1.  **Create Model_marginals_v2 class**
    *   Keep existing `marginal_array_smart_p` storage
    *   Add `as_tensor()` method returning borrowing Tensor view
    *   Implement new methods delegating to Math module
    *   Old methods continue working unchanged

#### **Phase 2: Replace Storage** (2-3 weeks)
*Structural change - high risk*

1.  **Migrate to structured storage**
    *   Replace giant array with `std::unordered_map<Rec_Event_name, Tensor<long double>>`
    *   Maintain offset maps for backward compatibility
    *   Update constructors and copy operations

#### **Phase 3: Refactor Operations** (3-4 weeks)
*Incremental replacement of logic*

1.  **Simple operations first**
    *   `operator+=` using `linalg::add`
    *   `operator-=` using `linalg::subtract`
    *   `add_pseudo_counts` using `linalg::scale`

2.  **Complex normalization**
    *   Replace recursive `iterate_normalize` with `linalg::normalize_axis`
    *   Handle parent-child event dependencies
    *   Maintain numerical equivalence

#### **Phase 4: Parallelization** (2-3 weeks)
*Performance enhancement*

1.  **Thread-local accumulation**
    *   Create thread-local `Model_marginals_v2` per thread
    *   Parallel inference loop with OpenMP
    *   Final reduction using `operator+=`

2.  **Broadcasting optimizations**
    *   Use `broadcast_multiply` for marginal scaling
    *   Parallel processing of independent events

#### **Phase 5: Validation** (1-2 weeks)
*Ensure correctness*

1.  **Numerical equivalence tests**
    *   Compare v1 vs v2 on benchmark datasets
    *   Verify marginals match within tolerance

2.  **Performance benchmarking**
    *   Measure inference throughput
    *   Assess parallelization scaling

**Total Estimated Timeline:** 3-4 months

**Critical Dependencies:**
- Phase 0 must complete before Phase 1
- Phase 2 is high-risk and requires extensive testing
- Phase 4 depends on Phase 3 completion
