# Comparison: `igor::math::Tensor` vs `xtensor::xarray`

This document compares the design of the newly implemented `igor::math::Tensor<T>` with `xtensor::xarray`, a popular C++ multi-dimensional array library found in the workspace.

## 1. Core Architecture: Dynamic Rank

Both implementations share the fundamental design choice of **handling dimensions at runtime**.

*   **`xtensor::xarray`**:
    *   Defines `rank = SIZE_MAX` (dynamic).
    *   Stores shape and strides in dynamic containers (e.g., `std::vector<size_t>`).
*   **`igor::math::Tensor`**:
    *   Defines dynamic rank via `std::vector<size_t> dims`.
    *   **Verdict**: Identical architectural approach. Both solve the problem of loading model parameters where the number of dimensions is not known at compile time.

## 2. Storage Strategy: Unification vs. Separation

*   **`xtensor::xarray`**:
    *   Separates ownership concepts into two distinct types:
        1.  `xarray_container` (Owning): Manages memory (like `std::vector`).
        2.  `xarray_adaptor` (Borrowing): Wraps existing memory buffers.
    *   **Pros**: Compile-time strictness.
    *   **Cons**: Requires templating code on the container type to handle both owning and borrowing arrays genericly.

*   **`igor::math::Tensor`**:
    *   **Unifies** both concepts in a single class using `HybridBuffer`.
    *   Uses a runtime boolean `is_owning_` to switch behavior.
    *   **Pros**: `Model_marginals` can store `std::map<string, Tensor>` containing a mix of owned and borrowed arrays without complex type erasure or variants.
    *   **Cons**: Slight runtime overhead (one branch in destructor).

## 3. Computation Model: Explicit Views vs. Expression Templates

*   **`xtensor::xarray`**:
    *   provides a "Rich API" (`operator+`, `broadcast`, etc.) on the container itself.
    *   Uses **Expression Templates** (Lazy Evaluation) to optimize composed operations.
    *   Supports runtime-rank indexing `arr(i, j, k)` (slower/generic).

*   **`igor::math::Tensor`**:
    *   Acts as a **"Dumb Container"**. It has no math operators.
    *   Forces the use of **Explicit Views** (`view<Rank>()`) to perform computation.
    *   **Why this fits IGoR**:
        *   IGoR's critical loops (marginalization) are deeply nested.
        *   By forcing `view<N>()`, we generate `mdspan<..., N>` which allows the compiler to unroll/vectorize loops as if variables were static arrays.
        *   It prevents accidental slow generic access.

## Summary

The `igor::math::Tensor` is essentially a **"Minimalist, Type-Eased xarray"**.

| Feature | `xtensor::xarray` | `igor::math::Tensor` |
| :--- | :--- | :--- |
| **Rank** | Dynamic | Dynamic |
| **Storage** | Template Policy (Owner/Adaptor) | Runtime Hybrid |
| **Philosophy** | Full Algebra Library | Raw Data Container |
| **Performance** | Expression Templates | Explicit `mdspan` Loops |

Your intuition was correct: we have re-implemented the core storage mechanics of `xarray`, but stripped away the complexity of expression templates in favor of standard C++23 `mdspan` views for optimization.
