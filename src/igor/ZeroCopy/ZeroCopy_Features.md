# ZeroCopy Module Features

This document summarizes the new features and architectural improvements introduced in the `ZeroCopy` layer of the IGoR project. The primary goal of this module is to enable efficient, zero-copy data access and manipulation using modern C++ standards (specifically `std::mdspan`).

## 1. Architectural Overview

*   **Location**: `src/igor/ZeroCopy/`
*   **Namespace**: `igor` (Directly integrated into the main namespace, removing the nested `zerocopy` namespace).
*   **Code Style**: Enforced **4-space indentation** following the Google C++ Style Guide, maintained via `.clang-format`.

## 2. Core Components & Improvements

### A. `Aligner` Refactoring
The `Aligner` class has been modernized to replace custom memory management with standard containers and views.

*   **`std::mdspan` Integration**: Internal dynamic programming matrices (Score, Memory, Tracker) now use `std::experimental::mdspan`.
*   **Typed Views**: Introduced specific view types for better type safety and clarity:
    ```cpp
    using IntMatrixView = std::experimental::mdspan<int, std::experimental::dextents<size_t, 2>>;
    using DoubleMatrixView = std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>>;
    ```
*   **Zero-Copy Internals**: The `sw_align` method now allocates data in flat `std::vector`s and passes lightweight `mdspan` views to the core alignment logic (`sw_align_common`). This eliminates the need for the legacy `Matrix` class and allows for potential future interoperability with other buffers.

### B. `Model_marginals` Enhancements
The `Model_marginals` class, responsible for storing event probabilities, now exposes its data via structured views.

*   **`MarginalView` (1D)**: Provides a direct, zero-copy view of the entire underlying marginals array.
    ```cpp
    using MarginalView = std::experimental::mdspan<long double, std::experimental::dextents<size_t, 1>>;
    MarginalView get_view() const;
    ```
*   **`EventView` (2D)**: Allows accessing the marginals for a specific event as a 2D matrix (realizations x dependencies), abstracting away the underlying flat array arithmetic.
    ```cpp
    using EventView = std::experimental::mdspan<long double, std::experimental::dextents<size_t, 2>>;
    EventView get_event_view(size_t offset, size_t rows, size_t cols) const;
    ```
*   **Refactored Logic**: Methods like `set_realization_proba` have been updated to use `EventView`, replacing error-prone manual pointer arithmetic with clean, multi-dimensional array access.

### C. Data Views (`ScoreView`, `SequenceView`, `AlignmentView`)
New view classes provide read-only, zero-copy access to structured data formats (e.g., Arrow/Parquet buffers).

*   **`ScoreView`**: Provides a 2D view over alignment scores.
*   **`SequenceView`**: Accesses sequence data without copying strings.
*   **`AlignmentView`**: Navigates alignment records efficiently.

## 3. Technical Benefits

1.  **Zero-Copy Efficiency**: Data is accessed in-place. Large arrays (like marginals or alignment matrices) are not copied when passed between functions; only lightweight views are created.
2.  **Modern C++ Standards**: Leverages `std::mdspan` (C++23 feature, currently via `<experimental/mdspan>`) for multi-dimensional array handling.
3.  **Simplified Memory Management**: Replaces custom, legacy container classes (like `Matrix`) with standard `std::vector` + Views pattern. This reduces technical debt and improves compatibility with standard algorithms.
4.  **Safety & Readability**: `mdspan` provides a standard interface for indexing, reducing the risk of off-by-one errors common in manual pointer arithmetic (e.g., `array[i * cols + j]`).

## 4. Verification

*   **Test Suite**: A comprehensive test suite (`zerocopy_tests`) covers all new components.
*   **Status**: All 45 tests are currently passing, verifying the correctness of the `mdspan` integration and the integrity of the data views.
