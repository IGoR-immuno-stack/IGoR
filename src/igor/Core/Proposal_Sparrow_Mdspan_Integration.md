# Proposal: Sparrow & mdspan Integration in Core

## 1. Executive Summary
This proposal outlines a strategy to integrate the **Sparrow** (Apache Arrow) library into the **Core** module of IGoR. It builds upon the **Streaming Architecture** defined in `Task_6_Streaming_Architecture.md`, extending it from "Streaming I/O" to **"Zero-Copy Processing"**.

While Task 6 focuses on loading data efficiently, this proposal focuses on **processing** that data without converting it back into expensive legacy C++ structures. We also propose the use of **`std::mdspan`** (C++23 feature, backported to C++20) to provide multi-dimensional views over contiguous Arrow buffers.

Finally, we propose a strategic shift towards **Native AIRR Support**, where IGoR's internal APIs and data models align directly with the AIRR (Adaptive Immune Receptor Repertoire) standard, eliminating the need for constant translation.

## 2. Analysis & Alignment with Task 6

### The "Core" Module (Current)
*   **Data Model:** Relies on **Array of Structures (AoS)** (e.g., `std::vector<Alignment_data>`).
*   **Memory Layout:** High fragmentation, pointer chasing.
*   **Task 6 Status:** Task 6 implements `row_to_sequence_data`, which reads from Arrow but **copies** data into these legacy structures. This incurs a per-row allocation cost (O(N) allocations) and limits the performance benefits of Arrow.

### The "Streaming" Module (Sparrow/Arrow)
*   **Data Model:** **Structure of Arrays (SoA)**.
*   **Memory Layout:** Columnar, contiguous buffers.
*   **Task 6 Status:** Successfully implements `sparrow::record_batch` for I/O.

### The Gap: "Copying" vs "Zero-Copy"
Task 6 achieves O(1) *memory usage* by streaming, but still suffers from O(N) *allocations* because it converts every row into a `SequenceData` struct.

**This Proposal:** Eliminates the `SequenceData` struct entirely for the hot path. Instead of converting `Arrow Row -> C++ Struct`, we use `Arrow Row -> C++ View`.

## 3. Proposed Architecture: "Zero-Copy Proxy Views"

We introduce **Proxy View Objects** in `Core`. These are lightweight C++ objects that *behave* like the current structs but *store no data*.

### Key Benefits
1.  **Readability:** Algorithms look like `sequence.get_gene_name()`.
2.  **Performance:** No data copying. Access is direct to Arrow buffers.
3.  **Thread Safety:** Inherits the read-only thread safety described in Task 6. The `std::atomic` row counter pattern works identically with Views.

### Conceptual Design

**Current (Task 6 Intermediate):**
```cpp
// Allocates memory for string and vector
SequenceData seq = row_to_sequence_data(batch[i], i);
process(seq);
```

**Proposed (View Based):**
```cpp
// No allocations. 'seq' is just a pointer and an index.
SequenceView seq(batch, i);
process(seq);
```

### The `SequenceView` Class
```cpp
class SequenceView {
    const sparrow::record_batch& batch;
    size_t row_index;

public:
    std::string_view get_sequence() const {
        // Zero-copy view into the arrow string buffer
        return batch.get_column("sequence")[row_index];
    }

    // Returns a view, not a vector
    AlignmentView get_alignment(Gene_class gene) const {
        return AlignmentView(batch, row_index, gene);
    }
};
```

## 4. Integration of `std::mdspan`

Task 6 highlights that "All V Alignments" take up **2000 MB** in the current format. `mdspan` solves this by treating the Arrow buffer as a matrix, requiring **0 bytes** of extra memory.

### Use Case A: Replacing `igor::Matrix`
The `Core` module currently uses a custom `Matrix<T>` class. This can be replaced by `mdspan` to wrap any data source (std::vector or Arrow buffer).

### Use Case B: The "Score Grid" (High Performance)
Instead of storing gene scores as separate columns (or a vector of structs), we store them in Arrow as a **FixedSizeList** (a column of arrays).

**In Memory (Arrow Buffer):**
```text
[ V1, D1, J1, V2, D2, J2, V3, D3, J3 ... ]
```

**The `mdspan` View:**
```cpp
// Hypothetical usage in Core
void process_scores(const sparrow::record_batch& batch) {
    // 1. Get raw pointer from Sparrow (Zero-Copy)
    auto* raw_data = batch.get_column("alignment_scores").data<double>();
    
    // 2. Create 2D view (NumRows x 3)
    // No memory allocation here. Just pointer arithmetic.
    auto score_matrix = std::mdspan(raw_data, batch.num_rows(), 3);
    
    // 3. Access
    double s = score_matrix[row_idx, 0]; // V-score
}
```

### Use Case C: `Model_marginals` Optimization
The `Model_marginals` class currently manages a large, flattened 1D array (`marginal_array_smart_p`) that represents a multi-dimensional probability space (events × realizations × dependencies).

*   **Current State:** Manual index arithmetic (`index = base + offset * stride`) is scattered throughout `Model_marginals.cpp`. This is error-prone and hard to read.
*   **Proposed Improvement:** Wrap the underlying `long double*` array with a `std::mdspan`.
    *   **Dimensions:** The dimensions (events, dependencies) are known at runtime based on `Model_Parms`.
    *   **Benefit:** We can replace complex manual indexing with `marginals[event_idx, realization_idx]`.
    *   **Zero-Overhead:** `mdspan` compiles down to the same pointer arithmetic but provides a safer, cleaner API.

## 5. Standardization: Native AIRR Types

To fully align with community standards, we propose that IGoR adopts **AIRR types and naming conventions** as its native internal representation, rather than maintaining a separate "IGoR format" that requires translation.

### 1. Schema Alignment
The in-memory Arrow schema will strictly follow the [AIRR Standards](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).
*   **Old:** `gene_name`, `offset`, `insertions`
*   **New:** `v_call`, `d_call`, `j_call`, `v_sequence_start`, `v_sequence_end`, `v_cigar`

### 2. API Alignment
The `SequenceView` interface will expose AIRR-compliant accessors. This makes the code self-documenting for anyone familiar with the standard.

```cpp
class SequenceView {
public:
    // AIRR Standard Naming
    std::string_view get_sequence_id() const;
    std::string_view get_sequence() const;
    
    // Direct access to AIRR fields
    std::string_view get_v_call() const;  // e.g., "IGHV1-69*01"
    std::string_view get_j_call() const;
    
    // Productivity (AIRR boolean)
    bool productive() const;
};
```

### 3. Type Aliases
We will introduce C++ strong types or aliases that match AIRR definitions to replace generic types.
*   `using SequenceID = std::string_view;`
*   `using GeneCall = std::string_view;` (Replacing internal `Gene_class` enums where appropriate for string data)

### 4. Benefit: "Zero-Translation"
By making the internal data model identical to the external standard:
1.  **Input:** Loading an AIRR TSV/Parquet file requires **zero field mapping**.
2.  **Output:** Writing results is a direct dump of the memory buffers.
3.  **Interoperability:** IGoR can directly link with other AIRR tools (like standard Python libraries) via Arrow C Data Interface.

## 6. Implementation Plan

### Phase 1: Foundation & Dependency
1.  **Update Build System:** Link `Core` to `sparrow`.
2.  **Vendor mdspan:** Add a single-header implementation of `mdspan` to `src/igor/Core/Vendor/mdspan.hpp`.

### Phase 2: The View Abstraction
1.  **`SequenceView` Class:** Implement the view wrapper with **AIRR-compliant naming**.
2.  **`AlignmentView` Class:** Implement the alignment view mapping to AIRR CIGAR/start/end fields.
3.  **Iterator Adapters:** Create iterators compatible with standard algorithms.

### Phase 3: Algorithm Adaptation
1.  **Refactor Aligner:** Modify `Aligner.h` to accept `mdspan<double, 2>` instead of `Matrix<double>`.
2.  **Refactor Model_marginals:** Update `Model_marginals` to use `mdspan` for internal array access.
3.  **Update Inference Loop:** Update `GenModel::infer_model_streaming` to use `SequenceView`.

### Phase 4: Deprecation
1.  Mark `SequenceData` struct as deprecated.
2.  Deprecate non-AIRR field names in public APIs.
