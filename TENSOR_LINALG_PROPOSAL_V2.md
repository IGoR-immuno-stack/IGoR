# Tensor and Linear Algebra Module Proposal v2

**Date:** 6 February 2026
**Target:** IGoR Bayesian inference infrastructure
**Scope:** Minimal viable linear algebra module with future extensibility

---

## 1. Design Philosophy

### Core Principles

1. **Modern C++ Standards**
   - C++23 `std::mdspan` for multi-dimensional views
   - C++20 concepts for compile-time interface validation
   - Standard library only (no external dependencies)

2. **Algorithm-Data Separation**
   - Algorithms operate on views (`std::mdspan`), not containers
   - Follows C++26 `std::linalg` proposal (P1673) design
   - Enables operations on owned buffers, borrowed memory, or sub-slices equally

3. **Zero-Copy Flexibility**
   - Hybrid storage model: owning vs. borrowing semantics
   - Solves `std::vector`'s ownership-only limitation
   - Critical for IGoR's event-based marginal arrays

4. **Minimal Surface Area**
   - 8 operations total (down from theoretical 20+)
   - Balanced: current needs (7 ops) + future-proofing (1 op)
   - Excludes: broadcasting, matrix algebra, multi-axis marginalization

---

## 2. Storage Layer

### 2.1 HybridBuffer<T>

A simple ownership-aware pointer wrapper solving the "borrow or own" problem.

```cpp
template<typename T>
class HybridBuffer {
    T* data_;
    size_t size_;
    bool is_owning_;

public:
    // Standard container type aliases
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;

    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    using iterator = pointer;
    using const_iterator = const_pointer;

    // Owning constructor: allocates memory
    explicit HybridBuffer(size_t n)
        : data_(new T[n]), size_(n), is_owning_(true) {}

    // Borrowing constructor: wraps existing memory
    HybridBuffer(T* ptr, size_t n)
        : data_(ptr), size_(n), is_owning_(false) {}

    ~HybridBuffer() {
        if (is_owning_) delete[] data_;
    }

    // Move semantics for ownership transfer
    HybridBuffer(HybridBuffer&& other) noexcept;
    HybridBuffer& operator=(HybridBuffer&& other) noexcept;

    // Delete copy (explicit ownership policy)
    HybridBuffer(const HybridBuffer&) = delete;
    HybridBuffer& operator=(const HybridBuffer&) = delete;

    // Data access
    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }

    // Iterator interface for STL compatibility
    iterator begin() { return data_; }
    iterator end() { return data_ + size_; }
    const_iterator begin() const { return data_; }
    const_iterator end() const { return data_ + size_; }
    const_iterator cbegin() const { return data_; }
    const_iterator cend() const { return data_ + size_; }

    // Element access
    reference operator[](size_type pos) { return data_[pos]; }
    const_reference operator[](size_type pos) const { return data_[pos]; }
};
```

**Design Rationale:**
- **Ownership flag** determines destructor behavior
- Prevents accidental double-free (move-only semantics)
- Zero overhead: same memory layout as raw pointer + size

---

### 2.2 Tensor<T, Rank>

An N-dimensional array container managing memory lifecycle and shape metadata.

```cpp
template<typename T, size_t Rank>
class Tensor {
    HybridBuffer<T> storage_;
    std::dextents<size_t, Rank> dims_;  // Critical: std::dextents, not std::array

public:
    // Standard container type aliases
    using self_type = Tensor<T, Rank>;
    using storage_type = HybridBuffer<T>;
    using value_type = T;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;

    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    using shape_type = std::dextents<size_type, Rank>;
    using extent_type = typename shape_type::index_type;

    // Iterator types (flat iteration over storage)
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    static constexpr size_type rank = Rank;

    // Owning constructor: allocate memory for given shape
    explicit Tensor(std::dextents<size_t, Rank> extents)
        : storage_(compute_size(extents)),
          dims_(extents) {}

    // Borrowing constructor: wrap existing memory
    Tensor(T* data, std::dextents<size_t, Rank> extents)
        : storage_(data, compute_size(extents)),
          dims_(extents) {}

    // View accessor: returns std::mdspan for indexing
    auto view() {
        return std::mdspan<T, std::dextents<size_t, Rank>>(
            storage_.data(), dims_
        );
    }

    auto view() const {
        return std::mdspan<const T, std::dextents<size_t, Rank>>(
            storage_.data(), dims_
        );
    }

    // Direct accessors for algorithm implementations
    T* data() { return storage_.data(); }
    const T* data() const { return storage_.data(); }
    const auto& shape() const { return dims_; }
    const auto& extents() const { return dims_; }

    // Size queries
    size_type size() const { return storage_.size(); }
    size_type dimension() const { return Rank; }

    // Flat iteration over storage (for STL algorithms)
    iterator begin() { return storage_.begin(); }
    iterator end() { return storage_.end(); }
    const_iterator begin() const { return storage_.begin(); }
    const_iterator end() const { return storage_.end(); }
    const_iterator cbegin() const { return storage_.cbegin(); }
    const_iterator cend() const { return storage_.cend(); }

    // Multi-dimensional indexing: forwards to view()
    template<typename... Indices>
    reference operator()(Indices... indices) {
        return view()(indices...);
    }

    template<typename... Indices>
    const_reference operator()(Indices... indices) const {
        return view()(indices...);
    }

private:
    static size_t compute_size(const std::dextents<size_t, Rank>& e) {
        size_t product = 1;
        for (size_t i = 0; i < Rank; ++i) {
            product *= e.extent(i);
        }
        return product;
    }
};
```

**Key Design Decisions:**

1. **std::dextents<size_t, Rank>** instead of `std::array<size_t, Rank>`
   - Supports `std::dynamic_extent` for runtime-determined rank
   - Example: `Tensor<double, std::dynamic_extent>` for flexible dimensions
   - Critical for `Model_marginals` where event counts vary

2. **view() returns mdspan**
   - Decouples storage from indexing semantics
   - Enables algorithms to work on sub-tensors via `std::submdspan`
   - Natural integration with future C++26 `std::linalg`

3. **operator() for convenience**
   - Direct indexing: `tensor(i, j, k)` instead of `tensor.view()(i, j, k)`
   - Forwards to `view()` - zero overhead
   - Matches xtensor convention for intuitive use

4. **No built-in indexing complexity**
   - Both interfaces coexist: `view()` for algorithms, `operator()` for convenience
   - Compiler inlines both in release builds
   - Explicit about storage vs. view separation

---

## 3. Algorithm Layer: `igor::linalg`

All algorithms operate on `std::mdspan` views, following P1673 design principles.

### 3.1 Operation Catalog

#### **Element-wise Operations (4 ops)**

```cpp
namespace igor::linalg {

// Addition: out[i] = a[i] + b[i]
template<typename InView1, typename InView2, typename OutView>
void element_wise_add(InView1 a, InView2 b, OutView out) {
    // Requires: compatible shapes
    // Supports: in-place (out may alias a or b)
}

// Multiplication: out[i] = a[i] * b[i] (Hadamard product)
template<typename InView1, typename InView2, typename OutView>
void element_wise_multiply(InView1 a, InView2 b, OutView out) {
    // Use case: joint probabilities P(A,B) = P(A) * P(B)
}

// Division: out[i] = a[i] / b[i]
template<typename InView1, typename InView2, typename OutView>
void element_wise_divide(InView1 a, InView2 b, OutView out) {
    // Use case: Bayesian inversion (Model_marginals.cpp:340)
    // Handles: zero division (sets result to 0.0)
}

// Scalar multiplication: out[i] = in[i] * alpha
template<typename InView, typename Scalar, typename OutView>
void scale(InView in, Scalar alpha, OutView out) {
    // Use case: normalization scaling
    // Supports: in-place
}

} // namespace igor::linalg
```

---

#### **Reductions (2 ops)**

```cpp
// Full reduction: returns sum of all elements
template<typename InView>
auto sum(InView in) -> typename InView::value_type {
    // Returns: scalar
    // Use case: compute normalization factors
}

// Tensor contraction: sum over trailing dimension
template<typename InView, typename OutView>
void marginalize(InView in, OutView out, size_t dim) {
    // Semantics: sums over dimension `dim` (eliminates it)
    // Example: (64, 4, 3) --[dim=2]--> (64, 4)
    // Use case: P(V, D) = Σ_J P(V, D, J)
}
```

**Design Note on Marginalize:**

Single-dimension only (not multi-axis). Rationale:
- IGoR's `Model_marginals::normalize()` already handles multi-dimensional via recursion
- Multi-axis can be composed: `marginalize(temp, out1, 2); marginalize(out1, out2, 1);`
- Simpler API, sufficient for current patterns

---

#### **Probability Operations (1 op)**

```cpp
// Normalization: enforces Σᵢ P(i) = 1 along dimension
template<typename InView, typename OutView>
void normalize(InView in, OutView out, size_t dim) {
    // Implementation:
    // 1. Compute sum along `dim`
    // 2. Divide each slice by its sum
    // Supports: in-place
    // Handles: zero sums (leaves slice unchanged)
}
```

**Mathematical Definition:**

$$\text{normalize}(P, \text{dim}=d) \implies \sum_{i_d} P[\ldots, i_d, \ldots] = 1$$

**Example:**
```cpp
// Marginals shape: (64 genes, 4 deletions)
normalize(marginals.view(), result.view(), 1);
// Result: each gene's deletion probabilities sum to 1
```

---

#### **Numerical Stability (1 op, future extension)**

```cpp
// Log-space addition: log(exp(a) + exp(b))
template<typename Scalar>
Scalar log_add_exp(Scalar log_a, Scalar log_b) {
    if (log_a > log_b) {
        return log_a + std::log1p(std::exp(log_b - log_a));
    } else {
        return log_b + std::log1p(std::exp(log_a - log_b));
    }
}
```

**Mathematical Definition:**

$$\text{log_add_exp}(a, b) = \max(a, b) + \log(1 + e^{-|a - b|})$$

**Rationale:**
- IGoR currently uses `long double` (80-bit precision), not log-space
- Future-proofing for extreme probabilities ($P < 10^{-100}$)
- Enables underflow-free accumulation: $\log(P_1 + P_2 + \ldots)$

**Status:** Optional for v1.0, include for extensibility

---

### 3.2 Design Specifications

#### **Memory Layout**

**Default:** Row-major (C-style, `std::layout_right`)

```cpp
// Tensor shape (64, 4): 64 genes × 4 deletions
// Memory layout: [gene₀_del₀, gene₀_del₁, gene₀_del₂, gene₀_del₃, gene₁_del₀, ...]
// Index formula: A[i,j] = data[i * 4 + j]
```

**Rationale:**
- Matches C++ convention
- Optimal for IGoR's access patterns (iterate genes, then deletions)
- Marginalization over trailing dimension is cache-friendly

**Extensibility:**
```cpp
// Future support for column-major if needed
template<typename T, size_t Rank, typename Layout = std::layout_right>
class Tensor { /* ... */ };
```

---

#### **In-Place Operations**

**All operations support in-place computation** when output aliases an input:

```cpp
element_wise_add(a.view(), b.view(), a.view());  // a += b
scale(data.view(), 2.0, data.view());            // data *= 2
normalize(x.view(), x.view(), 1);                // In-place normalization
```

**Implementation Requirement:**
- Read elements before writing (avoid read-after-write hazards)
- Example: `for (size_t i = 0; i < n; ++i) { out[i] = in[i] * scale; }`

---

#### **Dimension Semantics**

**Convention:** `dim` parameter specifies dimension to **eliminate** (sum over), following NumPy semantics.

```cpp
// P(V, D, J) has shape (64, 4, 3)
marginalize(joint.view(), marginal.view(), 2);
// Result: P(V, D) = Σ_J P(V, D, J), shape (64, 4)
```

**Rationale:** Matches `numpy.sum(axis=...)` behavior, widely understood convention.

---

## 4. Integration with IGoR

### 4.1 Model_marginals Refactoring

**Current (unsafe):**
```cpp
class Model_marginals {
    std::unique_ptr<long double[]> marginal_array_smart_p;  // Line 114
    size_t marginal_arr_size;

    // Unsafe raw pointer arithmetic everywhere
    void iterate(...) {
        marginal_array_smart_p[base_index + offset + iter] *= factor;
    }
};
```

**Proposed (safe):**
```cpp
class Model_marginals {
    std::map<Rec_Event_name, Tensor<long double, 2>> event_marginals;

    void add_to_marginals(Rec_Event& event, /* ... */) {
        auto view = event_marginals[event.name].view();
        linalg::element_wise_add(view, event_probs, view);  // Clear, safe
    }

    void normalize_event(Rec_Event_name name) {
        auto& tensor = event_marginals[name];
        linalg::normalize(tensor.view(), tensor.view(), 1);  // Explicit dimension
    }
};
```

**Benefits:**
- Type-safe indexing via `mdspan`
- Bounds checking in debug builds
- Clear algorithmic intent (`linalg::normalize` vs nested loops)
- Per-event tensors (instead of monolithic array)

---

### 4.2 Thread Safety Enhancement

**Current Problem:** Mutable state in `Rec_Event` prevents parallelization.

**Tensor-Based Solution:**
```cpp
// Phase 0: Extract state to InferenceContext
struct InferenceContext {
    Tensor<long double, 2> thread_local_marginals;  // One per thread
    double scenario_proba;
    Index_map base_index_map;
};

// Parallel inference
#pragma omp parallel for
for (size_t seq = 0; seq < sequences.size(); ++seq) {
    InferenceContext ctx = create_context();
    event.add_to_marginals(ctx, sequences[seq]);

    // Reduction phase
    #pragma omp critical
    linalg::element_wise_add(
        ctx.thread_local_marginals.view(),
        global_marginals.view(),
        global_marginals.view()
    );
}
```

**Impact:** Enables per-thread marginal accumulation without locks.

---

### 4.3 Numerical Precision Considerations

**IGoR's Current Choice:** `long double` (80-bit extended precision)

**Tensor Support:**
```cpp
using Probability = long double;
Tensor<Probability, 2> marginals({64, 4});  // Full precision
```

**Future Option:** Log-space for extreme probabilities
```cpp
Tensor<double, 2> log_marginals({64, 4});
// Accumulate: log_marginals[i] = log_add_exp(log_marginals[i], log_p);
```

**Recommendation:** Start with `long double`, add `log_add_exp` as extension point.

---

## 5. Implementation Roadmap

### Phase 1: Core Infrastructure (Week 1)
- [ ] Implement `HybridBuffer<T>`
- [ ] Implement `Tensor<T, Rank>` with `std::dextents`
- [ ] Unit tests: owning/borrowing semantics, move operations

### Phase 2: Element-wise Operations (Week 2)
- [ ] `element_wise_add`, `element_wise_multiply`, `element_wise_divide`
- [ ] `scale`
- [ ] Unit tests: in-place support, shape validation
- [ ] Benchmark: vs raw pointer arithmetic

### Phase 3: Reductions (Week 3)
- [ ] `sum` (full reduction)
- [ ] `marginalize` (single dimension)
- [ ] Unit tests: multi-dimensional tensors, edge cases
- [ ] Validate: matches `Model_marginals::iterate_normalize` behavior

### Phase 4: Probability Operations (Week 4)
- [ ] `normalize` (with dimension parameter)
- [ ] Unit tests: zero-sum handling, in-place normalization
- [ ] Integration test: replace one `Model_marginals` method

### Phase 5: Integration & Validation (Week 5)
- [ ] Create `TensorMarginals` class (parallel to `Model_marginals`)
- [ ] Port one event's `add_to_marginals` implementation
- [ ] Numerical validation: compare results with original implementation
- [ ] Performance benchmark: tensor API vs raw pointers

### Phase 6: Future Extensions (Post-v1.0)
- [ ] `log_add_exp` for numerical stability
- [ ] Column-major layout support (`std::layout_left`)
- [ ] SIMD optimizations for contiguous operations
- [ ] Multi-axis marginalization (if needed)

---

## 6. Performance Considerations

### 6.1 Contiguous Layout Optimization

For contiguous `mdspan` views, algorithms reduce to pointer arithmetic:

```cpp
template<typename InView, typename OutView>
void scale(InView in, double alpha, OutView out) {
    if constexpr (InView::is_always_contiguous()) {
        // Fast path: single loop
        const double* in_ptr = in.data_handle();
        double* out_ptr = out.data_handle();
        size_t n = in.size();
        for (size_t i = 0; i < n; ++i) {
            out_ptr[i] = in_ptr[i] * alpha;
        }
    } else {
        // Generic path: mdspan indexing
        for (size_t i = 0; i < in.extent(0); ++i) {
            for (size_t j = 0; j < in.extent(1); ++j) {
                out(i, j) = in(i, j) * alpha;
            }
        }
    }
}
```

**Expected Overhead:** <5% vs hand-written loops in release builds (`-O3`).

---

### 6.2 Cache Efficiency

**Row-major marginalization over trailing dimension:**
```cpp
// marginals[gene, deletion] has shape (64, 4)
// Sum over deletions (dim=1): cache-friendly stride-1 access
for (size_t gene = 0; gene < 64; ++gene) {
    double sum = 0.0;
    for (size_t del = 0; del < 4; ++del) {  // Contiguous in memory
        sum += marginals[gene, del];
    }
    result[gene] = sum;
}
```

**Cache Line Utilization:** ~100% for contiguous operations.

---

### 6.3 Compiler Vectorization

Modern compilers auto-vectorize simple loops on `mdspan`:

```cpp
// This loop vectorizes well (AVX2/AVX-512)
for (size_t i = 0; i < n; ++i) {
    out[i] = a[i] + b[i];
}
```

**Expected SIMD Speedup:** 2-4× on modern CPUs for element-wise operations.

---

## 7. Testing Strategy

### Unit Tests (per operation)
- **Shape validation:** Mismatched dimensions throw exceptions
- **Boundary conditions:** Empty tensors, single-element tensors
- **Numerical accuracy:** Compare with reference implementations
- **In-place correctness:** Verify aliased input/output works

### Integration Tests
- **Model_marginals equivalence:** Tensor-based results match original
- **Thread safety:** Parallel accumulation produces correct sums
- **Normalization:** Probability constraints satisfied ($\sum P = 1$)

### Benchmarks
- **Raw pointer baseline:** Measure overhead of `mdspan` abstraction
- **Memory bandwidth:** Verify cache-efficient access patterns
- **Scalability:** Performance with increasing tensor dimensions

---

## 8. What We're NOT Doing

To maintain minimal scope, these features are **explicitly excluded:**

### ❌ Broadcasting
- **Reason:** IGoR pre-aligns arrays to matching shapes
- **Example excluded:** `(64, 4) * (4,)` with shape broadcasting
- **Alternative:** Explicit reshaping if needed in future

### ❌ Multi-Axis Marginalization
- **Reason:** IGoR's recursive `normalize()` handles this
- **Example excluded:** `marginalize(In, Out, {1, 2})` summing multiple dims
- **Alternative:** Compose single-dim operations

### ❌ Matrix Operations
- **Reason:** IGoR has no matrix algebra (no GEMM, no transpose)
- **Example excluded:** `matrix_product(A, B, C)`
- **Alternative:** N/A (not needed)

### ❌ Advanced Reductions
- **Reason:** IGoR doesn't use `max`, `argmax`, `mean`, etc.
- **Example excluded:** `argmax_along_axis(In, Out, dim)`
- **Alternative:** Add on-demand if use cases emerge

### ❌ Complex Probability Operations
- **Reason:** Beyond current scope
- **Example excluded:** `outer_product`, `conditional_marginalize`
- **Alternative:** Implement in Phase 2 if Rec_Event refactoring requires

---

## 9. Summary

### Operation Count: 8 Total

| Category | Operations | Status |
|----------|-----------|--------|
| Element-wise | 4 (add, multiply, divide, scale) | Essential |
| Reductions | 2 (sum, marginalize) | Essential |
| Probability | 1 (normalize) | Essential |
| Numerical Stability | 1 (log_add_exp) | Future extension |

### Design Strengths
✅ Modern C++ (`std::mdspan`, `std::dextents`)
✅ Zero-copy flexibility (`HybridBuffer`)
✅ Algorithm-data separation
✅ Minimal dependencies (standard library only)
✅ In-place operation support
✅ Cache-efficient memory layout
✅ Future-proof without overdesign

### Integration Benefits
✅ Type-safe replacement for raw pointer arithmetic
✅ Enables thread-local marginal accumulation
✅ Clear algorithmic intent
✅ Bounds checking in debug builds
✅ Smooth migration path (parallel implementation)

### Performance Goals
- **Overhead:** <5% vs raw pointers in release builds
- **Vectorization:** Compiler auto-SIMD for element-wise ops
- **Cache efficiency:** 100% utilization for contiguous access
- **Scalability:** O(n) for all operations (no hidden complexity)

---

**Status:** Ready for implementation
**Next Step:** Phase 1 - Implement `HybridBuffer` and `Tensor` classes with unit tests
