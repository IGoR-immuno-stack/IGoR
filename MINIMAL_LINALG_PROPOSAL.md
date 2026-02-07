# Minimal Linear Algebra Module Proposal

## Design Philosophy
*   **Modern C++**: Utilizing C++23 `std::mdspan` for multi-dimensional views.
*   **Zero-Copy Flexibility**: Hybrid storage model (Owning vs. Borrowing) to solve the `std::vector` limitation.
*   **Algorithm-Data Separation**: Implementing algorithms similar to the C++26 `std::linalg` proposal (P1673), acting on views rather than containers.
*   **No Heavy Dependencies**: Pure C++ standard library usage.

## 1. Storage: The `Tensor` Class
An N-dimensional array container that manages memory lifecycle but delegating indexing to `mdspan`.

### Core Components
1.  **`HybridBuffer<T>`**:
    *   Simple wrapper managing a pointer.
    *   Flag `is_owning`: If true, `delete[]` in destructor. If false, do nothing.
    *   Supports construction from `size` (allocation) or `T*` (wrapping).

2.  **`Tensor<T, Rank>`**:
    *   Holds a `HybridBuffer<T>`.
    *   **Holds `std::extents<size_t, Rank>` (shape)**: Critically defines the dimensions of the data.
    *   Method `view()`: Returns `std::mdspan<double, std::dextents<size_t, Rank>>`.

```cpp
// Sketch
template<typename T, size_t Rank>
class Tensor {
    HybridBuffer<T> storage;
    std::dextents<size_t, Rank> dims; // Restored dimensions

public:
    // Calling with extents -> Allocates
    Tensor(std::dextents<size_t, Rank> e) : storage(product(e)), dims(e) {}

    // Calling with ptr + extents -> Wraps
    Tensor(T* ptr, std::dextents<size_t, Rank> e) : storage(ptr, product(e)), dims(e) {}

    auto view() { return std::mdspan(storage.data(), dims); }
};
```

## 2. Algorithms: `igor::linalg`
A set of templated functions operating strictly on `std::mdspan`. This ensures they work on whole Tensors, sub-slices, or external buffers equally well.

### Required Operations for IGoR
1.  **`element_wise_add(In1, In2, Out)`**:
    *   Equivalent to `Out = In1 + In2`.
    *   Optimized path for contiguous layouts (1D pointer loop).
    *   Used for: Accumulating thread-local marginals.

2.  **`scale(In, Scalar, Out)`**:
    *   Equivalent to `Out = In * Scalar`.
    *   Used for: Normalization.

3.  **`sum(In)`**:
    *   Full reduction (sum of all elements).
    *   Used for: Computing normalization factors.

4.  **`marginalize(In, Out, dim)`**:
    *   Tensor contraction: Sums values over a specific axis.
    *   Example: `(64, 4)` -> `sum(dim=1)` -> `(64)`.
    *   Used for: Bayesian Inversion (`invert_edge`), reducing joint probabilities.

5.  **`element_wise_multiply(A, B, C)`**:
    *   Hadamard product ($C_i = A_i \times B_i$).
    *   Used for: Joint probabilities, probability reweighting.

6.  **`element_wise_divide(A, B, C)`**:
    *   Element-wise division ($C_i = A_i / B_i$).
    *   Used for: Bayesian inversion (Model_marginals.cpp line 340).

7.  **`normalize(A, B, dim)`**:
    *   Probability constraint: $\sum_i P(i) = 1$ along dimension.
    *   Used for: Model_marginals.cpp line 943, fundamental requirement.

8.  **`log_add_exp(log_a, log_b)`** (future extension):
    *   Numerically stable: $\log(e^a + e^b)$.
    *   Future-proofing for extreme probabilities ($P < 10^{-100}$).

**Total: 8 operations** (minimal + future-proof, not overdesigned)

**Note:** Matrix operations excluded - IGoR has no matrix algebra needs.

## 3. Integration Plan

### Model_marginals Refactoring
Instead of:
```cpp
// Current
unique_ptr<double[]> giant_array;
double* get_event_block(name); // unsafe
```
We move to:
```cpp
// New
std::map<std::string, Tensor<double, Dynamic>> event_tables;
auto view = event_tables[name].view(); // safe mdspan
linalg::scale(view, 0.0); // clean API
```

### Safety & Performance
*   **Bounds Checking**: `mdspan` supports distinct implementations (debug vs release), allowing safety during dev and raw pointer speed during production.
*   **Vectorization**: By exposing contiguous pointers via `mdspan` accessors when layout permits, the compiler can auto-vectorize the `linalg` loops.

---

## Linear Algebra Quality Analysis

**Date:** 6 February 2026
**Review Focus:** Mathematical completeness and correctness for Bayesian inference

### ✅ Strengths

**1. Clean Abstraction Philosophy**
- Algorithm-data separation is mathematically sound
- Operations on views (not containers) matches modern numerical computing practice
- Matches `std::linalg` design principles (P1673)

**2. Core Operations are Well-Chosen**
The 5 proposed operations cover basic tensor arithmetic:
- `element_wise_add` - vector space addition
- `scale` - scalar multiplication
- `sum` - reduction operation
- `marginalize` - tensor contraction
- `matrix_product` - GEMM

---

### 🔴 Critical Mathematical Gaps

#### **Gap 1: Missing Tensor Contraction Generality**

Your `marginalize(In, Out, dim)` is too limited:

```cpp
// Your proposal: sum over single dimension
marginalize(In, Out, dim)  // (64, 4) --[dim=1]--> (64)

// Bayesian inference needs: sum over MULTIPLE dimensions
// P(A) = Σ_B Σ_C P(A,B,C)
marginalize(In, Out, {dim1, dim2})  // (64, 4, 3) --[{1,2}]--> (64)

// Also needs: CONDITIONAL marginalization
// P(A|B=b₀) = Σ_C P(A,B=b₀,C) / Σ_A Σ_C P(A,B=b₀,C)
conditional_marginalize(In, Out, condition_dim, condition_value)
```

**Why Critical:** Joint probability distributions have 3-7 dimensions in IGoR. Must sum over arbitrary subsets.

---

#### **Gap 2: No Numerical Stability Operations**

**Problem:** Bayesian inference multiplies many small probabilities:
```
P(sequence) = P(V) × P(D|V) × P(J|V,D) × ...
            = 0.01 × 0.05 × 0.02 × ...
            = 10⁻⁵⁰  ⟹ underflow!
```

**Missing Operations:**
```cpp
// Log-space arithmetic to prevent underflow
// Instead of: P = P1 × P2 × P3
// Compute: log(P) = log(P1) + log(P2) + log(P3)

log_add_exp(log_a, log_b)        // log(eᵃ + eᵇ) - numerically stable
log_space_multiply(log_a, log_b) // log(eᵃ × eᵇ) = a + b
log_normalize(In, Out, dims)     // exp(log_a - log_sum_exp(log_a))
```

**Mathematical Identity:**
```
log(eᵃ + eᵇ) = max(a,b) + log(1 + exp(-|a-b|))
```

This is **mandatory** for stable probability computation, not optional.

---

#### **Gap 3: Missing Probability-Specific Operations**

**Outer Product (Tensor Product):**
```cpp
// Build joint distribution from independent marginals
// P(A,B) = P(A) ⊗ P(B)
outer_product(A, B, Out)  // (n) × (m) --> (n,m)
```

**Element-wise Operations:**
```cpp
// Bayesian update: P(A|B) ∝ P(B|A) × P(A)
element_wise_multiply(likelihood, prior, posterior)

// Hadamard product for probability reweighting
hadamard_product(weights, probabilities, reweighted)
```

#### **2. Normalize Implementation**

**API:**
```cpp
void normalize(In input, Out result, size_t dim);
```

**Behavior:** Ensures probability constraint $\sum_i P(i) = 1$ along dimension `dim`.

**Example:**
```cpp
// Marginals shape: (64 genes, 4 deletions)
normalize(marginals, out, 1);  // Each gene's deletions sum to 1

// Implementation:
// 1. Sum along dimension: s = sum(input, dim)
// 2. Divide: result = input / s
```

**In-place safe:** `normalize(data, data, dim)` is valid.



#### **4. Memory Layout**

**Default:** Row-major (C-style, `std::layout_right`)

```cpp
// Tensor shape (64, 4): 64 genes × 4 deletions
// Memory: [gene0_del0, gene0_del1, ..., gene0_del3, gene1_del0, ...]
//         A[i,j] = data[i * 4 + j]
```

**Rationale:** Matches C++ convention, optimal for IGoR's access patterns.

**Extensibility:** Can template on `std::layout_left` later if needed:
```cpp
template<typename T, size_t Rank, typename Layout = std::layout_right>
class Tensor { /* ... */ };
```

**Performance note:** Marginalization over trailing dimension is cache-friendly.
Sum over dimension d:
- If d is contiguous: sequential memory access ✓
- If d is strided: cache-unfriendly ✗

**Required:** Explicitly specify default layout and provide stride-aware algorithms.

---

#### **Gap 6: No Reduction Operation Completeness**

Your `sum()` is one reduction. Need full suite:

```cpp
// Statistical reductions for posterior analysis
max(In)              // Maximum probability value
argmax(In)           // Index of maximum (most likely realization)
min(In)              // Minimum
mean(In)             // Average (for expectation values)
variance(In)         // For uncertainty quantification

// Conditional reductions
max_along_axis(In, Out, dim)     // Max over specified dimension
argmax_along_axis(In, Out, dim)  // Indices of maxima
```

**Use Case:**
```cpp
// Find most likely gene choice: argmax(P(gene))
auto most_likely_gene = argmax(gene_marginals);

// Compute expected deletion length: Σᵢ i × P(deletion=i)
auto expected_length = weighted_sum(deletion_marginals, indices);
```

---

#### **Gap 7: Linear Algebra Composability**

**Problem:** Operations don't compose mathematically.

**Example Workflow:**
```cpp
// Compute: P(A) = Σ_B,C P(A,B,C) / Σ_A,B,C P(A,B,C)
auto joint_sum = sum(joint);              // Scalar
auto marginal = marginalize(joint, {1,2}); // Vector
element_wise_divide(marginal, joint_sum, normalized); // Missing!
```

**Missing:**
```cpp
element_wise_divide(A, B, Out)  // A[i] / B[i]
reciprocal(In, Out)             // 1 / In[i]
clamp(In, Out, min, max)        // Bound values (avoid 0 or infinity)
```

---

### 🎯 Linear Algebra Correctness Issues

#### **Issue 1: `matrix_product` Ambiguity**

You write:
> `matrix_product(A, B, C)` - Classic GEMM (C = A × B)

**Mathematical Questions:**
1. What about α, β coefficients? (GEMM is actually `C = αAB + βC`)
2. Transpose options? (`A^T × B`, `A × B^T`, `A^T × B^T`)
3. Strided/non-contiguous input handling?

**Standard BLAS signature:**
```cpp
gemm(TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
```

**Minimal useful signature:**
```cpp
matrix_multiply(A, B, C, transpose_A=false, transpose_B=false)
```

---

#### **Issue 2: `marginalize` Dimension Semantics**

```cpp
marginalize(In, Out, dim)
```

**Ambiguity:**
- Is `dim` the dimension to **eliminate** (sum over)?
- Or the dimension to **keep**?

**NumPy convention:** `sum(A, axis=1)` eliminates dimension 1.

**Specify explicitly!**

---

#### **Issue 3: In-Place Operations Unspecified**

Many operations allow `In == Out` (in-place):
```cpp
scale(A, 2.0, A)  // A *= 2  (in-place)
```

**Critical for memory efficiency!**

But requires:
1. Mathematically correct ordering (no read-after-write hazards)
2. Explicit documentation of which operations support it

---

### 📊 Completeness Score for IGoR's Needs

| Category | Proposed | Required | Status |
|----------|----------|----------|--------|
| Basic Arithmetic | ✓ | ✓ | ✅ Complete |
| Tensor Contraction | Partial | Multi-dim | ⚠️ Insufficient |
| Log-Space Ops | ✗ | ✓ | ❌ Missing |
| Probability Ops | ✗ | ✓ | ❌ Missing |
| Broadcasting | ✗ | ✓ | ❌ Missing |
| Reductions | Partial | Full suite | ⚠️ Insufficient |
| Matrix Ops | Minimal | Transpose variants | ⚠️ Insufficient |

**Current Coverage: ~40% of required operations**

---

### 🔬 Required Mathematical Operations (Complete List)

#### **Tier 1: Critical (Cannot work without these)**
```cpp
// Multi-dimensional marginalization
marginalize(In, Out, dims_to_sum)

// Numerical stability
log_add_exp(log_a, log_b)
log_normalize(In, Out, dims)

// Probability normalization
normalize(In, Out, dims)

// Broadcasting multiply
broadcast_multiply(A, B, Out)

// Conditional slicing
conditional_marginalize(In, Out, condition_dim, condition_val)
```

#### **Tier 2: Essential (Needed for efficiency)**
```cpp
// Outer product
outer_product(A, B, Out)

// Element-wise operations
element_wise_multiply(A, B, Out)
element_wise_divide(A, B, Out)

// Statistical reductions
max_along_axis(In, Out, dim)
argmax_along_axis(In, Out, dim)
mean_along_axis(In, Out, dim)

// In-place variants
add_in_place(target, source)
scale_in_place(A, scalar)
```

#### **Tier 3: Optimization (Nice to have)**
```cpp
// Fused operations (fewer memory passes)
fused_multiply_add(A, B, C, Out)  // Out = A * B + C

// Vectorized reductions
weighted_sum(values, weights, Out)

// Matrix operations with options
gemm(A, B, C, trans_A, trans_B, alpha, beta)
```

---

### 📝 Revised Algorithm Priority

**Replace the list of 5 with this list of 8:**

1. **`marginalize(In, Out, dims)`** - Multi-dimensional tensor contraction ⭐
2. **`normalize(In, Out, dims)`** - Probability constraint enforcement ⭐
3. **`log_normalize(In, Out, dims)`** - Numerically stable normalization ⭐
4. **`broadcast_multiply(A, B, Out)`** - Vectorized probability scaling ⭐
5. **`outer_product(A, B, Out)`** - Joint distribution construction
6. **`element_wise_multiply(A, B, Out)`** - Bayesian update
7. **`conditional_marginalize(In, Out, cond_dim, cond_val)`** - Conditional probability
8. **`argmax_along_axis(In, Out, dim)`** - Most likely scenario finding

(⭐ = absolutely mandatory, cannot work without these)

---

### ✅ Quality Assessment Summary

**Mathematical Rigor: 6/10**
- Core concepts sound
- Missing essential probability operations
- No numerical stability considerations
- Broadcasting not addressed
- Layout/stride semantics unspecified

**Completeness: 4/10**
- ~40% of required operations
- Critical gaps in multi-dimensional operations
- Missing log-space arithmetic (mandatory for the domain)

**Correctness: 7/10**
- What's there is mathematically correct
- But insufficient specification (transpose options, dimension semantics, in-place support)

**Overall Verdict:** Solid foundation, needs ~8 additional operations and clarification of 3 ambiguities to be viable for Bayesian inference.

---

### 📋 Recommendations for Revision

1. **Add multi-dimensional marginalization** - `marginalize(In, Out, std::vector<size_t> dims)`
2. **Add log-space operations** - Critical for numerical stability
3. **Add normalize operation** - Probability constraint enforcement
4. **Specify broadcasting semantics** - How shapes combine
5. **Specify memory layout** - Row-major vs column-major default
6. **Add element-wise multiply/divide** - Bayesian update operations
7. **Clarify dimension semantics** - Which dims are eliminated vs kept
8. **Document in-place support** - Which operations allow it

With these additions, the proposal becomes mathematically complete for Bayesian probabilistic inference.
