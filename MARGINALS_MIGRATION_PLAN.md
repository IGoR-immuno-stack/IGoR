# Model_Marginals Migration Implementation Plan

**Date:** 9 February 2026
**Status:** In Progress (Phase 0.1 Complete)
**Timeline:** 3-4 months
**Risk Level:** High (Core inference component)

---

## Executive Summary

Modernize `Model_marginals` (956 lines of complex probability management) to use the new Math module, replacing raw pointer arithmetic with type-safe Tensor operations. Migration requires 5 phases with critical prerequisites.

**Key Constraint:** Cannot start migration until Phase 0 (Math module gaps) is complete.

---

## Phase 0: Math Module Extensions

**Duration:** 2 weeks
**Risk:** Low
**Blocking:** Must complete before Phase 1
**Status:** Phase 0.1 Complete (70 assertions), Phase 0.2 In Progress

### 0.1 Axis-Specific Reductions (✅ Complete - 70 assertions)

#### Implementation Tasks

```cpp
// In src/igor/Math/Linalg.h

/**
 * \brief Sum tensor along specific axis
 * @param in Input tensor of rank N
 * @param axis Axis to sum over (0-indexed)
 * @return Tensor of rank N-1 with sums
 */
template<typename In>
auto sum_axis(In in, size_t axis) -> Tensor<typename In::value_type>;

/**
 * \brief Normalize tensor along specific axis (sum to 1 per slice)
 * @param in Input tensor
 * @param out Output tensor (same shape as in)
 * @param axis Axis to normalize over
 */
template<typename In, typename Out>
void normalize_axis(In in, Out out, size_t axis);
```

#### Example Use Case
```cpp
// Normalize P(V,D,J) over J dimension only
// Each (V,D) slice sums to 1.0
Tensor<double> probs({64, 25, 12});
linalg::normalize_axis(probs, probs, 2);  // axis=2 is J
```

#### Test Requirements
- Test all axes (0, 1, 2 for 3D tensors)
- Test edge cases (rank-1, rank-5)
- Verify numerical stability
- Test in-place vs separate output

**Deliverable:** `test_Linalg_AxisOps.cpp` with 50+ assertions

### 0.2 Tensor Slicing (Week 2)

**Status:** ✅ Complete

**Approach:** Leverage existing Kokkos mdspan dependency (already used for GCC fallback) to provide `std::experimental::submdspan`.

#### Step 1: Extend MdspanCompat.h (Hybrid Approach) - **Done**

**Strategy:** Configure `MdspanCompat.h` to prefer Kokkos mdspan/submdspan when available for slicing support.

*Note:* Native `Tensor.slice()` method is **removed** as direct `submdspan` usage on views is preferred and sufficient. The `Tensor` class remains simple (contiguous).

**Deliverable:**
- Extended `MdspanCompat.h` with `submdspan` imports.

#### Step 2: Implement Tensor::Slice Helper (Member Function) - **Done**

**Strategy:** Implement `slice<Rank>(dim, index)` as member functions in `Tensor` class.
- Unlike native tensors, these functions return `mdspan` views (specifically `layout_stride`), avoiding the complexity of creating non-contiguous `Tensor` objects.
- Handles `const` and non-`const` overloads correctly.

**Deliverable:**
- `Tensor::slice` methods in `src/igor/Math/Tensor.h` / `Tensor.tpp`
- `tst/igor/Math/test_Tensor_Slicing.cpp` with usage examples.

```cpp
// In src/igor/Math/MdspanCompat.h

#pragma once

// Try standard <mdspan> first (Clang libc++, MSVC) - C++23
#if __has_include(<mdspan>) && !defined(ADJOINTCHECK_FORCE_KOKKOS_MDSPAN)
    #include <mdspan>

    // Standard mdspan is available, but submdspan is C++26 only
    // Import submdspan from Kokkos experimental
    #if __has_include(<experimental/mdspan>)
        #include <experimental/mdspan>
        namespace std {
            // Import only submdspan utilities from Kokkos
            using ::std::experimental::submdspan;
            using ::std::experimental::full_extent;
            using ::std::experimental::full_extent_t;
        }
    #else
        #warning "submdspan not available - slicing operations disabled"
        #define IGOR_NO_SUBMDSPAN
    #endif

// Fallback: Use Kokkos mdspan entirely (for GCC on Linux)
#elif __has_include(<experimental/mdspan>)
    #include <experimental/mdspan>
    namespace std {
        // Import everything from Kokkos
        using ::std::experimental::mdspan;
        using ::std::experimental::dextents;
        using ::std::experimental::extents;
        using ::std::experimental::layout_right;
        using ::std::experimental::layout_left;
        using ::std::experimental::layout_stride;
        using ::std::experimental::default_accessor;

        // Include submdspan
        using ::std::experimental::submdspan;
        using ::std::experimental::full_extent;
        using ::std::experimental::full_extent_t;
    }

#else
    #error "No mdspan implementation found."
#endif
```

**This approach:**
- ✅ Uses standard `std::mdspan` on Clang/MSVC (better performance)
- ✅ Falls back to Kokkos `submdspan` for slicing (available now)
- ✅ Minimal changes to existing code
- ✅ Prepares for C++26 when `std::submdspan` becomes standard

#### Step 2: Implement Tensor::slice()

```cpp
// In src/igor/Math/Tensor.h

#ifndef IGOR_NO_SUBMDSPAN
/**
 * \brief Extract slice along specific dimension (reduces rank by 1)
 * @tparam Rank Input tensor rank
 * @param dim Dimension to slice
 * @param index Index in that dimension
 * @return Tensor view (borrowing) of rank N-1
 *
 * Example: tensor.slice<3>(0, 5) extracts tensor[5, :, :]
 *
 * Note: Requires Kokkos mdspan submdspan support
 */
template<std::size_t Rank>
Tensor<T> slice(size_t dim, size_t index) const {
    if (ndim() != Rank) {
        throw std::invalid_argument("Tensor rank mismatch");
    }
    if (dim >= Rank || index >= shape()[dim]) {
        throw std::invalid_argument("Slice index out of bounds");
    }

    auto view_in = this->template view<Rank>();

    // Dispatch by rank and dimension
    if constexpr (Rank == 2) {
        if (dim == 0) {
            auto sub = std::submdspan(view_in, index, std::full_extent);
            return Tensor<T>(sub.data_handle(), {sub.extent(0)},
                           HybridBuffer<T>::Borrow{});
        } else {
            auto sub = std::submdspan(view_in, std::full_extent, index);
            return Tensor<T>(sub.data_handle(), {sub.extent(0)},
                           HybridBuffer<T>::Borrow{});
        }
    }
    else if constexpr (Rank == 3) {
        std::vector<size_t> out_shape;
        if (dim == 0) {
            auto sub = std::submdspan(view_in, index, std::full_extent, std::full_extent);
            return Tensor<T>(sub.data_handle(), {sub.extent(0), sub.extent(1)},
                           HybridBuffer<T>::Borrow{});
        } else if (dim == 1) {
            auto sub = std::submdspan(view_in, std::full_extent, index, std::full_extent);
            return Tensor<T>(sub.data_handle(), {sub.extent(0), sub.extent(1)},
                           HybridBuffer<T>::Borrow{});
        } else {
            auto sub = std::submdspan(view_in, std::full_extent, std::full_extent, index);
            return Tensor<T>(sub.data_handle(), {sub.extent(0), sub.extent(1)},
                           HybridBuffer<T>::Borrow{});
        }
    }
    // Similar for Rank 4, 5...
}
#endif // IGOR_NO_SUBMDSPAN
```

**Key points:**
- Uses `HybridBuffer<T>::Borrow{}` to create zero-copy views
- Returns new Tensor that borrows data (no allocation)
- Modifying slice modifies original tensor
- Guarded by `IGOR_NO_SUBMDSPAN` for platforms without support

```cpp
TEST_CASE("Tensor Slicing with submdspan", "[Math][Tensor][Slice]") {
    SECTION("Slice 2D -> 1D") {
        Tensor<double> t({3, 4});
        // Fill: [0, 1, 2, 3]
        //       [4, 5, 6, 7]
        //       [8, 9, 10, 11]
        for (size_t i = 0; i < t.size(); ++i) t.data()[i] = i;

        auto row1 = t.slice<2>(0, 1);  // Second row
        REQUIRE(row1.ndim() == 1);
        REQUIRE(row1.shape()[0] == 4);
        REQUIRE(row1(0) == 4.0);
        REQUIRE(row1(3) == 7.0);

        auto col2 = t.slice<2>(1, 2);  // Third column
        REQUIRE(col2.ndim() == 1);
        REQUIRE(col2.shape()[0] == 3);
        REQUIRE(col2(0) == 2.0);
        REQUIRE(col2(2) == 10.0);
    }
}
```
- Extended MdspanCompat.h with hybrid mdspan/submdspan approach
- Tensor::slice() method for ranks 2-5 (guarded by IGOR_NO_SUBMDSPAN)
- test_Tensor_Slicing.cpp with 40+ assertions
- Zero-copy semantics verified (modify slice → modifies original)        // ... test slicing along each dimension

**Deliverable:**
- Extended MdspanCompat.h with submdspan imports
- Tensor::slice() method (rank 2-5)
- test_Tensor_Slicing.cpp with 40+ assertions

### 0.3 Documentation

- [x] Update Math module README with axis operations
- [ ] Add axis operations and slicing examples to TENSOR_LINALG_PROPOSAL_V2.md
- [ ] Update MATH_MODULE_TASK_PLAN.md to mark Phase 0 complete
- [ ] Document submdspan usage patterns

**Milestone:** Phase 0 Complete - Math module ready for Model_marginals integration

---

## Phase 1: Borrowing Wrapper

**Duration:** 1 week
**Risk:** Low (no changes to existing code)
**Dependencies:** Phase 0

### 1.1 Create Model_marginals_v2 Class

#### File Structure
```
src/igor/Core/
  Model_marginals.h      # Keep existing
  Model_marginals.cpp    # Keep existing
  Model_marginals_v2.h   # NEW
  Model_marginals_v2.cpp # NEW
```

#### Implementation

```cpp
// Model_marginals_v2.h
#pragma once

#include <igor/Core/Model_marginals.h>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>

namespace igor {

/**
 * \brief Modernized Model_marginals using Math module
 *
 * Phase 1: Wraps legacy array with zero-copy Tensor views
 * Phase 2+: Replaces storage with structured Tensors
 */
class Model_marginals_v2 : public Model_marginals {
public:
    Model_marginals_v2();
    Model_marginals_v2(const Model_Parms&);

    // New tensor-based interface
    math::Tensor<long double> as_tensor();
    const math::Tensor<long double> as_tensor() const;

    // New operations using Math module
    void normalize_v2(const Model_Parms&);
    void add_to_marginals_v2(double event_proba,
                            std::list<std::shared_ptr<Rec_Event>> events,
                            const Model_Parms& parms);

    // Delegate to v1 for complex operations (temporary)
    using Model_marginals::normalize;
    using Model_marginals::add_to_marginals;
};

} // namespace igor
```

#### as_tensor() Implementation

```cpp
// Model_marginals_v2.cpp

math::Tensor<long double> Model_marginals_v2::as_tensor() {
    // Borrowing mode - wraps existing array without copying
    return math::Tensor<long double>(
        marginal_array_smart_p.get(),
        {marginal_arr_size}
    );
}

void Model_marginals_v2::normalize_v2(const Model_Parms& parms) {
    auto tensor = as_tensor();
    math::linalg::normalize(tensor, tensor);  // In-place
}
```

### 1.2 Testing Strategy

#### Test File: `tst/igor/Core/test_Model_marginals_v2.cpp`

```cpp
TEST_CASE("Model_marginals_v2: Borrowing semantics", "[marginals]") {
    Model_Parms parms = load_test_model();

    Model_marginals_v2 marginals(parms);

    // Verify zero-copy: modify via tensor, check array
    auto tensor = marginals.as_tensor();
    tensor.data()[0] = 42.0;

    REQUIRE(marginals.get_array()[0] == 42.0);  // Should match
}

TEST_CASE("Model_marginals_v2: Normalization equivalence", "[marginals]") {
    Model_Parms parms = load_test_model();

    Model_marginals v1(parms);
    Model_marginals_v2 v2(parms);

    // Initialize with same data
    v1.random_initialize(parms);
    v2.copy_from(v1);

    // Normalize both ways
    v1.normalize(...);
    v2.normalize_v2(parms);

    // Should produce identical results
    for (size_t i = 0; i < parms.get_array_size(); ++i) {
        REQUIRE_THAT(v1.get_array()[i],
                     WithinRel(v2.get_array()[i], 1e-10));
    }
}
```

### 1.3 Deliverables

- [ ] `Model_marginals_v2.h` with borrowing interface
- [ ] `Model_marginals_v2.cpp` with as_tensor() implementation
- [ ] Test suite with 20+ test cases
- [ ] Documentation comparing v1 vs v2 usage

**Milestone:** Phase 1 Complete - Wrapper proven with existing array

---

## Phase 2: Replace Storage

**Duration:** 2-3 weeks
**Risk:** High (structural change)
**Dependencies:** Phase 1

### 2.1 New Storage Architecture

#### Replace Giant Array

```cpp
class Model_marginals_v2 {
private:
    // OLD (delete after migration):
    // Marginal_array_p marginal_array_smart_p;
    // size_t marginal_arr_size;

    // NEW:
    struct EventMarginal {
        math::Tensor<long double> tensor;
        size_t offset;  // Keep for backward compat during transition
        std::vector<size_t> shape;
    };

    std::unordered_map<Rec_Event_name, EventMarginal> event_marginals_;

    // Navigation helpers (keep existing logic)
    std::unordered_map<Rec_Event_name, size_t> index_map_;
    std::unordered_map<Rec_Event_name,
                       std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>
        offset_map_;
};
```

### 2.2 Constructor Refactoring

```cpp
Model_marginals_v2::Model_marginals_v2(const Model_Parms& parms) {
    // Compute shapes for each event type
    auto event_list = parms.get_event_list();

    for (auto& event_ptr : event_list) {
        Rec_Event_name name = event_ptr->get_name();

        // Compute tensor shape for this event
        auto shape = compute_event_shape(event_ptr, parms);

        // Create tensor with computed shape
        event_marginals_[name] = EventMarginal{
            math::Tensor<long double>(shape),
            compute_offset(name, parms),  // Temporary
            shape
        };
    }

    // Build navigation maps (keep existing logic)
    index_map_ = get_index_map(parms);
    offset_map_ = get_offsets_map(parms);
}
```

### 2.3 Migration Tasks

#### Week 1: Storage Infrastructure
- [ ] Design EventMarginal structure
- [ ] Implement shape computation for each event type
- [ ] Update constructor to create Tensors
- [ ] Keep offset maps for compatibility

#### Week 2: Update Core Methods
- [ ] `operator=(const Model_marginals_v2&)` - deep copy tensors
- [ ] `operator+=(Model_marginals_v2)` - use linalg::add per event
- [ ] `operator-=(Model_marginals_v2)` - use linalg::subtract per event
- [ ] `compute_size()` - sum tensor sizes
- [ ] `null_initialize()` - zero all tensors

#### Week 3: Testing & Validation
- [ ] Unit tests for each updated method
- [ ] Numerical equivalence tests vs v1
- [ ] Memory leak checks with valgrind
- [ ] Performance regression tests

### 2.4 Rollback Strategy

Keep v1 implementation intact. If Phase 2 fails:
1. Revert to Phase 1 (borrowing wrapper)
2. Identify root cause
3. Fix and retry

**Milestone:** Phase 2 Complete - Structured storage working

---

## Phase 3: Refactor Operations

**Duration:** 3-4 weeks
**Risk:** Medium (logic changes)
**Dependencies:** Phase 2

### 3.1 Simple Operations (Week 1)

#### operator+= Refactoring

```cpp
// OLD: Model_marginals.cpp:721
Model_marginals& Model_marginals::operator+=(Model_marginals marginals) {
    if (this->marginal_arr_size != marginals.marginal_arr_size) {
        throw std::runtime_error("Size mismatch");
    }
    for (size_t i = 0; i < marginal_arr_size; ++i) {
        this->marginal_array_smart_p[i] += marginals.marginal_array_smart_p[i];
    }
    return *this;
}

// NEW:
Model_marginals_v2& Model_marginals_v2::operator+=(const Model_marginals_v2& other) {
    for (auto& [name, event_data] : event_marginals_) {
        if (!other.event_marginals_.count(name)) {
            throw std::runtime_error("Event missing: " + name);
        }

        auto& this_tensor = event_data.tensor;
        const auto& other_tensor = other.event_marginals_.at(name).tensor;

        // Use Math module
        math::linalg::add(this_tensor.view(),
                         other_tensor.view(),
                         this_tensor.view());
    }
    return *this;
}
```

#### Tasks
- [ ] `operator+=` using linalg::add
- [ ] `operator-=` using linalg::subtract
- [ ] `add_pseudo_counts` using linalg::scale
- [ ] Tests for each operation

### 3.2 Complex Normalization (Week 2-4)

#### Current Normalization Logic (Complex)

```cpp
// Model_marginals.cpp:895-956 (61 lines)
// Recursive iteration through parent-child event relationships
// Custom offset calculation
// Conditional normalization based on event types
```

#### Refactoring Strategy

**Step 1:** Identify normalization patterns
```cpp
// Pattern 1: Full event normalization
for each event {
    normalize_axis(event.tensor, event.tensor, last_axis);
}

// Pattern 2: Conditional normalization (if parent realized)
for each event with parent {
    if (parent realized) {
        normalize_axis(event.tensor, event.tensor, last_axis);
    }
}

// Pattern 3: Joint normalization (parent-child)
for each parent-child pair {
    normalize_axis(parent.tensor × child.tensor, out, joint_axis);
}
```

**Step 2:** Implement helper methods
```cpp
void normalize_event(const Rec_Event_name& name, size_t axis);
void normalize_conditional(const Rec_Event_name& parent,
                          const Rec_Event_name& child,
                          size_t axis);
void normalize_joint(const std::vector<Rec_Event_name>& events);
```

**Step 3:** Refactor main normalize() method
```cpp
void Model_marginals_v2::normalize(
    std::unordered_map<Rec_Event_name, std::list<...>> inverse_map,
    std::unordered_map<Rec_Event_name, int> index_map,
    std::queue<std::shared_ptr<Rec_Event>> queue
) {
    // Traverse event hierarchy
    while (!queue.empty()) {
        auto event = queue.front();
        queue.pop();

        auto name = event->get_name();
        auto& tensor = event_marginals_[name].tensor;

        // Determine normalization axis from event type
        size_t norm_axis = get_normalization_axis(event, index_map);

        // Use Math module
        math::linalg::normalize_axis(tensor.view(),
                                     tensor.view(),
                                     norm_axis);

        // Process children
        for (auto child : event->get_children()) {
            queue.push(child);
        }
    }
}
```

#### Tasks
- [ ] Map normalization patterns to Math ops
- [ ] Implement helper methods
- [ ] Refactor main normalize() (preserve exact logic)
- [ ] Extensive numerical equivalence tests
- [ ] Edge case testing (empty events, single realization, etc.)

### 3.3 add_to_marginals Refactoring (Week 4)

```cpp
// OLD: Manual pointer arithmetic
bool Model_marginals::add_to_marginals(
    double event_proba,
    std::list<std::shared_ptr<Rec_Event>> events,
    Model_Parms parms
) {
    // Complex offset calculation
    size_t index = 0;
    for (auto& event : events) {
        index += compute_offset(event, parms);
    }
    marginal_array_smart_p[index] += event_proba;
}

// NEW: Direct tensor access
bool Model_marginals_v2::add_to_marginals(
    double event_proba,
    const std::list<std::shared_ptr<Rec_Event>>& events,
    const Model_Parms& parms
) {
    for (auto& event : events) {
        auto name = event->get_name();
        auto& tensor = event_marginals_[name].tensor;

        // Compute multi-dimensional index
        auto indices = compute_tensor_indices(event, parms);

        // Direct access
        tensor(indices) += event_proba;
    }
    return true;
}
```

**Milestone:** Phase 3 Complete - All operations use Math module

---

## Phase 4: Parallelization

**Duration:** 2-3 weeks
**Risk:** Low (optional optimization)
**Dependencies:** Phase 3

### 4.1 Thread-Local Accumulation (Week 1-2)

#### Current Bottleneck
```cpp
// Single-threaded inference
for (auto& sequence : sequences) {
    double likelihood = compute_likelihood(sequence, parms);
    marginals.add_to_marginals(likelihood, ...);  // SERIAL
}
```

#### Parallel Strategy

```cpp
void parallel_inference(
    const std::vector<Sequence>& sequences,
    const Model_Parms& parms,
    Model_marginals_v2& global_marginals
) {
    #pragma omp parallel
    {
        // Each thread has private marginals
        Model_marginals_v2 thread_marginals(parms);
        thread_marginals.null_initialize();

        // Parallel loop - no synchronization needed
        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < sequences.size(); ++i) {
            double likelihood = compute_likelihood(sequences[i], parms);
            thread_marginals.add_to_marginals(likelihood, ...);
        }

        // Reduction: merge thread-local → global
        #pragma omp critical
        {
            global_marginals += thread_marginals;  // Uses linalg::add
        }
    }
}
```

#### Benefits
- Lock-free accumulation (no contention)
- Linear scaling with cores
- Minimal code changes

#### Tasks
- [ ] Implement parallel_inference wrapper
- [ ] Benchmark with 1, 2, 4, 8, 16 threads
- [ ] Verify numerical equivalence
- [ ] Document scalability results

### 4.2 Broadcasting Optimizations (Week 2-3)

#### Use Case: Marginal Scaling

```cpp
// Scale P(V,D,J) by P(J)
void scale_by_marginal(
    math::Tensor<long double>& p_vdj,
    const math::Tensor<long double>& p_j
) {
    // OLD: Triple loop
    for (size_t v = 0; v < n_v; ++v) {
        for (size_t d = 0; d < n_d; ++d) {
            for (size_t j = 0; j < n_j; ++j) {
                p_vdj(v, d, j) *= p_j(j);
            }
        }
    }

    // NEW: Zero-copy broadcasting
    math::linalg::broadcast_multiply(
        p_vdj.view<3>(),
        p_j.view<1>(),
        p_vdj.view<3>()
    );
}
```

#### Parallel Broadcasting

```cpp
// Process V slices in parallel
#pragma omp parallel for
for (size_t v = 0; v < n_v; ++v) {
    auto p_dj = p_vdj.slice(0, v);  // Extract [D, J] slice
    math::linalg::broadcast_multiply(
        p_dj.view<2>(),
        p_j.view<1>(),
        p_dj.view<2>()
    );
}
```

**Milestone:** Phase 4 Complete - Parallel inference working

---

## Phase 5: Validation & Deployment

**Duration:** 1-2 weeks
**Risk:** Low (verification only)
**Dependencies:** Phase 4

### 5.1 Numerical Validation (Week 1)

#### Test Suite

```cpp
TEST_CASE("Numerical equivalence: v1 vs v2", "[marginals][validation]") {
    // Load production test cases
    auto test_cases = load_validation_suite();

    for (auto& test_case : test_cases) {
        Model_marginals v1(test_case.parms);
        Model_marginals_v2 v2(test_case.parms);

        // Run same inference
        v1.run_inference(test_case.sequences);
        v2.run_inference(test_case.sequences);

        // Compare all marginals
        for (auto& event_name : test_case.events) {
            auto probs_v1 = v1.get_event_probs(event_name);
            auto probs_v2 = v2.get_event_probs(event_name);

            for (size_t i = 0; i < probs_v1.size(); ++i) {
                REQUIRE_THAT(probs_v1[i],
                            WithinRel(probs_v2[i], 1e-9));
            }
        }
    }
}
```

#### Test Cases
- [ ] Human TRB (100 sequences)
- [ ] Mouse TRB (100 sequences)
- [ ] Human BCR heavy chain (100 sequences)
- [ ] Edge cases (single sequence, all invalid, etc.)

### 5.2 Performance Benchmarking (Week 1-2)

#### Benchmark Suite

```cpp
BENCHMARK("Inference: v1 single-thread") {
    Model_marginals v1(parms);
    return v1.run_inference(benchmark_sequences);
};

BENCHMARK("Inference: v2 single-thread") {
    Model_marginals_v2 v2(parms);
    return v2.run_inference(benchmark_sequences);
};

BENCHMARK("Inference: v2 parallel (8 threads)") {
    Model_marginals_v2 v2(parms);
    return v2.parallel_inference(benchmark_sequences, 8);
};
```

#### Metrics to Measure
- Throughput (sequences/second)
- Memory usage (RSS, peak)
- Cache performance (L1/L2 hit rates)
- Scaling efficiency (speedup vs threads)

#### Performance Targets
- Single-thread: within 5% of v1
- 8-thread: 6-7x speedup
- Memory: ≤ 1.5x of v1

### 5.3 Documentation

- [ ] Migration guide for users
- [ ] API documentation (Doxygen)
- [ ] Performance comparison report
- [ ] Known issues and limitations

### 5.4 Deployment Plan

#### Gradual Rollout

**Stage 1: Internal Testing** (1 week)
- Deploy to test environment
- Run nightly validation suite
- Monitor for crashes/errors

**Stage 2: Beta Testing** (2 weeks)
- Selected users opt-in to v2
- Collect feedback
- Fix reported issues

**Stage 3: Production** (Phased)
- Week 1: 10% of production traffic
- Week 2: 50% of production traffic
- Week 3: 100% migration
- Week 4: Remove v1 code

**Rollback Triggers:**
- Numerical divergence > 1e-6
- Performance regression > 10%
- Crash rate increase

**Milestone:** Phase 5 Complete - v2 in production

---

## Risk Management

### High-Risk Items

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| Phase 0 takes longer | Blocks all work | Medium | Start immediately, parallel dev |
| Storage migration breaks logic | Critical failure | Medium | Extensive testing, keep v1 |
| Normalization equivalence fails | Wrong results | High | Numerical tests, gradual rollout |
| Performance regression | User complaints | Low | Benchmarking, optimization |
| Parallel bugs (race conditions) | Data corruption | Medium | Thread sanitizer, stress tests |

### Mitigation Strategies

1. **Keep v1 Intact**
   - No modifications to Model_marginals.cpp during migration
   - v2 is additive, not replacement initially

2. **Incremental Testing**
   - Test after each phase
   - Numerical validation at every step

3. **Parallel Development**
   - Phase 0 (Math) can proceed independently
   - Phase 1 wrapper can be prototyped early

4. **Performance Monitoring**
   - Continuous benchmarking
   - Early detection of regressions

---

## Success Criteria

### Phase 0
- [x] 70 tests passing for axis operations (exceeds 50+ target)
- [x] 40+ tests passing for slicing (using Kokkos submdspan)
- [x] Documentation complete

### Phase 1
- [ ] 20+ tests for borrowing wrapper
- [ ] Zero-copy verified
- [ ] v1/v2 equivalence on simple operations

### Phase 2
- [ ] All constructors working
- [ ] operator+= and operator-= pass tests
- [ ] Memory usage reasonable (≤2x of v1 temporarily)

### Phase 3
- [ ] All 20+ Model_marginals methods refactored
- [ ] Numerical equivalence on 100+ test cases
- [ ] No performance regression

### Phase 4
- [ ] 6-7x speedup on 8 cores
- [ ] Thread-safety verified (no data races)
- [ ] Broadcasting optimizations measured

### Phase 5
- [ ] All validation tests passing
- [ ] Performance targets met
- [ ] Production deployment successful

---

## Resources Required

### Personnel
- 1 Senior Engineer (lead implementation)
- 1 Engineer (testing & validation)
- 1 Part-time (code review & consultation)

### Timeline
- **Week 1:** Phase 0.1 (Axis operations) ✅ **COMPLETE**
- **Week 2:** Phase 0.2 (Slicing with Kokkos submdspan) ✅ **COMPLETE**
- **Week 3:** Phase 1 (Wrapper)
- **Weeks 4-6:** Phase 2 (Storage)
- **Weeks 7-10:** Phase 3 (Operations)
- **Weeks 11-12:** Phase 4 (Parallelization)
- **Weeks 13-14:** Phase 5 (Validation)

**Total:** 14 weeks (3.5 months)

### Infrastructure
- Dedicated test server for benchmarking
- CI/CD pipeline for automated testing
- Profiling tools (perf, valgrind, thread sanitizer)

---

## Appendix A: Code Examples

### Current vs Modernized Comparison

#### Example 1: Normalization

**Current (v1):**
```cpp
// 61 lines of complex recursive logic
void Model_marginals::normalize(...) {
    // Manual pointer arithmetic
    // Nested loops with offset calculation
    // Parent-child relationship handling
}
```

**Modernized (v2):**
```cpp
void Model_marginals_v2::normalize(...) {
    for (auto& [name, event_data] : event_marginals_) {
        size_t axis = get_normalization_axis(name);
        linalg::normalize_axis(event_data.tensor,
                              event_data.tensor,
                              axis);
    }
}
```

#### Example 2: Accumulation

**Current (v1):**
```cpp
// Manual index calculation
size_t index = base_offset;
for (auto& event : events) {
    index += event_offset * realization_index;
}
marginal_array[index] += probability;
```

**Modernized (v2):**
```cpp
// Direct tensor access
auto& tensor = event_marginals_[event_name].tensor;
tensor(v_index, d_index, j_index) += probability;
```

---

## Appendix B: Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-02-09 | Phase 0 mandatory before migration | Axis operations critical for normalization |
| 2026-02-09 | Keep v1 intact during migration | Enables rollback, A/B testing |
| 2026-02-09 | Use std::unordered_map for events | Clean separation, easier debugging |
| 2026-02-09 | Parallel by default in Phase 4 | OpenMP widely available, low effort |

---

*Last Updated: 9 February 2026*
*Status: Planning Complete - Awaiting Approval*
