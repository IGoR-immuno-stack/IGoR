# Tensor Layout Support: Complexity Analysis

**Question:** How much complexity does adding row-major vs column-major support add?

---

## Current Implementation (Row-Major Only)

**Lines of Code:** ~300 lines
**Complexity:** Simple - single stride computation, single mdspan layout

```cpp
// Current: Hard-coded row-major
static shape_type compute_strides(const shape_type& e) {
    // Row-major: strides = [cols*rows*..., cols, 1]
    for (std::size_t i = e.size(); i > 0; --i) {
        s[i - 1] = running_stride;
        running_stride *= e[i - 1];
    }
}
```

---

## Option 1: Template Parameter (Compile-Time Layout)

**Added Complexity:** ~50 lines, moderate complexity increase

### Changes Required:

```cpp
// 1. Add layout enum
enum class Layout { RowMajor, ColMajor };

// 2. Template on layout
template <typename T, Layout L = Layout::RowMajor>
class Tensor {
    // ...

    // 3. Conditional stride computation
    static shape_type compute_strides(const shape_type& e) {
        if constexpr (L == Layout::RowMajor) {
            // Current implementation
        } else {
            // Column-major: strides = [1, rows, rows*cols, ...]
            size_type running_stride = 1;
            for (std::size_t i = 0; i < e.size(); ++i) {
                s[i] = running_stride;
                running_stride *= e[i];
            }
        }
    }

    // 4. Conditional mdspan layout
    template<std::size_t Rank, std::size_t... Is>
    auto make_view_impl(std::index_sequence<Is...>) {
        using LayoutPolicy = std::conditional_t<
            L == Layout::RowMajor,
            std::layout_right,
            std::layout_left
        >;
        return std::mdspan<T, std::dextents<size_type, Rank>, LayoutPolicy>(
            storage_.data(), dims_[Is]...
        );
    }
};
```

**Pros:**
- ✅ Zero runtime overhead (compile-time dispatch)
- ✅ Type safety - layout is part of type
- ✅ Clean - no variant explosion

**Cons:**
- ❌ Different types: `Tensor<double, RowMajor>` ≠ `Tensor<double, ColMajor>`
- ❌ Cannot change layout at runtime
- ❌ Template bloat (two instantiations per value type)

**Complexity Score:** 6/10

---

## Option 2: Runtime Layout (Current Design + Runtime Flag)

**Added Complexity:** ~80 lines, higher complexity

### Changes Required:

```cpp
template <typename T>
class Tensor {
    // 1. Add layout storage
    Layout layout_ = Layout::RowMajor;

    // 2. Runtime stride computation
    static shape_type compute_strides(const shape_type& e, Layout layout) {
        if (layout == Layout::RowMajor) {
            // Row-major computation
        } else {
            // Column-major computation
        }
    }

    // 3. DOUBLE the variant types (one for each layout)
    using view_variant = std::variant<
        // Row-major variants
        std::mdspan<T, std::dextents<size_type, 1>, std::layout_right>,
        std::mdspan<T, std::dextents<size_type, 2>, std::layout_right>,
        // ... up to 5

        // Column-major variants
        std::mdspan<T, std::dextents<size_type, 1>, std::layout_left>,
        std::mdspan<T, std::dextents<size_type, 2>, std::layout_left>,
        // ... up to 5
    >;
    // Now 10 variant alternatives instead of 5!

    // 4. Dispatch in variant() factory
    view_variant variant() {
        if (layout_ == Layout::RowMajor) {
            switch (dims_.size()) {
                case 1: return view_row_major<1>();
                case 2: return view_row_major<2>();
                // ...
            }
        } else {
            switch (dims_.size()) {
                case 1: return view_col_major<1>();
                case 2: return view_col_major<2>();
                // ...
            }
        }
    }
};
```

**Pros:**
- ✅ Single type: `Tensor<double>` for both layouts
- ✅ Can change layout at runtime (via reshape/transpose)
- ✅ Flexible for model loading (detect layout from file)

**Cons:**
- ❌ Variant size doubles (10 alternatives instead of 5)
- ❌ More runtime checks
- ❌ Slightly slower compilation (larger variant)

**Complexity Score:** 8/10

---

## Option 3: Hybrid (Template + mdspan Auto-Detect)

**Added Complexity:** ~40 lines, minimal complexity

### Key Insight: Let mdspan handle it!

```cpp
template <typename T>
class Tensor {
    // Keep simple: always store row-major strides
    // But let mdspan views handle layout conversion

    template<std::size_t Rank, typename LayoutPolicy = std::layout_right>
    auto view() {
        // mdspan can reinterpret layout without changing storage
        return std::mdspan<T, std::dextents<size_type, Rank>, LayoutPolicy>(
            storage_.data(), dims_[Is]...
        );
    }

    // Usage:
    auto row_view = tensor.view<2, std::layout_right>();  // Row-major
    auto col_view = tensor.view<2, std::layout_left>();   // Col-major
};
```

**Pros:**
- ✅ Minimal code change
- ✅ Zero storage overhead
- ✅ mdspan handles layout conversion
- ✅ User chooses at call site

**Cons:**
- ❌ Confusing: storage layout ≠ view layout
- ❌ Potential performance trap (strided access)
- ❌ Not suitable if you need to **store** column-major data

**Complexity Score:** 4/10

---

## Recommendation

### For IGoR's Use Case:

**If models are always row-major:**
- ✅ **Keep current implementation** (0 added complexity)
- Add Option 3 (view-level layout) only if needed for interop

**If models can be either layout:**
- ✅ **Use Option 1** (template parameter) - best performance
- Instantiate both: `using RowTensor = Tensor<double, RowMajor>;`

**If layout determined at runtime (from file header):**
- ✅ **Use Option 2** (runtime layout) - most flexible
- Accept ~80 lines of added complexity

---

## Code Changes Estimate

| Option | Lines Added | Files Modified | Test Cases Added | Build Time Impact |
|--------|-------------|----------------|------------------|-------------------|
| Option 1 (Template) | ~50 | 1 (Tensor.h) | +3 sections | +10% (template instantiation) |
| Option 2 (Runtime) | ~80 | 2 (Tensor.h, tests) | +5 sections | +15% (larger variant) |
| Option 3 (Hybrid) | ~40 | 1 (Tensor.h) | +2 sections | +5% (minimal) |

---

## Concrete Implementation Preview

### Option 1 Implementation (Recommended):

```cpp
// In Tensor.h
enum class Layout { RowMajor, ColMajor };

template <typename T, Layout L = Layout::RowMajor>
class Tensor {
    // Stride computation becomes:
    static shape_type compute_strides(const shape_type& e) {
        shape_type s(e.size());
        if (e.empty()) return s;

        if constexpr (L == Layout::RowMajor) {
            // Current: strides decrease [N*M*..., M, 1]
            size_type running_stride = 1;
            for (std::size_t i = e.size(); i > 0; --i) {
                s[i - 1] = running_stride;
                running_stride *= e[i - 1];
            }
        } else {
            // Column-major: strides increase [1, N, N*M, ...]
            size_type running_stride = 1;
            for (std::size_t i = 0; i < e.size(); ++i) {
                s[i] = running_stride;
                running_stride *= e[i];
            }
        }
        return s;
    }

    // View creation becomes:
    template<std::size_t Rank, std::size_t... Is>
    auto make_view_impl(std::index_sequence<Is...>) {
        using LayoutPolicy = std::conditional_t<
            L == Layout::RowMajor,
            std::layout_right,
            std::layout_left
        >;
        return std::mdspan<T, std::dextents<size_type, Rank>, LayoutPolicy>(
            storage_.data(), dims_[Is]...
        );
    }
};

// Usage:
Tensor<double, Layout::RowMajor> row_tensor({10, 20});
Tensor<double, Layout::ColMajor> col_tensor({10, 20});
```

**Changes Summary:**
1. Add `Layout` enum (1 line)
2. Add template parameter (1 line)
3. Update `compute_strides` (8 lines)
4. Update `make_view_impl` (5 lines)
5. Update `make_view_const_impl` (5 lines)

**Total:** ~20 lines changed, ~50 lines if you count tests and type aliases.

---

## Question for You:

**What's your use case?**

1. **Models are always row-major** → Keep current (0 complexity)
2. **Need both, known at compile-time** → Option 1 (low complexity, best performance)
3. **Layout determined by file header at runtime** → Option 2 (moderate complexity)
4. **Just need different views on same data** → Option 3 (minimal complexity)

Let me know and I can implement the chosen option!
