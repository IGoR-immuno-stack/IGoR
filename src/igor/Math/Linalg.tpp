/*
 * Linalg.tpp
 *
 *  Created on: Feb 08, 2026
 *      Author: IGoR Agent
 */

#pragma once

#include <igor/Math/Tensor.h>

#include <algorithm>
// #include <numeric>
// #include <concepts>
#include <array>
#include <cmath>
#include <functional>

namespace igor::math::linalg {

namespace detail {

    // Recursive N-dim Iterator for Binary Ops
    template <size_t Rank, size_t CurrentDim>
    struct NDimIterator {
        template <typename In1, typename In2, typename Out, typename Func, typename... Indices>
        static void apply(In1 in1, In2 in2, Out out, Func f, Indices... idx) {
            for (size_t i = 0; i < out.extent(CurrentDim); ++i) {
                if constexpr (CurrentDim == Rank - 1) {
                     out[idx..., i] = f(in1[idx..., i], in2[idx..., i]);
                } else {
                     NDimIterator<Rank, CurrentDim + 1>::apply(in1, in2, out, f, idx..., i);
                }
            }
        }
    };

    // Recursive N-dim Iterator for Unary Ops
    template <size_t Rank, size_t CurrentDim>
    struct NDimIteratorUnary {
        template <typename In, typename Out, typename Func, typename... Indices>
        static void apply(In in, Out out, Func f, Indices... idx) {
            for (size_t i = 0; i < out.extent(CurrentDim); ++i) {
                if constexpr (CurrentDim == Rank - 1) {
                     out[idx..., i] = f(in[idx..., i]);
                } else {
                     NDimIteratorUnary<Rank, CurrentDim + 1>::apply(in, out, f, idx..., i);
                }
            }
        }
    };

    // Recursive N-dim Reducer
    template <size_t Rank, size_t CurrentDim>
    struct NDimReducer {
        template <typename In, typename T, typename Func, typename... Indices>
        static void apply(In in, T& acc, Func f, Indices... idx) {
            for (size_t i = 0; i < in.extent(CurrentDim); ++i) {
                if constexpr (CurrentDim == Rank - 1) {
                     acc = f(acc, in[idx..., i]);
                } else {
                     NDimReducer<Rank, CurrentDim + 1>::apply(in, acc, f, idx..., i);
                }
            }
        }
    };

    // Recursive N-dim ArgMax Reducer
    template <size_t Rank, size_t CurrentDim>
    struct ArgMaxReducer {
        template <typename In, typename T, typename... Indices>
        static void apply(In in, T& max_val, std::array<size_t, Rank>& max_idx, Indices... idx) {
            for (size_t i = 0; i < in.extent(CurrentDim); ++i) {
                if constexpr (CurrentDim == Rank - 1) {
                     T val = in[idx..., i];
                     // Use strict greater to get first occurrence (stable)
                     if (val > max_val) {
                         max_val = val;
                         max_idx = { static_cast<size_t>(idx)..., i };
                     }
                } else {
                     ArgMaxReducer<Rank, CurrentDim + 1>::apply(in, max_val, max_idx, idx..., i);
                }
            }
        }
    };

    // Helper to insert index at axis
    template <size_t InRank, typename... Indices>
    auto insert_axis_index(size_t axis, size_t val, Indices... idx) {
        std::array<size_t, InRank> result;
        std::array<size_t, InRank - 1> current = { static_cast<size_t>(idx)... };
        
        size_t src = 0;
        for (size_t dst = 0; dst < InRank; ++dst) {
            if (dst == axis) {
                result[dst] = val;
            } else {
                result[dst] = current[src++];
            }
        }
        return result;
    }

    // Generic Axis Operation Iterator
    // Iterates over R-1 dimensions (skipping axis)
    // Func: (indices...) -> void
    template <size_t OutRank, size_t CurrentDim>
    struct AxisIterator {
        template <typename Shape, typename Func, typename... Indices>
        static void apply(const Shape& shape, Func f, Indices... idx) {
            for (size_t i = 0; i < shape[CurrentDim]; ++i) {
                if constexpr (CurrentDim == OutRank - 1) {
                    f(idx..., i);
                } else {
                    AxisIterator<OutRank, CurrentDim + 1>::apply(shape, f, idx..., i);
                }
            }
        }
    };

} // namespace detail

template <typename In1, typename In2, typename Out, typename Func>
void apply_binary(In1 in1, In2 in2, Out out, Func f) {
    constexpr size_t R = In1::rank();
    static_assert(In2::rank() == R, "Ranks must match");
    static_assert(Out::rank() == R, "Ranks must match");

    detail::NDimIterator<R, 0>::apply(in1, in2, out, f);
}

template <typename In, typename Out, typename Func>
void apply_unary(In in, Out out, Func f) {
    constexpr size_t R = In::rank();
    static_assert(Out::rank() == R, "Ranks must match");

    detail::NDimIteratorUnary<R, 0>::apply(in, out, f);
}

template <typename In, typename T, typename Func>
void reduce(In in, T& acc, Func f) {
    constexpr size_t R = In::rank();
    detail::NDimReducer<R, 0>::apply(in, acc, f);
}

template <typename In1, typename In2, typename Out>
void add(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::plus<void>{});
}

template <typename In1, typename In2, typename Out>
void subtract(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::minus<void>{});
}

template <typename In1, typename In2, typename Out>
void multiply(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::multiplies<void>{});
}

template <typename In1, typename In2, typename Out>
void divide(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::divides<void>{});
}

template <typename In, typename Scalar, typename Out>
void scale(In in, Scalar s, Out out) {
    apply_unary(in, out, [s](const auto& val) { return val * s; });
}

template <typename In>
auto sum(In in) {
    using T = typename In::value_type;
    T total = 0;
    reduce(in, total, std::plus<void>{});
    return total;
}

template <typename In>
auto max(In in) {
    using T = typename In::value_type;
    T val = std::numeric_limits<T>::lowest();
    reduce(in, val, [](const T& a, const T& b) { return std::max(a, b); });
    return val;
}

template <typename In>
auto min(In in) {
    using T = typename In::value_type;
    T val = std::numeric_limits<T>::max();
    reduce(in, val, [](const T& a, const T& b) { return std::min(a, b); });
    return val;
}

template <typename In>
auto argmax(In in) {
    using T = typename In::value_type;
    constexpr size_t R = In::rank();

    T max_val = std::numeric_limits<T>::lowest();
    std::array<size_t, R> max_idx{};

    detail::ArgMaxReducer<R, 0>::apply(in, max_val, max_idx);

    return max_idx;
}

// -----------------------------------------------------------------------------
// Tensor Dispatch Helpers
// -----------------------------------------------------------------------------

namespace detail {

    template <typename T, typename Func>
    void dispatch_binary(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out, Func f) {
        if (in1.ndim() != in2.ndim() || in1.ndim() != out.ndim()) {
            throw std::invalid_argument("Tensor rank mismatch");
        }

#ifdef IGOR_MATH_DEBUG
        // Validate shape equality
        for (size_t i = 0; i < in1.ndim(); ++i) {
            if (in1.shape()[i] != in2.shape()[i] || in1.shape()[i] != out.shape()[i]) {
                throw std::invalid_argument("Tensor shape mismatch at dimension " + std::to_string(i) +
                                          ": in1[" + std::to_string(in1.shape()[i]) +
                                          "] vs in2[" + std::to_string(in2.shape()[i]) +
                                          "] vs out[" + std::to_string(out.shape()[i]) + "]");
            }
        }
#endif

        switch (in1.ndim()) {
            case 1: f(in1.template view<1>(), in2.template view<1>(), out.template view<1>()); break;
            case 2: f(in1.template view<2>(), in2.template view<2>(), out.template view<2>()); break;
            case 3: f(in1.template view<3>(), in2.template view<3>(), out.template view<3>()); break;
            case 4: f(in1.template view<4>(), in2.template view<4>(), out.template view<4>()); break;
            case 5: f(in1.template view<5>(), in2.template view<5>(), out.template view<5>()); break;
            default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
        }
    }

    template <typename T, typename Scalar, typename Func>
    void dispatch_unary_scalar(const Tensor<T>& in, Scalar s, Tensor<T>& out, Func f) {
        if (in.ndim() != out.ndim()) {
            throw std::invalid_argument("Tensor rank mismatch");
        }

#ifdef IGOR_MATH_DEBUG
        // Validate shape equality
        for (size_t i = 0; i < in.ndim(); ++i) {
            if (in.shape()[i] != out.shape()[i]) {
                throw std::invalid_argument("Tensor shape mismatch at dimension " + std::to_string(i) +
                                          ": in[" + std::to_string(in.shape()[i]) +
                                          "] vs out[" + std::to_string(out.shape()[i]) + "]");
            }
        }
#endif

        switch (in.ndim()) {
            case 1: f(in.template view<1>(), s, out.template view<1>()); break;
            case 2: f(in.template view<2>(), s, out.template view<2>()); break;
            case 3: f(in.template view<3>(), s, out.template view<3>()); break;
            case 4: f(in.template view<4>(), s, out.template view<4>()); break;
            case 5: f(in.template view<5>(), s, out.template view<5>()); break;
            default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
        }
    }

    template <typename T, typename Func>
    auto dispatch_reduce(const Tensor<T>& in, Func f) {
        switch (in.ndim()) {
            case 1: return f(in.template view<1>());
            case 2: return f(in.template view<2>());
            case 3: return f(in.template view<3>());
            case 4: return f(in.template view<4>());
            case 5: return f(in.template view<5>());
            default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
        }
    }
}

// -----------------------------------------------------------------------------
// Tensor Overload Implementations
// -----------------------------------------------------------------------------

template <typename T>
void add(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out) {
    detail::dispatch_binary(in1, in2, out, [](auto i1, auto i2, auto o) {
        add(i1, i2, o);
    });
}

template <typename T>
void subtract(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out) {
    detail::dispatch_binary(in1, in2, out, [](auto i1, auto i2, auto o) {
        subtract(i1, i2, o);
    });
}

template <typename T>
void multiply(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out) {
    detail::dispatch_binary(in1, in2, out, [](auto i1, auto i2, auto o) {
        multiply(i1, i2, o);
    });
}

template <typename T>
void divide(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out) {
    detail::dispatch_binary(in1, in2, out, [](auto i1, auto i2, auto o) {
        divide(i1, i2, o);
    });
}

template <typename T, typename Scalar>
void scale(const Tensor<T>& in, Scalar s, Tensor<T>& out) {
    detail::dispatch_unary_scalar(in, s, out, [](auto i, auto sc, auto o) {
        scale(i, sc, o);
    });
}

template <typename T>
T sum(const Tensor<T>& in) {
    return detail::dispatch_reduce(in, [](auto v) { return sum(v); });
}

template <typename T>
T max(const Tensor<T>& in) {
    return detail::dispatch_reduce(in, [](auto v) { return max(v); });
}

template <typename T>
T min(const Tensor<T>& in) {
    return detail::dispatch_reduce(in, [](auto v) { return min(v); });
}

template <typename T>
auto argmax(const Tensor<T>& in) {
    switch (in.ndim()) {
        case 1: { auto idx = argmax(in.template view<1>()); return std::vector<size_t>(idx.begin(), idx.end()); }
        case 2: { auto idx = argmax(in.template view<2>()); return std::vector<size_t>(idx.begin(), idx.end()); }
        case 3: { auto idx = argmax(in.template view<3>()); return std::vector<size_t>(idx.begin(), idx.end()); }
        case 4: { auto idx = argmax(in.template view<4>()); return std::vector<size_t>(idx.begin(), idx.end()); }
        case 5: { auto idx = argmax(in.template view<5>()); return std::vector<size_t>(idx.begin(), idx.end()); }
        default: throw std::runtime_error("Unsupported rank");
    }
}

// -----------------------------------------------------------------------------
// Axis-Specific Operations - Implementation
// -----------------------------------------------------------------------------

template <typename T>
Tensor<T> sum_axis(const Tensor<T>& in, size_t axis) {
    if (axis >= in.ndim()) {
        throw std::invalid_argument("Axis out of bounds");
    }
    
    // Calculate output shape
    std::vector<size_t> out_shape;
    out_shape.reserve(in.ndim() - 1);
    for (size_t i = 0; i < in.ndim(); ++i) {
        if (i != axis) {
            out_shape.push_back(in.shape()[i]);
        }
    }
    
    Tensor<T> out(out_shape);

    // Dispatch based on Input Rank
    auto dispatch_sum = [&]<size_t Rank>(auto in_view) {
        if constexpr (Rank == 1) {
             // Rank 1 -> Scalar (represented as 0-dim tensor or just computed)
             // But map_view<0> is not supported by standard variant switch usually
             // Special case: summing the only axis.
             // We'll treat out as 0-dim if possible, but Tensor<0> support is limited
             // For now assume Rank > 1 as per plan constraints or hack around it
             // Actually, sum(vector) -> scalar. 
             // We'll trust In::rank constraint
        }
        
        // Output Rank is Rank - 1
        constexpr size_t OutRank = Rank - 1;
        auto out_view = out.template view<OutRank>(); // View into output tensor
        
        // Shape for generic iterator (output shape)
        std::array<size_t, OutRank> shape_arr;
        for(size_t i=0; i<OutRank; ++i) shape_arr[i] = out_view.extent(i);

        // Iterate over output dimensions
        detail::AxisIterator<OutRank, 0>::apply(shape_arr, [&](auto... idx) {
            // For each output position, sum along axis
            T acc = 0;
            size_t axis_dim = in.shape()[axis];
            
            for (size_t k = 0; k < axis_dim; ++k) {
                // Construct input index
                auto full_idx = detail::insert_axis_index<Rank>(axis, k, idx...);
                // Apply input index
                // We need to apply std::array as arguments to mdspan
                // C++23 mdspan supports spans/arrays? 
                // Or we unpack. Unpacking array to variadic is annoying.
                // Alternative: access raw pointer with strides? No.
                // We'll usage a helper to apply array
                // But full_idx is std::array.
                // Let's rely on Tensor's templated operator which accepts indices... 
                // Wait, mdspan::operator[] takes indices.
                // We need `std::apply`.
                acc += std::apply([&](auto... args) { return in_view[args...]; }, full_idx);
            }
            out_view[idx...] = acc;
        });
    };

    switch (in.ndim()) {
        case 2: dispatch_sum.template operator()<2>(in.template view<2>()); break;
        case 3: dispatch_sum.template operator()<3>(in.template view<3>()); break;
        case 4: dispatch_sum.template operator()<4>(in.template view<4>()); break;
        case 5: dispatch_sum.template operator()<5>(in.template view<5>()); break;
        default: throw std::runtime_error("Unsupported rank for sum_axis (requires 2-5)");
    }

    return out;
}

template <typename T>
void normalize_axis(const Tensor<T>& in, Tensor<T>& out, size_t axis) {
    if (in.ndim() != out.ndim()) throw std::invalid_argument("Rank mismatch");
    if (axis >= in.ndim()) throw std::invalid_argument("Axis out of bounds");

    // Check shapes match
    for (size_t i=0; i<in.ndim(); ++i) {
        if (in.shape()[i] != out.shape()[i]) throw std::invalid_argument("Shape mismatch");
    }

    auto dispatch_norm = [&]<size_t Rank>(auto in_view, auto out_view) {
        // Output shape for iterator (same as input) but we want to iterate over "base" dims
        // We can reuse the logic: iterate over all dims EXCEPT axis
        constexpr size_t BaseRank = Rank - 1;
        
        std::array<size_t, BaseRank> base_shape;
        size_t j = 0;
        for (size_t i=0; i<Rank; ++i) {
            if (i != axis) base_shape[j++] = in_view.extent(i);
        }

        // Iterate over base dimensions
        detail::AxisIterator<BaseRank, 0>::apply(base_shape, [&](auto... idx) {
            // 1. Compute sum
            T acc = 0;
            size_t axis_dim = in_view.extent(axis);
            
            for (size_t k = 0; k < axis_dim; ++k) {
                auto full_idx = detail::insert_axis_index<Rank>(axis, k, idx...);
                acc += std::apply([&](auto... args) { return in_view[args...]; }, full_idx);
            }
            
            // 2. Normalize
            if (acc != 0) {
                T scale = T(1) / acc;
                for (size_t k = 0; k < axis_dim; ++k) {
                    auto full_idx = detail::insert_axis_index<Rank>(axis, k, idx...);
                    // We need to apply to OUT as well
                    T val = std::apply([&](auto... args) { return in_view[args...]; }, full_idx);
                    std::apply([&](auto... args) -> decltype(auto) { return out_view[args...]; }, full_idx) = val * scale;
                }
            } else {
                 // Zero sum: leave as 0? or NaN? 
                 // Usually 0 in probabilistic models
                 for (size_t k = 0; k < axis_dim; ++k) {
                     auto full_idx = detail::insert_axis_index<Rank>(axis, k, idx...);
                     std::apply([&](auto... args) -> decltype(auto) { return out_view[args...]; }, full_idx) = 0;
                 }
            }
        });
    };

    switch (in.ndim()) {
        case 1: {
            // Rank 1 special case: explicit loop easier? Or reuse logic with BaseRank=0 ??
            // If BaseRank=0, AxisIterator should run once.
            // Let's see if generic works.
            // BaseRank=0. Iterator::apply(empty_shape, f).
            // loop i=0 to shape[0] (which is invalid access).
            // Need special handling for Rank 1?
            // "AxisIterator" assumes OutRank > 0.
            auto v_in = in.template view<1>();
            auto v_out = out.template view<1>();
            T acc = 0;
            for(size_t i=0; i<v_in.extent(0); ++i) acc += v_in[i];
            if(acc != 0) scale(v_in, T(1)/acc, v_out);
            else scale(v_in, 0, v_out);
            break;
        }
        case 2: dispatch_norm.template operator()<2>(in.template view<2>(), out.template view<2>()); break;
        case 3: dispatch_norm.template operator()<3>(in.template view<3>(), out.template view<3>()); break;
        case 4: dispatch_norm.template operator()<4>(in.template view<4>(), out.template view<4>()); break;
        case 5: dispatch_norm.template operator()<5>(in.template view<5>(), out.template view<5>()); break;
        default: throw std::runtime_error("Unsupported rank");
    }
}

// -----------------------------------------------------------------------------
// Probabilistic Operations (Phase 2b) - Core Implementations
// -----------------------------------------------------------------------------

template <typename In1, typename In2>
auto dot(In1 in1, In2 in2) {
    static_assert(In1::rank() == 1 && In2::rank() == 1, "Dot product requires rank 1 tensors");
    if (in1.extent(0) != in2.extent(0)) {
        throw std::invalid_argument("Dot product dimension mismatch");
    }

    using T = typename In1::value_type;
    T sum = 0;
    for (size_t i = 0; i < in1.extent(0); ++i) {
        sum += in1[i] * in2[i];
    }
    return sum;
}

template <typename In1, typename In2, typename Out>
void matmul(In1 A, In2 B, Out C) {
    static_assert(In1::rank() == 2 && In2::rank() == 2 && Out::rank() == 2,
                  "Matmul requires rank 2 tensors");

    auto M = A.extent(0);
    auto K = A.extent(1);
    auto N = B.extent(1);

    if (B.extent(0) != K || C.extent(0) != M || C.extent(1) != N) {
        throw std::invalid_argument("Matmul dimension mismatch (Expected [M,K] x [K,N] -> [M,N])");
    }

    // Naive Triple Loop
    for (size_t m = 0; m < M; ++m) {
        for (size_t n = 0; n < N; ++n) {
            typename Out::value_type acc = 0;
            for (size_t k = 0; k < K; ++k) {
                acc += A[m, k] * B[k, n];
            }
            C[m, n] = acc;
        }
    }
}

template <typename In, typename Out>
void normalize(In in, Out out) {
    auto total = sum(in);
    if (total == typename In::value_type(0)) {
        throw std::invalid_argument("Cannot normalize: tensor sum is zero");
    }
    scale(in, typename In::value_type(1) / total, out);
}

// -----------------------------------------------------------------------------
// Probabilistic Operations - Tensor Overloads
// -----------------------------------------------------------------------------

template <typename T>
T dot(const Tensor<T>& in1, const Tensor<T>& in2) {
    if (in1.ndim() != 1 || in2.ndim() != 1) {
        throw std::invalid_argument("Dot product requires 1D tensors");
    }
    return dot(in1.template view<1>(), in2.template view<1>());
}

template <typename T>
void matmul(const Tensor<T>& A, const Tensor<T>& B, Tensor<T>& C) {
    if (A.ndim() != 2 || B.ndim() != 2 || C.ndim() != 2) {
        throw std::invalid_argument("Matmul requires 2D tensors");
    }
    matmul(A.template view<2>(), B.template view<2>(), C.template view<2>());
}

template <typename T>
void normalize(const Tensor<T>& in, Tensor<T>& out) {
    if (in.ndim() != out.ndim()) {
        throw std::invalid_argument("Tensor rank mismatch");
    }

    switch (in.ndim()) {
        case 1: normalize(in.template view<1>(), out.template view<1>()); break;
        case 2: normalize(in.template view<2>(), out.template view<2>()); break;
        case 3: normalize(in.template view<3>(), out.template view<3>()); break;
        case 4: normalize(in.template view<4>(), out.template view<4>()); break;
        case 5: normalize(in.template view<5>(), out.template view<5>()); break;
        default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
    }
}

// -----------------------------------------------------------------------------
// Log-Space & Stability (Phase 2c) - Core Implementations
// -----------------------------------------------------------------------------

template <typename T>
T log_add_exp(T a, T b) {
    if (a == -std::numeric_limits<T>::infinity()) return b;
    if (b == -std::numeric_limits<T>::infinity()) return a;
    if (a > b) return a + std::log1p(std::exp(b - a));
    return b + std::log1p(std::exp(a - b));
}

template <typename In, typename Out>
void center(In in, Out out) {
    auto m = max(in);
    apply_unary(in, out, [m](auto x) { return x - m; });
}

template <typename In>
auto log_sum_exp(In in) {
    using T = typename In::value_type;
    T m = max(in);
    if (m == -std::numeric_limits<T>::infinity()) return m;

    T sum_exp = 0;
    reduce(in, sum_exp, [m](T acc, T val) {
        return acc + std::exp(val - m);
    });
    return m + std::log(sum_exp);
}

template <typename In, typename Out>
void log_normalize(In in, Out out) {
    auto lse = log_sum_exp(in);
    apply_unary(in, out, [lse](auto x) { return x - lse; });
}


// -----------------------------------------------------------------------------
// Log-Space & Stability - Tensor Overloads
// -----------------------------------------------------------------------------

template <typename T>
void center(const Tensor<T>& in, Tensor<T>& out) {
    if (in.ndim() != out.ndim()) {
        throw std::invalid_argument("Tensor rank mismatch");
    }

    switch (in.ndim()) {
        case 1: center(in.template view<1>(), out.template view<1>()); break;
        case 2: center(in.template view<2>(), out.template view<2>()); break;
        case 3: center(in.template view<3>(), out.template view<3>()); break;
        case 4: center(in.template view<4>(), out.template view<4>()); break;
        case 5: center(in.template view<5>(), out.template view<5>()); break;
        default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
    }
}

template <typename T>
T log_sum_exp(const Tensor<T>& in) {
    return detail::dispatch_reduce(in, [](auto v) { return log_sum_exp(v); });
}

template <typename T>
void log_normalize(const Tensor<T>& in, Tensor<T>& out) {
    if (in.ndim() != out.ndim()) {
        throw std::invalid_argument("Tensor rank mismatch");
    }

    switch (in.ndim()) {
        case 1: log_normalize(in.template view<1>(), out.template view<1>()); break;
        case 2: log_normalize(in.template view<2>(), out.template view<2>()); break;
        case 3: log_normalize(in.template view<3>(), out.template view<3>()); break;
        case 4: log_normalize(in.template view<4>(), out.template view<4>()); break;
        case 5: log_normalize(in.template view<5>(), out.template view<5>()); break;
        default: throw std::runtime_error("Unsupported rank (supported: 1-5)");
    }
}

// -----------------------------------------------------------------------------
// Broadcasting Support (Phase 2d) - Implementations
// -----------------------------------------------------------------------------

namespace detail {
    template <std::size_t Rank, std::size_t... Is>
    auto make_dextents_from_array(const std::array<std::size_t, Rank>& arr, std::index_sequence<Is...>) {
        return std::dextents<std::size_t, Rank>(arr[Is]...);
    }
}

template<typename InMdspan, std::size_t OutRank>
auto broadcast_to(
    InMdspan in,
    std::array<std::size_t, OutRank> target_shape
) -> std::mdspan<typename InMdspan::element_type, std::dextents<std::size_t, OutRank>, std::layout_stride>
{
    using T = typename InMdspan::element_type;
    constexpr std::size_t InRank = InMdspan::extents_type::rank();

    static_assert(OutRank >= InRank, "Target rank must be >= Input rank");

    std::array<std::size_t, OutRank> new_strides;
    std::fill(new_strides.begin(), new_strides.end(), 0);

    int in_idx = static_cast<int>(InRank) - 1;
    int out_idx = static_cast<int>(OutRank) - 1;

    while (out_idx >= 0) {
        if (in_idx < 0) {
            new_strides[out_idx] = 0; // New dimension
        }
        else {
            size_t in_dim = in.extent(in_idx);
            size_t out_dim = target_shape[out_idx];

            if (in_dim == out_dim) {
                new_strides[out_idx] = in.stride(in_idx);
                in_idx--;
            }
            else if (in_dim == 1) {
                new_strides[out_idx] = 0; // Broadcast singleton
                in_idx--;
            }
            else {
                throw std::invalid_argument("Incompatible broadcast shapes");
            }
        }
        out_idx--;
    }

    using Extents = std::dextents<std::size_t, OutRank>;
    using Mapping = typename std::layout_stride::template mapping<Extents>;

    auto extents = detail::make_dextents_from_array(target_shape, std::make_index_sequence<OutRank>{});

    // Layout stride mapping constructor: mapping(extents, strides)
    Mapping map(extents, new_strides);

    return std::mdspan<T, Extents, std::layout_stride>(in.data_handle(), map);
}

template<typename LeftMdspan, typename RightMdspan, typename OutMdspan>
void broadcast_multiply(
    LeftMdspan left,
    RightMdspan right,
    OutMdspan out
) {
    constexpr std::size_t LeftRank = LeftMdspan::extents_type::rank();

    std::array<std::size_t, LeftRank> target_shape;
    for (size_t i = 0; i < LeftRank; ++i) {
        target_shape[i] = left.extent(i);
    }

    // Create broadcast view
    auto right_broadcasted = broadcast_to(right, target_shape);

    // multiply works generically on any mdspan with operator[]
    multiply(left, right_broadcasted, out);
}

} // namespace igor::math::linalg

