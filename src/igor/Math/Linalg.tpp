/*
 * Linalg.tpp
 *
 *  Created on: Feb 08, 2026
 *      Author: IGoR Agent
 */

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
    apply_binary(in1, in2, out, std::plus<>{});
}

template <typename In1, typename In2, typename Out>
void subtract(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::minus<>{});
}

template <typename In1, typename In2, typename Out>
void multiply(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::multiplies<>{});
}

template <typename In1, typename In2, typename Out>
void divide(In1 in1, In2 in2, Out out) {
    apply_binary(in1, in2, out, std::divides<>{});
}

template <typename In, typename Scalar, typename Out>
void scale(In in, Scalar s, Out out) {
    apply_unary(in, out, [s](const auto& val) { return val * s; });
}

template <typename In>
auto sum(In in) {
    using T = typename In::value_type;
    T total = 0;
    reduce(in, total, std::plus<>{});
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

} // namespace igor::math::linalg

