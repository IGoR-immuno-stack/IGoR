// TensorCreation.tpp ---

#pragma once

#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <string>

namespace igor::math::tensor {

template<typename T>
Tensor<T> zeros(std::initializer_list<std::size_t> shape) {
    Tensor<T> result{std::vector<std::size_t>(shape)};
    std::fill(result.begin(), result.end(), T(0));
    return result;
}

template<typename T>
Tensor<T> ones(std::initializer_list<std::size_t> shape) {
    Tensor<T> result{std::vector<std::size_t>(shape)};
    std::fill(result.begin(), result.end(), T(1));
    return result;
}

template<typename T>
Tensor<T> full(std::initializer_list<std::size_t> shape, T value) {
    Tensor<T> result{std::vector<std::size_t>(shape)};
    std::fill(result.begin(), result.end(), value);
    return result;
}

template<typename T>
Tensor<T> arange(T start, T stop, T step) {
    if (step == T(0)) {
        throw std::invalid_argument("arange: step cannot be zero");
    }

    if ((step > 0 && start >= stop) || (step < 0 && start <= stop)) {
        return Tensor<T>({0});
    }

    std::size_t n = static_cast<std::size_t>(std::ceil((stop - start) / step));
    if (n == 0) return Tensor<T>({0});

    Tensor<T> result({n});
    T value = start;
    for (std::size_t i = 0; i < n; ++i) {
        result.data()[i] = value;
        value += step;
    }
    return result;
}

template<typename T>
Tensor<T> arange(T stop) {
    return arange(T(0), stop, T(1));
}

template<typename T>
Tensor<T> linspace(T start, T stop, std::size_t num, bool endpoint) {
    if (num == 0) return Tensor<T>({0});

    if (num == 1) {
        Tensor<T> result({1});
        result.data()[0] = start;
        return result;
    }

    Tensor<T> result({num});
    T step = (stop - start) / (endpoint ? (num - 1) : num);

    for (std::size_t i = 0; i < num; ++i) {
        result.data()[i] = start + i * step;
    }
    return result;
}

template<typename T>
Tensor<T> eye(std::size_t n, std::size_t m) {
    if (m == 0) m = n;

    Tensor<T> result({n, m});
    std::fill(result.begin(), result.end(), T(0));

    std::size_t diag_len = std::min(n, m);
    auto view = result.template view<2>();
    for (std::size_t i = 0; i < diag_len; ++i) {
        view[i, i] = T(1);
    }
    return result;
}

template<typename T>
Tensor<T> identity(std::size_t n) {
    return eye<T>(n);
}

template<typename T>
Tensor<T> from_array(std::initializer_list<std::size_t> shape,
                     std::initializer_list<T> data) {
    std::size_t expected_size = std::accumulate(
        shape.begin(), shape.end(),
        std::size_t(1), std::multiplies<std::size_t>()
    );

    if (data.size() != expected_size) {
        throw std::invalid_argument(
            "from_array: data size (" + std::to_string(data.size()) +
            ") doesn't match shape size (" + std::to_string(expected_size) + ")"
        );
    }

    Tensor<T> result{std::vector<std::size_t>(shape)};
    std::copy(data.begin(), data.end(), result.begin());
    return result;
}

#ifndef IGOR_NO_SUBMDSPAN
template<std::size_t Rank, typename T>
auto slice(Tensor<T>& tensor, size_t dim, size_t index) {
    auto v = tensor.template view<Rank>();

    // Define the common return type (strided view of Rank-1)
    using return_type = std::mdspan<T, std::dextents<size_t, Rank-1>, std::layout_stride>;

    // Rank 2
    if constexpr (Rank == 2) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, index));
    }
    // Rank 3
    else if constexpr (Rank == 3) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent, std::full_extent));
        else if (dim == 1) return return_type(std::submdspan(v, std::full_extent, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, std::full_extent, index));
    }
    // Rank 4
    else if constexpr (Rank == 4) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent, std::full_extent, std::full_extent));
        else if (dim == 1) return return_type(std::submdspan(v, std::full_extent, index, std::full_extent, std::full_extent));
        else if (dim == 2) return return_type(std::submdspan(v, std::full_extent, std::full_extent, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, std::full_extent, std::full_extent, index));
    }
    else {
        static_assert(Rank <= 4, "Slicing implemented up to Rank 4");
        return return_type(std::submdspan(v, index)); // dummy
    }
}

template<std::size_t Rank, typename T>
auto slice(const Tensor<T>& tensor, size_t dim, size_t index) {
    auto v = tensor.template view<Rank>();

    // Define the common return type (strided view of Rank-1)
    using return_type = std::mdspan<const T, std::dextents<size_t, Rank-1>, std::layout_stride>;

    // Rank 2
    if constexpr (Rank == 2) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, index));
    }
    // Rank 3
    else if constexpr (Rank == 3) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent, std::full_extent));
        else if (dim == 1) return return_type(std::submdspan(v, std::full_extent, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, std::full_extent, index));
    }
    // Rank 4
    else if constexpr (Rank == 4) {
        if (dim == 0) return return_type(std::submdspan(v, index, std::full_extent, std::full_extent, std::full_extent));
        else if (dim == 1) return return_type(std::submdspan(v, std::full_extent, index, std::full_extent, std::full_extent));
        else if (dim == 2) return return_type(std::submdspan(v, std::full_extent, std::full_extent, index, std::full_extent));
        return return_type(std::submdspan(v, std::full_extent, std::full_extent, std::full_extent, index));
    }
    else {
        static_assert(Rank <= 4, "Slicing implemented up to Rank 4");
        return return_type(std::submdspan(v, index)); // dummy
    }
}
#endif

template<typename T>
Tensor<T> empty(std::initializer_list<std::size_t> shape) {
    // Just allocate without initialization
    return Tensor<T>{std::vector<std::size_t>(shape)};
}

template<typename T>
Tensor<T> zeros_like(const Tensor<T>& other) {
    Tensor<T> result{std::vector<std::size_t>(other.shape().begin(), other.shape().end())};
    std::fill(result.begin(), result.end(), T{0});
    return result;
}

template<typename T>
Tensor<T> ones_like(const Tensor<T>& other) {
    Tensor<T> result{std::vector<std::size_t>(other.shape().begin(), other.shape().end())};
    std::fill(result.begin(), result.end(), T{1});
    return result;
}

template<typename T>
Tensor<T> full_like(const Tensor<T>& other, T value) {
    Tensor<T> result{std::vector<std::size_t>(other.shape().begin(), other.shape().end())};
    std::fill(result.begin(), result.end(), value);
    return result;
}

template<typename T>
Tensor<T> diag(std::initializer_list<T> diagonal) {
    std::size_t n = diagonal.size();
    Tensor<T> result = zeros<T>({n, n});

    auto view = result.template view<2>();
    std::size_t idx = 0;
    for (const auto& val : diagonal) {
        view[idx, idx] = val;
        ++idx;
    }

    return result;
}

template<typename T>
Tensor<T> diag(const Tensor<T>& matrix) {
    if (matrix.ndim() != 2) {
        throw std::invalid_argument("diag: input must be a 2D tensor");
    }

    std::size_t n = std::min(matrix.shape()[0], matrix.shape()[1]);
    Tensor<T> result{std::vector<std::size_t>{n}};

    auto mat_view = matrix.template view<2>();
    auto diag_view = result.template view<1>();

    for (std::size_t i = 0; i < n; ++i) {
        diag_view[i] = mat_view[i, i];
    }

    return result;
}

} // namespace igor::math::tensor

// TensorCreation.tpp ends here
