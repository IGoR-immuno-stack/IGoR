#pragma once

#include <igor/Math/Tensor.h>
#include <igor/Math/TensorCreation.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <initializer_list>
#include <cstddef>
#include <cmath>

namespace igor::math::test {

/**
 * @brief Create a tensor filled with sequential values.
 *
 * Creates a tensor where data[i] = start + i * step.
 * Useful for creating predictable test data.
 *
 * @tparam T Element type
 * @param shape Shape of the tensor
 * @param start Starting value (default: 0)
 * @param step Step between consecutive elements (default: 1)
 * @return Tensor filled with sequential values
 *
 * @code
 * auto t = test::make_sequential<double>({2, 3});  // [0,1,2,3,4,5]
 * auto t2 = test::make_sequential<double>({2, 2}, 1.0);  // [1,2,3,4]
 * @endcode
 */
template<typename T = double>
Tensor<T> make_sequential(std::initializer_list<std::size_t> shape, T start = T{0}, T step = T{1}) {
    Tensor<T> result{std::vector<std::size_t>(shape)};
    T value = start;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result.data()[i] = value;
        value += step;
    }
    return result;
}

/**
 * @brief Create a tensor with a custom fill function.
 *
 * Fills tensor using a function f(linear_index) -> value.
 *
 * @tparam T Element type
 * @tparam Func Function type
 * @param shape Shape of the tensor
 * @param fill_func Function taking linear index and returning value
 * @return Tensor filled using the function
 *
 * @code
 * auto t = test::make_tensor<double>({3, 3}, [](size_t i) { return i * i; });
 * @endcode
 */
template<typename T = double, typename Func>
Tensor<T> make_tensor(std::initializer_list<std::size_t> shape, Func&& fill_func) {
    Tensor<T> result{std::vector<std::size_t>(shape)};
    for (std::size_t i = 0; i < result.size(); ++i) {
        result.data()[i] = fill_func(i);
    }
    return result;
}

/**
 * @brief Create an identity matrix using existing tensor::eye function.
 *
 * @tparam T Element type
 * @param n Size of the square matrix
 * @return n×n identity matrix
 */
template<typename T = double>
Tensor<T> make_identity(std::size_t n) {
    return tensor::eye<T>(n);
}

/**
 * @brief Verify that all elements in a tensor equal a specific value.
 *
 * @tparam T Element type
 * @param tensor The tensor to check
 * @param expected The expected value
 * @param tolerance Tolerance for floating point comparison (default: 0 for exact match)
 *
 * @code
 * auto zeros = tensor::zeros<double>({3, 4});
 * test::require_all_equal(zeros, 0.0);
 * @endcode
 */
template<typename T>
void require_all_equal(const Tensor<T>& tensor, T expected, T tolerance = T{0}) {
    for (std::size_t i = 0; i < tensor.size(); ++i) {
        if constexpr (std::is_floating_point_v<T>) {
            if (tolerance > T{0}) {
                REQUIRE(tensor.data()[i] == Catch::Approx(expected).margin(tolerance));
            } else {
                REQUIRE(tensor.data()[i] == expected);
            }
        } else {
            REQUIRE(tensor.data()[i] == expected);
        }
    }
}

/**
 * @brief Verify that a tensor has the expected shape.
 *
 * Checks both ndim and each dimension size.
 *
 * @tparam T Element type
 * @param tensor The tensor to check
 * @param expected_shape The expected shape
 *
 * @code
 * auto t = tensor::ones<double>({3, 4, 5});
 * test::require_shape(t, {3, 4, 5});
 * @endcode
 */
template<typename T>
void require_shape(const Tensor<T>& tensor, std::initializer_list<std::size_t> expected_shape) {
    REQUIRE(tensor.ndim() == expected_shape.size());

    std::size_t dim = 0;
    for (auto expected_extent : expected_shape) {
        REQUIRE(tensor.shape()[dim] == expected_extent);
        ++dim;
    }

    // Also verify total size matches
    std::size_t expected_size = 1;
    for (auto extent : expected_shape) {
        expected_size *= extent;
    }
    REQUIRE(tensor.size() == expected_size);
}

/**
 * @brief Verify that all elements satisfy a predicate.
 *
 * @tparam T Element type
 * @tparam Func Predicate function type
 * @param tensor The tensor to check
 * @param predicate Function taking T and returning bool
 * @param message Optional message for failure (default: "")
 *
 * @code
 * auto t = tensor::ones<double>({3, 3});
 * test::require_all(t, [](double x) { return x > 0.0; });
 * @endcode
 */
template<typename T, typename Func>
void require_all(const Tensor<T>& tensor, Func&& predicate, const char* message = "") {
    for (std::size_t i = 0; i < tensor.size(); ++i) {
        INFO("Element " << i << " failed predicate" << (message[0] ? ": " : "") << message);
        REQUIRE(predicate(tensor.data()[i]));
    }
}

/**
 * @brief Check if two tensors are equal within a tolerance.
 *
 * @tparam T Element type
 * @param a First tensor
 * @param b Second tensor
 * @param tolerance Tolerance for floating point comparison
 * @return true if tensors are equal within tolerance
 */
template<typename T>
bool tensors_equal(const Tensor<T>& a, const Tensor<T>& b, T tolerance = T{1e-10}) {
    if (a.ndim() != b.ndim() || a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.ndim(); ++i) {
        if (a.shape()[i] != b.shape()[i]) {
            return false;
        }
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T>) {
            if (std::abs(a.data()[i] - b.data()[i]) > tolerance) {
                return false;
            }
        } else {
            if (a.data()[i] != b.data()[i]) {
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Require that two tensors are equal within a tolerance.
 *
 * @tparam T Element type
 * @param a First tensor
 * @param b Second tensor
 * @param tolerance Tolerance for floating point comparison
 *
 * @code
 * auto a = tensor::ones<double>({3, 3});
 * auto b = tensor::full<double>({3, 3}, 1.0);
 * test::require_tensors_equal(a, b);
 * @endcode
 */
template<typename T>
void require_tensors_equal(const Tensor<T>& a, const Tensor<T>& b, T tolerance = T{1e-10}) {
    REQUIRE(a.ndim() == b.ndim());
    REQUIRE(a.size() == b.size());

    for (std::size_t i = 0; i < a.ndim(); ++i) {
        INFO("Shape mismatch at dimension " << i);
        REQUIRE(a.shape()[i] == b.shape()[i]);
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        INFO("Value mismatch at index " << i << ": " << a.data()[i] << " vs " << b.data()[i]);
        if constexpr (std::is_floating_point_v<T>) {
            REQUIRE(a.data()[i] == Catch::Approx(b.data()[i]).margin(tolerance));
        } else {
            REQUIRE(a.data()[i] == b.data()[i]);
        }
    }
}

/**
 * @brief Fill a 2D tensor with a pattern using indices.
 *
 * Fills tensor(i,j) = func(i, j).
 *
 * @tparam T Element type
 * @tparam Func Function type
 * @param tensor The 2D tensor to fill
 * @param func Function taking (i, j) and returning T
 *
 * @code
 * Tensor<double> t({3, 3});
 * test::fill_2d(t, [](size_t i, size_t j) { return i * 10 + j; });
 * @endcode
 */
template<typename T, typename Func>
void fill_2d(Tensor<T>& tensor, Func&& func) {
    REQUIRE(tensor.ndim() == 2);
    auto view = tensor.template view<2>();
    for (std::size_t i = 0; i < tensor.shape()[0]; ++i) {
        for (std::size_t j = 0; j < tensor.shape()[1]; ++j) {
            view[i, j] = func(i, j);
        }
    }
}

/**
 * @brief Fill a 3D tensor with a pattern using indices.
 *
 * @tparam T Element type
 * @tparam Func Function type
 * @param tensor The 3D tensor to fill
 * @param func Function taking (i, j, k) and returning T
 */
template<typename T, typename Func>
void fill_3d(Tensor<T>& tensor, Func&& func) {
    REQUIRE(tensor.ndim() == 3);
    auto view = tensor.template view<3>();
    for (std::size_t i = 0; i < tensor.shape()[0]; ++i) {
        for (std::size_t j = 0; j < tensor.shape()[1]; ++j) {
            for (std::size_t k = 0; k < tensor.shape()[2]; ++k) {
                view[i, j, k] = func(i, j, k);
            }
        }
    }
}

} // namespace igor::math::test
