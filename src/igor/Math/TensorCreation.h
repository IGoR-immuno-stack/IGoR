// TensorCreation.h ---

#pragma once

/**
 * @file TensorCreation.h
 * @brief Factory functions for creating Tensor objects.
 *
 * Provides numpy/xtensor-style creation functions:
 * - zeros, ones, full: Fill with constant values
 * - arange, linspace: Generate sequences
 * - eye, identity: Create identity matrices
 * - from_array: Create from initializer lists
 */

#include <igor/Math/Tensor.h>
#include <initializer_list>
#include <cstddef>

namespace igor::math::tensor {

/**
 * @brief Create tensor filled with zeros.
 * @tparam T Element type (default: double)
 * @param shape Dimensions of the tensor
 * @return Tensor filled with zeros
 *
 * Example:
 * @code
 * auto t = tensor::zeros({3, 4});  // 3x4 tensor of zeros
 * @endcode
 */
template<typename T = double>
Tensor<T> zeros(std::initializer_list<std::size_t> shape);

/**
 * @brief Create tensor filled with ones.
 * @tparam T Element type (default: double)
 * @param shape Dimensions of the tensor
 * @return Tensor filled with ones
 *
 * Example:
 * @code
 * auto t = tensor::ones({2, 2});  // 2x2 tensor of ones
 * @endcode
 */
template<typename T = double>
Tensor<T> ones(std::initializer_list<std::size_t> shape);

/**
 * @brief Create tensor filled with a constant value.
 * @tparam T Element type (default: double)
 * @param shape Dimensions of the tensor
 * @param value Value to fill the tensor with
 * @return Tensor filled with the specified value
 *
 * Example:
 * @code
 * auto t = tensor::full({3, 3}, 5.0);  // 3x3 tensor filled with 5.0
 * @endcode
 */
template<typename T = double>
Tensor<T> full(std::initializer_list<std::size_t> shape, T value);

/**
 * @brief Create 1D tensor with evenly spaced values.
 * @tparam T Element type (default: double)
 * @param start Start value (inclusive)
 * @param stop End value (exclusive)
 * @param step Step size
 * @return 1D tensor with values [start, start+step, ..., stop)
 *
 * Example:
 * @code
 * auto t = tensor::arange(0.0, 10.0, 2.0);  // [0, 2, 4, 6, 8]
 * @endcode
 */
template<typename T = double>
Tensor<T> arange(T start, T stop, T step = T(1));

/**
 * @brief Create 1D tensor with values from 0 to stop (exclusive).
 * @tparam T Element type (default: double)
 * @param stop End value (exclusive)
 * @return 1D tensor with values [0, 1, ..., stop)
 *
 * Example:
 * @code
 * auto t = tensor::arange(5.0);  // [0, 1, 2, 3, 4]
 * @endcode
 */
template<typename T = double>
Tensor<T> arange(T stop);

/**
 * @brief Create 1D tensor with linearly spaced values.
 * @tparam T Element type (default: double)
 * @param start Start value
 * @param stop End value
 * @param num Number of values to generate
 * @param endpoint If true, include stop value; if false, exclude it
 * @return 1D tensor with linearly spaced values
 *
 * Example:
 * @code
 * auto t = tensor::linspace(0.0, 1.0, 5);  // [0.0, 0.25, 0.5, 0.75, 1.0]
 * @endcode
 */
template<typename T = double>
Tensor<T> linspace(T start, T stop, std::size_t num, bool endpoint = true);

/**
 * @brief Create 2D identity matrix.
 * @tparam T Element type (default: double)
 * @param n Number of rows
 * @param m Number of columns (default: 0 means square matrix, m=n)
 * @return 2D tensor with ones on the diagonal and zeros elsewhere
 *
 * Example:
 * @code
 * auto I = tensor::eye(3);      // 3x3 identity
 * auto E = tensor::eye(3, 5);   // 3x5 matrix with diagonal ones
 * @endcode
 */
template<typename T = double>
Tensor<T> eye(std::size_t n, std::size_t m = 0);

/**
 * @brief Create square identity matrix.
 * @tparam T Element type (default: double)
 * @param n Size of the square matrix
 * @return Square identity matrix of size n×n
 *
 * Example:
 * @code
 * auto I = tensor::identity(4);  // 4x4 identity matrix
 * @endcode
 */
template<typename T = double>
Tensor<T> identity(std::size_t n);

/**
 * @brief Create tensor from initializer list data.
 * @tparam T Element type
 * @param shape Dimensions of the tensor
 * @param data Initializer list of data values (row-major order)
 * @return Tensor initialized with the provided data
 * @throws std::invalid_argument if data size doesn't match shape
 *
 * Example:
 * @code
 * auto t = tensor::from_array<double>(
 *     {2, 3},
 *     {1.0, 2.0, 3.0,
 *      4.0, 5.0, 6.0}
 * );
 * @endcode
 */
template<typename T>
Tensor<T> from_array(std::initializer_list<std::size_t> shape,
                     std::initializer_list<T> data);

#ifndef IGOR_NO_SUBMDSPAN
/**
 * @brief Get a slice of the tensor along a dimension.
 *
 * Returns a strided mdspan view of rank Rank-1 representing the slice.
 * This is a zero-copy operation that creates a view into the original tensor.
 *
 * @tparam Rank The rank of the input tensor (2, 3, or 4)
 * @tparam T The value type of the tensor
 * @param tensor The input tensor to slice
 * @param dim The dimension along which to slice (0-indexed)
 * @param index The index along the dimension to extract
 * @return A strided mdspan view of rank Rank-1
 *
 * @code
 * auto tensor = tn::ones<double>({3, 4, 5});
 * auto slice = tn::slice<3>(tensor, 0, 1);  // Get slice at index 1 along first dimension
 * @endcode
 */
template<std::size_t Rank, typename T>
auto slice(Tensor<T>& tensor, size_t dim, size_t index);

template<std::size_t Rank, typename T>
auto slice(const Tensor<T>& tensor, size_t dim, size_t index);
#endif

/**
 * @brief Create an uninitialized tensor with the given shape.
 *
 * Unlike zeros() or ones(), this function does not initialize the tensor data,
 * which can be more efficient when the values will be immediately overwritten.
 *
 * @tparam T The value type (default: double)
 * @param shape The shape of the tensor
 * @return A new tensor with uninitialized data
 *
 * @code
 * auto tensor = tn::empty<double>({100, 100});  // Fast allocation without initialization
 * @endcode
 */
template<typename T = double>
Tensor<T> empty(std::initializer_list<std::size_t> shape);

/**
 * @brief Create a tensor of zeros with the same shape as the input tensor.
 *
 * @tparam T The value type
 * @param other The tensor whose shape to match
 * @return A new tensor filled with zeros
 *
 * @code
 * auto original = tn::ones<double>({3, 4});
 * auto zeros = tn::zeros_like(original);  // Shape {3, 4}, filled with 0
 * @endcode
 */
template<typename T>
Tensor<T> zeros_like(const Tensor<T>& other);

/**
 * @brief Create a tensor of ones with the same shape as the input tensor.
 *
 * @tparam T The value type
 * @param other The tensor whose shape to match
 * @return A new tensor filled with ones
 *
 * @code
 * auto original = tn::zeros<double>({3, 4});
 * auto ones = tn::ones_like(original);  // Shape {3, 4}, filled with 1
 * @endcode
 */
template<typename T>
Tensor<T> ones_like(const Tensor<T>& other);

/**
 * @brief Create a tensor with the same shape as the input, filled with a value.
 *
 * @tparam T The value type
 * @param other The tensor whose shape to match
 * @param value The value to fill the tensor with
 * @return A new tensor filled with the specified value
 *
 * @code
 * auto original = tn::zeros<double>({3, 4});
 * auto filled = tn::full_like(original, 42.0);  // Shape {3, 4}, filled with 42.0
 * @endcode
 */
template<typename T>
Tensor<T> full_like(const Tensor<T>& other, T value);

/**
 * @brief Create a diagonal matrix from a 1D array of diagonal values.
 *
 * Creates a 2D square matrix with the given values on the diagonal and zeros elsewhere.
 *
 * @tparam T The value type (default: double)
 * @param diagonal The values to place on the diagonal
 * @return A 2D square tensor with the diagonal values
 *
 * @code
 * auto mat = tn::diag<double>({1.0, 2.0, 3.0});
 * // Result:
 * // 1 0 0
 * // 0 2 0
 * // 0 0 3
 * @endcode
 */
template<typename T = double>
Tensor<T> diag(std::initializer_list<T> diagonal);

/**
 * @brief Extract the diagonal from a 2D matrix.
 *
 * For a 2D tensor, extracts the main diagonal into a 1D tensor.
 *
 * @tparam T The value type
 * @param matrix The 2D tensor to extract the diagonal from
 * @return A 1D tensor containing the diagonal values
 *
 * @code
 * auto mat = tn::eye<double>(3);
 * auto diag_vals = tn::diag(mat);  // Result: {1.0, 1.0, 1.0}
 * @endcode
 */
template<typename T>
Tensor<T> diag(const Tensor<T>& matrix);

} // namespace igor::math::tensor

#include <igor/Math/TensorCreation.tpp>

// TensorCreation.h ends here
