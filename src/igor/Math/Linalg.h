/*
 * Linalg.h
 *
 *  Created on: Feb 08, 2026
 *      Author: IGoR Agent
 */

#pragma once

#include <igor/Math/Tensor.h>

#include <algorithm>
#include <numeric>
#include <concepts>
#include <array>

namespace igor::math::linalg {

// -----------------------------------------------------------------------------
// Helper Interface
// -----------------------------------------------------------------------------

template <typename In1, typename In2, typename Out, typename Func>
void apply_binary(In1 in1, In2 in2, Out out, Func f);

template <typename In, typename Out, typename Func>
void apply_unary(In in, Out out, Func f);

template <typename In, typename T, typename Func>
void reduce(In in, T& acc, Func f);

// -----------------------------------------------------------------------------
// Element-wise Operations (Core)
// -----------------------------------------------------------------------------

/**
 * \brief Element-wise Add: Out = In1 + In2 (mdspan)
 */
template <typename In1, typename In2, typename Out>
void add(In1 in1, In2 in2, Out out);

template <typename In1, typename In2, typename Out>
void subtract(In1 in1, In2 in2, Out out);

template <typename In1, typename In2, typename Out>
void multiply(In1 in1, In2 in2, Out out);

template <typename In1, typename In2, typename Out>
void divide(In1 in1, In2 in2, Out out);

/**
 * \brief Scale: Out = In * Scalar (mdspan)
 */
template <typename In, typename Scalar, typename Out>
void scale(In in, Scalar s, Out out);

// -----------------------------------------------------------------------------
// Element-wise Operations (Tensor Overloads)
// -----------------------------------------------------------------------------

template <typename T>
void add(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out);

template <typename T>
void subtract(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out);

template <typename T>
void multiply(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out);

template <typename T>
void divide(const Tensor<T>& in1, const Tensor<T>& in2, Tensor<T>& out);

template <typename T, typename Scalar>
void scale(const Tensor<T>& in, Scalar s, Tensor<T>& out);


// -----------------------------------------------------------------------------
// Reductions (Core)
// -----------------------------------------------------------------------------

template <typename In>
auto sum(In in);

template <typename In>
auto max(In in);

template <typename In>
auto min(In in);

template <typename In>
auto argmax(In in);

// -----------------------------------------------------------------------------
// Reductions (Tensor Overloads)
// -----------------------------------------------------------------------------

template <typename T>
T sum(const Tensor<T>& in);

template <typename T>
T max(const Tensor<T>& in);

template <typename T>
T min(const Tensor<T>& in);

template <typename T>
auto argmax(const Tensor<T>& in);

} // namespace igor::math::linalg

#include <igor/Math/Linalg.tpp>
