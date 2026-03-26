#pragma once

#include <igor/Math/TensorCreation.h>

namespace igor::model {

// ─── Constructor ───────────────────────────────────────────────────────

template <typename T>
MarkovInferenceHandler<T>::MarkovInferenceHandler(std::string name, igor::index_type uid, math::Tensor<T>& weights)
    : InferenceHandler<T>(std::move(name), uid)
    , m_weights(weights)
    , m_accumulator(math::tensor::zeros_like(weights))
{
}

// ─── Tensor access ─────────────────────────────────────────────────────

template <typename T>
const math::Tensor<T>& MarkovInferenceHandler<T>::weights(void) const
{
    return m_weights;
}

template <typename T>
math::Tensor<T>& MarkovInferenceHandler<T>::weights(void)
{
    return m_weights;
}

template <typename T>
const math::Tensor<T>& MarkovInferenceHandler<T>::accumulator(void) const
{
    return m_accumulator;
}

template <typename T>
math::Tensor<T>& MarkovInferenceHandler<T>::accumulator(void)
{
    return m_accumulator;
}

// ─── Accessors ─────────────────────────────────────────────────────────

template <typename T>
std::size_t MarkovInferenceHandler<T>::stateCount(void) const
{
    return m_weights.shape()[0];
}

template <typename T>
const std::vector<std::size_t>& MarkovInferenceHandler<T>::shape(void) const
{
    return m_weights.shape();
}

// ─── EM operations ─────────────────────────────────────────────────────

template <typename T>
void MarkovInferenceHandler<T>::resetAccumulator(void)
{
    std::fill(m_accumulator.begin(), m_accumulator.end(), T(0));
}

template <typename T>
void MarkovInferenceHandler<T>::maximizeLikelihood(void)
{
    if (m_weights.ndim() == 2) {
        // Fast path for unconditional Markov (no parents): row-wise normalize
        std::size_t n = m_weights.shape()[0];
        for (std::size_t from = 0; from < n; ++from) {
            T row_sum = T(0);
            for (std::size_t to = 0; to < n; ++to) {
                row_sum += m_accumulator(from, to);
            }
            if (row_sum > T(1e-15)) {
                for (std::size_t to = 0; to < n; ++to) {
                    m_weights(from, to) = m_accumulator(from, to) / row_sum;
                }
            } else {
                // Zero accumulator row → zero weights row (matches normalize_axis)
                for (std::size_t to = 0; to < n; ++to) {
                    m_weights(from, to) = T(0);
                }
            }
        }
    } else {
        // Multi-dimensional: normalize along LAST axis (to-state)
        // Shape is [parent1, parent2, ..., from_state, to_state]
        std::size_t to_state_axis = m_weights.ndim() - 1;
        math::linalg::normalize_axis(m_accumulator, m_weights, to_state_axis);
    }
}

} // namespace igor::model
