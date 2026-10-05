#pragma once

#include <igor/Math/TensorCreation.h>

#include <numeric>

#include <igor/Model/Forward.h>

namespace igor::inference {
using igor::model::RecombinationModel;
using igor::model::Navigator;
using igor::model::Topology;
using igor::model::SampledScenario;
using igor::model::SampledEvent;

// ─── Constructor ───────────────────────────────────────────────────────

template <typename T>
CategoricalInferenceHandler<T>::CategoricalInferenceHandler(std::string name, igor::core::index_type uid, math::Tensor<T>& weights)
    : InferenceHandler<T>(std::move(name), uid)
    , m_weights(weights)
    , m_accumulator(math::tensor::zeros_like(weights))
{
}

// ─── Tensor access ─────────────────────────────────────────────────────

template <typename T>
const math::Tensor<T>& CategoricalInferenceHandler<T>::weights(void) const
{
    return m_weights;
}

template <typename T>
math::Tensor<T>& CategoricalInferenceHandler<T>::weights(void)
{
    return m_weights;
}

template <typename T>
const math::Tensor<T>& CategoricalInferenceHandler<T>::accumulator(void) const
{
    return m_accumulator;
}

template <typename T>
math::Tensor<T>& CategoricalInferenceHandler<T>::accumulator(void)
{
    return m_accumulator;
}

// ─── Accessors ─────────────────────────────────────────────────────────

template <typename T>
std::size_t CategoricalInferenceHandler<T>::realizationCount(void) const
{
    // Shape is [parent1, parent2, ..., child], so the own realizations are the last axis,
    // the one maximizeLikelihood() normalises. Same convention as
    // CategoricalSamplingHandler, which reads shape().back() too.
    return m_weights.shape().back();
}

template <typename T>
const std::vector<std::size_t>& CategoricalInferenceHandler<T>::shape(void) const
{
    return m_weights.shape();
}

// ─── EM operations ─────────────────────────────────────────────────────

template <typename T>
void CategoricalInferenceHandler<T>::resetAccumulator(void)
{
    std::fill(m_accumulator.begin(), m_accumulator.end(), T(0));
}

template <typename T>
void CategoricalInferenceHandler<T>::maximizeLikelihood(void)
{
    if (m_weights.ndim() == 1) {
        // Fast path for unconditional events (no parents)
        T total = math::linalg::sum(m_accumulator);
        constexpr T tolerance = std::numeric_limits<T>::epsilon() * T(10);
        if (total > tolerance) {
            for (std::size_t i = 0; i < m_weights.size(); ++i) {
                m_weights.data()[i] = m_accumulator.data()[i] / total;
            }
        } else {
            // Zero accumulator → zero weights (matches normalize_axis behavior)
            std::fill(m_weights.begin(), m_weights.end(), T(0));
        }
    } else {
        // Multi-dimensional: normalize along LAST axis (child dimension)
        // Shape is [parent1, parent2, ..., child], normalize over child
        std::size_t child_axis = m_weights.ndim() - 1;
        math::linalg::normalize_axis(m_accumulator, m_weights, child_axis);
    }
}

} // namespace igor::inference
