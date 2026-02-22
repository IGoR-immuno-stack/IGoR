/*
 * CategoricalInferenceHandler.h
 *
 *  Created on: Feb 21, 2026
 *
 *  Handler for categorical distributions (Gene_choice, Deletion, Insertion).
 *
 *  The weight tensor has shape [n_realizations, parent1, parent2, ...]
 *  where the first dimension is the event's own realizations and the
 *  remaining dimensions correspond to parent event realizations.
 *
 *  The handler borrows a mutable reference to the probability tensor stored
 *  in RecombinationModel. It only owns the accumulator tensor needed for
 *  the E-step.
 *
 *  M-step: normalize along axis 0 (sum over own realizations = 1
 *          for each combination of parent values).
 */

#pragma once

#include <igor/Model/InferenceHandler.h>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>

#include <string>
#include <vector>

namespace igor::model {

template <typename T = double>
class CategoricalInferenceHandler : public InferenceHandler<T> {
public:
    /// Construct from name, uid, and a mutable reference to the model's weight tensor.
    CategoricalInferenceHandler(std::string name, igor::index_type uid, math::Tensor<T>& weights);

    // Tensor access
    const math::Tensor<T>& weights(void) const override;
          math::Tensor<T>& weights(void)       override;

    const math::Tensor<T>& accumulator(void) const override;
          math::Tensor<T>& accumulator(void)       override;

    // EM operations
    void resetAccumulator(void) override;
    void maximizeLikelihood(void) override;

    // Accessors
    std::size_t realizationCount(void) const;
    const std::vector<std::size_t>& shape(void) const;

private:
    math::Tensor<T>& m_weights;      // borrowed from RecombinationModel
    math::Tensor<T>  m_accumulator;   // owned, same shape as m_weights
};

} // namespace igor::model

#include <igor/Model/CategoricalInferenceHandler.tpp>
