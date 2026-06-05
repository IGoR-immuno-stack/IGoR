/*
 * MarkovInferenceHandler.h
 *
 *  Created on: Feb 21, 2026
 *
 *  Handler for Markov transition matrices (Dinucl_markov).
 *
 *  The weight tensor has shape [n_states, n_states, parent1, parent2, ...]
 *  where the first two dimensions are the transition matrix (from, to)
 *  and the remaining dimensions correspond to parent event realizations.
 *
 *  The handler borrows a mutable reference to the probability tensor stored
 *  in RecombinationModel. It only owns the accumulator tensor needed for
 *  the E-step.
 *
 *  M-step: normalize along axis 1 (sum over "to" states = 1
 *          for each "from" state and parent combination).
 */

#pragma once

#include <igor/Model/InferenceHandler.h>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>

#include <string>
#include <vector>

namespace igor::model {

template <typename T = double>
class MarkovInferenceHandler : public InferenceHandler<T> {
public:
    /// Construct from name, uid, and a mutable reference to the model's weight tensor.
    MarkovInferenceHandler(std::string name, igor::index_type uid, math::Tensor<T>& weights);

    // Tensor access
    const math::Tensor<T>& weights(void) const override;
          math::Tensor<T>& weights(void)       override;

    const math::Tensor<T>& accumulator(void) const override;
          math::Tensor<T>& accumulator(void)       override;

    // EM operations
    void resetAccumulator(void) override;
    void maximizeLikelihood(void) override;

    // Accessors
    std::size_t stateCount(void) const;
    const std::vector<std::size_t>& shape(void) const;

private:
    math::Tensor<T>& m_weights;      // borrowed from RecombinationModel
    math::Tensor<T>  m_accumulator;   // owned, same shape as m_weights
};

} // namespace igor::model

#include <igor/Model/MarkovInferenceHandler.tpp>
