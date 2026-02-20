/*
 * MarkovHandler.h
 *
 *  Created on: Feb 10, 2026
 *
 *  Handler for Markov transition matrices (Dinucl_markov).
 *
 *  The parameter tensor has shape [n_states, n_states, parent1, parent2, ...]
 *  where the first two dimensions are the transition matrix (from, to)
 *  and the remaining dimensions correspond to parent event realizations.
 *
 *  M-step: normalize along axis 1 (sum over "to" states = 1
 *          for each "from" state and parent combination).
 */

#pragma once

#include <igor/Model/MarginalHandler.h>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>

#include <algorithm>
#include <string>
#include <vector>

namespace igor::model {

template <typename T = double>
class MarkovHandler : public MarginalHandler<T> {
public:
    /// Construct from full shape: {n_states, n_states, parent1_size, ...}
    MarkovHandler(std::string name, std::vector<std::size_t> shape);

    /// Convenience: 2D (no parents) — shape = {n_states, n_states}
    MarkovHandler(std::string name, std::size_t n_states);

    // Tensor access
    const math::Tensor<T>& parameters() const override { return parameters_; }
    math::Tensor<T>& parameters() override { return parameters_; }
    const math::Tensor<T>& accumulator() const override { return accumulator_; }
    math::Tensor<T>& accumulator() override { return accumulator_; }

    // EM operations
    void reset_accumulator() override;
    void maximize_likelihood() override;

    // Sampling
    std::size_t sample(
        std::mt19937_64& generator,
        const std::vector<std::size_t>& parent_indices = {}) const override;

    // I/O
    void write_parameters(std::ostream& out) const override;
    void read_parameters(std::istream& in) override;

    // Accessors
    std::size_t n_states() const { return shape_[0]; }
    const std::vector<std::size_t>& shape() const { return shape_; }

private:
    std::vector<std::size_t> shape_;
    math::Tensor<T> parameters_;   // [n_states, n_states, parent1, ...]
    math::Tensor<T> accumulator_;  // [n_states, n_states, parent1, ...]
};

} // namespace igor::model

#include <igor/Model/MarkovHandler.tpp>
