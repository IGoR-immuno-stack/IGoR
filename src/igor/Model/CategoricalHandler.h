/*
 * CategoricalHandler.h
 *
 *  Created on: Feb 10, 2026
 *
 *  Handler for categorical distributions (Gene_choice, Deletion, Insertion).
 *
 *  The parameter tensor has shape [n_realizations, parent1, parent2, ...]
 *  where the first dimension is the event's own realizations and the
 *  remaining dimensions correspond to parent event realizations.
 *
 *  M-step: normalize along axis 0 (sum over own realizations = 1
 *          for each combination of parent values).
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
class CategoricalHandler : public MarginalHandler<T> {
public:
    /// Construct from full shape: {n_realizations, parent1_size, parent2_size, ...}
    CategoricalHandler(std::string name, std::vector<std::size_t> shape);

    /// Convenience: 1D (no parents) — shape = {n_realizations}
    CategoricalHandler(std::string name, std::size_t n_realizations);

    // Tensor access
    const math::Tensor<T>& parameters() const override { return parameters_; }
    math::Tensor<T>& parameters() override { return parameters_; }
    const math::Tensor<T>& accumulator() const override { return accumulator_; }
    math::Tensor<T>& accumulator() override { return accumulator_; }

    // EM operations
    void reset_accumulator() override;
    void maximize_likelihood() override;

    // I/O
    void write_parameters(std::ostream& out) const override;
    void read_parameters(std::istream& in) override;

    // Accessors
    std::size_t n_realizations() const { return parameters_.shape()[0]; }
    const std::vector<std::size_t>& shape() const { return shape_; }

private:
    std::vector<std::size_t> shape_;
    math::Tensor<T> parameters_;   // [n_realizations, parent1, ...]
    math::Tensor<T> accumulator_;  // [n_realizations, parent1, ...]
};

} // namespace igor::model

#include <igor/Model/CategoricalHandler.tpp>
