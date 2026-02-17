#pragma once

#include <numeric>

namespace igor::model {

// ─── Full-shape constructor ────────────────────────────────────────────

template <typename T>
CategoricalHandler<T>::CategoricalHandler(std::string name, std::vector<std::size_t> shape)
    : MarginalHandler<T>(std::move(name))
    , shape_(std::move(shape))
    , parameters_(shape_)
    , accumulator_(shape_)
{
    // Uniform initialization: each parent slice [:, p1, p2, ...] sums to 1
    // All elements set to 1/n_realizations ensures proper normalization
    const std::size_t n_realizations = shape_[0];
    const T uniform = T(1) / static_cast<T>(n_realizations);
    std::fill(parameters_.begin(), parameters_.end(), uniform);
    std::fill(accumulator_.begin(), accumulator_.end(), T(0));
}

// ─── 1D convenience constructor ────────────────────────────────────────

template <typename T>
CategoricalHandler<T>::CategoricalHandler(std::string name, std::size_t n_realizations)
    : CategoricalHandler(std::move(name), std::vector<std::size_t>{n_realizations})
{}

// ─── EM operations ─────────────────────────────────────────────────────

template <typename T>
void CategoricalHandler<T>::reset_accumulator() {
    std::fill(accumulator_.begin(), accumulator_.end(), T(0));
}

template <typename T>
void CategoricalHandler<T>::maximize_likelihood() {
    if (parameters_.ndim() == 1) {
        // Fast path for unconditional events (no parents)
        T total = math::linalg::sum(accumulator_);
        constexpr T tolerance = std::numeric_limits<T>::epsilon() * T(10);
        if (total > tolerance) {
            for (std::size_t i = 0; i < parameters_.size(); ++i) {
                parameters_.data()[i] = accumulator_.data()[i] / total;
            }
        }
    } else {
        // Multi-dimensional: normalize along axis 0
        // Each slice [*, p1, p2, ...] sums to 1
        math::linalg::normalize_axis(accumulator_, parameters_, 0);
    }
}

// ─── Sampling ──────────────────────────────────────────────────────────

template <typename T>
std::size_t CategoricalHandler<T>::sample(
    std::mt19937_64& generator,
    const std::vector<std::size_t>& parent_indices) const
{
    // Validate parent indices match tensor shape
    if (parent_indices.size() != parameters_.ndim() - 1) {
        throw std::invalid_argument("Parent indices count mismatch for " + this->name());
    }

    // Get the probability slice for this parent combination
    // For unconditional (ndim() == 1), parent_indices is empty
    const T* prob_slice = nullptr;
    std::size_t n_realizations = parameters_.shape()[0];

    if (parameters_.ndim() == 1) {
        // No parents: direct access
        prob_slice = parameters_.data();
    } else {
        // Multi-dimensional: slice access
        // Calculate linear index for [*, parent_indices...]
        std::size_t stride = 1;
        std::size_t offset = 0;
        for (int d = parameters_.ndim() - 1; d >= 1; --d) {
            offset += parent_indices[d - 1] * stride;
            stride *= parameters_.shape()[d];
        }
        prob_slice = parameters_.data() + offset * n_realizations;
    }

    // Sample using cumulative sum
    // CRITICAL: RNG call order must match legacy implementation
    std::uniform_real_distribution<T> dist(T(0), T(1));
    T random_val = dist(generator);
    T cumsum = T(0);

    for (std::size_t i = 0; i < n_realizations; ++i) {
        cumsum += prob_slice[i];
        if (random_val <= cumsum) {
            return i;  // Return realization index
        }
    }

    // Safety fallback: return last realization if rounding errors occur
    return n_realizations - 1;
}

// ─── I/O ───────────────────────────────────────────────────────────────

template <typename T>
void CategoricalHandler<T>::write_parameters(std::ostream& out) const {
    out << this->name() << "\n";
    for (std::size_t i = 0; i < parameters_.size(); ++i) {
        out << parameters_.data()[i];
        if (i + 1 < parameters_.size()) out << ",";
    }
    out << "\n";
}

template <typename T>
void CategoricalHandler<T>::read_parameters(std::istream& in) {
    std::string line;
    std::getline(in, line); // skip name line
    std::getline(in, line); // values line

    if (!in.good()) {
        throw std::runtime_error("Failed to read parameters for " + this->name_);
    }

    std::size_t idx = 0;
    std::size_t pos = 0;
    while (pos < line.size() && idx < parameters_.size()) {
        auto comma = line.find(',', pos);
        if (comma == std::string::npos) comma = line.size();
        parameters_.data()[idx] = static_cast<T>(std::stold(line.substr(pos, comma - pos)));
        pos = comma + 1;
        ++idx;
    }

    if (idx != parameters_.size()) {
        throw std::runtime_error("Incorrect number of parameters for " + this->name_);
    }
}

} // namespace igor::model
