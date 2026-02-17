#pragma once

namespace igor::model {

// ─── Full-shape constructor ────────────────────────────────────────────

template <typename T>
MarkovHandler<T>::MarkovHandler(std::string name, std::vector<std::size_t> shape)
    : MarginalHandler<T>(std::move(name))
    , shape_(std::move(shape))
    , parameters_(shape_)
    , accumulator_(shape_)
{
    // Uniform initialization: each row [from, :, p1, p2, ...] sums to 1
    // All elements set to 1/n_to_states ensures proper row-wise normalization
    const std::size_t n_to_states = shape_[1];  // Number of "to" states per row
    const T uniform = T(1) / static_cast<T>(n_to_states);
    std::fill(parameters_.begin(), parameters_.end(), uniform);
    std::fill(accumulator_.begin(), accumulator_.end(), T(0));
}

// ─── 2D convenience constructor ────────────────────────────────────────

template <typename T>
MarkovHandler<T>::MarkovHandler(std::string name, std::size_t n_states)
    : MarkovHandler(std::move(name), std::vector<std::size_t>{n_states, n_states})
{}

// ─── EM operations ─────────────────────────────────────────────────────

template <typename T>
void MarkovHandler<T>::reset_accumulator() {
    std::fill(accumulator_.begin(), accumulator_.end(), T(0));
}

template <typename T>
void MarkovHandler<T>::maximize_likelihood() {
    if (parameters_.ndim() == 2) {
        // Fast path for unconditional Markov (no parents): row-wise normalize
        std::size_t n = shape_[0];
        for (std::size_t from = 0; from < n; ++from) {
            T row_sum = T(0);
            for (std::size_t to = 0; to < n; ++to) {
                row_sum += accumulator_(from, to);
            }
            if (row_sum > T(1e-15)) {
                for (std::size_t to = 0; to < n; ++to) {
                    parameters_(from, to) = accumulator_(from, to) / row_sum;
                }
            }
        }
    } else {
        // Multi-dimensional: normalize along axis 1 (to-state)
        // Each slice [from, *, parent1, ...] sums to 1
        math::linalg::normalize_axis(accumulator_, parameters_, 1);
    }
}

// ─── Sampling ──────────────────────────────────────────────────────────

template <typename T>
std::size_t MarkovHandler<T>::sample(
    std::mt19937_64& generator,
    const std::vector<std::size_t>& parent_indices) const
{
    // For Markov handlers, parent_indices should contain:
    // parent_indices[0] = from_state (previous state in the sequence)
    // parent_indices[1..n] = parent event realizations (optional)

    if (parent_indices.empty()) {
        // First nucleotide: sample from the marginal distribution
        // computed by summing each row of the transition matrix
        std::size_t n_states = parameters_.shape()[0];
        std::size_t n_to_states = parameters_.shape()[1];

        // Compute marginal: for each 'from' state, sum the row
        std::vector<T> marginal(n_states, T(0));
        for (std::size_t s = 0; s < n_states; ++s) {
            const T* row = parameters_.data() + s * n_to_states;
            for (std::size_t t = 0; t < n_to_states; ++t) {
                marginal[s] += row[t];
            }
        }

        // Normalize
        T total = T(0);
        for (auto v : marginal) total += v;
        if (total > T(0)) {
            for (auto& v : marginal) v /= total;
        }

        // Sample
        std::uniform_real_distribution<T> dist(T(0), T(1));
        T random_val = dist(generator);
        T cumsum = T(0);
        for (std::size_t s = 0; s < n_states; ++s) {
            cumsum += marginal[s];
            if (random_val <= cumsum) {
                return s;
            }
        }
        return n_states - 1;
    }

    std::size_t from_state = parent_indices[0];
    std::size_t n_to_states = parameters_.shape()[1];

    if (from_state >= shape_[0]) {
        throw std::out_of_range("from_state index out of range for " + this->name());
    }

    // Get the transition probability row for this from_state
    // If there are more parent indices (conditional Markov with upstream dependencies),
    // we need to slice appropriately
    const T* transition_probs = nullptr;

    if (parameters_.ndim() == 2) {
        // No additional parents: direct row access
        transition_probs = parameters_.data() + from_state * n_to_states;
    } else {
        // Multi-dimensional: slice for [from_state, :, parent1, ...]
        // Calculate linear offset for parent indices at dimensions [2:]
        std::size_t stride = n_to_states;
        std::size_t offset = from_state * stride;  // from_state offset

        for (int d = parameters_.ndim() - 1; d >= 2; --d) {
            offset += parent_indices[d - 1] * stride;
            stride *= parameters_.shape()[d];
        }
        transition_probs = parameters_.data() + offset;
    }

    // Sample next state using cumulative sum
    // CRITICAL: RNG call order must match legacy implementation
    std::uniform_real_distribution<T> dist(T(0), T(1));
    T random_val = dist(generator);
    T cumsum = T(0);

    for (std::size_t to = 0; to < n_to_states; ++to) {
        cumsum += transition_probs[to];
        if (random_val <= cumsum) {
            return to;  // Return "to" state index
        }
    }

    // Safety fallback: return last state if rounding errors occur
    return n_to_states - 1;
}

// ─── I/O ───────────────────────────────────────────────────────────────

template <typename T>
void MarkovHandler<T>::write_parameters(std::ostream& out) const {
    out << this->name() << "\n";
    for (std::size_t i = 0; i < parameters_.size(); ++i) {
        out << parameters_.data()[i];
        if (i + 1 < parameters_.size()) out << ",";
    }
    out << "\n";
}

template <typename T>
void MarkovHandler<T>::read_parameters(std::istream& in) {
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
