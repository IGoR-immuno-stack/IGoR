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
