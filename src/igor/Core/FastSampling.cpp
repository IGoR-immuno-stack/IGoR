/*
 * FastSampling.cpp
 *
 *  Created on: Feb 4, 2026
 *      Performance optimization for sequence generation
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <igor/Core/FastSampling.h>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <queue>

namespace igor {
namespace fast {

//==============================================================================
// CategoricalSampler Implementation
//==============================================================================

void CategoricalSampler::initialize(const std::vector<double>& probabilities, bool use_alias) {
    initialize(probabilities.data(), probabilities.size(), use_alias);
}

void CategoricalSampler::initialize(const double* probs, size_t n, bool use_alias) {
    if (n == 0) {
        probabilities_.clear();
        cdf_.clear();
        return;
    }

    probabilities_.resize(n);
    std::copy(probs, probs + n, probabilities_.begin());

    // Normalize probabilities
    double sum = std::accumulate(probabilities_.begin(), probabilities_.end(), 0.0);
    if (sum > 0) {
        for (auto& p : probabilities_) {
            p /= sum;
        }
    } else {
        // Uniform distribution if all zeros
        double uniform = 1.0 / static_cast<double>(n);
        std::fill(probabilities_.begin(), probabilities_.end(), uniform);
    }

    use_alias_method_ = use_alias;
    build_cdf();

    if (use_alias) {
        build_alias();
    }
}

void CategoricalSampler::initialize(const long double* probs, size_t n, bool use_alias) {
    if (n == 0) {
        probabilities_.clear();
        cdf_.clear();
        return;
    }

    probabilities_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        probabilities_[i] = static_cast<double>(probs[i]);
    }

    // Normalize probabilities
    double sum = std::accumulate(probabilities_.begin(), probabilities_.end(), 0.0);
    if (sum > 0) {
        for (auto& p : probabilities_) {
            p /= sum;
        }
    } else {
        // Uniform distribution if all zeros
        double uniform = 1.0 / static_cast<double>(n);
        std::fill(probabilities_.begin(), probabilities_.end(), uniform);
    }

    use_alias_method_ = use_alias;
    build_cdf();

    if (use_alias) {
        build_alias();
    }
}

void CategoricalSampler::build_cdf() {
    size_t n = probabilities_.size();
    cdf_.resize(n);

    if (n == 0) return;

    cdf_[0] = probabilities_[0];
    for (size_t i = 1; i < n; ++i) {
        cdf_[i] = cdf_[i - 1] + probabilities_[i];
    }

    // Ensure CDF ends at exactly 1.0 (numerical stability)
    cdf_.back() = 1.0;
}

void CategoricalSampler::build_alias() {
    // Vose's Alias Method for O(1) sampling
    size_t n = probabilities_.size();
    if (n == 0) return;

    alias_.resize(n);
    prob_alias_.resize(n);

    // Scale probabilities by n
    std::vector<double> scaled(n);
    for (size_t i = 0; i < n; ++i) {
        scaled[i] = probabilities_[i] * static_cast<double>(n);
    }

    // Partition into small and large
    std::queue<size_t> small, large;
    for (size_t i = 0; i < n; ++i) {
        if (scaled[i] < 1.0) {
            small.push(i);
        } else {
            large.push(i);
        }
    }

    // Build alias table
    while (!small.empty() && !large.empty()) {
        size_t s = small.front(); small.pop();
        size_t l = large.front(); large.pop();

        prob_alias_[s] = scaled[s];
        alias_[s] = l;

        scaled[l] = scaled[l] + scaled[s] - 1.0;

        if (scaled[l] < 1.0) {
            small.push(l);
        } else {
            large.push(l);
        }
    }

    // Handle remaining (due to numerical issues)
    while (!large.empty()) {
        size_t l = large.front(); large.pop();
        prob_alias_[l] = 1.0;
        alias_[l] = l;
    }

    while (!small.empty()) {
        size_t s = small.front(); small.pop();
        prob_alias_[s] = 1.0;
        alias_[s] = s;
    }
}


//==============================================================================
// ConditionalSampler Implementation
//==============================================================================

void ConditionalSampler::initialize(const std::vector<std::vector<double>>& probs, bool use_alias) {
    samplers_.clear();
    samplers_.reserve(probs.size());

    for (const auto& row : probs) {
        CategoricalSampler sampler;
        sampler.initialize(row, use_alias);
        samplers_.push_back(std::move(sampler));
    }
}

void ConditionalSampler::initialize(const double* probs, size_t num_conditions,
                                    size_t num_outcomes, bool use_alias) {
    samplers_.clear();
    samplers_.reserve(num_conditions);

    for (size_t i = 0; i < num_conditions; ++i) {
        CategoricalSampler sampler;
        sampler.initialize(probs + i * num_outcomes, num_outcomes, use_alias);
        samplers_.push_back(std::move(sampler));
    }
}

void ConditionalSampler::initialize(const long double* probs, size_t num_conditions,
                                    size_t num_outcomes, bool use_alias) {
    samplers_.clear();
    samplers_.reserve(num_conditions);

    for (size_t i = 0; i < num_conditions; ++i) {
        CategoricalSampler sampler;
        sampler.initialize(probs + i * num_outcomes, num_outcomes, use_alias);
        samplers_.push_back(std::move(sampler));
    }
}


//==============================================================================
// DinucleotideMarkovSampler Implementation
//==============================================================================

void DinucleotideMarkovSampler::initialize(const double* dinuc_probs, size_t size) {
    if (size == 0) return;

    // Build conditional probability matrix
    // dinuc_probs[i * size + j] = P(next=j | prev=i)
    // But IGoR uses dinuc_probs[i + rows * j] (column-major like Fortran)
    // We need to handle this properly

    std::vector<std::vector<double>> cond_probs(size, std::vector<double>(NUM_STANDARD_NT, 0.0));

    for (size_t prev = 0; prev < size; ++prev) {
        double sum = 0.0;
        // Only sample standard nucleotides (0-3), but condition can be ambiguous
        for (size_t next = 0; next < NUM_STANDARD_NT; ++next) {
            // IGoR's Matrix uses column-major: matrix(i,j) = array[i + rows*j]
            double p = dinuc_probs[prev + size * next];
            cond_probs[prev][next] = p;
            sum += p;
        }
        // Normalize if needed
        if (sum > 0 && std::abs(sum - 1.0) > 1e-10) {
            for (size_t next = 0; next < NUM_STANDARD_NT; ++next) {
                cond_probs[prev][next] /= sum;
            }
        }
    }

    conditional_sampler_.initialize(cond_probs, true);  // Use alias for O(1)
}


}  // namespace fast
}  // namespace igor
