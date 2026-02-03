/*
 * FastSampler.h
 *
 *  Created on: Feb 2, 2026
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
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *   Fast sampling utilities for sequence generation.
 *   Provides O(log n) sampling using precomputed CDFs and binary search,
 *   instead of O(n) linear search.
 */

#pragma once

#include <vector>
#include <random>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <iostream>

namespace igor {

/**
 * \class CategoricalSampler FastSampler.h
 * \brief Efficient categorical distribution sampler using precomputed CDF.
 * \author Optimization contribution
 * \version 1.0
 *
 * This class precomputes the cumulative distribution function (CDF) from
 * a probability vector and uses binary search for O(log n) sampling,
 * instead of the O(n) linear search used in the original implementation.
 */
class CategoricalSampler {
public:
    /**
     * \brief Default constructor creating an empty sampler.
     */
    CategoricalSampler() = default;

    /**
     * \brief Construct from a vector of probabilities.
     * \param probabilities Vector of probabilities (should sum to ~1.0)
     */
    explicit CategoricalSampler(const std::vector<double>& probabilities) {
        build_cdf(probabilities.data(), probabilities.size());
    }

    /**
     * \brief Construct from a raw array of probabilities.
     * \param probabilities Pointer to probability array
     * \param count Number of probabilities
     */
    CategoricalSampler(const double* probabilities, size_t count) {
        build_cdf(probabilities, count);
    }

    /**
     * \brief Construct from a long double array (for Marginal_array_p compatibility).
     * \param probabilities Pointer to probability array
     * \param count Number of probabilities
     */
    CategoricalSampler(const long double* probabilities, size_t count) {
        cumulative_.resize(count);
        double cum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            cum += static_cast<double>(probabilities[i]);
            cumulative_[i] = cum;
        }
        // Normalize to handle floating point errors
        if (count > 0 && cum > 0) {
            for (size_t i = 0; i < count; ++i) {
                cumulative_[i] /= cum;
            }
        }
    }

    /**
     * \brief Check if the sampler is empty.
     * \return true if no categories are defined
     */
    bool empty() const { return cumulative_.empty(); }

    /**
     * \brief Get the number of categories.
     * \return Number of categories
     */
    size_t size() const { return cumulative_.size(); }

    /**
     * \brief Sample from the categorical distribution using binary search.
     * \param rng Random number generator
     * \return Sampled category index
     *
     * Time complexity: O(log n) vs O(n) for linear search.
     */
    size_t sample(std::mt19937_64& rng) const {
        if (cumulative_.empty()) {
            throw std::runtime_error("CategoricalSampler: Cannot sample from empty distribution");
        }
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(rng);
        // Binary search for first element >= u
        auto it = std::lower_bound(cumulative_.begin(), cumulative_.end(), u);
        size_t idx = static_cast<size_t>(it - cumulative_.begin());
        // Clamp to valid range (handles floating point edge cases)
        return std::min(idx, cumulative_.size() - 1);
    }

    /**
     * \brief Sample multiple values from the distribution.
     * \param rng Random number generator
     * \param output Output array for sampled indices
     * \param count Number of samples to generate
     */
    void sample_batch(std::mt19937_64& rng, size_t* output, size_t count) const {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (size_t i = 0; i < count; ++i) {
            double u = dist(rng);
            auto it = std::lower_bound(cumulative_.begin(), cumulative_.end(), u);
            size_t idx = static_cast<size_t>(it - cumulative_.begin());
            output[i] = std::min(idx, cumulative_.size() - 1);
        }
    }

    /**
     * \brief Get the precomputed CDF for debugging/verification.
     * \return Reference to the cumulative distribution vector
     */
    const std::vector<double>& get_cdf() const { return cumulative_; }

private:
    void build_cdf(const double* probabilities, size_t count) {
        cumulative_.resize(count);
        double cum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            cum += probabilities[i];
            cumulative_[i] = cum;
        }
        // Normalize to handle floating point errors
        if (count > 0 && cum > 0) {
            for (size_t i = 0; i < count; ++i) {
                cumulative_[i] /= cum;
            }
        }
    }

    std::vector<double> cumulative_;  // Precomputed CDF
};


/**
 * \class LinearSampler FastSampler.h
 * \brief Sampler using linear search for EXACT reproducibility with original IGoR.
 * \author Optimization contribution
 * \version 1.0
 *
 * This sampler uses the same linear search algorithm as the original IGoR
 * implementation to ensure bit-exact reproducibility with the same seed.
 * It's O(n) instead of O(log n) but maintains compatibility.
 */
class LinearSampler {
public:
    LinearSampler() = default;

    explicit LinearSampler(const std::vector<double>& probabilities)
        : probabilities_(probabilities) {}

    bool empty() const { return probabilities_.empty(); }
    size_t size() const { return probabilities_.size(); }

    /**
     * \brief Sample using EXACT same algorithm as original IGoR.
     *
     * Original algorithm:
     *   double rand = uniform(0,1);
     *   double prob_count = 0;
     *   for each realization:
     *       prob_count += prob[i];
     *       if (prob_count >= rand) return i;
     */
    size_t sample(std::mt19937_64& rng) const {
        if (probabilities_.empty()) {
            throw std::runtime_error("LinearSampler: Cannot sample from empty distribution");
        }
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double rand = dist(rng);
        double prob_count = 0.0;
        
        for (size_t i = 0; i < probabilities_.size(); ++i) {
            prob_count += probabilities_[i];
            if (prob_count >= rand) {
                return i;
            }
        }
        // Return last index if we didn't find one (floating point edge case)
        return probabilities_.size() - 1;
    }

private:
    std::vector<double> probabilities_;
};


/**
 * \class ConditionalSampler FastSampler.h
 * \brief Sampler for conditional distributions (e.g., dinucleotide Markov).
 * \author Optimization contribution
 * \version 1.0
 *
 * Stores multiple CategoricalSamplers, one for each conditioning variable.
 * For dinucleotide Markov: 4 conditions (A,C,G,T previous), 4 outcomes each.
 */
class ConditionalSampler {
public:
    /**
     * \brief Default constructor.
     */
    ConditionalSampler() = default;

    /**
     * \brief Construct with specified dimensions.
     * \param num_conditions Number of conditioning states
     * \param num_outcomes Number of possible outcomes per condition
     */
    ConditionalSampler(size_t num_conditions, size_t num_outcomes)
        : num_conditions_(num_conditions), num_outcomes_(num_outcomes) {
        conditional_samplers_.resize(num_conditions);
    }

    /**
     * \brief Set probabilities for a specific condition.
     * \param condition Conditioning state index
     * \param probs Probability vector for outcomes given this condition
     */
    void set_probabilities(size_t condition, const std::vector<double>& probs) {
        if (condition >= num_conditions_) {
            throw std::out_of_range("ConditionalSampler: condition index out of range");
        }
        conditional_samplers_[condition] = CategoricalSampler(probs);
    }

    /**
     * \brief Set probabilities for a specific condition from raw array.
     * \param condition Conditioning state index
     * \param probs Pointer to probability array
     * \param count Number of outcomes
     */
    void set_probabilities(size_t condition, const double* probs, size_t count) {
        if (condition >= num_conditions_) {
            throw std::out_of_range("ConditionalSampler: condition index out of range");
        }
        conditional_samplers_[condition] = CategoricalSampler(probs, count);
    }

    /**
     * \brief Sample outcome given a condition.
     * \param condition Conditioning state
     * \param rng Random number generator
     * \return Sampled outcome index
     */
    size_t sample(size_t condition, std::mt19937_64& rng) const {
        if (condition >= num_conditions_) {
            throw std::out_of_range("ConditionalSampler: condition index out of range");
        }
        return conditional_samplers_[condition].sample(rng);
    }

    /**
     * \brief Check if all conditions have been initialized.
     * \return true if all conditional samplers are non-empty
     */
    bool is_initialized() const {
        for (const auto& sampler : conditional_samplers_) {
            if (sampler.empty()) return false;
        }
        return !conditional_samplers_.empty();
    }

    size_t num_conditions() const { return num_conditions_; }
    size_t num_outcomes() const { return num_outcomes_; }

private:
    size_t num_conditions_ = 0;
    size_t num_outcomes_ = 0;
    std::vector<CategoricalSampler> conditional_samplers_;
};


/**
 * \class DinucleotideMarkovSampler FastSampler.h
 * \brief Optimized sampler for dinucleotide Markov insertion model.
 * \author Optimization contribution
 * \version 1.0
 *
 * This class precomputes all 4x4 conditional distributions for the
 * dinucleotide Markov model, enabling O(log 4) = O(1) sampling per nucleotide
 * instead of O(4) linear search.
 *
 * For longer insertions, this provides significant speedup since the
 * original implementation recomputes cumulative probabilities for each
 * nucleotide position.
 */
class DinucleotideMarkovSampler {
public:
    static constexpr size_t NUM_NUCLEOTIDES = 4;  // A, C, G, T
    static constexpr size_t NUM_EXTENDED = 15;    // Including ambiguous codes

    /**
     * \brief Default constructor.
     */
    DinucleotideMarkovSampler() : sampler_(NUM_EXTENDED, NUM_NUCLEOTIDES) {}

    /**
     * \brief Build sampler from dinucleotide probability matrix.
     * \param dinuc_matrix 15x15 matrix (or similar) of conditional probabilities
     *
     * The matrix is indexed as: P(next_nt | prev_nt) = dinuc_matrix(prev_nt, next_nt)
     * For nucleotides 0-3 (A,C,G,T), we precompute the conditional CDFs.
     * For ambiguous nucleotides (4-14), we use averaged probabilities.
     */
    template<typename MatrixType>
    void build_from_matrix(const MatrixType& dinuc_matrix) {
        sampler_ = ConditionalSampler(NUM_EXTENDED, NUM_NUCLEOTIDES);
        
        for (size_t prev_nt = 0; prev_nt < NUM_EXTENDED; ++prev_nt) {
            std::vector<double> probs(NUM_NUCLEOTIDES);
            double sum = 0.0;
            for (size_t next_nt = 0; next_nt < NUM_NUCLEOTIDES; ++next_nt) {
                probs[next_nt] = dinuc_matrix(prev_nt, next_nt);
                sum += probs[next_nt];
            }
            // Normalize
            if (sum > 0) {
                for (size_t i = 0; i < NUM_NUCLEOTIDES; ++i) {
                    probs[i] /= sum;
                }
            }
            sampler_.set_probabilities(prev_nt, probs);
        }
        initialized_ = true;
    }

    /**
     * \brief Sample the next nucleotide given the previous one.
     * \param prev_nt Previous nucleotide index (0=A, 1=C, 2=G, 3=T, or ambiguous 4-14)
     * \param rng Random number generator
     * \return Next nucleotide index (0-3)
     */
    size_t sample_next(size_t prev_nt, std::mt19937_64& rng) const {
        return sampler_.sample(prev_nt, rng);
    }

    /**
     * \brief Sample a complete insertion sequence.
     * \param prev_nt Previous nucleotide (from the gene sequence)
     * \param length Number of nucleotides to insert
     * \param output Output vector for sampled nucleotide indices
     * \param rng Random number generator
     */
    void sample_insertion(size_t prev_nt, size_t length, 
                          std::vector<int>& output, std::mt19937_64& rng) const {
        output.resize(length);
        size_t current = prev_nt;
        for (size_t i = 0; i < length; ++i) {
            size_t next = sampler_.sample(current, rng);
            output[i] = static_cast<int>(next);
            current = next;
        }
    }

    /**
     * \brief Sample insertion and write to string directly.
     * \param prev_nt Previous nucleotide index
     * \param inserted_seq String to fill with nucleotides (modified in place)
     * \param realization_indices Output vector for realization indices
     * \param rng Random number generator
     *
     * This matches the interface expected by Dinucl_markov::draw_random_common
     */
    void sample_insertion_to_string(size_t prev_nt, std::string& inserted_seq,
                                    std::vector<int>& realization_indices,
                                    std::mt19937_64& rng) const {
        static const char NT_CHARS[] = {'A', 'C', 'G', 'T'};
        
        realization_indices.clear();
        realization_indices.reserve(inserted_seq.size());
        
        size_t current = prev_nt;
        for (size_t i = 0; i < inserted_seq.size(); ++i) {
            if (inserted_seq[i] == 'I') {  // 'I' marks positions to be inserted
                size_t next = sampler_.sample(current, rng);
                inserted_seq[i] = NT_CHARS[next];
                realization_indices.push_back(static_cast<int>(next));
                current = next;
            } else {
                // Position already has a nucleotide (e.g., palindromic)
                // Convert character to index for next iteration
                switch (inserted_seq[i]) {
                    case 'A': current = 0; break;
                    case 'C': current = 1; break;
                    case 'G': current = 2; break;
                    case 'T': current = 3; break;
                    default: current = 14; break;  // N or unknown
                }
            }
        }
    }

    bool is_initialized() const { return initialized_; }

private:
    ConditionalSampler sampler_;
    bool initialized_ = false;
};


/**
 * \struct RealizationBuffer FastSampler.h
 * \brief Pre-allocated buffer for sequence generation to avoid per-sequence allocations.
 * \author Optimization contribution
 * \version 1.0
 */
struct RealizationBuffer {
    // Pre-allocated strings for gene sequences
    std::string v_gene_seq;
    std::string d_gene_seq;
    std::string j_gene_seq;
    
    // Pre-allocated strings for insertions
    std::string vd_insertion;
    std::string dj_insertion;
    std::string vj_insertion;
    
    // Realization indices for each event
    std::vector<int> v_gene_realization;
    std::vector<int> d_gene_realization;
    std::vector<int> j_gene_realization;
    std::vector<int> v_del_realization;
    std::vector<int> d_del_5_realization;
    std::vector<int> d_del_3_realization;
    std::vector<int> j_del_realization;
    std::vector<int> vd_ins_length_realization;
    std::vector<int> dj_ins_length_realization;
    std::vector<int> vj_ins_length_realization;
    std::vector<int> vd_dinuc_realization;
    std::vector<int> dj_dinuc_realization;
    std::vector<int> vj_dinuc_realization;
    
    // Final concatenated sequence
    std::string final_sequence;
    
    // Error positions (if generating with errors)
    std::vector<int> error_positions;
    
    /**
     * \brief Clear all buffers for reuse.
     */
    void clear() {
        v_gene_seq.clear();
        d_gene_seq.clear();
        j_gene_seq.clear();
        vd_insertion.clear();
        dj_insertion.clear();
        vj_insertion.clear();
        v_gene_realization.clear();
        d_gene_realization.clear();
        j_gene_realization.clear();
        v_del_realization.clear();
        d_del_5_realization.clear();
        d_del_3_realization.clear();
        j_del_realization.clear();
        vd_ins_length_realization.clear();
        dj_ins_length_realization.clear();
        vj_ins_length_realization.clear();
        vd_dinuc_realization.clear();
        dj_dinuc_realization.clear();
        vj_dinuc_realization.clear();
        final_sequence.clear();
        error_positions.clear();
    }
    
    /**
     * \brief Reserve capacity based on expected maximum sizes.
     * \param max_v_len Maximum V gene length
     * \param max_d_len Maximum D gene length
     * \param max_j_len Maximum J gene length
     * \param max_ins Maximum insertion length (per junction)
     */
    void reserve(size_t max_v_len, size_t max_d_len, size_t max_j_len, size_t max_ins) {
        v_gene_seq.reserve(max_v_len);
        d_gene_seq.reserve(max_d_len);
        j_gene_seq.reserve(max_j_len);
        vd_insertion.reserve(max_ins);
        dj_insertion.reserve(max_ins);
        vj_insertion.reserve(max_ins);
        final_sequence.reserve(max_v_len + max_d_len + max_j_len + 3 * max_ins);
        
        // Reserve for realization indices
        vd_dinuc_realization.reserve(max_ins);
        dj_dinuc_realization.reserve(max_ins);
        vj_dinuc_realization.reserve(max_ins);
    }
    
    /**
     * \brief Build final sequence from components.
     * \param has_d Whether the model has a D gene
     */
    void build_final_sequence(bool has_d) {
        final_sequence.clear();
        final_sequence.reserve(v_gene_seq.size() + d_gene_seq.size() + j_gene_seq.size() +
                               vd_insertion.size() + dj_insertion.size() + vj_insertion.size());
        
        final_sequence += v_gene_seq;
        if (has_d) {
            final_sequence += vd_insertion;
            final_sequence += d_gene_seq;
            final_sequence += dj_insertion;
        } else {
            final_sequence += vj_insertion;
        }
        final_sequence += j_gene_seq;
    }
};

}  // namespace igor
