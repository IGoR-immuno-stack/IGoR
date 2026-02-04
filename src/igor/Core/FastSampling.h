/*
 * FastSampling.h
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

#pragma once

#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <mutex>
#include <atomic>
#include <thread>
#include <functional>
#include <fstream>
#include <sstream>

namespace igor {
namespace fast {

/**
 * \class CategoricalSampler
 * \brief Fast categorical sampling using precomputed CDF and binary search.
 *
 * Precomputes cumulative distribution function (CDF) for O(log n) sampling
 * instead of O(n) linear search. For very frequent sampling, also supports
 * the alias method for O(1) sampling.
 */
class CategoricalSampler {
public:
    CategoricalSampler() = default;

    /**
     * \brief Initialize sampler with probabilities
     * \param probabilities Vector of probabilities (need not be normalized)
     * \param use_alias Whether to use alias method (O(1)) vs binary search (O(log n))
     */
    void initialize(const std::vector<double>& probabilities, bool use_alias = false);

    /**
     * \brief Initialize from raw array
     * \param probs Pointer to probability array
     * \param n Number of elements
     * \param use_alias Whether to use alias method
     */
    void initialize(const double* probs, size_t n, bool use_alias = false);

    /**
     * \brief Initialize from long double array (model marginals)
     * \param probs Pointer to probability array
     * \param n Number of elements
     * \param use_alias Whether to use alias method
     */
    void initialize(const long double* probs, size_t n, bool use_alias = false);

    /**
     * \brief Sample an index from the distribution
     * \param rng Random number generator
     * \return Sampled index
     */
    template<typename RNG>
    size_t sample(RNG& rng) const {
        if (use_alias_method_) {
            return sample_alias(rng);
        }
        return sample_binary_search(rng);
    }

    /**
     * \brief Get the probability of a specific outcome
     */
    double probability(size_t index) const {
        return index < probabilities_.size() ? probabilities_[index] : 0.0;
    }

    /**
     * \brief Get the number of categories
     */
    size_t size() const { return probabilities_.size(); }

    /**
     * \brief Check if sampler is initialized
     */
    bool is_initialized() const { return !cdf_.empty(); }

private:
    template<typename RNG>
    size_t sample_binary_search(RNG& rng) const {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(rng);
        // Binary search for the first CDF value >= u
        auto it = std::lower_bound(cdf_.begin(), cdf_.end(), u);
        return std::min(static_cast<size_t>(it - cdf_.begin()), cdf_.size() - 1);
    }

    template<typename RNG>
    size_t sample_alias(RNG& rng) const {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::uniform_int_distribution<size_t> int_dist(0, alias_.size() - 1);
        size_t i = int_dist(rng);
        double u = dist(rng);
        return (u < prob_alias_[i]) ? i : alias_[i];
    }

    void build_cdf();
    void build_alias();

    std::vector<double> probabilities_;
    std::vector<double> cdf_;

    // Alias method structures
    bool use_alias_method_ = false;
    std::vector<size_t> alias_;
    std::vector<double> prob_alias_;
};


/**
 * \class ConditionalSampler
 * \brief Sampler for conditional distributions P(X|Y).
 *
 * Stores multiple CategoricalSamplers, one for each conditioning value.
 * Optimized memory layout for cache efficiency.
 */
class ConditionalSampler {
public:
    ConditionalSampler() = default;

    /**
     * \brief Initialize with conditional probabilities
     * \param probs 2D probability matrix [condition][outcome]
     * \param use_alias Whether to use alias method for sampling
     */
    void initialize(const std::vector<std::vector<double>>& probs, bool use_alias = false);

    /**
     * \brief Initialize from flat array with dimensions
     * \param probs Flat array of probabilities
     * \param num_conditions Number of conditioning values
     * \param num_outcomes Number of possible outcomes per condition
     * \param use_alias Whether to use alias method
     */
    void initialize(const double* probs, size_t num_conditions, size_t num_outcomes, bool use_alias = false);

    /**
     * \brief Initialize from long double flat array
     */
    void initialize(const long double* probs, size_t num_conditions, size_t num_outcomes, bool use_alias = false);

    /**
     * \brief Sample given a conditioning value
     * \param condition The conditioning value
     * \param rng Random number generator
     * \return Sampled outcome
     */
    template<typename RNG>
    size_t sample(size_t condition, RNG& rng) const {
        return samplers_[condition].sample(rng);
    }

    /**
     * \brief Get the number of conditioning values
     */
    size_t num_conditions() const { return samplers_.size(); }

    /**
     * \brief Check if sampler is initialized
     */
    bool is_initialized() const { return !samplers_.empty(); }

private:
    std::vector<CategoricalSampler> samplers_;
};


/**
 * \class DinucleotideMarkovSampler
 * \brief Optimized sampler for dinucleotide Markov model.
 *
 * Precomputes all 4x4 (or 15x4 with ambiguous) conditional distributions
 * for efficient sequence generation. Replaces per-nucleotide probability
 * computation with a single lookup and sample.
 */
class DinucleotideMarkovSampler {
public:
    // Standard nucleotide indices (matching IGoR's nt2int)
    static constexpr size_t NT_A = 0;
    static constexpr size_t NT_C = 1;
    static constexpr size_t NT_G = 2;
    static constexpr size_t NT_T = 3;
    static constexpr size_t NUM_STANDARD_NT = 4;
    static constexpr size_t NUM_ALL_NT = 15;  // Including ambiguous

    DinucleotideMarkovSampler() = default;

    /**
     * \brief Initialize from dinucleotide probability matrix
     * \param dinuc_probs 15x15 or 4x4 probability matrix
     * \param size Size of the matrix (4 or 15)
     */
    void initialize(const double* dinuc_probs, size_t size);

    /**
     * \brief Initialize from Matrix<double> (IGoR's matrix type)
     */
    template<typename MatrixT>
    void initialize_from_matrix(const MatrixT& matrix) {
        size_t n = matrix.get_n_rows();
        std::vector<double> probs(n * n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                probs[i * n + j] = static_cast<double>(matrix(i, j));
            }
        }
        initialize(probs.data(), n);
    }

    /**
     * \brief Sample next nucleotide given previous nucleotide
     * \param prev_nt Previous nucleotide (integer code)
     * \param rng Random number generator
     * \return Next nucleotide (integer code 0-3)
     */
    template<typename RNG>
    int sample_next(int prev_nt, RNG& rng) const {
        // Handle ambiguous nucleotides by sampling from uniform distribution
        if (prev_nt < 0 || prev_nt >= static_cast<int>(conditional_sampler_.num_conditions())) {
            std::uniform_int_distribution<int> dist(0, 3);
            return dist(rng);
        }
        return static_cast<int>(conditional_sampler_.sample(static_cast<size_t>(prev_nt), rng));
    }

    /**
     * \brief Generate entire inserted sequence
     * \param first_nt First nucleotide (from previous gene)
     * \param length Number of nucleotides to generate
     * \param output Output buffer for nucleotide sequence (integer codes)
     * \param rng Random number generator
     */
    template<typename RNG>
    void generate_sequence(int first_nt, size_t length, int* output, RNG& rng) const {
        if (length == 0) return;

        int prev = first_nt;
        for (size_t i = 0; i < length; ++i) {
            output[i] = sample_next(prev, rng);
            prev = output[i];
        }
    }

    /**
     * \brief Generate sequence as string
     */
    template<typename RNG>
    std::string generate_sequence_str(int first_nt, size_t length, RNG& rng) const {
        static const char NT_CHARS[] = {'A', 'C', 'G', 'T'};
        std::string result(length, 'N');
        int prev = first_nt;
        for (size_t i = 0; i < length; ++i) {
            int nt = sample_next(prev, rng);
            result[i] = NT_CHARS[nt];
            prev = nt;
        }
        return result;
    }

    /**
     * \brief Check if sampler is initialized
     */
    bool is_initialized() const { return conditional_sampler_.is_initialized(); }

private:
    ConditionalSampler conditional_sampler_;
};


/**
 * \class RealizationBuffer
 * \brief Pre-allocated buffer for storing sequence generation results.
 *
 * Avoids repeated allocations during sequence generation.
 */
class RealizationBuffer {
public:
    explicit RealizationBuffer(size_t max_seq_length = 1000)
        : max_length_(max_seq_length) {
        v_gene_seq_.reserve(max_length_);
        d_gene_seq_.reserve(max_length_);
        j_gene_seq_.reserve(max_length_);
        vd_ins_seq_.reserve(max_length_ / 4);
        dj_ins_seq_.reserve(max_length_ / 4);
        vj_ins_seq_.reserve(max_length_ / 4);
        final_seq_.reserve(max_length_);
        realizations_.reserve(20);  // Typical number of events
    }

    void clear() {
        v_gene_seq_.clear();
        d_gene_seq_.clear();
        j_gene_seq_.clear();
        vd_ins_seq_.clear();
        dj_ins_seq_.clear();
        vj_ins_seq_.clear();
        final_seq_.clear();
        realizations_.clear();
    }

    // Gene sequences
    std::string v_gene_seq_;
    std::string d_gene_seq_;
    std::string j_gene_seq_;

    // Insertion sequences
    std::string vd_ins_seq_;
    std::string dj_ins_seq_;
    std::string vj_ins_seq_;

    // Final assembled sequence
    std::string final_seq_;

    // Event realizations (index for each event)
    std::vector<std::vector<int>> realizations_;

private:
    size_t max_length_;
};


/**
 * \class BufferPool
 * \brief Thread-safe pool of RealizationBuffers.
 *
 * Reduces memory allocation overhead by reusing buffers across
 * sequence generation calls.
 */
class BufferPool {
public:
    explicit BufferPool(size_t pool_size = 0, size_t max_seq_length = 1000)
        : max_seq_length_(max_seq_length) {
        for (size_t i = 0; i < pool_size; ++i) {
            pool_.push_back(std::make_unique<RealizationBuffer>(max_seq_length_));
        }
    }

    /**
     * \brief Acquire a buffer from the pool
     * \return Unique pointer to buffer (returns to pool on destruction)
     */
    std::unique_ptr<RealizationBuffer> acquire() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (pool_.empty()) {
            return std::make_unique<RealizationBuffer>(max_seq_length_);
        }
        auto buffer = std::move(pool_.back());
        pool_.pop_back();
        buffer->clear();
        return buffer;
    }

    /**
     * \brief Release a buffer back to the pool
     */
    void release(std::unique_ptr<RealizationBuffer> buffer) {
        std::lock_guard<std::mutex> lock(mutex_);
        pool_.push_back(std::move(buffer));
    }

    /**
     * \brief Get pool statistics
     */
    size_t available() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return pool_.size();
    }

private:
    mutable std::mutex mutex_;
    std::vector<std::unique_ptr<RealizationBuffer>> pool_;
    size_t max_seq_length_;
};


/**
 * \class BufferedWriter
 * \brief Buffered I/O for sequence output.
 *
 * Batches sequence writes to reduce I/O overhead.
 * Thread-safe with batched writes to minimize contention.
 */
class BufferedWriter {
public:
    explicit BufferedWriter(const std::string& filename,
                           size_t buffer_size = 1024 * 1024)  // 1MB default
        : filename_(filename), buffer_size_(buffer_size) {
        buffer_.reserve(buffer_size_);
    }

    ~BufferedWriter() {
        flush();
    }

    /**
     * \brief Write a line to the buffer (thread-safe)
     */
    void write_line(const std::string& line) {
        std::lock_guard<std::mutex> lock(mutex_);
        buffer_ += line;
        buffer_ += '\n';
        if (buffer_.size() >= buffer_size_) {
            flush_internal();
        }
    }

    /**
     * \brief Write multiple lines at once (more efficient for batched writes)
     */
    void write_lines(const std::vector<std::string>& lines) {
        if (lines.empty()) return;

        // Calculate total size
        size_t total_size = 0;
        for (const auto& line : lines) {
            total_size += line.size() + 1;  // +1 for newline
        }

        // Build batch string
        std::string batch;
        batch.reserve(total_size);
        for (const auto& line : lines) {
            batch += line;
            batch += '\n';
        }

        // Single lock for entire batch
        std::lock_guard<std::mutex> lock(mutex_);
        buffer_ += batch;
        if (buffer_.size() >= buffer_size_) {
            flush_internal();
        }
    }

    /**
     * \brief Write formatted sequence data
     */
    void write_sequence(size_t index, const std::string& sequence) {
        std::ostringstream oss;
        oss << index << ";" << sequence;
        write_line(oss.str());
    }

    /**
     * \brief Flush buffer to disk
     */
    void flush() {
        std::lock_guard<std::mutex> lock(mutex_);
        flush_internal();
    }

    /**
     * \brief Get number of bytes written
     */
    size_t bytes_written() const {
        return total_bytes_written_.load();
    }

private:
    void flush_internal() {
        if (buffer_.empty()) return;

        std::ofstream file(filename_, std::ios::app | std::ios::binary);
        if (file) {
            file.write(buffer_.data(), buffer_.size());
            total_bytes_written_ += buffer_.size();
        }
        buffer_.clear();
    }

    std::string filename_;
    size_t buffer_size_;
    std::string buffer_;
    std::mutex mutex_;
    std::atomic<size_t> total_bytes_written_{0};
};


/**
 * \struct ThreadContext
 * \brief Per-thread context for parallel sequence generation.
 */
struct ThreadContext {
    std::mt19937_64 rng;
    std::unique_ptr<RealizationBuffer> buffer;

    explicit ThreadContext(uint64_t seed, size_t max_seq_length = 1000)
        : rng(seed), buffer(std::make_unique<RealizationBuffer>(max_seq_length)) {}
};


/**
 * \class FastGeneratorConfig
 * \brief Configuration for fast sequence generation.
 */
struct FastGeneratorConfig {
    size_t num_threads = 0;  // 0 = auto-detect
    size_t batch_size = 10000;  // Sequences per batch
    size_t buffer_pool_size = 0;  // 0 = num_threads * 2
    size_t io_buffer_size = 1024 * 1024;  // 1MB
    bool show_progress = true;
    uint64_t base_seed = 0;  // 0 = random seed
};


/**
 * \brief Callback type for generated sequences
 * \param seq_index Index of the generated sequence
 * \param sequence The generated sequence string
 * \param realizations Event realizations
 */
using SequenceCallback = std::function<void(
    size_t seq_index,
    const std::string& sequence,
    const std::vector<std::vector<int>>& realizations
)>;


/**
 * \brief Progress callback type
 * \param completed Number of sequences completed
 * \param total Total number of sequences
 */
using ProgressCallback = std::function<void(size_t completed, size_t total)>;


}  // namespace fast
}  // namespace igor
