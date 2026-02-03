/*
 * FastGenerator.h
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
 *   Fast parallel sequence generator with optimized sampling and I/O.
 */

#pragma once

#include <igor/Core/FastSampler.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <fstream>
#include <sstream>
#include <mutex>
#include <atomic>
#include <thread>
#include <queue>

namespace igor {

/**
 * \enum GenerationMode
 * \brief Mode for sequence generation, trading off speed vs reproducibility.
 */
enum class GenerationMode {
    ExactMatch,     ///< Produces sequences identical to original IGoR (slower)
    MaxSpeed        ///< Maximum performance using precomputed samplers (faster, different sequences)
};

/**
 * \struct GenerationConfig FastGenerator.h
 * \brief Configuration for fast sequence generation.
 */
struct GenerationConfig {
    size_t num_threads = 0;            ///< Number of threads (0 = auto-detect)
    size_t batch_size = 10000;         ///< Sequences per batch for I/O
    size_t io_buffer_size = 1 << 20;   ///< 1 MB I/O buffer
    bool generate_errors = false;      ///< Whether to apply error model
    uint64_t seed = 0;                 ///< Random seed (0 = random seed)
    GenerationMode mode = GenerationMode::ExactMatch;  ///< Generation mode
    
    // Output options
    bool output_realizations = true;   ///< Output realization indices
    bool output_sequences_only = false; ///< Only output sequences (faster)
    
    /**
     * \brief Get effective number of threads.
     * \return Number of threads to use
     */
    size_t effective_threads() const {
        if (num_threads == 0) {
            return std::max(1u, std::thread::hardware_concurrency());
        }
        return num_threads;
    }
};


/**
 * \class BufferedWriter FastGenerator.h
 * \brief Thread-safe buffered file writer for high-throughput I/O.
 */
class BufferedWriter {
public:
    /**
     * \brief Construct a buffered writer.
     * \param path Output file path
     * \param buffer_size Buffer size in bytes
     */
    BufferedWriter(const std::string& path, size_t buffer_size = 1 << 20)
        : buffer_size_(buffer_size) {
        file_.open(path, std::ios::out | std::ios::binary);
        if (!file_.is_open()) {
            throw std::runtime_error("BufferedWriter: Cannot open file: " + path);
        }
        buffer_.reserve(buffer_size);
    }

    ~BufferedWriter() {
        flush();
        file_.close();
    }

    /**
     * \brief Write a line to the buffer.
     * \param line Text to write
     */
    void write(const std::string& line) {
        std::lock_guard<std::mutex> lock(mutex_);
        buffer_ += line;
        if (buffer_.size() >= buffer_size_) {
            flush_internal();
        }
    }

    /**
     * \brief Write multiple lines at once (more efficient).
     * \param lines Vector of lines to write
     */
    void write_batch(const std::vector<std::string>& lines) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (const auto& line : lines) {
            buffer_ += line;
        }
        if (buffer_.size() >= buffer_size_) {
            flush_internal();
        }
    }

    /**
     * \brief Flush buffer to file.
     */
    void flush() {
        std::lock_guard<std::mutex> lock(mutex_);
        flush_internal();
    }

private:
    void flush_internal() {
        if (!buffer_.empty()) {
            file_.write(buffer_.data(), buffer_.size());
            buffer_.clear();
        }
    }

    std::ofstream file_;
    std::string buffer_;
    size_t buffer_size_;
    std::mutex mutex_;
};


/**
 * \class FastGenerator FastGenerator.h
 * \brief High-performance parallel sequence generator.
 * \author Optimization contribution
 * \version 1.0
 *
 * This class provides optimized sequence generation with:
 * - Precomputed sampling distributions (CDF + binary search)
 * - Parallel generation using multiple threads
 * - Buffered I/O to reduce disk access overhead
 * - Memory pooling to avoid per-sequence allocations
 *
 * Target: 100x speedup over original implementation (~200,000 seq/sec)
 */
class FastGenerator {
public:
    /**
     * \brief Construct a fast generator from model parameters.
     * \param model_parms Model parameters defining the recombination graph
     * \param marginals Marginal probabilities for each event
     * \param config Generation configuration
     */
    FastGenerator(const Model_Parms& model_parms,
                  const Model_marginals& marginals,
                  const GenerationConfig& config = GenerationConfig());

    /**
     * \brief Generate sequences to files.
     * \param num_sequences Number of sequences to generate
     * \param sequence_file Output file for sequences
     * \param realization_file Output file for realizations
     */
    void generate(size_t num_sequences,
                  const std::string& sequence_file,
                  const std::string& realization_file);

    /**
     * \brief Generate sequences with custom callback.
     * \param num_sequences Number of sequences to generate
     * \param callback Function called for each generated sequence
     */
    void generate(size_t num_sequences,
                  std::function<void(size_t idx, const std::string& seq,
                                     const std::vector<std::vector<int>>& realizations)> callback);

    /**
     * \brief Generate sequences to memory (for small samples).
     * \param num_sequences Number of sequences to generate
     * \return Vector of (sequence, realizations) pairs
     */
    std::vector<std::pair<std::string, std::vector<std::vector<int>>>>
    generate_to_memory(size_t num_sequences);

    /**
     * \brief Set progress callback for monitoring generation.
     * \param callback Function called periodically with (done, total)
     */
    void set_progress_callback(std::function<void(size_t done, size_t total)> callback) {
        progress_callback_ = callback;
    }

    /**
     * \brief Get generation statistics.
     */
    struct Stats {
        double sequences_per_second;
        double total_time_seconds;
        size_t total_sequences;
    };
    Stats get_stats() const { return stats_; }

private:
    /**
     * \struct ThreadContext
     * \brief Per-thread resources to avoid contention.
     */
    struct ThreadContext {
        std::mt19937_64 rng;
        RealizationBuffer buffer;
        std::stringstream output_seq;
        std::stringstream output_real;
        std::vector<std::string> seq_batch;
        std::vector<std::string> real_batch;
        
        // Per-thread model data to avoid race conditions
        // Each thread needs its own copy because events have mutable internal state
        std::unique_ptr<Model_Parms> model_parms;
        std::queue<std::shared_ptr<Rec_Event>> model_queue;
        std::unordered_map<Rec_Event_name, int> index_map;
        std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map;
    };

    /**
     * \brief Prepare samplers from model marginals.
     */
    void prepare_samplers();

    /**
     * \brief Generate a single sequence using ExactMatch mode (calls original event methods).
     * \param ctx Thread context with RNG and buffers
     * \param seq_idx Sequence index for output
     */
    void generate_single_sequence(ThreadContext& ctx, size_t seq_idx);
    
    /**
     * \brief Generate a single sequence using MaxSpeed mode (precomputed samplers).
     * \param ctx Thread context with RNG and buffers
     * \param seq_idx Sequence index for output
     * 
     * This method uses precomputed CDF-based samplers for maximum performance.
     * Sequences are statistically equivalent but not identical to the original.
     */
    void generate_single_sequence_fast(ThreadContext& ctx, size_t seq_idx);

    /**
     * \brief Worker function for parallel generation.
     * \param thread_id Thread identifier
     * \param start_idx Starting sequence index
     * \param count Number of sequences to generate
     * \param seq_writer Buffered writer for sequences
     * \param real_writer Buffered writer for realizations (may be null)
     */
    void generate_worker(size_t thread_id, size_t start_idx, size_t count,
                         BufferedWriter* seq_writer, BufferedWriter* real_writer);

    /**
     * \brief Format realizations for output.
     * \param buffer Realization buffer with all event realizations
     * \return Formatted string for CSV output
     */
    std::string format_realizations(const RealizationBuffer& buffer) const;

    // Model data - stored as copies to ensure consistent internal state
    Model_Parms model_parms_;
    Model_marginals marginals_;
    GenerationConfig config_;
    
    // Precomputed samplers for each event type
    // LinearSampler: O(n) for exact reproducibility with original IGoR (ExactMatch mode)
    // CategoricalSampler: O(log n) binary search for maximum speed (MaxSpeed mode)
    std::vector<LinearSampler> gene_choice_samplers_;
    std::vector<LinearSampler> deletion_samplers_;
    std::vector<LinearSampler> insertion_length_samplers_;
    std::vector<DinucleotideMarkovSampler> dinuc_samplers_;
    
    // Fast samplers using binary search (O(log n)) for MaxSpeed mode
    std::vector<CategoricalSampler> fast_gene_choice_samplers_;
    std::vector<CategoricalSampler> fast_deletion_samplers_;
    std::vector<CategoricalSampler> fast_insertion_length_samplers_;
    
    // Named fast gene choice samplers - conditional on parent gene choices
    // V gene is unconditional
    CategoricalSampler fast_v_gene_sampler_;
    // J gene is conditional on V gene: fast_j_gene_samplers_[v_gene_idx]
    std::vector<CategoricalSampler> fast_j_gene_samplers_;
    // D gene is conditional on V gene AND J gene: 
    // fast_d_gene_samplers_[v_gene_idx * num_j_genes_ + j_gene_idx]
    std::vector<CategoricalSampler> fast_d_gene_samplers_;
    
    // Named fast deletion samplers - ALL are conditional on gene choice
    // V deletion is conditional on V gene choice
    std::vector<CategoricalSampler> fast_v_del_samplers_;
    // D deletion samplers indexed by D gene choice (conditional on D gene)
    // fast_d_del_5_samplers_[d_gene_index] gives sampler for that D gene
    std::vector<CategoricalSampler> fast_d_del_5_samplers_;
    // D_3' deletion is conditional on BOTH D gene choice AND D_5' deletion
    // fast_d_del_3_samplers_[d_gene_idx * num_d_del_5 + d_del_5_idx]
    std::vector<CategoricalSampler> fast_d_del_3_samplers_;
    // J deletion is conditional on J gene choice
    std::vector<CategoricalSampler> fast_j_del_samplers_;
    
    // Number of gene choices and deletion options (for conditional sampler indexing)
    size_t num_v_genes_ = 0;
    size_t num_d_genes_ = 0;
    size_t num_j_genes_ = 0;
    size_t num_d_del_5_ = 0;
    
    // Cached gene sequences for fast lookup (indexed by sampler result)
    std::vector<std::string> v_gene_sequences_;
    std::vector<std::string> d_gene_sequences_;
    std::vector<std::string> j_gene_sequences_;
    
    // Realization index mappings (sampler index -> original realization index)
    std::vector<int> v_realization_indices_;
    std::vector<int> d_realization_indices_;
    std::vector<int> j_realization_indices_;
    
    // Cached deletion values (indexed by sampler result)
    std::vector<int> v_del_values_;
    std::vector<int> d_del_5_values_;
    std::vector<int> d_del_3_values_;
    std::vector<int> j_del_values_;
    
    // Deletion index mappings
    std::vector<int> v_del_indices_;
    std::vector<int> d_del_5_indices_;
    std::vector<int> d_del_3_indices_;
    std::vector<int> j_del_indices_;
    
    // Cached insertion length values
    std::vector<int> vd_ins_values_;
    std::vector<int> dj_ins_values_;
    std::vector<int> vj_ins_values_;
    
    // Insertion index mappings
    std::vector<int> vd_ins_indices_;
    std::vector<int> dj_ins_indices_;
    std::vector<int> vj_ins_indices_;
    
    // Event order and indices
    std::vector<std::tuple<Event_type, Gene_class, Seq_side>> event_order_;
    std::unordered_map<Rec_Event_name, int> index_map_;
    std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map_;
    std::queue<std::shared_ptr<Rec_Event>> model_queue_;  // Cached model queue
    
    // Thread contexts
    std::vector<ThreadContext> thread_contexts_;
    
    // Progress tracking
    std::atomic<size_t> sequences_generated_{0};
    std::function<void(size_t, size_t)> progress_callback_;
    
    // Statistics
    Stats stats_;
    
    // Model characteristics
    bool has_d_gene_ = true;
};

}  // namespace igor
