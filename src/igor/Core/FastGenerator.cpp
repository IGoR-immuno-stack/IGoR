/*
 * FastGenerator.cpp
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
 */

#include <igor/Core/FastGenerator.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Dinuclmarkov.h>

#include <chrono>
#include <iostream>
#include <iomanip>

namespace igor {

FastGenerator::FastGenerator(const Model_Parms& model_parms,
                             const Model_marginals& marginals,
                             const GenerationConfig& config)
    : model_parms_(model_parms)
    , marginals_(marginals)
    , config_(config)
{
    // Determine if model has D gene
    auto event_list = model_parms_.get_event_list();
    has_d_gene_ = false;
    for (const auto& event : event_list) {
        if (event->get_class() == D_gene || event->get_class() == VD_genes ||
            event->get_class() == DJ_genes || event->get_class() == VDJ_genes) {
            has_d_gene_ = true;
            break;
        }
    }
    
    // Initialize thread contexts
    size_t num_threads = config_.effective_threads();
    thread_contexts_.resize(num_threads);
    
    // Seed RNGs for each thread
    uint64_t base_seed = config_.seed;
    if (base_seed == 0) {
        // Generate random seed from high-resolution clock
        auto now = std::chrono::high_resolution_clock::now();
        base_seed = now.time_since_epoch().count();
    }
    
    // Use different seeds for each thread for independent streams
    for (size_t i = 0; i < num_threads; ++i) {
        thread_contexts_[i].rng.seed(base_seed + i * 1000000007ULL);
        thread_contexts_[i].buffer.reserve(500, 50, 100, 50);  // Typical max sizes
        
        // Only create per-thread model copies for ExactMatch mode
        // MaxSpeed mode uses precomputed samplers and doesn't need model copies
        if (config_.mode == GenerationMode::ExactMatch) {
            thread_contexts_[i].model_parms = std::make_unique<Model_Parms>(model_parms_);
        }
    }
    
    // Prepare samplers from marginals
    prepare_samplers();
}

void FastGenerator::prepare_samplers() {
    // Get model queue, index map, and offset map - computed ONCE and cached
    model_queue_ = model_parms_.get_model_queue();
    index_map_ = marginals_.get_index_map(model_parms_, model_queue_);
    offset_map_ = marginals_.get_offsets_map(model_parms_, model_queue_);
    
    // Get the marginal array
    const auto& marginal_array = marginals_.marginal_array_smart_p;
    
    // CRITICAL: Call update_event_internal_probas for each event to populate
    // internal probability structures like dinuc_proba_matrix in Dinucl_markov
    // This mirrors what generate_unique_sequence does in the original code
    std::queue<std::shared_ptr<Rec_Event>> init_queue = model_queue_;
    while (!init_queue.empty()) {
        init_queue.front()->update_event_internal_probas(marginal_array, index_map_);
        init_queue.pop();
    }
    
    // Process events in MODEL QUEUE ORDER (same as original generate_unique_sequence)
    std::queue<std::shared_ptr<Rec_Event>> queue_copy = model_queue_;
    
    while (!queue_copy.empty()) {
        auto event = queue_copy.front();
        queue_copy.pop();
        
        Event_type type = event->get_type();
        Gene_class gene_class = event->get_class();
        Seq_side side = event->get_side();
        Rec_Event_name name = event->get_name();
        
        event_order_.push_back(std::make_tuple(type, gene_class, side));
        
        int base_index = index_map_.at(name);
        const auto& realizations = event->get_realizations_map();
        
        switch (type) {
            case GeneChoice_t: {
                // Build arrays sorted by INDEX for reproducibility
                size_t num_realizations = realizations.size();
                
                // Build sorted index mapping
                std::vector<std::pair<int, std::string>> sorted_reals;
                sorted_reals.reserve(num_realizations);
                for (const auto& [key, real] : realizations) {
                    sorted_reals.emplace_back(real.index, real.value_str);
                }
                std::sort(sorted_reals.begin(), sorted_reals.end(),
                    [](const auto& a, const auto& b) { return a.first < b.first; });
                
                std::vector<int> index_mapping;
                std::vector<std::string> seq_mapping;
                index_mapping.reserve(num_realizations);
                seq_mapping.reserve(num_realizations);
                for (const auto& [idx, seq] : sorted_reals) {
                    index_mapping.push_back(idx);
                    seq_mapping.push_back(seq);
                }
                
                // Build unconditional sampler (for LinearSampler/ExactMatch mode)
                std::vector<double> probs(num_realizations);
                for (size_t i = 0; i < num_realizations; ++i) {
                    probs[i] = static_cast<double>(marginal_array[base_index + index_mapping[i]]);
                }
                gene_choice_samplers_.push_back(LinearSampler(probs));
                fast_gene_choice_samplers_.push_back(CategoricalSampler(probs));
                
                // Store mappings and build CONDITIONAL fast samplers
                if (gene_class == V_gene) {
                    v_gene_sequences_ = std::move(seq_mapping);
                    v_realization_indices_ = std::move(index_mapping);
                    num_v_genes_ = num_realizations;
                    // V gene is unconditional - use the same sampler
                    fast_v_gene_sampler_ = CategoricalSampler(probs);
                } else if (gene_class == J_gene) {
                    j_gene_sequences_ = std::move(seq_mapping);
                    j_realization_indices_ = index_mapping;
                    num_j_genes_ = num_realizations;
                    
                    // J gene depends on V gene - build conditional samplers
                    // Check if there's an offset from V gene
                    int v_offset = 0;
                    if (offset_map_.count(name) > 0) {
                        for (const auto& [parent_event, offset] : offset_map_.at(name)) {
                            if (parent_event->get_class() == V_gene) {
                                v_offset = offset;
                                break;
                            }
                        }
                    }
                    
                    if (v_offset > 0 && num_v_genes_ > 0) {
                        // Build conditional samplers for each V gene choice
                        // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                        fast_j_gene_samplers_.resize(num_v_genes_);
                        for (size_t v_idx = 0; v_idx < num_v_genes_; ++v_idx) {
                            std::vector<double> cond_probs(num_realizations);
                            // Use v_idx (position) directly, not v_realization_indices_[v_idx]
                            int cond_base = base_index + static_cast<int>(v_idx) * v_offset;
                            for (size_t i = 0; i < num_realizations; ++i) {
                                cond_probs[i] = static_cast<double>(marginal_array[cond_base + index_mapping[i]]);
                            }
                            fast_j_gene_samplers_[v_idx] = CategoricalSampler(cond_probs);
                        }
                    } else {
                        // J gene is unconditional - use single sampler
                        fast_j_gene_samplers_.resize(1);
                        fast_j_gene_samplers_[0] = CategoricalSampler(probs);
                    }
                } else if (gene_class == D_gene) {
                    d_gene_sequences_ = std::move(seq_mapping);
                    d_realization_indices_ = index_mapping;
                    num_d_genes_ = num_realizations;
                    
                    // D gene depends on V gene AND J gene - build conditional samplers
                    int v_offset = 0, j_offset = 0;
                    if (offset_map_.count(name) > 0) {
                        for (const auto& [parent_event, offset] : offset_map_.at(name)) {
                            if (parent_event->get_class() == V_gene) {
                                v_offset = offset;
                            } else if (parent_event->get_class() == J_gene) {
                                j_offset = offset;
                            }
                        }
                    }
                    
                    if ((v_offset > 0 || j_offset > 0) && num_v_genes_ > 0 && num_j_genes_ > 0) {
                        // Build conditional samplers for each (V, J) combination
                        // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                        fast_d_gene_samplers_.resize(num_v_genes_ * num_j_genes_);
                        for (size_t v_idx = 0; v_idx < num_v_genes_; ++v_idx) {
                            for (size_t j_idx = 0; j_idx < num_j_genes_; ++j_idx) {
                                std::vector<double> cond_probs(num_realizations);
                                // Use position indices directly
                                int cond_base = base_index 
                                    + static_cast<int>(v_idx) * v_offset
                                    + static_cast<int>(j_idx) * j_offset;
                                for (size_t i = 0; i < num_realizations; ++i) {
                                    cond_probs[i] = static_cast<double>(marginal_array[cond_base + index_mapping[i]]);
                                }
                                size_t sampler_idx = v_idx * num_j_genes_ + j_idx;
                                fast_d_gene_samplers_[sampler_idx] = CategoricalSampler(cond_probs);
                            }
                        }
                    } else {
                        // D gene is unconditional - use single sampler
                        fast_d_gene_samplers_.resize(1);
                        fast_d_gene_samplers_[0] = CategoricalSampler(probs);
                    }
                }
                break;
            }
            case Deletion_t: {
                // Build arrays sorted by index for deterministic iteration
                std::vector<std::tuple<int, int, double>> sorted_reals;
                sorted_reals.reserve(realizations.size());
                
                for (const auto& [key, real] : realizations) {
                    double p = static_cast<double>(marginal_array[base_index + real.index]);
                    sorted_reals.emplace_back(real.index, real.value_int, p);
                }
                
                std::sort(sorted_reals.begin(), sorted_reals.end(),
                    [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
                
                std::vector<double> probs;
                std::vector<int> index_mapping;
                std::vector<int> value_mapping;
                
                probs.reserve(realizations.size());
                index_mapping.reserve(realizations.size());
                value_mapping.reserve(realizations.size());
                
                for (const auto& [idx, val, p] : sorted_reals) {
                    probs.push_back(p);
                    index_mapping.push_back(idx);
                    value_mapping.push_back(val);
                }
                
                LinearSampler sampler(probs);
                deletion_samplers_.push_back(std::move(sampler));
                
                // Also build CategoricalSampler for MaxSpeed mode
                CategoricalSampler fast_sampler(probs);
                fast_deletion_samplers_.push_back(std::move(fast_sampler));
                
                if (gene_class == V_gene) {
                    // Build CONDITIONAL samplers for each V gene choice
                    // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                    v_del_values_ = value_mapping;
                    v_del_indices_ = index_mapping;
                    size_t del_size = realizations.size();
                    fast_v_del_samplers_.resize(num_v_genes_);
                    for (size_t v_idx = 0; v_idx < num_v_genes_; ++v_idx) {
                        std::vector<double> cond_probs(del_size);
                        for (size_t i = 0; i < del_size; ++i) {
                            // Use v_idx (position) directly
                            cond_probs[i] = static_cast<double>(
                                marginal_array[base_index + static_cast<int>(v_idx) * static_cast<int>(del_size) + index_mapping[i]]);
                        }
                        fast_v_del_samplers_[v_idx] = CategoricalSampler(cond_probs);
                    }
                } else if (gene_class == D_gene && side == Five_prime) {
                    // Build CONDITIONAL samplers for each D gene choice
                    // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                    d_del_5_values_ = value_mapping;
                    d_del_5_indices_ = index_mapping;
                    num_d_del_5_ = realizations.size();  // Store for D_3' conditional indexing
                    size_t del_size = realizations.size();
                    fast_d_del_5_samplers_.resize(num_d_genes_);
                    for (size_t d_idx = 0; d_idx < num_d_genes_; ++d_idx) {
                        std::vector<double> cond_probs(del_size);
                        for (size_t i = 0; i < del_size; ++i) {
                            // Use d_idx (position) directly
                            cond_probs[i] = static_cast<double>(
                                marginal_array[base_index + static_cast<int>(d_idx) * static_cast<int>(del_size) + index_mapping[i]]);
                        }
                        fast_d_del_5_samplers_[d_idx] = CategoricalSampler(cond_probs);
                    }
                } else if (gene_class == D_gene && side == Three_prime) {
                    // Build CONDITIONAL samplers for each (D gene, D_5' deletion) pair
                    // D_3' deletion depends on BOTH D gene choice AND D_5' deletion
                    // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                    d_del_3_values_ = value_mapping;
                    d_del_3_indices_ = index_mapping;
                    size_t d_del_3_size = realizations.size();
                    
                    // Total combinations: num_d_genes_ * num_d_del_5_
                    fast_d_del_3_samplers_.resize(num_d_genes_ * num_d_del_5_);
                    
                    for (size_t d_idx = 0; d_idx < num_d_genes_; ++d_idx) {
                        for (size_t d5_idx = 0; d5_idx < num_d_del_5_; ++d5_idx) {
                            std::vector<double> cond_probs(d_del_3_size);
                            for (size_t i = 0; i < d_del_3_size; ++i) {
                                // Use position indices directly
                                size_t array_idx = base_index 
                                    + d_idx * (num_d_del_5_ * d_del_3_size)
                                    + d5_idx * d_del_3_size
                                    + static_cast<size_t>(index_mapping[i]);
                                cond_probs[i] = static_cast<double>(marginal_array[array_idx]);
                            }
                            size_t sampler_idx = d_idx * num_d_del_5_ + d5_idx;
                            fast_d_del_3_samplers_[sampler_idx] = CategoricalSampler(cond_probs);
                        }
                    }
                } else if (gene_class == J_gene) {
                    // Build CONDITIONAL samplers for each J gene choice
                    // NOTE: Marginal array uses SEQUENTIAL position (0,1,2...) not realization index
                    j_del_values_ = value_mapping;
                    j_del_indices_ = index_mapping;
                    size_t del_size = realizations.size();
                    fast_j_del_samplers_.resize(num_j_genes_);
                    for (size_t j_idx = 0; j_idx < num_j_genes_; ++j_idx) {
                        std::vector<double> cond_probs(del_size);
                        for (size_t i = 0; i < del_size; ++i) {
                            // Use j_idx (position) directly
                            cond_probs[i] = static_cast<double>(
                                marginal_array[base_index + static_cast<int>(j_idx) * static_cast<int>(del_size) + index_mapping[i]]);
                        }
                        fast_j_del_samplers_[j_idx] = CategoricalSampler(cond_probs);
                    }
                }
                break;
            }
            case Insertion_t: {
                // Build arrays sorted by index for deterministic iteration
                std::vector<std::tuple<int, int, double>> sorted_reals;
                sorted_reals.reserve(realizations.size());
                
                for (const auto& [key, real] : realizations) {
                    double p = static_cast<double>(marginal_array[base_index + real.index]);
                    sorted_reals.emplace_back(real.index, real.value_int, p);
                }
                
                std::sort(sorted_reals.begin(), sorted_reals.end(),
                    [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
                
                std::vector<double> probs;
                std::vector<int> index_mapping;
                std::vector<int> value_mapping;
                
                probs.reserve(realizations.size());
                index_mapping.reserve(realizations.size());
                value_mapping.reserve(realizations.size());
                
                for (const auto& [idx, val, p] : sorted_reals) {
                    probs.push_back(p);
                    index_mapping.push_back(idx);
                    value_mapping.push_back(val);
                }
                
                LinearSampler sampler(probs);
                insertion_length_samplers_.push_back(std::move(sampler));
                
                // Also build CategoricalSampler for MaxSpeed mode
                CategoricalSampler fast_sampler(probs);
                fast_insertion_length_samplers_.push_back(std::move(fast_sampler));
                
                if (gene_class == VD_genes) {
                    vd_ins_values_ = std::move(value_mapping);
                    vd_ins_indices_ = std::move(index_mapping);
                } else if (gene_class == DJ_genes) {
                    dj_ins_values_ = std::move(value_mapping);
                    dj_ins_indices_ = std::move(index_mapping);
                } else if (gene_class == VJ_genes) {
                    vj_ins_values_ = std::move(value_mapping);
                    vj_ins_indices_ = std::move(index_mapping);
                }
                break;
            }
            case Dinuclmarkov_t: {
                // Create dinucleotide Markov sampler
                DinucleotideMarkovSampler dinuc_sampler;
                
                // Build the conditional probability matrix from marginals
                Matrix<double> dinuc_matrix(15, 15);
                for (int prev = 0; prev < 15; ++prev) {
                    double row_sum = 0.0;
                    for (int next = 0; next < 4; ++next) {
                        double prob = marginal_array[base_index + prev * 4 + next];
                        dinuc_matrix(prev, next) = prob;
                        row_sum += prob;
                    }
                    if (row_sum > 0) {
                        for (int next = 0; next < 4; ++next) {
                            dinuc_matrix(prev, next) /= row_sum;
                        }
                    }
                }
                dinuc_sampler.build_from_matrix(dinuc_matrix);
                dinuc_samplers_.push_back(std::move(dinuc_sampler));
                break;
            }
            default:
                break;
        }
    }
    
    // Initialize per-thread model data for ExactMatch mode only
    // MaxSpeed mode uses precomputed samplers and doesn't need this
    if (config_.mode == GenerationMode::ExactMatch) {
        for (auto& ctx : thread_contexts_) {
            ctx.model_queue = ctx.model_parms->get_model_queue();
            ctx.index_map = marginals_.get_index_map(*ctx.model_parms, ctx.model_queue);
            ctx.offset_map = marginals_.get_offsets_map(*ctx.model_parms, ctx.model_queue);
            
            // Call update_event_internal_probas for each event in this thread's model
            std::queue<std::shared_ptr<Rec_Event>> init_queue = ctx.model_queue;
            while (!init_queue.empty()) {
                init_queue.front()->update_event_internal_probas(marginal_array, ctx.index_map);
                init_queue.pop();
            }
        }
    }
}

void FastGenerator::generate_single_sequence(ThreadContext& ctx, size_t seq_idx) {
    ctx.buffer.clear();
    
    // Use the EXACT same sampling logic as GenModel::generate_unique_sequence
    // This ensures perfect reproducibility with conditional dependencies
    
    // Use THREAD-LOCAL model_queue, and COPY index_map (it gets modified during sampling)
    // Each thread has its own model copy to avoid race conditions
    auto model_queue = ctx.model_queue;  // Copy the queue (lightweight - just shared_ptrs to thread's own events)
    auto index_map = ctx.index_map;      // Copy the index_map (gets modified by offset updates)
    
    std::unordered_map<Seq_type, std::string> constructed_sequences;
    std::queue<std::queue<int>> realizations;
    
    // Process events in queue order, exactly like the original
    while (!model_queue.empty()) {
        auto event = model_queue.front();
        model_queue.pop();
        
        // Call the original event's draw_random_realization method
        auto event_realizations = event->draw_random_realization(
            marginals_.marginal_array_smart_p,
            index_map,
            ctx.offset_map,  // Use thread's offset_map
            constructed_sequences,
            ctx.rng
        );
        realizations.push(event_realizations);
    }
    
    // Build final sequence from constructed_sequences
    ctx.buffer.final_sequence = constructed_sequences[V_gene_seq] 
                              + constructed_sequences[VJ_ins_seq]
                              + constructed_sequences[VD_ins_seq] 
                              + constructed_sequences[D_gene_seq] 
                              + constructed_sequences[DJ_ins_seq] 
                              + constructed_sequences[J_gene_seq];
    
    // Convert realizations queue to our buffer format
    // Order in queue matches model_queue order: V, J, D, V_del, D_del5, D_del3, J_del, VD_ins, VD_dinuc, DJ_ins, DJ_dinuc
    auto queue_to_vector = [](std::queue<int>& q) {
        std::vector<int> v;
        while (!q.empty()) {
            v.push_back(q.front());
            q.pop();
        }
        return v;
    };
    
    // V gene
    if (!realizations.empty()) {
        ctx.buffer.v_gene_realization = queue_to_vector(realizations.front());
        realizations.pop();
    }
    // J gene
    if (!realizations.empty()) {
        ctx.buffer.j_gene_realization = queue_to_vector(realizations.front());
        realizations.pop();
    }
    
    if (has_d_gene_) {
        // D gene
        if (!realizations.empty()) {
            ctx.buffer.d_gene_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // V deletion
        if (!realizations.empty()) {
            ctx.buffer.v_del_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // D deletion 5'
        if (!realizations.empty()) {
            ctx.buffer.d_del_5_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // D deletion 3'
        if (!realizations.empty()) {
            ctx.buffer.d_del_3_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // J deletion
        if (!realizations.empty()) {
            ctx.buffer.j_del_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // VD insertion length
        if (!realizations.empty()) {
            ctx.buffer.vd_ins_length_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // VD dinucleotide
        if (!realizations.empty()) {
            ctx.buffer.vd_dinuc_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // DJ insertion length
        if (!realizations.empty()) {
            ctx.buffer.dj_ins_length_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // DJ dinucleotide
        if (!realizations.empty()) {
            ctx.buffer.dj_dinuc_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
    } else {
        // VJ recombination (no D gene)
        // V deletion
        if (!realizations.empty()) {
            ctx.buffer.v_del_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // J deletion
        if (!realizations.empty()) {
            ctx.buffer.j_del_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // VJ insertion length
        if (!realizations.empty()) {
            ctx.buffer.vj_ins_length_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
        // VJ dinucleotide
        if (!realizations.empty()) {
            ctx.buffer.vj_dinuc_realization = queue_to_vector(realizations.front());
            realizations.pop();
        }
    }
    
    // Apply errors if configured
    if (config_.generate_errors) {
        auto err_rate = model_parms_.get_err_rate_p();
        if (err_rate) {
            auto error_queue = err_rate->generate_errors(ctx.buffer.final_sequence, ctx.rng);
            while (!error_queue.empty()) {
                ctx.buffer.error_positions.push_back(error_queue.front());
                error_queue.pop();
            }
        }
    }
}

void FastGenerator::generate_single_sequence_fast(ThreadContext& ctx, size_t seq_idx) {
    ctx.buffer.clear();
    
    // MaxSpeed mode: Use precomputed CONDITIONAL samplers for all events
    // This properly captures the dependency structure in the model
    static const char NT_CHARS[] = {'A', 'C', 'G', 'T'};
    
    // Sample gene choices using CONDITIONAL fast categorical samplers
    
    // V gene - unconditional
    size_t v_gene_idx = fast_v_gene_sampler_.sample(ctx.rng);
    ctx.buffer.v_gene_realization.push_back(v_realization_indices_[v_gene_idx]);
    const std::string& v_seq = v_gene_sequences_[v_gene_idx];
    
    // J gene - CONDITIONAL on V gene
    size_t j_gene_idx;
    if (fast_j_gene_samplers_.size() > 1) {
        // Use conditional sampler based on V gene
        j_gene_idx = fast_j_gene_samplers_[v_gene_idx].sample(ctx.rng);
    } else {
        // Unconditional
        j_gene_idx = fast_j_gene_samplers_[0].sample(ctx.rng);
    }
    ctx.buffer.j_gene_realization.push_back(j_realization_indices_[j_gene_idx]);
    const std::string& j_seq = j_gene_sequences_[j_gene_idx];
    
    // D gene (if present) - CONDITIONAL on V gene AND J gene
    const std::string* d_seq_ptr = nullptr;
    size_t d_gene_idx = 0;
    if (has_d_gene_) {
        if (fast_d_gene_samplers_.size() > 1) {
            // Use conditional sampler based on (V, J) combination
            size_t d_sampler_idx = v_gene_idx * num_j_genes_ + j_gene_idx;
            d_gene_idx = fast_d_gene_samplers_[d_sampler_idx].sample(ctx.rng);
        } else {
            // Unconditional
            d_gene_idx = fast_d_gene_samplers_[0].sample(ctx.rng);
        }
        ctx.buffer.d_gene_realization.push_back(d_realization_indices_[d_gene_idx]);
        d_seq_ptr = &d_gene_sequences_[d_gene_idx];
    }
    
    // Sample deletions using CONDITIONAL fast categorical samplers
    // ALL deletions are conditional on their respective gene choice
    
    // V deletion - CONDITIONAL on V gene choice
    size_t v_del_result = fast_v_del_samplers_[v_gene_idx].sample(ctx.rng);
    int v_del = v_del_values_[v_del_result];
    ctx.buffer.v_del_realization.push_back(v_del_indices_[v_del_result]);
    
    // Pre-reserve final sequence capacity
    ctx.buffer.final_sequence.clear();
    ctx.buffer.final_sequence.reserve(500);  // Typical max sequence length
    
    if (has_d_gene_) {
        // D deletion 5' - CONDITIONAL on D gene choice
        size_t d_del_5_result = fast_d_del_5_samplers_[d_gene_idx].sample(ctx.rng);
        int d_del_5 = d_del_5_values_[d_del_5_result];
        ctx.buffer.d_del_5_realization.push_back(d_del_5_indices_[d_del_5_result]);
        
        // D deletion 3' - CONDITIONAL on D gene choice AND D_5' deletion
        // Sampler index = d_gene_idx * num_d_del_5_ + d_del_5_result
        size_t d3_sampler_idx = d_gene_idx * num_d_del_5_ + d_del_5_result;
        size_t d_del_3_result = fast_d_del_3_samplers_[d3_sampler_idx].sample(ctx.rng);
        int d_del_3 = d_del_3_values_[d_del_3_result];
        ctx.buffer.d_del_3_realization.push_back(d_del_3_indices_[d_del_3_result]);
        
        // J deletion - CONDITIONAL on J gene choice
        size_t j_del_result = fast_j_del_samplers_[j_gene_idx].sample(ctx.rng);
        int j_del = j_del_values_[j_del_result];
        ctx.buffer.j_del_realization.push_back(j_del_indices_[j_del_result]);
        
        // Sample insertions using fast categorical samplers
        size_t ins_sampler_idx = 0;
        
        // VD insertion
        size_t vd_ins_result = fast_insertion_length_samplers_[ins_sampler_idx++].sample(ctx.rng);
        int vd_ins_len = vd_ins_values_[vd_ins_result];
        ctx.buffer.vd_ins_length_realization.push_back(vd_ins_indices_[vd_ins_result]);
        
        // DJ insertion
        size_t dj_ins_result = fast_insertion_length_samplers_[ins_sampler_idx++].sample(ctx.rng);
        int dj_ins_len = dj_ins_values_[dj_ins_result];
        ctx.buffer.dj_ins_length_realization.push_back(dj_ins_indices_[dj_ins_result]);
        
        // Build final sequence directly without intermediate strings
        // V gene (with 3' deletion)
        size_t v_len = v_seq.size();
        size_t v_effective_len = (v_del >= 0 && static_cast<size_t>(v_del) < v_len) 
                                ? v_len - v_del : v_len;
        ctx.buffer.final_sequence.append(v_seq, 0, v_effective_len);
        
        // VD insertion - generate directly into final sequence
        if (vd_ins_len > 0 && !dinuc_samplers_.empty()) {
            size_t prev_nt = v_seq.empty() ? 0 : static_cast<size_t>(
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'A') ? 0 :
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'C') ? 1 :
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'G') ? 2 : 3);
            ctx.buffer.vd_dinuc_realization.resize(vd_ins_len);
            for (int i = 0; i < vd_ins_len; ++i) {
                size_t next_nt = dinuc_samplers_[0].sample_next(prev_nt, ctx.rng);
                ctx.buffer.vd_dinuc_realization[i] = static_cast<int>(next_nt);
                ctx.buffer.final_sequence.push_back(NT_CHARS[next_nt]);
                prev_nt = next_nt;
            }
        }
        
        // D gene (with 5' and 3' deletions)
        if (d_seq_ptr && !d_seq_ptr->empty()) {
            const std::string& d_seq = *d_seq_ptr;
            size_t d_len = d_seq.size();
            size_t start = std::min(static_cast<size_t>(d_del_5), d_len);
            size_t end_trim = std::min(static_cast<size_t>(d_del_3), d_len - start);
            if (start + end_trim < d_len) {
                ctx.buffer.final_sequence.append(d_seq, start, d_len - start - end_trim);
            }
        }
        
        // DJ insertion - generate directly into final sequence (reverse order)
        if (dj_ins_len > 0 && dinuc_samplers_.size() > 1) {
            size_t prev_nt = j_seq.empty() ? 0 : static_cast<size_t>(
                (j_seq[0] == 'A') ? 0 : (j_seq[0] == 'C') ? 1 : (j_seq[0] == 'G') ? 2 : 3);
            ctx.buffer.dj_dinuc_realization.resize(dj_ins_len);
            // Generate in forward order, then reverse
            size_t insert_pos = ctx.buffer.final_sequence.size();
            for (int i = 0; i < dj_ins_len; ++i) {
                size_t next_nt = dinuc_samplers_[1].sample_next(prev_nt, ctx.rng);
                ctx.buffer.dj_dinuc_realization[dj_ins_len - 1 - i] = static_cast<int>(next_nt);
                ctx.buffer.final_sequence.push_back(NT_CHARS[next_nt]);
                prev_nt = next_nt;
            }
            // Reverse the inserted portion
            std::reverse(ctx.buffer.final_sequence.begin() + insert_pos, ctx.buffer.final_sequence.end());
        }
        
        // J gene (with 5' deletion)
        size_t j_start = (j_del >= 0 && static_cast<size_t>(j_del) < j_seq.size()) 
                        ? static_cast<size_t>(j_del) : 0;
        ctx.buffer.final_sequence.append(j_seq, j_start, std::string::npos);
    } else {
        // VJ recombination (no D gene)
        // J deletion - CONDITIONAL on J gene choice
        size_t j_del_result = fast_j_del_samplers_[j_gene_idx].sample(ctx.rng);
        int j_del = j_del_values_[j_del_result];
        ctx.buffer.j_del_realization.push_back(j_del_indices_[j_del_result]);
        
        // VJ insertion
        size_t vj_ins_result = fast_insertion_length_samplers_[0].sample(ctx.rng);
        int vj_ins_len = vj_ins_values_[vj_ins_result];
        ctx.buffer.vj_ins_length_realization.push_back(vj_ins_indices_[vj_ins_result]);
        
        // Build final sequence directly
        // V gene (with 3' deletion)
        size_t v_len = v_seq.size();
        size_t v_effective_len = (v_del >= 0 && static_cast<size_t>(v_del) < v_len)
                                ? v_len - v_del : v_len;
        ctx.buffer.final_sequence.append(v_seq, 0, v_effective_len);
        
        // VJ insertion
        if (vj_ins_len > 0 && !dinuc_samplers_.empty()) {
            size_t prev_nt = v_seq.empty() ? 0 : static_cast<size_t>(
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'A') ? 0 :
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'C') ? 1 :
                (v_seq[v_effective_len > 0 ? v_effective_len - 1 : 0] == 'G') ? 2 : 3);
            ctx.buffer.vj_dinuc_realization.resize(vj_ins_len);
            for (int i = 0; i < vj_ins_len; ++i) {
                size_t next_nt = dinuc_samplers_[0].sample_next(prev_nt, ctx.rng);
                ctx.buffer.vj_dinuc_realization[i] = static_cast<int>(next_nt);
                ctx.buffer.final_sequence.push_back(NT_CHARS[next_nt]);
                prev_nt = next_nt;
            }
        }
        
        // J gene (with 5' deletion)
        size_t j_start = (j_del >= 0 && static_cast<size_t>(j_del) < j_seq.size())
                        ? static_cast<size_t>(j_del) : 0;
        ctx.buffer.final_sequence.append(j_seq, j_start, std::string::npos);
    }
    
    // Apply errors if configured
    if (config_.generate_errors) {
        auto err_rate = model_parms_.get_err_rate_p();
        if (err_rate) {
            auto error_queue = err_rate->generate_errors(ctx.buffer.final_sequence, ctx.rng);
            while (!error_queue.empty()) {
                ctx.buffer.error_positions.push_back(error_queue.front());
                error_queue.pop();
            }
        }
    }
}

std::string FastGenerator::format_realizations(const RealizationBuffer& buffer) const {
    std::ostringstream oss;
    
    // Format must match original model queue order:
    // V_gene, J_gene, D_gene, V_del, D_del_5, D_del_3, J_del, VD_ins, VD_dinuc, DJ_ins, DJ_dinuc
    
    auto format_vec = [&oss](const std::vector<int>& vec) {
        oss << "(";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i > 0) oss << ",";
            oss << vec[i];
        }
        oss << ")";
    };
    
    // Gene choices first: V, J, D
    format_vec(buffer.v_gene_realization);
    oss << ";";
    format_vec(buffer.j_gene_realization);
    oss << ";";
    
    if (has_d_gene_) {
        format_vec(buffer.d_gene_realization);
        oss << ";";
    }
    
    // Deletions: V_del, D_del_5, D_del_3, J_del
    format_vec(buffer.v_del_realization);
    oss << ";";
    
    if (has_d_gene_) {
        format_vec(buffer.d_del_5_realization);
        oss << ";";
        format_vec(buffer.d_del_3_realization);
        oss << ";";
    }
    
    format_vec(buffer.j_del_realization);
    oss << ";";
    
    // Insertions and dinucleotides
    if (has_d_gene_) {
        format_vec(buffer.vd_ins_length_realization);
        oss << ";";
        format_vec(buffer.vd_dinuc_realization);
        oss << ";";
        format_vec(buffer.dj_ins_length_realization);
        oss << ";";
        format_vec(buffer.dj_dinuc_realization);
    } else {
        format_vec(buffer.vj_ins_length_realization);
        oss << ";";
        format_vec(buffer.vj_dinuc_realization);
    }
    
    if (config_.generate_errors) {
        oss << ";";
        format_vec(buffer.error_positions);
    }
    
    return oss.str();
}

void FastGenerator::generate_worker(size_t thread_id, size_t start_idx, size_t count,
                                    BufferedWriter* seq_writer, BufferedWriter* real_writer) {
    auto& ctx = thread_contexts_[thread_id];
    ctx.seq_batch.clear();
    ctx.real_batch.clear();
    ctx.seq_batch.reserve(config_.batch_size);
    if (real_writer) {
        ctx.real_batch.reserve(config_.batch_size);
    }
    
    // Choose generation method based on mode
    const bool use_fast_mode = (config_.mode == GenerationMode::MaxSpeed);
    
    for (size_t i = 0; i < count; ++i) {
        size_t seq_idx = start_idx + i;
        
        // Generate sequence using appropriate method
        if (use_fast_mode) {
            generate_single_sequence_fast(ctx, seq_idx);
        } else {
            generate_single_sequence(ctx, seq_idx);
        }
        
        // Format output
        ctx.seq_batch.push_back(std::to_string(seq_idx) + ";" + ctx.buffer.final_sequence + "\n");
        
        if (real_writer && config_.output_realizations) {
            ctx.real_batch.push_back(std::to_string(seq_idx) + ";" + format_realizations(ctx.buffer) + "\n");
        }
        
        // Write batch when full
        if (ctx.seq_batch.size() >= config_.batch_size) {
            seq_writer->write_batch(ctx.seq_batch);
            ctx.seq_batch.clear();
            
            if (real_writer) {
                real_writer->write_batch(ctx.real_batch);
                ctx.real_batch.clear();
            }
        }
        
        // Update progress
        ++sequences_generated_;
    }
    
    // Write remaining
    if (!ctx.seq_batch.empty()) {
        seq_writer->write_batch(ctx.seq_batch);
    }
    if (real_writer && !ctx.real_batch.empty()) {
        real_writer->write_batch(ctx.real_batch);
    }
}

void FastGenerator::generate(size_t num_sequences,
                             const std::string& sequence_file,
                             const std::string& realization_file) {
    auto start_time = std::chrono::high_resolution_clock::now();
    sequences_generated_ = 0;
    
    // Create writers
    BufferedWriter seq_writer(sequence_file, config_.io_buffer_size);
    std::unique_ptr<BufferedWriter> real_writer;
    
    if (config_.output_realizations && !realization_file.empty()) {
        real_writer = std::make_unique<BufferedWriter>(realization_file, config_.io_buffer_size);
    }
    
    // Write headers
    seq_writer.write("seq_index;nt_sequence\n");
    if (real_writer) {
        // Build header from model events
        std::string header = "seq_index";
        auto model_queue = model_parms_.get_model_queue();
        while (!model_queue.empty()) {
            header += ";" + model_queue.front()->get_name();
            model_queue.pop();
        }
        if (config_.generate_errors) {
            header += ";Errors";
        }
        header += "\n";
        real_writer->write(header);
    }
    
    size_t num_threads = config_.effective_threads();
    
    if (num_threads == 1) {
        // Single-threaded generation
        generate_worker(0, 0, num_sequences, &seq_writer, real_writer.get());
    } else {
        // Multi-threaded generation
        std::vector<std::thread> threads;
        size_t sequences_per_thread = num_sequences / num_threads;
        size_t remaining = num_sequences % num_threads;
        
        size_t current_start = 0;
        for (size_t t = 0; t < num_threads; ++t) {
            size_t count = sequences_per_thread + (t < remaining ? 1 : 0);
            threads.emplace_back(&FastGenerator::generate_worker, this,
                                 t, current_start, count,
                                 &seq_writer, real_writer.get());
            current_start += count;
        }
        
        // Wait for all threads
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    // Flush writers
    seq_writer.flush();
    if (real_writer) {
        real_writer->flush();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    
    stats_.total_sequences = num_sequences;
    stats_.total_time_seconds = elapsed.count();
    stats_.sequences_per_second = num_sequences / elapsed.count();
}

void FastGenerator::generate(size_t num_sequences,
                             std::function<void(size_t idx, const std::string& seq,
                                                const std::vector<std::vector<int>>& realizations)> callback) {
    sequences_generated_ = 0;
    
    auto& ctx = thread_contexts_[0];
    
    for (size_t i = 0; i < num_sequences; ++i) {
        generate_single_sequence(ctx, i);
        
        // Collect realizations
        std::vector<std::vector<int>> realizations;
        realizations.push_back(ctx.buffer.v_gene_realization);
        realizations.push_back(ctx.buffer.v_del_realization);
        if (has_d_gene_) {
            realizations.push_back(ctx.buffer.d_gene_realization);
            realizations.push_back(ctx.buffer.d_del_5_realization);
            realizations.push_back(ctx.buffer.d_del_3_realization);
            realizations.push_back(ctx.buffer.vd_ins_length_realization);
            realizations.push_back(ctx.buffer.vd_dinuc_realization);
        } else {
            realizations.push_back(ctx.buffer.vj_ins_length_realization);
            realizations.push_back(ctx.buffer.vj_dinuc_realization);
        }
        realizations.push_back(ctx.buffer.j_gene_realization);
        realizations.push_back(ctx.buffer.j_del_realization);
        if (has_d_gene_) {
            realizations.push_back(ctx.buffer.dj_ins_length_realization);
            realizations.push_back(ctx.buffer.dj_dinuc_realization);
        }
        if (config_.generate_errors) {
            realizations.push_back(ctx.buffer.error_positions);
        }
        
        callback(i, ctx.buffer.final_sequence, realizations);
        
        ++sequences_generated_;
        if (progress_callback_ && (i % 10000 == 0)) {
            progress_callback_(i, num_sequences);
        }
    }
}

std::vector<std::pair<std::string, std::vector<std::vector<int>>>>
FastGenerator::generate_to_memory(size_t num_sequences) {
    std::vector<std::pair<std::string, std::vector<std::vector<int>>>> results;
    results.reserve(num_sequences);
    
    generate(num_sequences, [&results](size_t, const std::string& seq,
                                        const std::vector<std::vector<int>>& realizations) {
        results.emplace_back(seq, realizations);
    });
    
    return results;
}

}  // namespace igor
