/*
 * FastGenerator.cpp
 *
 *  Created on: Feb 4, 2026
 *      Fast parallel sequence generation
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

#include <igor/Core/FastGenerator.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace igor {
namespace fast {


void FastGenerator::initialize(const Model_Parms& model_parms,
                               const Model_marginals& model_marginals) {
    event_samplers_.clear();

    // Get model queue and maps for conditional sampling
    std::queue<std::shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    std::unordered_map<Rec_Event_name, int> index_map =
        model_marginals.get_index_map(model_parms, model_queue);

    // Get offset map for conditional dependencies
    model_queue = model_parms.get_model_queue();
    auto offset_map = model_marginals.get_offsets_map(model_parms, model_queue);

    // Reset queue
    model_queue = model_parms.get_model_queue();

    // Check for D gene and count genes
    has_d_gene_ = false;
    auto events_map = model_parms.get_events_map();

    // Get gene counts from events
    if (events_map.count(std::make_tuple(GeneChoice_t, V_gene, Undefined_side)) > 0) {
        auto v_event = events_map.at(std::make_tuple(GeneChoice_t, V_gene, Undefined_side));
        num_v_genes_ = v_event->get_realizations_map().size();
    }
    if (events_map.count(std::make_tuple(GeneChoice_t, J_gene, Undefined_side)) > 0) {
        auto j_event = events_map.at(std::make_tuple(GeneChoice_t, J_gene, Undefined_side));
        num_j_genes_ = j_event->get_realizations_map().size();
    }
    if (events_map.count(std::make_tuple(GeneChoice_t, D_gene, Undefined_side)) > 0) {
        has_d_gene_ = true;
        auto d_event = events_map.at(std::make_tuple(GeneChoice_t, D_gene, Undefined_side));
        num_d_genes_ = d_event->get_realizations_map().size();
    }

    // Initialize sampler for each event in order
    while (!model_queue.empty()) {
        auto event = model_queue.front();
        model_queue.pop();

        FastEventSampler sampler;
        sampler.type = event->get_type();
        sampler.gene_class = event->get_class();
        sampler.side = event->get_side();
        sampler.name = event->get_name();

        int base_index = index_map.at(event->get_name());

        // Determine parent dependencies from offset_map
        if (offset_map.count(event->get_name()) > 0) {
            for (const auto& parent_offset : offset_map.at(event->get_name())) {
                Gene_class parent_class = parent_offset.first->get_class();
                if (sampler.parent_gene_class == Undefined_gene) {
                    sampler.parent_gene_class = parent_class;
                } else if (sampler.parent_gene_class_2 == Undefined_gene) {
                    sampler.parent_gene_class_2 = parent_class;
                }
            }
        }

        switch (sampler.type) {
            case GeneChoice_t:
                initialize_gene_choice_sampler(sampler, event,
                    model_marginals.marginal_array_smart_p, base_index);
                break;

            case Insertion_t:
                initialize_insertion_sampler(sampler, event,
                    model_marginals.marginal_array_smart_p, base_index);
                break;

            case Deletion_t:
                initialize_deletion_sampler(sampler, event,
                    model_marginals.marginal_array_smart_p, base_index);
                break;

            case Dinuclmarkov_t:
                initialize_dinucl_sampler(sampler, event,
                    model_marginals.marginal_array_smart_p, base_index);
                break;

            default:
                break;
        }

        sampler.is_initialized = true;
        event_samplers_.push_back(std::move(sampler));
    }

    initialized_ = true;
}


void FastGenerator::initialize_gene_choice_sampler(
    FastEventSampler& sampler,
    const std::shared_ptr<Rec_Event>& event,
    const Marginal_array_p& marginals,
    int base_index
) {
    auto realizations = event->get_realizations_map();
    size_t n = realizations.size();

    // Store gene sequences
    sampler.gene_sequences.resize(n);
    sampler.gene_sequences_int.resize(n);

    for (const auto& pair : realizations) {
        const Event_realization& real = pair.second;
        sampler.gene_sequences[real.index] = real.value_str;
        sampler.gene_sequences_int[real.index] = real.value_str_int;
    }

    // Determine if this is conditional
    Gene_class gc = sampler.gene_class;

    if (gc == V_gene) {
        // V gene is unconditional (root of the model)
        std::vector<double> probs(n);
        for (const auto& pair : realizations) {
            const Event_realization& real = pair.second;
            probs[real.index] = static_cast<double>(marginals[base_index + real.index]);
        }
        sampler.gene_sampler.initialize(probs, true);
        sampler.is_conditional = false;
    }
    else if (gc == J_gene) {
        // J gene is conditional on V: P(J|V)
        // Marginals are stored as: marginals[base_index + v_idx * num_j + j_idx]
        sampler.is_conditional = true;
        sampler.parent_gene_class = V_gene;

        std::vector<std::vector<double>> cond_probs(num_v_genes_);
        for (size_t v = 0; v < num_v_genes_; ++v) {
            cond_probs[v].resize(n);
            double sum = 0.0;
            for (const auto& pair : realizations) {
                const Event_realization& real = pair.second;
                double p = static_cast<double>(marginals[base_index + v * n + real.index]);
                cond_probs[v][real.index] = p;
                sum += p;
            }
            // Normalize
            if (sum > 0) {
                for (size_t j = 0; j < n; ++j) {
                    cond_probs[v][j] /= sum;
                }
            }
        }
        sampler.conditional_gene_sampler.initialize(cond_probs, true);
    }
    else if (gc == D_gene) {
        // D gene is conditional on V and J: P(D|V,J)
        // Marginals: marginals[base_index + v_idx * num_j * num_d + j_idx * num_d + d_idx]
        sampler.is_conditional = true;
        sampler.is_conditional_2d = true;
        sampler.parent_gene_class = V_gene;
        sampler.parent_gene_class_2 = J_gene;

        // Build 2D conditional: samplers[v][j] -> P(D|V=v, J=j)
        sampler.conditional_gene_sampler_2d.resize(num_v_genes_);
        for (size_t v = 0; v < num_v_genes_; ++v) {
            std::vector<std::vector<double>> cond_probs_vj(num_j_genes_);
            for (size_t j = 0; j < num_j_genes_; ++j) {
                cond_probs_vj[j].resize(n);
                double sum = 0.0;
                for (const auto& pair : realizations) {
                    const Event_realization& real = pair.second;
                    double p = static_cast<double>(marginals[base_index + v * num_j_genes_ * n + j * n + real.index]);
                    cond_probs_vj[j][real.index] = p;
                    sum += p;
                }
                // Normalize
                if (sum > 0) {
                    for (size_t d = 0; d < n; ++d) {
                        cond_probs_vj[j][d] /= sum;
                    }
                }
            }
            sampler.conditional_gene_sampler_2d[v].initialize(cond_probs_vj, true);
        }
    }
    else {
        // Fallback: treat as unconditional
        std::vector<double> probs(n);
        for (const auto& pair : realizations) {
            const Event_realization& real = pair.second;
            probs[real.index] = static_cast<double>(marginals[base_index + real.index]);
        }
        sampler.gene_sampler.initialize(probs, true);
        sampler.is_conditional = false;
    }
}


void FastGenerator::initialize_insertion_sampler(
    FastEventSampler& sampler,
    const std::shared_ptr<Rec_Event>& event,
    const Marginal_array_p& marginals,
    int base_index
) {
    auto realizations = event->get_realizations_map();

    // Find max insertion value and build probability vector
    sampler.max_insertions = 0;
    for (const auto& pair : realizations) {
        if (pair.second.value_int > sampler.max_insertions) {
            sampler.max_insertions = pair.second.value_int;
        }
    }

    // Build probability vector indexed by insertion count
    std::vector<double> probs(sampler.max_insertions + 1, 0.0);
    for (const auto& pair : realizations) {
        const Event_realization& real = pair.second;
        if (real.value_int >= 0 && real.value_int <= sampler.max_insertions) {
            probs[real.value_int] = static_cast<double>(marginals[base_index + real.index]);
        }
    }

    sampler.insertion_sampler.initialize(probs, true);
}


void FastGenerator::initialize_deletion_sampler(
    FastEventSampler& sampler,
    const std::shared_ptr<Rec_Event>& event,
    const Marginal_array_p& marginals,
    int base_index
) {
    auto realizations = event->get_realizations_map();
    size_t num_realizations = realizations.size();

    // Find deletion range and build index-to-value mapping
    sampler.min_deletion = INT_MAX;
    sampler.max_deletion = INT_MIN;
    sampler.deletion_idx_to_value.resize(num_realizations);

    for (const auto& pair : realizations) {
        const Event_realization& real = pair.second;
        int val = real.value_int;
        if (val < sampler.min_deletion) sampler.min_deletion = val;
        if (val > sampler.max_deletion) sampler.max_deletion = val;
        if (real.index >= 0 && static_cast<size_t>(real.index) < num_realizations) {
            sampler.deletion_idx_to_value[real.index] = val;
        }
    }

    Gene_class gc = sampler.gene_class;
    Seq_side side = sampler.side;

    // Special case: D 3' deletion is conditional on both D gene and D 5' deletion
    // Dim[num_d, num_d5_del, num_d3_del]
    if (gc == D_gene && side == Three_prime) {
        sampler.is_conditional = true;
        sampler.is_conditional_2d = true;
        sampler.parent_gene_class = D_gene;
        sampler.parent_gene_class_2 = D_gene;  // Second parent is D 5' deletion (same gene class)

        // Build 2D conditional: samplers[d_gene][d5_del] -> P(d3_del|D, d5_del)
        // Assuming D 5' deletion has 21 realizations (same as D 3')
        size_t num_d5_del = num_realizations;  // Same number of realizations

        sampler.deletion_sampler_2d.resize(num_d_genes_);
        for (size_t d = 0; d < num_d_genes_; ++d) {
            std::vector<std::vector<double>> cond_probs_d(num_d5_del);
            for (size_t d5 = 0; d5 < num_d5_del; ++d5) {
                cond_probs_d[d5].resize(num_realizations, 0.0);
                double sum = 0.0;

                for (const auto& pair : realizations) {
                    const Event_realization& real = pair.second;
                    // Marginals: base_index + d * num_d5 * num_d3 + d5 * num_d3 + d3_idx
                    double p = static_cast<double>(marginals[base_index + d * num_d5_del * num_realizations +
                                                             d5 * num_realizations + real.index]);
                    cond_probs_d[d5][real.index] = p;
                    sum += p;
                }

                // Normalize
                if (sum > 0) {
                    for (size_t i = 0; i < num_realizations; ++i) {
                        cond_probs_d[d5][i] /= sum;
                    }
                }
            }
            sampler.deletion_sampler_2d[d].initialize(cond_probs_d, true);
        }
        return;
    }

    // Standard 1D conditional case
    size_t num_conditions = 1;
    if (gc == V_gene) {
        num_conditions = num_v_genes_;
        sampler.parent_gene_class = V_gene;
    } else if (gc == J_gene) {
        num_conditions = num_j_genes_;
        sampler.parent_gene_class = J_gene;
    } else if (gc == D_gene) {
        num_conditions = num_d_genes_;
        sampler.parent_gene_class = D_gene;
    }

    // Build conditional probability matrix P(del|gene)
    std::vector<std::vector<double>> cond_probs(num_conditions);

    for (size_t g = 0; g < num_conditions; ++g) {
        cond_probs[g].resize(num_realizations, 0.0);
        double sum = 0.0;

        for (const auto& pair : realizations) {
            const Event_realization& real = pair.second;
            // Marginals stored as: base_index + gene_idx * num_realizations + real_idx
            double p = static_cast<double>(marginals[base_index + g * num_realizations + real.index]);
            cond_probs[g][real.index] = p;
            sum += p;
        }

        // Normalize
        if (sum > 0) {
            for (size_t i = 0; i < num_realizations; ++i) {
                cond_probs[g][i] /= sum;
            }
        }
    }

    sampler.deletion_sampler.initialize(cond_probs, true);
    sampler.is_conditional = (num_conditions > 1);
}


void FastGenerator::initialize_dinucl_sampler(
    FastEventSampler& sampler,
    const std::shared_ptr<Rec_Event>& /* event */,
    const Marginal_array_p& marginals,
    int base_index
) {
    // The dinucleotide model has a 15x15 or 4x4 transition matrix
    // We need to extract it from the marginals

    // For a 4x4 standard nucleotide matrix:
    // P(next=j | prev=i) is stored at marginals[base_index + i*4 + j]
    // or using IGoR's convention: marginals[base_index + i + 4*j]

    const size_t n = 4;  // Standard nucleotides only for generation
    std::vector<double> dinuc_probs(n * n);

    // Extract the transition probabilities
    // IGoR stores these as: marginals[base_index + prev_nt + n_nucleotides * next_nt]
    for (size_t prev = 0; prev < n; ++prev) {
        double row_sum = 0.0;
        for (size_t next = 0; next < n; ++next) {
            double p = static_cast<double>(marginals[base_index + prev + n * next]);
            dinuc_probs[prev * n + next] = p;
            row_sum += p;
        }
        // Normalize each row
        if (row_sum > 0) {
            for (size_t next = 0; next < n; ++next) {
                dinuc_probs[prev * n + next] /= row_sum;
            }
        }
    }

    sampler.dinuc_sampler.initialize(dinuc_probs.data(), n);
}


int FastGenerator::get_parent_gene_index(
    const FastEventSampler& sampler,
    const SamplingContext& context
) const {
    switch (sampler.parent_gene_class) {
        case V_gene: return context.v_gene_idx;
        case J_gene: return context.j_gene_idx;
        case D_gene: return context.d_gene_idx;
        default: return 0;
    }
}


void FastGenerator::sample_gene_choice(
    const FastEventSampler& sampler,
    std::mt19937_64& rng,
    std::unordered_map<Seq_type, std::string>& sequences,
    std::vector<int>& realization,
    SamplingContext& context
) const {
    size_t gene_idx;
    Gene_class gc = sampler.gene_class;

    if (gc == V_gene) {
        // V gene: unconditional marginal
        gene_idx = sampler.gene_sampler.sample(rng);
        context.v_gene_idx = static_cast<int>(gene_idx);
    }
    else if (gc == J_gene && sampler.is_conditional) {
        // J gene: conditional on V
        gene_idx = sampler.conditional_gene_sampler.sample(
            static_cast<size_t>(context.v_gene_idx), rng);
        context.j_gene_idx = static_cast<int>(gene_idx);
    }
    else if (gc == D_gene && sampler.is_conditional_2d) {
        // D gene: conditional on V and J
        gene_idx = sampler.conditional_gene_sampler_2d[context.v_gene_idx].sample(
            static_cast<size_t>(context.j_gene_idx), rng);
        context.d_gene_idx = static_cast<int>(gene_idx);
    }
    else {
        // Fallback to unconditional
        gene_idx = sampler.gene_sampler.sample(rng);
    }

    realization.push_back(static_cast<int>(gene_idx));

    Seq_type seq_type;
    switch (gc) {
        case V_gene: seq_type = V_gene_seq; break;
        case D_gene: seq_type = D_gene_seq; break;
        case J_gene: seq_type = J_gene_seq; break;
        default: return;
    }

    sequences[seq_type] = sampler.gene_sequences[gene_idx];
}


void FastGenerator::sample_insertion(
    const FastEventSampler& sampler,
    std::mt19937_64& rng,
    std::unordered_map<Seq_type, std::string>& sequences,
    std::vector<int>& realization,
    SamplingContext& /* context */
) const {
    size_t num_ins = sampler.insertion_sampler.sample(rng);
    realization.push_back(static_cast<int>(num_ins));

    // Create placeholder string for insertions (will be filled by dinucleotide model)
    std::string ins_seq(num_ins, 'I');

    Seq_type seq_type;
    switch (sampler.gene_class) {
        case VD_genes: seq_type = VD_ins_seq; break;
        case DJ_genes: seq_type = DJ_ins_seq; break;
        case VJ_genes: seq_type = VJ_ins_seq; break;
        default: return;
    }

    sequences[seq_type] = ins_seq;
}


void FastGenerator::sample_deletion(
    const FastEventSampler& sampler,
    std::mt19937_64& rng,
    std::unordered_map<Seq_type, std::string>& sequences,
    std::vector<int>& realization,
    SamplingContext& context
) const {
    size_t del_idx;

    // Handle 2D conditional (D 3' deletion depends on D gene and D 5' deletion)
    if (sampler.is_conditional_2d && sampler.gene_class == D_gene && sampler.side == Three_prime) {
        int d_gene_idx = context.d_gene_idx;
        int d5_del_idx = context.d_5_del_idx;

        if (d_gene_idx < 0) d_gene_idx = 0;
        if (d5_del_idx < 0) d5_del_idx = 0;

        // Sample from P(d3_del|D=d_gene_idx, d5_del=d5_del_idx)
        del_idx = sampler.deletion_sampler_2d[d_gene_idx].sample(static_cast<size_t>(d5_del_idx), rng);
    } else {
        // Standard 1D conditional case
        int parent_idx = get_parent_gene_index(sampler, context);
        if (parent_idx < 0) parent_idx = 0;

        // Sample deletion realization index conditional on gene
        del_idx = sampler.deletion_sampler.sample(static_cast<size_t>(parent_idx), rng);

        // Track D 5' deletion index for D 3' deletion conditioning
        if (sampler.gene_class == D_gene && sampler.side == Five_prime) {
            context.d_5_del_idx = static_cast<int>(del_idx);
        }
    }

    // Convert realization index to actual deletion value
    int num_del = sampler.deletion_idx_to_value[del_idx];
    realization.push_back(static_cast<int>(del_idx));  // Store the realization index for consistency

    // Apply deletion to the appropriate sequence
    Seq_type seq_type;
    switch (sampler.gene_class) {
        case V_gene: seq_type = V_gene_seq; break;
        case D_gene: seq_type = D_gene_seq; break;
        case J_gene: seq_type = J_gene_seq; break;
        default: return;
    }

    std::string& seq = sequences[seq_type];

    if (num_del >= 0) {
        // Positive deletion: remove bases
        if (sampler.side == Three_prime) {
            // Remove from end
            if (static_cast<size_t>(num_del) < seq.size()) {
                seq.erase(seq.size() - num_del);
            }
        } else if (sampler.side == Five_prime) {
            // Remove from beginning
            if (static_cast<size_t>(num_del) < seq.size()) {
                seq.erase(0, num_del);
            }
        }
    } else {
        // Negative deletion = palindromic insertion
        int num_palindrome = -num_del;
        if (sampler.side == Three_prime) {
            apply_palindrome(seq, num_palindrome, true);
        } else if (sampler.side == Five_prime) {
            apply_palindrome(seq, num_palindrome, false);
        }
    }
}


void FastGenerator::sample_dinucl_markov(
    const FastEventSampler& sampler,
    std::mt19937_64& rng,
    std::unordered_map<Seq_type, std::string>& sequences,
    std::vector<int>& realization,
    SamplingContext& /* context */
) const {
    static const char NT_CHARS[] = {'A', 'C', 'G', 'T'};
    static const int CHAR_TO_INT[] = {
        0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3
    };  // A=0, C=2, G=6, T=19 offset from 'A'

    auto get_last_nt = [&](const std::string& s) -> int {
        if (s.empty()) return 0;  // Default to A
        char c = s.back();
        int idx = c - 'A';
        if (idx >= 0 && idx < 20 && CHAR_TO_INT[idx] >= 0) {
            return CHAR_TO_INT[idx];
        }
        return 0;
    };

    auto get_first_nt = [&](const std::string& s) -> int {
        if (s.empty()) return 0;  // Default to A
        char c = s.front();
        int idx = c - 'A';
        if (idx >= 0 && idx < 20 && CHAR_TO_INT[idx] >= 0) {
            return CHAR_TO_INT[idx];
        }
        return 0;
    };

    // Process based on gene class
    if (sampler.gene_class == VD_genes || sampler.gene_class == VDJ_genes) {
        std::string& ins_seq = sequences[VD_ins_seq];
        if (!ins_seq.empty()) {
            int prev_nt = get_last_nt(sequences[V_gene_seq]);
            for (size_t i = 0; i < ins_seq.size(); ++i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }

    if (sampler.gene_class == DJ_genes || sampler.gene_class == VDJ_genes) {
        std::string& ins_seq = sequences[DJ_ins_seq];
        if (!ins_seq.empty()) {
            // DJ insertions are generated from J side going towards D
            int prev_nt = get_first_nt(sequences[J_gene_seq]);
            // Generate in reverse order
            for (int i = static_cast<int>(ins_seq.size()) - 1; i >= 0; --i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }

    if (sampler.gene_class == VJ_genes) {
        std::string& ins_seq = sequences[VJ_ins_seq];
        if (!ins_seq.empty()) {
            int prev_nt = get_last_nt(sequences[V_gene_seq]);
            for (size_t i = 0; i < ins_seq.size(); ++i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }
}


void FastGenerator::apply_palindrome(std::string& seq, int num_bases, bool is_3prime) {
    auto complement_char = [](char c) -> char {
        switch (c) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default: return c;
        }
    };

    if (num_bases <= 0 || seq.empty()) return;

    std::string palindrome;
    if (is_3prime) {
        // Take last num_bases from seq, reverse complement
        size_t start = (num_bases <= static_cast<int>(seq.size()))
                       ? seq.size() - num_bases : 0;
        size_t len = std::min(static_cast<size_t>(num_bases), seq.size());
        palindrome = seq.substr(start, len);
        std::reverse(palindrome.begin(), palindrome.end());
        for (char& c : palindrome) {
            c = complement_char(c);
        }
        seq += palindrome;
    } else {
        // Take first num_bases from seq, reverse complement, prepend
        size_t len = std::min(static_cast<size_t>(num_bases), seq.size());
        palindrome = seq.substr(0, len);
        std::reverse(palindrome.begin(), palindrome.end());
        for (char& c : palindrome) {
            c = complement_char(c);
        }
        seq = palindrome + seq;
    }
}


std::string FastGenerator::assemble_sequence(
    const std::unordered_map<Seq_type, std::string>& sequences
) const {
    std::string result;
    result.reserve(500);  // Typical sequence length

    // Assemble in order: V + VD_ins + D + DJ_ins + J (for VDJ)
    // or: V + VJ_ins + J (for VJ)

    auto get_seq = [&sequences](Seq_type t) -> const std::string& {
        static const std::string empty;
        auto it = sequences.find(t);
        return (it != sequences.end()) ? it->second : empty;
    };

    result += get_seq(V_gene_seq);

    if (has_d_gene_) {
        result += get_seq(VD_ins_seq);
        result += get_seq(D_gene_seq);
        result += get_seq(DJ_ins_seq);
    } else {
        result += get_seq(VJ_ins_seq);
    }

    result += get_seq(J_gene_seq);

    return result;
}


void FastGenerator::generate_single(std::mt19937_64& rng, GeneratedSequence& result) const {
    result.clear();

    std::unordered_map<Seq_type, std::string> sequences;
    sequences.reserve(6);

    result.realizations.reserve(event_samplers_.size());

    // Create sampling context for conditional dependencies
    SamplingContext context;

    // Process each event in order
    for (const auto& sampler : event_samplers_) {
        std::vector<int> event_realization;
        event_realization.reserve(4);

        switch (sampler.type) {
            case GeneChoice_t:
                sample_gene_choice(sampler, rng, sequences, event_realization, context);
                break;

            case Insertion_t:
                sample_insertion(sampler, rng, sequences, event_realization, context);
                break;

            case Deletion_t:
                sample_deletion(sampler, rng, sequences, event_realization, context);
                break;

            case Dinuclmarkov_t:
                sample_dinucl_markov(sampler, rng, sequences, event_realization, context);
                break;

            default:
                break;
        }

        result.realizations.push_back(std::move(event_realization));
    }

    result.sequence = assemble_sequence(sequences);
}


std::vector<GeneratedSequence> FastGenerator::generate(
    size_t num_sequences,
    const FastGeneratorConfig& config,
    SequenceCallback callback,
    ProgressCallback progress
) {
    if (!initialized_) {
        throw std::runtime_error("FastGenerator not initialized. Call initialize() first.");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Determine thread count
    size_t num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = get_optimal_thread_count();
    }

    // Get base seed
    uint64_t base_seed = config.base_seed;
    if (base_seed == 0) {
        base_seed = draw_random_seed();
    }

    std::vector<GeneratedSequence> results;
    if (!callback) {
        results.resize(num_sequences);
    }

    std::atomic<size_t> completed{0};
    std::atomic<size_t> next_progress_report{0};
    const size_t progress_interval = std::max(size_t(1000), num_sequences / 100);

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::mt19937_64 rng(base_seed + thread_id * 1000000ULL);
        GeneratedSequence local_result;

        #pragma omp for schedule(dynamic, config.batch_size)
        for (size_t i = 0; i < num_sequences; ++i) {
            generate_single(rng, local_result);

            if (callback) {
                callback(i, local_result.sequence, local_result.realizations);
            } else {
                results[i] = local_result;
            }

            size_t done = ++completed;

            // Progress reporting
            if (progress && config.show_progress) {
                size_t expected = next_progress_report.load();
                if (done >= expected) {
                    if (next_progress_report.compare_exchange_strong(expected, expected + progress_interval)) {
                        progress(done, num_sequences);
                    }
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();

    stats_.sequences_generated = num_sequences;
    stats_.total_time_seconds = elapsed;
    stats_.sequences_per_second = num_sequences / elapsed;

    return results;
}


void FastGenerator::generate_to_files(
    size_t num_sequences,
    const std::string& seq_filename,
    const std::string& real_filename,
    const FastGeneratorConfig& config,
    ProgressCallback progress
) {
    if (!initialized_) {
        throw std::runtime_error("FastGenerator not initialized. Call initialize() first.");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Determine thread count
    size_t num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = get_optimal_thread_count();
    }

    // Get base seed
    uint64_t base_seed = config.base_seed;
    if (base_seed == 0) {
        base_seed = draw_random_seed();
    }

    // Phase 1: Generate all sequences in parallel WITHOUT any I/O
    std::vector<GeneratedSequence> results(num_sequences);

    std::atomic<size_t> completed{0};
    std::atomic<size_t> next_progress_report{0};
    const size_t progress_interval = std::max(size_t(1000), num_sequences / 100);

    auto gen_start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::mt19937_64 rng(base_seed + thread_id * 1000000ULL);
        GeneratedSequence local_result;

        #pragma omp for schedule(static)
        for (size_t i = 0; i < num_sequences; ++i) {
            generate_single(rng, local_result);
            results[i] = std::move(local_result);

            size_t done = ++completed;

            // Progress reporting (generation phase: 0-80%)
            if (progress && config.show_progress) {
                size_t expected = next_progress_report.load();
                if (done >= expected) {
                    if (next_progress_report.compare_exchange_strong(expected, expected + progress_interval)) {
                        progress(done * 80 / 100, num_sequences);
                    }
                }
            }
        }
    }

    auto gen_end = std::chrono::high_resolution_clock::now();
    double gen_time = std::chrono::duration<double>(gen_end - gen_start).count();

    // Phase 2: Write to files sequentially with large buffers
    auto io_start = std::chrono::high_resolution_clock::now();

    const size_t write_buffer_size = 4 * 1024 * 1024;  // 4MB buffer

    std::ofstream seq_file(seq_filename, std::ios::binary);
    std::ofstream real_file(real_filename, std::ios::binary);

    if (!seq_file || !real_file) {
        throw std::runtime_error("Failed to open output files");
    }

    // Write headers
    seq_file << "seq_index;nt_sequence\n";

    real_file << "seq_index";
    for (const auto& sampler : event_samplers_) {
        real_file << ";" << sampler.name;
    }
    real_file << "\n";

    // Write data with buffering
    std::string seq_buffer, real_buffer;
    seq_buffer.reserve(write_buffer_size);
    real_buffer.reserve(write_buffer_size);

    size_t total_bytes = 0;

    for (size_t i = 0; i < num_sequences; ++i) {
        const auto& result = results[i];

        // Format sequence line
        seq_buffer += std::to_string(i);
        seq_buffer += ';';
        seq_buffer += result.sequence;
        seq_buffer += '\n';

        // Format realizations line
        real_buffer += std::to_string(i);
        for (const auto& event_real : result.realizations) {
            real_buffer += ";(";
            for (size_t j = 0; j < event_real.size(); ++j) {
                if (j > 0) real_buffer += ',';
                real_buffer += std::to_string(event_real[j]);
            }
            real_buffer += ')';
        }
        real_buffer += '\n';

        // Flush buffers when they're large enough
        if (seq_buffer.size() >= write_buffer_size) {
            seq_file.write(seq_buffer.data(), seq_buffer.size());
            total_bytes += seq_buffer.size();
            seq_buffer.clear();
        }
        if (real_buffer.size() >= write_buffer_size) {
            real_file.write(real_buffer.data(), real_buffer.size());
            total_bytes += real_buffer.size();
            real_buffer.clear();
        }

        // Progress reporting (write phase: 80-100%)
        if (progress && config.show_progress && (i % progress_interval == 0)) {
            progress(num_sequences * 80 / 100 + i * 20 / 100, num_sequences);
        }
    }

    // Final flush
    if (!seq_buffer.empty()) {
        seq_file.write(seq_buffer.data(), seq_buffer.size());
        total_bytes += seq_buffer.size();
    }
    if (!real_buffer.empty()) {
        real_file.write(real_buffer.data(), real_buffer.size());
        total_bytes += real_buffer.size();
    }

    auto io_end = std::chrono::high_resolution_clock::now();
    double io_time = std::chrono::duration<double>(io_end - io_start).count();

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();

    // Print timing breakdown
    std::cerr << "  [Timing] Generation: " << gen_time << "s, I/O: " << io_time << "s, Total: " << elapsed << "s" << std::endl;

    stats_.sequences_generated = num_sequences;
    stats_.total_time_seconds = elapsed;
    stats_.sequences_per_second = num_sequences / elapsed;
    stats_.bytes_written = total_bytes;
}


}  // namespace fast
}  // namespace igor
