/*
 * Deletion.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: Quentin Marcou
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 analyze and model immune receptors
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

#include <igor/Core/Aligner.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/SequenceTypes.h>
#include <igor/Model/InferenceEngine.h>

using namespace std;

Deletion::Deletion() : Deletion(Undefined_gene, Undefined_side) { }

Deletion::Deletion(int gene, Seq_side side, std::pair<int, int> bounds) : Deletion(gene, side)
{
    int min_del = bounds.first;
    int max_del = bounds.second;
    for (int i = min_del; i <= max_del; ++i) {
        string del_name = to_string(i);
        Event_realization er(del_name, i, "", Int_Str(), this->size());
        this->event_realizations.insert({ del_name, er });
    }
}

Deletion::Deletion(int gene, Seq_side side)
    : Rec_Event(gene, side),
      memory_layer_off_threep(-1),
      memory_layer_off_fivep(-1),
      memory_layer_safety_upstream(-1),
      memory_layer_safety_downstream(-1),
      memory_layer_off_check1(-1),
      memory_layer_off_check2(-1),
      memory_layer_cs(-1),
      upstream_seq_type(-1),
      downstream_seq_type(-1),
      upstream_exists(false),
      downstream_exists(false),
      upstream_chosen(false),
      downstream_chosen(false),
      upstream_ins_type(-1),
      downstream_ins_type(-1)
{
    this->type = Event_type::Deletion_t;

    // Set sequence_type_id based on event_class for generation
    if (gene == V_gene)
        this->sequence_type_id = SequenceTypeRegistry::V_GENE_SEQ;
    else if (gene == D_gene)
        this->sequence_type_id = SequenceTypeRegistry::D_GENE_SEQ;
    else if (gene == J_gene)
        this->sequence_type_id = SequenceTypeRegistry::J_GENE_SEQ;

    this->update_event_name();
}

Deletion::Deletion(int gene, Seq_side side, unordered_map<string, Event_realization> &realizations)
    : Deletion(gene, side)
{
    this->event_realizations = realizations;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        int val = (*iter).second.value_int;
        if (val > this->len_max) {
            this->len_max = val;
        } else if (val < this->len_min) {
            this->len_min = val;
        }
    }
}

Deletion::~Deletion()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Deletion::copy()
{
    shared_ptr<Deletion> new_deletion_p =
            shared_ptr<Deletion>(new Deletion(this->event_class, this->event_side, this->event_realizations));
    new_deletion_p->priority = this->priority;
    new_deletion_p->nickname = this->nickname;
    new_deletion_p->fixed = this->fixed;

    // Copy tandem D support fields
    new_deletion_p->sequence_type_id = this->sequence_type_id;
    new_deletion_p->upstream_seq_type = this->upstream_seq_type;
    new_deletion_p->downstream_seq_type = this->downstream_seq_type;
    new_deletion_p->upstream_exists = this->upstream_exists;
    new_deletion_p->downstream_exists = this->downstream_exists;
    new_deletion_p->upstream_chosen = this->upstream_chosen;
    new_deletion_p->downstream_chosen = this->downstream_chosen;
    new_deletion_p->upstream_ins_type = this->upstream_ins_type;
    new_deletion_p->downstream_ins_type = this->downstream_ins_type;

    new_deletion_p->name = this->name; // Copy name directly
    new_deletion_p->set_event_identifier(this->event_index);
    return new_deletion_p;
}

void Deletion::iterate(double &scenario_proba, Downstream_scenario_proba_bound_map &downstream_proba_map,
                       const string &sequence, const Int_Str &int_sequence, Index_map &base_index_map,
                       const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
                       shared_ptr<Next_event_ptr> &next_event_ptr_arr, Marginal_array_p &updated_marginals_pointer,
                       const Marginal_array_p &model_parameters_pointer,
                       const unordered_map<int, vector<Alignment_data>> &allowed_realizations,
                       Seq_type_str_p_map &constructed_sequences, Seq_offsets_map &seq_offsets,
                       shared_ptr<Error_rate> &error_rate_p, map<size_t, shared_ptr<Counter>> &counters_list,
                       const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
                       Safety_bool_map &safety_set, Mismatch_vectors_map &mismatches_lists,
                       double &seq_max_prob_scenario, double &proba_threshold_factor)
{
    int base_index = base_index_map.at(this->event_index);
    int type_id = this->sequence_type_id;

    const Int_Str &previous_str = constructed_sequences.at(type_id, memory_layer_cs - 1);

    int prev_5_off = seq_offsets.at(type_id, Five_prime, memory_layer_off_fivep);
    int prev_3_off = seq_offsets.at(type_id, Three_prime, memory_layer_off_threep);

    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        int del_len = (*iter).second.value_int;

        if (del_len > (int)previous_str.size())
            continue;

        Int_Str current_str;
        int current_5_off = prev_5_off;
        int current_3_off = prev_3_off;

        if (this->event_side == Three_prime) {
            current_str.assign(previous_str.begin(), previous_str.end() - del_len);
            current_3_off = prev_3_off - del_len;
        } else {
            current_str.assign(previous_str.begin() + del_len, previous_str.end());
            current_5_off = prev_5_off + del_len;
        }

        // Safety Checks
        bool is_safe = true;
        if (event_side == Three_prime && downstream_exists && downstream_chosen) {
            int ds_5_off = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_off_check2);
            if (current_3_off >= ds_5_off) {
                is_safe = false;
            } else {
                safety_set.set_value(memory_layer_safety_downstream, true, memory_layer_safety_downstream);
                // Junction length check
                int junction_len = ds_5_off - current_3_off - 1;
                if (junction_length_best_proba_maps.count(downstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(downstream_ins_type).count(junction_len) == 0) {
                    is_safe = false;
                } else {
                    downstream_proba_map.set_value(
                            downstream_ins_type,
                            junction_length_best_proba_maps.at(downstream_ins_type).at(junction_len),
                            memory_layer_proba_map_junction);
                }
            }
        } else if (event_side == Five_prime && upstream_exists && upstream_chosen) {
            int us_3_off = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_off_check1);
            if (current_5_off <= us_3_off) {
                is_safe = false;
            } else {
                safety_set.set_value(memory_layer_safety_upstream, true, memory_layer_safety_upstream);
                // Junction length check
                int junction_len = current_5_off - us_3_off - 1;
                if (junction_length_best_proba_maps.count(upstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(upstream_ins_type).count(junction_len) == 0) {
                    is_safe = false;
                } else {
                    downstream_proba_map.set_value(
                            upstream_ins_type, junction_length_best_proba_maps.at(upstream_ins_type).at(junction_len),
                            memory_layer_proba_map_junction_upstream);
                }
            }
        }

        if (!is_safe)
            continue;

        double proba_contribution =
                Rec_Event::iterate_common((*iter).second.index, base_index, base_index_map, model_parameters_pointer);
        double new_scenario_proba = scenario_proba * proba_contribution;

        constructed_sequences.set_value(type_id, current_str, memory_layer_cs);
        seq_offsets.set_value(type_id, Five_prime, current_5_off, memory_layer_off_fivep);
        seq_offsets.set_value(type_id, Three_prime, current_3_off, memory_layer_off_threep);

        double scenario_upper_bound_proba = new_scenario_proba;
        downstream_proba_map.multiply_all(scenario_upper_bound_proba, current_downstream_proba_memory_layers.data());

        if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
            Rec_Event::iterate_wrap_up(new_scenario_proba, downstream_proba_map, sequence, int_sequence, base_index_map,
                                       offset_map, next_event_ptr_arr, updated_marginals_pointer,
                                       model_parameters_pointer, allowed_realizations, constructed_sequences,
                                       seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                                       mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
        }
    }
}

queue<int> Deletion::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        std::unordered_map<int, std::string> &constructed_sequences, std::mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;

    int type_id = this->sequence_type_id;

    // Check if type_id is valid and exists in constructed_sequences
    if (type_id < 0) {
        return realization_queue;
    }
    if (constructed_sequences.count(type_id) == 0) {
        return realization_queue;
    }

    string &seq = constructed_sequences[type_id];

    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        // Check if event name is in index_map
        if (index_map.count(this->get_name()) == 0) {
            break;
        }
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            int del_len = (*iter).second.value_int;
            // Take absolute value - deletion length should be positive
            if (del_len < 0) del_len = -del_len;

            if (del_len > (int)seq.length()) {
                del_len = seq.length();
            }

            if (this->event_side == Three_prime) {
				seq.resize(seq.length() - del_len);
      } else {
                seq = seq.substr(del_len);
            }

            realization_queue.push((*iter).second.index);
            if (offset_map.count(this->get_name()) != 0) {
                for (vector<pair<shared_ptr<const Rec_Event>, int>>::const_iterator jiter =
                             offset_map.at(this->get_name()).begin();
                     jiter != offset_map.at(this->get_name()).end(); ++jiter) {
                    index_map.at((*jiter).first->get_name()) += (*iter).second.index * (*jiter).second;
                }
            }
            break;
        }
    }
    return realization_queue;
}

// Phase 2: InferenceEngine-based sampling for Deletion
queue<int> Deletion::draw_random_realization_with_engine(
    const igor::model::InferenceEngine<long double> &engine,
    std::unordered_map<int, std::string> &constructed_sequences,
    std::unordered_map<std::string, std::size_t> &sampled_indices,
    const Model_Parms &model_parms,
    std::mt19937_64 &generator) const
{
    queue<int> realization_queue;

    try {
        // Check if the sequence type exists
        int type_id = this->sequence_type_id;
        if (type_id < 0 || constructed_sequences.count(type_id) == 0) {
            return realization_queue;  // No upstream sequence to delete from
        }

        // Get the handler for this event
        auto& handler = engine.handler(this->get_nickname());

        // Resolve parent indices from sampled_indices
        std::vector<std::size_t> parent_indices;
        auto parents = model_parms.get_parents(this->get_name());
        for (const auto& parent : parents) {
            auto it = sampled_indices.find(parent->get_nickname());
            if (it == sampled_indices.end()) {
                throw std::runtime_error(
                    "Parent event '" + parent->get_nickname() + "' not yet sampled for " + this->get_nickname());
            }
            parent_indices.push_back(it->second);
        }

        // Sample a realization index
        std::size_t realization_idx = handler.sample(generator, parent_indices);

        // Record in sampled_indices for downstream events
        sampled_indices[this->get_nickname()] = realization_idx;

        // Find the realization with this index and apply deletion
        string &seq = constructed_sequences[type_id];
        for (const auto& [name, real] : this->event_realizations) {
            if (real.index == static_cast<int>(realization_idx)) {
                int del_len = real.value_int;
                // Take absolute value - deletion length should be positive
                if (del_len < 0) del_len = -del_len;

                if (del_len > (int)seq.length()) {
                    del_len = seq.length();
                }

                // Apply deletion from Five_prime or Three_prime based on event_side
                if (this->event_side == Three_prime) {
                    seq.resize(seq.length() - del_len);
                } else {
                    seq = seq.substr(del_len);
                }

                realization_queue.push(real.index);
                break;
            }
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "Failed to sample from InferenceEngine for " + this->get_name() + ": " + e.what());
    }

    return realization_queue;
}

void Deletion::write2txt(ofstream &outfile)
{
    outfile << "#Deletion;" << (Gene_class)event_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.name << ";" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

std::vector<std::size_t> Deletion::inherent_shape() const {
    return { event_realizations.size() };
}

void Deletion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    // 1. Determine Sequence Type ID
    if (this->sequence_type_id == -1) {
        auto &registry = get_sequence_type_registry();
        if (!this->nickname.empty()) {
            int type_id = registry.try_get_type_id(this->nickname + "_seq");
            if (type_id >= 0) {
                this->sequence_type_id = type_id;
            }
        }
        if (this->sequence_type_id == -1) {
            if (this->event_class == V_gene)
                this->sequence_type_id = SequenceTypeRegistry::V_GENE_SEQ;
            else if (this->event_class == D_gene)
                this->sequence_type_id = SequenceTypeRegistry::D_GENE_SEQ;
            else if (this->event_class == J_gene)
                this->sequence_type_id = SequenceTypeRegistry::J_GENE_SEQ;
        }
    }

    // 2. Neighbors
    auto &registry = get_sequence_type_registry();
    auto upstream_neighbors = registry.get_upstream_neighbors((SequenceTypeRegistry::TypeId)this->sequence_type_id);
    auto downstream_neighbors = registry.get_downstream_neighbors((SequenceTypeRegistry::TypeId)this->sequence_type_id);

    // Reset
    upstream_seq_type = -1;
    downstream_seq_type = -1;
    upstream_exists = false;
    downstream_exists = false;
    upstream_chosen = false;
    downstream_chosen = false;
    upstream_ins_type = -1;
    downstream_ins_type = -1;

    for (const auto &neighbor : upstream_neighbors) {
        auto status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map, processed_events);
        if (status.exists) {
            upstream_exists = true;
            upstream_seq_type = neighbor.neighbor_type;
            upstream_ins_type = neighbor.junction_type;
            if (status.chosen)
                upstream_chosen = true;
            break;
        }
    }

    for (const auto &neighbor : downstream_neighbors) {
        auto status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map, processed_events);
        if (status.exists) {
            downstream_exists = true;
            downstream_seq_type = neighbor.neighbor_type;
            downstream_ins_type = neighbor.junction_type;
            if (status.chosen)
                downstream_chosen = true;
            break;
        }
    }

    // 3. Request Memory Layers
    seq_offsets.request_memory_layer(this->sequence_type_id, Three_prime);
    this->memory_layer_off_threep = seq_offsets.get_current_memory_layer(this->sequence_type_id, Three_prime);

    seq_offsets.request_memory_layer(this->sequence_type_id, Five_prime);
    this->memory_layer_off_fivep = seq_offsets.get_current_memory_layer(this->sequence_type_id, Five_prime);

    constructed_sequences.request_memory_layer(this->sequence_type_id);
    this->memory_layer_cs = constructed_sequences.get_current_memory_layer(this->sequence_type_id);

    if (upstream_exists) {
        memory_layer_off_check1 = seq_offsets.get_current_memory_layer(upstream_seq_type, Three_prime);
        safety_set.request_memory_layer(upstream_ins_type);
        memory_layer_safety_upstream = safety_set.get_current_memory_layer(upstream_ins_type);

        downstream_proba_map.request_memory_layer(upstream_ins_type);
        memory_layer_proba_map_junction_upstream = downstream_proba_map.get_current_memory_layer(upstream_ins_type);
    }

    if (downstream_exists) {
        memory_layer_off_check2 = seq_offsets.get_current_memory_layer(downstream_seq_type, Five_prime);
        safety_set.request_memory_layer(downstream_ins_type);
        memory_layer_safety_downstream = safety_set.get_current_memory_layer(downstream_ins_type);

        downstream_proba_map.request_memory_layer(downstream_ins_type);
        memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(downstream_ins_type);
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list, seq_offsets,
                                      index_map);
}

void Deletion::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                          const Marginal_array_p &model_parameters_point, Index_map &base_index_map)
{
    Seq_type_str_p_map constructed_sequences(SequenceTypeRegistry::get_instance().size());
    junction_length_best_proba_maps.clear();

    if (downstream_exists) {
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(
                downstream_ins_type, junction_length_best_proba_maps[downstream_ins_type], model_queue, init_proba,
                model_parameters_point, base_index_map, constructed_sequences);
    }
    if (upstream_exists) {
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(
                upstream_ins_type, junction_length_best_proba_maps[upstream_ins_type], model_queue, init_proba,
                model_parameters_point, base_index_map, constructed_sequences);
    }
}
void Deletion::set_nickname(string name)
{
    this->nickname = name;
    auto &registry = get_sequence_type_registry();

    // Try to get type_id directly first
    int type_id = registry.try_get_type_id(this->nickname);

    if (type_id >= 0) {
        // Type already registered, use it
        this->sequence_type_id = type_id;
    } else {
        // Check for tandem D junction deletion patterns
        // Junction deletions have names like V_3_del (standard) or VD1_3_del (tandem D)
        if ((this->event_class == VD_genes || this->event_class == DJ_genes)
            && (this->nickname.find("VD") == 0 || this->nickname.find("DJ") == 0)) {
            // Register as a junction-specific deletion
            this->sequence_type_id = registry.register_junction_type(this->nickname);
        } else {
            // Use event_class as fallback for standard deletions
            this->sequence_type_id = this->event_class;
        }
    }
}

bool Deletion::has_effect_on(int) const
{
    return false; // deletion events affect something?
}

void Deletion::add_to_marginals(long double, Marginal_array_p &) const
{
    // TODO: implement this if needed
}

void Deletion::iterate_initialize_Len_proba(int considered_junction, std::map<int, double> &length_best_proba_map,
                                            std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
                                            const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                            Seq_type_str_p_map &constructed_sequences, int &seq_len) const
{
    // Deletion events don't directly affect junction length probabilities
    // They influence safety constraints which are handled elsewhere
    Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map, model_queue,
                                                    scenario_proba, model_parameters_point, base_index_map,
                                                    constructed_sequences, seq_len);
}
