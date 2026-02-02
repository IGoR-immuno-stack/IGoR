/*
 * Insertion.cpp
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
#include <igor/Core/EventUtils.h>
#include <igor/Core/Insertion.h>

using namespace std;

Insertion::Insertion() : Insertion(Undefined_gene, Undefined_side) { }

Insertion::Insertion(Gene_class gene, std::pair<int, int> bounds)
    : Insertion(gene, Undefined_side)
{
    int min_ins = bounds.first;
    int max_ins = bounds.second;
    for (int i = min_ins; i <= max_ins; ++i) {
        string ins_name = to_string(i);
        Event_realization er(ins_name, i, "", Int_Str(), this->size());
        this->event_realizations.insert({ins_name, er});
    }
}

Insertion::Insertion(Gene_class gene) : Insertion(gene, Undefined_side) { }

Insertion::Insertion(Gene_class gene, unordered_map<string, Event_realization> &realizations)
    : Insertion(gene, Undefined_side, realizations)
{
}

Insertion::Insertion(Gene_class gene, Seq_side side)
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
    this->type = Event_type::Insertion_t;
    this->update_event_name();
}

Insertion::Insertion(Gene_class gene, Seq_side side,
                     unordered_map<string, Event_realization> &realizations)
    : Insertion(gene, side)
{
    this->event_realizations = realizations;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        int val = (*iter).second.value_int;
        if (val > this->len_max) {
            this->len_max = val;
        } else if (val < this->len_min) {
            this->len_min = val;
        }
    }
}

Insertion::~Insertion()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Insertion::copy()
{
    shared_ptr<Insertion> new_insertion_p =
            shared_ptr<Insertion>(new Insertion(this->event_class, this->event_side, this->event_realizations));
    new_insertion_p->priority = this->priority;
    new_insertion_p->nickname = this->nickname;
    new_insertion_p->fixed = this->fixed;

    // Tandem D
    new_insertion_p->sequence_type_id = this->sequence_type_id;
    new_insertion_p->upstream_seq_type = this->upstream_seq_type;
    new_insertion_p->downstream_seq_type = this->downstream_seq_type;
    new_insertion_p->upstream_exists = this->upstream_exists;
    new_insertion_p->downstream_exists = this->downstream_exists;
    new_insertion_p->upstream_chosen = this->upstream_chosen;
    new_insertion_p->downstream_chosen = this->downstream_chosen;
    new_insertion_p->upstream_ins_type = this->upstream_ins_type;
    new_insertion_p->downstream_ins_type = this->downstream_ins_type;

    new_insertion_p->name = this->name; // Copy name directly
    new_insertion_p->set_event_identifier(this->event_index);
    return new_insertion_p;
}

void Insertion::iterate(
        double &scenario_proba, Downstream_scenario_proba_bound_map &downstream_proba_map,
        const string &sequence, const Int_Str &int_sequence, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        shared_ptr<Next_event_ptr> &next_event_ptr_arr, Marginal_array_p &updated_marginals_pointer,
        const Marginal_array_p &model_parameters_pointer,
        const unordered_map<Gene_class, vector<Alignment_data>> &allowed_realizations,
        Seq_type_str_p_map &constructed_sequences, Seq_offsets_map &seq_offsets,
        shared_ptr<Error_rate> &error_rate_p, map<size_t, shared_ptr<Counter>> &counters_list,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        Safety_bool_map &safety_set, Mismatch_vectors_map &mismatches_lists,
        double &seq_max_prob_scenario, double &proba_threshold_factor)
{
    int base_index = base_index_map.at(this->event_index);
    int type_id = this->sequence_type_id;

    if (!upstream_exists || !downstream_exists || !upstream_chosen || !downstream_chosen)
        return;

    int us_3_off = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_off_check1);
    int ds_5_off = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_off_check2);

    int ins_len = ds_5_off - us_3_off - 1;

    if (ins_len < 0)
        return;

    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second .value_int!= ins_len)
            continue;

        double proba_contribution =
                Rec_Event::iterate_common((*iter).second.index, base_index, base_index_map,
                                          model_parameters_pointer);
        double new_scenario_proba = scenario_proba * proba_contribution;

        // Fill insertions with 'N' (14) if they haven't been inferred by Dinuclmarkov yet
        Int_Str ins_str(ins_len, 14);
        constructed_sequences.set_value(type_id, ins_str, memory_layer_cs);

        double scenario_upper_bound_proba = new_scenario_proba;
        downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                          current_downstream_proba_memory_layers.data());

        if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
            Rec_Event::iterate_wrap_up(
                    new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                    base_index_map, offset_map, next_event_ptr_arr, updated_marginals_pointer,
                    model_parameters_pointer, allowed_realizations, constructed_sequences,
                    seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                    mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
        }
    }
}

queue<int> Insertion::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        std::unordered_map<int, std::string> &constructed_sequences,
        std::mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;

    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            int ins_len = (*iter).second.value_int;
            string ins_seq(ins_len, 'N'); // Placeholder
            constructed_sequences[this->sequence_type_id] = ins_seq;

            realization_queue.push((*iter).second.index);
            if (offset_map.count(this->get_name()) != 0) {
                for (vector<pair<shared_ptr<const Rec_Event>, int>>::const_iterator jiter =
                             offset_map.at(this->get_name()).begin();
                     jiter != offset_map.at(this->get_name()).end(); ++jiter) {
                    index_map.at((*jiter).first->get_name()) +=
                            (*iter).second.index * (*jiter).second;
                }
            }
            break;
        }
    }
    return realization_queue;
}

void Insertion::write2txt(ofstream &outfile)
{
    outfile << "#Insertion;" << event_class << ";" << event_side << ";" << priority << ";"
            << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.name << ";" << (*iter).second .value_int<< ";"
                << (*iter).second.index << endl;
    }
}

void Insertion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map,
        Seq_type_str_p_map &constructed_sequences, Safety_bool_map &safety_set,
        shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    // 1. Determine Sequence Type ID
    if (this->sequence_type_id == -1) {
        auto &registry = get_sequence_type_registry();
        if (!this->nickname.empty()) {
            int type_id = registry.try_get_type_id(this->nickname);
            if (type_id >= 0)
                this->sequence_type_id = type_id;
        }
        if (this->sequence_type_id == -1) {
            if (this->event_class == VD_genes)
                this->sequence_type_id = SequenceTypeRegistry::VD_INS_SEQ;
            else if (this->event_class == DJ_genes)
                this->sequence_type_id = SequenceTypeRegistry::DJ_INS_SEQ;
            else if (this->event_class == VJ_genes)
                this->sequence_type_id = SequenceTypeRegistry::VJ_INS_SEQ;
        }
    }

    // 2. Neighbors
    auto &registry = get_sequence_type_registry();
    auto neighbors = registry.get_neighbors_for_junction(
            (SequenceTypeRegistry::TypeId)this->sequence_type_id);

    upstream_exists = false;
    downstream_exists = false;
    upstream_chosen = false;
    downstream_chosen = false;

    for (const auto &neighbor : neighbors) {
        auto status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map,
                                                    processed_events);
        if (status.exists) {
            if (neighbor.is_upstream) {
                upstream_exists = true;
                upstream_seq_type = neighbor.neighbor_type;
                if (status.chosen)
                    upstream_chosen = true;
            } else {
                downstream_exists = true;
                downstream_seq_type = neighbor.neighbor_type;
                if (status.chosen)
                    downstream_chosen = true;
            }
        }
    }

    // 3. Request Memory Layers
    constructed_sequences.request_memory_layer(this->sequence_type_id);
    this->memory_layer_cs = constructed_sequences.get_current_memory_layer(this->sequence_type_id);

    if (upstream_exists) {
        memory_layer_off_check1 =
                seq_offsets.get_current_memory_layer(upstream_seq_type, Three_prime);
    }
    if (downstream_exists) {
        memory_layer_off_check2 =
                seq_offsets.get_current_memory_layer(downstream_seq_type, Five_prime);
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map,
                                      downstream_proba_map, constructed_sequences, safety_set,
                                      error_rate_p, mismatches_list, seq_offsets, index_map);
}

void Insertion::iterate_initialize_Len_proba(
        int considered_junction, std::map<int, double> &length_best_proba_map,
        std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
        const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
        Seq_type_str_p_map &constructed_sequences, int &seq_len /*=0*/) const
{
    if (considered_junction == this->sequence_type_id) {
        int base_index = base_index_map.at(this->event_index, 0);
        for (unordered_map<string, Event_realization>::const_iterator iter =
                     this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()]
                    > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index
                                                            + i * this->size()];
                }
            }
            Rec_Event::iterate_initialize_Len_proba_wrap_up(
                    considered_junction, length_best_proba_map, model_queue,
                    scenario_proba * real_max_proba, model_parameters_point, base_index_map,
                    constructed_sequences, seq_len + (*iter).second.value_int);
        }
    } else {
        Rec_Event::iterate_initialize_Len_proba_wrap_up(
                considered_junction, length_best_proba_map, model_queue, scenario_proba,
                model_parameters_point, base_index_map, constructed_sequences, seq_len);
    }
}
void Insertion::set_nickname(string name)
{
    this->nickname = name;
    auto &registry = get_sequence_type_registry();
    int type_id = registry.try_get_type_id(this->nickname);
    
    if (type_id >= 0) {
        // Type already registered, use it
        this->sequence_type_id = type_id;
    } else {
        // Check for tandem D junction patterns: VD1_ins, V1D2_ins, D2J_ins
        // Also handle VD_genes and DJ_genes standard cases
        if (this->nickname.find("VD") == 0 || this->nickname.find("DJ") == 0) {
            // Register as a tandem D junction
            this->sequence_type_id = registry.register_junction_type(this->nickname);
        } else if (this->event_class == VD_genes) {
            this->sequence_type_id = SequenceTypeRegistry::VD_INS_SEQ;
        } else if (this->event_class == DJ_genes) {
            this->sequence_type_id = SequenceTypeRegistry::DJ_INS_SEQ;
        } else if (this->event_class == VJ_genes) {
            this->sequence_type_id = SequenceTypeRegistry::VJ_INS_SEQ;
        } else {
            this->sequence_type_id = -1;
        }
    }
}

bool Insertion::has_effect_on(int) const
{
    return false; // to be implemented if needed
}

void Insertion::add_to_marginals(long double, Marginal_array_p &) const
{
    // TODO: implement
}

void Insertion::set_crude_upper_bound_proba(size_t, size_t, Marginal_array_p&)
{
    // TODO: implement
}

void Insertion::update_event_internal_probas(const Marginal_array_p &,
                                         const std::unordered_map<std::string, int> &)
{
    // TODO: implement
}

void Insertion::initialize_crude_scenario_proba_bound(
        double &,
        std::forward_list<double *> &,
        const std::unordered_map<std::tuple<Event_type, int, Seq_side>,
                                std::shared_ptr<Rec_Event>> &)
{
    // TODO: implement
}
