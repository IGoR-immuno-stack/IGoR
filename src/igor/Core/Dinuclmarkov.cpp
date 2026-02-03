/*
 * Dinuclmarkov.cpp
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
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/SequenceTypes.h>

using namespace std;

Dinucl_markov::Dinucl_markov() : Dinucl_markov(Undefined_gene, Undefined_side) { }

Dinucl_markov::Dinucl_markov(Gene_class gene) : Dinucl_markov(gene, Undefined_side) { }

Dinucl_markov::Dinucl_markov(Gene_class gene, Seq_side side)
    : Rec_Event(gene, side),
      memory_layer_off_threep(-1),
      memory_layer_off_fivep(-1),
      memory_layer_cs(-1),
      upstream_seq_type(-1),
      downstream_seq_type(-1),
      upstream_exists(false),
      downstream_exists(false)
{
    this->type = Event_type::Dinuclmarkov_t;
    this->update_event_name();
}

Dinucl_markov::Dinucl_markov(Gene_class gene, Seq_side side,
                             map<string, Event_realization> &realizations)
    : Dinucl_markov(gene, side)
{
    this->event_realizations = realizations;
}

Dinucl_markov::~Dinucl_markov()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Dinucl_markov::copy()
{
    shared_ptr<Dinucl_markov> new_dinucl_markov_p =
            shared_ptr<Dinucl_markov>(new Dinucl_markov(this->event_class, this->event_side, this->event_realizations));
    new_dinucl_markov_p->priority = this->priority;
    new_dinucl_markov_p->nickname = this->nickname;
    new_dinucl_markov_p->fixed = this->fixed;

    new_dinucl_markov_p->sequence_type_id = this->sequence_type_id;
    new_dinucl_markov_p->upstream_seq_type = this->upstream_seq_type;
    new_dinucl_markov_p->downstream_seq_type = this->downstream_seq_type;
    new_dinucl_markov_p->upstream_exists = this->upstream_exists;
    new_dinucl_markov_p->downstream_exists = this->downstream_exists;

    new_dinucl_markov_p->name = this->name; // Copy name directly
    new_dinucl_markov_p->set_event_identifier(this->event_index);
    return new_dinucl_markov_p;
}

void Dinucl_markov::iterate(
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

    const Int_Str &ins_str = constructed_sequences.at(type_id, memory_layer_cs);
    int ins_len = ins_str.size();

    int us_3_off = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_off_threep);
    int ds_5_off = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_off_fivep);

    if (ins_len != (ds_5_off - us_3_off - 1))
        return;

    double proba = 1.0;
    int prev_nt = int_sequence.at(us_3_off);

    for (int i = 0; i < ins_len; ++i) {
        int current_nt = int_sequence.at(us_3_off + 1 + i);
        string realization_name = "";
        realization_name += int2nt(prev_nt);
        realization_name += int2nt(current_nt);

        if (this->event_realizations.count(realization_name) == 0) {
            proba = 0;
            break;
        }

        proba *= model_parameters_pointer[base_index + this->event_realizations.at(realization_name).index];
        prev_nt = current_nt;
    }

    if (proba > 0) {
        double new_scenario_proba = scenario_proba * proba;
        double scenario_upper_bound_proba = new_scenario_proba;

        downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                          current_downstream_proba_memory_layers.data());

        if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
            // Update constructed sequence with actual nucleotides
            Int_Str real_ins_str(ins_len);
            for (int i = 0; i < ins_len; ++i) {
                real_ins_str[i] = int_sequence.at(us_3_off + 1 + i);
            }
            constructed_sequences.set_value(type_id, real_ins_str, memory_layer_cs);

            Rec_Event::iterate_wrap_up(
                    new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                    base_index_map, offset_map, next_event_ptr_arr, updated_marginals_pointer,
                    model_parameters_pointer, allowed_realizations, constructed_sequences,
                    seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                    mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
        }
    }
}

void Dinucl_markov::add_to_marginals(long double scenario_proba,
                                    Marginal_array_p &updated_marginals, const Int_Str &int_sequence,
                                    Seq_offsets_map &seq_offsets) const
{
    int us_3_off = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_off_threep);
    int ds_5_off = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_off_fivep);
    int ins_len = ds_5_off - us_3_off - 1;

    int prev_nt = int_sequence.at(us_3_off);
    for (int i = 0; i < ins_len; ++i) {
        int current_nt = int_sequence.at(us_3_off + 1 + i);
        string realization_name = "";
        realization_name += int2nt(prev_nt);
        realization_name += int2nt(current_nt);

        int real_index = this->event_realizations.at(realization_name).index;
        if (viterbi_run) {
            updated_marginals[this->new_index + real_index] = scenario_proba;
        } else {
            updated_marginals[this->new_index + real_index] += scenario_proba;
        }
        prev_nt = current_nt;
    }
}

queue<int> Dinucl_markov::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        std::unordered_map<int, std::string> &constructed_sequences,
        std::mt19937_64 &generator) const
{
    queue<int> realization_queue;
    int type_id = this->sequence_type_id;
    string &ins_seq = constructed_sequences[type_id];
    int ins_len = ins_seq.length();

    // Need to find the nucleotide upstream of the insertion
    // For now, let's assume it's the last nucleotide of the upstream sequence
    // This is a bit simplified for Tandem D
    int prev_nt_int = 0; // Default to A
    
    // In generation, constructed_sequences should have the upstream sequence
    // But we need to know which one is upstream.
    // This part of IGoR's generation is usually handled by knowing the order of events.

    for (int i = 0; i < ins_len; ++i) {
        uniform_real_distribution<double> distribution(0.0, 1.0);
        double rand = distribution(generator);
        double prob_count = 0;
        bool found = false;

        for (int current_nt_int = 0; current_nt_int < 4; ++current_nt_int) {
            string realization_name = "";
            realization_name += int2nt(prev_nt_int);
            realization_name += int2nt(current_nt_int);

            prob_count += model_marginals_p[index_map.at(this->get_name()) + this->event_realizations.at(realization_name).index];
            if (prob_count >= rand) {
                ins_seq[i] = int2nt(current_nt_int);
                realization_queue.push(this->event_realizations.at(realization_name).index);
                prev_nt_int = current_nt_int;
                found = true;
                break;
            }
        }
        if (!found) {
            // Fallback
             ins_seq[i] = 'A';
             prev_nt_int = 0;
        }
    }

    return realization_queue;
}

void Dinucl_markov::write2txt(ofstream &outfile)
{
    outfile << "#DinucMarkov;" << event_class << ";" << event_side << ";" << priority << ";"
            << nickname << endl;
    // DinucMarkov doesn't usually write realizations to txt in this way in IGoR, 
    // it's usually a matrix. But let's keep it consistent.
}

void Dinucl_markov::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map,
        Seq_type_str_p_map &constructed_sequences, Safety_bool_map &safety_set,
        shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    // 1. Determine Sequence Type ID (the junction)
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

    for (const auto &neighbor : neighbors) {
        auto status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map,
                                                    processed_events);
        if (status.exists) {
            if (neighbor.is_upstream) {
                upstream_exists = true;
                upstream_seq_type = neighbor.neighbor_type;
            } else {
                downstream_exists = true;
                downstream_seq_type = neighbor.neighbor_type;
            }
        }
    }

    // 3. Memory layers
    constructed_sequences.request_memory_layer(this->sequence_type_id);
    this->memory_layer_cs = constructed_sequences.get_current_memory_layer(this->sequence_type_id);

    if (upstream_exists) {
        memory_layer_off_threep =
                seq_offsets.get_current_memory_layer(upstream_seq_type, Three_prime);
    }
    if (downstream_exists) {
        memory_layer_off_fivep =
                seq_offsets.get_current_memory_layer(downstream_seq_type, Five_prime);
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map,
                                      downstream_proba_map, constructed_sequences, safety_set,
                                      error_rate_p, mismatches_list, seq_offsets, index_map);
}
void Dinucl_markov::set_nickname(string name)
{
    this->nickname = name;
    auto &registry = get_sequence_type_registry();
    int type_id = registry.try_get_type_id(this->nickname);
    if (type_id >= 0) {
        this->sequence_type_id = type_id;
    }
}

void Dinucl_markov::iterate_initialize_Len_proba(
        int considered_junction,
        std::map<int, double> &length_best_proba_map,
        std::queue<std::shared_ptr<Rec_Event>> &model_queue,
        double &scenario_proba,
        const Marginal_array_p &model_parameters_point,
        Index_map &base_index_map,
        Seq_type_str_p_map &constructed_sequences,
        int &seq_len) const
{
    // TODO: implement for Tandem D
    // For now, just call wrap_up with current parameters
    Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map,
                model_queue, scenario_proba, model_parameters_point,
                base_index_map, constructed_sequences, seq_len);
}
