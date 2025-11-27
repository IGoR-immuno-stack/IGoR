/*
 * Deletion.cpp
 *
 *  Created on: Dec 9, 2014
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

#include <igor/Core/Deletion.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/SequenceTypes.h>

using namespace std;
Deletion::Deletion() : Deletion(Undefined_gene, Undefined_side)
{
    this->type = Event_type::Deletion_t;
    this->update_event_name();
}

Deletion::Deletion(Gene_class gene, Seq_side del_side, pair<int, int> del_range)
    : Deletion(gene, del_side)
{
    // TODO prevent undefined_side entry(throw exception)

    int min_del = min(del_range.first, del_range.second);
    int max_del = max(del_range.first, del_range.second);

    this->type = Event_type::Deletion_t;
    this->len_max = -min_del;
    this->len_min = -max_del;

    for (int i = min_del; i != max_del + 1; ++i) {
        this->add_realization(i);
    }
    this->update_event_name();
}

Deletion::Deletion(Gene_class gene, Seq_side side)
    : Rec_Event(gene, side),
      new_scenario_proba(-1),
      d_3_max_del(INT16_MAX),
      v_3_new_offset(INT16_MAX),
      memory_layer_mismatches(-1),
      memory_layer_offset_del(-1),
      v_3_min_del(INT16_MAX),
      memory_layer_offset_check2(-1),
      d_5_min_offset(INT16_MAX),
      new_index(-1),
      j_5_max_offset(INT16_MAX),
      d_3_max_offset(INT16_MAX),
      dj_check(true),
      j_5_offset(INT16_MAX),
      end_reached(false),
      d_chosen(false),
      memory_layer_safety_1(-1),
      new_tmp_err_w_proba(-1),
      d_5_max_del(INT16_MAX),
      v_3_min_offset(INT16_MAX),
      j_chosen(false),
      memory_layer_safety_2(-1),
      err_rate_upper_bound(-1),
      v_chosen(false),
      vd_check(true),
      j_5_min_offset(INT16_MAX),
      v_3_max_offset(INT16_MAX),
      v_3_max_del(INT16_MAX),
      d_3_min_offset(INT16_MAX),
      d_5_max_offset(INT16_MAX),
      j_5_min_del(INT16_MAX),
      d_5_new_offset(INT16_MAX),
      vj_check(true),
      d_3_min_del(INT16_MAX),
      base_index(INT16_MAX),
      memory_layer_offset_check1(-1),
      d_5_offset(INT16_MAX),
      d_del_opposite_side_processed(false),
      v_3_offset(INT16_MAX),
      d_5_min_del(INT16_MAX),
      previous_marginal_index(INT16_MAX),
      deletion_value(INT16_MAX),
      j_5_max_del(INT16_MAX),
      d_3_offset(INT16_MAX),
      proba_contribution(-1),
      d_3_new_offset(INT16_MAX),
      memory_layer_cs(-1),
      j_5_new_offset(INT16_MAX)
{
    this->type = Event_type::Deletion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > (-this->len_min)) {
            this->len_min = -(*iter).second.value_int;
        } else if ((*iter).second.value_int < (-this->len_max)) {
            this->len_max = -(*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Deletion::Deletion(Gene_class gene, Seq_side side,
                   unordered_map<string, Event_realization> &realizations)
    : Deletion(gene, side)
{
    this->event_realizations = realizations;

    this->type = Event_type::Deletion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > (-this->len_min)) {
            this->len_min = -(*iter).second.value_int;
        } else if ((*iter).second.value_int < (-this->len_max)) {
            this->len_max = -(*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Deletion::~Deletion()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Deletion::copy()
{
    shared_ptr<Deletion> new_deletion_p = shared_ptr<Deletion>(
            new Deletion(this->event_class, this->event_side,
                         this->event_realizations)); // FIXME remove this new for all
    // events and for the error rates
    new_deletion_p->priority = this->priority;
    new_deletion_p->nickname = this->nickname;
    new_deletion_p->fixed = this->is_fixed();
    new_deletion_p->update_event_name();
    new_deletion_p->set_event_identifier(this->event_index);
    return new_deletion_p;
}

void Deletion::add_realization(int del_number)
{

    this->Rec_Event::add_realization(
            Event_realization(to_string(del_number), del_number, "", Int_Str(), this->size()));

    if (del_number > (-this->len_min)) {
        this->len_min = (-del_number);
    } else if (del_number < (-this->len_max)) {
        this->len_max = (-del_number);
    }
    this->update_event_name();
}
/**
 * General: Loop over all possible number of deletions for a given gene on a
 * given sequence side
 *
 * Specific:
 * -First check whether any of these number of deletions is possible given the
 * current position and number of deletions on other genes -Loop over # of
 * deletions in decreasing order
 */
void Deletion::iterate(
        double &scenario_proba, Downstream_scenario_proba_bound_map &downstream_proba_map,
        const string &sequence, const Int_Str &int_sequence, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        shared_ptr<Next_event_ptr> &next_event_ptr_arr, Marginal_array_p &updated_marginals_point,
        const Marginal_array_p &model_parameters_point,
        const unordered_map<Gene_class, vector<Alignment_data>> &allowed_realizations,
        Seq_type_str_p_map &constructed_sequences, Seq_offsets_map &seq_offsets,
        shared_ptr<Error_rate> &error_rate_p, map<size_t, shared_ptr<Counter>> &counters_list,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        Safety_bool_map &safety_set, Mismatch_vectors_map &mismatches_lists,
        double &seq_max_prob_scenario, double &proba_threshold_factor)
{
    int deletion_value;
    double new_scenario_proba;
    double proba_contribution;
    double scenario_upper_bound_proba;
    int new_index;
    Int_Str new_str;
    Int_Str tmp_str;
    vector<int> mismatches_vector;
    vector<int>::iterator mis_iter;
    bool end_reached;
    int endogeneous_mismatches;

    // Helper lambda for transversions on Int_Str
    auto make_transversions = [](Int_Str &seq) {
        for (auto &val : seq) {
            if (val == 0)
                val = 3;
            else if (val == 1)
                val = 2;
            else if (val == 2)
                val = 1;
            else if (val == 3)
                val = 0;
            else if (val == 4)
                val = 5;
            else if (val == 5)
                val = 4;
        }
    };

    // Get base index
    int base_index = base_index_map.at(this->event_index);

    // Get previous string and mismatches
    Int_Str &previous_str =
            (*constructed_sequences.at(this->sequence_type_id, memory_layer_cs - 1));
    vector<int> &current_mismatch_list =
            *mismatches_lists.at(this->sequence_type_id, memory_layer_mismatches - 1);

    // Iterate over realizations
    for (forward_list<Event_realization>::const_iterator iter = (*this).int_value_and_index.begin();
         iter != (*this).int_value_and_index.end(); ++iter) {

        deletion_value = (*iter).value_int;

        // Check if deletion is possible (length check)
        if (deletion_value > 0 && (int)previous_str.size() < deletion_value) {
            continue;
        }
        if (deletion_value < 0 && -deletion_value > (int)previous_str.size()) {
            continue;
        }

        int current_offset;
        int new_offset;
        bool unsafe = false;

        // Calculate new offset and check safety
        if (this->event_side == Three_prime) {
            current_offset = seq_offsets.at(this->sequence_type_id, Three_prime,
                                            memory_layer_offset_del - 1);
            new_offset = current_offset - deletion_value;

            // Safety check against downstream neighbors
            if (downstream_chosen) {
                for (const auto &neighbor : active_downstream_junctions) {
                    int neighbor_5_offset =
                            seq_offsets.at(neighbor.first, Five_prime, memory_layer_offset_check1);
                    if (new_offset - neighbor_5_offset >= neighbor_min_del) {
                        unsafe = true;
                        break;
                    }
                }
            }
        } else { // Five_prime
            current_offset =
                    seq_offsets.at(this->sequence_type_id, Five_prime, memory_layer_offset_del - 1);
            new_offset = current_offset + deletion_value;

            if (new_offset >= (int)int_sequence.size()) {
                unsafe = true;
            } else if (upstream_chosen) {
                for (const auto &neighbor : active_upstream_junctions) {
                    int neighbor_3_offset =
                            seq_offsets.at(neighbor.first, Three_prime, memory_layer_offset_check1);
                    if (neighbor_3_offset - new_offset >= neighbor_min_del) {
                        unsafe = true;
                        break;
                    }
                }
            }
        }

        if (unsafe)
            continue;

        // Store deletion value and index
        current_realizations_index_vec[0] = (*iter).index;
        new_index = base_index + (*iter).index;
        new_scenario_proba = scenario_proba;
        proba_contribution = 1;

        this->iterate_common(iter, base_index_map, offset_map, model_parameters_point);

        // Apply Deletion / Palindrome
        if (deletion_value >= 0) {
            if (this->event_side == Three_prime) {
                new_str = previous_str.substr(0, previous_str.size() - deletion_value);

                // Update mismatches (keep <= new_offset)
                if (!current_mismatch_list.empty()) {
                    mis_iter = current_mismatch_list.end();
                    mis_iter--;
                    end_reached = false;
                    while ((*mis_iter) > new_offset) {
                        if (mis_iter == current_mismatch_list.begin()) {
                            end_reached = true;
                            break;
                        } else {
                            --mis_iter;
                        }
                    }
                    if (end_reached) {
                        mismatches_vector.clear();
                    } else {
                        ++mis_iter;
                        mismatches_vector.assign(current_mismatch_list.begin(), mis_iter);
                    }
                } else {
                    mismatches_vector.clear();
                }
            } else { // Five_prime
                new_str = previous_str.substr(deletion_value);

                // Update mismatches (keep >= new_offset)
                if (!current_mismatch_list.empty()) {
                    mis_iter = current_mismatch_list.begin();
                    end_reached = false;
                    while ((*mis_iter) < new_offset) {
                        ++mis_iter;
                        if (mis_iter == current_mismatch_list.end()) {
                            end_reached = true;
                            break;
                        }
                    }
                    if (end_reached) {
                        mismatches_vector.clear();
                    } else {
                        mismatches_vector.assign(mis_iter, current_mismatch_list.end());
                    }
                } else {
                    mismatches_vector.clear();
                }
            }
        } else { // Negative (Palindrome)
            if (this->event_side == Three_prime) {
                if (new_offset < (int)sequence.size()) {
                    tmp_str = previous_str.substr(previous_str.size() + deletion_value);
                    reverse(tmp_str.begin(), tmp_str.end());
                    make_transversions(tmp_str);
                    new_str = previous_str + tmp_str;

                    mismatches_vector = current_mismatch_list;
                    for (int i = 0; i != -deletion_value; ++i) {
                        if (current_offset + 1 + i >= 0
                            && current_offset + 1 + i < int_sequence.size()
                            && !comp_nt_int(tmp_str[i], int_sequence.at(current_offset + 1 + i))) {
                            mismatches_vector.push_back(current_offset + 1 + i);
                        }
                    }
                } else {
                    continue;
                }
            } else { // Five_prime
                if (new_offset >= 0) {
                    tmp_str = previous_str.substr(0, -deletion_value);
                    reverse(tmp_str.begin(), tmp_str.end());
                    make_transversions(tmp_str);
                    new_str = tmp_str + previous_str;

                    mismatches_vector = current_mismatch_list;
                    for (int i = 0; i != -deletion_value; ++i) {
                        if (new_offset + i < int_sequence.size()
                            && !comp_nt_int(tmp_str[i], int_sequence.at(new_offset + i))) {
                            mismatches_vector.push_back(new_offset + i);
                        }
                    }
                    sort(mismatches_vector.begin(), mismatches_vector.end());
                } else {
                    continue;
                }
            }
        }

        // Update State
        constructed_sequences.set_value(this->sequence_type_id, &new_str, memory_layer_cs);
        seq_offsets.set_value(this->sequence_type_id, this->event_side, new_offset,
                              memory_layer_offset_del);
        mismatches_lists.set_value(this->sequence_type_id, &mismatches_vector,
                                   memory_layer_mismatches);

        // Downstream Probability
        scenario_upper_bound_proba = new_scenario_proba;

        // Junctions
        if (this->event_side == Three_prime) {
            if (downstream_chosen) {
                for (const auto &neighbor : active_downstream_junctions) {
                    int neighbor_5_offset =
                            seq_offsets.at(neighbor.first, Five_prime, memory_layer_offset_check1);
                    int junction_len = neighbor_5_offset - new_offset - 1;

                    if (junction_length_best_proba_maps[neighbor.second].count(junction_len) <= 0) {
                        unsafe = true;
                        break;
                    }
                    downstream_proba_map.set_value(
                            neighbor.second,
                            junction_length_best_proba_maps[neighbor.second].at(junction_len),
                            memory_layer_proba_map_junction);
                }
            }
        } else { // Five_prime
            if (upstream_chosen) {
                for (const auto &neighbor : active_upstream_junctions) {
                    int neighbor_3_offset =
                            seq_offsets.at(neighbor.first, Three_prime, memory_layer_offset_check1);
                    int junction_len = new_offset - neighbor_3_offset - 1;

                    if (junction_length_best_proba_maps[neighbor.second].count(junction_len) <= 0) {
                        unsafe = true;
                        break;
                    }
                    downstream_proba_map.set_value(
                            neighbor.second,
                            junction_length_best_proba_maps[neighbor.second].at(junction_len),
                            memory_layer_proba_map_junction);
                }
            }
        }

        if (unsafe)
            continue;

        // Sequence Mismatches
        if (opposite_side_processed) {
            endogeneous_mismatches = mismatches_vector.size();
            downstream_proba_map.set_value(
                    this->sequence_type_id,
                    error_rate_p->get_err_rate_upper_bound(endogeneous_mismatches,
                                                           new_str.size() - endogeneous_mismatches),
                    memory_layer_proba_map_seq);
        } else {
            downstream_proba_map.set_value(this->sequence_type_id, 1.0, memory_layer_proba_map_seq);
        }

        // Multiply all downstream probas
        downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                          current_downstream_proba_memory_layers);

        if (scenario_upper_bound_proba < (seq_max_prob_scenario * proba_threshold_factor)) {
            if (deletion_value >= 0)
                break; // Optimization for positive deletions
            else
                continue;
        }

        new_scenario_proba *= proba_contribution;
        scenario_upper_bound_proba = new_scenario_proba;

        // Re-calculate with length proba
        if (this->event_side == Three_prime) {
            if (downstream_chosen) {
                for (const auto &neighbor : active_downstream_junctions) {
                    int neighbor_5_offset =
                            seq_offsets.at(neighbor.first, Five_prime, memory_layer_offset_check1);
                    int junction_len = neighbor_5_offset - new_offset - 1;
                    downstream_proba_map.set_value(
                            neighbor.second,
                            junction_length_best_proba_maps[neighbor.second].at(junction_len),
                            memory_layer_proba_map_junction);
                }
            }
        } else {
            if (upstream_chosen) {
                for (const auto &neighbor : active_upstream_junctions) {
                    int neighbor_3_offset =
                            seq_offsets.at(neighbor.first, Three_prime, memory_layer_offset_check1);
                    int junction_len = new_offset - neighbor_3_offset - 1;
                    downstream_proba_map.set_value(
                            neighbor.second,
                            junction_length_best_proba_maps[neighbor.second].at(junction_len),
                            memory_layer_proba_map_junction);
                }
            }
        }

        downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                          current_downstream_proba_memory_layers);

        if (scenario_upper_bound_proba < (seq_max_prob_scenario * proba_threshold_factor)) {
            continue;
        }

        Rec_Event::iterate_wrap_up(new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                                   base_index_map, offset_map, next_event_ptr_arr,
                                   updated_marginals_point, model_parameters_point,
                                   allowed_realizations, constructed_sequences, seq_offsets,
                                   error_rate_p, counters_list, events_map, safety_set,
                                   mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
    }
}

/*
 * This short method performs the iterate operations common to all Rec_event
 * (modify index map and fetch realization probability)
 */
void Deletion::iterate_common(
        forward_list<Event_realization>::const_iterator &iter, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        const Marginal_array_p &model_parameters_point)
{

    proba_contribution = Rec_Event::iterate_common((*iter).index, base_index, base_index_map,
                                                   model_parameters_point);
}

queue<int> Deletion::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        unordered_map<int, string> &constructed_sequences, mt19937_64 &generator) const
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

            if (constructed_sequences.count(this->sequence_type_id)) {
                string &seq = constructed_sequences.at(this->sequence_type_id);
                int val = (*iter).second.value_int;

                if (val >= 0) {
                    if (this->event_side == Three_prime) {
                        if (seq.size() >= val)
                            seq.erase(seq.size() - val);
                    } else { // Five_prime
                        if (seq.size() >= val)
                            seq.erase(0, val);
                    }
                } else { // Negative (Palindrome)
                    string gen_tmp_str;
                    if (this->event_side == Three_prime) {
                        gen_tmp_str = seq.substr(seq.size() + val);
                        reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                        make_transversions(gen_tmp_str, false);
                        seq += gen_tmp_str;
                    } else { // Five_prime
                        gen_tmp_str = seq.substr(0, -val);
                        reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                        make_transversions(gen_tmp_str, false);
                        seq = gen_tmp_str + seq;
                    }
                }
            }

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

void Deletion::write2txt(ofstream &outfile)
{
    outfile << "#Deletion;" << event_class << ";" << event_side << ";" << priority << ";"
            << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Deletion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map,
        Seq_type_str_p_map &constructed_sequences, Safety_bool_map &safety_set,
        shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    // Initialize realizations
    int_value_and_index.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 (*this).event_realizations.begin();
         iter != (*this).event_realizations.end(); ++iter) {
        int_value_and_index.push_front((*iter).second);
    }
    int_value_and_index.sort(del_numb_compare);

    // Determine sequence_type_id
    // Try to find type ID from nickname first (e.g. "D1_gene")
    string name_to_check = this->nickname;
    if (name_to_check.empty()) {
        // Fallback to class-based naming if nickname is empty
        if (this->event_class == V_gene)
            name_to_check = "V_gene";
        else if (this->event_class == D_gene)
            name_to_check = "D_gene"; // Default D
        else if (this->event_class == J_gene)
            name_to_check = "J_gene";
    }

    this->sequence_type_id = SequenceTypeRegistry::get_instance().get_type_id(name_to_check);

    // Fallback for backward compatibility if not found by name
    if (this->sequence_type_id == -1) {
        if (this->event_class == V_gene)
            this->sequence_type_id = SequenceTypeRegistry::V_GENE_SEQ;
        else if (this->event_class == J_gene)
            this->sequence_type_id = SequenceTypeRegistry::J_GENE_SEQ;
        else if (this->event_class == D_gene)
            this->sequence_type_id = SequenceTypeRegistry::D_GENE_SEQ;
    }

    // Find neighbors and populate active junctions
    auto upstream_neighbors =
            SequenceTypeRegistry::get_instance().get_upstream_neighbors(this->sequence_type_id);
    auto downstream_neighbors =
            SequenceTypeRegistry::get_instance().get_downstream_neighbors(this->sequence_type_id);

    active_upstream_junctions.clear();
    upstream_exists = false;
    upstream_chosen = false;

    for (const auto &neighbor : upstream_neighbors) {
        if (neighbor.exists) {
            upstream_exists = true;
            active_upstream_junctions.push_back({ neighbor.neighbor_type, neighbor.junction_type });

            // Check if neighbor is chosen
            if (events_map.count(tuple<Event_type, int, Seq_side>(
                        GeneChoice_t, neighbor.neighbor_type, Undefined_side))) {
                auto event = events_map.at(tuple<Event_type, int, Seq_side>(
                        GeneChoice_t, neighbor.neighbor_type, Undefined_side));
                if (processed_events.count(event->get_name())) {
                    upstream_chosen = true;
                }
            }
        }
    }

    active_downstream_junctions.clear();
    downstream_exists = false;
    downstream_chosen = false;

    for (const auto &neighbor : downstream_neighbors) {
        if (neighbor.exists) {
            downstream_exists = true;
            active_downstream_junctions.push_back(
                    { neighbor.neighbor_type, neighbor.junction_type });

            // Check if neighbor is chosen
            if (events_map.count(tuple<Event_type, int, Seq_side>(
                        GeneChoice_t, neighbor.neighbor_type, Undefined_side))) {
                auto event = events_map.at(tuple<Event_type, int, Seq_side>(
                        GeneChoice_t, neighbor.neighbor_type, Undefined_side));
                if (processed_events.count(event->get_name())) {
                    downstream_chosen = true;
                }
            }
        }
    }

    // Request Memory Layers
    // 1. Offsets
    if (this->event_side == Three_prime) {
        seq_offsets.request_memory_layer(this->sequence_type_id, Three_prime);
        memory_layer_offset_del =
                seq_offsets.get_current_memory_layer(this->sequence_type_id, Three_prime);
    } else if (this->event_side == Five_prime) {
        seq_offsets.request_memory_layer(this->sequence_type_id, Five_prime);
        memory_layer_offset_del =
                seq_offsets.get_current_memory_layer(this->sequence_type_id, Five_prime);
    }

    // 2. Mismatches and Constructed Sequences
    mismatches_list.request_memory_layer(this->sequence_type_id);
    this->memory_layer_mismatches =
            mismatches_list.get_current_memory_layer(this->sequence_type_id);

    constructed_sequences.request_memory_layer(this->sequence_type_id);
    this->memory_layer_cs = constructed_sequences.get_current_memory_layer(this->sequence_type_id);

    // 3. Downstream Proba Map
    downstream_proba_map.request_memory_layer(this->sequence_type_id);
    this->memory_layer_proba_map_seq =
            downstream_proba_map.get_current_memory_layer(this->sequence_type_id);

    // 4. Junctions
    // Request memory for all active downstream junctions
    for (const auto &junction : active_downstream_junctions) {
        downstream_proba_map.request_memory_layer(junction.second);
        // We store the memory layer for the first one found for now, or use a map if needed.
        // The original code used memory_layer_proba_map_junction.
        // We should probably store a map of memory layers if we have multiple downstream junctions.
        // For now, assuming one active path or using the last one (which is risky).
        // But wait, iterate() uses memory_layer_proba_map_junction.
        // If we have multiple downstream junctions (e.g. D1->D2 and D1->J?), we need to handle this.
        // In linear topology, there is only one downstream neighbor.
        // So we can just use the first one.
        memory_layer_proba_map_junction =
                downstream_proba_map.get_current_memory_layer(junction.second);
    }

    // 5. Neighbor Offsets (for overlap checks)
    // We need to request memory layers for neighbor offsets to check overlaps.
    // Original code used memory_layer_offset_check1/2.
    // We can iterate neighbors and request.
    // But we need to store them to access in iterate().
    // Since we don't have a map for memory layers in the class, we might be limited.
    // However, we can use the fact that we only check overlap if the neighbor is chosen.
    // And usually there is only one relevant neighbor for overlap check (the one we are deleting towards).
    // If 3' del, we check downstream neighbor 5' offset.
    // If 5' del, we check upstream neighbor 3' offset.

    if (this->event_side == Three_prime) {
        for (const auto &neighbor : active_downstream_junctions) {
            if (downstream_chosen) { // Should check specific neighbor chosen status if multiple
                seq_offsets.request_memory_layer(neighbor.first, Five_prime);
                memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(
                        neighbor.first, Five_prime); // Reusing check1
            }
        }
    } else if (this->event_side == Five_prime) {
        for (const auto &neighbor : active_upstream_junctions) {
            if (upstream_chosen) {
                seq_offsets.request_memory_layer(neighbor.first, Three_prime);
                memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(
                        neighbor.first, Three_prime); // Reusing check1
            }
        }
    }

    // Note: The original code had check1 and check2 for D gene (checking both V and J).
    // But D gene has 5' del (checks V) and 3' del (checks J).
    // They are separate events (Deletion 5' and Deletion 3').
    // So for a single Deletion event, we only check one side.
    // The original code for D_gene 5' del checked V (check1).
    // The original code for D_gene 3' del checked J (check2).
    // So we only need one check variable per Deletion event.

    // 6. Safety
    // Skipping safety optimization for generic case for now.
    this->Rec_Event::initialize_event(processed_events, events_map, offset_map,
                                      downstream_proba_map, constructed_sequences, safety_set,
                                      error_rate_p, mismatches_list, seq_offsets, index_map);

    // 7. Get Neighbor Deletion Ranges
    neighbor_min_del = 0;
    neighbor_max_del = 0;
    opposite_side_processed = true; // Default to true if no opposite event exists

    if (this->event_side == Three_prime) {
        // Check for opposite side (Five_prime) deletion on same sequence
        if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id,
                                                              Five_prime))) {
            auto opp_event = events_map.at(tuple<Event_type, int, Seq_side>(
                    Deletion_t, this->sequence_type_id, Five_prime));
            if (processed_events.count(opp_event->get_name()) == 0) {
                opposite_side_processed = false;
            }
        }

        // Check downstream neighbors (usually one)
        for (const auto &neighbor : active_downstream_junctions) {
            // We need to find the Deletion event for this neighbor at Five_prime side
            if (events_map.count(
                        tuple<Event_type, int, Seq_side>(Deletion_t, neighbor.first, Five_prime))) {
                auto del_event = events_map.at(
                        tuple<Event_type, int, Seq_side>(Deletion_t, neighbor.first, Five_prime));
                if (processed_events.count(del_event->get_name()) == 0) {
                    neighbor_min_del = del_event->get_len_max();
                    neighbor_max_del = del_event->get_len_min();
                }
            }
        }
    } else if (this->event_side == Five_prime) {
        // Check for opposite side (Three_prime) deletion on same sequence
        if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id,
                                                              Three_prime))) {
            auto opp_event = events_map.at(tuple<Event_type, int, Seq_side>(
                    Deletion_t, this->sequence_type_id, Three_prime));
            if (processed_events.count(opp_event->get_name()) == 0) {
                opposite_side_processed = false;
            }
        }

        // Check upstream neighbors
        for (const auto &neighbor : active_upstream_junctions) {
            if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, neighbor.first,
                                                                  Three_prime))) {
                auto del_event = events_map.at(
                        tuple<Event_type, int, Seq_side>(Deletion_t, neighbor.first, Three_prime));
                if (processed_events.count(del_event->get_name()) == 0) {
                    neighbor_min_del = del_event->get_len_max();
                    neighbor_max_del = del_event->get_len_min();
                }
            }
        }
    }
}

void Deletion::add_to_marginals(long double scenario_proba,
                                Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        updated_marginals[this->new_index] = scenario_proba;
    } else {
        updated_marginals[this->new_index] += scenario_proba;
    }
}

string &make_transversions(string &init_sequence, bool is_int_seq)
{
    if (is_int_seq) {
        for (string::iterator iter = init_sequence.begin(); iter != init_sequence.end(); ++iter) {
            if ((*iter) == '0') {
                (*iter) = '3';
            } else if ((*iter) == '1') {
                (*iter) = '2';
            } else if ((*iter) == '2') {
                (*iter) = '1';
            } else if ((*iter) == '3') {
                (*iter) = '0';
            } else if ((*iter) == '4') {
                (*iter) = '5';
            } else if ((*iter) == '5') {
                (*iter) = '4';
            } else if ((*iter) == '8') {
                // Nothing to do
            } else if ((*iter) == '9') {
                // Nothing to do
            } else if ((*iter) == '14') {
                // Nothing to do
            } else {
                throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq "
                                    + init_sequence + " in make_transversions()");
            }
        }
    } else {
        for (string::iterator iter = init_sequence.begin(); iter != init_sequence.end(); ++iter) {
            if ((*iter) == 'A') {
                (*iter) = 'T';
            } else if ((*iter) == 'C') {
                (*iter) = 'G';
            } else if ((*iter) == 'G') {
                (*iter) = 'C';
            } else if ((*iter) == 'T') {
                (*iter) = 'A';
            } else {
                throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq "
                                    + init_sequence + " in make_transversions()");
            }
        }
    }
    return init_sequence;
}

Int_Str &make_transversions(Int_Str &init_sequence)
{

    for (Int_Str::iterator iter = init_sequence.begin(); iter != init_sequence.end(); ++iter) {
        if ((*iter) == 0) {
            (*iter) = 3;
        } else if ((*iter) == 1) {
            (*iter) = 2;
        } else if ((*iter) == 2) {
            (*iter) = 1;
        } else if ((*iter) == 3) {
            (*iter) = 0;
        } else if ((*iter) == 4) {
            (*iter) = 5;
        } else if ((*iter) == 5) {
            (*iter) = 4;
        } else if ((*iter) == 8) {
            // Nothing to do
        } else if ((*iter) == 9) {
            // Nothing to do
        } else if ((*iter) == 14) {
            // Nothing to do
        } else {

            string error_str("Unknown int nucleotide " + to_string((*iter)) + " in seq ");
            for (Int_Str::iterator jiter = init_sequence.begin(); jiter != init_sequence.end();
                 ++jiter) {
                error_str += to_string((*jiter));
            }
            error_str += " in make_transversions()";
            throw runtime_error(error_str);
        }
    }

    return init_sequence;
}

bool del_numb_compare(const Event_realization &real1, const Event_realization &real2)
{
    return real1.value_int > real2.value_int;
}

bool Deletion::has_effect_on(int seq_type) const
{
    // Use topology from SequenceTypeRegistry to determine if this deletion affects
    // the given junction type. This supports both standard V-D-J and tandem D models.
    //
    // A deletion affects a junction if:
    // - 3' deletion: the junction connects this gene to a downstream neighbor
    // - 5' deletion: the junction connects this gene to an upstream neighbor

    const auto &registry = SequenceTypeRegistry::get_instance();

    if (this->event_side == Three_prime) {
        // 3' deletion affects downstream junctions
        auto downstream_neighbors = registry.get_downstream_neighbors(this->sequence_type_id);
        for (const auto &neighbor : downstream_neighbors) {
            if (neighbor.exists && neighbor.junction_type == seq_type) {
                return true;
            }
        }
    } else if (this->event_side == Five_prime) {
        // 5' deletion affects upstream junctions
        auto upstream_neighbors = registry.get_upstream_neighbors(this->sequence_type_id);
        for (const auto &neighbor : upstream_neighbors) {
            if (neighbor.exists && neighbor.junction_type == seq_type) {
                return true;
            }
        }
    }

    return false;
}

void Deletion::iterate_initialize_Len_proba(
        int considered_junction, std::map<int, double> &length_best_proba_map,
        std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
        const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
        Seq_type_str_p_map &constructed_sequences, int &seq_len /*=0*/) const
{

    if (this->has_effect_on(considered_junction)) {
        base_index = base_index_map.at(this->event_index, 0);
        for (unordered_map<string, Event_realization>::const_iterator iter =
                     this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {

            /*		//Update base index map
                      for(forward_list<tuple<int,int,int>>::const_iterator jiter
         = memory_and_offsets.begin() ; jiter!=memory_and_offsets.end() ;
         ++jiter){
                              //Get previous index for the considered event
                              int previous_index =
         base_index_map.at(get<0>(*jiter),get<1>(*jiter)-1);
                              //Update the index given the realization and the
         offset previous_index += iter->second.index *get<2>(*jiter);
                              //Set the value
                              base_index_map.set_value(get<0>(*jiter) ,
         previous_index , get<1>(*jiter));
                      }*/

            // Get the max proba for this realization (in case the event is child of
            // another)
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()]
                    > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index
                                                            + i * this->size()];
                }
            }
            // Update the length and the probability in the recursive call
            Rec_Event::iterate_initialize_Len_proba_wrap_up(
                    considered_junction, length_best_proba_map, model_queue,
                    scenario_proba * real_max_proba, model_parameters_point, base_index_map,
                    constructed_sequences, seq_len - (*iter).second.value_int);
        }
    } else {
        Rec_Event::iterate_initialize_Len_proba_wrap_up(
                considered_junction, length_best_proba_map, model_queue, scenario_proba,
                model_parameters_point, base_index_map, constructed_sequences, seq_len);
    }
}

void Deletion::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                          const Marginal_array_p &model_parameters_point,
                                          Index_map &base_index_map)
{
    Seq_type_str_p_map constructed_sequences(SequenceTypeRegistry::get_instance().size());
    junction_length_best_proba_maps.clear();

    // Iterate over active downstream junctions (only if downstream neighbor is chosen)
    if (downstream_chosen) {
        for (const auto &junction : active_downstream_junctions) {
            int junction_type = junction.second;
            junction_length_best_proba_maps[junction_type].clear();

            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(
                    junction_type, junction_length_best_proba_maps[junction_type], model_queue,
                    init_proba, model_parameters_point, base_index_map, constructed_sequences);
        }
    }

    // Iterate over active upstream junctions (only if upstream neighbor is chosen)
    if (upstream_chosen) {
        for (const auto &junction : active_upstream_junctions) {
            int junction_type = junction.second;
            junction_length_best_proba_maps[junction_type].clear();

            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(
                    junction_type, junction_length_best_proba_maps[junction_type], model_queue,
                    init_proba, model_parameters_point, base_index_map, constructed_sequences);
        }
    }
}
