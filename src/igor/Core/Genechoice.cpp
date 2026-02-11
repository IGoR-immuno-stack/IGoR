/*
 * Genechoice.cpp
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

#include <igor/Core/Aligner.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/SequenceTypes.h>

using namespace std;

Gene_choice::Gene_choice() : Gene_choice(Undefined_gene)
{
    this->type = Event_type::GeneChoice_t;
    this->update_event_name();
}

Gene_choice::Gene_choice(int gene)
    : Rec_Event(gene, Undefined_side),
      vd_check(true),
      vj_check(true),
      dj_check(true),
      d_5_min_offset(INT16_MAX),
      d_5_max_offset(INT16_MAX),
      j_5_min_offset(INT16_MAX),
      j_5_max_offset(INT16_MAX),
      v_5_off(INT16_MAX),
      v_3_off(INT16_MAX),
      d_offset(INT16_MAX),
      j_offset(INT16_MAX),
      v_offset(INT16_MAX),
      v_3_min_offset(INT16_MAX),
      v_3_max_offset(INT16_MAX),
      d_3_off(INT16_MAX),
      d_5_off(INT16_MAX),
      d_3_min_offset(INT16_MAX),
      d_3_max_offset(INT16_MAX),
      j_5_off(INT16_MAX),
      no_d_align(true),
      d_size(INT16_MAX),
      d_full_3_offset(INT16_MAX),
      base_index(-1),
      new_scenario_proba(-1),
      new_tmp_err_w_proba(-1),
      proba_contribution(-1),
      new_index(-1),
      alignment_offset_p(NULL),
      memory_layer_cs(-1),
      memory_layer_mismatches(-1),
      memory_layer_safety_1(-1),
      memory_layer_safety_2(-1),
      memory_layer_off_threep(-1),
      memory_layer_off_fivep(-1),
      memory_layer_offset_check1(-1),
      memory_layer_offset_check2(-1),
      v_chosen(false),
      v_choice_exist(true),
      d_chosen(false),
      d_choice_exist(false),
      j_chosen(false),
      j_choice_exist(true),
      d_5_max_del(INT16_MIN),
      d_5_min_del(INT16_MAX),
      d_5_real_max_del(INT16_MIN),
      j_5_max_del(INT16_MIN),
      j_5_min_del(INT16_MIN),
      v_3_max_del(INT16_MIN),
      v_3_min_del(INT16_MAX),
      d_3_max_del(INT16_MIN),
      d_3_min_del(INT16_MAX)
{
    this->type = Event_type::GeneChoice_t;
    this->sequence_type_id = gene; // Set to gene for proper events_map lookups
    this->update_event_name();
}

/*
 * Should probably be avoided, default initialize the event and use
 * add_event_realization() instead unless you're sure about the index fields in
 * the Event_realizations instances
 */
Gene_choice::Gene_choice(int gene, unordered_map<string, Event_realization> &realizations) : Gene_choice(gene)
{
    this->event_realizations = realizations;

    this->type = Event_type::GeneChoice_t;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        int str_len = (*iter).second.value_str.length();
        if (str_len > this->len_max) {
            this->len_max = str_len;
        } else if (str_len < this->len_min) {
            this->len_min = str_len;
        }
    }
    this->update_event_name();
}

Gene_choice::Gene_choice(int gene, vector<pair<string, string>> genomic_sequences) : Gene_choice(gene)
{
    this->type = Event_type::GeneChoice_t;
    for (vector<pair<string, string>>::const_iterator seq_it = genomic_sequences.begin();
         seq_it != genomic_sequences.end(); ++seq_it) {
        int str_len = (*seq_it).second.length();
        if (str_len > this->len_max) {
            this->len_max = str_len;
        } else if (str_len < this->len_min) {
            this->len_min = str_len;
        }
        this->add_realization((*seq_it).first, (*seq_it).second);
    }
    this->update_event_name();
}

Gene_choice::~Gene_choice()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Gene_choice::copy()
{
    shared_ptr<Gene_choice> new_gene_choice_p =
            shared_ptr<Gene_choice>(new Gene_choice(this->event_class, this->event_realizations));
    new_gene_choice_p->priority = this->priority;
    new_gene_choice_p->nickname = this->nickname;
    new_gene_choice_p->fixed = this->fixed;
    // Copy tandem D support fields
    new_gene_choice_p->sequence_type_id = this->sequence_type_id;
    new_gene_choice_p->upstream_seq_type = this->upstream_seq_type;
    new_gene_choice_p->downstream_seq_type = this->downstream_seq_type;
    new_gene_choice_p->upstream_exists = this->upstream_exists;
    new_gene_choice_p->downstream_exists = this->downstream_exists;
    new_gene_choice_p->upstream_chosen = this->upstream_chosen;
    new_gene_choice_p->downstream_chosen = this->downstream_chosen;
    new_gene_choice_p->upstream_ins_type = this->upstream_ins_type;
    new_gene_choice_p->downstream_ins_type = this->downstream_ins_type;

    new_gene_choice_p->name = this->name; // Copy name directly to ensure stability
    new_gene_choice_p->set_event_identifier(this->event_index);
    return new_gene_choice_p;
}

bool Gene_choice::add_realization(string gene_name, string gene_sequence)
{
    int str_len = gene_sequence.length();
    if (str_len > this->len_max) {
        this->len_max = str_len;
    } else if (str_len < this->len_min) {
        this->len_min = str_len;
    }
    this->Rec_Event::add_realization(Event_realization(gene_name, INT16_MAX, gene_sequence, nt2int(gene_sequence),
                                                       this->event_realizations.size())); //FIXME nonsense new
    this->update_event_name();
    return 1;
}

void Gene_choice::set_genomic_templates(const vector<pair<string, string>> &genomic_templates)
{
    // First remove previous realizations
    this->event_realizations.clear();
    for (vector<pair<string, string>>::const_iterator iter = genomic_templates.begin(); iter != genomic_templates.end();
         ++iter) {
        this->add_realization((*iter).first, (*iter).second);
    }
}

void Gene_choice::iterate(
        double &scenario_proba, Downstream_scenario_proba_bound_map &downstream_proba_map, const string &sequence,
        const Int_Str &int_sequence, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        shared_ptr<Next_event_ptr> &next_event_ptr_arr, Marginal_array_p &updated_marginals_pointer,
        const Marginal_array_p &model_parameters_pointer,
        const unordered_map<int, vector<Alignment_data>> &allowed_realizations,
        Seq_type_str_p_map &constructed_sequences, Seq_offsets_map &seq_offsets, shared_ptr<Error_rate> &error_rate_p,
        map<size_t, shared_ptr<Counter>> &counters_list,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        Safety_bool_map &safety_set, Mismatch_vectors_map &mismatches_lists, double &seq_max_prob_scenario,
        double &proba_threshold_factor)
{
    base_index = base_index_map.at(this->event_index);

    // Generic Iteration Logic

    // Determine if we are V, D, or J like based on neighbors
    bool is_v_like = (upstream_seq_type == -1); // No upstream
    bool is_j_like = (downstream_seq_type == -1); // No downstream
    bool is_d_like = (!is_v_like && !is_j_like); // Middle

    if (allowed_realizations.count(this->event_class) == 0) {
        if (this->event_class == D_gene) {
            return;
        }
        return;
    }

    // --- V-like Logic (Upstream-most) ---
    if (is_v_like) {
        // Iterate over allowed realizations
        for (const auto &iter : allowed_realizations.at(this->event_class)) {
            if (this->event_realizations.count(iter.gene_name) == 0) {
                continue;
            }
            const Int_Str &gene_seq = this->event_realizations.at(iter.gene_name).value_str_int;

            constructed_sequences.set_value(this->sequence_type_id, gene_seq, memory_layer_cs);

            int v_3_off = iter.offset + gene_seq.size() - 1;
            int v_5_off = iter.offset;

            // Check downstream safety (e.g. VD or VJ)
            if (downstream_exists && downstream_chosen) {
                // Check against downstream neighbor
                int ds_5_offset = seq_offsets.at(downstream_seq_type, Five_prime,
                                                 memory_layer_offset_check2); // Using check2 for downstream

                // Max deletion check (Impossible to fit)
                if ((v_3_off + v_3_max_del) >= (ds_5_offset - j_5_max_del)) {
                    // Overlap even with max deletions -> Bad
                    continue;
                }

                // Min deletion check
                if ((v_3_off - v_3_min_del) < (ds_5_offset - j_5_min_del)) {
                    // No overlap with min deletions -> Safe
                    safety_set.set_value(memory_layer_safety_downstream, true, memory_layer_safety_downstream);
                } else {
                    // In deletion range -> Unsafe (potential overlap)
                    safety_set.set_value(memory_layer_safety_downstream, false, memory_layer_safety_downstream);
                }
            } else if (downstream_exists) {
                if (memory_layer_safety_downstream != -1)
                    safety_set.set_value(memory_layer_safety_downstream, true, memory_layer_safety_downstream);
            }

            // Update indices and probas
            current_realizations_index_vec[0] = this->event_realizations.at(iter.gene_name).index;
            double proba_contribution = iterate_common(1.0, current_realizations_index_vec[0], base_index,
                                                       base_index_map, offset_map, model_parameters_pointer);
            double new_scenario_proba = scenario_proba * proba_contribution;

            // Set offsets
            seq_offsets.set_value(this->sequence_type_id, Three_prime, v_3_off, memory_layer_off_threep);
            seq_offsets.set_value(this->sequence_type_id, Five_prime, v_5_off, memory_layer_off_fivep);
            mismatches_lists.set_value(this->sequence_type_id, &iter.mismatches, memory_layer_mismatches);

            // Downstream proba
            double scenario_upper_bound_proba = new_scenario_proba;

            if (downstream_exists && downstream_chosen) {
                int ds_5_offset = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_offset_check2);
                int junction_len = ds_5_offset - v_3_off - 1;

                if (junction_length_best_proba_maps.count(downstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(downstream_ins_type).count(junction_len) == 0) {
                    continue;
                }

                downstream_proba_map.set_value(downstream_ins_type,
                                               junction_length_best_proba_maps.at(downstream_ins_type).at(junction_len),
                                               memory_layer_proba_map_junction);
            }

            // Mismatches
            downstream_proba_map.set_value(this->sequence_type_id, 1.0,
                                           memory_layer_proba_map_seq); // Simplified

            // Recurse
            downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                              current_downstream_proba_memory_layers.data());

            if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
                Rec_Event::iterate_wrap_up(new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                                           base_index_map, offset_map, next_event_ptr_arr, updated_marginals_pointer,
                                           model_parameters_pointer, allowed_realizations, constructed_sequences,
                                           seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                                           mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
            }
        }
    }
    // --- D-like Logic (Middle) ---
    else if (is_d_like) {
        // Check Upstream Safety
        bool us_check = false;
        int us_3_offset = 0;
        if (upstream_exists && upstream_chosen) {
            us_3_offset = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_offset_check1);
            us_check = true;
        } else {
            if (upstream_exists && memory_layer_safety_upstream != -1) {
                safety_set.set_value(memory_layer_safety_upstream, true, memory_layer_safety_upstream);
            }
        }

        // Check Downstream Safety
        bool ds_check = false;
        int ds_5_offset = 0;
        if (downstream_exists && downstream_chosen) {
            ds_5_offset = seq_offsets.at(downstream_seq_type, Five_prime, memory_layer_offset_check2);
            ds_check = true;
        } else {
            if (downstream_exists) {
                if (memory_layer_safety_downstream != -1)
                    safety_set.set_value(memory_layer_safety_downstream, true, memory_layer_safety_downstream);
            }
            if (!downstream_chosen && downstream_exists) { // Assuming J-like behavior
                ds_5_offset = sequence.size() - 1;
            }
        }

        for (const auto &iter : allowed_realizations.at(this->event_class)) {
            if (this->event_realizations.count(iter.gene_name) == 0) {
                continue;
            }
            const Int_Str &gene_seq = this->event_realizations.at(iter.gene_name).value_str_int;
            constructed_sequences.set_value(this->sequence_type_id, gene_seq, memory_layer_cs);

            int d_5_off = iter.offset;
            int d_3_off = iter.offset + gene_seq.size() - 1;

            // Upstream Check
            if (us_check) {
                if ((d_5_off - d_5_max_del) <= (us_3_offset + v_3_max_del)) {
                    continue; // Overlap
                }
                if ((d_5_off - d_5_min_del) > (us_3_offset - v_3_max_del)) {
                    safety_set.set_value(memory_layer_safety_upstream, true, memory_layer_safety_upstream);
                } else {
                    safety_set.set_value(memory_layer_safety_upstream, false, memory_layer_safety_upstream);
                }

                // Check upstream junction length
                int junction_len = d_5_off - us_3_offset - 1;
                if (junction_length_best_proba_maps.count(upstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(upstream_ins_type).count(junction_len) == 0) {
                    continue;
                }
                downstream_proba_map.set_value(upstream_ins_type,
                                               junction_length_best_proba_maps.at(upstream_ins_type).at(junction_len),
                                               memory_layer_proba_map_junction_upstream);
            }

            // Downstream Check
            if (ds_check) {
                if ((d_3_off + d_3_max_del) >= (ds_5_offset - j_5_max_del)) {
                    continue; // Overlap
                }
                if ((d_3_off + d_3_min_del) < (ds_5_offset - j_5_min_del)) {
                    safety_set.set_value(memory_layer_safety_downstream, true, memory_layer_safety_downstream);
                } else {
                    safety_set.set_value(memory_layer_safety_downstream, false, memory_layer_safety_downstream);
                }

                // Check downstream junction length
                int junction_len = ds_5_offset - d_3_off - 1;
                if (junction_length_best_proba_maps.count(downstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(downstream_ins_type).count(junction_len) == 0) {
                    continue;
                }
                downstream_proba_map.set_value(downstream_ins_type,
                                               junction_length_best_proba_maps.at(downstream_ins_type).at(junction_len),
                                               memory_layer_proba_map_junction);
            }

            // Update indices and probas
            current_realizations_index_vec[0] = this->event_realizations.at(iter.gene_name).index;
            double proba_contribution = iterate_common(1.0, current_realizations_index_vec[0], base_index,
                                                       base_index_map, offset_map, model_parameters_pointer);
            double new_scenario_proba = scenario_proba * proba_contribution;

            // Set offsets
            seq_offsets.set_value(this->sequence_type_id, Five_prime, d_5_off, memory_layer_off_fivep);
            seq_offsets.set_value(this->sequence_type_id, Three_prime, d_3_off, memory_layer_off_threep);
            mismatches_lists.set_value(this->sequence_type_id, &iter.mismatches, memory_layer_mismatches);

            // Downstream proba & Mismatches (Simplified)
            double scenario_upper_bound_proba = new_scenario_proba;
            downstream_proba_map.set_value(this->sequence_type_id, 1.0, memory_layer_proba_map_seq);

            // Recurse
            downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                              current_downstream_proba_memory_layers.data());

            if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
                Rec_Event::iterate_wrap_up(new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                                           base_index_map, offset_map, next_event_ptr_arr, updated_marginals_pointer,
                                           model_parameters_pointer, allowed_realizations, constructed_sequences,
                                           seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                                           mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
            }
        }
    }
    // --- J-like Logic (Downstream-most) ---
    else if (is_j_like) {
        // Check Upstream Safety
        bool us_check = false;
        int us_3_offset = 0;
        if (upstream_exists && upstream_chosen) {
            us_3_offset = seq_offsets.at(upstream_seq_type, Three_prime, memory_layer_offset_check1);
            us_check = true;
        } else {
            if (upstream_exists && memory_layer_safety_upstream != -1) {
                safety_set.set_value(memory_layer_safety_upstream, true, memory_layer_safety_upstream);
            }
        }

        for (const auto &iter : allowed_realizations.at(this->event_class)) {
            if (this->event_realizations.count(iter.gene_name) == 0) {
                continue;
            }
            const Int_Str &gene_seq = this->event_realizations.at(iter.gene_name).value_str_int;
            constructed_sequences.set_value(this->sequence_type_id, gene_seq, memory_layer_cs);

            int j_5_off = iter.offset;
            int j_3_off = iter.offset + gene_seq.size() - 1;

            // Upstream Check
            if (us_check) {
                if ((j_5_off - j_5_max_del) <= (us_3_offset + v_3_max_del)) {
                    continue; // Overlap
                }
                if ((j_5_off - j_5_min_del) > (us_3_offset - v_3_max_del)) {
                    safety_set.set_value(memory_layer_safety_upstream, true, memory_layer_safety_upstream);
                } else {
                    safety_set.set_value(memory_layer_safety_upstream, false, memory_layer_safety_upstream);
                }

                // Check upstream junction length
                int junction_len = j_5_off - us_3_offset - 1;
                if (junction_length_best_proba_maps.count(upstream_ins_type) == 0
                    || junction_length_best_proba_maps.at(upstream_ins_type).count(junction_len) == 0) {
                    continue;
                }
                downstream_proba_map.set_value(upstream_ins_type,
                                               junction_length_best_proba_maps.at(upstream_ins_type).at(junction_len),
                                               memory_layer_proba_map_junction_upstream);
            }

            // Update indices and probas
            current_realizations_index_vec[0] = this->event_realizations.at(iter.gene_name).index;
            double proba_contribution = iterate_common(1.0, current_realizations_index_vec[0], base_index,
                                                       base_index_map, offset_map, model_parameters_pointer);
            double new_scenario_proba = scenario_proba * proba_contribution;

            // Set offsets
            seq_offsets.set_value(this->sequence_type_id, Five_prime, j_5_off, memory_layer_off_fivep);
            seq_offsets.set_value(this->sequence_type_id, Three_prime, j_3_off, memory_layer_off_threep);
            mismatches_lists.set_value(this->sequence_type_id, &iter.mismatches, memory_layer_mismatches);

            // Downstream proba & Mismatches (Simplified)
            double scenario_upper_bound_proba = new_scenario_proba;
            downstream_proba_map.set_value(this->sequence_type_id, 1.0, memory_layer_proba_map_seq);

            // Recurse
            downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                              current_downstream_proba_memory_layers.data());

            if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
                Rec_Event::iterate_wrap_up(new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                                           base_index_map, offset_map, next_event_ptr_arr, updated_marginals_pointer,
                                           model_parameters_pointer, allowed_realizations, constructed_sequences,
                                           seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                                           mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
            }
        }
    }
};

double Gene_choice::iterate_common(
        double scenario_proba, const int &gene_index, int base_index, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        const Marginal_array_p &model_parameters)
{
    return scenario_proba * Rec_Event::iterate_common(gene_index, base_index, base_index_map, model_parameters);
}

queue<int> Gene_choice::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        std::unordered_map<int, std::string> &constructed_sequences, std::mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;

    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            constructed_sequences[this->sequence_type_id] = (*iter).second.value_str;
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
void Gene_choice::write2txt(ofstream &outfile)
{
    outfile << "#GeneChoice;" << (Gene_class)event_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.name << ";" << (*iter).second.value_str << ";" << (*iter).second.index << endl;
    }
}

void Gene_choice::initialize_event(
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

        // Try nickname-based lookup first (for tandem D: "D1_gene" -> "D1_gene_seq")
        if (!this->nickname.empty()) {
            int type_id = registry.try_get_type_id(this->nickname + "_seq");
            if (type_id >= 0) {
                this->sequence_type_id = type_id;
            }
        }

        // Fall back to event_class-based assignment
        if (this->sequence_type_id == -1) {
            if (this->event_class == V_gene) {
                this->sequence_type_id = SequenceTypeRegistry::V_GENE_SEQ;
            } else if (this->event_class == D_gene) {
                this->sequence_type_id = SequenceTypeRegistry::D_GENE_SEQ;
            } else if (this->event_class == J_gene) {
                this->sequence_type_id = SequenceTypeRegistry::J_GENE_SEQ;
            }
        }
    }

    // 2. Determine Neighbors using Registry
    auto &registry = get_sequence_type_registry();
    auto upstream_neighbors = registry.get_upstream_neighbors((SequenceTypeRegistry::TypeId)this->sequence_type_id);
    auto downstream_neighbors = registry.get_downstream_neighbors((SequenceTypeRegistry::TypeId)this->sequence_type_id);

    active_upstream_junctions.clear();
    active_downstream_junctions.clear();

    // Reset
    upstream_seq_type = -1;
    downstream_seq_type = -1;
    upstream_exists = false;
    downstream_exists = false;
    upstream_chosen = false;
    downstream_chosen = false;
    upstream_ins_type = -1;
    downstream_ins_type = -1;

    // Find active upstream neighbor
    for (const auto &neighbor : upstream_neighbors) {
        auto v_status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map, processed_events);
        if (v_status.exists) {
            active_upstream_junctions.push_back({ neighbor.neighbor_type, neighbor.junction_type });

            if (!upstream_exists) {
                upstream_exists = true;
                upstream_seq_type = neighbor.neighbor_type;
                upstream_ins_type = neighbor.junction_type;

                if (v_status.chosen) {
                    upstream_chosen = true;
                }
            }
        }
    }

    // Find active downstream neighbor
    for (const auto &neighbor : downstream_neighbors) {
        auto j_status = EventUtils::check_gene_choice((Gene_class)neighbor.neighbor_type, events_map, processed_events);
        if (j_status.exists) {
            active_downstream_junctions.push_back({ neighbor.neighbor_type, neighbor.junction_type });

            if (!downstream_exists) {
                downstream_exists = true;
                downstream_seq_type = neighbor.neighbor_type;
                downstream_ins_type = neighbor.junction_type;

                if (j_status.chosen) {
                    downstream_chosen = true;
                }
            }
        }
    }

    // 3. Request Memory Layers
    seq_offsets.request_memory_layer(this->sequence_type_id, Three_prime);
    this->memory_layer_off_threep = seq_offsets.get_current_memory_layer(this->sequence_type_id, Three_prime);

    seq_offsets.request_memory_layer(this->sequence_type_id, Five_prime);
    this->memory_layer_off_fivep = seq_offsets.get_current_memory_layer(this->sequence_type_id, Five_prime);

    mismatches_list.request_memory_layer(this->sequence_type_id);
    this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(this->sequence_type_id);

    constructed_sequences.request_memory_layer(this->sequence_type_id);
    this->memory_layer_cs = constructed_sequences.get_current_memory_layer(this->sequence_type_id);

    downstream_proba_map.request_memory_layer(this->sequence_type_id);
    memory_layer_proba_map_seq = downstream_proba_map.get_current_memory_layer(this->sequence_type_id);

    // Safety and Junctions
    if (upstream_exists) {
        // Use upstream_ins_type as key for safety
        safety_set.request_memory_layer(upstream_ins_type);
        memory_layer_safety_upstream = safety_set.get_current_memory_layer(upstream_ins_type);

        memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(upstream_seq_type, Three_prime);

        // Map to old variables for compatibility if needed
        memory_layer_safety_1 = memory_layer_safety_upstream;

        downstream_proba_map.request_memory_layer(upstream_ins_type);
        memory_layer_proba_map_junction_upstream = downstream_proba_map.get_current_memory_layer(upstream_ins_type);
    }

    if (downstream_exists) {
        safety_set.request_memory_layer(downstream_ins_type);
        memory_layer_safety_downstream = safety_set.get_current_memory_layer(downstream_ins_type);

        downstream_proba_map.request_memory_layer(downstream_ins_type);
        memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(downstream_ins_type);

        memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(downstream_seq_type, Five_prime);

        // Map to old variables
        memory_layer_safety_2 = memory_layer_safety_downstream;
    }

    // Deletions (Generic lookup)
    v_3_min_del = 0;
    v_3_max_del = 0;
    if (upstream_exists) {
        if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, upstream_seq_type, Three_prime))) {
            auto del_p = events_map.at(tuple<Event_type, int, Seq_side>(Deletion_t, upstream_seq_type, Three_prime));
            if (processed_events.count(del_p->get_name()) == 0) {
                v_3_min_del = del_p->get_len_max();
                v_3_max_del = del_p->get_len_min();
            }
        }
    }

    j_5_min_del = 0;
    j_5_max_del = 0;
    if (downstream_exists) {
        if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, downstream_seq_type, Five_prime))) {
            auto del_p = events_map.at(tuple<Event_type, int, Seq_side>(Deletion_t, downstream_seq_type, Five_prime));
            if (processed_events.count(del_p->get_name()) == 0) {
                j_5_min_del = del_p->get_len_max();
                j_5_max_del = del_p->get_len_min();
            }
        }
    }

    d_5_min_del = 0;
    d_5_max_del = 0;
    if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id, Five_prime))) {
        auto del_p = events_map.at(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id, Five_prime));
        if (processed_events.count(del_p->get_name()) == 0) {
            d_5_min_del = del_p->get_len_max();
            d_5_max_del = del_p->get_len_min();
        }
    }

    d_3_min_del = 0;
    d_3_max_del = 0;
    if (events_map.count(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id, Three_prime))) {
        auto del_p = events_map.at(tuple<Event_type, int, Seq_side>(Deletion_t, this->sequence_type_id, Three_prime));
        if (processed_events.count(del_p->get_name()) == 0) {
            d_3_min_del = del_p->get_len_max();
            d_3_max_del = del_p->get_len_min();
        }
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list, seq_offsets,
                                      index_map);
}

void Gene_choice::set_nickname(string name)
{
    this->nickname = name;
    // For backward compatibility, sequence_type_id is set in constructor
}

void Gene_choice::update_event_internal_probas(const Marginal_array_p &marginal_array,
                                               const unordered_map<Rec_Event_name, int> &index_map)
{
}

/**
 * All add_to_marginals should take into account the possibility to perform
 * viterbi runs(take only the most likely scenario into account)
 */
void Gene_choice::add_to_marginals(long double scenario_proba, Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        updated_marginals[this->new_index] = scenario_proba;
    } else {
        updated_marginals[this->new_index] += scenario_proba;
    }
}

bool Gene_choice::has_effect_on(int seq_type) const
{
    for (const auto &junction : active_upstream_junctions) {
        if (junction.second == seq_type)
            return true;
    }
    for (const auto &junction : active_downstream_junctions) {
        if (junction.second == seq_type)
            return true;
    }
    return false;
}

void Gene_choice::iterate_initialize_Len_proba(int considered_junction, std::map<int, double> &length_best_proba_map,
                                               std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                               double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                               Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                               int &seq_len /*=0*/) const
{

    if (this->has_effect_on(considered_junction)) {
        base_index = base_index_map.at(this->event_index, 0);
        for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {

            // Get the max proba for this realization (in case the event is child of
            // another)
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()] > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index + i * this->size()];
                }
            }
            // Update the length within and probability in the recursive call
            Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map, model_queue,
                                                            scenario_proba * real_max_proba, model_parameters_point,
                                                            base_index_map, constructed_sequences,
                                                            seq_len + (*iter).second.value_str.length());
        }
    } else {
        Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map, model_queue,
                                                        scenario_proba, model_parameters_point, base_index_map,
                                                        constructed_sequences, seq_len);
    }
}

void Gene_choice::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                             const Marginal_array_p &model_parameters_point, Index_map &base_index_map)
{

    Seq_type_str_p_map constructed_sequences(SequenceTypeRegistry::get_instance().size());
    junction_length_best_proba_maps.clear();
    vj_length_d_position_proba.clear();

    // Iterate over active downstream junctions
    for (const auto &junction : active_downstream_junctions) {
        int ins_type = junction.second;
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(ins_type, junction_length_best_proba_maps[ins_type], model_queue,
                                                      init_proba, model_parameters_point, base_index_map,
                                                      constructed_sequences);
    }

    // Iterate over active upstream junctions
    for (const auto &junction : active_upstream_junctions) {
        int ins_type = junction.second;
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(ins_type, junction_length_best_proba_maps[ins_type], model_queue,
                                                      init_proba, model_parameters_point, base_index_map,
                                                      constructed_sequences);
    }
}
