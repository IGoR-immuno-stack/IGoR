/*
 * Genechoice.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Quentin Marcou
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

#include <igor/Core/EventUtils.h>
#include <igor/Core/Genechoice.h>

using namespace std;

namespace {

using LegacyEventsMap = unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>>;
using SeqEventsMap = unordered_map<tuple<Event_type, Seq_type, Seq_side>, shared_ptr<Rec_Event>>;

const LegacyEventsMap &empty_legacy_events_map()
{
    static const LegacyEventsMap kEmptyLegacyEventsMap;
    return kEmptyLegacyEventsMap;
}

SeqEventsMap build_seq_events_map(const LegacyEventsMap &events_map)
{
    SeqEventsMap seq_events_map;
    seq_events_map.reserve(events_map.size());
    for (const auto &entry : events_map) {
        tuple<Event_type, Seq_type, Seq_side> seq_key;
        if (!EventUtils::try_event_key_to_seq_key(get<0>(entry.first), get<1>(entry.first), get<2>(entry.first),
                                                   seq_key)) {
            continue;
        }
        seq_events_map.emplace(seq_key, entry.second);
    }
    return seq_events_map;
}

EventUtils::GeneChoiceStatus check_gene_choice_seq_type(
        Gene_class gene, const SeqEventsMap &events_map, const unordered_set<Rec_Event_name> &processed_events)
{
    EventUtils::GeneChoiceStatus status{ false, false, nullptr };
    Seq_type gene_seq = V_gene_seq;
    if (!EventUtils::try_gene_class_to_gene_seq_type(gene, gene_seq)) {
        return status;
    }

    shared_ptr<Rec_Event> event_ptr;
    if (!EventUtils::try_get_event(events_map, GeneChoice_t, gene_seq, Undefined_side, event_ptr)) {
        return status;
    }

    status.exists = true;
    status.chosen = processed_events.count(event_ptr->get_name()) != 0;
    status.event_ptr = event_ptr;
    return status;
}

} // namespace

Gene_choice::Gene_choice() : Gene_choice(Undefined_gene)
{
    this->type = Event_type::GeneChoice_t;
    this->update_event_name();
}

Gene_choice::Gene_choice(Gene_class gene)
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
    this->update_event_name();
}

/*
 * Should probably be avoided, default initialize the event and use add_event_realization() instead
 * unless you're sure about the index fields in the Event_realizations instances
 */
Gene_choice::Gene_choice(Gene_class gene, unordered_map<string, Event_realization> &realizations) : Gene_choice(gene)
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

Gene_choice::Gene_choice(Gene_class gene, vector<pair<string, string>> genomic_sequences) : Gene_choice(gene)
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
    new_gene_choice_p->update_event_name();
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
    //First remove previous realizations
    this->event_realizations.clear();
    for (vector<pair<string, string>>::const_iterator iter = genomic_templates.begin(); iter != genomic_templates.end();
         ++iter) {
        this->add_realization((*iter).first, (*iter).second);
    }
}

/**
 * @brief Context-based iterate() implementation
 *
 * Unpacks 5 context objects into legacy parameters and delegates
 * to the existing iterate() implementation.
 *
 */
void Gene_choice::iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    base_index = exploration.index_map.at(this->event_index);
    const double base_scenario_proba = scenario.scenario_proba;

    switch (this->event_class) {

        //TODO take into account in-dels and construct them in the constructed sequences
        //TODO withdraw the assignment body to iterate  so that other same kind of functions can be constructed
        //such as one to iterate one time and generate ""counters""

    case V_gene: {

        //Check D choice
        if (d_chosen) {
            //If D chosen need to check V safety
            d_offset = scenario.get_offset(D_gene_seq, Five_prime, memory_layer_offset_check1);
            d_5_min_offset = d_offset - d_5_min_del;
            d_5_max_offset = d_offset - d_5_max_del;

            vd_check = true; //Further check needed
        } else {
            vd_check = false;
            if (d_choice_exist) {
                exploration.set_overlap_safety(Event_safety::VD_safe, false, memory_layer_safety_1);

            } else {
                //If no D choice V choice is safe
                exploration.set_overlap_safety(Event_safety::VD_safe, true, memory_layer_safety_1);
            }
        }

        //Check J choice
        if (j_chosen) {
            //If J chosen need to check V safety
            j_offset = scenario.get_offset(J_gene_seq, Five_prime, memory_layer_offset_check2);
            j_5_min_offset = j_offset - j_5_min_del;
            j_5_max_offset = j_offset - j_5_max_del;

            vj_check = true; //Further check needed
        } else {
            vj_check = false;
            if (j_choice_exist) {
                exploration.set_overlap_safety(Event_safety::VJ_safe, false, memory_layer_safety_2);
            } else {
                //If no J choice V choice is safe
                exploration.set_overlap_safety(Event_safety::VJ_safe, true, memory_layer_safety_2);
            }
        }

        //Iterate over possible realizations (alignments provided for the V gene)
        for (vector<Alignment_data>::const_iterator iter = query.gene_alignments.at(V_gene).begin();
             iter != query.gene_alignments.at(V_gene).end(); ++iter) {

            if ((*iter).offset >= 0) {
                //gene_seq = this->event_realizations.at((*iter).gene_name).value_str ;
                //Use integer sequence (allow indexing on nucleotide identity)
                gene_seq = this->event_realizations.at((*iter).gene_name).value_str_int;
                v_5_off = (*iter).offset;
            } else {
                //If the offset is negative then the whole V-gene is not visible in the sequence thus only the aligned part of the gene is used.
                //gene_seq = this->event_realizations.at((*iter).gene_name).value_str.substr( -(*iter).offset ) ;
                //Use integer sequence (allow indexing on nucleotide identity)
                gene_seq = this->event_realizations.at((*iter).gene_name).value_str_int.substr(-(*iter).offset);
                v_5_off = 0;
            }
            //Insert the gene sequence as the constructed V gene sequence
            scenario.set_sequence_segment(V_gene_seq, &gene_seq, memory_layer_cs);

            //Compute v_3_offset
            v_3_off = v_5_off + gene_seq.size() - 1;

            //Check VD if needed
            if (vd_check) {

                if ((v_3_off + v_3_max_del) >= (d_5_max_offset)) {
                    //Even with maximum number of deletions on each side the V and D overlap => bad alignments
                    continue;
                }
                if ((v_3_off + v_3_min_del) < (d_5_min_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::VD_safe, true, memory_layer_safety_1);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::VD_safe, false, memory_layer_safety_1);
                }
            }

            //Check VJ if needed
            if (vj_check) {

                if ((v_3_off + v_3_max_del) >= (j_5_max_offset)) {
                    //Even with maximum number of deletions on each side the V and J overlap => bad alignments
                    continue;
                }
                if ((v_3_off + v_3_min_del) < (j_5_min_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::VJ_safe, true, memory_layer_safety_2);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::VJ_safe, false, memory_layer_safety_2);
                }
            }

            //Compute gene choice realization index
            current_realizations_index_vec[0] = this->event_realizations.at((*iter).gene_name).index;
            new_index = base_index + current_realizations_index_vec[0];
            new_scenario_proba = base_scenario_proba;
            //new_tmp_err_w_proba = tmp_err_w_proba;
            proba_contribution = 1;

            //State pointers
            current_realization_index = &this->event_realizations.at((*iter).gene_name).index;
            alignment_offset_p = &(*iter).offset;

            proba_contribution = iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                                exploration.index_map, model.offset_map, model.model_parameters);

            //Update scenario probability
            new_scenario_proba *= proba_contribution;
            new_tmp_err_w_proba *= proba_contribution;

            //Set seq offsets
            scenario.set_offset(V_gene_seq, Five_prime, v_5_off, memory_layer_off_fivep);
            scenario.set_offset(V_gene_seq, Three_prime, v_3_off, memory_layer_off_threep);

            //Set the the V mismatch list using the mismatch list computed during the alignment
            scenario.set_mismatches(V_gene_seq, &(*iter).mismatches, memory_layer_mismatches);


            //Get VD or VJ junction upper bound proba
            if (d_chosen) {
                if (vd_length_best_proba_map.count(d_offset - v_3_off - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(VD_ins_seq, vd_length_best_proba_map.at(d_offset - v_3_off - 1),
                                               memory_layer_proba_map_junction);
            } else if (j_chosen) {
                if (vj_length_best_proba_map.count(j_offset - v_3_off - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(VJ_ins_seq, vj_length_best_proba_map.at(j_offset - v_3_off - 1),
                                               memory_layer_proba_map_junction);
            }

            //Count the number of mismatches that will not go away whatever the number off deletions
            endogeneous_mismatches = 0;
            mism_iter = iter->mismatches.begin();
            while ((mism_iter != iter->mismatches.end()) and ((*mism_iter) <= v_3_off + v_3_max_del)) {
                //Count one mismatch
                ++endogeneous_mismatches;
                ++mism_iter;
            }
            exploration.downstream_proba_map.set_value(
                    V_gene_seq,
                    accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                           gene_seq.size() - v_3_max_del - endogeneous_mismatches),
                    memory_layer_proba_map_seq);

            //Multiply all downstream probas
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );

            if (exploration.should_prune(scenario_upper_bound_proba)) {
                continue;
            }

            // Update context with new probability before proceeding
            scenario.scenario_proba = new_scenario_proba;

            Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
        }
    } break;

    case D_gene: {

        //Check V choice
        if (v_chosen) {
            //If V chosen need to check D safety
            v_offset = scenario.get_offset(V_gene_seq, Three_prime, memory_layer_offset_check1);
            v_3_max_offset = v_offset + v_3_min_del;
            v_3_min_offset = v_offset + v_3_max_del;

            vd_check = true; //Further check needed

        } else {
            vd_check = false;
            if (v_choice_exist) {
                exploration.set_overlap_safety(Event_safety::VD_safe, false, memory_layer_safety_1);
            } else {
                //If no V choice D choice is safe
                exploration.set_overlap_safety(Event_safety::VD_safe, true, memory_layer_safety_1);
            }
        }

        //Check J choice
        if (j_chosen) {
            //If J chosen need to check D safety
            j_offset = scenario.get_offset(J_gene_seq, Five_prime, memory_layer_offset_check2);
            j_5_min_offset = j_offset - j_5_min_del;
            j_5_max_offset = j_offset - j_5_max_del;

            dj_check = true; //Further check needed
        } else {
            dj_check = false;

            //Useful in case of no D, however for speed purpose it might be better to process J choice first
            j_5_min_offset = query.sequence.size() - 1;
            j_5_max_offset = j_5_min_offset;

            if (j_choice_exist) {
                //safety_set.emplace(Event_safety::DJ_unsafe);
                exploration.set_overlap_safety(Event_safety::DJ_safe, false, memory_layer_safety_2);
            } else {
                //If no J choice V choice is safe
                //safety_set.emplace(Event_safety::DJ_safe);
                exploration.set_overlap_safety(Event_safety::DJ_safe, true, memory_layer_safety_2);
            }
        }

        no_d_align = true;

        //Iterate over possible realizations (alignments provided for the D gene)
        for (vector<Alignment_data>::const_iterator iter = query.gene_alignments.at(D_gene).begin();
             iter != query.gene_alignments.at(D_gene).end(); ++iter) {

            //gene_seq = this->event_realizations.at((*iter).gene_name).value_str;
            gene_seq = this->event_realizations.at((*iter).gene_name).value_str_int;

            scenario.set_sequence_segment(D_gene_seq, &gene_seq, memory_layer_cs);

            d_5_off = (*iter).offset;
            d_3_off = (*iter).offset + gene_seq.size() - 1;

            if (vd_check) {
                if ((d_5_off - d_5_max_del) <= (v_3_min_offset)) {
                    //Even with maximum number of deletions on each side the V and D overlap => bad alignments
                    continue;
                }
                if ((d_5_off - d_5_min_del) > (v_3_max_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::VD_safe, true, memory_layer_safety_1);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::VD_safe, false, memory_layer_safety_1);
                }
            }
            if (dj_check) {
                if ((d_3_off + d_3_max_del) >= (j_5_max_offset)) {
                    //Even with maximum number of deletions on each side the D and J overlap => bad alignments
                    continue;
                }
                if ((d_3_off + d_3_min_del) < (j_5_min_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::DJ_safe, true, memory_layer_safety_2);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::DJ_safe, false, memory_layer_safety_2);
                }
            }

            //FIXME deal with state pointers for D

            current_realizations_index_vec[0] = this->event_realizations.at((*iter).gene_name).index;
            new_index = base_index + current_realizations_index_vec[0];
            new_scenario_proba = base_scenario_proba;
            //new_tmp_err_w_proba = tmp_err_w_proba;
            proba_contribution = 1;

            proba_contribution = iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                                exploration.index_map, model.offset_map, model.model_parameters);

            new_scenario_proba *= proba_contribution;
            //new_tmp_err_w_proba*=proba_contribution;

            //Assume that the whole D is in the sequence and add the D sequence to the constructed sequences
            scenario.set_offset(D_gene_seq, Five_prime, (*iter).offset, memory_layer_off_fivep);
            scenario.set_offset(D_gene_seq, Three_prime, (*iter).offset + gene_seq.size() - 1,
                                  memory_layer_off_threep);

            scenario.set_mismatches(D_gene_seq, &(*iter).mismatches, memory_layer_mismatches);


            //Get DJ or VJ junction upper bound proba
            if (v_chosen and j_chosen) {
                if (vd_length_best_proba_map.count(d_5_off - v_offset - 1) <= 0
                    or dj_length_best_proba_map.count(j_offset - d_3_off - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(VJ_ins_seq, 1.0, memory_layer_proba_map_junction);
                exploration.downstream_proba_map.set_value(VD_ins_seq, vd_length_best_proba_map.at(d_5_off - v_offset - 1),
                                               memory_layer_proba_map_junction_d2);
                exploration.downstream_proba_map.set_value(DJ_ins_seq, dj_length_best_proba_map.at(j_offset - d_3_off - 1),
                                               memory_layer_proba_map_junction_d3);
            } else if (v_chosen) {
                if (vd_length_best_proba_map.count(d_5_off - v_offset - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(VD_ins_seq, vd_length_best_proba_map.at(d_5_off - v_offset - 1),
                                               memory_layer_proba_map_junction_d2);
            } else if (j_chosen) {
                if (dj_length_best_proba_map.count(j_offset - d_3_off - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(DJ_ins_seq, dj_length_best_proba_map.at(j_offset - d_3_off - 1),
                                               memory_layer_proba_map_junction_d3);
            }

            //Count the number of mismatches that will not go away even with maximum number of deletions
            endogeneous_mismatches = 0;
            if ((d_5_off - d_5_max_del) < (d_3_off + d_3_max_del)) {
                mism_iter = iter->mismatches.begin();
                while (mism_iter != iter->mismatches.end()) {
                    if ((*mism_iter) >= (d_5_off - d_5_max_del) and (*mism_iter) <= (d_3_off + d_3_max_del)) {
                        //Count one mismatch
                        ++endogeneous_mismatches;
                    }
                    ++mism_iter;
                }
                exploration.downstream_proba_map.set_value(D_gene_seq,
                                               accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                                                      (d_3_off + d_3_max_del)
                                                                                              - (d_5_off - d_5_max_del)
                                                                                              - endogeneous_mismatches),
                                               memory_layer_proba_map_seq);
            } else {
                exploration.downstream_proba_map.set_value(D_gene_seq, 1.0, memory_layer_proba_map_seq);
            }

            //Multiply all downstream probas
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );

            if (exploration.should_prune(scenario_upper_bound_proba)) {
                continue;
            }
            no_d_align = false;

            // Update context with new probability before proceeding
            scenario.scenario_proba = new_scenario_proba;

            Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
        }

        if (no_d_align) {
            //int test = 0;

            //Pass the mismatch vector pointer to the memory map once (will be updated in the next loop)
            scenario.set_mismatches(D_gene_seq, &no_d_mismatches, memory_layer_mismatches);

            if (v_chosen and j_chosen) {
                int vj_len = j_offset - v_offset - 1;
                if (vj_length_d_position_proba.count(vj_len) != 0) {
                    const vector<tuple<string, int, int, double>> &d_positions_vector =
                            vj_length_d_position_proba.at(vj_len);
                    for (vector<tuple<string, int, int, double>>::const_iterator d_position_iter =
                                 d_positions_vector.begin();
                         d_position_iter != d_positions_vector.end(); ++d_position_iter) {

                        const Event_realization &d_real = this->event_realizations.at(get<0>(*d_position_iter));

                        //d_5_off is v 3' offset + vd junction length
                        d_5_off = v_offset + get<1>(*d_position_iter);

                        if (d_5_off - d_5_min_del >= j_5_max_offset) {
                            continue;
                        }

                        d_size = d_real.value_str.size();

                        d_full_3_offset = d_5_off + d_size - 1;
                        d_3_max_offset = d_full_3_offset + d_3_min_del;

                        if (d_3_max_offset <= v_3_min_offset) {
                            continue;
                        }

                        gene_seq = d_real.value_str_int;
                        scenario.set_sequence_segment(D_gene_seq, &gene_seq, memory_layer_cs);

                        current_realizations_index_vec[0] = d_real.index;
                        new_index = base_index + current_realizations_index_vec[0];

                        //Proba contribution is the same wherever is the gene
                        proba_contribution = 1;
                        proba_contribution =
                                iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                               exploration.index_map, model.offset_map, model.model_parameters);

                        new_scenario_proba = scenario.scenario_proba * proba_contribution;


                        //Get DJ or VJ junction upper bound proba
                        /*							if(v_chosen and j_chosen){
									if(vd_length_best_proba_map.count(d_5_off - v_offset -1)<=0 or dj_length_best_proba_map.count(j_offset - d_full_3_offset  -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}*/
                        exploration.downstream_proba_map.set_value(VJ_ins_seq, 1.0, memory_layer_proba_map_junction);
                        exploration.downstream_proba_map.set_value(VD_ins_seq,
                                                       vd_length_best_proba_map.at(get<1>(*d_position_iter)),
                                                       memory_layer_proba_map_junction_d2);
                        exploration.downstream_proba_map.set_value(DJ_ins_seq,
                                                       dj_length_best_proba_map.at(get<2>(*d_position_iter)),
                                                       memory_layer_proba_map_junction_d3);
                        exploration.downstream_proba_map.set_value(D_gene_seq, 1.0,
                                                       memory_layer_proba_map_seq); //Lift the penalty on D gene seq

                        /*							}
								else if(v_chosen){
									if(vd_length_best_proba_map.count(d_5_off - v_offset -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}
									downstream_proba_map.set_value(VD_ins_seq , vd_length_best_proba_map.at(d_5_off - v_offset -1) , memory_layer_proba_map_junction_d2);
								}
								else if(j_chosen){
									if(dj_length_best_proba_map.count(j_offset - d_full_3_offset  -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}
									downstream_proba_map.set_value(DJ_ins_seq , dj_length_best_proba_map.at(j_offset - d_full_3_offset  -1) , memory_layer_proba_map_junction_d3);
								}*/

                        //Multiply all downstream probas
                        scenario_upper_bound_proba = exploration.compute_upper_bound(
                            new_scenario_proba,
                            current_downstream_proba_memory_layers
                        );

                        //If even without taking the weight of errors into account not good, then any lower one not good
                        if (exploration.should_prune(scenario_upper_bound_proba)) {
                            break;
                        }

                        //Get mismatches between D gene and sequence
                        no_d_mismatches.clear();
                        for (int i = 0; i != d_size; ++i) {
                            if (((d_5_off + i) >= 0) && (d_5_off + i) < query.int_sequence.size()) {
                                if (gene_seq[i] != query.int_sequence[d_5_off + i]) {
                                    no_d_mismatches.push_back(d_5_off + i);
                                }
                            }
                        }

                        //Count the number of mismatches that will not go away even with maximum number of deletions
                        endogeneous_mismatches = 0;
                        if ((d_5_off - d_5_max_del) < (d_full_3_offset + d_3_max_del)) {
                            mism_iter = no_d_mismatches.begin();
                            while (mism_iter != no_d_mismatches.end()) {
                                if ((*mism_iter) >= (d_5_off - d_5_max_del)
                                    and (*mism_iter) <= (d_full_3_offset + d_3_max_del)) {
                                    //Count one mismatch
                                    ++endogeneous_mismatches;
                                }
                                ++mism_iter;
                            }
                            //Weigh D_gene_seq accordingly
                            exploration.downstream_proba_map.set_value(
                                    D_gene_seq,
                                    accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                                           (d_full_3_offset + d_3_max_del)
                                                                                   - (d_5_off - d_5_max_del)
                                                                                   - endogeneous_mismatches),
                                    memory_layer_proba_map_seq);

                            //Multiply all downstream probas
                            scenario_upper_bound_proba = exploration.compute_upper_bound(
                                new_scenario_proba,
                                current_downstream_proba_memory_layers
                            );

                            if (exploration.should_prune(scenario_upper_bound_proba)) {
                                continue;
                            }

                        } else {
                            exploration.downstream_proba_map.set_value(D_gene_seq, 1.0, memory_layer_proba_map_seq);
                        }

                        //Assume that the whole D is in the sequence and add the D sequence to the constructed sequences
                        scenario.set_offset(D_gene_seq, Five_prime, d_5_off, memory_layer_off_fivep);
                        scenario.set_offset(D_gene_seq, Three_prime, d_5_off + d_size - 1, memory_layer_off_threep);

                        // Update context with new probability before proceeding
                        scenario.scenario_proba = new_scenario_proba;

                        Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
                    }
                }
            } else {
                for (unordered_map<string, Event_realization>::const_iterator d_gene_iter =
                             this->event_realizations.begin();
                     d_gene_iter != this->event_realizations.end(); ++d_gene_iter) {

                    //Starts the D one nucleotide after v_3_min offset(-1 if no V chosen) given max deletions on the 5' of the D
                    //FIXME v_min offset set to -1 if no V chosen
                    d_size = (*d_gene_iter).second.value_str.size();

                    //Take care of the fact that not all D have the same length
                    // and that the maximum number of deletions might be greater than the D itself
                    if ((-d_5_max_del) > d_size) {
                        d_5_real_max_del = -d_size;
                    } else {
                        d_5_real_max_del = d_5_max_del;
                    }

                    if (v_3_min_offset > 0) {
                        d_5_off = v_3_min_offset + d_5_real_max_del + 1;
                    } else {
                        d_5_off = 1 + d_5_real_max_del
                                + 1; //Consider that V cannot be absent from the read, at least one nucleotide is present
                    }

                    d_full_3_offset = d_5_off + d_size - 1;
                    d_3_max_offset = d_full_3_offset + d_3_min_del; //Useless?
                    if (abs(d_3_max_del) < d_size) {
                        d_3_min_offset = d_full_3_offset + d_3_max_del;
                    } else {
                        d_3_min_offset = d_5_off;
                    }

                    //Always the same sequence for the given D
                    gene_seq = (*d_gene_iter).second.value_str_int;
                    scenario.set_sequence_segment(D_gene_seq, &gene_seq, memory_layer_cs);

                    current_realizations_index_vec[0] = d_gene_iter->second.index;
                    new_index = base_index + current_realizations_index_vec[0];

                    //Proba contribution is the same wherever is the gene
                    proba_contribution = 1;
                    proba_contribution =
                            iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                           exploration.index_map, model.offset_map, model.model_parameters);
                    //new_tmp_err_w_proba = tmp_err_w_proba*proba_contribution;
                    /*					compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
						if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
							continue;
						}*/

                    while (d_3_min_offset < j_5_min_offset) {
                        //Slides the D one nucleotide at a time towards 3', updating the mismatch list,offsets

                        //Get mismatches between D gene and sequence at the 5' most position
                        no_d_mismatches.clear();
                        for (int i = 0; i != d_size; ++i) {
                            if (((d_5_off + i) >= 0) && (d_5_off + i) < query.int_sequence.size()) {
                                if (gene_seq[i] != query.int_sequence[d_5_off + i]) {
                                    no_d_mismatches.push_back(d_5_off + i);
                                }
                            }
                        }

                        new_scenario_proba = scenario.scenario_proba * proba_contribution;
                        //new_tmp_err_w_proba = tmp_err_w_proba*proba_contribution;

                        /*if( (d_full_3_offset<0)){
								cout<<"problem in gene choice"<<endl;
								cout<<d_full_3_offset<<endl;
								cout<<v_3_min_offset<<endl;
								cout<<d_5_max_del<<endl;
							}*/

                        //Assume that the whole D is in the sequence and add the D sequence to the constructed sequences
                        scenario.set_offset(D_gene_seq, Five_prime, d_5_off, memory_layer_off_fivep);
                        scenario.set_offset(D_gene_seq, Three_prime, d_full_3_offset, memory_layer_off_threep);


                        //Get DJ or VJ junction upper bound proba
                        if (v_chosen and j_chosen) {
                            if (vd_length_best_proba_map.count(d_5_off - v_offset - 1) <= 0
                                or dj_length_best_proba_map.count(j_offset - d_full_3_offset - 1) <= 0) {
                                continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                            }
                            exploration.downstream_proba_map.set_value(VJ_ins_seq, 1.0, memory_layer_proba_map_junction);
                            exploration.downstream_proba_map.set_value(VD_ins_seq,
                                                           vd_length_best_proba_map.at(d_5_off - v_offset - 1),
                                                           memory_layer_proba_map_junction_d2);
                            exploration.downstream_proba_map.set_value(DJ_ins_seq,
                                                           dj_length_best_proba_map.at(j_offset - d_full_3_offset - 1),
                                                           memory_layer_proba_map_junction_d3);
                        } else if (v_chosen) {
                            if (vd_length_best_proba_map.count(d_5_off - v_offset - 1) <= 0) {
                                continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                            }
                            exploration.downstream_proba_map.set_value(VD_ins_seq,
                                                           vd_length_best_proba_map.at(d_5_off - v_offset - 1),
                                                           memory_layer_proba_map_junction_d2);
                        } else if (j_chosen) {
                            if (dj_length_best_proba_map.count(j_offset - d_full_3_offset - 1) <= 0) {
                                continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                            }
                            exploration.downstream_proba_map.set_value(DJ_ins_seq,
                                                           dj_length_best_proba_map.at(j_offset - d_full_3_offset - 1),
                                                           memory_layer_proba_map_junction_d3);
                        }

                        //Count the number of mismatches that will not go away even with maximum number of deletions
                        endogeneous_mismatches = 0;
                        if ((d_5_off - d_5_max_del) < (d_full_3_offset + d_3_max_del)) {
                            mism_iter = no_d_mismatches.begin();
                            while (mism_iter != no_d_mismatches.end()) {
                                if ((*mism_iter) >= (d_5_off - d_5_max_del)
                                    and (*mism_iter) <= (d_full_3_offset + d_3_max_del)) {
                                    //Count one mismatch
                                    ++endogeneous_mismatches;
                                }
                                ++mism_iter;
                            }
                            exploration.downstream_proba_map.set_value(
                                    D_gene_seq,
                                    accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                                           (d_full_3_offset + d_3_max_del)
                                                                                   - (d_5_off - d_5_max_del)
                                                                                   - endogeneous_mismatches),
                                    memory_layer_proba_map_seq);
                        } else {
                            exploration.downstream_proba_map.set_value(D_gene_seq, 1.0, memory_layer_proba_map_seq);
                        }

                        //Multiply all downstream probas
                        scenario_upper_bound_proba = exploration.compute_upper_bound(
                            new_scenario_proba,
                            current_downstream_proba_memory_layers
                        );

                        if (exploration.should_prune(scenario_upper_bound_proba)) {
                            continue;
                        }

                        // Update context with new probability before proceeding
                        scenario.scenario_proba = new_scenario_proba;

                        Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);

                        //test++;

                        //Slide the D from 1 nucleotide
                        ++d_5_off;
                        ++d_full_3_offset;
                        ++d_3_min_offset;
                        ++d_3_max_offset;

                        /*//Adapt D mismatches if needed
							if(!no_d_mismacthes.empty()){
								if(no_d_mismacthes[0]<d_5_off) {
									no_d_mismacthes.erase(no_d_mismacthes.begin());
								}
							}
							if(gene_seq[d_size-1] != int_sequence[d_full_3_offset]) no_d_mismacthes.push_back(d_full_3_offset);
	*/
                    }
                }
                //cout<<"Seq "<<sequence<<"; #Ds made up" <<test<<endl;
            }
        }

    } break;

    case J_gene: {

        //Check D choice
        if (d_chosen) {
            //If D chosen need to check J safety
            d_offset = scenario.get_offset(D_gene_seq, Three_prime, memory_layer_offset_check2);
            d_3_min_offset = d_offset + d_3_max_del;
            d_3_max_offset = d_offset + d_3_min_del;

            dj_check = true; //Further check needed
        } else {
            dj_check = false;
            if (d_choice_exist) {
                exploration.set_overlap_safety(Event_safety::DJ_safe, false, memory_layer_safety_2);
            } else {
                //If no D choice V choice is safe
                exploration.set_overlap_safety(Event_safety::DJ_safe, true, memory_layer_safety_2);
            }
        }

        //Check V choice
        if (v_chosen) {
            //If V chosen need to check J safety
            v_offset = scenario.get_offset(V_gene_seq, Three_prime, memory_layer_offset_check1);
            v_3_min_offset = v_offset + v_3_max_del;
            v_3_max_offset = v_offset + v_3_min_del;

            vj_check = true; //Further check needed
        } else {
            vj_check = false;
            if (v_choice_exist) {
                exploration.set_overlap_safety(Event_safety::VJ_safe, false, memory_layer_safety_1);
            } else {
                //If no V choice J choice is safe
                exploration.set_overlap_safety(Event_safety::VJ_safe, true, memory_layer_safety_1);
            }
        }

        //Iterate over possible realizations (J gene alignments)
        for (vector<Alignment_data>::const_iterator iter = query.gene_alignments.at(J_gene).begin();
             iter != query.gene_alignments.at(J_gene).end(); ++iter) {

            j_5_off = (*iter).offset;

            if (vj_check) {

                if ((j_5_off - j_5_max_del) <= (v_3_min_offset)) {
                    //Even with maximum number of deletions on each side the V and D overlap => bad alignments
                    continue;
                }
                if ((j_5_off - j_5_min_del) > (v_3_max_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::VJ_safe, true, memory_layer_safety_1);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::VJ_safe, false, memory_layer_safety_1);
                }
            }
            if (dj_check) {
                if ((j_5_off - j_5_max_del) <= (d_3_min_offset)) {
                    //Even with maximum number of deletions on each side the D and J overlap => bad alignments
                    continue;
                }
                if ((j_5_off - j_5_min_del) > (d_3_max_offset)) {
                    //Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
                    exploration.set_overlap_safety(Event_safety::DJ_safe, true, memory_layer_safety_2);
                } else {
                    //In the deletion range => some number of deletion won't be allowed and will be discarded in the deletion process
                    exploration.set_overlap_safety(Event_safety::DJ_safe, false, memory_layer_safety_2);
                }
            }

            current_realizations_index_vec[0] = this->event_realizations.at((*iter).gene_name).index;
            new_index = base_index + current_realizations_index_vec[0];
            new_scenario_proba = base_scenario_proba;
            //new_tmp_err_w_proba = tmp_err_w_proba;
            proba_contribution = 1;

            //State pointers
            current_realization_index = &this->event_realizations.at((*iter).gene_name).index;
            alignment_offset_p = &(*iter).offset;

            proba_contribution = iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                                exploration.index_map, model.offset_map, model.model_parameters);

            new_scenario_proba *= proba_contribution;
            new_tmp_err_w_proba *= proba_contribution;

            //Compute the number of nucleotides at the end of the sequence that are not aligned with the J-gene and remove them
            //gene_seq = this->event_realizations.at((*iter).gene_name).value_str.substr(0,sequence.size() - (*iter).offset);
            gene_seq = this->event_realizations.at((*iter).gene_name)
                               .value_str_int.substr(0, query.sequence.size() - (*iter).offset);

            scenario.set_sequence_segment(J_gene_seq, &gene_seq, memory_layer_cs);

            scenario.set_offset(J_gene_seq, Five_prime, (*iter).offset, memory_layer_off_fivep);
            scenario.set_offset(J_gene_seq, Three_prime, (*iter).offset + gene_seq.size() - 1,
                                  memory_layer_off_threep);

            //Mismatches list computed during alignment
            scenario.set_mismatches(J_gene_seq, &(*iter).mismatches, memory_layer_mismatches);


            //Get DJ or VJ junction upper bound proba
            if (d_chosen) {
                if (dj_length_best_proba_map.count(j_5_off - d_offset - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(DJ_ins_seq, dj_length_best_proba_map.at(j_5_off - d_offset - 1),
                                               memory_layer_proba_map_junction);
            } else if (v_chosen) {
                if (vj_length_best_proba_map.count(j_5_off - v_offset - 1) <= 0) {
                    continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                }
                exploration.downstream_proba_map.set_value(VJ_ins_seq, vj_length_best_proba_map.at(j_5_off - v_offset - 1),
                                               memory_layer_proba_map_junction);
            }

            //Count the number of mismatches that will not go away even with maximum number of deletions
            endogeneous_mismatches = 0;
            rev_mism_iter = iter->mismatches.rbegin();
            while ((rev_mism_iter != iter->mismatches.rend()) and ((*rev_mism_iter) >= j_5_off - j_5_max_del)) {
                //Count one mismatch
                ++endogeneous_mismatches;
                ++rev_mism_iter;
            }
            exploration.downstream_proba_map.set_value(
                    J_gene_seq,
                    accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                           gene_seq.size() - j_5_max_del - endogeneous_mismatches),
                    memory_layer_proba_map_seq);

            //Multiply all downstream probas
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );

            if (exploration.should_prune(scenario_upper_bound_proba)) {
                continue;
            }

            // Update context with new probability before proceeding
            scenario.scenario_proba = new_scenario_proba;

            Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
        }
    } break;

    default:
        throw invalid_argument(std::string("Unknown gene_class for GeneChoice: ") + this->event_class);
        break;
    }
}

/*
 *This short method performs the iterate operations common to all Rec_event (modify index map and fetch realization probability)
 *
 */
double Gene_choice::iterate_common(
        double scenario_proba, const int &gene_index, int base_index, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        const Marginal_array_p &model_parameters)
{
    return scenario_proba * Rec_Event::iterate_common(gene_index, base_index,
                                                      base_index_map,
                                                      model_parameters);
}

queue<int> Gene_choice::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        unordered_map<Seq_type, string> &constructed_sequences, mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;

    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            switch (this->event_class) {
            case V_gene:
                constructed_sequences[V_gene_seq] = (*iter).second.value_str;
                break;
            case D_gene:
                constructed_sequences[D_gene_seq] = (*iter).second.value_str;
                break;
            case J_gene:
                constructed_sequences[J_gene_seq] = (*iter).second.value_str;
                break;
            default:
                break;
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
void Gene_choice::write2txt(ofstream &outfile)
{
    outfile << "#GeneChoice;" << event_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.name << ";" << (*iter).second.value_str << ";" << (*iter).second.index << endl;
    }
}

void Gene_choice::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
    const unordered_map<tuple<Event_type, Seq_type, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{

    //Check V choice
    auto v_status = check_gene_choice_seq_type(V_gene, events_map, processed_events);
    v_choice_exist = v_status.exists;
    v_chosen = v_status.chosen;

    //Check D choice
    auto d_status = check_gene_choice_seq_type(D_gene, events_map, processed_events);
    d_choice_exist = d_status.exists;
    d_chosen = d_status.chosen;

    //Check J choice
    auto j_status = check_gene_choice_seq_type(J_gene, events_map, processed_events);
    j_choice_exist = j_status.exists;
    j_chosen = j_status.chosen;

    switch (this->event_class) {
    case V_gene:
        seq_offsets.request_memory_layer(V_gene_seq, Three_prime);
        this->memory_layer_off_threep = seq_offsets.get_current_memory_layer(V_gene_seq, Three_prime);
        seq_offsets.request_memory_layer(V_gene_seq, Five_prime);
        this->memory_layer_off_fivep = seq_offsets.get_current_memory_layer(V_gene_seq, Five_prime);
        mismatches_list.request_memory_layer(V_gene_seq);
        this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(V_gene_seq);
        constructed_sequences.request_memory_layer(V_gene_seq);
        this->memory_layer_cs = constructed_sequences.get_current_memory_layer(V_gene_seq);
        //if(d_chosen){
        safety_set.request_memory_layer(VD_safe);
        memory_layer_safety_1 = safety_set.get_current_memory_layer(VD_safe);
        //cout<<"V_choice 1: "<<memory_layer_safety_1<<endl;
        //}
        //if(j_chosen){
        safety_set.request_memory_layer(VJ_safe);
        memory_layer_safety_2 = safety_set.get_current_memory_layer(VJ_safe);
        //cout<<"V_choice 2: "<<memory_layer_safety_2<<endl;
        //}

        downstream_proba_map.request_memory_layer(V_gene_seq);
        memory_layer_proba_map_seq = downstream_proba_map.get_current_memory_layer(V_gene_seq);
        if (d_chosen) {
            downstream_proba_map.request_memory_layer(VD_ins_seq);
            memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(VD_ins_seq);
        } else if (j_chosen) {
            downstream_proba_map.request_memory_layer(VJ_ins_seq);
            memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(VJ_ins_seq);
        }

        if (d_chosen) {
            memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(D_gene_seq, Five_prime);
        }
        if (j_chosen) {
            memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(J_gene_seq, Five_prime);
        }

        break;
    case D_gene:
        seq_offsets.request_memory_layer(D_gene_seq, Three_prime);
        this->memory_layer_off_threep = seq_offsets.get_current_memory_layer(D_gene_seq, Three_prime);
        seq_offsets.request_memory_layer(D_gene_seq, Five_prime);
        this->memory_layer_off_fivep = seq_offsets.get_current_memory_layer(D_gene_seq, Five_prime);
        mismatches_list.request_memory_layer(D_gene_seq);
        this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(D_gene_seq);
        constructed_sequences.request_memory_layer(D_gene_seq);
        this->memory_layer_cs = constructed_sequences.get_current_memory_layer(D_gene_seq);
        //if(v_chosen){
        safety_set.request_memory_layer(VD_safe);
        memory_layer_safety_1 = safety_set.get_current_memory_layer(VD_safe);
        //cout<<"D_choice 1: "<<memory_layer_safety_1<<endl;
        //}
        //if(j_chosen){
        safety_set.request_memory_layer(DJ_safe);
        memory_layer_safety_2 = safety_set.get_current_memory_layer(DJ_safe);
        //cout<<"D_choice 2: "<<memory_layer_safety_2<<endl;
        //}

        downstream_proba_map.request_memory_layer(D_gene_seq);
        memory_layer_proba_map_seq = downstream_proba_map.get_current_memory_layer(D_gene_seq);
        if (v_chosen) {
            downstream_proba_map.request_memory_layer(VD_ins_seq);
            memory_layer_proba_map_junction_d2 = downstream_proba_map.get_current_memory_layer(VD_ins_seq);
        }
        if (j_chosen) {
            downstream_proba_map.request_memory_layer(DJ_ins_seq);
            memory_layer_proba_map_junction_d3 = downstream_proba_map.get_current_memory_layer(DJ_ins_seq);
        }
        if (v_chosen and j_chosen) {
            downstream_proba_map.request_memory_layer(VJ_ins_seq);
            memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(VJ_ins_seq);
        }

        if (v_chosen) {
            memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(V_gene_seq, Three_prime);
        } else {
            v_3_min_offset = 0;
            v_3_max_offset = 0;
        }
        if (j_chosen) {
            memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(J_gene_seq, Five_prime);
        }

        break;
    case J_gene:
        seq_offsets.request_memory_layer(J_gene_seq, Three_prime);
        this->memory_layer_off_threep = seq_offsets.get_current_memory_layer(J_gene_seq, Three_prime);
        seq_offsets.request_memory_layer(J_gene_seq, Five_prime);
        this->memory_layer_off_fivep = seq_offsets.get_current_memory_layer(J_gene_seq, Five_prime);
        mismatches_list.request_memory_layer(J_gene_seq);
        this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(J_gene_seq);
        constructed_sequences.request_memory_layer(J_gene_seq);
        this->memory_layer_cs = constructed_sequences.get_current_memory_layer(J_gene_seq);
        //if(v_chosen){
        safety_set.request_memory_layer(VJ_safe);
        memory_layer_safety_1 = safety_set.get_current_memory_layer(VJ_safe);
        //cout<<"j_choice 1: "<<memory_layer_safety_1<<endl;
        //}
        //if(d_chosen){
        safety_set.request_memory_layer(DJ_safe);
        memory_layer_safety_2 = safety_set.get_current_memory_layer(DJ_safe);
        //cout<<"j_choice 2: "<<memory_layer_safety_2<<endl;
        //}

        downstream_proba_map.request_memory_layer(J_gene_seq);
        memory_layer_proba_map_seq = downstream_proba_map.get_current_memory_layer(J_gene_seq);
        if (d_chosen) {
            downstream_proba_map.request_memory_layer(DJ_ins_seq);
            memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(DJ_ins_seq);
        } else if (v_chosen) {
            downstream_proba_map.request_memory_layer(VJ_ins_seq);
            memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(VJ_ins_seq);
        }

        if (v_chosen) {
            memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(V_gene_seq, Three_prime);
        }
        if (d_chosen) {
            memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(D_gene_seq, Three_prime);
        }

        break;
    default:
        break;
    }

    //downstream_proba_map.get_all_current_memory_layer(current_downstream_proba_memory_layers);

    //Get V 3' deletion
    shared_ptr<Rec_Event> del_v_p;
    if (EventUtils::try_get_event(events_map, Deletion_t, V_gene_seq, Three_prime, del_v_p)) {
        if (processed_events.count(del_v_p->get_name()) != 0) {
            v_3_min_del = 0;
            v_3_max_del = 0;
        } else {
            v_3_min_del = del_v_p->get_len_max();
            v_3_max_del = del_v_p->get_len_min();
        }
    } else {
        v_3_min_del = 0;
        v_3_max_del = 0;
    }

    //Get D 5' deletion range
    shared_ptr<Rec_Event> del_d_p;
    if (EventUtils::try_get_event(events_map, Deletion_t, D_gene_seq, Five_prime, del_d_p)) {
        if (processed_events.count(del_d_p->get_name()) != 0) {
            d_5_min_del = 0;
            d_5_max_del = 0;
        } else {
            d_5_min_del = del_d_p->get_len_max();
            d_5_max_del = del_d_p->get_len_min();
        }
    } else {
        d_5_min_del = 0;
        d_5_max_del = 0;
    }

    //Get D 3' deletion
    shared_ptr<Rec_Event> del_d_3_p;
    if (EventUtils::try_get_event(events_map, Deletion_t, D_gene_seq, Three_prime, del_d_3_p)) {
        if (processed_events.count(del_d_3_p->get_name()) != 0) {
            d_3_min_del = 0;
            d_3_max_del = 0;
        } else {
            d_3_min_del = del_d_3_p->get_len_max();
            d_3_max_del = del_d_3_p->get_len_min();
        }
    } else {
        d_3_min_del = 0;
        d_3_max_del = 0;
    }

    //Get J 5' deletion range
    shared_ptr<Rec_Event> del_j_p;
    if (EventUtils::try_get_event(events_map, Deletion_t, J_gene_seq, Five_prime, del_j_p)) {
        if (processed_events.count(del_j_p->get_name()) != 0) {
            j_5_min_del = 0;
            j_5_max_del = 0;
        } else {
            j_5_min_del = del_j_p->get_len_max();
            j_5_max_del = del_j_p->get_len_min();
        }
    } else {
        j_5_min_del = 0;
        j_5_max_del = 0;
    }
    this->Rec_Event::initialize_event(processed_events, empty_legacy_events_map(), offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list, seq_offsets,
                                      index_map);
}

/**
 * All add_to_marginals should take into account the possibility to perform viterbi runs(take only the most likely scenario into account)
 */
void Gene_choice::add_to_marginals(long double scenario_proba, Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        updated_marginals[this->new_index] = scenario_proba;
    } else {
        updated_marginals[this->new_index] += scenario_proba;
    }
}

bool Gene_choice::has_effect_on(Seq_type seq_type) const
{
    switch (this->event_class) {
    case V_gene:
        return false;
        break;

    case D_gene:
        if (seq_type == VJ_ins_seq) {
            return true;
        } else {
            return false;
        }
        break;

    case J_gene:
        return false;
        break;
    default:
        return false;
    }
}

void Gene_choice::iterate_initialize_Len_proba(Seq_type considered_junction,
                                               std::map<int, double> &length_best_proba_map,
                                               std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                               double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                               Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                               int &seq_len /*=0*/) const
{

    if (this->has_effect_on(considered_junction)) {
        base_index = base_index_map.at(this->event_index, 0);
        for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {

            /*		//Update base index map
			for(forward_list<tuple<int,int,int>>::const_iterator jiter = memory_and_offsets.begin() ; jiter!=memory_and_offsets.end() ; ++jiter){
				//Get previous index for the considered event
				int previous_index = base_index_map.at(get<0>(*jiter),get<1>(*jiter)-1);
				//Update the index given the realization and the offset
				previous_index += iter->second.index *get<2>(*jiter);
				//Set the value
				base_index_map.set_value(get<0>(*jiter) , previous_index , get<1>(*jiter));
			}*/

            //Get the max proba for this realization (in case the event is child of another)
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()] > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index + i * this->size()];
                }
            }
            //Update the length within and probability in the recursive call
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

    Seq_type_str_p_map constructed_sequences(6);
    switch (this->event_class) {
    case V_gene:
        vd_length_best_proba_map.clear();
        vj_length_best_proba_map.clear();

        if (d_chosen) {
            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(VD_ins_seq, vd_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        } else if (j_chosen) {
            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(VJ_ins_seq, vj_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        }
        break;
    case D_gene:
        vd_length_best_proba_map.clear();
        dj_length_best_proba_map.clear();
        vj_length_d_position_proba.clear();

        if (v_chosen) {
            double init_proba = 1.0;
            constructed_sequences.reset();
            this->Rec_Event::iterate_initialize_Len_proba(VD_ins_seq, vd_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        }
        if (j_chosen) {
            double init_proba = 1.0;
            constructed_sequences.reset();
            this->Rec_Event::iterate_initialize_Len_proba(DJ_ins_seq, dj_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        }

        if (v_chosen and j_chosen) {
            int junction_len;

            //Loop over D gene choices
            for (unordered_map<string, Event_realization>::const_iterator d_gene_iter =
                         this->event_realizations.begin();
                 d_gene_iter != this->event_realizations.end(); ++d_gene_iter) {
                //Get considered D gene best proba
                double d_gene_max_proba = 0;
                base_index = base_index_map.at(this->event_index, 0);
                for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                    if (model_parameters_point[base_index + d_gene_iter->second.index + i * this->size()]
                        > d_gene_max_proba) {
                        d_gene_max_proba =
                                model_parameters_point[base_index + d_gene_iter->second.index + i * this->size()];
                    }
                }
                //Loop over possible VD junction lengths
                for (map<int, double>::const_iterator vd_len_iter = vd_length_best_proba_map.begin();
                     vd_len_iter != vd_length_best_proba_map.end(); ++vd_len_iter) {
                    //Loop over possible DJ junction lengths
                    for (map<int, double>::const_iterator dj_len_iter = dj_length_best_proba_map.begin();
                         dj_len_iter != dj_length_best_proba_map.end(); ++dj_len_iter) {
                        junction_len = d_gene_iter->second.value_str.size() + vd_len_iter->first + dj_len_iter->first;

                        if (vj_length_d_position_proba.count(junction_len) != 0) {
                            vj_length_d_position_proba.at(junction_len)
                                    .emplace_back(d_gene_iter->first, vd_len_iter->first, dj_len_iter->first,
                                                  (d_gene_max_proba * vd_len_iter->second * dj_len_iter->second));
                        } else {
                            vj_length_d_position_proba.emplace(
                                    piecewise_construct, make_tuple(junction_len),
                                    make_tuple(1,
                                               make_tuple(d_gene_iter->first, vd_len_iter->first, dj_len_iter->first,
                                                          (d_gene_max_proba * vd_len_iter->second
                                                           * dj_len_iter->second))));
                        }
                    }
                }
            }

            //Now sort each vector in the map in decreasing order of probability (according to the model)
            for (map<int, vector<tuple<string, int, int, double>>>::iterator d_position_map_iter =
                         vj_length_d_position_proba.begin();
                 d_position_map_iter != vj_length_d_position_proba.end(); ++d_position_map_iter) {
                sort(d_position_map_iter->second.begin(), d_position_map_iter->second.end(), D_position_tuple);
            }
        }

        break;
    case J_gene:
        dj_length_best_proba_map.clear();
        vj_length_best_proba_map.clear();

        if (d_chosen) {
            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(DJ_ins_seq, dj_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        } else if (v_chosen) {
            double init_proba = 1.0;
            this->Rec_Event::iterate_initialize_Len_proba(VJ_ins_seq, vj_length_best_proba_map, model_queue, init_proba,
                                                          model_parameters_point, base_index_map,
                                                          constructed_sequences);
        }
        break;
    default:
        break;
    }
}
