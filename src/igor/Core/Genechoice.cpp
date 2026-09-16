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
#include <igor/Core/gene_to_seqtype_migr.h>

using namespace std;



Gene_choice::Gene_choice() : Gene_choice(Undefined_gene)
{
    this->type = Event_type::GeneChoice_t;
    this->update_event_name();
}

Gene_choice::Gene_choice(Gene_class gene)
    : Rec_Event(gene, Undefined_side),
      d_5_off(INT16_MAX),
      d_3_min_offset(INT16_MAX),
      d_3_max_offset(INT16_MAX),
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
      memory_layer_off_threep(-1),
      memory_layer_off_fivep(-1),
      v_chosen(false),
      v_choice_exist(true),
      d_chosen(false),
      d_choice_exist(false),
      j_chosen(false),
      j_choice_exist(true),
      d_5_max_del(INT16_MIN),
      d_5_min_del(INT16_MAX),
      d_5_real_max_del(INT16_MIN),
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
    new_gene_choice_p->set_seq_type(this->get_seq_type());
    new_gene_choice_p->set_seq_type_id(this->get_seq_type_id());
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
/*
 * One body for V, D and J. Once initialize_event() has settled the topology, the three-way
 * switch over event_class this replaces differed in exactly four things:
 *
 *   - which neighbouring ends the placement is checked against -- flank_checks_;
 *   - whether the template may overhang the read, and at which end;
 *   - which junctions the placement bounds -- junction_bounds_, resolved in S4c;
 *   - how the surviving error-free core is credited -- section 7.1's defect, carried verbatim.
 *
 * None of those is a gene class, which is why a tandem D needs no fourth case here. `event_class`
 * stays, and stays the alignment-strategy key: query.gene_alignments is keyed by Gene_class, and
 * D1/D2 correctly share one alignment set (parent plan B0).
 */
void Gene_choice::iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    base_index = exploration.index_map.get(this->event_index);
    const double base_scenario_proba = scenario.scenario_proba;

    //Where each checked neighbour sits, read once per scenario. A neighbour that has not been
    //placed cannot be checked against, and the verdict recorded for it says only whether it will
    //ever need checking -- which is what the six preambles this replaces each said.
    for (FlankCheck &check : flank_checks_) {
        if (check.partner_chosen) {
            const Seq_Offset partner_offset =
                    scenario.get_offset(static_cast<Seq_type>(check.partner_id), check.partner_side,
                                        check.partner_offset_layer);
            neighbour_offset_[check.partner_id] = partner_offset;
            neighbour_reach_[check.partner_id] =
                    pending_.reachable(check.partner_id, check.partner_side, partner_offset);
            check.active = true;
        } else {
            check.active = false;
            exploration.set_overlap_safety(check.safety_slot, not check.partner_exists, check.safety_layer);
            if (exhaustive_position_fallback_ and check.this_is_left) {
                //Nothing placed to the right, so the exhaustive scan's right wall is the end of
                //the read. Read only by the no_d_align path below, which is 5b's.
                const Seq_Offset read_end = static_cast<Seq_Offset>(query.sequence.size()) - 1;
                neighbour_reach_[check.partner_id] = {read_end, read_end};
            }
        }
    }

    no_d_align = true;

    for (const Alignment_data &alignment : query.gene_alignments.at(this->event_class)) {
        const Event_realization &realization = this->event_realizations.at(alignment.gene_name);

        //Clip the template to the read. Only a gene at an end of the ordering can overhang it,
        //and the two ends clip on opposite sides -- an asymmetry T0 pins in two sections, and
        //which B3 turns into a flank segment of its own rather than a shortening of this one.
        const int clipped_prefix = clip_template_before_read_ ? std::max(0, -alignment.offset) : 0;
        my_5_off = alignment.offset + clipped_prefix;
        if (clip_template_after_read_) {
            gene_seq = realization.value_str_int.substr(clipped_prefix,
                                                        query.sequence.size() - alignment.offset);
        } else {
            gene_seq = realization.value_str_int.substr(clipped_prefix);
        }
        my_3_off = my_5_off + static_cast<Seq_Offset>(gene_seq.size()) - 1;

        scenario.set_sequence_segment(static_cast<Seq_type>(this->seq_type_id), &gene_seq, memory_layer_cs);

        //Can this placement still avoid colliding with each placed neighbour? One predicate
        //(section 2.2) for what was twelve hand-written comparisons, differing only in which
        //side of the pair this event sits on. The gap is 0 at every site today: every segment
        //that can sit between two checked ends can also be empty (section 7.2).
        bool infeasible = false;
        for (const FlankCheck &check : flank_checks_) {
            if (not check.active) {
                continue;
            }
            const JunctionGeometry::OffsetInterval mine =
                    pending_.reachable(this->seq_type_id, check.this_is_left ? Three_prime : Five_prime,
                                       check.this_is_left ? my_3_off : my_5_off);
            const JunctionGeometry::OffsetInterval &theirs = neighbour_reach_[check.partner_id];
            const JunctionGeometry::Overlap verdict =
                    check.this_is_left ? JunctionGeometry::check_overlap(mine, theirs, 0)
                                       : JunctionGeometry::check_overlap(theirs, mine, 0);
            if (verdict == JunctionGeometry::Overlap::Infeasible) {
                //No combination of the pending deletions separates them: a bad alignment.
                infeasible = true;
                break;
            }
            exploration.set_overlap_safety(check.safety_slot,
                                           verdict == JunctionGeometry::Overlap::Safe,
                                           check.safety_layer);
        }
        if (infeasible) {
            continue;
        }

        current_realizations_index_vec[0] = realization.index;
        new_index = base_index + current_realizations_index_vec[0];
        new_scenario_proba = base_scenario_proba;
        proba_contribution = 1;

        if (publishes_alignment_state_) {
            //The D arm never published these -- the "FIXME deal with state pointers for D" this
            //body replaces -- and the hypermutation error rates read them off the D event, so
            //the omission is observable. Carried verbatim: it is a defect for phase R.
            current_realization_index = &realization.index;
            alignment_offset_p = &alignment.offset;
        }

        proba_contribution = iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                            exploration.index_map, model.offset_map, model.model_parameters);

        new_scenario_proba *= proba_contribution;
        if (publishes_alignment_state_) {
            new_tmp_err_w_proba *= proba_contribution;
        }

        scenario.set_offset(static_cast<Seq_type>(this->seq_type_id), Five_prime, my_5_off,
                            memory_layer_off_fivep);
        scenario.set_offset(static_cast<Seq_type>(this->seq_type_id), Three_prime, my_3_off,
                            memory_layer_off_threep);
        scenario.set_mismatches(static_cast<Seq_type>(this->seq_type_id), &alignment.get_all_mismatches(),
                                memory_layer_mismatches);

        //Bound every junction this placement fixes. Which ones, in which slot and at which layer
        //was settled in initialize_event() (S4c), so there is no junction to identify here.
        if (not write_junction_bounds(exploration.downstream_proba_map, my_5_off, my_3_off)) {
            continue; //No scenario reaches one of those distances, would need to be changed for error models with in/dels
        }

        //The mismatches that survive whatever the pending deletions do: the ones inside the core
        //neither end can retract past. One window (section 2.8) for what was three, each written
        //out in terms of the one deletion its gene class happened to have.
        const Seq_Offset core_5 = pending_.reachable(this->seq_type_id, Five_prime, my_5_off).hi;
        const Seq_Offset core_3 = pending_.reachable(this->seq_type_id, Three_prime, my_3_off).lo;

        endogeneous_mismatches = 0;
        if (endogenous_core_ == EndogenousCore::Truncated and core_5 >= core_3) {
            //Nothing of this template is guaranteed to survive, so it constrains nothing.
            exploration.downstream_proba_map.set(this->seq_type_id, 1.0, memory_layer_proba_map_seq);
        } else {
            for (const size_t mismatch_position : alignment.mismatches) {
                const Seq_Offset position = static_cast<Seq_Offset>(mismatch_position);
                if (position >= core_5 and position <= core_3) {
                    ++endogeneous_mismatches;
                }
            }
            exploration.downstream_proba_map.set(
                    this->seq_type_id,
                    accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                                      credited_core_length(core_5, core_3)),
                    memory_layer_proba_map_seq);
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

    //A gene with a read end on one side is anchored by it; one with neighbours on both sides is
    //not, so when no alignment survives there is still a position range to scan. That is the
    //whole of the V/J-versus-D asymmetry the `case D_gene` this replaces carried (decision O6),
    //and the body below is B11b's, unchanged here.
    if (exhaustive_position_fallback_ and no_d_align) {
        //int test = 0;

        //Pass the mismatch vector pointer to the memory map once (will be updated in the next loop)
        scenario.set_mismatches(D_gene_seq, &no_d_mismatches, memory_layer_mismatches);

        //Record an overlap verdict at this event's layer, which the exhaustive path
        //below never does. The preamble already wrote one for a neighbour that has not
        //been chosen; when it has been, the alignment loop is what writes it, and this
        //path does not run that loop. A downstream Deletion reads
        //memory_layer_safety - 1, i.e. exactly this layer, and LayeredArray::get()
        //refuses a layer above the last one written -- so leaving it unwritten aborts
        //the run. Before the B8 container port the same read returned uninitialized
        //storage instead, which is why this went unnoticed.
        //
        //`false` means "not established safe", so the downstream deletion performs its
        //own check rather than skipping it: the conservative direction, and the same
        //value the alignment loop writes whenever the verdict is undetermined. B11
        //should compute the real verdict here, since d_5_off and d_full_3_offset are
        //known per position. See docs/ITERATE_GENERIC_REWRITE_PLAN.md section 7.9.
        if (v_chosen) {
            exploration.set_overlap_safety(flank_checks_[0].safety_slot, false, flank_checks_[0].safety_layer);
        }
        if (j_chosen) {
            exploration.set_overlap_safety(flank_checks_[1].safety_slot, false, flank_checks_[1].safety_layer);
        }

        if (v_chosen and j_chosen) {
            int vj_len = neighbour_offset_[J_gene_seq] - neighbour_offset_[V_gene_seq] - 1;
            if (vj_length_d_position_proba.count(vj_len) != 0) {
                const vector<tuple<string, int, int, double>> &d_positions_vector =
                        vj_length_d_position_proba.at(vj_len);
                for (vector<tuple<string, int, int, double>>::const_iterator d_position_iter =
                             d_positions_vector.begin();
                     d_position_iter != d_positions_vector.end(); ++d_position_iter) {

                    const Event_realization &d_real = this->event_realizations.at(get<0>(*d_position_iter));

                    //d_5_off is v 3' offset + vd junction length
                    d_5_off = neighbour_offset_[V_gene_seq] + get<1>(*d_position_iter);

                    if (d_5_off - d_5_min_del >= neighbour_reach_[J_gene_seq].hi) {
                        continue;
                    }

                    d_size = d_real.value_str.size();

                    d_full_3_offset = d_5_off + d_size - 1;
                    d_3_max_offset = d_full_3_offset + d_3_min_del;

                    if (d_3_max_offset <= neighbour_reach_[V_gene_seq].lo) {
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
									if(vd_length_best_proba_map.count(d_5_off - neighbour_offset_[V_gene_seq] -1)<=0 or dj_length_best_proba_map.count(neighbour_offset_[J_gene_seq] - d_full_3_offset  -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}*/
                    //The two distances come out of the retained decomposition, so they are
                    //entries of the profiles by construction and need no guard.
                    exploration.downstream_proba_map.set(junction_bound(kEnclosingJunction).proba_key(), 1.0,
                                                   junction_bound(kEnclosingJunction).memory_layer());
                    exploration.downstream_proba_map.set(
                            junction_bound(kLeftJunction).proba_key(),
                            *junction_bound(kLeftJunction).profile().best_for(get<1>(*d_position_iter)),
                            junction_bound(kLeftJunction).memory_layer());
                    exploration.downstream_proba_map.set(
                            junction_bound(kRightJunction).proba_key(),
                            *junction_bound(kRightJunction).profile().best_for(get<2>(*d_position_iter)),
                            junction_bound(kRightJunction).memory_layer());
                    exploration.downstream_proba_map.set(D_gene_seq, 1.0,
                                                   memory_layer_proba_map_seq); //Lift the penalty on D gene seq

                    /*							}
								else if(v_chosen){
									if(vd_length_best_proba_map.count(d_5_off - neighbour_offset_[V_gene_seq] -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}
									downstream_proba_map.set(VD_ins_seq , vd_length_best_proba_map.at(d_5_off - neighbour_offset_[V_gene_seq] -1) , memory_layer_proba_map_junction_d2);
								}
								else if(j_chosen){
									if(dj_length_best_proba_map.count(neighbour_offset_[J_gene_seq] - d_full_3_offset  -1)<=0){
										continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
									}
									downstream_proba_map.set(DJ_ins_seq , dj_length_best_proba_map.at(neighbour_offset_[J_gene_seq] - d_full_3_offset  -1) , memory_layer_proba_map_junction_d3);
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
                        exploration.downstream_proba_map.set(
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
                        exploration.downstream_proba_map.set(D_gene_seq, 1.0, memory_layer_proba_map_seq);
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

                if (neighbour_reach_[V_gene_seq].lo > 0) {
                    d_5_off = neighbour_reach_[V_gene_seq].lo + d_5_real_max_del + 1;
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

                while (d_3_min_offset < neighbour_reach_[J_gene_seq].lo) {
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
								cout<<neighbour_reach_[V_gene_seq].lo<<endl;
								cout<<d_5_max_del<<endl;
							}*/

                    //Assume that the whole D is in the sequence and add the D sequence to the constructed sequences
                    scenario.set_offset(D_gene_seq, Five_prime, d_5_off, memory_layer_off_fivep);
                    scenario.set_offset(D_gene_seq, Three_prime, d_full_3_offset, memory_layer_off_threep);


                    //Get the upper bound probas for the junctions this D placement creates
                    if (not write_junction_bounds(exploration.downstream_proba_map, d_5_off,
                                                  d_full_3_offset)) {
                        continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
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
                        exploration.downstream_proba_map.set(
                                D_gene_seq,
                                accumulation.error_rate->get_err_rate_upper_bound(endogeneous_mismatches,
                                                                       (d_full_3_offset + d_3_max_del)
                                                                               - (d_5_off - d_5_max_del)
                                                                               - endogeneous_mismatches),
                                memory_layer_proba_map_seq);
                    } else {
                        exploration.downstream_proba_map.set(D_gene_seq, 1.0, memory_layer_proba_map_seq);
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
        }    }
}

/*
 * The bound each junction this placement fixes contributes. A gene at an end of the ordering
 * fixes one; an internal one fixes both its flanks and replaces the estimate on the junction it
 * splits with the neutral element, the refinement living in vj_length_d_position_proba instead.
 * Which is which was resolved in initialize_event(), so this asks no span questions.
 *
 * Both distances are computed before either is written, so a placement that fails on the second
 * leaves nothing behind -- which is what the D arm did and the V and J arms got for free by
 * having only one.
 */
bool Gene_choice::write_junction_bounds(Downstream_scenario_proba_bound_map &proba_map,
                                       Seq_Offset five_off, Seq_Offset three_off) const
{
    const JunctionBound &left = junction_bound(kLeftJunction);
    const JunctionBound &right = junction_bound(kRightJunction);

    std::optional<double> left_proba;
    if (left.resolved()) {
        left_proba = left.profile().best_for(five_off - neighbour_offset_[left.span().left.id] - 1);
        if (not left_proba) {
            return false;
        }
    }
    std::optional<double> right_proba;
    if (right.resolved()) {
        right_proba = right.profile().best_for(neighbour_offset_[right.span().right.id] - three_off - 1);
        if (not right_proba) {
            return false;
        }
    }

    const JunctionBound &enclosing = junction_bound(kEnclosingJunction);
    if (enclosing.resolved()) {
        proba_map.set(enclosing.proba_key(), 1.0, enclosing.memory_layer());
    }
    if (left.resolved()) {
        proba_map.set(left.proba_key(), *left_proba, left.memory_layer());
    }
    if (right.resolved()) {
        proba_map.set(right.proba_key(), *right_proba, right.memory_layer());
    }
    return true;
}

Event_safety Gene_choice::legacy_safety_slot(std::size_t left_position, std::size_t right_position)
{
    //Positions in the legacy 5'->3' gene order: 0 = V, 1 = D, 2 = J.
    if (left_position == 0 and right_position == 1) {
        return VD_safe;
    }
    if (left_position == 1 and right_position == 2) {
        return DJ_safe;
    }
    if (left_position == 0 and right_position == 2) {
        return VJ_safe;
    }
    throw std::logic_error("Gene_choice: the Event_safety enum cannot name the gene pair at "
                           "ordering positions " + std::to_string(left_position) + " and "
                           + std::to_string(right_position) + "; that is what S5's safety row "
                           "bitmask replaces it with");
}

int Gene_choice::credited_core_length(Seq_Offset core_5, Seq_Offset core_3) const
{
    if (endogenous_core_ == EndogenousCore::Truncated) {
        //The core between the two movable ends. Right in shape, one short of the inclusive
        //count -- section 7.1's off-by-one, which R5 fixes together with the other arm.
        return static_cast<int>(core_3 - core_5) - static_cast<int>(endogeneous_mismatches);
    }
    //Section 7.1's sign inversion, reproduced rather than derived (decision O4). `travel` is how
    //far the one movable end can retract; the surviving core is `size - travel`, and this credits
    //`size + travel`, so the bound comes out too small and prunes harder than the model justifies.
    const int travel = pending_.offset_delta(this->seq_type_id, Five_prime).max
                     - pending_.offset_delta(this->seq_type_id, Three_prime).min;
    return static_cast<int>(gene_seq.size()) + travel - static_cast<int>(endogeneous_mismatches);
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
    write2txt_legacy(outfile);
}

void Gene_choice::write2txt_legacy(ofstream &outfile)
{
    outfile << "#GeneChoice;" << event_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (const auto &iter : event_realizations) {
        outfile << "%" << iter.second.name << ";" << iter.second.value_str << ";" << iter.second.index << endl;
    }
}

void Gene_choice::write2txt_v2(ofstream &outfile)
{
    outfile << "#GeneChoice;" << event_class << ";" << seq_type << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (const auto &iter : event_realizations) {
        outfile << "%" << iter.second.name << ";" << iter.second.value_str << ";" << iter.second.index << endl;
    }
}

void Gene_choice::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const Events_map &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    //Check V choice
    auto v_status = EventUtils::check_gene_choice("V_gene_seq", events_map, processed_events);
    v_choice_exist = v_status.exists;
    v_chosen = v_status.chosen;

    //Check D choice
    auto d_status = EventUtils::check_gene_choice("D_gene_seq", events_map, processed_events);
    d_choice_exist = d_status.exists;
    d_chosen = d_status.chosen;

    //Check J choice
    auto j_status = EventUtils::check_gene_choice("J_gene_seq", events_map, processed_events);
    j_choice_exist = j_status.exists;
    j_chosen = j_status.chosen;

    //Everything topological is settled here, so iterate() reads answers instead of asking
    //event_class again. The three switch arms of layer requests collapse because they only ever
    //differed in which seq_type they named.
    const SeqTypeRegistry &registry = constructed_sequences.registry();
    pending_.rebuild(registry, events_map, processed_events);
    neighbour_reach_.assign(registry.total_count(), JunctionGeometry::OffsetInterval{});
    neighbour_offset_.assign(registry.total_count(), 0);

    //A gene at an end of the constructed sequence is anchored by the read rather than by a
    //neighbour, and that single fact is the whole V/J-versus-D asymmetry: its template can
    //overhang the read on that side, it has no position range to scan when the aligner finds
    //nothing (decision O6), and it is the arm that credits its error-free core the wrong way
    //round (section 7.1). Read off the ordering, so a tandem D1/D2 pair gets the internal
    //behaviour without either of them being named anywhere.
    const bool at_left_end = registry.left_neighbor(this->seq_type_id) == kNoSeqType;
    const bool at_right_end = registry.right_neighbor(this->seq_type_id) == kNoSeqType;
    clip_template_before_read_ = at_left_end;
    clip_template_after_read_ = at_right_end;
    exhaustive_position_fallback_ = not at_left_end and not at_right_end;
    publishes_alignment_state_ = at_left_end or at_right_end;
    endogenous_core_ = (at_left_end or at_right_end) ? EndogenousCore::Inflated
                                                     : EndogenousCore::Truncated;

    //The neighbours this placement is checked against: the other gene segments, 5' to 3', which
    //is the order the two claimed safety layers were in. The candidate list is still the legacy
    //three by name -- parent plan B9 step 3, which puts the whole topology in the registry, is
    //what replaces this literal with registry.ordering(); until then a VJ model must still claim
    //the VD safety slot it never uses, exactly as the preambles this replaces did.
    static const std::array<Seq_type_String, 3> kGeneSegments = {"V_gene_seq", "D_gene_seq",
                                                                 "J_gene_seq"};
    const std::size_t my_position =
            std::find(kGeneSegments.begin(), kGeneSegments.end(), registry.name(this->seq_type_id))
            - kGeneSegments.begin();
    flank_checks_.clear();
    for (std::size_t position = 0; position != kGeneSegments.size(); ++position) {
        if (position == my_position or not registry.contains(kGeneSegments[position])) {
            continue;
        }
        const Seq_type_String &partner_name = kGeneSegments[position];
        const igor::migration::GeneChoiceStatus status =
                EventUtils::check_gene_choice(partner_name, events_map, processed_events);
        FlankCheck check;
        check.partner_id = registry.id(partner_name);
        check.this_is_left = position > my_position;
        check.partner_side = check.this_is_left ? Five_prime : Three_prime;
        check.safety_slot = legacy_safety_slot(std::min(position, my_position),
                                               std::max(position, my_position));
        check.partner_exists = status.exists;
        check.partner_chosen = status.chosen;
        flank_checks_.push_back(check);
    }

    seq_offsets.request_layer(this->seq_type_id, Three_prime);
    this->memory_layer_off_threep = seq_offsets.claimed_layer(this->seq_type_id, Three_prime);
    seq_offsets.request_layer(this->seq_type_id, Five_prime);
    this->memory_layer_off_fivep = seq_offsets.claimed_layer(this->seq_type_id, Five_prime);
    mismatches_list.request_layer(this->seq_type_id);
    this->memory_layer_mismatches = mismatches_list.claimed_layer(this->seq_type_id);
    constructed_sequences.request_layer(this->seq_type_id);
    this->memory_layer_cs = constructed_sequences.claimed_layer(this->seq_type_id);

    //Claimed whether or not the neighbour is there, which is what the commented-out `if`s in the
    //three arms this replaces were about: a downstream Deletion reads this layer unconditionally.
    for (FlankCheck &check : flank_checks_) {
        safety_set.request_layer(check.safety_slot);
        check.safety_layer = safety_set.claimed_layer(check.safety_slot);
    }

    downstream_proba_map.request_layer(this->seq_type_id);
    memory_layer_proba_map_seq = downstream_proba_map.claimed_layer(this->seq_type_id);

    //The junctions this event's placement fixes: one against the nearest placed neighbour on
    //each side, plus -- when there is one on both -- the junction it sits inside and splits.
    //Three slots is enough in every topology (see JunctionBound), which is what removed the enum
    //ceiling in S4c; here is where they are filled.
    const FlankCheck *left_partner = nullptr;
    const FlankCheck *right_partner = nullptr;
    for (const FlankCheck &check : flank_checks_) {
        if (not check.partner_chosen) {
            continue;
        }
        if (check.this_is_left) {
            if (right_partner == nullptr) {
                right_partner = &check; //first placed neighbour to the 3' side
            }
        } else {
            left_partner = &check;      //last placed neighbour to the 5' side
        }
    }

    auto resolve_junction = [&](JunctionSlot slot, SeqTypeId left_id, SeqTypeId right_id,
                                JunctionBound::Fold fold) {
        const SegmentSpan span = SegmentSpan::gap(left_id, right_id);
        const SeqTypeId junction_key = legacy_junction_of(span);
        downstream_proba_map.request_layer(junction_key);
        junction_bound(slot).resolve(span, junction_key, downstream_proba_map.claimed_layer(junction_key),
                                     fold);
    };

    if (left_partner != nullptr) {
        resolve_junction(kLeftJunction, left_partner->partner_id, this->seq_type_id,
                         JunctionBound::Fold::Yes);
    }
    if (right_partner != nullptr) {
        resolve_junction(kRightJunction, this->seq_type_id, right_partner->partner_id,
                         JunctionBound::Fold::Yes);
    }
    if (left_partner != nullptr and right_partner != nullptr) {
        //Placing this segment refines the enclosing junction into the two halves above rather
        //than measuring it, so it carries no folded profile of its own: iterate() writes the
        //neutral 1.0 there and the decomposition lives in vj_length_d_position_proba.
        resolve_junction(kEnclosingJunction, left_partner->partner_id, right_partner->partner_id,
                         JunctionBound::Fold::No);
    }

    for (FlankCheck &check : flank_checks_) {
        if (check.partner_chosen) {
            check.partner_offset_layer = seq_offsets.claimed_layer(check.partner_id, check.partner_side);
        }
    }

    //downstream_proba_map.get_all_current_memory_layer(current_downstream_proba_memory_layers);


    //The two D deletion ranges survive only because the no_d_align path below still reads the
    //scalars; `pending_` answers the same question for everyone else, and these go with that
    //path in 5b. The V 3' and J 5' lookups they used to sit between are already gone.
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

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
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

OffsetDelta Gene_choice::get_offset_delta_bounds(SeqTypeId, Seq_side) const
{
    //A gene choice *sets* both offsets from the alignment; it never shifts an existing one,
    //so it contributes no travel to any end. What its templates make possible is a length,
    //reported by get_length_contribution().
    return {};
}

LengthContribution Gene_choice::get_length_contribution(SeqTypeId type_id) const
{
    if (type_id != this->seq_type_id || this->event_realizations.empty()) {
        return {};
    }
    int shortest = std::numeric_limits<int>::max();
    int longest = std::numeric_limits<int>::min();
    for (const auto &[name, realization] : this->event_realizations) {
        (void)name;
        const int length = static_cast<int>(realization.value_str.size());
        shortest = std::min(shortest, length);
        longest = std::max(longest, length);
    }
    //The whole template. The portion clipped for an overhanging alignment is not subtracted
    //here: it depends on the query, not on the model, and task B3 turns it into a flank
    //segment of its own rather than a shortening of this one.
    return {shortest, longest};
}

SeqConstructionRole Gene_choice::get_seq_construction_role(SeqTypeId type_id) const
{
    return type_id == this->seq_type_id ? SeqConstructionRole::Creates : SeqConstructionRole::None;
}

OffsetRole Gene_choice::get_offset_role(SeqTypeId type_id, Seq_side) const
{
    //Both ends, since the alignment fixes the segment's position outright.
    return type_id == this->seq_type_id ? OffsetRole::Creates : OffsetRole::None;
}


bool Gene_choice::affects_length_of(SegmentSpan span) const
{
    //A gene's template contributes length to a span only when the gene sits strictly *inside*
    //it: at either end it is the anchor the span is measured from, so its own length is outside
    //the frame. Over the legacy junctions that leaves exactly one case, D within V->J.
    //S4b generalises this off the enum, once the traversal carries the registry ordering.
    switch (this->event_class) {
    case V_gene:
        return false;

    case D_gene:
        return legacy_junction_of(span) == VJ_ins_seq;

    case J_gene:
        return false;

    default:
        return false;
    }
}

int Gene_choice::length_delta(const Event_realization &realization) const
{
    //The chosen template's length, before any deletion trims it: the gene creates the segment.
    return static_cast<int>(realization.value_str.length());
}

/*
 * The retained decomposition of the junction this D splits: for each achievable V->J distance,
 * every (D gene, VD distance, DJ distance) that reaches it, best first. Section 2.5's composition
 * operator in the variant that keeps its arguments rather than max-folding them -- which is why it
 * cannot be rebuilt from the two profiles after the fact, and why it is built here rather than
 * derived on demand.
 *
 * Called from Rec_Event::initialize_Len_proba_bound() once every folded junction of this event
 * exists, so the profiles it reads are the current EM iteration's. The three-way enum switch it
 * replaces is gone entirely: a V or a J gene choice has nothing to add beyond its folded profile,
 * which the base class already built.
 */
void Gene_choice::finalize_Len_proba_bound(const Marginal_array_p &model_parameters_point, Index_map &base_index_map)
{
    if (this->event_class != D_gene) {
        return;
    }

    vj_length_d_position_proba.clear();

    const JunctionBound &left = junction_bound(kLeftJunction);
    const JunctionBound &right = junction_bound(kRightJunction);
    if (not left.resolved() or not right.resolved()) {
        return;
    }

    //Loop over D gene choices
    for (unordered_map<string, Event_realization>::const_iterator d_gene_iter = this->event_realizations.begin();
         d_gene_iter != this->event_realizations.end(); ++d_gene_iter) {
        //Get considered D gene best proba
        double d_gene_max_proba = 0;
        base_index_map.set_current_layer(this->event_index, 0);
        base_index = base_index_map.get(this->event_index);
        for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
            if (model_parameters_point[base_index + d_gene_iter->second.index + i * this->size()] > d_gene_max_proba) {
                d_gene_max_proba = model_parameters_point[base_index + d_gene_iter->second.index + i * this->size()];
            }
        }
        //Loop over possible VD junction lengths
        for (const SpanProfile::Entry vd : left.profile()) {
            //Loop over possible DJ junction lengths
            for (const SpanProfile::Entry dj : right.profile()) {
                const int junction_len = d_gene_iter->second.value_str.size() + vd.distance + dj.distance;
                const double position_proba = d_gene_max_proba * vd.proba * dj.proba;

                if (vj_length_d_position_proba.count(junction_len) != 0) {
                    vj_length_d_position_proba.at(junction_len)
                            .emplace_back(d_gene_iter->first, vd.distance, dj.distance, position_proba);
                } else {
                    vj_length_d_position_proba.emplace(
                            piecewise_construct, make_tuple(junction_len),
                            make_tuple(1, make_tuple(d_gene_iter->first, vd.distance, dj.distance, position_proba)));
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

/*
 * The adopting half of the above. vj_length_d_position_proba is not in any JunctionBound, so the
 * base class cannot move it; both this and finalize_Len_proba_bound() disappear when 5b turns the
 * retained decomposition into JunctionBound::Fold::Retain.
 */
void Gene_choice::adopt_finalized_Len_proba_bound(const Rec_Event &source)
{
    if (this->event_class != D_gene) {
        return;
    }

    //The caller passes the same event of another thread's model copy, so this is a Gene_choice(D).
    const Gene_choice &gene_source = dynamic_cast<const Gene_choice &>(source);
    vj_length_d_position_proba = gene_source.vj_length_d_position_proba;
}
