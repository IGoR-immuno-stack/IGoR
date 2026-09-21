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
      placement_5_off(INT16_MAX),
      no_alignment_survived(true),
      template_size(INT16_MAX),
      placement_3_off(INT16_MAX),
      base_index(-1),
      new_scenario_proba(-1),
      proba_contribution(-1),
      new_index(-1),
      alignment_offset_p(NULL),
      memory_layer_cs(-1),
      memory_layer_mismatches(-1),
      memory_layer_off_threep(-1),
      memory_layer_off_fivep(-1)
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
            exploration.set_overlap_safety(check.safety_cell, not check.partner_exists, check.safety_layer);
            if (exhaustive_position_fallback_ and check.this_is_left) {
                //Nothing placed to the right, so the exhaustive scan's right wall is the end of
                //the read. Read only by the no_alignment_survived path below, which is 5b's.
                const Seq_Offset read_end = static_cast<Seq_Offset>(query.sequence.size()) - 1;
                neighbour_reach_[check.partner_id] = {read_end, read_end};
            }
        }
    }

    no_alignment_survived = true;

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
            exploration.set_overlap_safety(check.safety_cell,
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
        no_alignment_survived = false;

        // Update context with new probability before proceeding
        scenario.scenario_proba = new_scenario_proba;

        Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
    }

    //A gene with a read end on one side is anchored by it; one with neighbours on both sides is
    //not, so when no alignment survives there is still a position range to scan. That is the
    //whole of the V/J-versus-D asymmetry the `case D_gene` this replaces carried (decision O6).
    //
    //Two sub-branches, and what separates them is whether the segment has a placed neighbour on
    //*both* sides. With both, the span between them is known and the placements come out of its
    //retained decomposition, ordered by decreasing bound. With fewer, there is nothing to
    //decompose and the window slides against the read.
    if (exhaustive_position_fallback_ and no_alignment_survived) {
        const Seq_type my_seq_type = static_cast<Seq_type>(this->seq_type_id);

        //Pass the mismatch vector pointer to the memory map once (will be updated in the next loop)
        scenario.set_mismatches(my_seq_type, &placement_mismatches, memory_layer_mismatches);

        //Record an overlap verdict at this event's layer, which the scans below never do. The
        //preamble already wrote one for a neighbour that has not been placed; when it has, the
        //alignment loop is what writes it, and this path does not run that loop. A downstream
        //Deletion reads `memory_layer_safety - 1`, i.e. exactly this layer, and a claimed but
        //unwritten layer is unreadable -- so leaving it unwritten aborts the run (§7.9).
        //
        //`false` means "not established safe", so the downstream deletion performs its own check
        //rather than skipping it: the conservative direction, and the same value the alignment
        //loop writes whenever the verdict is undetermined. The real verdict is computable here,
        //since both ends of the placement are known -- a tightening, and therefore not this
        //step's (§1).
        for (const FlankCheck &check : flank_checks_) {
            if (check.partner_chosen) {
                exploration.set_overlap_safety(check.safety_cell, false, check.safety_layer);
            }
        }

        const JunctionBound &enclosing = junction_bound(kEnclosingJunction);
        const JunctionBound &left = junction_bound(kLeftJunction);
        const JunctionBound &right = junction_bound(kRightJunction);

        if (enclosing.resolved()) {
            //Both flanks placed. Every way this segment can sit between them, at the distance
            //they actually leave, best first -- §2.5's `⊗ᵉⁿᵘᵐ`, built at initialization by the
            //`Fold::Retain` mode and read here by total distance.
            const SeqTypeId left_id = enclosing.span().left.id;
            const SeqTypeId right_id = enclosing.span().right.id;
            const int span_len = neighbour_offset_[right_id] - neighbour_offset_[left_id] - 1;

            for (const SpanDecomposition::Placement &placement : enclosing.decomposition().at(span_len)) {
                const Event_realization &realization = *realizations_by_index_[placement.realization_index];

                //The 5' end sits one junction-length past the left neighbour's 3' end.
                placement_5_off = neighbour_offset_[left_id] + placement.left_distance;
                template_size = realization.value_str.size();
                placement_3_off = placement_5_off + template_size - 1;

                //Can this placement still clear each neighbour once its own deletions have had
                //their say? Same predicate as everywhere else (§2.2), with this segment's end
                //carrying the travel rather than the neighbour's.
                const JunctionGeometry::OffsetInterval my_five_reach = own_five_prime_reach(placement_5_off);
                const JunctionGeometry::OffsetInterval my_three_reach = own_three_prime_reach(placement_3_off);
                if (JunctionGeometry::check_overlap(my_five_reach, neighbour_reach_[right_id], 0)
                    == JunctionGeometry::Overlap::Infeasible) {
                    continue;
                }
                if (JunctionGeometry::check_overlap(neighbour_reach_[left_id], my_three_reach, 0)
                    == JunctionGeometry::Overlap::Infeasible) {
                    continue;
                }

                gene_seq = realization.value_str_int;
                scenario.set_sequence_segment(my_seq_type, &gene_seq, memory_layer_cs);

                current_realizations_index_vec[0] = realization.index;
                new_index = base_index + current_realizations_index_vec[0];

                //Proba contribution is the same wherever the segment sits
                proba_contribution = 1;
                proba_contribution =
                        iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                       exploration.index_map, model.offset_map, model.model_parameters);

                //§7.16, carried: this reads the *live* scenario probability, which the hand-off
                //below has already multiplied for the previous placement, so the placements
                //compound. Every one of them is the same realization of the same event and they
                //should all carry the value this event inherited. R7 repairs it, and moves the
                //no_alignment_survived golden data with it.
                new_scenario_proba = scenario.scenario_proba * proba_contribution;

                //The two distances come out of the retained decomposition, so they are entries
                //of the profiles by construction and need no guard.
                exploration.downstream_proba_map.set(enclosing.proba_key(), 1.0, enclosing.memory_layer());
                exploration.downstream_proba_map.set(left.proba_key(),
                                                     *left.profile().best_for(placement.left_distance),
                                                     left.memory_layer());
                exploration.downstream_proba_map.set(right.proba_key(),
                                                     *right.profile().best_for(placement.right_distance),
                                                     right.memory_layer());
                //Lift the penalty on this segment
                exploration.downstream_proba_map.set(my_seq_type, 1.0, memory_layer_proba_map_seq);

                //Multiply all downstream probas
                scenario_upper_bound_proba = exploration.compute_upper_bound(
                    new_scenario_proba,
                    current_downstream_proba_memory_layers
                );

                //If even without taking the weight of errors into account not good, then any
                //lower one not good -- exact only because the decomposition is sorted by
                //decreasing bound, which is the property 5a pinned (§6.15).
                if (exploration.should_prune(scenario_upper_bound_proba)) {
                    break;
                }

                score_placement_against_read(query, accumulation, exploration, my_five_reach,
                                             my_three_reach);

                if (my_five_reach.hi < my_three_reach.lo) {
                    //Something of this template is guaranteed to survive, so the error bound it
                    //carries is worth re-testing against.
                    scenario_upper_bound_proba = exploration.compute_upper_bound(
                        new_scenario_proba,
                        current_downstream_proba_memory_layers
                    );

                    if (exploration.should_prune(scenario_upper_bound_proba)) {
                        continue;
                    }
                }

                //Assume that the whole template is in the sequence and record where it sits
                scenario.set_offset(my_seq_type, Five_prime, placement_5_off, memory_layer_off_fivep);
                scenario.set_offset(my_seq_type, Three_prime, placement_3_off, memory_layer_off_threep);

                // Update context with new probability before proceeding
                scenario.scenario_proba = new_scenario_proba;

                Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
            }
        } else if (nearest_left_check_ >= 0 and nearest_right_check_ >= 0) {
            //Fewer than two placed neighbours, so there is no span to decompose: slide the
            //template through the window the read and whichever neighbour *is* placed leave.
            //It needs a candidate on each side to bound the window, placed or not -- an
            //unplaced one contributes the read's end, which the preamble already recorded.
            const SeqTypeId left_id = flank_checks_[nearest_left_check_].partner_id;
            const SeqTypeId right_id = flank_checks_[nearest_right_check_].partner_id;

            for (unordered_map<string, Event_realization>::const_iterator gene_iter =
                         this->event_realizations.begin();
                 gene_iter != this->event_realizations.end(); ++gene_iter) {

                //Take care of the fact that not all templates have the same length, and that
                //the maximum number of deletions might be greater than the template itself
                template_size = gene_iter->second.value_str.size();
                //The unsigned comparison is the one the scalar pair this replaces made, and
                //it is kept: a template shorter than the deletion range cannot be trimmed past
                //its own start.
                const int five_prime_travel =
                        static_cast<std::size_t>(own_five_prime_travel_.max) > template_size
                                ? -static_cast<int>(template_size)
                                : -own_five_prime_travel_.max;

                //Start one nucleotide after the left neighbour's furthest reach, given maximum
                //deletions on this segment's 5' end.
                //FIXME left reach set to -1 when the neighbour is not placed
                if (neighbour_reach_[left_id].lo > 0) {
                    placement_5_off = neighbour_reach_[left_id].lo + five_prime_travel + 1;
                } else {
                    //Consider that the left neighbour cannot be absent from the read: at least
                    //one nucleotide of it is present
                    placement_5_off = 1 + five_prime_travel + 1;
                }

                placement_3_off = placement_5_off + template_size - 1;
                //Likewise on the 3' side: where the window's trailing edge starts from.
                Seq_Offset three_prime_floor =
                        std::abs(own_three_prime_travel_.min) < static_cast<int>(template_size)
                                ? placement_3_off + own_three_prime_travel_.min
                                : placement_5_off;

                //Always the same sequence for the given template
                gene_seq = gene_iter->second.value_str_int;
                scenario.set_sequence_segment(my_seq_type, &gene_seq, memory_layer_cs);

                current_realizations_index_vec[0] = gene_iter->second.index;
                new_index = base_index + current_realizations_index_vec[0];

                //Proba contribution is the same wherever the segment sits
                proba_contribution = 1;
                proba_contribution =
                        iterate_common(proba_contribution, current_realizations_index_vec[0], base_index,
                                       exploration.index_map, model.offset_map, model.model_parameters);

                //Slide one nucleotide at a time towards 3'. The advance is in the header, which
                //is what makes §7.17 -- a `continue` that skipped the increments and hung --
                //unexpressible rather than fixed.
                for (; three_prime_floor < neighbour_reach_[right_id].lo;
                     ++placement_5_off, ++placement_3_off, ++three_prime_floor) {

                    const JunctionGeometry::OffsetInterval my_five_reach = own_five_prime_reach(placement_5_off);
                    const JunctionGeometry::OffsetInterval my_three_reach =
                            own_three_prime_reach(placement_3_off);

                    //§7.16 again: the live read, and the same compounding. See above.
                    new_scenario_proba = scenario.scenario_proba * proba_contribution;

                    //Assume that the whole template is in the sequence and record where it sits
                    scenario.set_offset(my_seq_type, Five_prime, placement_5_off, memory_layer_off_fivep);
                    scenario.set_offset(my_seq_type, Three_prime, placement_3_off, memory_layer_off_threep);

                    //Get the upper bound probas for the junctions this placement creates
                    if (not write_junction_bounds(exploration.downstream_proba_map, placement_5_off,
                                                  placement_3_off)) {
                        continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
                    }

                    score_placement_against_read(query, accumulation, exploration, my_five_reach,
                                                 my_three_reach);

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
            }
        }
    }
}

/*
 * The mismatches a placement scores against the read, and the error bound they earn it.
 *
 * Shared by the two exhaustive sub-branches, which wrote it out twice identically. The window is
 * the part of the template neither end can retract past -- the same "surviving core" the
 * alignment path computes, measured here from the placement's own ends because there is no
 * alignment to read it off.
 *
 * When the two reaches cross, nothing of the template is guaranteed to survive and it constrains
 * nothing: the slot goes neutral. The alignment path spells the same case `EndogenousCore::
 * Truncated and core_5 >= core_3`.
 */
void Gene_choice::score_placement_against_read(const QuerySequenceContext &query,
                                               AccumulationContext &accumulation,
                                               ExplorationContext &exploration,
                                               JunctionGeometry::OffsetInterval five_reach,
                                               JunctionGeometry::OffsetInterval three_reach)
{
    placement_mismatches.clear();
    for (std::size_t i = 0; i != template_size; ++i) {
        const Seq_Offset position = placement_5_off + static_cast<Seq_Offset>(i);
        if (position >= 0 and static_cast<std::size_t>(position) < query.int_sequence.size()) {
            if (gene_seq[i] != query.int_sequence[position]) {
                placement_mismatches.push_back(position);
            }
        }
    }

    //Count the number of mismatches that will not go away even with maximum deletions
    endogeneous_mismatches = 0;
    if (five_reach.hi < three_reach.lo) {
        for (placement_mism_iter = placement_mismatches.begin(); placement_mism_iter != placement_mismatches.end(); ++placement_mism_iter) {
            if (static_cast<Seq_Offset>(*placement_mism_iter) >= five_reach.hi
                and static_cast<Seq_Offset>(*placement_mism_iter) <= three_reach.lo) {
                //Count one mismatch
                ++endogeneous_mismatches;
            }
        }
        exploration.downstream_proba_map.set(
                static_cast<Seq_type>(this->seq_type_id),
                accumulation.error_rate->get_err_rate_upper_bound(
                        endogeneous_mismatches,
                        three_reach.lo - five_reach.hi - endogeneous_mismatches),
                memory_layer_proba_map_seq);
    } else {
        exploration.downstream_proba_map.set(static_cast<Seq_type>(this->seq_type_id), 1.0,
                                             memory_layer_proba_map_seq);
    }
}

/*
 * The bound each junction this placement fixes contributes. A gene at an end of the ordering
 * fixes one; an internal one fixes both its flanks and replaces the estimate on the junction it
 * splits with the neutral element, the refinement living in that junction's retained
 * decomposition instead.
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
        SafetyMatrix &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
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
        check.safety_cell = safety_set.cell(this->seq_type_id, check.partner_id);
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

    //The nearest candidate on each side, placed or not: what bounds the sliding window when
    //there is no span to decompose. flank_checks_ is in 5'->3' order, so the last partner to
    //the 5' side and the first to the 3' side are the nearest ones.
    nearest_left_check_ = -1;
    nearest_right_check_ = -1;
    for (std::size_t i = 0; i != flank_checks_.size(); ++i) {
        if (flank_checks_[i].this_is_left) {
            if (nearest_right_check_ < 0) {
                nearest_right_check_ = static_cast<int>(i);
            }
        } else {
            nearest_left_check_ = static_cast<int>(i);
        }
    }

    //This event's realizations by index, so the retained decomposition can name one without a
    //string hash per candidate inside the enumeration (§2.5).
    realizations_by_index_.assign(this->size(), nullptr);
    for (const auto &[realization_name, realization] : this->event_realizations) {
        (void)realization_name;
        realizations_by_index_[realization.index] = &realization;
    }

    //Claimed whether or not the neighbour is there, which is what the commented-out `if`s in the
    //three arms this replaces were about: a downstream Deletion reads this layer unconditionally.
    //One claim per *row*, not per check: two neighbours on the same side of this gene name two
    //cells of one row, which share a word and therefore a layer. Claiming twice would leave the
    //second layer unwritten, and an unwritten layer is unreadable by design (section 7.9).
    for (std::size_t i = 0; i != flank_checks_.size(); ++i) {
        FlankCheck &check = flank_checks_[i];
        std::size_t earlier = 0;
        while (earlier != i and flank_checks_[earlier].safety_cell.row != check.safety_cell.row) {
            ++earlier;
        }
        if (earlier != i) {
            check.safety_layer = flank_checks_[earlier].safety_layer;
            continue;
        }
        safety_set.request_layer(check.safety_cell);
        check.safety_layer = safety_set.claimed_layer(check.safety_cell);
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
        //neutral 1.0 there. A gene that will scan positions when the aligner finds nothing needs
        //to know *which* composition reaches each total and not only the best one, so it asks
        //for the decomposition instead -- §2.5's `⊗ᵉⁿᵘᵐ`, and the last thing in the sweep that
        //used to be a subclass hook.
        resolve_junction(kEnclosingJunction, left_partner->partner_id, right_partner->partner_id,
                         exhaustive_position_fallback_ ? JunctionBound::Fold::Retain
                                                       : JunctionBound::Fold::No);
    }

    for (FlankCheck &check : flank_checks_) {
        if (check.partner_chosen) {
            check.partner_offset_layer = seq_offsets.claimed_layer(check.partner_id, check.partner_side);
        }
    }

    //downstream_proba_map.get_all_current_memory_layer(current_downstream_proba_memory_layers);


    //How far this segment's own two ends can still be trimmed, for the exhaustive scan below --
    //the only consumer left of a travel range on *this* event's segment; the alignment path
    //asks `pending_`. The two hand-written lookups of "the D deletion on side X" are gone, so a
    //tandem D1/D2 pair needs no new code here either.
    //
    //It is still `legacy_offset_delta` and not `pending_`: §7.4, the same short bound 4b had to
    //reproduce, and the same one-line swap in R10.
    own_five_prime_travel_ = JunctionGeometry::legacy_offset_delta(this->seq_type_id, Five_prime,
                                                                   events_map, processed_events);
    own_three_prime_travel_ = JunctionGeometry::legacy_offset_delta(this->seq_type_id, Three_prime,
                                                                    events_map, processed_events);

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

