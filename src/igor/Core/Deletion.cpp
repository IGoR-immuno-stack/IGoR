/*
 * Deletion.cpp
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

#include <igor/Core/Deletion.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/gene_to_seqtype_migr.h>

#include <algorithm>
#include <limits>
#include <utility>

using namespace std;

namespace {

Seq_type get_deletion_target_seq_type(Gene_class gene_class)
{
    switch (gene_class) {
    case V_gene:
        return V_gene_seq;
    case D_gene:
        return D_gene_seq;
    case J_gene:
        return J_gene_seq;
    default:
        throw invalid_argument(std::string("Unknown gene for deletions : ") + to_string(gene_class));
    }
}

Gene_class get_deletion_gene_class(Seq_type target_seq)
{
    switch (target_seq) {
    case V_gene_seq:
        return V_gene;
    case D_gene_seq:
        return D_gene;
    case J_gene_seq:
        return J_gene;
    default:
        throw invalid_argument(std::string("Unknown seq_type for deletions"));
    }
}

vector<Seq_type> get_deletion_effective_junctions(Seq_type target_seq_type, Seq_side event_side)
{
    switch (target_seq_type) {
    case V_gene_seq:
        return { VD_ins_seq, VJ_ins_seq };
    case D_gene_seq:
        if (event_side == Five_prime) {
            return { VD_ins_seq };
        }
        if (event_side == Three_prime) {
            return { DJ_ins_seq };
        }
        return {};
    case J_gene_seq:
        return { VJ_ins_seq, DJ_ins_seq };
    default:
        return {};
    }
}

} // namespace

Deletion::Deletion() : Deletion(V_gene_seq, Undefined_side)
{
    this->event_class = Undefined_gene;
    this->type = Event_type::Deletion_t;
    this->update_event_name();
}

Deletion::Deletion(Seq_type target_seq, Seq_side del_side, pair<int, int> del_range) : Deletion(target_seq, del_side)
{
    //TODO prevent undefined_side entry(throw exception)

    int min_del = std::min(del_range.first, del_range.second);
    int max_del = std::max(del_range.first, del_range.second);

    this->type = Event_type::Deletion_t;
    this->len_max = -min_del;
    this->len_min = -max_del;

    for (int i = min_del; i != max_del + 1; ++i) {
        this->add_realization(i);
    }
    this->update_event_name();
}

Deletion::Deletion(Seq_type target_seq, Seq_side side)
    : Rec_Event(get_deletion_gene_class(target_seq), side),
      target_seq_type(target_seq),
      base_index(INT16_MAX),
      new_scenario_proba(-1),
      proba_contribution(-1),
      new_index(-1),
      end_reached(false),
      deletion_value(INT16_MAX),
      memory_layer_cs(-1),
      memory_layer_mismatches(-1),
      memory_layer_offset_del(-1),
      memory_layer_proba_map_seq(-1)
{
    this->type = Event_type::Deletion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > (-this->len_min)) {
            this->len_min = -(*iter).second.value_int;
        } else if ((*iter).second.value_int < (-this->len_max)) {
            this->len_max = -(*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Deletion::Deletion(Seq_type target_seq, Seq_side side, unordered_map<string, Event_realization> &realizations)
    : Deletion(target_seq, side)
{
    this->event_realizations = realizations;

    this->type = Event_type::Deletion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
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
            new Deletion(this->target_seq_type, this->event_side,
                         this->event_realizations)); //FIXME remove this new for all events and for the error rates
    new_deletion_p->priority = this->priority;
    new_deletion_p->nickname = this->nickname;
    new_deletion_p->fixed = this->is_fixed();
    new_deletion_p->update_event_name();
    new_deletion_p->set_event_identifier(this->event_index);
    new_deletion_p->set_seq_type(this->get_seq_type());
    new_deletion_p->set_seq_type_id(this->get_seq_type_id());
    return new_deletion_p;
}

void Deletion::add_realization(int del_number)
{

    this->Rec_Event::add_realization(Event_realization(to_string(del_number), del_number, "", Int_Str(), this->size()));

    if (del_number > (-this->len_min)) {
        this->len_min = (-del_number);
    } else if (del_number < (-this->len_max)) {
        this->len_max = (-del_number);
    }
    this->update_event_name();
}

/**
 * @brief Context-based iterate() implementation
 *
 * General: Loop over all possible number of deletions for a given gene on a given sequence side
 *
 * Specific:
 * -First check whether any of these number of deletions is possible given the current position and number of deletions on other genes
 * -Loop over # of deletions in decreasing order
 */
/*
 * One body for V 3', D 5', D 3' and J 5'. Once initialize_event() has settled the topology, the
 * four arms this replaces differed in exactly three things:
 *
 *   - **which side of its segment the deletion moves** -- and therefore which neighbours it is
 *     compared against, which end of the template it trims, which junction it widens, which way
 *     the mismatch list is cut and whether a palindrome's new positions need re-sorting. All of
 *     that is `event_side`, and it is read here rather than switched on;
 *   - **whether the segment is anchored on an end of the read** -- which decides whether the
 *     template may be deleted away entirely (§2.7), whether the enumeration carries the
 *     dominated early prune stage (§6.14), and how the palindrome is bounded against the read
 *     (§7.15). Read off the ordering, exactly as `Gene_choice` does since B11a;
 *   - **whether anything still moves the opposite end** -- which decides whether this segment's
 *     error bound can be written at all, and which used to be spelled `d_del_opposite_side_processed`.
 *
 * None of those is a gene class, which is why a tandem D needs no fifth case.
 */
void Deletion::iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    base_index = exploration.index_map.get(this->event_index);
    const double base_scenario_proba = scenario.scenario_proba;
    const Seq_type my_seq_type = static_cast<Seq_type>(this->seq_type_id);

    my_offset = scenario.get_offset(my_seq_type, this->event_side, memory_layer_offset_del - 1);

    //Where each checked neighbour sits, read once per scenario, and whether it still has to be
    //compared against at all. A pair the enclosing depth already established as separated needs
    //no comparison here; re-asserting the verdict at this event's layer is what carries it
    //down to the next reader, and is what the four `else { check = false; set safe }` arms did.
    for (FlankCheck &check : flank_checks_) {
        check.active = false;
        if (not check.partner_chosen) {
            continue;
        }
        //Read whether or not the comparison runs. The V arm this replaces read it only when it
        //did, and then used it to measure the junction -- so a V deletion in a model with no D
        //measured the V->J junction against a stale, and on the first scenario uninitialized,
        //offset. See §7.21.
        check.offset = scenario.get_offset(static_cast<Seq_type>(check.partner_id),
                                           check.partner_side, check.partner_offset_layer);
        if (exploration.is_overlap_safe(check.safety_cell, check.safety_layer - 1)) {
            exploration.set_overlap_safety(check.safety_cell, true, check.safety_layer);
        } else {
            check.reach = {static_cast<Seq_Offset>(check.offset + check.partner_delta.min),
                           static_cast<Seq_Offset>(check.offset + check.partner_delta.max)};
            check.active = true;
        }
    }

    Int_Str &previous_str = *scenario.get_sequence_segment(my_seq_type, memory_layer_cs - 1);
    const vector<size_t> &previous_mismatches =
            *scenario.get_mismatches(my_seq_type, memory_layer_mismatches - 1);

    //The junction on the side this deletion trims. Which span that is, where its value goes and
    //at which layer were settled in initialize_event() (S4c); all iterate() does is measure the
    //distance.
    const JunctionBound &junction = junction_bound(trims_three_prime_ ? kRightJunction : kLeftJunction);

    for (forward_list<Event_realization>::const_iterator iter = int_value_and_index.begin();
         iter != int_value_and_index.end(); ++iter) {
        const int deletions = iter->value_int;

        //A segment anchored on the read keeps at least one nucleotide of its template; an
        //internal one may delete itself away entirely, leaving a written-but-empty segment at
        //the degenerate offsets `three_prime == five_prime - 1`. §2.7 calls this a modelling
        //decision rather than an accident, and 4a pins it in all four arms.
        if (keep_one_nucleotide_ ? static_cast<int>(previous_str.size()) <= deletions
                                 : static_cast<int>(previous_str.size()) < deletions) {
            continue;
        }

        my_new_offset = trims_three_prime_ ? my_offset - deletions : my_offset + deletions;

        //§7.3's surviving `//FIXME`, carried verbatim *including* its unsigned comparison: a
        //negative offset converts to a huge size_t and is rejected here too, which is the
        //second job 4a found it doing (§6.14). Only an internal segment's 5' deletion has it.
        if (discard_offset_outside_read_
            and static_cast<std::size_t>(my_new_offset) >= query.int_sequence.size()) {
            continue;
        }

        //Can this deletion still avoid colliding with each checked neighbour? One predicate
        //(§2.2) for what was eight hand-written comparisons, differing only in which side of the
        //pair this segment sits on. The moving end is a point here -- the realization pins it --
        //which is the difference §2.3's `reachable()` absorbs.
        bool infeasible = false;
        for (const FlankCheck &check : flank_checks_) {
            if (not check.active) {
                continue;
            }
            const JunctionGeometry::OffsetInterval mine{my_new_offset, my_new_offset};
            const JunctionGeometry::Overlap verdict =
                    trims_three_prime_ ? JunctionGeometry::check_overlap(mine, check.reach, 0)
                                       : JunctionGeometry::check_overlap(check.reach, mine, 0);
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

        //Store the deletion value in a private variable: the error rates and the coverage
        //counter read it off this event by pointer.
        deletion_value = deletions;

        current_realizations_index_vec[0] = iter->index;
        new_index = base_index + iter->index;
        new_scenario_proba = base_scenario_proba;
        proba_contribution = 1;

        this->iterate_common(iter, exploration.index_map, model.offset_map, model.model_parameters);

        //Positive or negative deletion (palindroms) mechanism
        if (deletions >= 0) {
            if (trims_three_prime_) {
                previous_str.substr(new_str, 0, previous_str.size() - deletions);
            } else {
                previous_str.substr(new_str, deletions, string::npos);
            }
            trim_mismatches(previous_mismatches);
        } else {
            //A negative deletion is a palindromic insertion: the |k| nucleotides nearest the
            //trimmed end are reversed, complemented and put back on the far side of it.
            //
            //Only a 3' trim checks that its new end is still inside the read before scoring the
            //palindrome against it. That asymmetry is §7.15, carried as observed: the four arms
            //disagreed about *what* to do when the comparison runs off the read, so choosing one
            //is a modelling decision rather than a refactor.
            if (trims_three_prime_ and my_new_offset >= static_cast<int>(query.sequence.size())) {
                continue;
            }
            //The palindrome cannot be longer than the template it is mirrored from.
            if (-deletions > static_cast<int>(previous_str.size())) {
                continue;
            }

            if (trims_three_prime_) {
                previous_str.substr(tmp_str, previous_str.size() + deletions, string::npos);
            } else {
                previous_str.substr(tmp_str, 0, -deletions);
            }
            reverse(tmp_str.begin(), tmp_str.end());
            make_transversions(tmp_str);
            new_str = trims_three_prime_ ? previous_str + tmp_str : tmp_str + previous_str;

            //Count mismatches and add them to the mismatches list. The |k| read positions the
            //palindrome newly occupies run upwards from just past the old 3' end, or from the
            //new 5' end -- the same set, measured from whichever end moved.
            mismatches_vector = previous_mismatches;
            const Seq_Offset first_new_position = trims_three_prime_ ? my_offset + 1 : my_new_offset;
            for (int i = 0; i != -deletions; ++i) {
                const Seq_Offset position = first_new_position + i;
                //Half a guard each, and different halves -- §7.15's table. An internal segment's
                //3' palindrome can run below the read and its 5' one past the end; the arms
                //anchored on the read have neither, which is the defect that table records.
                if (guard_palindrome_positions_
                    and (trims_three_prime_
                                 ? position < 0
                                 : static_cast<std::size_t>(position) >= query.int_sequence.size())) {
                    continue;
                }
                if (not comp_nt_int(tmp_str[i], query.int_sequence.at(position))) {
                    mismatches_vector.push_back(position);
                }
            }
            if (not trims_three_prime_) {
                //A 5' palindrome prepends positions below everything already in the list; a 3'
                //one appends above it, and needs no re-sort.
                sort(mismatches_vector.begin(), mismatches_vector.end());
            }
        }

        scenario.set_sequence_segment(my_seq_type, &new_str, memory_layer_cs);
        scenario.set_offset(my_seq_type, this->event_side, my_new_offset, memory_layer_offset_del);
        scenario.set_mismatches(my_seq_type, &mismatches_vector, memory_layer_mismatches);

        //Get the junction upper bound proba for the span this deletion widens.
        std::optional<double> junction_bound_proba;
        if (junction.resolved()) {
            const Seq_Offset partner_offset = flank_checks_[junction_partner_].offset;
            const int junction_length = trims_three_prime_ ? partner_offset - my_new_offset - 1
                                                           : my_new_offset - partner_offset - 1;
            junction_bound_proba = junction.profile().best_for(junction_length);
            if (not junction_bound_proba) {
                continue; //This means no scenario can lead to a correct solution, would need to be changed for Error models with in/dels
            }
            //With the early stage the slot starts neutral and takes its value below, once that
            //stage has had its chance to stop the enumeration; without it there is one write.
            exploration.downstream_proba_map.set(junction.proba_key(),
                                                 early_prune_stage_ ? 1.0 : *junction_bound_proba,
                                                 junction.memory_layer());
        }

        //Update the mismatches penalty. While another event can still move the opposite end,
        //nothing of this template is settled, so it constrains nothing -- the `//TODO finish
        //this part (compute endogeneous mismatches)` the two D arms carried stays uncomputed.
        if (opposite_end_still_moves_) {
            exploration.downstream_proba_map.set(my_seq_type, 1.0, memory_layer_proba_map_seq);
        } else {
            exploration.downstream_proba_map.set(
                    my_seq_type,
                    accumulation.error_rate->get_err_rate_upper_bound(
                            mismatches_vector.size(), new_str.size() - mismatches_vector.size()),
                    memory_layer_proba_map_seq);
        }

        if (early_prune_stage_) {
            //Stage one: the bound without this realization's own marginal and with the junction
            //still neutral. 4a showed it can only fire where stage two fires too, and that both
            //bounds are monotone along the enumeration -- which is decreasing deletion count --
            //so the `break` is an optimisation rather than a behaviour (§6.14). It is kept
            //because dropping it costs work, and it stays conditional because unifying the two
            //shapes is not free: the arms without it multiply the contribution into the *bound*,
            //which rounds differently from multiplying it into the probability first.
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );
            if (exploration.should_prune(scenario_upper_bound_proba)) {
                break;
            }

            new_scenario_proba *= proba_contribution;
            //Same junction, same distance -- the map was not touched in between, so the value
            //read above still stands.
            if (junction.resolved()) {
                exploration.downstream_proba_map.set(junction.proba_key(), *junction_bound_proba,
                                                     junction.memory_layer());
            }
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );
        } else {
            scenario_upper_bound_proba = exploration.compute_upper_bound(
                new_scenario_proba,
                current_downstream_proba_memory_layers
            );
            new_scenario_proba *= proba_contribution;
            scenario_upper_bound_proba *= proba_contribution;
        }

        if (exploration.should_prune(scenario_upper_bound_proba)) {
            continue;
        }

        // Update context with new probability before proceeding
        scenario.scenario_proba = new_scenario_proba;

        Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
    }
}

/*
 * Keep the mismatches that survive the trim, as a contiguous subrange of an ordered list.
 *
 * The direction is the opposite of the side's name: a 3' deletion keeps the *prefix* of the
 * list and a 5' deletion the *suffix*, because the positions are read coordinates and the end
 * that moved is the one whose side they fall on. Both walks are carried verbatim from the arms
 * they came from, including their signed/unsigned comparison against the new offset.
 */
void Deletion::trim_mismatches(const vector<size_t> &previous_mismatches)
{
    if (previous_mismatches.empty()) {
        mismatches_vector.clear();
        return;
    }

    end_reached = false;
    if (trims_three_prime_) {
        mis_iter = previous_mismatches.end();
        //Get an iterator referring to an actual integer (end() is not dereferenceable)
        --mis_iter;
        while ((*mis_iter) > my_new_offset) {
            if (mis_iter == previous_mismatches.begin()) {
                end_reached = true;
                break;
            }
            --mis_iter;
        }
        if (end_reached) {
            //Clear the vector but the capacity remains the same
            mismatches_vector.clear();
        } else {
            ++mis_iter;
            mismatches_vector.assign(previous_mismatches.begin(), mis_iter);
        }
    } else {
        mis_iter = previous_mismatches.begin();
        while ((*mis_iter) < my_new_offset) {
            ++mis_iter;
            if (mis_iter == previous_mismatches.end()) {
                end_reached = true;
                break;
            }
        }
        if (end_reached) {
            mismatches_vector.clear();
        } else {
            mismatches_vector.assign(mis_iter, previous_mismatches.end());
        }
    }
}

/*
 * This short method performs the iterate operations common to all Rec_event (modify index map and fetch realization probability)
 */
void Deletion::iterate_common(
        forward_list<Event_realization>::const_iterator &iter, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        const Marginal_array_p &model_parameters_point)
{
    proba_contribution = Rec_Event::iterate_common(
        (*iter).index, base_index, base_index_map, model_parameters_point);
}

queue<int> Deletion::draw_random_realization(
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
            const Seq_type target_seq_type = this->target_seq_type;

            switch (target_seq_type) {

            case V_gene_seq:
                if ((*iter).second.value_int >= 0) {
                    constructed_sequences.at(V_gene_seq)
                            .erase(constructed_sequences.at(V_gene_seq).size() - (*iter).second.value_int);
                } else {
                    string &v_gene_seq = constructed_sequences.at(V_gene_seq);
                    gen_tmp_str = v_gene_seq.substr(v_gene_seq.size() + (*iter).second.value_int, string::npos);
                    reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                    make_transversions(gen_tmp_str, false);
                    v_gene_seq += gen_tmp_str;
                }

                break;

            case D_gene_seq:
                switch (this->event_side) {

                case Five_prime:
                    if ((*iter).second.value_int >= 0) {
                        constructed_sequences.at(D_gene_seq).erase(0, (*iter).second.value_int);
                    } else {
                        string &d_gene_seq = constructed_sequences.at(D_gene_seq);
                        gen_tmp_str = d_gene_seq.substr(0, -(*iter).second.value_int);
                        reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                        make_transversions(gen_tmp_str, false);
                        gen_new_str = gen_tmp_str + d_gene_seq;
                        d_gene_seq = gen_new_str;
                    }

                    break;

                case Three_prime:
                    if ((*iter).second.value_int >= 0) {
                        constructed_sequences.at(D_gene_seq)
                                .erase(constructed_sequences.at(D_gene_seq).size() - (*iter).second.value_int);
                    } else {
                        string &d_gene_seq = constructed_sequences.at(D_gene_seq);
                        gen_tmp_str = d_gene_seq.substr(d_gene_seq.size() + (*iter).second.value_int, string::npos);
                        reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                        make_transversions(gen_tmp_str, false);
                        d_gene_seq += gen_tmp_str;
                    }

                    break;

                default:
                    break;
                }
                break;
            case J_gene_seq:
                if ((*iter).second.value_int >= 0) {
                    constructed_sequences.at(J_gene_seq).erase(0, (*iter).second.value_int);
                } else {
                    string &j_gene_seq = constructed_sequences.at(J_gene_seq);
                    gen_tmp_str = j_gene_seq.substr(0, -(*iter).second.value_int);
                    reverse(gen_tmp_str.begin(), gen_tmp_str.end());
                    make_transversions(gen_tmp_str, false);
                    gen_new_str = gen_tmp_str + j_gene_seq;
                    j_gene_seq = gen_new_str;
                }

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

void Deletion::write2txt(ofstream &outfile)
{
    write2txt_legacy(outfile);
}

void Deletion::write2txt_legacy(ofstream &outfile)
{
    outfile << "#Deletion;" << event_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Deletion::write2txt_v2(ofstream &outfile)
{
    outfile << "#Deletion;Undefined_gene;" << seq_type << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Deletion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const Events_map &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        SafetyMatrix &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    //A deletion trims one end of one segment, so a model that does not say which end is not a
    //model this event can run in. Rejected here rather than at the scenario node the four-arm
    //body rejected it at: the generic body has no side to switch on, and a topology error
    //belongs before inference rather than a million nodes into it.
    if (this->event_side != Five_prime and this->event_side != Three_prime) {
        throw invalid_argument("Deletion " + this->get_name()
                               + ": a deletion must trim the 5' or the 3' end of its segment, "
                                 "and this one names neither");
    }

    //TODO change this and the usage of int_value_and_index
    int_value_and_index.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter = (*this).event_realizations.begin();
         iter != (*this).event_realizations.end(); ++iter) {
        int_value_and_index.push_front((*iter).second);
    }
    int_value_and_index.sort(del_numb_compare);

    //Everything topological is settled here, so iterate() reads answers instead of switching on
    //the gene again. The topology has to be there to be read: a segment the model registered
    //but left out of the ordering has no neighbours, no read ends and no pairs, so every flag
    //below would be answered by silence rather than by the model.
    if (not constructed_sequences.registry().contains(
                constructed_sequences.registry().name(this->seq_type_id))) {
        throw invalid_argument("Deletion " + this->get_name() + ": its segment "
                               + constructed_sequences.registry().name(this->seq_type_id)
                               + " is not in the model's 5'->3' ordering, so there is no side "
                                 "for it to trim towards");
    }

    //The eight `*_min_del` / `*_max_del` scalars and the sixty lines that filled them by
    //looking up the V 3', D 5', D 3' and J 5' deletion events by hand collapse to one query per
    //checked partner (§2.1); a topology with more than one deletion per end needs no new code
    //for it, because the deltas sum.
    const SeqTypeRegistry &registry = constructed_sequences.registry();

    trims_three_prime_ = (this->event_side == Three_prime);

    //A segment anchored on an end of the read behaves differently in three places, and all
    //three used to be spelled "V or J". Read off the ordering, so a tandem D1/D2 pair gets the
    //internal behaviour without either of them being named -- the same resolution B11a gave
    //`Gene_choice`.
    const bool at_read_end = registry.left_neighbor(this->seq_type_id) == kNoSeqType
                          or registry.right_neighbor(this->seq_type_id) == kNoSeqType;
    keep_one_nucleotide_ = at_read_end;
    early_prune_stage_ = at_read_end;
    guard_palindrome_positions_ = not at_read_end;
    discard_offset_outside_read_ = not at_read_end and this->event_side == Five_prime;

    //Is anything still going to move the other end of my segment? Asked through A0 rather than
    //by looking up "the D deletion on the opposite side", which is what makes it answerable in
    //a topology with two Ds. The gene choice that created the segment is always processed by
    //now -- a deletion runs after its own gene choice by construction -- so the only events
    //that can answer `Modifies` here are deletions still to come.
    const Seq_side opposite_side = trims_three_prime_ ? Five_prime : Three_prime;
    opposite_end_still_moves_ = false;
    for (const auto &[key, event] : events_map) {
        (void)key;
        if (not event or processed_events.count(event->get_name()) != 0) {
            continue;
        }
        if (event->get_name() == this->get_name()) {
            continue;
        }
        if (event->get_offset_role(this->seq_type_id, opposite_side) == OffsetRole::Modifies) {
            opposite_end_still_moves_ = true;
        }
    }

    seq_offsets.request_layer(this->seq_type_id, this->event_side);
    memory_layer_offset_del = seq_offsets.claimed_layer(this->seq_type_id, this->event_side);
    mismatches_list.request_layer(this->seq_type_id);
    this->memory_layer_mismatches = mismatches_list.claimed_layer(this->seq_type_id);
    constructed_sequences.request_layer(this->seq_type_id);
    this->memory_layer_cs = constructed_sequences.claimed_layer(this->seq_type_id);

    //The neighbours this deletion is checked against: the other gene segments **on the side it
    //trims**, 5' to 3'. The other side of the segment is not moving, so nothing there can newly
    //collide -- which is why each of the four arms checked one or two partners and never all
    //three. The candidate list is still the legacy three by name; parent plan B9 step 3, which
    //puts the whole topology in the registry, is what replaces this literal with
    //registry.ordering().
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
        if ((position > my_position) != trims_three_prime_) {
            continue;
        }
        const Seq_type_String &partner_name = kGeneSegments[position];
        FlankCheck check;
        check.partner_id = registry.id(partner_name);
        check.partner_side = trims_three_prime_ ? Five_prime : Three_prime;
        check.safety_cell = safety_set.cell(this->seq_type_id, check.partner_id);
        check.partner_chosen =
                EventUtils::check_gene_choice(partner_name, events_map, processed_events).chosen;
        check.partner_delta = JunctionGeometry::legacy_offset_delta(
                check.partner_id, check.partner_side, events_map, processed_events);
        flank_checks_.push_back(check);
    }

    //Claimed only for a neighbour that has been placed, which is what the four arms' `if
    //(x_chosen)` guards said: an unplaced neighbour has no offset to compare against and no
    //verdict to record. One claim per *row*, not per check -- a 3' deletion's partners all sit
    //in the row of its own segment and share a word, so a second claim would leave a layer
    //unwritten, and an unwritten layer is unreadable by design (§7.9).
    for (std::size_t i = 0; i != flank_checks_.size(); ++i) {
        FlankCheck &check = flank_checks_[i];
        if (not check.partner_chosen) {
            continue;
        }
        std::size_t earlier = 0;
        while (earlier != i
               and (not flank_checks_[earlier].partner_chosen
                    or flank_checks_[earlier].safety_cell.row != check.safety_cell.row)) {
            ++earlier;
        }
        if (earlier != i) {
            check.safety_layer = flank_checks_[earlier].safety_layer;
        } else {
            safety_set.request_layer(check.safety_cell);
            check.safety_layer = safety_set.claimed_layer(check.safety_cell);
        }
        check.partner_offset_layer = seq_offsets.claimed_layer(check.partner_id, check.partner_side);
    }

    downstream_proba_map.request_layer(this->seq_type_id);
    memory_layer_proba_map_seq = downstream_proba_map.claimed_layer(this->seq_type_id);

    //The junction this deletion widens: the one between its trimmed end and the **nearest
    //placed** neighbour on that side. Nearest, because a nearer placed segment makes the
    //further junction somebody else's to bound -- which is what the `if (d_chosen) … else if
    //(j_chosen)` pairs said, once per arm.
    junction_partner_ = -1;
    for (std::size_t i = 0; i != flank_checks_.size(); ++i) {
        if (not flank_checks_[i].partner_chosen) {
            continue;
        }
        junction_partner_ = static_cast<int>(i);
        if (trims_three_prime_) {
            break; //the first placed neighbour to the 3' side
        }
        //...and to the 5' side the last one, which the loop reaches by not breaking.
    }
    if (junction_partner_ >= 0) {
        const SeqTypeId partner_id = flank_checks_[junction_partner_].partner_id;
        const SegmentSpan span = trims_three_prime_
                                         ? SegmentSpan::gap(this->seq_type_id, partner_id)
                                         : SegmentSpan::gap(partner_id, this->seq_type_id);
        const SeqTypeId junction_key = legacy_junction_of(span);
        downstream_proba_map.request_layer(junction_key);
        resolve_junction(span, junction_key, downstream_proba_map.claimed_layer(junction_key));
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list,
                                      seq_offsets, index_map);
}

namespace {

/// The inclusive deletion range, taken straight from the realization set.
///
/// Deliberately not read off len_min / len_max: those are accumulated by an `if / else if`
/// over an unordered map, so a strictly ascending iteration order never reaches the second
/// branch and leaves one bound at its INT16 sentinel. Latent today -- the hash order for
/// realistic key sets happens to cooperate -- but the failure would be silent and severe.
/// See docs/ITERATE_GENERIC_REWRITE_PLAN.md section 7.4.
std::pair<int, int> deletion_range(const std::unordered_map<std::string, Event_realization> &realizations)
{
    if (realizations.empty()) {
        return {0, 0};
    }
    int min_del = std::numeric_limits<int>::max();
    int max_del = std::numeric_limits<int>::min();
    for (const auto &[name, realization] : realizations) {
        (void)name;
        min_del = std::min(min_del, realization.value_int);
        max_del = std::max(max_del, realization.value_int);
    }
    return {min_del, max_del};
}

} // namespace

OffsetDelta Deletion::get_offset_delta_bounds(SeqTypeId type_id, Seq_side side) const
{
    if (type_id != this->seq_type_id || side != this->event_side) {
        return {};
    }
    const auto [min_del, max_del] = deletion_range(this->event_realizations);

    //A 3' end retreats as nucleotides are removed, a 5' end advances; a negative deletion
    //(a palindromic insertion) moves it the other way. This sign flip is the one the
    //current code spells out at each of its eight comparison sites.
    if (side == Three_prime) {
        return {-max_del, -min_del};
    }
    return {min_del, max_del};
}

LengthContribution Deletion::get_length_contribution(SeqTypeId type_id) const
{
    if (type_id != this->seq_type_id) {
        return {};
    }
    const auto [min_del, max_del] = deletion_range(this->event_realizations);
    //Removing `k` nucleotides shortens the segment by `k`, whichever end they come from.
    return {-max_del, -min_del};
}

SeqConstructionRole Deletion::get_seq_construction_role(SeqTypeId type_id) const
{
    return type_id == this->seq_type_id ? SeqConstructionRole::Modifies : SeqConstructionRole::None;
}

OffsetRole Deletion::get_offset_role(SeqTypeId type_id, Seq_side side) const
{
    if (type_id != this->seq_type_id || side != this->event_side) {
        return OffsetRole::None;
    }
    return OffsetRole::Modifies;
}

void Deletion::add_to_marginals(long double scenario_proba, Marginal_array_p &updated_marginals) const
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
                //Nothing to do
            } else if ((*iter) == '9') {
                //Nothing to do
            } else if ((*iter) == '14') {
                //Nothing to do
            } else {
                throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq " + init_sequence
                                    + " in make_transversions()");
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
                throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq " + init_sequence
                                    + " in make_transversions()");
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
            //Nothing to do
        } else if ((*iter) == 9) {
            //Nothing to do
        } else if ((*iter) == 14) {
            //Nothing to do
        } else {

            string error_str("Unknown int nucleotide " + to_string((*iter)) + " in seq ");
            for (Int_Str::iterator jiter = init_sequence.begin(); jiter != init_sequence.end(); ++jiter) {
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

bool Deletion::affects_length_of(SegmentSpan span) const
{
    //A deletion moves an anchor's boundary away from where it was created, so it widens the
    //span that boundary bounds rather than shortening the anchor -- see the frame in section
    //2.5 of docs/ITERATE_GENERIC_REWRITE_PLAN.md.
    //
    //The generic form of this table is "(target, side) is an inward-facing endpoint of the
    //span", but get_deletion_effective_junctions() is side-insensitive for V and J, so the two
    //differ for a hypothetical V 5' or J 3' deletion. Keeping the table verbatim here holds S4a
    //bitwise; reconciling them is S4b's.
    const Seq_type junction = legacy_junction_of(span);
    const auto effective_junctions =
            get_deletion_effective_junctions(this->target_seq_type, this->event_side);
    return find(effective_junctions.begin(), effective_junctions.end(), junction) != effective_junctions.end();
}

int Deletion::length_delta(const Event_realization &realization) const
{
    //Negative: a deletion moves an anchor's boundary away from where it was created, which
    //widens the span it bounds rather than shortening the anchor (§2.5's frame).
    return -realization.value_int;
}

/*
 * A deletion reads exactly one junction: the one on the side it trims. The table this replaces
 * named two for a V or a J deletion -- {VD, VJ} and {VJ, DJ} -- of which one was always dead,
 * because the consumption site picks by `d_chosen` and the other branch is never taken. On the
 * TRB corpus that dead V->J fold was 176 ms of the 706 ms sweep, about a quarter (§6.10
 * finding 3). Resolving the junction in initialize_event() removes the choice, so nothing is
 * built that nothing reads.
 */
void Deletion::resolve_junction(SegmentSpan span, SeqTypeId proba_key, int memory_layer)
{
    junction_bound(this->event_side == Five_prime ? kLeftJunction : kRightJunction)
            .resolve(span, proba_key, memory_layer, JunctionBound::Fold::Yes);
}
