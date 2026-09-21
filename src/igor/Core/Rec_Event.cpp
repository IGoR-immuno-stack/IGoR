/*
 * Rec_Event.cpp
 *
 *  Created on: 3 nov. 2014
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
 *
 */

#include <igor/Core/Rec_Event.h>
#include <igor/Core/gene_to_seqtype_migr.h>
#include <igor/Core/Counter.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/Scenario.h>  // For Scenario view construction
#include <igor/Core/BoundTightness.h>

#include <cassert>
#include <iostream>

using namespace std;


Rec_Event::Rec_Event(Gene_class gene, Seq_side side)
    : priority(0),
      event_class(gene),
      event_side(side),
      name("Undefined_event_name"),
      len_min(INT16_MAX),
      len_max(INT16_MIN),
      type(Undefined_t),
      event_index(INT16_MIN),
      updated(false),
      fixed(false),
      current_realizations_index_vec(vector<int>()),
      scenario_downstream_upper_bound_proba(-1),
      event_upper_bound_proba(-1),
      scenario_upper_bound_proba(-1),
      current_realization_index(nullptr)
{
} //FIXME why does this exist? anyway fix initilization

Rec_Event::Rec_Event(Gene_class gene, Seq_side side, unordered_map<string, Event_realization> &realizations)
    : Rec_Event(gene, side)
{
    this->event_realizations = realizations;
}

Rec_Event::Rec_Event() : Rec_Event(Undefined_gene, Undefined_side) { }

//TODO see this later
/*
Rec_Event::Rec_Event(list<Event_realization> realizations_list){
	for(list<Event_realization>::const_iterator iter = realizations_list.begin() ; iter!=realizations_list.end() ; iter++){
		this->add_realization((*iter));
	}
}


Rec_Event::Rec_Event(list<Event_realization> realization_list, int new_priority ) : Rec_Event(realization_list) {
	this->priority = new_priority;
}
*/

Rec_Event::~Rec_Event()
{
    // TODO Auto-generated destructor stub
}

bool Rec_Event::operator==(const Rec_Event &other) const
{
    if (this->get_type() != other.get_type())
        return 0;
    if (this->event_class != other.event_class)
        return 0;
    if (this->event_side != other.event_side)
        return 0;
    if (this->priority != other.priority)
        return 0;
    if (this->event_realizations.size() != other.event_realizations.size())
        return 0;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if (other.event_realizations.count((*iter).first) != 1)
            return 0;
    }
    return 1;
}

void Rec_Event::update_event_name()
{
    this->name = string() + this->type + string("_") + this->event_class + string("_") + this->event_side
            + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

Rec_Event_name Rec_Event::get_v2_name() const
{
    if (seq_type.empty())
        return name;
    return string() + this->type + string("_") + this->event_class + string("_") + seq_type
           + string("_") + this->event_side
           + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

Rec_Event_name Rec_Event::get_legacy_name() const
{
    return string() + this->type + string("_") + this->event_class + string("_") + this->event_side
           + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

void Rec_Event::add_realization(const Event_realization &realization)
{
    this->event_realizations.insert(make_pair((realization).name, (realization)));
    this->update_event_name();
}

bool Rec_Event::set_priority(int new_priority)
{
    this->priority = new_priority;
    this->update_event_name();
    return 1;
}

int Rec_Event::size() const
{
    return event_realizations.size();
}

void Rec_Event::set_event_identifier(size_t identifier)
{
    this->event_index = identifier;
}

int Rec_Event::get_event_identifier() const
{
    return event_index;
}

double Rec_Event::iterate_common(int realization_index, int base_index,
                                 Index_map &base_index_map,
                                 const Marginal_array_p &model_parameters) {
    // Update parent realization tracking for marginal array indexing of child events
    update_parent_tracking(realization_index, base_index_map);

    // Return marginal probability for this realization
    return model_parameters[base_index + realization_index];
}

/**
 * @brief Context-based iterate_wrap_up with Counter interface
 *
 * Handles leaf-node scenario completion:
 * - Computes error-weighted probability
 * - Creates Scenario view from ScenarioContext
 * - Calls counters and marginal accumulation
 *
 * This replaces the legacy adapter pattern with direct context usage.
 */
void Rec_Event::iterate_wrap_up(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    if (exploration.next_event_ptr_arr.get()[this->event_index]) {
        //One node of the scenario tree: this event's realization, and everything the walk
        //explores under it. The pair brackets the descent so the instrumentation can compare the
        //bound that admitted this node against the best any leaf below it turned out to reach --
        //which is the quantity pruning acts on, and the one the leaf ratio cannot show. Both
        //calls compile to nothing unless IGOR_BOUND_INSTRUMENTATION is defined.
        //The name labels this depth's row in the per-event decomposition. Taken by value and
        //copied straight into a static table: a pointer into a per-thread event copy would
        //dangle long before the report runs.
        BoundTightness::enter(this->scenario_upper_bound_proba, this->get_name().c_str());

        // Not a leaf node - recursively call next event's iterate_wrap_up
        exploration.next_event_ptr_arr.get()[this->event_index]->iterate(
            query,
            model,
            scenario,
            exploration,
            accumulation
        );

        BoundTightness::leave();
    } else {
        // Leaf node - complete scenario and accumulate

#ifndef NDEBUG
        // Every position of every segment must be determined by now: int_undefined is a
        // placeholder that e.g. Insertion leaves for its Dinucl_markov, never a value a consumer can
        // interpret. Checked here because this is the one boundary where the invariant has to
        // hold, and only under assertions -- the walk is linear in the scenario's length, and
        // the default build defines NDEBUG, so release pays nothing.
        if (const SeqTypeId unfilled = first_unfilled_segment(scenario.constructed_sequences);
            unfilled != kNoSeqType) {
            std::cerr << "Scenario leaf reached with unfilled nucleotides in "
                      << scenario.constructed_sequences.registry().name(unfilled)
                      << std::endl;
            assert(false && "unfilled nucleotide in a completed scenario");
        }
#endif

        // Compute error-weighted probability using error_rate
        // TODO (Future): Consider changing compute_scenario_error_probability() to void return
        //                Method already assigns scenario.scenario_error_w_proba internally,
        //                making this assignment redundant. Void return would match
        //                Rec_Event::iterate() pattern (contexts mutated, no return value).
        //                See docs/ITERATE_REFACTORING_PLAN.md "Future Refactoring Opportunities"
        scenario.scenario_error_w_proba = accumulation.error_rate->compute_scenario_error_probability(
            query,
            model,
            scenario,
            exploration
        );

        //How far above the realized probability the bound that let this scenario through sits.
        //`scenario_upper_bound_proba` is the value *this* event computed before handing off, and
        //at a leaf this event is the last one in the queue -- so it is the bound the whole
        //descent ended on. Compiled out unless IGOR_BOUND_INSTRUMENTATION is defined; see
        //BoundTightness.h and section 6.10 of docs/ITERATE_GENERIC_REWRITE_PLAN.md.
        //
        //Recorded before the threshold test, not after: a scenario the threshold discards is
        //still one the bound admitted, and leaving those out would measure only the tail that
        //survived.
        BoundTightness::record(this->scenario_upper_bound_proba, scenario.scenario_error_w_proba,
                               this->get_name().c_str());

        // Check pruning threshold. `is_below_threshold` rather than `should_prune`: the value
        // tested here is what the scenario realized, not a bound on what it might still become,
        // and the decision discards one completed scenario rather than a subtree. The
        // instrumentation counts subtree prunings to tell a barren node that a tighter bound
        // could have deleted from one that died on geometry, and this test belongs to neither.
        if (not exploration.is_below_threshold(scenario.scenario_error_w_proba)) {
            // Update best scenario probability if needed
            exploration.update_max_prob(scenario.scenario_error_w_proba);

            // Create Scenario view and call counters with context interface
            Scenario scenario_view(scenario);

            for (auto& [counter_id, counter] : accumulation.counters) {
                counter->count_scenario(scenario_view, query, model);
            }

            // Accumulate marginals for non-fixed events
            for (const auto& [event_key, event_ptr] : model.events_map) {
                if (!event_ptr->is_fixed()) {
                    event_ptr->add_to_marginals(scenario.scenario_error_w_proba, accumulation.updated_marginals);
                }
            }
        }
    }
}

void Rec_Event::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const Events_map &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        SafetyMatrix &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    initialize_event_common(processed_events, offset_map, downstream_proba_map, index_map);
}

void Rec_Event::initialize_event_common(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map,
        Index_map &index_map)
{
    //No action performed on the event by default if the method is not overloaded
    //Need to call Rec_Event::initialize_event() to apply these common actions when the method is overloaded
    current_realizations_index_vec.push_back(-1);

    if (offset_map.count(this->get_name()) != 0) {
        const vector<pair<shared_ptr<const Rec_Event>, int>> &offset_vector = offset_map.at(this->get_name());
        for (vector<pair<shared_ptr<const Rec_Event>, int>>::const_iterator iter = offset_vector.begin();
             iter != offset_vector.end(); ++iter) {
            //Request memory layer
            int event_identitfier = (*iter).first->get_event_identifier();
            index_map.request_layer(event_identitfier);
            memory_and_offsets.emplace_front(event_identitfier, index_map.claimed_layer(event_identitfier),
                                             (*iter).second);
        }
    }

    //Snapshot, not a view: the map's layers move as the traversal proceeds.
    const auto downstream_layers = downstream_proba_map.claimed_layers();
    current_downstream_proba_memory_layers.assign(downstream_layers.begin(), downstream_layers.end());

    processed_events.emplace(this->name);
}

void Rec_Event::ind_normalize(Marginal_array_p &marginal_array_p, size_t base_index) const
{
    long double sum_marginals = 0;
    for (int i = 0; i != this->size(); ++i) {
        sum_marginals += marginal_array_p[base_index + i];
    }
    if (sum_marginals != 0) {
        for (int i = 0; i != this->size(); ++i) {
            marginal_array_p[base_index + i] /= sum_marginals;
        }
    }
}

void Rec_Event::set_crude_upper_bound_proba(size_t base_index, size_t event_size, Marginal_array_p &marginal_array_p)
{

    double max_proba = 0;
    for (size_t i = 0; i != event_size; ++i) {
        if (marginal_array_p[base_index + i] > max_proba) {
            max_proba = marginal_array_p[base_index + i];
        }
    }
    this->event_upper_bound_proba = max_proba;
}

void Rec_Event::set_upper_bound_proba(double proba)
{
    this->event_upper_bound_proba = proba;
}

/**
 * Does nothing since in general events will not need to perform any operation on the marginal probabilities
 */
void Rec_Event::update_event_internal_probas(const Marginal_array_p &marginal_array,
                                             const unordered_map<Rec_Event_name, int> &index_map)
{
    //Do nothing
}

/*
 * This method initialize the scenario probability upper bound for each event
 * The point is to compute the upper bound probability (given the model) of the scenario for each event
 * This allows to discard scenarios with too low probability at early stages
 */
void Rec_Event::initialize_crude_scenario_proba_bound(
        double &downstream_proba_bound, forward_list<double *> &updated_proba_list,
        const Events_map &events_map)
{
    this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
    this->updated_proba_bounds_list = updated_proba_list;
    if (!this->is_updated()) {
        downstream_proba_bound *= this->event_upper_bound_proba;
    } else {
        throw logic_error("Updated events should overload Rec_event::initialize_scenario_proba_bound()");
    }
}

/*
 * Description??
 */
double *Rec_Event::get_updated_ptr()
{
    throw logic_error("Updated events should overload Rec_event::get_updated_ptr()");
}

/*
 * Updates the value of scenario_upper_bound_proba according to the error weighted scenario and the upper bound of downstream scenarios
 */
void Rec_Event::compute_crude_upper_bound_scenario_proba(double &tmp_err_w_proba)
{
    scenario_upper_bound_proba = tmp_err_w_proba * scenario_downstream_upper_bound_proba;
    for (forward_list<double *>::const_iterator iter = updated_proba_bounds_list.begin();
         iter != updated_proba_bounds_list.end(); ++iter) {
        scenario_upper_bound_proba *= (*(*iter));
    }
}

void Rec_Event::iterate_initialize_Len_proba(SegmentSpan span, SpanProfile &profile,
                                             const SpanParticipants &participants,
                                             double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                             Index_map &base_index_map, SpanAccumulator &lengths) const
{
    int seq_len = 0;
    //This overload is the traversal's entry point, and `this` is not in `participants` -- that
    //list holds the events *after* this one -- so it has to apply the same test to itself that
    //built the list. It is not vacuous: Gene_choice(V) opens the VD span traversal but
    //contributes nothing to it.
    if (this->participates_in_span(span)) {
        this->iterate_initialize_Len_proba(span, profile, participants, 0, scenario_proba,
                                           model_parameters_point, base_index_map, lengths, seq_len);
    } else {
        this->iterate_initialize_Len_proba_wrap_up(span, profile, participants, 0, scenario_proba,
                                                   model_parameters_point, base_index_map, lengths, seq_len);
    }
}

/*
 * One body for all four event kinds. Called only for an event the filter has already found to
 * participate in `span`, so the two branches below are "enumerates its realizations" and
 * "contributes a probability factor without enumerating" -- which is exactly the
 * affects_length_of / affects_proba_of split.
 */
void Rec_Event::iterate_initialize_Len_proba(SegmentSpan span, SpanProfile &profile,
                                             const SpanParticipants &participants, std::size_t cursor,
                                             double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                             Index_map &base_index_map, SpanAccumulator &lengths, int &seq_len) const
{
    //A local, not the subclasses' `mutable int base_index`: this body is shared, and each of the
    //four declares its own. Safe because every reader of that member sets it first in the same
    //call chain -- iterate() and initialize_event() both do -- so the fold never had to publish it.
    base_index_map.set_current_layer(this->event_index, 0);
    const int span_base_index = base_index_map.get(this->event_index);

    if (not this->affects_length_of(span)) {
        //Dinucl_markov: no realization of its own contributes length, and its p^L factor reads a
        //length some upstream creator already published. One factor, one recursive call.
        double contributed_proba = scenario_proba * this->span_proba_factor(span, lengths);
        this->iterate_initialize_Len_proba_wrap_up(span, profile, participants, cursor, contributed_proba,
                                                   model_parameters_point, base_index_map, lengths, seq_len);
        return;
    }

    //Only a segment's creator publishes its length, so each key has one writer per path and a
    //published value is a real segment size. A Deletion contributes its negative delta to the
    //span total without touching the accumulator.
    const bool publishes_length =
            this->get_seq_construction_role(this->seq_type_id) == SeqConstructionRole::Creates;

    for (std::unordered_map<std::string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        const Event_realization &realization = iter->second;

        //Get the max proba for this realization (in case the event is child of another).
        //This maxᵢ is what R6 replaces with a max taken jointly over a conditioned clique.
        double real_max_proba = 0;
        for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
            if (model_parameters_point[span_base_index + realization.index + i * this->size()] > real_max_proba) {
                real_max_proba = model_parameters_point[span_base_index + realization.index + i * this->size()];
            }
        }

        const int delta = this->length_delta(realization);
        if (publishes_length) {
            lengths.set(this->seq_type_id, delta);
        }

        this->iterate_initialize_Len_proba_wrap_up(span, profile, participants, cursor,
                                                   scenario_proba * real_max_proba, model_parameters_point,
                                                   base_index_map, lengths, seq_len + delta);
    }
}

void Rec_Event::iterate_initialize_Len_proba_wrap_up(SegmentSpan span, SpanProfile &profile,
                                                     const SpanParticipants &participants, std::size_t cursor,
                                                     double scenario_proba,
                                                     const Marginal_array_p &model_parameters_point,
                                                     Index_map &base_index_map, SpanAccumulator &lengths,
                                                     int seq_len) const
{
    //The events that neither change this span's length nor contribute a probability factor to it
    //were dropped when `participants` was built, so descending is one index step: the depth stays
    //proportional to the number of contributing events rather than to the model size, and nothing
    //is copied on the way down.
    if (cursor < participants.size()) {
        // Explore realizations of this event
        participants[cursor]->iterate_initialize_Len_proba(span, profile, participants, cursor + 1, scenario_proba,
                                                           model_parameters_point, base_index_map, lengths, seq_len);
    } else {
        // Every event contributing to this span has chosen, so this path reaches `seq_len` with
        // `scenario_proba`. SpanProfile::record() keeps the better of that and what is already
        // stored, in one descent rather than the count/at/operator[] trio this replaced.
        profile.record(seq_len, scenario_proba);
    }
}

/*
 * The four per-subclass overrides, minus their enum-to-member switches. Everything topological
 * was settled in initialize_event(), so all that is left is "fold each junction I read".
 */
void Rec_Event::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                           const Marginal_array_p &model_parameters_point, Index_map &base_index_map)
{
    //Scoped to one fold and reset between junctions: an entry is a length published by the
    //segment's creator on the current path, and the paths of two junctions share nothing.
    SpanAccumulator lengths(legacy_seq_type_registry().total_count());

    //Flatten the queue of downstream events once, here, instead of copying it at every node of
    //every fold. `model_queue` is left as the caller gave it: the junction loop below reads the
    //same suffix for each of this event's junctions.
    SpanParticipants downstream;
    downstream.reserve(model_queue.size());
    for (queue<shared_ptr<Rec_Event>> remaining = model_queue; not remaining.empty(); remaining.pop()) {
        downstream.push_back(remaining.front().get());
    }

    for (JunctionBound &bound : junction_bounds_) {
        if (not bound.resolved() or not bound.folded()) {
            continue;
        }

        //Which events participate depends only on the span and on the events themselves, both
        //fixed for the whole fold, so the filter runs once per junction rather than per node.
        SpanParticipants participants;
        participants.reserve(downstream.size());
        for (const Rec_Event *const downstream_event : downstream) {
            if (downstream_event->participates_in_span(bound.span())) {
                participants.push_back(downstream_event);
            }
        }

        bound.mutable_profile().clear();
        lengths.reset();
        double init_proba = 1.0;
        this->iterate_initialize_Len_proba(bound.span(), bound.mutable_profile(), participants, init_proba,
                                           model_parameters_point, base_index_map, lengths);
    }

    this->finalize_Len_proba_bound(model_parameters_point, base_index_map);
}

void Rec_Event::adopt_Len_proba_bound(const Rec_Event &source)
{
    for (std::size_t slot = 0; slot != kJunctionSlotCount; ++slot) {
        JunctionBound &bound = junction_bounds_[slot];
        if (not bound.resolved() or not bound.folded()) {
            continue;
        }
        bound.mutable_profile() = source.junction_bounds_[slot].profile();
    }

    this->adopt_finalized_Len_proba_bound(source);
}
