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
        // Not a leaf node - recursively call next event's iterate_wrap_up
        exploration.next_event_ptr_arr.get()[this->event_index]->iterate(
            query,
            model,
            scenario,
            exploration,
            accumulation
        );
    } else {
        // Leaf node - complete scenario and accumulate

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

        // Check pruning threshold
        if (not exploration.should_prune(scenario.scenario_error_w_proba)) {
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
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
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
            memory_and_offsets.emplace_front(event_identitfier, index_map.current_layer(event_identitfier),
                                             (*iter).second);
        }
    }

    //Snapshot, not a view: the map's layers move as the traversal proceeds.
    const auto downstream_layers = downstream_proba_map.current_layers();
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

void Rec_Event::iterate_initialize_Len_proba(Seq_type considered_junction, std::map<int, double> &length_best_proba_map,
                                             std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                             double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                             Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences) const
{
    int seq_len = 0;
    this->iterate_initialize_Len_proba(considered_junction, length_best_proba_map, model_queue, scenario_proba,
                                       model_parameters_point, base_index_map, constructed_sequences, seq_len);
}

/*
 * Called when iterating over all possible scenarios during initialization
 * Fills up the length-max_proba_bound for a given junction , and links the call to iterate_initialize_len_proba for two events
 *
 * TODO constructed sequences should not be used but it is useful to compute the dinucl contribution
 */
void Rec_Event::iterate_initialize_Len_proba_wrap_up(Seq_type considered_junction,
                                                     std::map<int, double> &length_best_proba_map,
                                                     std::queue<std::shared_ptr<Rec_Event>> model_queue,
                                                     double scenario_proba,
                                                     const Marginal_array_p &model_parameters_point,
                                                     Index_map &base_index_map,
                                                     Seq_type_str_p_map &constructed_sequences, int seq_len) const
{

    if (not model_queue.empty()) {
        std::shared_ptr<Rec_Event> next_event_p = model_queue.front();
        model_queue.pop();
        //TODO fix this and find a way not to loop over all events
        //if(next_event_p->has_effect_on(considered_junction)){
        // Explore realizations of this event
        next_event_p->iterate_initialize_Len_proba(considered_junction, length_best_proba_map, model_queue,
                                                   scenario_proba, model_parameters_point, base_index_map,
                                                   constructed_sequences, seq_len);
        //}
        //else{
        // If this event has no effect on the junction skip it using a recursive call
        //next_event_p->iterate_initialize_Len_proba_wrap_up(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
        //}
    } else {
        // When all events with an effect on the junction have been processed update the length-proba map
        if (length_best_proba_map.count(seq_len) > 0) {
            if (scenario_proba > length_best_proba_map.at(seq_len)) {
                //Keep the best proba for each length
                length_best_proba_map.at(seq_len) = scenario_proba;
            }
        } else {
            length_best_proba_map[seq_len] = scenario_proba;
        }
    }
}
