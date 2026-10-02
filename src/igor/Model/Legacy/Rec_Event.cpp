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

#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Counter.h>
#include <igor/Model/Legacy/EventUtils.h>
#include <igor/Model/Legacy/Scenario.h>  // For Scenario view construction
#include <igor/Model/Legacy/BoundTightness.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>


namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

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
      fixed(false),
      current_realizations_index_vec(vector<int>()),
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

/**
 * The side token a generated name carries. Shared, because Insertion and Dinucl_markov override
 * update_event_name() to derive the class token from their ins_seq_type.
 *
 * It is the event's own side, except for DinucMarkov, where it is always Undefined_side. Two
 * reasons, one of substance and one of robustness.
 *
 * Substance: the side does not identify a DinucMarkov event. Its seq_type already says which
 * junction it belongs to, and event_side only carries a direction. Model_Parms is already
 * explicit about this — get_events_map() keys DinucMarkov with Undefined_side, and
 * write2txt_legacy() writes Undefined_side in the event line whatever the in-memory side is.
 * Names now follow the same rule instead of contradicting it.
 *
 * Robustness: the name used to depend on the ORDER of the setters. set_priority() refreshes
 * the name, and read_model_parms() sets the side after the priority, so a DinucMarkov kept a
 * name saying Undefined_side while its side was Three_prime — the truth by accident. Anything
 * that called update_event_name() once more produced a different name for the same event, and
 * that name is the key of Model_marginals::get_index_map(), of Model_Parms::edges, and a
 * column header of the scenario and generation outputs. Making the rule explicit here is what
 * lets set_event_side() refresh the name like every other setter.
 */
Seq_side Rec_Event::name_side(Event_type type, Seq_side side)
{
    return (type == Event_type::Dinuclmarkov_t) ? Undefined_side : side;
}

void Rec_Event::update_event_name()
{
    this->name = string() + this->type + string("_") + this->event_class + string("_")
            + name_side(this->type, this->event_side)
            + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

nlohmann::json Rec_Event::to_json() const
{
    nlohmann::json out;
    out["type"] = string() + this->type;   // "GeneChoice", "Deletion", "Insertion", "DinucMarkov"
    out["gene_class"] = to_string(this->event_class);
    out["seq_type"] = this->seq_type;
    out["side"] = to_string(this->event_side);
    out["priority"] = this->priority;
    out["nickname"] = this->nickname;

    nlohmann::json realizations = nlohmann::json::array();
    for (const auto &entry : this->event_realizations) {
        const Event_realization &real = entry.second;
        nlohmann::json r;
        r["index"] = real.index;
        r["name"] = real.name;
        // INT16_MAX is the "no integer value" sentinel the gene-choice events carry, and an
        // empty value_str is the "no sequence" one. Emitting either would not round-trip.
        if (real.value_int != INT16_MAX)
            r["value_int"] = real.value_int;
        if (!real.value_str.empty())
            r["value_str"] = real.value_str;
        realizations.push_back(std::move(r));
    }
    std::sort(realizations.begin(), realizations.end(),
              [](const nlohmann::json &a, const nlohmann::json &b) {
                  return a.at("index").get<int>() < b.at("index").get<int>();
              });
    out["realizations"] = std::move(realizations);
    return out;
}

Rec_Event_name Rec_Event::get_v2_name() const
{
    if (seq_type.empty())
        return name;
    return string() + this->type + string("_") + this->event_class + string("_") + seq_type
           + string("_") + name_side(this->type, this->event_side)
           + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

Rec_Event_name Rec_Event::get_legacy_name() const
{
    return string() + this->type + string("_") + this->event_class + string("_")
           + name_side(this->type, this->event_side)
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
        // One invariant with two halves: by the time a scenario reaches a leaf, every seq_type
        // in the model must have both a sequence and its offsets -- not necessarily written by
        // the same event. Checked here because this is the one boundary where it has to hold,
        // and only under assertions -- each walk is linear in the scenario's length, and the
        // default build defines NDEBUG, so release pays nothing.
        //
        // Content. int_undefined is a placeholder inside one scenario, never a value a consumer
        // can interpret. Since O12 (a') no event hands one on at all -- Dinucl_markov creates
        // the junction it fills -- so this has gone from a leaf property to a global one, and
        // the leaf is now the last place it could possibly fail rather than the only place it
        // was ever true.
        if (const SeqTypeId unfilled = first_unfilled_segment(scenario.constructed_sequences);
            unfilled != kNoSeqType) {
            std::cerr << "Scenario leaf reached with unfilled nucleotides in "
                      << scenario.constructed_sequences.registry().name(unfilled)
                      << std::endl;
            assert(false && "unfilled nucleotide in a completed scenario");
        }

        // Placement, which became checkable with R3: until Insertion wrote its junction's
        // offsets this fired on every scenario that has one. A segment nobody placed is one no
        // consumer can locate on the read, and the only reason that did not bite is that every
        // consumer today knows an insertion's span runs between its neighbours -- exactly the
        // coupling the capability queries exist to remove.
        if (const auto [unplaced, side] = first_unplaced_segment_end(scenario.seq_offsets);
            unplaced != kNoSeqType) {
            std::cerr << "Scenario leaf reached with no "
                      << (side == Five_prime ? "5'" : "3'") << " offset for "
                      << scenario.seq_offsets.five_prime.registry().name(unplaced) << std::endl;
            assert(false && "unplaced segment end in a completed scenario");
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

namespace {

//What an event leaves in the accumulator for one realization: how many nucleotides of its
//segment its offsets imply that nobody has chosen yet (see the body below for the three cases).
//Shared by the body, for a participant, and by initialize_Len_proba_bound(), for the reader: the
//reader publishes exactly as any participant would -- it is only its probability that stays out.
void publish_unfilled_length(SeqTypeId seq_type_id, int delta, bool creates_offsets, bool creates_sequence,
                             UnfilledSegmentLengths &lengths)
{
    if (creates_offsets) {
        lengths.set(seq_type_id, creates_sequence ? 0 : delta);
    }
}

} // namespace

void Rec_Event::SpanFold::run(FoldState start) const
{
    //The events that neither change this span's length nor contribute a probability factor to it
    //were dropped when `participants` was built, so there is one step per contributing event
    //rather than per event in the model.
    FoldFrontier frontier;
    frontier.offer(std::move(start), 1.0);
    for (std::size_t cursor = 0; cursor != participants.size(); ++cursor) {
        FoldFrontier next;
        for (const auto &[state, proba] : frontier) {
            participants[cursor]->fold_step(*this, cursor, state, proba, next);
        }
        frontier = std::move(next);
    }

    //Every event contributing to this span has chosen. Several states can share a length --
    //they differ in what a participant that has now run would have read -- and
    //SpanProfile::record() keeps the best of them.
    for (const auto &[state, proba] : frontier) {
        profile.record(state.length, proba);
    }
}

/*
 * One step for all four event kinds. Called only for an event the filter has already found to
 * participate in the span, so the two branches below are "enumerates its realizations" and
 * "contributes a probability factor without enumerating" -- which is exactly the
 * affects_length_of / affects_proba_of split.
 */
void Rec_Event::fold_step(const SpanFold &fold, std::size_t cursor, const FoldState &state, double proba,
                          FoldFrontier &next) const
{
    if (not this->affects_length_of(fold.span)) {
        //Dinucl_markov: no realization of its own contributes length, and its p^L factor reads a
        //length some upstream creator already published. One factor, one state on.
        FoldState successor = state;
        successor.parent_offsets[cursor] = 0;
        next.offer(std::move(successor), proba * this->span_proba_factor(fold.span, state.lengths));
        return;
    }

    //What the accumulator carries is *how many nucleotides a segment's offsets imply that
    //nobody has chosen yet* -- how many are still to be chosen. So the publisher is whoever
    //places the offsets, and the value is the length they imply only while no sequence exists
    //for them yet (O12, plan section 6.10 finding 8).
    //
    //Three cases, and each key still has exactly one writer per path:
    //  - Insertion  -- creates the offsets, not the sequence: publishes the whole length,
    //                  which is what Dinucl_markov::span_proba_factor raises p to.
    //  - Gene_choice -- creates both: publishes 0. Its nucleotides are fixed by the template,
    //                  so none of them are still to be chosen. It used to publish the template
    //                  length, which nothing reads, so this corrects a silent wrong answer
    //                  rather than changing one anybody sees.
    //  - Deletion   -- modifies an offset rather than creating one: publishes nothing, and
    //                  contributes its negative delta to the span total instead.
    const bool creates_offsets = this->creates_own_offsets();
    const bool creates_sequence =
            this->get_seq_construction_role(this->seq_type_id) == SeqConstructionRole::Creates;

    //Read at layer 0, where every event's base index sits before the walk writes any other: a
    //pure read, where the fold used to rewind the key's current layer to 0 first -- a write
    //from a const traversal, and a no-op at initialization, when nothing stands above layer 0.
    //What the path has already chosen of this event's parents is added on top, exactly as the
    //walk's update_parent_tracking() would have added it (R6).
    const int base_index = fold.base_index_map.get(this->event_index, 0) + state.parent_offsets[cursor];
    const SpanConditioning::Participant &conditioning = fold.conditioning.participants[cursor];

    for (std::unordered_map<std::string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        const Event_realization &realization = iter->second;
        const int delta = this->length_delta(realization);
        const double bound = this->realization_bound(realization, fold.model_parameters, base_index,
                                                     conditioning.free_parent_offsets);

        FoldState successor = state;
        successor.length += delta;
        //Read above, and by nobody after: clearing it is what lets two paths that reached this
        //participant from different parents merge from here on.
        successor.parent_offsets[cursor] = 0;
        publish_unfilled_length(this->seq_type_id, delta, creates_offsets, creates_sequence, successor.lengths);
        for (const auto &[child, stride] : conditioning.children) {
            successor.parent_offsets[child] += realization.index * stride;
        }
        next.offer(std::move(successor), proba * bound);
    }
}

double Rec_Event::realization_bound(const Event_realization &realization, const Marginal_array_p &model_parameters,
                                    int base_index, const std::vector<int> &parent_offsets) const
{
    //The realization's marginal is stored once per realization of the conditioning parents, at a
    //stride of this event's size, and the bound has to hold whichever of those the scenario
    //turns out to carry -- among the ones still open.
    double best = 0;
    for (const int parent_offset : parent_offsets) {
        if (model_parameters[base_index + realization.index + parent_offset] > best) {
            best = model_parameters[base_index + realization.index + parent_offset];
        }
    }
    return best;
}

std::vector<int> Rec_Event::parent_offsets_free_of(const std::vector<std::pair<int, int>> &fixed) const
{
    //The parents' realizations are a mixed-radix number, each parent's digit weighted by the
    //stride it was given in memory_and_offsets, so a configuration leaves a parent at its first
    //realization exactly when that digit is 0.
    std::vector<int> offsets;
    offsets.reserve(this->event_marginal_size / this->size());
    for (std::size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
        const int offset = static_cast<int>(i * this->size());
        bool free = true;
        for (const auto &[stride, parent_size] : fixed) {
            if ((offset / stride) % parent_size != 0) {
                free = false;
                break;
            }
        }
        if (free) {
            offsets.push_back(offset);
        }
    }
    return offsets;
}

Rec_Event::SpanConditioning Rec_Event::conditioning_within(const SpanParticipants &participants, SegmentSpan span,
                                                           bool reader_conditioned) const
{
    //Who indexes whom, read off the parents' side: an event's memory_and_offsets lists its
    //children and the stride its realization moves each of them by. Only an event the fold
    //enumerates has a realization on the path -- one that changes the span's length, or the
    //reader of a table conditioned on its own -- so only those fix anything.
    SpanConditioning conditioning;
    conditioning.participants.resize(participants.size());
    std::vector<std::vector<std::pair<int, int>>> fixed(participants.size());

    const auto link_children = [&](const Rec_Event &parent, std::size_t first_cursor,
                                   std::vector<std::pair<std::size_t, int>> &children) {
        for (const auto &[child_id, layer, stride] : parent.memory_and_offsets) {
            for (std::size_t cursor = first_cursor; cursor != participants.size(); ++cursor) {
                if (participants[cursor]->event_index == child_id) {
                    children.emplace_back(cursor, stride);
                    fixed[cursor].emplace_back(stride, parent.size());
                }
            }
        }
    };

    if (reader_conditioned) {
        link_children(*this, 0, conditioning.reader_children);
    }
    for (std::size_t cursor = 0; cursor != participants.size(); ++cursor) {
        if (participants[cursor]->affects_length_of(span)) {
            link_children(*participants[cursor], cursor + 1, conditioning.participants[cursor].children);
        }
    }
    for (std::size_t cursor = 0; cursor != participants.size(); ++cursor) {
        conditioning.participants[cursor].free_parent_offsets = participants[cursor]->parent_offsets_free_of(fixed[cursor]);
    }
    return conditioning;
}

bool Rec_Event::creates_own_offsets() const
{
    return this->get_offset_role(this->seq_type_id, Five_prime) == OffsetRole::Creates
           || this->get_offset_role(this->seq_type_id, Three_prime) == OffsetRole::Creates;
}

/*
 * The four per-subclass overrides, minus their enum-to-member switches. Everything topological
 * was settled in initialize_event(), so all that is left is "fold each junction I read".
 */
void Rec_Event::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                           const Marginal_array_p &model_parameters_point,
                                           const Index_map &base_index_map)
{
    //Every fold starts from a state of its own, with no length published: an entry is a length
    //published by the segment's creator on one path, and the paths of two folds share nothing.
    const std::size_t seq_type_count = legacy_seq_type_registry().total_count();

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

        //A table never contains its reader: what the reader reads is a bound on what is still to
        //be realized once it has chosen, and its own choice is realized by then (§7.19, R8). Its
        //probability therefore stays out -- it is already in the scenario the bound multiplies.
        //
        //Its realization must not be maxed over either, and that is the half R8 first got wrong
        //by role: if the reader changes this span's length, the table is **conditioned** on its
        //realization -- one profile each, folded from that realization at weight 1, with its
        //length counted in the key and published into the accumulator exactly as a
        //participant's would be. So the key is always the summed `length_delta` of everything
        //enumerated, and a reader looks up its own realization's profile at the gap as it
        //stood before its choice. Nothing here depends on what kind of event the reader is or
        //on what follows it (R13).
        //
        //A reader that does not change the span's length -- the gene choice the span is
        //measured from -- has nothing to condition on: one profile, from the events after it.
        //
        //Within the table, a participant conditioned on an event the fold has already chosen on
        //the current path -- an earlier participant, or the reader when the table is conditioned
        //on it -- reads its marginals at that realization rather than at the best of them (R6).
        const SegmentSpan span = bound.span();
        const bool conditioned = this->affects_length_of(span);
        bound.reset_profiles(conditioned, this->size());
        const SpanConditioning conditioning = this->conditioning_within(participants, span, conditioned);

        if (not conditioned) {
            SpanFold{span, participants, conditioning, model_parameters_point, base_index_map,
                     bound.mutable_profile()}
                    .run(FoldState{0, std::vector<int>(participants.size(), 0),
                                   UnfilledSegmentLengths(seq_type_count)});
            continue;
        }

        const bool creates_offsets = this->creates_own_offsets();
        const bool creates_sequence =
                this->get_seq_construction_role(this->seq_type_id) == SeqConstructionRole::Creates;
        for (std::unordered_map<std::string, Event_realization>::const_iterator iter = this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {
            const Event_realization &realization = iter->second;
            const int delta = this->length_delta(realization);
            FoldState start{delta, std::vector<int>(participants.size(), 0), UnfilledSegmentLengths(seq_type_count)};
            publish_unfilled_length(this->seq_type_id, delta, creates_offsets, creates_sequence, start.lengths);
            for (const auto &[child, stride] : conditioning.reader_children) {
                start.parent_offsets[child] += realization.index * stride;
            }
            SpanFold{span, participants, conditioning, model_parameters_point, base_index_map,
                     bound.mutable_profile_for(realization.index)}
                    .run(std::move(start));
        }
    }

    this->build_retained_decomposition(model_parameters_point, base_index_map);
}

/*
 * The `Retain` mode of the fold: not "how good can this span get at each length", but "which
 * ways can this event sit inside it, at each length, best first". Section 2.5's composition
 * operator in the variant that keeps its arguments -- which is why it is built here, next to the
 * halves it composes, and not derived from them afterwards.
 *
 * `Retain` sits on the enclosing junction because that is the one this event splits: its two
 * halves are the left and right slots, already folded above, and the third factor is this
 * event's own realization set. Nothing here is a gene: the template's length is
 * `length_delta()`, the per-realization bound is the same maxᵢ the fold takes, and which
 * junctions are involved was settled in initialize_event().
 */
void Rec_Event::build_retained_decomposition(const Marginal_array_p &model_parameters_point,
                                             const Index_map &base_index_map)
{
    JunctionBound &enclosing = junction_bounds_[kEnclosingJunction];
    if (not enclosing.resolved() or not enclosing.retained()) {
        return;
    }

    SpanDecomposition &decomposition = enclosing.mutable_decomposition();
    decomposition.clear();

    const JunctionBound &left = junction_bounds_[kLeftJunction];
    const JunctionBound &right = junction_bounds_[kRightJunction];
    if (not left.resolved() or not right.resolved()) {
        //Only one flank placed, so there is no enclosed span to decompose. The consumer falls
        //back to scanning the read, which needs no decomposition.
        return;
    }

    const int base_index = base_index_map.get(this->event_index, 0);
    //Every parent of the reader is upstream of it, so none is known when the table is built.
    const std::vector<int> every_parent_offset = this->parent_offsets_free_of({});

    for (std::unordered_map<std::string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        const Event_realization &realization = iter->second;

        //The same bound the fold takes, for the same reason: this event may be a child, and the
        //bound has to hold whichever parent realization the scenario turns out to carry.
        const double real_max_proba =
                this->realization_bound(realization, model_parameters_point, base_index, every_parent_offset);

        const int own_length = this->length_delta(realization);

        for (const SpanProfile::Entry near_ : left.profile()) {
            for (const SpanProfile::Entry far_ : right.profile()) {
                decomposition.record(own_length + near_.distance + far_.distance,
                                     SpanDecomposition::Placement{realization.index, near_.distance,
                                                                  far_.distance,
                                                                  real_max_proba * near_.proba * far_.proba});
            }
        }
    }

    //The consumer `break`s on the first prune, which is exact only in this order.
    decomposition.sort_by_decreasing_proba();
}

void Rec_Event::adopt_Len_proba_bound(const Rec_Event &source)
{
    for (std::size_t slot = 0; slot != kJunctionSlotCount; ++slot) {
        JunctionBound &bound = junction_bounds_[slot];
        if (not bound.resolved() or not bound.folded()) {
            continue;
        }
        bound.adopt_profiles(source.junction_bounds_[slot]);
    }

    JunctionBound &enclosing = junction_bounds_[kEnclosingJunction];
    if (enclosing.resolved() and enclosing.retained()) {
        enclosing.mutable_decomposition() = source.junction_bounds_[kEnclosingJunction].decomposition();
    }
}

} // namespace igor::model::legacy
