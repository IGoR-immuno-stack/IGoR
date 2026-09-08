/*
 * test_utils.cpp
 *
 *  Created on: Jan 21, 2026
 *      Author: IGoR Test Suite
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

#include "test_utils.h"

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/EventUtils.h>
#include <igor/Core/gene_to_seqtype_migr.h>
#include <algorithm>
#include <cmath>
#include <forward_list>
#include <functional>
#include <numeric>
#include <queue>
#include <stack>
#include <stdexcept>

namespace IgorTestUtils {

Alignment_data create_mock_alignment_data(
    const std::string& gene_name,
    int offset,
    size_t five_p_offset,
    size_t three_p_offset,
    const std::vector<size_t>& mismatches,
    double score
) {
    // Calculate alignment length
    size_t align_length = three_p_offset - five_p_offset;
    
    // Create empty vectors for insertions and deletions
    std::vector<size_t> empty_insertions;
    std::vector<size_t> empty_deletions;
    
    // Use the appropriate constructor
    Alignment_data align_data(
        gene_name,
        offset,
        five_p_offset,
        three_p_offset,
        align_length,
        empty_insertions,
        empty_deletions,
        mismatches,
        score
    );
    
    return align_data;
}

// ============================================================================
// iterate() test harness
// ============================================================================

const SeqTypeRegistry &vdj_seq_type_registry()
{
    static const SeqTypeRegistry registry = [] {
        SeqTypeRegistry built;
        built.register_legacy_seq_types();
        built.set_ordered_types({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
        built.freeze();
        return built;
    }();
    return registry;
}

const SeqTypeRegistry &vj_seq_type_registry()
{
    static const SeqTypeRegistry registry = [] {
        SeqTypeRegistry built;
        built.register_legacy_seq_types();
        built.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
        built.freeze();
        return built;
    }();
    return registry;
}

IterateTestState create_iterate_state(const std::string &sequence, std::size_t marginal_array_size,
                                      std::size_t max_events, const SeqTypeRegistry &registry)
{
    return IterateTestState(sequence, marginal_array_size, max_events, registry);
}

std::string LayerViolation::describe() const
{
    return map_name + " key " + std::to_string(key) + ": claimed layer "
           + std::to_string(claimed_layer) + " but its data stands at layer "
           + std::to_string(current_layer)
           + " -- a layer was requested and never written on this path";
}

namespace {

template <typename Map>
std::vector<int> claimed_layers_of(const Map &map)
{
    std::vector<int> layers;
    layers.reserve(map.count());
    for (std::size_t key = 0; key != map.count(); ++key) {
        layers.push_back(map.claimed_layer(key));
    }
    return layers;
}

template <typename Map>
std::vector<int> current_layers_of(const Map &map)
{
    std::vector<int> layers;
    layers.reserve(map.count());
    for (std::size_t key = 0; key != map.count(); ++key) {
        layers.push_back(map.current_layer(key));
    }
    return layers;
}

} // namespace

LayerSnapshot capture_layers(const IterateTestState &state)
{
    LayerSnapshot snapshot;
    auto add = [&](const std::string &name, auto &&map) {
        snapshot.claimed.emplace(name, claimed_layers_of(map));
        snapshot.current.emplace(name, current_layers_of(map));
    };
    add("constructed_sequences", state.scenario.constructed_sequences);
    add("seq_offsets.five_prime", state.scenario.seq_offsets.five_prime);
    add("seq_offsets.three_prime", state.scenario.seq_offsets.three_prime);
    add("mismatches_lists", state.scenario.mismatches_lists);
    add("downstream_proba_map", state.exploration.downstream_proba_map);
    add("safety_set", state.exploration.safety_set);
    add("pruning_mismatch_floor", state.exploration.pruning_mismatch_floor);
    //index_map is deliberately excluded: its layering is driven by parent-realization
    //tracking through offset_map, which single-event tests do not populate, so it carries no
    //contract here.
    return snapshot;
}

std::string int_str_to_nt(const Int_Str &seq)
{
    static const char kBases[] = {'A', 'C', 'G', 'T'};
    std::string out;
    out.reserve(seq.size());
    for (int value : seq) {
        if (value >= 0 && value < 4) {
            out.push_back(kBases[value]);
        } else if (value == int_undefined) {
            //Rendered apart from 'N' on purpose: "not filled yet" and "filled, but the read is
            //ambiguous here" are different states, and a failure message that prints both the
            //same way hides exactly the confusion this notation exists to prevent.
            out.push_back('.');
        } else {
            out.push_back('N');
        }
    }
    return out;
}

RecordingEvent::RecordingEvent(int event_id) : Rec_Event()
{
    this->type = Event_type::Undefined_t;
    this->set_event_identifier(event_id);
    this->set_seq_type("Undefined_seq");
    this->fix(true);
}

std::shared_ptr<Rec_Event> RecordingEvent::copy()
{
    return std::make_shared<RecordingEvent>(this->get_event_identifier());
}

std::queue<int> RecordingEvent::draw_random_realization(
        const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
        const std::unordered_map<Rec_Event_name,
                                 std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
        std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const
{
    return std::queue<int>();
}

void RecordingEvent::iterate(QuerySequenceContext &, const ModelContext &, ScenarioContext &scenario,
                             ExplorationContext &exploration, AccumulationContext &)
{
    //Record what the real next event would read, then stop. No recursion, no error rate.
    static const Seq_type kAllSeqTypes[] = {V_gene_seq, VD_ins_seq,  D_gene_seq,
                                            DJ_ins_seq, J_gene_seq,  VJ_ins_seq};

    ScenarioSnapshot snapshot;
    snapshot.scenario_proba = scenario.scenario_proba;

    for (Seq_type seq_type : kAllSeqTypes) {
        if (scenario.seq_offsets.exists(seq_type, Five_prime)
            && scenario.seq_offsets.exists(seq_type, Three_prime)) {
            snapshot.offsets.emplace(seq_type,
                                     std::make_pair(scenario.get_offset(seq_type, Five_prime),
                                                    scenario.get_offset(seq_type, Three_prime)));
        }
        if (scenario.constructed_sequences.exists(seq_type)) {
            const Int_Str *segment = scenario.get_sequence_segment(seq_type);
            snapshot.sequences.emplace(seq_type,
                                       segment ? int_str_to_nt(*segment) : std::string());
        }
        if (scenario.mismatches_lists.exists(seq_type)) {
            const std::vector<std::size_t> *mismatches = scenario.get_mismatches(seq_type);
            snapshot.mismatches.emplace(seq_type, mismatches ? *mismatches
                                                             : std::vector<std::size_t>{});
        }
        if (exploration.downstream_proba_map.exists(seq_type)) {
            snapshot.downstream_bounds.emplace(seq_type,
                                               exploration.downstream_proba_map.get(seq_type));
        }
    }

    for (Event_safety safety : {Event_safety::VD_safe, Event_safety::VJ_safe, Event_safety::DJ_safe}) {
        if (exploration.safety_set.exists(safety)) {
            snapshot.safety.emplace(safety, exploration.safety_set.get(safety));
        }
    }

    //Layer contract: every layer this event requested must have been written before it
    //hands off. Checked here rather than in each test, so all sections get it for free.
    auto check_map = [&](const std::string &name, const std::vector<int> &current_now) {
        const auto claimed = layer_baseline.claimed.find(name);
        const auto before = layer_before_init.claimed.find(name);
        if (claimed == layer_baseline.claimed.end() || before == layer_before_init.claimed.end()) {
            return;
        }
        for (std::size_t key = 0; key != current_now.size(); ++key) {
            if (key >= claimed->second.size() || key >= before->second.size()) {
                break;
            }
            //Only keys this event claimed a layer for carry the promise.
            if (claimed->second[key] <= before->second[key]) {
                continue;
            }
            if (current_now[key] != claimed->second[key]) {
                layer_violations.push_back(
                        LayerViolation{calls.size(), name, key, claimed->second[key],
                                       current_now[key]});
            }
        }
    };
    check_map("constructed_sequences", current_layers_of(scenario.constructed_sequences));
    check_map("seq_offsets.five_prime", current_layers_of(scenario.seq_offsets.five_prime));
    check_map("seq_offsets.three_prime", current_layers_of(scenario.seq_offsets.three_prime));
    check_map("mismatches_lists", current_layers_of(scenario.mismatches_lists));
    check_map("downstream_proba_map", current_layers_of(exploration.downstream_proba_map));
    check_map("safety_set", current_layers_of(exploration.safety_set));
    check_map("pruning_mismatch_floor", current_layers_of(exploration.pruning_mismatch_floor));

    calls.push_back(std::move(snapshot));
}

std::shared_ptr<RecordingEvent> call_iterate_recording(const std::shared_ptr<Rec_Event> &event,
                                                       IterateTestState &state)
{
    auto recorder = std::make_shared<RecordingEvent>(31);
    call_iterate(event, state, recorder);

    //Every section gets the layer contract checked, without asking for it.
    for (const LayerViolation &violation : recorder->layer_violations) {
        UNSCOPED_INFO("layer contract violated at hand-off " << violation.call_index << ": "
                                                             << violation.describe());
    }
    CHECK(recorder->layer_violations.empty());

    return recorder;
}

void call_iterate(const std::shared_ptr<Rec_Event> &event, IterateTestState &state,
                  const std::shared_ptr<RecordingEvent> &next)
{
    //The order below mirrors GenModel's setup. It is not incidental: the length-proba
    //bounds are built by a reverse pass over the queue *after* every event is initialized,
    //and omitting it leaves the *_length_best_proba_map members empty, which makes the
    //junction-length guard discard every scenario. That is the single easiest way to write
    //a test that passes for the wrong reason.

    if (next) {
        next->layer_before_init = capture_layers(state);
    }

    auto &events_map = const_cast<Events_map &>(state.model.events_map);
    auto &offset_map = const_cast<std::unordered_map<
            Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &>(
            state.model.offset_map);

    //The event under test must be reachable through events_map like any other.
    events_map[std::make_tuple(event->get_type(), event->get_seq_type(), event->get_side())] = event;

    //Step 1: every event gets a base index of 0 at layer 0, plus a marginal size and a
    //crude upper bound. iterate_common() and add_to_marginals() both read these.
    for (const auto &[key, ev] : events_map) {
        (void)key;
        const int event_index = ev->get_event_identifier();
        state.exploration.index_map.request_layer(event_index);
        state.exploration.index_map.set(event_index, state.base_index_for(event_index), 0);
        ev->set_event_marginal_size(ev->size());
        ev->set_crude_upper_bound_proba(0, ev->size(),
                                        const_cast<Marginal_array_p &>(state.model.model_parameters));
        ev->set_viterbi_run(false);
    }

    //Step 2: the queue holds the event under test alone. Neighbours marked with
    //mark_chosen() are deliberately *not* initialized: initialize_event() requests a memory
    //layer per seq_type it touches, so running it on a neighbour would shift that
    //neighbour's current layer to 1 while preset_segment() writes at 0, and the event under
    //test would read an unwritten layer. Registering them in events_map and in
    //processed_events reproduces everything the code under test actually queries about them
    //-- existence, chosen-ness, and their offsets -- with the layer numbering pinned.
    std::queue<std::shared_ptr<Rec_Event>> queue;
    std::stack<std::shared_ptr<Rec_Event>> init_stack;
    queue.push(event);
    for (const auto &downstream : state.downstream_events()) {
        queue.push(downstream);
    }
    state.model_queue() = queue;

    //Step 3: initialize_event() in queue order. Note the aliasing: processed_events is the
    //state's own set, so the chosen events stay marked while the event under test
    //initializes, which is exactly what drives its *_chosen flags.
    {
        std::queue<std::shared_ptr<Rec_Event>> init_queue = queue;
        while (!init_queue.empty()) {
            std::shared_ptr<Rec_Event> ev = init_queue.front();
            init_queue.pop();
            init_stack.push(ev);
            ev->initialize_event(state.processed_events(), events_map, offset_map,
                                 state.exploration.downstream_proba_map,
                                 state.scenario.constructed_sequences, state.exploration.safety_set,
                                 state.accumulation.error_rate, state.scenario.mismatches_lists,
                                 state.scenario.seq_offsets, state.exploration.index_map);

            //The layer contract baseline, captured the instant the event under test has
            //requested its layers and before any downstream event requests more. Downstream
            //events are initialized but never iterate -- the recorder intercepts first -- so
            //their requested layers are legitimately unwritten and must not be in the
            //baseline. See LayerContract in test_utils.h.
            if (next && ev == event) {
                next->layer_baseline = capture_layers(state);
            }
        }
    }

    //Step 4: next-event chain. With a recorder, the event under test hands off to it and
    //stops; without one, iterate_wrap_up() takes the leaf path, which pulls in Error_rate
    //and therefore needs a *complete* scenario (V, D and J all constructed).
    for (const auto &[key, ev] : events_map) {
        (void)key;
        state.exploration.next_event_ptr_arr.get()[ev->get_event_identifier()] = nullptr;
    }
    if (next) {
        state.exploration.next_event_ptr_arr.get()[event->get_event_identifier()] = next.get();
    }

    //Step 5: probability bounds, reverse queue order.
    {
        double downstream_proba_bound = 1.0;
        std::forward_list<double *> updated_proba_list;
        while (!init_stack.empty()) {
            std::shared_ptr<Rec_Event> ev = init_stack.top();
            init_stack.pop();

            std::queue<std::shared_ptr<Rec_Event>> remaining = queue;
            while (!remaining.empty() && remaining.front() != ev) {
                remaining.pop();
            }
            if (!remaining.empty()) {
                remaining.pop();
            }

            ev->initialize_crude_scenario_proba_bound(downstream_proba_bound, updated_proba_list,
                                                      events_map);
            ev->initialize_Len_proba_bound(remaining,
                                           const_cast<Marginal_array_p &>(state.model.model_parameters),
                                           state.exploration.index_map);
        }
    }

    event->iterate(state.query, state.model, state.scenario, state.exploration, state.accumulation);
}

// ============================================================================
// State inspection
// ============================================================================

Seq_Offset get_seq_offset(const IterateTestState &state, Seq_type seq_type, Seq_side side,
                          std::size_t layer)
{
    return state.scenario.seq_offsets.get(seq_type, side, layer);
}

bool has_seq_offset(const IterateTestState &state, Seq_type seq_type, Seq_side side)
{
    return state.scenario.seq_offsets.exists(seq_type, side);
}

const Int_Str *get_constructed_sequence(const IterateTestState &state, Seq_type seq_type,
                                        std::size_t layer)
{
    return state.scenario.constructed_sequences.get(seq_type, layer);
}

bool has_constructed_sequence(const IterateTestState &state, Seq_type seq_type)
{
    return state.scenario.constructed_sequences.exists(seq_type);
}

std::vector<std::size_t> get_mismatches(const IterateTestState &state, Seq_type seq_type,
                                        std::size_t layer)
{
    const std::vector<std::size_t> *v = state.scenario.mismatches_lists.get(seq_type, layer);
    return v ? *v : std::vector<std::size_t>{};
}

bool is_safe(const IterateTestState &state, Event_safety safety_type, std::size_t layer)
{
    return state.exploration.safety_set.get(safety_type, layer);
}

bool has_safety(const IterateTestState &state, Event_safety safety_type)
{
    return state.exploration.safety_set.exists(safety_type);
}

int safety_current_layer(const IterateTestState &state, Event_safety safety_type)
{
    return static_cast<int>(state.exploration.safety_set.claimed_layer(safety_type));
}

double get_downstream_bound(const IterateTestState &state, Seq_type seq_type, std::size_t layer)
{
    return state.exploration.downstream_proba_map.get(seq_type, layer);
}

long double total_marginal_mass(const IterateTestState &state, std::size_t marginal_array_size)
{
    long double total = 0.0;
    for (std::size_t i = 0; i != marginal_array_size; ++i) {
        total += state.accumulation.updated_marginals[i];
    }
    return total;
}

// ============================================================================
// Event builders
// ============================================================================

std::shared_ptr<Gene_choice> make_gene_choice(Gene_class gene_class,
                                              const std::vector<std::pair<std::string, std::string>> &genes,
                                              int event_id, bool fixed)
{
    auto event = std::make_shared<Gene_choice>(gene_class);
    for (const auto &[name, sequence] : genes) {
        event->add_realization(name, sequence);
    }
    event->set_event_identifier(event_id);
    event->set_priority(1);
    Seq_type target = V_gene_seq;
    if (!igor::migration::try_gene_class_to_gene_seq_type(gene_class, target)) {
        throw std::invalid_argument("make_gene_choice: gene class has no gene seq_type");
    }
    const Seq_type_String seq_type = EventUtils::seq_type_to_string(target);
    event->set_seq_type(seq_type);
    event->set_seq_type_id(legacy_seq_type_registry().id(seq_type));
    event->update_event_name();
    event->fix(fixed);
    return event;
}

std::shared_ptr<Deletion> make_deletion(Seq_type target, Seq_side side, int min_del, int max_del,
                                        int event_id)
{
    auto event = std::make_shared<Deletion>(target, side, std::make_pair(min_del, max_del));
    event->set_event_identifier(event_id);
    event->set_priority(1);
    const Seq_type_String seq_type = EventUtils::seq_type_to_string(target);
    event->set_seq_type(seq_type);
    event->set_seq_type_id(legacy_seq_type_registry().id(seq_type));
    event->update_event_name();
    event->fix(true);
    return event;
}

std::shared_ptr<Insertion> make_insertion(Seq_type target, int min_ins, int max_ins, int event_id)
{
    auto event = std::make_shared<Insertion>(target, std::make_pair(min_ins, max_ins));
    event->set_event_identifier(event_id);
    event->set_priority(1);
    const Seq_type_String seq_type = EventUtils::seq_type_to_string(target);
    event->set_seq_type(seq_type);
    event->set_seq_type_id(legacy_seq_type_registry().id(seq_type));
    event->update_event_name();
    event->fix(true);
    return event;
}

std::string segment_run(Seq_Offset five_prime, Seq_Offset three_prime)
{
    if (three_prime < five_prime - 1) {
        throw std::invalid_argument("segment_run: three_prime is more than one before five_prime; "
                                    "the empty-segment convention is off(3') == off(5') - 1");
    }
    static const std::string kPattern = "ACGT";
    const auto length = static_cast<std::size_t>(three_prime - five_prime + 1);
    std::string run;
    run.reserve(length);
    for (std::size_t i = 0; i != length; ++i) {
        run += kPattern[i % kPattern.size()];
    }
    return run;
}

std::shared_ptr<Dinucl_markov> make_dinucl_markov(Seq_type target, int event_id)
{
    auto event = std::make_shared<Dinucl_markov>(target);
    event->set_event_identifier(event_id);
    event->set_priority(1);
    const Seq_type_String seq_type = EventUtils::seq_type_to_string(target);
    event->set_seq_type(seq_type);
    event->set_seq_type_id(legacy_seq_type_registry().id(seq_type));
    event->update_event_name();
    event->fix(true);
    return event;
}

Alignment_data create_perfect_alignment(const std::string &gene_name, int offset, int gene_length)
{
    return create_mock_alignment_data(gene_name, offset, offset >= 0 ? offset : 0,
                                      offset + gene_length - 1, {}, 100.0);
}

Alignment_data create_alignment_with_mismatches(const std::string &gene_name, int offset,
                                                int gene_length,
                                                const std::vector<std::size_t> &mismatch_positions)
{
    return create_mock_alignment_data(gene_name, offset, offset >= 0 ? offset : 0,
                                      offset + gene_length - 1, mismatch_positions,
                                      100.0 - 5.0 * static_cast<double>(mismatch_positions.size()));
}

} // namespace IgorTestUtils
