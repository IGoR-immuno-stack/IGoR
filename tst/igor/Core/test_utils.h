/*
 * test_utils.h
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

#pragma once

#include <igor/Core/AccumulationContext.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/ExplorationContext.h>
#include <igor/Core/GenModel.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/ModelContext.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/QuerySequenceContext.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/ScenarioContext.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Singleerrorrate.h>
#include <igor/Core/Utils.h>

#include <deque>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace IgorTestUtils {

/**
 * @brief Create mock alignment data for testing
 * 
 * @param gene_name Name of the gene (e.g., "TRBV1", "TRBJ1-1")
 * @param offset Alignment position on target sequence
 * @param five_p_offset 5' alignment boundary
 * @param three_p_offset 3' alignment boundary
 * @param mismatches Vector of mismatch positions
 * @param score Alignment score
 * @return Alignment_data Mock alignment data structure
 */
Alignment_data create_mock_alignment_data(
    const std::string& gene_name,
    int offset,
    size_t five_p_offset,
    size_t three_p_offset,
    const std::vector<size_t>& mismatches = {},
    double score = 100.0
);

/**
 * @brief iterate() test harness
 *
 * Ported and adapted from the sketch on feature/2_unittests; see
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md section 6.1. The Storage/Context split is that
 * branch's design. Everything below the Storage structs was rewritten for the
 * registry-sized containers of B2/B8, and for the ability to place a neighbouring event in
 * the "already chosen" state -- which the sketch could not express, and which is why every
 * safety-check section there was left an empty stub.
 *
 * Design contract:
 *  1. Storage objects own the data and guarantee the lifetime of the Context references.
 *  2. Context objects are the only surface tests touch.
 *  3. When a Context grows a field, update the corresponding Storage struct here.
 */

/// Layers preallocated in every layered container. request_layer() grows on demand, so
/// this is a starting point rather than a cap.
inline constexpr std::size_t kTestLayers = 8;

/// Storage for QuerySequenceContext.
struct QueryStorage {
    std::string sequence;
    Int_Str int_sequence;
    std::unordered_map<Gene_class, std::vector<Alignment_data>> gene_alignments;

    explicit QueryStorage(const std::string &seq)
        : sequence(seq), int_sequence(nt2int(seq)), gene_alignments{}
    {
        //Gene_choice indexes gene_alignments.at(gene_class) unconditionally, so every class
        //needs an entry even when a test supplies no alignment for it.
        gene_alignments[V_gene];
        gene_alignments[D_gene];
        gene_alignments[J_gene];
    }
};

/// Storage for ModelContext.
struct ModelStorage {
    Marginal_array_p model_marginals;
    std::unordered_map<Rec_Event_name,
                       std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map;
    Events_map events_map;
    std::queue<std::shared_ptr<Rec_Event>> model_queue;

    explicit ModelStorage(std::size_t marginal_array_size)
        : model_marginals(new long double[marginal_array_size]()),
          offset_map{}, events_map{}, model_queue{}
    {}
};

/// Storage for ScenarioContext. All three maps are sized from the frozen legacy registry,
/// so SeqTypeId == Seq_type for the six standard types and enum-keyed access is valid.
struct ScenarioStorage {
    double scenario_proba;
    Seq_type_str_p_map constructed_sequences;
    Seq_offsets_map seq_offsets;
    Mismatch_vectors_map mismatches_lists;

    ScenarioStorage()
        : scenario_proba(1.0),
          constructed_sequences(legacy_seq_type_registry(), kTestLayers),
          seq_offsets(legacy_seq_type_registry(), kTestLayers),
          mismatches_lists(legacy_seq_type_registry(), kTestLayers)
    {}
};

/// Storage for ExplorationContext.
struct ExplorationStorage {
    Downstream_scenario_proba_bound_map downstream_proba_map;
    double seq_max_prob;
    double proba_threshold;
    Index_map index_map;
    std::shared_ptr<Next_event_ptr> next_event_ptr_arr;
    Safety_bool_map safety_set;
    Pruning_mismatch_floor_map pruning_mismatch_floor;

    explicit ExplorationStorage(std::size_t max_events)
        : downstream_proba_map(legacy_seq_type_registry(), kTestLayers),
          //ExplorationContext copies proba_threshold_factor by value, so it can only be set
          //here. Fixing it at 1 makes seq_max_prob itself the pruning threshold, and it is a
          //reference, so set_pruning_threshold() can move it afterwards. Starting at 0 means
          //should_prune() is false for any non-negative bound: nothing is dropped
          //incidentally, only by an explicit geometric or range check.
          seq_max_prob(0.0),
          proba_threshold(1.0),
          index_map(max_events, kTestLayers),
          next_event_ptr_arr(new Next_event_ptr[max_events](),
                             std::default_delete<Next_event_ptr[]>()),
          safety_set(3, kTestLayers),
          pruning_mismatch_floor(legacy_seq_type_registry(), kTestLayers)
    {
        //Mirrors GenModel: downstream bounds start at 1 so multiply_all() is neutral.
        downstream_proba_map.init_first_layer(1.0);
    }
};

/// Storage for AccumulationContext.
struct AccumulationStorage {
    Marginal_array_p updated_marginals;
    std::map<std::size_t, std::shared_ptr<Counter>> counters_list;
    std::shared_ptr<Error_rate> error_rate;

    explicit AccumulationStorage(std::size_t marginal_array_size)
        : updated_marginals(new long double[marginal_array_size]()),
          counters_list{},
          error_rate(std::make_shared<Single_error_rate>(0.0))
    {}
};

/**
 * @brief Owns every input iterate() needs, and exposes the five contexts over it.
 *
 * Neither copyable nor movable: the contexts hold references into this object's own
 * storage, so a move would leave them dangling. C++17 guaranteed copy elision is what lets
 * create_iterate_state() return one by value regardless.
 */
struct IterateTestState {
private:
    QueryStorage query_storage;
    ModelStorage model_storage;
    ScenarioStorage scenario_storage;
    ExplorationStorage exploration_storage;
    AccumulationStorage accumulation_storage;

    /// Passed to initialize_event(). An event listed here reads to the code under test as
    /// "already chosen"/"already processed", which is what drives every *_chosen flag.
    std::unordered_set<Rec_Event_name> processed_events_;

    /// Stable backing store for sequences handed to the scenario by preset_segment().
    /// std::deque so that previously returned addresses stay valid as it grows.
    std::deque<Int_Str> preset_sequences_;
    std::deque<std::vector<std::size_t>> preset_mismatches_;

    /// Events after the event under test in the model queue, in order.
    std::vector<std::shared_ptr<Rec_Event>> downstream_;

    std::unordered_map<int, std::size_t> base_index_overrides_;

public:
    QuerySequenceContext query;
    ModelContext model;
    ScenarioContext scenario;
    ExplorationContext exploration;
    AccumulationContext accumulation;

    IterateTestState(const std::string &seq, std::size_t marginal_array_size, std::size_t max_events)
        : query_storage(seq),
          model_storage(marginal_array_size),
          scenario_storage(),
          exploration_storage(max_events),
          accumulation_storage(marginal_array_size),
          processed_events_{},
          preset_sequences_{},
          preset_mismatches_{},
          downstream_{},
          base_index_overrides_{},
          query(query_storage.sequence, query_storage.int_sequence, query_storage.gene_alignments),
          model(model_storage.model_marginals, model_storage.offset_map, model_storage.events_map,
                model_storage.model_queue),
          scenario(scenario_storage.scenario_proba, scenario_storage.constructed_sequences,
                   scenario_storage.seq_offsets, scenario_storage.mismatches_lists),
          exploration(exploration_storage.downstream_proba_map, exploration_storage.seq_max_prob,
                      exploration_storage.proba_threshold, exploration_storage.index_map,
                      exploration_storage.next_event_ptr_arr, exploration_storage.safety_set,
                      exploration_storage.pruning_mismatch_floor),
          accumulation(accumulation_storage.updated_marginals, accumulation_storage.counters_list,
                       accumulation_storage.error_rate)
    {}

    IterateTestState(const IterateTestState &) = delete;
    IterateTestState &operator=(const IterateTestState &) = delete;
    IterateTestState(IterateTestState &&) = delete;
    IterateTestState &operator=(IterateTestState &&) = delete;

    /// Alignments the gene choice under test will enumerate.
    void set_alignments(Gene_class gene_class, std::vector<Alignment_data> alignments)
    {
        query_storage.gene_alignments[gene_class] = std::move(alignments);
    }

    /**
     * Register an event so the code under test can find it. The key is derived from the
     * event itself, so the event must already carry its seq_type -- see
     * make_gene_choice()/make_deletion(), which set it.
     */
    void add_event(const std::shared_ptr<Rec_Event> &event)
    {
        model_storage.events_map[std::make_tuple(event->get_type(), event->get_seq_type(),
                                                 event->get_side())] = event;
    }

    /// Mark an event as already processed, so *_chosen is true for the event under test.
    void mark_chosen(const std::shared_ptr<Rec_Event> &event)
    {
        processed_events_.insert(event->get_name());
    }

    /**
     * Write a segment as an upstream event would have left it: offsets on both ends, a
     * constructed sequence and an (empty by default) mismatch list, all at layer 0.
     *
     * Pair this with mark_chosen() on the corresponding Gene_choice; mark_chosen() alone
     * makes the code *look* for the offsets, and they have to be there.
     */
    void preset_segment(Seq_type seq_type, Seq_Offset five_prime, Seq_Offset three_prime,
                        const std::string &segment = "",
                        const std::vector<std::size_t> &mismatches = {})
    {
        scenario_storage.seq_offsets.set(seq_type, Five_prime, five_prime, 0);
        scenario_storage.seq_offsets.set(seq_type, Three_prime, three_prime, 0);
        preset_sequences_.push_back(nt2int(segment));
        scenario_storage.constructed_sequences.set(seq_type, &preset_sequences_.back(), 0);
        preset_mismatches_.push_back(mismatches);
        scenario_storage.mismatches_lists.set(seq_type, &preset_mismatches_.back(), 0);
    }

    /// The probability the event under test inherits from upstream. Defaults to 1.
    void set_scenario_proba(double proba) { scenario_storage.scenario_proba = proba; }

    /// Override the base index an event reads its marginals from. call_iterate() puts every
    /// event at 0 unless told otherwise; a non-zero value pins that the marginal read is
    /// base_index + realization_index rather than realization_index alone.
    void set_base_index(int event_id, std::size_t base_index)
    {
        base_index_overrides_[event_id] = base_index;
    }

    std::size_t base_index_for(int event_id) const
    {
        const auto found = base_index_overrides_.find(event_id);
        return found == base_index_overrides_.end() ? 0 : found->second;
    }

    /// One entry of the model marginal array, by flat index.
    void set_marginal(std::size_t index, long double value)
    {
        model_storage.model_marginals[index] = value;
    }

    /**
     * Register an event that sits *after* the event under test in the model queue.
     *
     * These are never iterated -- the recorder intercepts first -- but they are what the
     * reverse initialize_Len_proba_bound() pass walks to populate the junction-length
     * maps. Without any, the map collapses to {0: 1.0} and the junction-length guard
     * discards every scenario whose neighbours are not exactly adjacent, which silently
     * turns a safety-check test into a test of the junction guard.
     */
    void add_downstream_event(const std::shared_ptr<Rec_Event> &event)
    {
        add_event(event);
        downstream_.push_back(event);
    }

    const std::vector<std::shared_ptr<Rec_Event>> &downstream_events() const { return downstream_; }

    /// Replace the error model. Defaults to Single_error_rate(0.0); a non-zero rate is
    /// what makes the downstream error bound a number worth asserting on.
    void set_error_rate(double rate)
    {
        accumulation_storage.error_rate = std::make_shared<Single_error_rate>(rate);
    }

    /// Prune any scenario whose upper bound falls below `threshold`. Off by default.
    /// The threshold *factor* cannot be changed after construction (ExplorationContext
    /// holds it by value), so this moves seq_max_prob_scenario instead, which the context
    /// holds by reference; with the factor pinned at 1 the two are the same knob.
    void set_pruning_threshold(double threshold) { exploration_storage.seq_max_prob = threshold; }

    std::unordered_set<Rec_Event_name> &processed_events() { return processed_events_; }
    std::queue<std::shared_ptr<Rec_Event>> &model_queue() { return model_storage.model_queue; }
};

/// Factory. marginal_array_size and max_events are generous defaults; raise them for a
/// test that needs more.
IterateTestState create_iterate_state(const std::string &sequence,
                                      std::size_t marginal_array_size = 1000,
                                      std::size_t max_events = 32);

/// Decode an Int_Str back to letters. Placeholders (-1, written by Insertion) become 'N'.
std::string int_str_to_nt(const Int_Str &seq);

/**
 * @brief What the next event in the chain would see, captured per surviving realization.
 *
 * Only seq_types that have actually been written are present in the maps, so a missing key
 * is itself an observation.
 */
struct ScenarioSnapshot {
    double scenario_proba = 0.0;
    std::map<Seq_type, std::pair<Seq_Offset, Seq_Offset>> offsets;
    std::map<Seq_type, std::string> sequences; ///< decoded back to ACGT
    std::map<Seq_type, std::vector<std::size_t>> mismatches;
    std::map<Event_safety, bool> safety;
    std::map<Seq_type, double> downstream_bounds;

    /// Convenience: the 5' offset of a segment, or throws if it was never written.
    Seq_Offset five_prime(Seq_type seq_type) const { return offsets.at(seq_type).first; }
    Seq_Offset three_prime(Seq_type seq_type) const { return offsets.at(seq_type).second; }
};

/**
 * @brief A stand-in for "the next event", which records instead of recursing.
 *
 * This is the primary observation instrument. Testing a single event through the leaf path
 * instead would drag in Error_rate, which needs a *complete* scenario (it concatenates V,
 * D and J) and therefore says nothing about the event under test. Recording at the
 * hand-off point gives one snapshot per realization that survived every check, which is
 * exactly the quantity these branches decide.
 */
class RecordingEvent : public Rec_Event {
public:
    explicit RecordingEvent(int event_id);

    std::vector<ScenarioSnapshot> calls;

    /// Number of realizations that reached the next event.
    std::size_t call_count() const { return calls.size(); }

    void iterate(QuerySequenceContext &query, const ModelContext &model, ScenarioContext &scenario,
                 ExplorationContext &exploration, AccumulationContext &accumulation) override;

    //Unused surface, present only because Rec_Event declares it pure virtual.
    std::shared_ptr<Rec_Event> copy() override;
    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name,
                                     std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const override;
    void write2txt(std::ofstream &) override {}
    void write2txt_legacy(std::ofstream &) override {}
    void write2txt_v2(std::ofstream &) override {}
    void add_to_marginals(long double, Marginal_array_p &) const override {}
    bool has_effect_on(Seq_type) const override { return false; }
    void iterate_initialize_Len_proba(Seq_type, std::map<int, double> &,
                                      std::queue<std::shared_ptr<Rec_Event>> &, double &,
                                      const Marginal_array_p &, Index_map &, Seq_type_str_p_map &,
                                      int &) const override
    {}
    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &, const Marginal_array_p &,
                                    Index_map &) override
    {}
};

/**
 * @brief Run the production initialization sequence, then call iterate() on `event`.
 *
 * Mirrors GenModel::infer_model / compute_Pgen ordering, which matters: the length-proba
 * bounds are built by a *reverse* pass over the queue after every event is initialized,
 * and skipping it leaves vd/dj/vj_length_best_proba_map empty, which silently discards
 * every scenario at the junction-length guard.
 *
 * `event` is initialized last, so anything marked with mark_chosen() is already in
 * processed_events when its initialize_event() runs.
 */
void call_iterate(const std::shared_ptr<Rec_Event> &event, IterateTestState &state,
                  const std::shared_ptr<RecordingEvent> &next = nullptr);

/// Build a recorder, wire it as `event`'s successor, run call_iterate(), return it.
std::shared_ptr<RecordingEvent> call_iterate_recording(const std::shared_ptr<Rec_Event> &event,
                                                       IterateTestState &state);

// ============================================================================
// State inspection
// ============================================================================

Seq_Offset get_seq_offset(const IterateTestState &state, Seq_type seq_type, Seq_side side,
                          std::size_t layer = 0);
bool has_seq_offset(const IterateTestState &state, Seq_type seq_type, Seq_side side);

const Int_Str *get_constructed_sequence(const IterateTestState &state, Seq_type seq_type,
                                        std::size_t layer = 0);
bool has_constructed_sequence(const IterateTestState &state, Seq_type seq_type);

std::vector<std::size_t> get_mismatches(const IterateTestState &state, Seq_type seq_type,
                                        std::size_t layer = 0);

bool is_safe(const IterateTestState &state, Event_safety safety_type, std::size_t layer = 0);
bool has_safety(const IterateTestState &state, Event_safety safety_type);

double get_downstream_bound(const IterateTestState &state, Seq_type seq_type,
                            std::size_t layer = 0);

/// Sum of the whole updated-marginals array. Zero means no scenario reached a leaf, which
/// is the cheapest positive/negative control a section can assert on.
long double total_marginal_mass(const IterateTestState &state, std::size_t marginal_array_size = 1000);

// ============================================================================
// Event builders
// ============================================================================

/**
 * Build a Gene_choice carrying one realization per (name, sequence) pair.
 *
 * `fixed` defaults to true for a reason: at a leaf, iterate_wrap_up() calls
 * add_to_marginals() on **every** non-fixed event in events_map, and a stub that never
 * iterated still has new_index == -1, so an unfixed stub writes out of bounds. Only the
 * event under test should be unfixed.
 */
std::shared_ptr<Gene_choice> make_gene_choice(Gene_class gene_class,
                                              const std::vector<std::pair<std::string, std::string>> &genes,
                                              int event_id, bool fixed = true);

/// Build a Deletion event over the inclusive realization range [min_del, max_del].
/// Present purely so a neighbour has non-zero deletion bounds; always fixed.
std::shared_ptr<Deletion> make_deletion(Seq_type target, Seq_side side, int min_del, int max_del,
                                        int event_id);

Alignment_data create_perfect_alignment(const std::string &gene_name, int offset, int gene_length);

Alignment_data create_alignment_with_mismatches(const std::string &gene_name, int offset,
                                                int gene_length,
                                                const std::vector<std::size_t> &mismatch_positions);

} // namespace IgorTestUtils

#include <igor/Core/SeqTypeRegistry.h>

/**
 * A frozen registry holding the six legacy seq_types at their Seq_type enum ids, with no
 * ordering set. This is what the scenario maps are sized from in tests that only need
 * enum-keyed access; tests exercising the ordered traversal build their own registry with
 * an explicit ordering.
 */
//legacy_seq_type_registry() now lives in SeqTypeRegistry.h; kept included here so that
//tests including test_utils.h continue to see it.
