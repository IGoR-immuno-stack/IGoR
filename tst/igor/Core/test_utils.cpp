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
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <functional>

namespace IgorTestUtils {

std::pair<std::shared_ptr<Model_Parms>, std::shared_ptr<Model_marginals>>
load_simple_model(bool with_D) {
    // Determine model file path based on with_D flag
    // Tests run with WORKING_DIRECTORY set to tst/ source directory
    std::string model_file = with_D ? 
        "models/monoscenario_wD_parms.txt" : 
        "models/monoscenario_noD_parms.txt";
    
    // Load model parameters
    auto model_parms = std::make_shared<Model_Parms>();
    model_parms->read_model_parms(model_file);
    
    // Create marginals from the model structure
    auto model_marginals = std::make_shared<Model_marginals>(*model_parms);
    
    // Initialize with uniform distribution
    model_marginals->uniform_initialize(*model_parms);
    
    return {model_parms, model_marginals};
}

Alignment_data create_mock_alignment_data(
    const std::string& gene_name,
    int offset,
    size_t five_p_offset,
    size_t three_p_offset,
    const std::vector<int>& mismatches,
    double score
) {
    // Calculate alignment length
    size_t align_length = three_p_offset - five_p_offset;
    
    // Create empty forward_lists for insertions and deletions
    std::forward_list<int> empty_insertions;
    std::forward_list<int> empty_deletions;
    
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

bool validate_marginals(const Model_marginals& marginals,
                       const Model_Parms& model_parms,
                       double tolerance) {
    // Get the marginal array
    const auto& marginal_array = marginals.marginal_array_smart_p;
    size_t array_size = marginals.get_length();
    
    if (array_size == 0) {
        return false;
    }
    
    // Check that all probabilities are non-negative
    for (size_t i = 0; i < array_size; ++i) {
        if (marginal_array[i] < 0.0) {
            return false;
        }
    }
    
    // For a more thorough check, we would need to iterate through each event
    // and verify that marginals sum to 1.0 for each realization combination
    // This is a simplified check that the array exists and has valid values
    
    return true;
}

bool compare_sequences(const std::string& seq1, const std::string& seq2) {
    return seq1 == seq2;
}

bool compare_sequences_with_tolerance(const std::string& seq1, 
                                      const std::string& seq2,
                                      int max_mismatches) {
    if (seq1.length() != seq2.length()) {
        return false;
    }
    
    int mismatches = count_mismatches(seq1, seq2);
    return mismatches <= max_mismatches;
}

int count_mismatches(const std::string& seq1, const std::string& seq2) {
    if (seq1.length() != seq2.length()) {
        throw std::invalid_argument("Sequences must have equal length");
    }
    
    int mismatches = 0;
    for (size_t i = 0; i < seq1.length(); ++i) {
        if (seq1[i] != seq2[i]) {
            ++mismatches;
        }
    }
    
    return mismatches;
}

bool is_normalized(const std::vector<double>& probabilities, double tolerance) {
    double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    return std::abs(sum - 1.0) < tolerance;
}

// MockRandomGenerator implementation
MockRandomGenerator::MockRandomGenerator(const std::vector<double>& values)
    : values_(values), current_index_(0) {
    if (values_.empty()) {
        throw std::invalid_argument("MockRandomGenerator requires at least one value");
    }
    
    // Validate that all values are in [0, 1]
    for (double val : values_) {
        if (val < 0.0 || val > 1.0) {
            throw std::invalid_argument("All values must be in range [0.0, 1.0]");
        }
    }
}

double MockRandomGenerator::next() {
    double value = values_[current_index_];
    current_index_ = (current_index_ + 1) % values_.size();
    return value;
}

void MockRandomGenerator::reset() {
    current_index_ = 0;
}

std::vector<double> extract_event_marginals(
    const Model_marginals& marginals,
    const std::string& event_name,
    const Model_Parms& model_parms
) {
    // Get event pointer
    auto event_ptr = model_parms.get_event_pointer(event_name, false);
    if (!event_ptr) {
        throw std::runtime_error("Event not found: " + event_name);
    }
    
    // Compute marginal probability for this event
    auto marginal_result = marginals.compute_event_marginal_probability(
        event_ptr->get_name(), model_parms);
    
    // Extract the probability array
    const auto& prob_array = marginal_result.second;
    
    // Get the number of realizations for this event
    // Use len_max - len_min + 1 as the size (for events like Insertion)
    // For Gene_choice events, this would be the number of genes
    size_t num_realizations = static_cast<size_t>(event_ptr->get_len_max() - event_ptr->get_len_min() + 1);
    
    // Extract probabilities into vector
    std::vector<double> probs(num_realizations);
    for (size_t i = 0; i < num_realizations; ++i) {
        probs[i] = static_cast<double>(prob_array.get()[i]);
    }
    
    return probs;
}

std::string create_test_sequence(
    const std::string& v_gene, int v_del,
    const std::string& vd_ins,
    const std::string& d_gene, int d5_del, int d3_del,
    const std::string& dj_ins,
    int j_del, const std::string& j_gene
) {
    std::string sequence;
    
    // Add V gene (with 3' deletions)
    if (v_del < static_cast<int>(v_gene.length())) {
        sequence += v_gene.substr(0, v_gene.length() - v_del);
    }
    
    // Add VD insertion
    sequence += vd_ins;
    
    // Add D gene (with 5' and 3' deletions)
    if (d5_del + d3_del < static_cast<int>(d_gene.length())) {
        sequence += d_gene.substr(d5_del, d_gene.length() - d5_del - d3_del);
    }
    
    // Add DJ insertion
    sequence += dj_ins;
    
    // Add J gene (with 5' deletions)
    if (j_del < static_cast<int>(j_gene.length())) {
        sequence += j_gene.substr(j_del);
    }
    
    return sequence;
}

bool compare_sequences_with_mismatches(const std::string& seq1,
                                      const std::string& seq2,
                                      int max_mismatches) {
    if (seq1.length() != seq2.length()) {
        return false;
    }
    
    int mismatches = 0;
    for (size_t i = 0; i < seq1.length(); ++i) {
        if (seq1[i] != seq2[i]) {
            ++mismatches;
            if (mismatches > max_mismatches) {
                return false;
            }
        }
    }
    
    return true;
}

// IterateTestState implementation - uses context factory functions

// ============================================================================
// Test State Factory
// ============================================================================

IterateTestState create_iterate_state(
    const std::string& sequence,
    size_t marginal_array_size
) {
    // Constructor handles all initialization
    return IterateTestState(sequence, marginal_array_size);
}

void call_iterate(
    std::shared_ptr<Rec_Event> event,
    IterateTestState& state
) {
    // Set base_index_map for ALL events in events_map (required by IGoR)
    // Access through const_cast since model context provides const reference
    auto& events_map = const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(
        state.model.events_map
    );
    for (const auto& event_pair : events_map) {
        const auto& ev = event_pair.second;
        if (ev) {
            // Access through const_cast since exploration context provides non-const reference
            const_cast<Index_map&>(state.exploration.index_map)[ev->get_event_identifier()] = 0;
        }
    }
    
    // Initialize the event - this will call initialize_event() which requests memory layers
    std::unordered_set<Rec_Event_name> initialized_events;
    event->initialize_event(
        initialized_events,
        events_map,
        const_cast<std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>&>(
            state.model.offset_map
        ),
        const_cast<Downstream_scenario_proba_bound_map&>(state.exploration.downstream_proba_map),
        const_cast<Seq_type_str_p_map&>(state.scenario.constructed_sequences),
        const_cast<Safety_bool_map&>(state.exploration.safety_set),
        const_cast<std::shared_ptr<Error_rate>&>(state.accumulation.error_rate),
        const_cast<Mismatch_vectors_map&>(state.scenario.mismatches_lists),
        const_cast<Seq_offsets_map&>(state.scenario.seq_offsets),
        const_cast<Index_map&>(state.exploration.index_map)
    );
    
    // POST-INITIALIZATION: Initialize empty sequences and mismatch lists for genes 
    // not initialized by the event. This prevents null pointer dereferencing in 
    // error rate calculation (compare_sequences) when testing events in isolation.
    // Use static storage to ensure pointers remain valid throughout test execution.
    static Int_Str empty_seq = nt2int("");
    static std::vector<int> empty_mismatches;
    
    auto& constructed_seqs = const_cast<Seq_type_str_p_map&>(state.scenario.constructed_sequences);
    auto& mismatch_lists = const_cast<Mismatch_vectors_map&>(state.scenario.mismatches_lists);
    
    if (!constructed_seqs.exist(V_gene_seq)) {
        constructed_seqs.set_value(V_gene_seq, &empty_seq, 0);
        mismatch_lists.set_value(V_gene_seq, &empty_mismatches, 0);
    }
    if (!constructed_seqs.exist(J_gene_seq)) {
        constructed_seqs.set_value(J_gene_seq, &empty_seq, 0);
        mismatch_lists.set_value(J_gene_seq, &empty_mismatches, 0);
    }
    
    // Create a next_event_ptr array for the event chain
    const size_t max_events = 100;
    auto& next_ptr = const_cast<std::shared_ptr<Next_event_ptr>&>(state.exploration.next_event_ptr_arr);
    next_ptr = std::shared_ptr<Next_event_ptr>(
        new Next_event_ptr[max_events](),
        std::default_delete<Next_event_ptr[]>()
    );
    
    // Call iterate with context objects (already members of state)
    event->iterate(state.query, state.model, state.scenario, state.exploration, state.accumulation);
}

// ============================================================================
// State Inspection Helpers
// ============================================================================

int get_seq_offset(
    const IterateTestState& state,
    Seq_type seq_type,
    Seq_side side,
    int layer
) {
    return state.scenario.seq_offsets.at(seq_type, side, layer);
}

bool has_seq_offset(
    const IterateTestState& state,
    Seq_type seq_type,
    Seq_side side,
    int layer
) {
    try {
        state.scenario.seq_offsets.at(seq_type, side, layer);
        return true;
    } catch (const std::out_of_range&) {
        return false;
    }
}

const Int_Str* get_constructed_sequence(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer
) {
    return state.scenario.constructed_sequences.at(seq_type, layer);
}

bool has_constructed_sequence(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer
) {
    try {
        state.scenario.constructed_sequences.at(seq_type, layer);
        return true;
    } catch (const std::out_of_range&) {
        return false;
    }
}

std::vector<int> get_mismatches(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer
) {
    auto* mismatch_vec = state.scenario.mismatches_lists.at(seq_type, layer);
    if (mismatch_vec) {
        return *mismatch_vec;
    }
    return {};
}

bool is_safe(
    const IterateTestState& state,
    Event_safety safety_type,
    int layer
) {
    return state.exploration.safety_set.at(safety_type, layer);
}

bool has_safety(
    const IterateTestState& state,
    Event_safety safety_type,
    int layer
) {
    try {
        state.exploration.safety_set.at(safety_type, layer);
        return true;
    } catch (const std::out_of_range&) {
        return false;
    }
}

// ============================================================================
// Stub Event Creators
// ============================================================================

std::shared_ptr<Gene_choice> create_stub_gene_choice(
    Gene_class gene_class,
    const std::string& gene_name,
    const std::string& gene_seq,
    int event_id
) {
    auto stub = std::make_shared<Gene_choice>(gene_class);
    stub->add_realization(gene_name, gene_seq);
    stub->set_event_identifier(event_id);
    stub->set_priority(1);
    return stub;
}

std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>
create_events_map(
    std::shared_ptr<Gene_choice> v_event,
    std::shared_ptr<Gene_choice> d_stub,
    std::shared_ptr<Gene_choice> j_stub
) {
    std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> events_map;
    
    if (v_event) {
        events_map[std::make_tuple(GeneChoice_t, V_gene, Undefined_side)] = v_event;
    }
    if (d_stub) {
        events_map[std::make_tuple(GeneChoice_t, D_gene, Undefined_side)] = d_stub;
    }
    if (j_stub) {
        events_map[std::make_tuple(GeneChoice_t, J_gene, Undefined_side)] = j_stub;
    }
    
    return events_map;
}

// ============================================================================
// Alignment Helpers
// ============================================================================

Alignment_data create_perfect_alignment(
    const std::string& gene_name,
    int gene_length
) {
    return create_mock_alignment_data(
        gene_name,
        0,              // offset
        0,              // 5' offset
        gene_length - 1,// 3' offset (0-indexed)
        {},             // no mismatches
        100.0           // perfect score
    );
}

Alignment_data create_alignment_with_mismatches(
    const std::string& gene_name,
    int offset,
    int gene_length,
    const std::vector<int>& mismatch_positions
) {
    return create_mock_alignment_data(
        gene_name,
        offset,
        offset >= 0 ? offset : 0,
        offset + gene_length - 1,
        mismatch_positions,
        100.0 - mismatch_positions.size() * 5.0  // Rough score
    );
}

// Note: Next_event_ptr is just a typedef for Rec_Event*
// For unit tests, we can simply pass nullptr since we're testing isolated events

} // namespace IgorTestUtils
