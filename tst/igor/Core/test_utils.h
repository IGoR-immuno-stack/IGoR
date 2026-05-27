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

#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <memory>
#include <string>
#include <vector>
#include <random>

namespace IgorTestUtils {

/**
 * @brief Load a simple model from test model files
 * 
 * Loads the appropriate model file (with or without D gene) from tst/models/
 * and creates marginals from the model structure.
 * 
 * @param with_D If true, loads monoscenario_wD_parms.txt, otherwise monoscenario_noD_parms.txt
 * @return std::pair<std::shared_ptr<Model_Parms>, std::shared_ptr<Model_marginals>>
 */
std::pair<std::shared_ptr<Model_Parms>, std::shared_ptr<Model_marginals>>
load_simple_model(bool with_D = true);

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
    const std::vector<int>& mismatches = {},
    double score = 100.0
);

/**
 * @brief Validate that marginal probabilities are properly normalized
 * 
 * @param marginals Model marginals to validate
 * @param model_parms Associated model parameters
 * @param tolerance Tolerance for floating point comparison (default 1e-9)
 * @return true if marginals are valid, false otherwise
 */
bool validate_marginals(const Model_marginals& marginals,
                       const Model_Parms& model_parms,
                       double tolerance = 1e-9);

/**
 * @brief Check if two sequences are equal
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @return true if sequences match, false otherwise
 */
bool compare_sequences(const std::string& seq1, const std::string& seq2);

/**
 * @brief Check if two sequences match with allowed mismatches
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param max_mismatches Maximum allowed mismatches
 * @return true if sequences match within tolerance, false otherwise
 */
bool compare_sequences_with_mismatches(const std::string& seq1,
                                      const std::string& seq2,
                                      int max_mismatches);

// ============================================================================
// Storage Objects (Private implementation - own data, guarantee Context lifetime)
// ============================================================================

/**
 * @brief Storage for QuerySequenceContext fields
 * 
 * Private storage that owns query data and guarantees lifetime for QuerySequenceContext references.
 * Update this struct when QuerySequenceContext evolves.
 */
struct QueryStorage {
    std::string sequence;
    Int_Str int_sequence;
    std::unordered_map<Gene_class, std::vector<Alignment_data>> gene_alignments;
    
    explicit QueryStorage(const std::string& seq)
        : sequence(seq), int_sequence(nt2int(seq)), gene_alignments{} {}
};

/**
 * @brief Storage for ModelContext fields
 * 
 * Private storage that owns model data and guarantees lifetime for ModelContext references.
 * Update this struct when ModelContext evolves.
 */
struct ModelStorage {
    std::unique_ptr<long double[]> model_marginals;
    std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map;
    std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>> events_map;
    std::queue<std::shared_ptr<Rec_Event>> model_queue;
    
    explicit ModelStorage(size_t marginal_array_size)
        : model_marginals(std::make_unique<long double[]>(marginal_array_size)),
          offset_map{}, events_map{}, model_queue{} {
        // Initialize marginal array to 0
        for (size_t i = 0; i < marginal_array_size; ++i) {
            model_marginals[i] = 0.0;
        }
    }
};

/**
 * @brief Storage for ScenarioContext fields
 * 
 * Private storage that owns scenario data and guarantees lifetime for ScenarioContext references.
 * Update this struct when ScenarioContext evolves.
 */
struct ScenarioStorage {
    double scenario_proba;
    Seq_type_str_p_map constructed_sequences;
    Seq_offsets_map seq_offsets;
    Mismatch_vectors_map mismatches_lists;
    
    ScenarioStorage()
        : scenario_proba(1.0),
          constructed_sequences(Seq_type_str_p_map(6)),
          seq_offsets(Seq_offsets_map(6, 3)),
          mismatches_lists(Mismatch_vectors_map(6)) {}
};

/**
 * @brief Storage for ExplorationContext fields
 * 
 * Private storage that owns exploration data and guarantees lifetime for ExplorationContext references.
 * Update this struct when ExplorationContext evolves.
 */
struct ExplorationStorage {
    Downstream_scenario_proba_bound_map downstream_proba_map;
    double seq_max_prob;
    double proba_threshold;
    Index_map base_index_map;
    std::shared_ptr<Next_event_ptr> next_event_ptr_arr;
    Safety_bool_map safety_set;
    
    ExplorationStorage()
        : downstream_proba_map(Downstream_scenario_proba_bound_map(6)),
          seq_max_prob(1.0),
          proba_threshold(1e-10),
          base_index_map(Index_map(100)),
          next_event_ptr_arr(nullptr),
          safety_set(Safety_bool_map(3)) {}
};

/**
 * @brief Storage for AccumulationContext fields
 * 
 * Private storage that owns accumulation data and guarantees lifetime for AccumulationContext references.
 * Update this struct when AccumulationContext evolves.
 */
struct AccumulationStorage {
    std::unique_ptr<long double[]> updated_marginals;
    std::map<size_t, std::shared_ptr<Counter>> counters_list;
    std::shared_ptr<Error_rate> error_rate;
    
    explicit AccumulationStorage(size_t marginal_array_size)
        : updated_marginals(std::make_unique<long double[]>(marginal_array_size)),
          counters_list{},
          error_rate(std::make_shared<Single_error_rate>(0.0)) {
        // Initialize marginal array to 0
        for (size_t i = 0; i < marginal_array_size; ++i) {
            updated_marginals[i] = 0.0;
        }
    }
};

/**
 * @brief State container for iterate() method testing (storage-backed context design)
 * 
 * Clean separation of concerns:
 * - Private Storage objects own data and guarantee Context reference lifetime
 * - Public Context objects provide the API interface
 * - Tests interact exclusively through Context objects
 * - When a Context evolves, update the corresponding Storage struct
 * 
 * Design contract:
 * 1. Storage objects are private, guarantee lifetime, have constructors
 * 2. Context objects are public, reference Storage data
 * 3. All test access goes through Context API (state.query.*, state.model.*, etc.)
 * 4. Constructor: construct Storage → wire references to Contexts
 * 
 * Example usage:
 *   auto state = create_iterate_state("ACGTACGT");
 *   state.query.gene_alignments[V_gene] = v_alignments;
 *   state.model.events_map[{GeneChoice_t, V_gene, Undefined_side}] = v_event;
 *   call_iterate(event, state);
 */
struct IterateTestState {
private:
    // ========================================================================
    // Storage objects (own data, guarantee Context lifetime)
    // These are implementation details - tests never access them directly
    // ========================================================================
    QueryStorage query_storage;
    ModelStorage model_storage;
    ScenarioStorage scenario_storage;
    ExplorationStorage exploration_storage;
    AccumulationStorage accumulation_storage;
    
public:
    // ========================================================================
    // Context objects (public interface - reference Storage data)
    // Tests access everything through these
    // ========================================================================
    QuerySequenceContext query;
    ModelContext model;
    ScenarioContext scenario;
    ExplorationContext exploration;
    AccumulationContext accumulation;
    
    /**
     * @brief Constructor - initialize Storage objects, then wire Contexts to them
     * 
     * Construction order:
     * 1. Construct Storage objects (with initialization logic in their constructors)
     * 2. Construct Context objects that reference Storage data
     * 
     * This encapsulates all initialization complexity in Storage constructors,
     * keeping this constructor clean and focused on wiring.
     */
    IterateTestState(
        const std::string& seq,
        size_t marginal_array_size
    ) : // Step 1: Construct Storage objects
        query_storage(seq),
        model_storage(marginal_array_size),
        scenario_storage(),
        exploration_storage(),
        accumulation_storage(marginal_array_size),
        
        // Step 2: Wire Contexts to Storage references
        query(query_storage.sequence, query_storage.int_sequence, query_storage.gene_alignments),
        model(model_storage.model_marginals, model_storage.offset_map, model_storage.events_map, model_storage.model_queue),
        scenario(scenario_storage.scenario_proba, scenario_storage.constructed_sequences, 
                scenario_storage.seq_offsets, scenario_storage.mismatches_lists),
        exploration(exploration_storage.downstream_proba_map, exploration_storage.seq_max_prob,
                   exploration_storage.proba_threshold, exploration_storage.base_index_map,
                   exploration_storage.next_event_ptr_arr, exploration_storage.safety_set),
        accumulation(accumulation_storage.updated_marginals, accumulation_storage.counters_list,
                    accumulation_storage.error_rate)
    {}
    
    // Prevent copying (contexts hold references)
    IterateTestState(const IterateTestState&) = delete;
    IterateTestState& operator=(const IterateTestState&) = delete;
    
    // Allow moving
    IterateTestState(IterateTestState&&) = default;
    IterateTestState& operator=(IterateTestState&&) = delete;
};

// ============================================================================
// Test State Factory
// ============================================================================

/**
 * @brief Create an IterateTestState with reasonable defaults
 * 
 * Creates a test state initialized for iterate() method testing.
 * All data is owned by the returned state, and tests interact with it
 * through the Context object members (query, model, scenario, etc.).
 * 
 * After creation, populate test-specific data through the Context objects:
 *   state.query.gene_alignments[V_gene] = v_alignments;
 *   state.model.events_map[{GeneChoice_t, V_gene, Undefined_side}] = v_event;
 * 
 * @param sequence The test sequence to evaluate
 * @param marginal_array_size Size of marginal array to allocate (default 1000)
 * @return IterateTestState Initialized state for testing
 */
IterateTestState create_iterate_state(
    const std::string& sequence,
    size_t marginal_array_size = 1000
);

/**
 * @brief Wrapper function to call iterate() on a Rec_Event with context-based interface
 * 
 * This helper function simplifies calling iterate() by taking an IterateTestState struct,
 * initializing the event, setting up post-initialization for uninitialized genes, and
 * dispatching context objects to the iterate() method.
 * 
 * Post-initialization pattern: After initialize_event(), this function initializes empty
 * sequences/mismatch lists for genes not set by the event. This prevents null pointer
 * dereferencing in error rate calculation when testing events in isolation.
 * 
 * @param event The recombination event to iterate
 * @param state The test state containing all iterate parameters and storage
 */
void call_iterate(
    std::shared_ptr<Rec_Event> event,
    IterateTestState& state
);

// ============================================================================
// State Inspection Helpers
// ============================================================================

/**
 * @brief Get sequence offset from IterateTestState
 * 
 * Helper that accesses scenario context data.
 * 
 * @param state Test state
 * @param seq_type Type of sequence (V_gene_seq, D_gene_seq, etc.)
 * @param side Five_prime or Three_prime
 * @param layer Memory layer (default 0)
 * @return Offset value
 */
int get_seq_offset(
    const IterateTestState& state,
    Seq_type seq_type,
    Seq_side side,
    int layer = 0
);

/**
 * @brief Check if a seq_offset exists in the map
 * 
 * Helper that accesses scenario context data.
 * 
 * @param state Test state
 * @param seq_type Sequence type (V_gene_seq, D_gene_seq, J_gene_seq)
 * @param side Sequence side (Five_prime, Three_prime)
 * @param layer Memory layer (default 0)
 * @return true if the offset exists, false otherwise
 */
bool has_seq_offset(
    const IterateTestState& state,
    Seq_type seq_type,
    Seq_side side,
    int layer = 0
);

/**
 * @brief Get constructed sequence from IterateTestState
 * 
 * Helper that accesses scenario context data.
 * 
 * @param state Test state
 * @param seq_type Type of sequence
 * @param layer Memory layer (default 0)
 * @return Pointer to constructed sequence, or nullptr if not set
 */
const Int_Str* get_constructed_sequence(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer = 0
);

/**
 * @brief Check if a constructed sequence exists in the map
 * 
 * Helper that accesses scenario context data.
 * 
 * @param state Test state
 * @param seq_type Sequence type (V_gene_seq, D_gene_seq, J_gene_seq)
 * @param layer Memory layer (default 0)
 * @return true if the constructed sequence exists, false otherwise
 */
bool has_constructed_sequence(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer = 0
);

/**
 * @brief Get mismatch positions for a sequence type
 * 
 * Helper that accesses scenario context data.
 * 
 * @param state Test state
 * @param seq_type Type of sequence
 * @param layer Memory layer (default 0)
 * @return Vector of mismatch positions, empty if none
 */
std::vector<int> get_mismatches(
    const IterateTestState& state,
    Seq_type seq_type,
    int layer = 0
);

/**
 * @brief Check if a safety flag is set
 * 
 * Helper that accesses exploration context data.
 * 
 * @param state Test state
 * @param safety_type VD_safe, VJ_safe, or DJ_safe
 * @param layer Memory layer (default 0)
 * @return true if safety condition is met
 */
bool is_safe(
    const IterateTestState& state,
    Event_safety safety_type,
    int layer = 0
);

/**
 * @brief Check if a safety flag exists in the map
 * 
 * Helper that accesses exploration context data.
 * 
 * @param state Test state
 * @param safety_type VD_safe, VJ_safe, or DJ_safe
 * @param layer Memory layer (default 0)
 * @return true if the safety flag exists, false otherwise
 */
bool has_safety(
    const IterateTestState& state,
    Event_safety safety_type,
    int layer = 0
);

// ============================================================================
// Stub Event Creators
// ============================================================================

/**
 * @brief Create a stub Gene_choice event for testing
 * 
 * Creates a minimal Gene_choice with single realization that sets offsets
 * but doesn't perform full iterate logic. Useful for testing events that
 * depend on other genes being present.
 * 
 * @param gene_class V_gene, D_gene, or J_gene
 * @param gene_name Name of the gene (e.g., "TRBV1")
 * @param gene_seq Gene sequence
 * @param event_id Event identifier
 * @return Shared pointer to stub Gene_choice event
 */
std::shared_ptr<Gene_choice> create_stub_gene_choice(
    Gene_class gene_class,
    const std::string& gene_name,
    const std::string& gene_seq,
    int event_id
);

/**
 * @brief Create a minimal events_map with V, D, and J stubs
 * 
 * @param v_event V gene choice event (can be real or stub)
 * @param d_stub D gene choice stub (optional, can be nullptr)
 * @param j_stub J gene choice stub (optional, can be nullptr)
 * @return events_map suitable for testing
 */
std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>
create_events_map(
    std::shared_ptr<Gene_choice> v_event,
    std::shared_ptr<Gene_choice> d_stub = nullptr,
    std::shared_ptr<Gene_choice> j_stub = nullptr
);

// ============================================================================
// Alignment Helpers
// ============================================================================

/**
 * @brief Create a perfect alignment (no mismatches) at offset 0
 * 
 * @param gene_name Name of the gene
 * @param gene_length Length of the gene sequence
 * @return Alignment_data for perfect match
 */
Alignment_data create_perfect_alignment(
    const std::string& gene_name,
    int gene_length
);

/**
 * @brief Create alignment with specific mismatch positions
 * 
 * @param gene_name Name of the gene
 * @param offset Alignment offset
 * @param gene_length Length of the gene
 * @param mismatch_positions Vector of 0-based mismatch positions
 * @return Alignment_data with specified mismatches
 */
Alignment_data create_alignment_with_mismatches(
    const std::string& gene_name,
    int offset,
    int gene_length,
    const std::vector<int>& mismatch_positions
);

/**
 * @brief Extract marginals for a specific event from marginal array
 * 
 * @param marginals Full marginal array
 * @param event Event to extract marginals for
 * @return std::vector<double> Event-specific marginal probabilities
 */
bool compare_sequences_with_tolerance(const std::string& seq1, 
                                      const std::string& seq2,
                                      int max_mismatches);

/**
 * @brief Count mismatches between two sequences
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @return Number of mismatches
 */
int count_mismatches(const std::string& seq1, const std::string& seq2);

/**
 * @brief Validate that a probability distribution sums to 1.0
 * 
 * @param probabilities Vector of probabilities
 * @param tolerance Tolerance for floating point comparison
 * @return true if distribution is normalized, false otherwise
 */
bool is_normalized(const std::vector<double>& probabilities, double tolerance = 1e-9);

/**
 * @brief Mock random generator for deterministic testing
 * 
 * Allows controlling random number generation by providing predefined values
 */
class MockRandomGenerator {
public:
    /**
     * @brief Construct mock generator with predefined values
     * 
     * @param values Vector of predetermined random values (0.0 to 1.0)
     */
    explicit MockRandomGenerator(const std::vector<double>& values);
    
    /**
     * @brief Get next random value
     * 
     * @return double Next value in the sequence (cycles if exhausted)
     */
    double next();
    
    /**
     * @brief Reset to beginning of sequence
     */
    void reset();
    
    /**
     * @brief Get current index in sequence
     */
    size_t get_index() const { return current_index_; }
    
private:
    std::vector<double> values_;
    size_t current_index_;
};

/**
 * @brief Extract marginal probabilities for a specific event
 * 
 * @param marginals Model marginals
 * @param event_name Name of the event
 * @param model_parms Model parameters
 * @return std::vector<double> Vector of marginal probabilities for each realization
 */
std::vector<double> extract_event_marginals(
    const Model_marginals& marginals,
    const std::string& event_name,
    const Model_Parms& model_parms
);

/**
 * @brief Create a simple sequence with known structure for testing
 * 
 * @param v_gene V gene sequence
 * @param v_del Number of V gene 3' deletions
 * @param vd_ins VD insertion sequence
 * @param d_gene D gene sequence
 * @param d5_del Number of D gene 5' deletions
 * @param d3_del Number of D gene 3' deletions
 * @param dj_ins DJ insertion sequence
 * @param j_del Number of J gene 5' deletions
 * @param j_gene J gene sequence
 * @return std::string Constructed sequence
 */
std::string create_test_sequence(
    const std::string& v_gene, int v_del,
    const std::string& vd_ins,
    const std::string& d_gene, int d5_del, int d3_del,
    const std::string& dj_ins,
    int j_del, const std::string& j_gene
);

} // namespace IgorTestUtils
