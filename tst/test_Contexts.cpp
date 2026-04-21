#include <catch2/catch_test_macros.hpp>

#include <igor/Core/QuerySequenceContext.h>
#include <igor/Core/ModelContext.h>
#include <igor/Core/ScenarioContext.h>
#include <igor/Core/ExplorationContext.h>
#include <igor/Core/AccumulationContext.h>

#include <igor/Core/Utils.h>
#include <igor/Core/IntStr.h>
#include <igor/Core/Aligner.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Counter.h>
#include <igor/Core/Errorrate.h>
#include <igor/Core/Singleerrorrate.h>

#include "test_utils.h"

#include <string>
#include <memory>

using namespace std;
using namespace IgorTestUtils;

// ============================================================================
// QuerySequenceContext Tests
// ============================================================================

TEST_CASE("QuerySequenceContext construction and reference binding", "[Context][QuerySequenceContext]") {
    string seq = "ACGTACGT";
    Int_Str int_seq = nt2int(seq);
    
    // Mock alignment data
    unordered_map<Gene_class, vector<Alignment_data>> alignments;
    Alignment_data v_align = create_mock_alignment_data(
            "TRBV7-2*01", -250, 0, 16, {7, 9, 12}, 95.0
        );
    Alignment_data j_align = create_mock_alignment_data(
            "TRBJ1-1*01", 60, 60, 70, {}, 30.0
        );
    alignments[V_gene] = {v_align};
    alignments[J_gene] = {j_align};
    
    SECTION("Basic construction") {
        QuerySequenceContext query(seq, int_seq, alignments);
        
        // Verify references are bound correctly
        REQUIRE(query.sequence == seq);
        REQUIRE(&query.sequence == &seq);  // Same address
        REQUIRE(&query.int_sequence == &int_seq);
        REQUIRE(&query.gene_alignments == &alignments);
    }
    
    SECTION("Access alignment data") {
        QuerySequenceContext query(seq, int_seq, alignments);
        
        REQUIRE(query.gene_alignments.at(V_gene).size() == 1);
        REQUIRE(query.gene_alignments.at(V_gene)[0].gene_name == "TRBV7-2*01");
        REQUIRE(query.gene_alignments.at(J_gene).size() == 1);
        REQUIRE(query.gene_alignments.at(J_gene)[0].gene_name == "TRBJ1-1*01");
    }
    
    SECTION("Multiple alignments per gene") {
        Alignment_data v_align2 = create_mock_alignment_data(
            "TRBV7-3*01", -250, 0, 16, {7, 9, 12}, 95.0
        );

        alignments[V_gene].push_back(v_align2);
        
        QuerySequenceContext query(seq, int_seq, alignments);
        
        REQUIRE(query.gene_alignments.at(V_gene).size() == 2);
        REQUIRE(query.gene_alignments.at(V_gene)[1].gene_name == "TRBV7-3*01");
    }
    
    SECTION("Empty alignments") {
        unordered_map<Gene_class, vector<Alignment_data>> empty_alignments;
        
        QuerySequenceContext query(seq, int_seq, empty_alignments);
        
        REQUIRE(query.gene_alignments.empty());
    }
}

// ============================================================================
// ModelContext Tests
// ============================================================================

TEST_CASE("ModelContext construction and immutability", "[Context][ModelContext]") {
    // Create mock model data
    auto model_params = make_unique<long double[]>(100);
    for (int i = 0; i < 100; ++i) {
        model_params[i] = static_cast<long double>(i) / 100.0;
    }
    
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map;
    unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>> events_map;
    queue<shared_ptr<Rec_Event>> model_queue;  // Empty queue for tests
    
    SECTION("Basic construction") {
        ModelContext model(model_params, offset_map, events_map, model_queue);
        
        // Verify const references are bound correctly
        REQUIRE(&model.model_parameters == &model_params);
        REQUIRE(&model.offset_map == &offset_map);
        REQUIRE(&model.events_map == &events_map);
        REQUIRE(&model.model_queue == &model_queue);
    }
    
    SECTION("Access model parameters") {
        ModelContext model(model_params, offset_map, events_map, model_queue);
        
        // Can read through const reference
        REQUIRE(model.model_parameters[0] == 0.0L);
        REQUIRE(model.model_parameters[50] == 0.5L);
        REQUIRE(model.model_parameters[99] == 0.99L);
    }
    
    SECTION("Empty events map") {
        ModelContext model(model_params, offset_map, events_map, model_queue);
        
        REQUIRE(model.events_map.empty());
        REQUIRE(model.offset_map.empty());
    }
}

// ============================================================================
// ScenarioContext Tests
// ============================================================================

TEST_CASE("ScenarioContext construction and mutability", "[Context][ScenarioContext]") {
    double proba = 1.0;
    Seq_type_str_p_map constructed_sequences(6);
    Seq_offsets_map seq_offsets(6, 3);
    Mismatch_vectors_map mismatches(6);
    
    SECTION("Basic construction") {
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        
        // Verify references are bound correctly
        REQUIRE(&scenario.scenario_proba == &proba);
        REQUIRE(&scenario.constructed_sequences == &constructed_sequences);
        REQUIRE(&scenario.seq_offsets == &seq_offsets);
        REQUIRE(&scenario.mismatches_lists == &mismatches);
    }
    
    SECTION("Modify probability through context") {
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        
        scenario.scenario_proba *= 0.5;
        
        REQUIRE(proba == 0.5);  // Original modified
        REQUIRE(scenario.scenario_proba == 0.5);
    }
    
    SECTION("Modify constructed sequences") {
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        
        auto v_seq = make_shared<Int_Str>(nt2int("ACGT"));
        scenario.constructed_sequences.set_value(V_gene_seq, v_seq.get(), 0);
        
        REQUIRE(scenario.constructed_sequences.at(V_gene_seq, 0) == v_seq.get());
        REQUIRE(scenario.constructed_sequences.at(V_gene_seq, 0)->size() == 4);
    }
    
    SECTION("Hot-path unified access") {
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        
        // Simulate hot-path operations: probability and sequence together
        scenario.scenario_proba *= 0.8;
        
        auto j_seq = make_shared<Int_Str>(nt2int("GC"));
        scenario.constructed_sequences.set_value(J_gene_seq, j_seq.get(), 0);
        
        scenario.seq_offsets.set_value(J_gene_seq, Five_prime, 60, 0);
        
        // All modifications visible
        REQUIRE(proba == 0.8);
        REQUIRE(scenario.constructed_sequences.at(J_gene_seq, 0)->size() == 2);
        REQUIRE(scenario.seq_offsets.at(J_gene_seq, Five_prime, 0) == 60);
    }
}

// ============================================================================
// ExplorationContext Tests
// ============================================================================

TEST_CASE("ExplorationContext construction and pruning", "[Context][ExplorationContext]") {
    Downstream_scenario_proba_bound_map proba_map(6);
    double max_prob = 1e-5;
    double threshold = 0.001;
    Index_map index_map(10);
    std::shared_ptr<Next_event_ptr> next_event_ptr (
        new Next_event_ptr[1](),  // () initializes all to nullptr
        default_delete<Next_event_ptr[]>()
    );
    Safety_bool_map safety_set(3);  // For exploration decisions
    
    SECTION("Basic construction") {
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        
        // Verify references are bound correctly
        REQUIRE(&exploration.downstream_proba_map == &proba_map);
        REQUIRE(&exploration.seq_max_prob_scenario == &max_prob);
        REQUIRE(&exploration.index_map == &index_map);
        REQUIRE(&exploration.next_event_ptr_arr == &next_event_ptr);
        REQUIRE(&exploration.safety_set == &safety_set);
        REQUIRE(exploration.proba_threshold_factor == threshold);
    }
    
    SECTION("should_prune method - aggressive pruning") {
        double threshold_factor = 0.001;
        ExplorationContext exploration(proba_map, max_prob, threshold_factor, index_map, next_event_ptr, safety_set);
        
        // max_prob = 1e-5, threshold = 0.001
        // Prune if scenario_prob < 1e-5 * 0.001 = 1e-8
        REQUIRE(exploration.should_prune(1e-9) == true);   // Below threshold
        REQUIRE(exploration.should_prune(1e-7) == false);  // Above threshold
        REQUIRE(exploration.should_prune(1e-8) == false);  // At threshold
    }
    
    SECTION("should_prune method - conservative pruning") {
        double threshold_factor = 1e-6;
        ExplorationContext exploration(proba_map, max_prob, threshold_factor, index_map, next_event_ptr, safety_set);
        
        // max_prob = 1e-5, threshold = 1e-6
        // Prune if scenario_prob < 1e-5 * 1e-6 = 1e-11
        REQUIRE(exploration.should_prune(1e-12) == true);   // Below threshold
        REQUIRE(exploration.should_prune(1e-10) == false);  // Above threshold
        REQUIRE(exploration.should_prune(1e-9) == false);   // Well above
    }
    
    SECTION("update_max_prob method") {
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        
        REQUIRE(max_prob == 1e-5);
        
        // Update with lower probability - no change
        exploration.update_max_prob(1e-6);
        REQUIRE(max_prob == 1e-5);
        REQUIRE(exploration.seq_max_prob_scenario == 1e-5);
        
        // Update with higher probability - changes
        exploration.update_max_prob(1e-4);
        REQUIRE(max_prob == 1e-4);
        REQUIRE(exploration.seq_max_prob_scenario == 1e-4);
        
        // Update again
        exploration.update_max_prob(1e-3);
        REQUIRE(max_prob == 1e-3);
    }
    
    SECTION("index_map tracking") {
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        
        // Track parent realizations
        exploration.index_map.set_value(0, 5, 0);  // Event 0, realization 5, layer 0
        exploration.index_map.set_value(1, 3, 0);  // Event 1, realization 3, layer 0
        
        REQUIRE(exploration.index_map.at(0, 0) == 5);
        REQUIRE(exploration.index_map.at(1, 0) == 3);
    }
    
    SECTION("Adaptive pruning with threshold updates") {
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        
        // Initially, max_prob = 1e-5, threshold = 0.001
        // Cutoff = 1e-8
        REQUIRE(exploration.should_prune(1e-7) == false);
        
        // Update max probability
        exploration.update_max_prob(1e-4);
        
        // Now cutoff = 1e-4 * 0.001 = 1e-7
        REQUIRE(exploration.should_prune(1e-7 + 1e-8) == false);  // At boundary
        REQUIRE(exploration.should_prune(1e-8) == true);   // Now pruned
    }
}

// ============================================================================
// AccumulationContext Tests
// ============================================================================

TEST_CASE("AccumulationContext construction", "[Context][AccumulationContext]") {
    auto marginals = make_unique<long double[]>(100);
    map<size_t, shared_ptr<Counter>> counters;
    shared_ptr<Error_rate> error_rate = make_shared<Single_error_rate>();
    
    SECTION("Basic construction") {
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        // Verify references are bound correctly
        REQUIRE(&accumulation.updated_marginals == &marginals);
        REQUIRE(&accumulation.counters == &counters);
        REQUIRE(&accumulation.error_rate == &error_rate);
    }
    
    SECTION("Access marginals") {
        // Initialize some marginal values
        marginals[0] = 0.1L;
        marginals[50] = 0.5L;
        
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        // Can access through context
        REQUIRE(accumulation.updated_marginals[0] == 0.1L);
        REQUIRE(accumulation.updated_marginals[50] == 0.5L);
        
        // Can modify through context
        accumulation.updated_marginals[0] += 0.2L;
        
        REQUIRE(marginals[0] == 0.3L);  // Original modified
        REQUIRE(accumulation.updated_marginals[0] == 0.3L);
    }
    
    SECTION("Empty counters map") {
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        REQUIRE(accumulation.counters.empty());
    }
    
    SECTION("Error rate reference") {
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        REQUIRE(accumulation.error_rate != nullptr);
        REQUIRE(accumulation.error_rate.get() == error_rate.get());
    }
}

// ============================================================================
// Cross-Context Integration Tests
// ============================================================================

TEST_CASE("Multiple contexts work together", "[Context][Integration]") {
    // Setup query
    string seq = "ACGTACGTACGT";
    Int_Str int_seq = nt2int(seq);
    unordered_map<Gene_class, vector<Alignment_data>> alignments;
    
    // Setup model
    auto model_params = make_unique<long double[]>(100);
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map;
    unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>> events_map;
    queue<shared_ptr<Rec_Event>> model_queue;  // Empty queue for tests
    
    // Setup scenario
    double proba = 1.0;
    Seq_type_str_p_map constructed_sequences(6);
    Seq_offsets_map seq_offsets(6, 3);
    Mismatch_vectors_map mismatches(6);
    
    // Setup exploration
    Downstream_scenario_proba_bound_map proba_map(6);
    double max_prob = 1e-5;
    double threshold = 0.001;
    Index_map index_map(10);
    std::shared_ptr<Next_event_ptr> next_event_ptr (
        new Next_event_ptr[1](),  // () initializes all to nullptr
        default_delete<Next_event_ptr[]>()
    );
    Safety_bool_map safety_set(3);  // For exploration decisions
    
    // Setup accumulation
    auto marginals = make_unique<long double[]>(100);
    map<size_t, shared_ptr<Counter>> counters;
    shared_ptr<Error_rate> error_rate = make_shared<Single_error_rate>();
    
    SECTION("Create all five contexts") {
        QuerySequenceContext query(seq, int_seq, alignments);
        ModelContext model(model_params, offset_map, events_map, model_queue);
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        // All contexts exist and are independent
        REQUIRE(&query.sequence == &seq);
        REQUIRE(&model.model_parameters == &model_params);
        REQUIRE(&scenario.scenario_proba == &proba);
        REQUIRE(&exploration.seq_max_prob_scenario == &max_prob);
        REQUIRE(&exploration.next_event_ptr_arr == &next_event_ptr);
        REQUIRE(&exploration.safety_set == &safety_set);  // safety_set in ExplorationContext
        REQUIRE(&accumulation.updated_marginals == &marginals);
    }
    
    SECTION("Contexts don't interfere with each other") {
        QuerySequenceContext query(seq, int_seq, alignments);
        ModelContext model(model_params, offset_map, events_map, model_queue);
        ScenarioContext scenario(proba, constructed_sequences, seq_offsets, mismatches);
        ExplorationContext exploration(proba_map, max_prob, threshold, index_map, next_event_ptr, safety_set);
        AccumulationContext accumulation(marginals, counters, error_rate);
        
        // Modify scenario
        scenario.scenario_proba *= 0.5;
        
        // Modify exploration
        exploration.update_max_prob(1e-4);
        
        // Accumulation marginals
        accumulation.updated_marginals[0] = 0.1L;
        
        // Verify changes are independent
        REQUIRE(proba == 0.5);
        REQUIRE(max_prob == 1e-4);
        REQUIRE(marginals[0] == 0.1L);
        
        // Query and model unchanged (const)
        REQUIRE(query.sequence == seq);
        REQUIRE(&model.model_parameters == &model_params);
    }
}
