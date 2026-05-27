/*
 * test_gene_choice_iterate.cpp
 *
 *  Created on: Jan 22, 2026
 *      Author: IGoR Test Suite
 *
 *  Tests for Gene_choice event - focusing on testable aspects without full iterate()
 *  Testing iterate() requires complex mocking, so we focus on:
 *  - Event construction and initialization
 *  - Realizations management
 *  - Marginal extraction
 *  - draw_random_realization (generation)
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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "test_utils.h"
#include <igor/Core/Genechoice.h>
#include <igor/Core/Errorrate.h>

using namespace IgorTestUtils;
using namespace Catch::Matchers;


// =============================================================================
// TEMPORARY TESTS USED FOR TEST MOCKING ARCHITECTURE DEVELOPMENT
// =============================================================================

TEST_CASE("Gene_choice event construction and realizations", "[gene_choice][basic]") {
    SECTION("Create V gene choice with single realization") {
        Gene_choice v_choice(V_gene);
        std::string gene_seq = "ACGTACGTACGT";  // 12bp V gene
        
        REQUIRE(v_choice.add_realization("TRBV1", gene_seq));
        
        // Verify event properties
        REQUIRE(v_choice.get_class() == V_gene);
        REQUIRE(v_choice.get_type() == GeneChoice_t);
        REQUIRE(v_choice.size() == 1);
        
        // Verify realization was added
        auto realizations = v_choice.get_realizations_map();
        REQUIRE(realizations.size() == 1);
        REQUIRE(realizations.count("TRBV1") == 1);
    }
    
    SECTION("Add multiple realizations") {
        Gene_choice v_choice(V_gene);
        
        REQUIRE(v_choice.add_realization("TRBV1", "ACGTACGTACGT"));
        REQUIRE(v_choice.add_realization("TRBV2", "TGCATGCATGCA"));
        REQUIRE(v_choice.add_realization("TRBV3", "GGCCGGCCGGCC"));
        
        REQUIRE(v_choice.size() == 3);
        
        auto realizations = v_choice.get_realizations_map();
        REQUIRE(realizations.size() == 3);
        REQUIRE(realizations.count("TRBV1") == 1);
        REQUIRE(realizations.count("TRBV2") == 1);
        REQUIRE(realizations.count("TRBV3") == 1);
    }
}

TEST_CASE("Gene_choice draw_random_realization", "[gene_choice][generation]") {
    // Note: draw_random_realization is the generation part of iterate
    // It picks a gene based on marginals and returns indices
    // Testing this requires complex setup, so we focus on simpler tests
    
    SECTION("Create gene choice with realizations") {
        Gene_choice v_choice(V_gene);
        v_choice.add_realization("TRBV1", "ACGTACGTACGT");
        v_choice.add_realization("TRBV2", "TGCATGCATGCA");
        v_choice.set_event_identifier(0);
        
        // Verify realizations were added
        REQUIRE(v_choice.size() == 2);
        auto realizations = v_choice.get_realizations_map();
        REQUIRE(realizations.size() == 2);
    }
}

TEST_CASE("Alignment_data construction", "[gene_choice][alignment]") {
    SECTION("Create alignment with no mismatches") {
        auto align = create_mock_alignment_data(
            "TRBV1",
            0,      // offset
            0,      // 5' offset
            11,     // 3' offset
            {},     // no mismatches
            100.0   // score
        );
        
        REQUIRE(align.gene_name == "TRBV1");
        REQUIRE(align.offset == 0);
        REQUIRE(align.mismatches.empty());
    }
    
    SECTION("Create alignment with mismatches") {
        auto align = create_mock_alignment_data(
            "TRBV2",
            5,
            5,
            16,
            {7, 9, 12},  // mismatches at positions 7, 9, 12
            95.0
        );
        
        REQUIRE(align.gene_name == "TRBV2");
        REQUIRE(align.offset == 5);
        REQUIRE(align.mismatches.size() == 3);
        REQUIRE(align.mismatches[0] == 7);
        REQUIRE(align.mismatches[1] == 9);
        REQUIRE(align.mismatches[2] == 12);
    }
}

TEST_CASE("IterateTestState creation", "[gene_choice][test_utils]") {
    SECTION("Create state with default size") {
        std::string test_seq = "ACGTACGTACGTNNNNNN";
        auto state = create_iterate_state(test_seq);
        
        REQUIRE(state.query.sequence == test_seq);
        REQUIRE(state.query.int_sequence.size() == test_seq.size());
        REQUIRE(state.accumulation.updated_marginals != nullptr);
        REQUIRE(state.model.model_parameters != nullptr);
        REQUIRE(state.accumulation.error_rate != nullptr);
    }
    
    SECTION("Create state with custom marginal size") {
        std::string test_seq = "ACGT";
        auto state = create_iterate_state(test_seq, 500);
        
        // Verify marginal arrays are initialized to 0
        bool all_zero = true;
        for (size_t i = 0; i < 500; ++i) {
            if (state.accumulation.updated_marginals[i] != 0.0 || state.model.model_parameters[i] != 0.0) {
                all_zero = false;
                break;
            }
        }
        REQUIRE(all_zero);
    }
}

TEST_CASE("Gene_choice::iterate basic call", "[gene_choice][iterate]") {
    // Test that we can successfully call iterate() using the wrapper utility
    // No assertions yet - just verifying the call completes without errors
    
    SECTION("Single V gene with single alignment - using helper functions") {
        // Create a V gene choice event
        Gene_choice v_choice(V_gene);
        std::string gene_seq = "ACGTACGTACGT";  // 12bp V gene
        v_choice.add_realization("TRBV1", gene_seq);
        v_choice.set_priority(1);
        v_choice.set_event_identifier(0);
        
        // Create test sequence containing the V gene
        std::string test_sequence = "ACGTACGTACGTNNNNNN";  // V gene + junction
        auto state = create_iterate_state(test_sequence);
        
        // Set model marginals - P(TRBV1) = 1.0
        const_cast<long double*>(state.model.model_parameters.get())[0] = 1.0;
        
        // Create alignment using helper
        std::vector<Alignment_data> v_alignments;
        v_alignments.push_back(create_perfect_alignment("TRBV1", 12));
        
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = v_alignments;
        
        // Create events map using helper
        auto v_choice_ptr = std::make_shared<Gene_choice>(v_choice);
        auto events_map = create_events_map(v_choice_ptr);
        
        // Set base_index_map for this event
        const_cast<Index_map&>(state.exploration.index_map)[v_choice_ptr->get_event_identifier()] = 0;
        
        // Call iterate using the wrapper utility (handles initialization internally)
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_choice_ptr, state);
        
        // Verify using helper functions
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime) == 0);
        REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime) == 11);
        REQUIRE(get_constructed_sequence(state, V_gene_seq) != nullptr);
        REQUIRE(get_mismatches(state, V_gene_seq).empty());
    }
    
    SECTION("V gene with multiple alignments") {
        // Create V gene choice with multiple genes
        Gene_choice v_choice(V_gene);
        v_choice.add_realization("TRBV1", "ACGTACGTACGT");  // 12bp
        v_choice.add_realization("TRBV2", "TGCATGCATGCA");  // 12bp, different
        v_choice.set_priority(1);
        v_choice.set_event_identifier(0);
        
        std::string test_sequence = "ACGTACGTACGTNNNNNN";
        auto state = create_iterate_state(test_sequence);
        
        // Set marginals: P(TRBV1) = 0.7, P(TRBV2) = 0.3
        const_cast<long double*>(state.model.model_parameters.get())[0] = 0.7;
        const_cast<long double*>(state.model.model_parameters.get())[1] = 0.3;
        
        // Both genes align (but only TRBV1 matches perfectly)
        std::vector<Alignment_data> v_alignments;
        v_alignments.push_back(create_mock_alignment_data("TRBV1", 0, 0, 11, {}));
        v_alignments.push_back(create_mock_alignment_data("TRBV2", 0, 0, 11, {0, 1, 2})); // mismatches
        
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = v_alignments;
        
        std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> events_map;
        auto v_choice_ptr = std::make_shared<Gene_choice>(v_choice);
        events_map[std::make_tuple(GeneChoice_t, V_gene, Undefined_side)] = v_choice_ptr;
        
        const_cast<Index_map&>(state.exploration.index_map)[v_choice_ptr->get_event_identifier()] = 0;
        
        // Call iterate - should explore both alignments
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_choice_ptr, state);
    }
    
    SECTION("V gene with negative offset") {
        // Negative offset means gene starts before the read
        Gene_choice v_choice(V_gene);
        std::string full_gene = "AAAACCCGGGTTT";  // 13bp
        v_choice.add_realization("TRBV1", full_gene);
        v_choice.set_priority(1);
        v_choice.set_event_identifier(0);
        
        // Read starts at position 4 of the gene
        std::string test_sequence = "CCCGGGTTTNNNNNN";  // Last 9bp of gene visible
        auto state = create_iterate_state(test_sequence);
        
        const_cast<long double*>(state.model.model_parameters.get())[0] = 1.0;
        
        // Alignment with negative offset
        std::vector<Alignment_data> v_alignments;
        v_alignments.push_back(create_mock_alignment_data(
            "TRBV1",
            -4,    // gene starts 4bp before read
            0,     // 5' offset in read
            8,     // 3' offset in read
            {}
        ));
        
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = v_alignments;
        
        std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> events_map;
        auto v_choice_ptr = std::make_shared<Gene_choice>(v_choice);
        events_map[std::make_tuple(GeneChoice_t, V_gene, Undefined_side)] = v_choice_ptr;
        
        const_cast<Index_map&>(state.exploration.index_map)[v_choice_ptr->get_event_identifier()] = 0;
        
        // Call iterate with negative offset
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_choice_ptr, state);
    }
}

// =============================================================================
// COMPREHENSIVE TEST COVERAGE - V_GENE CODE PATHS
// =============================================================================

TEST_CASE("V_gene iterate - comprehensive code path coverage", "[gene_choice][iterate][v_gene]") {
    
    SECTION("Basic V alignment - no other genes") {
        // Simplest case: Only V gene, no D or J
        // Tests core V alignment processing without safety checks
        
        std::string gene_seq = "ACGTACGTACGT";
        std::string test_sequence = "ACGTACGTACGTNNNNNN";
        double realization_proba = .5; 
        
        auto v_event = create_stub_gene_choice(V_gene, "TRBV1", gene_seq, 0);
        v_event->fix(false);  // Ensure event is not fixed (required for marginal updates)
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("TRBV1", gene_seq.length())
        };
        
        auto state = create_iterate_state(test_sequence);
        const_cast<long double*>(state.model.model_parameters.get())[0] = realization_proba;
        
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        
        // For isolated V testing without D or J genes, create events_map with only V
        auto events_map = create_events_map(v_event);
        
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        
        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        
        call_iterate(v_event, state);
        
        // Verify V gene offsets and sequence
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime) == 0);
        REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime) == gene_seq.size() -1);
        REQUIRE(get_constructed_sequence(state, V_gene_seq) != nullptr);
        REQUIRE(get_mismatches(state, V_gene_seq).empty());
        
        // When testing a single event in isolation without downstream events,
        // scenario_proba isn't updated by iterate(). Instead, check updated_marginals
        // which contains the computed probability for this realization.
        REQUIRE(state.accumulation.updated_marginals[0] == realization_proba);
    }
    
    SECTION("VD safety check: overlap detected") {
        // Code path: lines 196-209 in v1.4.0 Genechoice.cpp
        // When D chosen and (v_3_off + v_3_max_del) >= (d_5_max_offset) => continue (skip)
        
        // TODO: Set up scenario where V and D would overlap even with max deletions
        std::string seq = "ACGTACGTACGTACGT"; // TODO: Tune sequence
        std::string v_gene_seq = "ACGTACGTACGT"; // TODO: Tune V gene
        std::string d_gene_seq = "TGCATGCA"; // TODO: Tune D gene
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene_seq, 0);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene_seq, 1);
        
        // TODO: Set v_3_max_del, d_5_max_offset such that overlap occurs
        // Need to initialize D offset first via events_map
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene_seq.length())
        };
        
        auto events_map = create_events_map(v_event, d_stub, nullptr);
        auto state = create_iterate_state(seq);
        
        // TODO: Pre-set D offset in state to trigger overlap check
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_event, state);
        
        // Assertion: iterate should skip this realization (no V offset set)
        // TODO: Verify that when overlap occurs, the realization is skipped
        // REQUIRE(state.scenario_proba == 1.0); // Unchanged
        // OR check that V offsets were not set (indicating skip)
    }
    
    SECTION("VD safety: safe with min deletions") {
        // Code path: lines 206-209 in Genechoice.cpp
        // When (v_3_off + v_3_min_del) < (d_5_min_offset) => VD_safe = true
        
        // TODO: Configure scenario with large gap between V and D
        std::string seq = "ACGTACGTNNNNNNNNNNNTGCATGCA"; // TODO: Large insertion
        std::string v_gene = "ACGTACGT"; // TODO: Tune
        std::string d_gene = "TGCATGCA"; // TODO: Tune
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene.length())
        };
        
        auto events_map = create_events_map(v_event, d_stub, nullptr);
        auto state = create_iterate_state(seq);
        
        // TODO: Set D offset far from V to ensure safety
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_event, state);
        
        // When gap is large enough, even with min deletions no overlap occurs
        REQUIRE(has_safety(state, VD_safe, 0));
        REQUIRE(is_safe(state, VD_safe, 0) == true);
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) == 0);
        REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime, 0) == v_gene.length() - 1);
    }
    
    SECTION("VD safety: unsafe in deletion range") {
        // Code path: lines 209-212 in Genechoice.cpp
        // When in deletion range => VD_safe = false (some deletions cause overlap)
        
        // TODO: Configure intermediate gap requiring careful deletion control
        std::string seq = "ACGTACGTNNNNTGCA"; // TODO: Small junction
        std::string v_gene = "ACGTACGT"; // TODO: Tune
        std::string d_gene = "TGCATGCA"; // TODO: Tune
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene.length())
        };
        
        auto events_map = create_events_map(v_event, d_stub, nullptr);
        auto state = create_iterate_state(seq);
        
        // TODO: Set D offset to create ambiguous deletion scenario
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_event, state);
        
        // When in deletion range, some deletions would cause overlap
        REQUIRE(has_safety(state, VD_safe, 0));
        REQUIRE(is_safe(state, VD_safe, 0) == false);
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) == 0);
        REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime, 0) == v_gene.length() - 1);
    }
    
    SECTION("VJ safety checks (no D gene)") {
        // Code path: lines 214-236 in Genechoice.cpp
        // Similar logic to VD but for VJ junction (when D not present)
        
        // TODO: Configure VJ model (no D gene)
        std::string seq = "ACGTACGTNNNNNNGGGGTTTT"; // TODO: VJ sequence
        std::string v_gene = "ACGTACGT"; // TODO: V gene
        std::string j_gene = "GGGGTTTT"; // TODO: J gene
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto j_stub = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 1);
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene.length())
        };
        
        auto events_map = create_events_map(v_event, nullptr, j_stub);
        auto state = create_iterate_state(seq);
        
        // TODO: Set J offset appropriately
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_event, state);
        
        // Verify VJ safety flag is set correctly (depends on gap size)
        // REQUIRE(is_safe(state, VJ_safe, 0) == true); // or false, depending on configuration
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) == 0);
        REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime, 0) == v_gene.length() - 1);
    }
    
    SECTION("Multiple alignments with different positions") {
        // Code path: lines 176-289 (iterates over all alignments)
        // Tests that iterate processes all alignment candidates
        
        std::string seq = "ACGTACGTACGTACGT"; // TODO: Tune
        std::string v_gene = "ACGTACGTACGT"; // TODO: Tune
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        
        // Multiple alignments at different positions
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene.length()),
            create_mock_alignment_data("IGHV1-1*01", 2, 2, 2 + v_gene.length() - 1, {}, 100.0),
            create_mock_alignment_data("IGHV1-1*01", 4, 4, 4 + v_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(v_event, state);
        
        // Verify that at least one alignment was processed
        REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) >= 0);
        REQUIRE(get_constructed_sequence(state, V_gene_seq, 0) != nullptr);
        // TODO: Verify all three alignments contributed to marginals
        // (Would need to track how many times iterate_wrap_up was called)
    }
    
    SECTION("Alignment with mismatches") {
        // Code path: lines 250-257 (endogenous mismatch counting)
        // Tests mismatch handling and error rate computation
        
        std::string seq = "ACGTACGTACGT"; // TODO: Tune
        std::string v_gene = "ACTTACGTACGT"; // TODO: Mismatches at specific positions
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        
        // Alignment with mismatches at specific positions
        std::vector<int> mismatch_positions = {2, 5}; // TODO: Tune positions
        auto alignments = std::vector<Alignment_data>{
            create_alignment_with_mismatches("IGHV1-1*01", 0, v_gene.length(), mismatch_positions)
        };
        
        auto state = create_iterate_state(seq);
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(v_event, state);
        
        // Verify mismatches are tracked
        REQUIRE(get_mismatches(state, V_gene_seq, 0).size() == 2);
        // TODO: Verify specific mismatch positions match expected values
        // REQUIRE(get_mismatches(state, V_gene_seq, 0)[0] == 2);
        // REQUIRE(get_mismatches(state, V_gene_seq, 0)[1] == 5);
    }
    
    SECTION("Junction length probability check") {
        // Code path: lines 240-247 (junction length probability map lookup)
        // When junction length not in map => continue (skip scenario)
        
        // TODO: Configure scenario with impossible junction length
        std::string seq = "ACGTNNNNNNNNNNNNNNNNNNNNNTGCA"; // TODO: Very long junction
        std::string v_gene = "ACGT"; // TODO: Short V
        std::string d_gene = "TGCA"; // TODO: Short D
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        auto alignments = std::vector<Alignment_data>{
            create_perfect_alignment("IGHV1-1*01", v_gene.length())
        };
        
        auto events_map = create_events_map(v_event, d_stub, nullptr);
        auto state = create_iterate_state(seq);
        
        // TODO: Configure D offset to create junction length not in probability map
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(v_event, state);
        
        // When junction length not in probability map, scenario should be skipped
        // TODO: Verify V offsets were not set or scenario_proba unchanged
        // REQUIRE(state.scenario_proba == 1.0); // Should remain at initial value
    }
    
    SECTION("Probability threshold filtering") {
        // Code path: lines 261-263 (threshold check)
        // When scenario_upper_bound_proba < (seq_max_prob_scenario * threshold_factor) => continue
        
        // TODO: Configure low-probability scenario
        std::string seq = "ACGTACGTACGT"; // TODO: Tune
        std::string v_gene = "TTTTTTTTTTT"; // TODO: Mismatched gene (low probability)
        
        auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        
        // Many mismatches => low probability
        std::vector<int> many_mismatches = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // TODO: Tune
        auto alignments = std::vector<Alignment_data>{
            create_alignment_with_mismatches("IGHV1-1*01", 0, v_gene.length(), many_mismatches)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Set seq_max_prob_scenario and threshold to trigger filtering
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[V_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(v_event, state);
        
        // When probability below threshold, scenario should be filtered
        // TODO: Verify realization was skipped due to low probability
        // REQUIRE(state.scenario_proba == 1.0); // Should remain unchanged
        // OR verify that iterate_wrap_up was not called (no offsets set)
    }
}

// =============================================================================
// COMPREHENSIVE TEST COVERAGE - D_GENE CODE PATHS
// =============================================================================

TEST_CASE("D_gene iterate - comprehensive code path coverage", "[gene_choice][iterate][d_gene]") {
    
    SECTION("Basic D alignment - no other genes") {
        // Simplest case: Only D gene, no V or J
        // Tests core D alignment processing without safety checks
        
        std::string seq = "NNNNNNNTGCANNNNNN"; // TODO: Tune sequence
        std::string d_gene = "TGCA"; // TODO: Tune D gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        int d_offset = 7; // TODO: Tune D alignment position
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHD1-1*01", d_offset, d_offset, d_offset + d_gene.length() - 1, {}, 100.0)
        };
        
        auto events_map = create_events_map(nullptr, d_event, nullptr);
        auto state = create_iterate_state(seq);
        state.model.model_parameters.get()[0] = 1.0; // P(IGHD1-1*01) = 1.0
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(d_event, state);
        
        // Basic assertions: D offsets should be set correctly
        REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) == d_offset);
        REQUIRE(get_seq_offset(state, D_gene_seq, Three_prime, 0) == d_offset + d_gene.length() - 1);
        REQUIRE(get_constructed_sequence(state, D_gene_seq, 0) != nullptr);
        // TODO: Verify mismatches if alignment has them
        // No safety checks needed - no other genes
    }
    
    SECTION("Basic D exhaustive search - no other genes") {
        // Simplest exhaustive case: Only D gene, no alignments provided
        // Tests D sliding through all positions without V/J constraints
        
        std::string seq = "NNNNNNNNNNNNNNNN"; // TODO: Tune sequence
        std::string d_gene = "TGCA"; // TODO: Tune D gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        auto alignments = std::vector<Alignment_data>{}; // Empty - triggers exhaustive search
        
        auto events_map = create_events_map(nullptr, d_event, nullptr);
        auto state = create_iterate_state(seq);
        const_cast<long double*>(state.model.model_parameters.get())[0] = 1.0; // P(IGHD1-1*01) = 1.0
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(d_event, state);
        
        // D should be positioned somewhere in the sequence
        // TODO: Verify D offset is within valid range
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) >= 0);
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) < seq.length());
        // REQUIRE(get_constructed_sequence(state, D_gene_seq, 0) != nullptr);
    }
    
    SECTION("Alignment-based with VD safety check") {
        // Code path: lines 349-381 (D gene with alignments, VD overlap check)
        
        // TODO: Configure VDJ scenario with D alignments
        std::string seq = "ACGTNNNNNNNTGCANNNNGGGG"; // TODO: Full VDJ sequence
        std::string v_gene = "ACGT"; // TODO: V gene
        std::string d_gene = "TGCA"; // TODO: D gene
        std::string j_gene = "GGGG"; // TODO: J gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        auto v_stub = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto j_stub = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        
        int d_offset = 10; // TODO: Tune D alignment position
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHD1-1*01", d_offset, d_offset, d_offset + d_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Pre-set V and J offsets in state
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) == d_offset);
        REQUIRE(get_seq_offset(state, D_gene_seq, Three_prime, 0) == d_offset + d_gene.length() - 1);
        // TODO: Verify VD_safe flag based on V offset configuration
        // REQUIRE(is_safe(state, VD_safe, 0) == true); // or false
    }
    
    SECTION("DJ safety overlap detection") {
        // Code path: lines 382-395 (DJ overlap check)
        // When (d_3_off + d_3_max_del) >= (j_5_max_offset) => continue (skip)
        
        // TODO: Configure D and J that would overlap
        std::string seq = "NNNNTGCAGGGG"; // TODO: Short DJ junction
        std::string d_gene = "TGCATGCA"; // TODO: Long D gene
        std::string j_gene = "GGGGTTTT"; // TODO: J gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        auto j_stub = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        
        int d_offset = 4; // TODO: D too close to J
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHD1-1*01", d_offset, d_offset, d_offset + d_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Set J offset to trigger overlap
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        // When D and J would overlap, realization should be skipped
        // TODO: Verify D offsets were not set (indicating skip)
        // REQUIRE(state.scenario_proba == 1.0); // Should remain unchanged
    }
    
    SECTION("Exhaustive search with V and J chosen") {
        // Code path: lines 470-587 (no_d_align=true, exhaustive D positioning)
        // When no alignments provided, slides D through all valid positions
        
        // TODO: Configure scenario without D alignments
        std::string seq = "ACGTNNNNNNNNGGGG"; // TODO: VJ sequence with unknown D position
        std::string v_gene = "ACGT"; // TODO: V gene
        std::string d_gene = "TGCA"; // TODO: D gene (will be positioned exhaustively)
        std::string j_gene = "GGGG"; // TODO: J gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        auto v_stub = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        auto j_stub = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        
        // Empty alignments => triggers exhaustive search
        auto alignments = std::vector<Alignment_data>{};
        
        auto state = create_iterate_state(seq);
        // TODO: Pre-set V and J offsets to define search range
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        // In exhaustive search, D should be positioned somewhere valid
        // TODO: Verify D offset is in expected range
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) >= 0);
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) < seq.length());
        // Mismatches should be computed for the positioned D
        // TODO: Verify mismatch vector has expected size
    }
    
    SECTION("Exhaustive search sliding D position") {
        // Code path: lines 598-746 (while loop sliding D one nucleotide at a time)
        // Tests exhaustive D positioning when only V or only J chosen
        
        // TODO: Configure VD model (no J) or DJ model (no V)
        std::string seq = "ACGTNNNNNNNNNNNN"; // TODO: V + insertion region
        std::string v_gene = "ACGT"; // TODO: V gene
        std::string d_gene = "TGCA"; // TODO: D gene to slide
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        auto v_stub = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        
        auto alignments = std::vector<Alignment_data>{}; // No alignments
        
        auto state = create_iterate_state(seq);
        // TODO: Set V offset, configure deletion bounds for sliding range
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        // D should slide through valid positions between V and sequence end
        // TODO: Verify D was positioned (at least one valid position found)
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) >= 0);
        // REQUIRE(get_constructed_sequence(state, D_gene_seq, 0) != nullptr);
        // Mismatches should be updated for each sliding position
        // TODO: Verify mismatch vector is populated
    }
    
    SECTION("Both ends truncated (max deletions)") {
        // Code path: lines 463-468 (D fully deleted case check)
        // When (d_5_off - d_5_max_del) >= (d_3_off + d_3_max_del) => D invisible
        
        // TODO: Configure scenario with large deletions on D
        std::string seq = "ACGTNNNNGGGG"; // TODO: VJ with no visible D
        std::string d_gene = "TGCA"; // TODO: Short D gene
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        int d_offset = 8; // TODO: Position with heavy deletions
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHD1-1*01", d_offset, d_offset, d_offset + d_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Set deletion bounds such that D fully truncated
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        // When D fully truncated, error rate should be 1.0 (no penalty)
        REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) == d_offset);
        // TODO: Verify downstream_proba_map[D_gene_seq] == 1.0
        // (Would need to expose downstream_proba_map in test state)
    }
    
    SECTION("Mismatch computation in exhaustive search") {
        // Code path: lines 541-553 (mismatch detection in exhaustive positioning)
        // Tests mismatch list construction when D positioned manually
        
        // TODO: Configure sequence with known mismatches to D gene
        std::string seq = "ACGTNNNNNNNNGGGG"; // TODO: Sequence with mismatches
        std::string d_gene = "TGCATGCA"; // TODO: D that won't perfectly match
        
        auto d_event = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        auto v_stub = create_stub_gene_choice(V_gene, "IGHV1-1*01", "ACGT", 0);
        auto j_stub = create_stub_gene_choice(J_gene, "IGHJ1*01", "GGGG", 2);
        
        auto alignments = std::vector<Alignment_data>{}; // Exhaustive search
        
        auto state = create_iterate_state(seq);
        // TODO: Position D to have specific mismatches
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[D_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(d_event, state);
        
        // Mismatches should be computed between D gene and sequence
        // TODO: Verify mismatch positions are correct for positioned D
        // REQUIRE(get_mismatches(state, D_gene_seq, 0).size() > 0);
        // REQUIRE(get_seq_offset(state, D_gene_seq, Five_prime, 0) >= 0);
    }
}

// =============================================================================
// COMPREHENSIVE TEST COVERAGE - J_GENE CODE PATHS  
// =============================================================================

TEST_CASE("J_gene iterate - comprehensive code path coverage", "[gene_choice][iterate][j_gene]") {
    
    SECTION("Basic J alignment - no other genes") {
        // Simplest case: Only J gene, no V or D
        // Tests core J alignment processing without safety checks
        
        std::string seq = "NNNNNNGGGGTTTT"; // TODO: Tune sequence
        std::string j_gene = "GGGGTTTT"; // TODO: Tune J gene
        
        auto j_event = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        
        int j_offset = 6; // TODO: Tune J alignment position
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHJ1*01", j_offset, j_offset, j_offset + j_gene.length() - 1, {}, 100.0)
        };
        
        auto events_map = create_events_map(nullptr, nullptr, j_event);
        auto state = create_iterate_state(seq);
        const_cast<long double*>(state.model.model_parameters.get())[0] = 1.0; // P(IGHJ1*01) = 1.0
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[J_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        const_cast<std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&>(state.model.events_map) = events_map;

        call_iterate(j_event, state);
        
        // Basic assertions: J offsets should be set correctly
        REQUIRE(get_seq_offset(state, J_gene_seq, Five_prime, 0) == j_offset);
        REQUIRE(get_seq_offset(state, J_gene_seq, Three_prime, 0) == j_offset + j_gene.length() - 1);
        REQUIRE(get_constructed_sequence(state, J_gene_seq, 0) != nullptr);
        // TODO: Verify mismatches if alignment has them
        // No safety checks needed - no other genes
    }
    
    SECTION("DJ safety check") {
        // Code path: lines 858-878 (J gene with D chosen, DJ overlap check)
        
        // TODO: Configure DJ junction scenario
        std::string seq = "NNNNTGCANNNNGGGGTTTT"; // TODO: D and J genes
        std::string d_gene = "TGCA"; // TODO: D gene
        std::string j_gene = "GGGGTTTT"; // TODO: J gene
        
        auto j_event = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        int j_offset = 12; // TODO: J alignment position
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHJ1*01", j_offset, j_offset, j_offset + j_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Pre-set D offset for safety check
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[J_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(j_event, state);
        
        REQUIRE(get_seq_offset(state, J_gene_seq, Five_prime, 0) == j_offset);
        REQUIRE(get_seq_offset(state, J_gene_seq, Three_prime, 0) == j_offset + j_gene.length() - 1);
        // TODO: Verify DJ_safe flag based on D offset configuration
        // REQUIRE(is_safe(state, DJ_safe, 0) == true); // or false
    }
    
    SECTION("VJ safety check (no D)") {
        // Code path: lines 880-900 (J with V chosen, VJ overlap check)
        
        // TODO: Configure VJ model
        std::string seq = "ACGTNNNNGGGGTTTT"; // TODO: VJ sequence
        std::string v_gene = "ACGT"; // TODO: V gene
        std::string j_gene = "GGGGTTTT"; // TODO: J gene
        
        auto j_event = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        auto v_stub = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
        
        int j_offset = 8; // TODO: J position
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHJ1*01", j_offset, j_offset, j_offset + j_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Set V offset for VJ check
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[J_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(j_event, state);
        
        // Verify VJ safety flag is set correctly
        REQUIRE(get_seq_offset(state, J_gene_seq, Five_prime, 0) == j_offset);
        // TODO: Verify VJ_safe flag based on V offset and gap configuration
        // REQUIRE(is_safe(state, VJ_safe, 0) == true); // or false
    }
    
    SECTION("Overlap causes skip") {
        // Code path: lines 910-913 or 925-928 (overlap detection => continue)
        
        // TODO: Configure J too close to D or V
        std::string seq = "NNNNTGCAGGGG"; // TODO: Short junction
        std::string d_gene = "TGCATGCA"; // TODO: Long D
        std::string j_gene = "GGGGTTTT"; // TODO: J gene
        
        auto j_event = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
        
        int j_offset = 8; // TODO: Too close to D
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHJ1*01", j_offset, j_offset, j_offset + j_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        // TODO: Set D offset to trigger overlap
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[J_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(j_event, state);
        
        // When J overlaps with D or V, realization should be skipped
        // TODO: Verify J offsets were not set (indicating skip)
        // REQUIRE(state.scenario_proba == 1.0); // Should remain unchanged
    }
    
    SECTION("Multiple J alignments") {
        // Code path: lines 903-990+ (iterate over all J alignments)
        
        // TODO: Multiple possible J alignments
        std::string seq = "NNNNGGGGGGGGGGGGTTTT"; // TODO: Ambiguous J region
        std::string j_gene = "GGGGTTTT"; // TODO: J with repetitive start
        
        auto j_event = create_stub_gene_choice(J_gene, "IGHJ1*01", j_gene, 2);
        
        // Multiple J positions possible
        auto alignments = std::vector<Alignment_data>{
            create_mock_alignment_data("IGHJ1*01", 4, 4, 4 + j_gene.length() - 1, {}, 100.0),
            create_mock_alignment_data("IGHJ1*01", 8, 8, 8 + j_gene.length() - 1, {}, 100.0),
            create_mock_alignment_data("IGHJ1*01", 12, 12, 12 + j_gene.length() - 1, {}, 100.0)
        };
        
        auto state = create_iterate_state(seq);
        std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
        allowed_realizations[J_gene] = alignments;
        const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

        call_iterate(j_event, state);
        
        // Verify at least one J alignment was processed
        REQUIRE(get_seq_offset(state, J_gene_seq, Five_prime, 0) >= 0);
        REQUIRE(get_constructed_sequence(state, J_gene_seq, 0) != nullptr);
        // TODO: Verify all three alignments contributed to marginals
    }
}

// =============================================================================
// EDGE CASES AND INTEGRATION SCENARIOS
// =============================================================================

TEST_CASE("Gene_choice iterate - empty alignment list", "[gene_choice][iterate][edge]") {
    // Edge case: No alignments provided (except for D exhaustive search)
    
    std::string seq = "ACGTACGTACGT";
    std::string v_gene = "ACGTACGT";
    
    auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
    auto alignments = std::vector<Alignment_data>{}; // Empty
    
    auto state = create_iterate_state(seq);
    std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
    allowed_realizations[V_gene] = alignments;
    const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

    call_iterate(v_event, state);
    
    // With empty alignments (non-D gene), nothing should be processed
    // TODO: Verify state remains unchanged
    // REQUIRE(state.scenario_proba == 1.0);
    // REQUIRE(state.updated_marginals[0] == 0.0);
}

TEST_CASE("Gene_choice iterate - sequence boundary conditions", "[gene_choice][iterate][edge]") {
    // Edge case: Gene extends beyond sequence boundaries
    
    // TODO: Configure gene longer than sequence
    std::string seq = "ACGT"; // TODO: Very short sequence
    std::string v_gene = "ACGTACGTACGTACGT"; // TODO: Long gene
    
    auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
    
    int offset = -8; // TODO: Heavy negative offset
    auto alignments = std::vector<Alignment_data>{
        create_mock_alignment_data("IGHV1-1*01", offset, offset, offset + v_gene.length() - 1, {}, 100.0)
    };
    
    auto state = create_iterate_state(seq);
    std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
    allowed_realizations[V_gene] = alignments;
    const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

    call_iterate(v_event, state);
    
    // With heavy negative offset, only visible part of gene should be used
    // V 5' offset should be normalized to 0
    // TODO: Verify offset normalization and substring usage
    // REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) == 0);
    // REQUIRE(get_constructed_sequence(state, V_gene_seq, 0) != nullptr);
}

TEST_CASE("Gene_choice iterate - zero-length junction", "[gene_choice][iterate][edge]") {
    // Edge case: Genes directly adjacent (junction length = 0)
    
    // TODO: V and D/J immediately adjacent
    std::string seq = "ACGTTGCA"; // TODO: No insertion
    std::string v_gene = "ACGT";
    std::string d_gene = "TGCA";
    
    auto v_event = create_stub_gene_choice(V_gene, "IGHV1-1*01", v_gene, 0);
    auto d_stub = create_stub_gene_choice(D_gene, "IGHD1-1*01", d_gene, 1);
    
    auto alignments = std::vector<Alignment_data>{
        create_perfect_alignment("IGHV1-1*01", v_gene.length())
    };
    
    auto state = create_iterate_state(seq);
    // TODO: Set D offset = v_3_off + 1 (zero junction)
    std::unordered_map<Gene_class, std::vector<Alignment_data>> allowed_realizations;
    allowed_realizations[V_gene] = alignments;
    const_cast<std::unordered_map<Gene_class, std::vector<Alignment_data>>&>(state.query.gene_alignments) = allowed_realizations;

    call_iterate(v_event, state);
    
    // With zero-length junction, V and D are directly adjacent
    REQUIRE(get_seq_offset(state, V_gene_seq, Five_prime, 0) == 0);
    REQUIRE(get_seq_offset(state, V_gene_seq, Three_prime, 0) == v_gene.length() - 1);
    // Junction length should be calculated as: d_offset - v_3_off - 1 = -1
    // TODO: Verify junction length probability map is consulted with -1
}

TEST_CASE("Gene_choice iterate - all three genes with complex junctions", "[gene_choice][iterate][integration]") {
    // Integration test: Full VDJ with realistic parameters
    
    // TODO: Complete VDJ recombination scenario
    std::string seq = "ACGTACGTNNNNNNNTGCATGCANNNNNNNNGGGGTTTTCCCC"; // TODO: Realistic BCR/TCR sequence
    std::string v_gene = "ACGTACGT"; // TODO: Real V gene
    std::string d_gene = "TGCATGCA"; // TODO: Real D gene
    std::string j_gene = "GGGGTTTTCCCC"; // TODO: Real J gene
    
    // TODO: Test each gene in sequence: V -> D -> J
    // TODO: Verify safety flags transition correctly
    // TODO: Check all offsets consistent
    // TODO: Validate mismatch handling across all genes
}
