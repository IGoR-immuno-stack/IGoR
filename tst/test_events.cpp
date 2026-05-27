/*
 * test_events.cpp
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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "test_utils.h"
#include <igor/Core/GenModel.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Dinuclmarkov.h>

using namespace IgorTestUtils;

TEST_CASE("Gene_choice basic functionality", "[events][gene_choice]") {
    // Load simple model with D gene
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    SECTION("Model loads successfully") {
        REQUIRE(model_parms != nullptr);
        REQUIRE(model_marginals != nullptr);
        REQUIRE(model_marginals->get_length() > 0);
    }
    
    SECTION("Gene choice events exist") {
        auto v_event = model_parms->get_event_pointer("v_choice", true);
        auto d_event = model_parms->get_event_pointer("d_gene", true);
        auto j_event = model_parms->get_event_pointer("j_choice", true);
        
        REQUIRE(v_event != nullptr);
        REQUIRE(d_event != nullptr);
        REQUIRE(j_event != nullptr);
        
        // Verify event types
        REQUIRE(v_event->get_type() == GeneChoice_t);
        REQUIRE(d_event->get_type() == GeneChoice_t);
        REQUIRE(j_event->get_type() == GeneChoice_t);
    }
}

TEST_CASE("Gene_choice generation with single realization", "[events][gene_choice][generation]") {
    // Load simple model
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    // Get gene choice events
    auto v_event = std::dynamic_pointer_cast<Gene_choice>(
        model_parms->get_event_pointer("v_choice", true));
    auto d_event = std::dynamic_pointer_cast<Gene_choice>(
        model_parms->get_event_pointer("d_gene", true));
    auto j_event = std::dynamic_pointer_cast<Gene_choice>(
        model_parms->get_event_pointer("j_choice", true));
    
    REQUIRE(v_event != nullptr);
    REQUIRE(d_event != nullptr);
    REQUIRE(j_event != nullptr);
    
    SECTION("Draw random realization from marginals") {
        // Get model queue for proper ordering
        auto model_queue = model_parms->get_model_queue();
        auto index_map = model_marginals->get_index_map(*model_parms);
        auto offset_map = model_marginals->get_offsets_map(*model_parms);
        
        // Map to store generated sequences
        std::unordered_map<Seq_type, std::string> constructed_sequences;
        
        // Random number generator
        std::mt19937_64 rng(12345);
        
        // Draw V gene realization
        auto v_realizations = v_event->draw_random_realization(
            model_marginals->marginal_array_smart_p,
            index_map,
            offset_map,
            constructed_sequences,
            rng
        );
        
        // With single realization, should be deterministic
        REQUIRE(!v_realizations.empty());
        
        // Draw D gene realization
        auto d_realizations = d_event->draw_random_realization(
            model_marginals->marginal_array_smart_p,
            index_map,
            offset_map,
            constructed_sequences,
            rng
        );
        
        REQUIRE(!d_realizations.empty());
        
        // Draw J gene realization
        auto j_realizations = j_event->draw_random_realization(
            model_marginals->marginal_array_smart_p,
            index_map,
            offset_map,
            constructed_sequences,
            rng
        );
        
        REQUIRE(!j_realizations.empty());
        
        // Verify that sequences were constructed
        REQUIRE(constructed_sequences.find(V_gene_seq) != constructed_sequences.end());
        REQUIRE(constructed_sequences.find(D_gene_seq) != constructed_sequences.end());
        REQUIRE(constructed_sequences.find(J_gene_seq) != constructed_sequences.end());
        
        // Verify sequences are not empty
        REQUIRE(!constructed_sequences[V_gene_seq].empty());
        REQUIRE(!constructed_sequences[D_gene_seq].empty());
        REQUIRE(!constructed_sequences[J_gene_seq].empty());
    }
}

TEST_CASE("Gene_choice with marginal extraction", "[events][gene_choice][marginals]") {
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    SECTION("Extract V gene marginals") {
        auto v_marginals = extract_event_marginals(
            *model_marginals, "v_choice", *model_parms);
        
        // Should have at least one realization
        REQUIRE(!v_marginals.empty());
        
        // With single realization and uniform initialization, 
        // probability should be 1.0
        if (v_marginals.size() == 1) {
            REQUIRE_THAT(v_marginals[0], 
                Catch::Matchers::WithinAbs(1.0, 1e-6));
        }
        
        // All probabilities should be non-negative
        for (double prob : v_marginals) {
            REQUIRE(prob >= 0.0);
        }
        
        // Marginals should sum to approximately 1.0
        double sum = 0.0;
        for (double prob : v_marginals) {
            sum += prob;
        }
        REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-6));
    }
    
    SECTION("Extract D gene marginals") {
        auto d_marginals = extract_event_marginals(
            *model_marginals, "d_gene", *model_parms);
        
        REQUIRE(!d_marginals.empty());
        
        // Verify normalization
        REQUIRE(is_normalized(d_marginals, 1e-6));
    }
    
    SECTION("Extract J gene marginals") {
        auto j_marginals = extract_event_marginals(
            *model_marginals, "j_choice", *model_parms);
        
        REQUIRE(!j_marginals.empty());
        
        // Verify normalization
        REQUIRE(is_normalized(j_marginals, 1e-6));
    }
}

TEST_CASE("Deletion event basic tests", "[events][deletion]") {
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    SECTION("Deletion events exist") {
        // V 3' deletion
        auto v_3_del = model_parms->get_event_pointer("v_3_del", true);
        REQUIRE(v_3_del != nullptr);
        REQUIRE(v_3_del->get_type() == Deletion_t);
        
        // D 5' and 3' deletions
        auto d_5_del = model_parms->get_event_pointer("d_5_del", true);
        auto d_3_del = model_parms->get_event_pointer("d_3_del", true);
        REQUIRE(d_5_del != nullptr);
        REQUIRE(d_3_del != nullptr);
        
        // J 5' deletion
        auto j_5_del = model_parms->get_event_pointer("j_5_del", true);
        REQUIRE(j_5_del != nullptr);
        REQUIRE(j_5_del->get_type() == Deletion_t);
    }
    
    SECTION("Deletion marginals exist and are normalized") {
        auto v_3_del_marginals = extract_event_marginals(
            *model_marginals, "v_3_del", *model_parms);
        
        REQUIRE(!v_3_del_marginals.empty());
        REQUIRE(is_normalized(v_3_del_marginals, 1e-6));
    }
}

TEST_CASE("Insertion event basic tests", "[events][insertion]") {
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    SECTION("Insertion events exist") {
        // VD insertion
        auto vd_ins = model_parms->get_event_pointer("vd_ins", true);
        REQUIRE(vd_ins != nullptr);
        REQUIRE(vd_ins->get_type() == Insertion_t);
        
        // DJ insertion
        auto dj_ins = model_parms->get_event_pointer("dj_ins", true);
        REQUIRE(dj_ins != nullptr);
        REQUIRE(dj_ins->get_type() == Insertion_t);
    }
    
    SECTION("Insertion marginals are normalized") {
        auto vd_ins_marginals = extract_event_marginals(
            *model_marginals, "vd_ins", *model_parms);
        
        REQUIRE(!vd_ins_marginals.empty());
        REQUIRE(is_normalized(vd_ins_marginals, 1e-6));
        
        auto dj_ins_marginals = extract_event_marginals(
            *model_marginals, "dj_ins", *model_parms);
        
        REQUIRE(!dj_ins_marginals.empty());
        REQUIRE(is_normalized(dj_ins_marginals, 1e-6));
    }
}

TEST_CASE("Model without D gene", "[events][no_d_gene]") {
    auto [model_parms, model_marginals] = load_simple_model(false);
    
    SECTION("Model loads successfully") {
        REQUIRE(model_parms != nullptr);
        REQUIRE(model_marginals != nullptr);
    }
    
    SECTION("V and J events exist but not D") {
        auto v_event = model_parms->get_event_pointer("v_choice", true);
        auto j_event = model_parms->get_event_pointer("j_choice", true);
        
        REQUIRE(v_event != nullptr);
        REQUIRE(j_event != nullptr);
        
        REQUIRE_THROWS(model_parms->get_event_pointer("d_gene", true));
    }
    
    SECTION("VJ insertion exists instead of VD and DJ") {
        auto vj_ins = model_parms->get_event_pointer("vj_ins", true);
        REQUIRE(vj_ins != nullptr);
        REQUIRE(vj_ins->get_type() == Insertion_t);
    }
}

TEST_CASE("Marginal validation", "[marginals][validation]") {
    auto [model_parms, model_marginals] = load_simple_model(true);
    
    SECTION("Marginals pass validation") {
        bool is_valid = validate_marginals(*model_marginals, *model_parms);
        REQUIRE(is_valid);
    }
    
    SECTION("Marginal array is properly sized") {
        size_t array_size = model_marginals->get_length();
        REQUIRE(array_size > 0);
        
        // Array should be accessible
        const auto& marginal_array = model_marginals->marginal_array_smart_p;
        REQUIRE(marginal_array != nullptr);
        
        // All values should be non-negative (after uniform initialization)
        for (size_t i = 0; i < array_size; ++i) {
            REQUIRE(marginal_array[i] >= 0.0);
        }
    }
}
