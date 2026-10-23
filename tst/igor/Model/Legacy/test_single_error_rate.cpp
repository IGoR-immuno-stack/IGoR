/*
 * test_single_error_rate.cpp
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
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * The segments Single_error_rate charges at a leaf (capability plan T7a; ERROR_MODEL_REWRITE_PLAN
 * X1 and X2).
 *
 * The leaf multiplies the scenario by (r/3)^k (1-r)^(n-k), k the mismatches and n the nucleotides
 * of the segments it reads. It read V, D and J at the ids the Seq_type enum pins. A tandem model's
 * D1 and D2 have ids of their own, so their mismatches were free (X1), while their gene choices
 * still bounded their errors below 1: the product of the bounds could fall below the leaf value,
 * and a scenario above the threshold could be pruned (X2).
 *
 * Each case presets the segments as the events leave them at a leaf: a gene with its sequence and
 * its mismatch list, a junction with its ends and, once its Markov chain has run, a sequence but
 * no list. The VDJ and VJ values are those the enum version gave.
 */

#include "test_utils.h"

#include <igor/Core/Legacy/SeqTypeRegistry.h>
#include <igor/Model/Legacy/Singleerrorrate.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "LayerClaim.h"

using namespace IgorTestUtils;

namespace {

constexpr double kRate = 0.1;
constexpr double kScenarioProba = 0.25;
const std::string kRead = "ACGTACGTACGTACGTACGTACGT"; // 24 nt

/// The tandem fixture's layout (tst/test_data/model_formats/tandem_d_model.json), over the
/// pinned legacy ids, as Model_Parms::finalize() registers it.
const SeqTypeRegistry &tandem_registry()
{
    static const SeqTypeRegistry registry = [] {
        SeqTypeRegistry built;
        built.register_legacy_seq_types();
        built.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq", "D2_gene_seq",
                                 "DJ_ins_seq", "J_gene_seq"});
        built.freeze();
        return built;
    }();
    return registry;
}

/// A gene segment at the leaf: where it sits in the read, and its mismatches.
struct Gene {
    const char *seq_type;
    Seq_Offset five_prime;
    Seq_Offset three_prime;
    std::vector<std::size_t> mismatches;

    int length() const { return three_prime - five_prime + 1; }
};

/// P(scenario) (r/3)^k (1-r)^(n-k), the factor the leaf applies to k mismatches over n nucleotides.
double leaf_value(int mismatches, int nucleotides)
{
    return kScenarioProba * std::pow(kRate / 3.0, mismatches) * std::pow(1.0 - kRate, nucleotides - mismatches);
}

/// A leaf over `genes`: each preset with its sequence and mismatch list, then the error model asked.
struct Leaf {
    IterateTestState state;
    const SeqTypeRegistry &registry;

    explicit Leaf(const SeqTypeRegistry &registry_)
        : state(create_iterate_state(kRead, 1000, 32, registry_)), registry(registry_)
    {
        state.set_error_rate(kRate);
        state.set_scenario_proba(kScenarioProba);
    }

    void preset(const Gene &gene)
    {
        state.preset_segment(registry.id(gene.seq_type), gene.five_prime, gene.three_prime,
                             kRead.substr(gene.five_prime, gene.length()), gene.mismatches);
    }

    double charge()
    {
        return state.accumulation.error_rate->compute_scenario_error_probability(state.query, state.model,
                                                                                state.scenario, state.exploration);
    }

    Error_rate &error_rate() { return *state.accumulation.error_rate; }
};

}  // namespace

TEST_CASE("Single_error_rate charges V, D and J in a VDJ scenario", "[error_rate][single]")
{
    Leaf leaf(vdj_seq_type_registry());
    leaf.preset({"V_gene_seq", 0, 7, {3}});
    leaf.preset({"D_gene_seq", 10, 13, {11}});
    leaf.preset({"J_gene_seq", 16, 23, {}});
    leaf.state.preset_junction(VD_ins_seq, 8, 2);
    leaf.state.preset_junction(DJ_ins_seq, 14, 2);

    CHECK_THAT(leaf.charge(), Catch::Matchers::WithinRel(leaf_value(2, 20), 1e-12));
    CHECK_THAT(leaf.error_rate().get_seq_mean_error_number(), Catch::Matchers::WithinRel(2.0, 1e-12));
}

TEST_CASE("Single_error_rate charges V and J in a VJ scenario", "[error_rate][single]")
{
    Leaf leaf(vj_seq_type_registry());
    leaf.preset({"V_gene_seq", 0, 9, {1, 4}});
    leaf.preset({"J_gene_seq", 14, 23, {20}});
    leaf.state.preset_junction(VJ_ins_seq, 10, 4);

    CHECK_THAT(leaf.charge(), Catch::Matchers::WithinRel(leaf_value(3, 20), 1e-12));
}

TEST_CASE("Single_error_rate leaves a junction's nucleotides uncharged", "[error_rate][single]")
{
    // The Markov chain has built the junction: it has a sequence, read off the read, and no
    // mismatch list. Only the genes count.
    Leaf leaf(vdj_seq_type_registry());
    leaf.preset({"V_gene_seq", 0, 7, {3}});
    leaf.preset({"D_gene_seq", 10, 13, {}});
    leaf.preset({"J_gene_seq", 16, 23, {}});
    leaf.state.preset_junction(DJ_ins_seq, 14, 2);
    leaf.state.preset_junction(VD_ins_seq, 8, 2);
    Int_Str junction = nt2int(kRead.substr(8, 2));
    claim_layer_zero(leaf.state.scenario.constructed_sequences, VD_ins_seq);
    leaf.state.scenario.constructed_sequences.set(VD_ins_seq, &junction, 0);

    CHECK_THAT(leaf.charge(), Catch::Matchers::WithinRel(leaf_value(1, 20), 1e-12));
}

TEST_CASE("Single_error_rate charges both Ds of a tandem scenario", "[error_rate][single][tandem_d]")
{
    const std::vector<Gene> genes = {
            {"V_gene_seq", 0, 7, {3}},
            {"D1_gene_seq", 10, 13, {11}},
            {"D2_gene_seq", 16, 19, {16, 19}},
            {"J_gene_seq", 22, 23, {}},
    };
    Leaf leaf(tandem_registry());
    for (const Gene &gene : genes) {
        leaf.preset(gene);
    }

    // 4 mismatches over 18 nucleotides: V's one, D1's one, D2's two.
    const double charged = leaf.charge();
    CHECK_THAT(charged, Catch::Matchers::WithinRel(leaf_value(4, 18), 1e-12));

    SECTION("so the M-step learns from both Ds")
    {
        leaf.error_rate().add_to_norm_counter();
        leaf.error_rate().update();
        CHECK_THAT(dynamic_cast<const Single_error_rate &>(leaf.error_rate()).get_model_rate(),
                   Catch::Matchers::WithinRel(4.0 / 18.0, 1e-12));
    }

    SECTION("and the bounds the gene choices write for each segment bound the leaf (X2)")
    {
        // A gene choice and its deletions bound the segment's error factor at its own key; the
        // pruning multiplies those bounds, so it is sound only if their product is no less than
        // what the leaf applies. With one rate, the two are equal.
        double bound_product = kScenarioProba;
        for (const Gene &gene : genes) {
            const int k = static_cast<int>(gene.mismatches.size());
            bound_product *= leaf.error_rate().get_err_rate_upper_bound(k, gene.length() - k);
        }
        CHECK(charged <= bound_product * (1.0 + 1e-12));
    }
}
