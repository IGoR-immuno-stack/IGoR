/*
 * test_tandem_guards.cpp
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
 * The error models that refuse a model they cannot read (capability plan T7a).
 *
 * The two hypermutation models find V, D and J by their legacy names and assemble the read from
 * them, so on a tandem model they would answer without an error, and wrongly (X17 of
 * ERROR_MODEL_REWRITE_PLAN). Until T7b moves them off the enum, they refuse any ordering but V-D-J
 * and V-J, naming where the ordering leaves both. Single_error_rate reads any ordering since T7a.
 * The counters' half is in tst/igor/Inference/Legacy/test_tandem_guards.cpp.
 */

#include "test_utils.h"

#include <igor/Model/Legacy/vdj_ordering_migr.h>
#include <igor/Model/Legacy/HypermutationfullNmererrorrate.h>
#include <igor/Model/Legacy/Hypermutationglobalerrorrate.h>
#include <igor/Model/Legacy/Singleerrorrate.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace IgorTestUtils;
using Catch::Matchers::ContainsSubstring;

namespace {

SeqTypeRegistry ordered(const std::vector<std::string> &segments)
{
    SeqTypeRegistry registry;
    registry.register_legacy_seq_types();
    registry.set_ordered_types(segments);
    registry.freeze();
    return registry;
}

/// The tandem fixture's layout (tst/test_data/model_formats/tandem_d_model.json).
const SeqTypeRegistry &tandem_registry()
{
    static const SeqTypeRegistry registry = ordered({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq",
                                                     "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"});
    return registry;
}

constexpr const char *kTandemOrdering =
        "V_gene_seq, VD1_ins_seq, D1_gene_seq, D1D2_ins_seq, D2_gene_seq, DJ_ins_seq, J_gene_seq";

}  // namespace

TEST_CASE("A V-D-J or V-J reader names where an ordering leaves both shapes", "[tandem_guard][tandem_d]")
{
    using igor::model::legacy::migration::require_vdj_or_vj_ordering;

    CHECK_NOTHROW(require_vdj_or_vj_ordering(vdj_seq_type_registry(), "reader"));
    CHECK_NOTHROW(require_vdj_or_vj_ordering(vj_seq_type_registry(), "reader"));

    SECTION("a tandem D, at its first junction")
    {
        CHECK_THROWS_WITH(require_vdj_or_vj_ordering(tandem_registry(), "reader"),
                          std::string("reader reads a V-D-J or a V-J model only, by the legacy segment names; this "
                                      "model's ordering is ")
                                  + kTandemOrdering + ", where it cannot place VD1_ins_seq");
    }
    SECTION("a VDJ model with a D named otherwise, at that D")
    {
        CHECK_THROWS_WITH(require_vdj_or_vj_ordering(
                                  ordered({"V_gene_seq", "VD_ins_seq", "D1_gene_seq", "DJ_ins_seq", "J_gene_seq"}),
                                  "reader"),
                          ContainsSubstring("where it cannot place D1_gene_seq"));
    }
    SECTION("a segment past a whole V-D-J")
    {
        CHECK_THROWS_WITH(require_vdj_or_vj_ordering(ordered({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq",
                                                              "J_gene_seq", "C_gene_seq"}),
                                                     "reader"),
                          ContainsSubstring("where it cannot place C_gene_seq"));
    }
    SECTION("an ordering that stops short, at what it lacks")
    {
        CHECK_THROWS_WITH(require_vdj_or_vj_ordering(ordered({"V_gene_seq", "VJ_ins_seq"}), "reader"),
                          ContainsSubstring("ordering is V_gene_seq, VJ_ins_seq, which lacks J_gene_seq"));
    }
    SECTION("the error is an invalid_argument, as the TraversalPlan's")
    {
        CHECK_THROWS_AS(require_vdj_or_vj_ordering(tandem_registry(), "reader"), std::invalid_argument);
    }
}

TEST_CASE("The global hypermutation model refuses a tandem model", "[tandem_guard][tandem_d][error_rate]")
{
    const Hypermutation_global_errorrate model(3, VDJ_genes, VDJ_genes, 0.1);
    CHECK_NOTHROW(model.check_admissible(vdj_seq_type_registry()));
    CHECK_NOTHROW(model.check_admissible(vj_seq_type_registry()));
    CHECK_THROWS_WITH(model.check_admissible(tandem_registry()),
                      ContainsSubstring("HypermutationGlobalErrorRate reads a V-D-J or a V-J model only")
                              && ContainsSubstring("where it cannot place VD1_ins_seq"));
}

TEST_CASE("The full N-mer hypermutation model refuses a tandem model", "[tandem_guard][tandem_d][error_rate]")
{
    const Hypermutation_full_Nmer_errorrate model(3, VDJ_genes, VDJ_genes, 0.1);
    CHECK_NOTHROW(model.check_admissible(vdj_seq_type_registry()));
    CHECK_NOTHROW(model.check_admissible(vj_seq_type_registry()));
    CHECK_THROWS_WITH(model.check_admissible(tandem_registry()),
                      ContainsSubstring("HypermutationFullNmerErrorrate reads a V-D-J or a V-J model only")
                              && ContainsSubstring("where it cannot place VD1_ins_seq"));
}

TEST_CASE("Single_error_rate reads a tandem model", "[tandem_guard][tandem_d][error_rate]")
{
    CHECK_NOTHROW(Single_error_rate(0.1).check_admissible(tandem_registry()));
}
