/*
 * test_generation_state.cpp
 *
 *  Unit tests for GenerationState (generation plan, G2).
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

#include <igor/Core/Legacy/StdTypedefs.h>
#include <igor/Model/Legacy/GenerationState.h>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <string>

using namespace igor::core::legacy;
using namespace igor::model::legacy;

namespace {

SeqTypeRegistry frozen(const std::vector<Seq_type_String> &ordering)
{
    SeqTypeRegistry registry;
    registry.register_legacy_seq_types();
    registry.set_ordered_types(ordering);
    registry.freeze();
    return registry;
}

} // namespace

TEST_CASE("GenerationState: a segment no event built does not exist", "[generation]")
{
    const SeqTypeRegistry registry = frozen({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    GenerationState state(registry);
    const SeqTypeId v = registry.id("V_gene_seq");

    CHECK_FALSE(state.built(v));
    CHECK_THROWS_AS(state.read(v), std::out_of_range);
    CHECK_THROWS_AS(state.modify(v), std::out_of_range);
}

TEST_CASE("GenerationState: built and empty is not the same as not built", "[generation]")
{
    const SeqTypeRegistry registry = frozen({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    GenerationState state(registry);
    const SeqTypeId vj = registry.id("VJ_ins_seq");

    state.create(vj, "");
    CHECK(state.built(vj));
    CHECK(state.read(vj).empty());
}

TEST_CASE("GenerationState: a modifier edits the segment in place", "[generation]")
{
    const SeqTypeRegistry registry = frozen({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    GenerationState state(registry);
    const SeqTypeId v = registry.id("V_gene_seq");

    state.create(v, "AACCGG");
    state.modify(v).erase(4);
    CHECK(state.read(v) == "AACC");

    SECTION("creating again replaces the segment")
    {
        state.create(v, "TT");
        CHECK(state.read(v) == "TT");
    }
}

TEST_CASE("GenerationState: assembly follows the registry's ordering", "[generation]")
{
    SECTION("VDJ, segments created out of order")
    {
        const SeqTypeRegistry registry =
                frozen({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
        GenerationState state(registry);
        state.create(registry.id("J_gene_seq"), "JJ");
        state.create(registry.id("D_gene_seq"), "DD");
        state.create(registry.id("V_gene_seq"), "VV");
        state.create(registry.id("DJ_ins_seq"), "dj");
        state.create(registry.id("VD_ins_seq"), "vd");
        CHECK(state.assemble() == "VVvdDDdjJJ");
    }

    SECTION("a segment nobody built is skipped")
    {
        const SeqTypeRegistry registry = frozen({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
        GenerationState state(registry);
        state.create(registry.id("V_gene_seq"), "VV");
        state.create(registry.id("J_gene_seq"), "JJ");
        CHECK(state.assemble() == "VVJJ");
    }

    SECTION("a segment outside the ordering is never assembled")
    {
        //register_legacy_seq_types() registers all six legacy names in every model, so a VJ
        //registry has an id for D_gene_seq. Nothing should build it; if something did, it is
        //not part of the sequence.
        const SeqTypeRegistry registry = frozen({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
        GenerationState state(registry);
        state.create(registry.id("V_gene_seq"), "VV");
        state.create(registry.id("D_gene_seq"), "DD");
        CHECK(state.assemble() == "VV");
    }

    SECTION("a tandem-D layout")
    {
        const SeqTypeRegistry registry = frozen({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq",
                                                 "D2_gene_seq", "D2J_ins_seq", "J_gene_seq"});
        GenerationState state(registry);
        for (const auto &name : registry.get_ordered_types()) {
            state.create(registry.id(name), name.substr(0, 2));
        }
        CHECK(state.assemble() == "V_VDD1D1D2D2J_");
    }
}

TEST_CASE("GenerationState: the registry must be frozen", "[generation]")
{
    SeqTypeRegistry registry;
    registry.register_legacy_seq_types();
    CHECK_THROWS_AS(GenerationState(registry), std::logic_error);
}
