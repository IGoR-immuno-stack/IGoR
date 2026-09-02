/*
 * test_event_capabilities.cpp
 *
 *  Tests for the reduced Phase A capability queries (task A0).
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

/**
 * The four queries added by task A0 replace get_len_min() / get_len_max(), whose meaning
 * differs per subclass. These tests pin each subclass's answers, and in particular pin the
 * *sign* conventions -- the thing every consumer currently re-derives by hand.
 *
 * See docs/ITERATE_GENERIC_REWRITE_PLAN.md sections 2.1, 2.2 and 7.4.
 */

#include "test_utils.h"

#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <limits>

using namespace IgorTestUtils;

namespace {

SeqTypeId id_of(Seq_type seq_type)
{
    return static_cast<SeqTypeId>(seq_type);
}

/// The id of a seq_type this event does not target, for the "everything else is zero" checks.
constexpr SeqTypeId kOtherId = static_cast<SeqTypeId>(J_gene_seq);

} // namespace

TEST_CASE("Deletion capability queries", "[capabilities][deletion]")
{
    SECTION("A 3' deletion retreats the 3' end and shortens the segment")
    {
        // Deletions 0..16. The 3' end moves left by the deletion, so the delta is negative
        // and its *minimum* corresponds to the largest deletion.
        auto del = make_deletion(V_gene_seq, Three_prime, 0, 16, /*id=*/0);

        CHECK(del->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{-16, 0});
        CHECK(del->get_length_contribution(id_of(V_gene_seq)) == LengthContribution{-16, 0});
        CHECK(del->get_seq_construction_role(id_of(V_gene_seq)) == SeqConstructionRole::Modifies);
        CHECK(del->get_offset_role(id_of(V_gene_seq), Three_prime) == OffsetRole::Modifies);

        // The other end of the same segment is untouched.
        CHECK(del->get_offset_delta_bounds(id_of(V_gene_seq), Five_prime) == OffsetDelta{});
        CHECK(del->get_offset_role(id_of(V_gene_seq), Five_prime) == OffsetRole::None);
        // So is every other segment.
        CHECK(del->get_offset_delta_bounds(kOtherId, Three_prime) == OffsetDelta{});
        CHECK(del->get_length_contribution(kOtherId) == LengthContribution{});
        CHECK(del->get_seq_construction_role(kOtherId) == SeqConstructionRole::None);
    }

    SECTION("A 5' deletion advances the 5' end: the sign flips, the length does not")
    {
        // The asymmetry the eight comparison sites in Deletion.cpp and Genechoice.cpp each
        // spell out by hand. A 5' end moves *right* as nucleotides are removed.
        auto del = make_deletion(D_gene_seq, Five_prime, 0, 16, /*id=*/0);

        CHECK(del->get_offset_delta_bounds(id_of(D_gene_seq), Five_prime) == OffsetDelta{0, 16});
        // Removing nucleotides shortens the segment whichever end they come from, so the
        // length contribution keeps the sign the offset delta just flipped.
        CHECK(del->get_length_contribution(id_of(D_gene_seq)) == LengthContribution{-16, 0});
        CHECK(del->get_offset_role(id_of(D_gene_seq), Five_prime) == OffsetRole::Modifies);
        CHECK(del->get_offset_delta_bounds(id_of(D_gene_seq), Three_prime) == OffsetDelta{});
    }

    SECTION("Palindromic insertions move the end the other way")
    {
        // Range [-4, 16]: a negative deletion adds nucleotides, so a 3' end can travel four
        // positions to the right and a 5' end four to the left.
        auto three_prime = make_deletion(V_gene_seq, Three_prime, -4, 16, /*id=*/0);
        CHECK(three_prime->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime)
              == OffsetDelta{-16, 4});
        CHECK(three_prime->get_length_contribution(id_of(V_gene_seq)) == LengthContribution{-16, 4});

        auto five_prime = make_deletion(D_gene_seq, Five_prime, -4, 16, /*id=*/0);
        CHECK(five_prime->get_offset_delta_bounds(id_of(D_gene_seq), Five_prime)
              == OffsetDelta{-4, 16});
        CHECK(five_prime->get_length_contribution(id_of(D_gene_seq)) == LengthContribution{-16, 4});
    }

    SECTION("The bounds are right across a range of shapes")
    {
        for (const auto &[min_del, max_del] :
             std::vector<std::pair<int, int>>{{0, 1}, {0, 16}, {-4, 16}, {-4, -1}, {3, 3}}) {
            auto del = make_deletion(V_gene_seq, Three_prime, min_del, max_del, /*id=*/0);
            INFO("deletion range [" << min_del << ", " << max_del << "]");
            CHECK(del->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime)
                  == OffsetDelta{-max_del, -min_del});
        }
    }

    SECTION("The bounds survive an event built by add_realization in ascending order")
    {
        // Plan section 7.4, made reproducible. Both accumulate len_min / len_max with an
        // `if / else if`:
        //
        //     if (v > -len_min)        len_min = -v;
        //     else if (v < -len_max)   len_max = -v;   // only when v is not a new maximum
        //
        // so a strictly ascending sequence never reaches the second branch. The map-based
        // constructor escapes this only because libstdc++'s hash order for the realistic
        // key sets happens not to be ascending -- an accident, not a guarantee. Calling
        // add_realization() in order reproduces it every time.
        //
        // This is exactly why A0 re-derives from the realization set: get_len_max() is read
        // as a *minimum* deletion count by every overlap check, and a sentinel there makes
        // the reachable interval absurdly wide, marking every pair safe.
        Deletion del(V_gene_seq, Three_prime);
        del.add_realization(0);
        del.add_realization(1);
        del.add_realization(2);
        del.set_seq_type("V_gene_seq");
        del.set_seq_type_id(legacy_seq_type_registry().id("V_gene_seq"));

        // The legacy accessor is left at its sentinel.
        REQUIRE(del.get_len_max() == std::numeric_limits<int16_t>::min());
        // The capability query is not.
        CHECK(del.get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{-2, 0});
        CHECK(del.get_length_contribution(id_of(V_gene_seq)) == LengthContribution{-2, 0});
    }

    SECTION("A single-realization deletion gives a degenerate, non-empty range")
    {
        auto del = make_deletion(V_gene_seq, Three_prime, 5, 5, /*id=*/0);
        CHECK(del->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{-5, -5});
    }
}

TEST_CASE("Gene_choice capability queries", "[capabilities][gene_choice]")
{
    SECTION("Creates a segment and both its offsets, and shifts nothing")
    {
        auto gc = make_gene_choice(V_gene, {{"V1", "ACGTACGT"}, {"V2", "ACGTACGTACGT"}}, 0);

        CHECK(gc->get_seq_construction_role(id_of(V_gene_seq)) == SeqConstructionRole::Creates);
        CHECK(gc->get_offset_role(id_of(V_gene_seq), Five_prime) == OffsetRole::Creates);
        CHECK(gc->get_offset_role(id_of(V_gene_seq), Three_prime) == OffsetRole::Creates);
        // An alignment fixes the position outright, so there is no travel to report.
        CHECK(gc->get_offset_delta_bounds(id_of(V_gene_seq), Five_prime) == OffsetDelta{});
        CHECK(gc->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{});
    }

    SECTION("The length contribution spans the shortest and longest templates")
    {
        auto gc = make_gene_choice(
                D_gene, {{"D1", "TTTT"}, {"D2", "TTTTTTTTTTTT"}, {"D3", "TTTTTT"}}, 0);
        CHECK(gc->get_length_contribution(id_of(D_gene_seq)) == LengthContribution{4, 12});
        CHECK(gc->get_length_contribution(kOtherId) == LengthContribution{});
    }

    SECTION("A gene choice with no realizations contributes nothing")
    {
        auto gc = make_gene_choice(V_gene, {}, 0);
        CHECK(gc->get_length_contribution(id_of(V_gene_seq)) == LengthContribution{});
    }
}

TEST_CASE("Insertion capability queries", "[capabilities][insertion]")
{
    auto ins = std::make_shared<Insertion>(VD_ins_seq, std::make_pair(0, 8));
    ins->set_seq_type("VD_ins_seq");
    ins->set_seq_type_id(legacy_seq_type_registry().id("VD_ins_seq"));

    SECTION("Creates a placeholder segment and writes no offsets")
    {
        CHECK(ins->get_seq_construction_role(id_of(VD_ins_seq)) == SeqConstructionRole::Creates);
        // The span is derived from where the neighbours already sit -- which is exactly what
        // makes the generic B6 rule possible.
        CHECK(ins->get_offset_role(id_of(VD_ins_seq), Five_prime) == OffsetRole::None);
        CHECK(ins->get_offset_role(id_of(VD_ins_seq), Three_prime) == OffsetRole::None);
        CHECK(ins->get_offset_delta_bounds(id_of(VD_ins_seq), Five_prime) == OffsetDelta{});
    }

    SECTION("The length contribution is the insertion count range, and is non-negative")
    {
        CHECK(ins->get_length_contribution(id_of(VD_ins_seq)) == LengthContribution{0, 8});
        CHECK(ins->get_length_contribution(kOtherId) == LengthContribution{});
    }
}

TEST_CASE("Dinucl_markov capability queries", "[capabilities][dinuclmarkov]")
{
    auto dinuc = std::make_shared<Dinucl_markov>(VD_ins_seq);
    dinuc->set_seq_type("VD_ins_seq");
    dinuc->set_seq_type_id(legacy_seq_type_registry().id("VD_ins_seq"));

    SECTION("Fills an existing segment: no length, no offsets")
    {
        CHECK(dinuc->get_seq_construction_role(id_of(VD_ins_seq)) == SeqConstructionRole::Fills);
        CHECK(dinuc->get_length_contribution(id_of(VD_ins_seq)) == LengthContribution{});
        CHECK(dinuc->get_offset_role(id_of(VD_ins_seq), Five_prime) == OffsetRole::None);
        CHECK(dinuc->get_offset_delta_bounds(id_of(VD_ins_seq), Three_prime) == OffsetDelta{});
    }

    SECTION("This is the subclass whose len_min / len_max were never set")
    {
        // Rec_Event initializes them to INT16_MAX / INT16_MIN and Dinucl_markov never
        // touches them, so a consumer reading the old accessors here gets sentinels. The
        // new queries answer honestly.
        CHECK(dinuc->get_len_min() == std::numeric_limits<int16_t>::max());
        CHECK(dinuc->get_len_max() == std::numeric_limits<int16_t>::min());
        CHECK(dinuc->get_length_contribution(id_of(VD_ins_seq)) == LengthContribution{});
    }
}

TEST_CASE("Capability queries agree across the subclasses that share a segment",
          "[capabilities][integration]")
{
    // What a generic body will actually do: ask every event in the model how far one end can
    // still travel, and how long the segment between two points can be. The answers have to
    // compose without the caller knowing which subclass gave them.
    auto v_choice = make_gene_choice(V_gene, {{"V1", "ACGTACGTACGT"}}, 0);
    auto v_del = make_deletion(V_gene_seq, Three_prime, 0, 4, 1);
    auto vd_ins = std::make_shared<Insertion>(VD_ins_seq, std::make_pair(0, 6));
    vd_ins->set_seq_type("VD_ins_seq");
    vd_ins->set_seq_type_id(legacy_seq_type_registry().id("VD_ins_seq"));

    SECTION("Exactly one event reports a role for each segment and end")
    {
        const std::vector<std::shared_ptr<Rec_Event>> events{v_choice, v_del, vd_ins};

        int creates_v_three_prime = 0;
        int modifies_v_three_prime = 0;
        for (const auto &event : events) {
            const OffsetRole role = event->get_offset_role(id_of(V_gene_seq), Three_prime);
            creates_v_three_prime += (role == OffsetRole::Creates);
            modifies_v_three_prime += (role == OffsetRole::Modifies);
        }
        CHECK(creates_v_three_prime == 1);
        CHECK(modifies_v_three_prime == 1);
    }

    SECTION("Length contributions sum to the reachable segment length")
    {
        // V is 12 nt and can lose up to 4, so the constructed V segment spans 8..12.
        const LengthContribution from_choice = v_choice->get_length_contribution(id_of(V_gene_seq));
        const LengthContribution from_del = v_del->get_length_contribution(id_of(V_gene_seq));
        CHECK(from_choice.min + from_del.min == 8);
        CHECK(from_choice.max + from_del.max == 12);
    }

    SECTION("Only the deletion reports travel for V's 3' end")
    {
        CHECK(v_choice->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{});
        CHECK(vd_ins->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{});
        CHECK(v_del->get_offset_delta_bounds(id_of(V_gene_seq), Three_prime) == OffsetDelta{-4, 0});
    }
}
