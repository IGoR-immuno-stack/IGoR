/*
 * test_span_profile.cpp
 *
 *  Unit tests for SpanProfile and JunctionBound -- the junction-length bound's storage and
 *  addressing after S4c. The behaviour that matters is narrow: keep the best probability per
 *  distance, answer value-or-absent in one lookup, and start out unresolved so that an event
 *  with no junction on a side reads as such rather than as "distance zero, probability zero".
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

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/SegmentSpan.h>
#include <igor/Core/SpanProfile.h>

#include <vector>

TEST_CASE("SpanProfile keeps the best probability per distance", "[unit][span_profile]")
{
    SpanProfile profile;
    REQUIRE(profile.empty());

    SECTION("a first record is stored as is")
    {
        profile.record(4, 0.25);
        REQUIRE(profile.best_for(4) == 0.25);
        REQUIRE(profile.size() == 1);
    }

    SECTION("a better probability replaces the stored one, whichever order they arrive in")
    {
        profile.record(4, 0.25);
        profile.record(4, 0.75);
        REQUIRE(profile.best_for(4) == 0.75);

        SpanProfile reversed;
        reversed.record(4, 0.75);
        reversed.record(4, 0.25);
        REQUIRE(reversed.best_for(4) == 0.75);

        //Only the best survives: the map has one entry per distance, not one per path.
        REQUIRE(profile.size() == 1);
        REQUIRE(reversed.size() == 1);
    }

    SECTION("distances are independent")
    {
        profile.record(3, 0.5);
        profile.record(4, 0.1);
        REQUIRE(profile.best_for(3) == 0.5);
        REQUIRE(profile.best_for(4) == 0.1);
    }
}

TEST_CASE("SpanProfile answers absent rather than zero", "[unit][span_profile]")
{
    SpanProfile profile;

    //The distinction the consumers turn on: a distance no scenario reaches means *discard this
    //branch*, which is not the same as a branch whose bound happens to be very small. Returning
    //0.0 for both would silently keep the first and prune it later on probability.
    REQUIRE_FALSE(profile.best_for(7).has_value());

    profile.record(7, 0.0);
    REQUIRE(profile.best_for(7).has_value());
    REQUIRE(*profile.best_for(7) == 0.0);

    profile.clear();
    REQUIRE_FALSE(profile.best_for(7).has_value());
    REQUIRE(profile.empty());
}

TEST_CASE("SpanProfile iterates in distance order", "[unit][span_profile]")
{
    //Gene_choice(D) walks two profiles to build the retained decomposition, and relies on the
    //walk covering every achievable distance exactly once.
    SpanProfile profile;
    profile.record(5, 0.1);
    profile.record(1, 0.2);
    profile.record(3, 0.3);

    std::vector<int> distances;
    std::vector<double> probas;
    for (const SpanProfile::Entry entry : profile) {
        distances.push_back(entry.distance);
        probas.push_back(entry.proba);
    }
    //1 and 5 are the extremes and 3 sits between them, so the walk has to skip the four slots
    //the dense array holds for distances 2 and 4 without ever reporting them.
    REQUIRE(distances == std::vector<int>{1, 3, 5});
    REQUIRE(probas == std::vector<double>{0.2, 0.3, 0.1});
    REQUIRE(profile.size() == 3);
}

TEST_CASE("SpanProfile grows to whichever distances arrive", "[unit][span_profile]")
{
    //The fold records leaves in no particular order, so the array has to extend at either end.
    //Nothing else pins this: every other test happens to record in a convenient order.
    SpanProfile profile;
    profile.record(10, 0.5);
    profile.record(-5, 0.4);
    profile.record(30, 0.3);

    REQUIRE(profile.best_for(10) == 0.5);
    REQUIRE(profile.best_for(-5) == 0.4);
    REQUIRE(profile.best_for(30) == 0.3);
    REQUIRE(profile.size() == 3);

    //Distances inside the extent but never recorded, and distances outside it, are alike absent.
    REQUIRE_FALSE(profile.best_for(0).has_value());
    REQUIRE_FALSE(profile.best_for(-6).has_value());
    REQUIRE_FALSE(profile.best_for(31).has_value());
    REQUIRE_FALSE(profile.best_for(-100000).has_value());
    REQUIRE_FALSE(profile.best_for(100000).has_value());

    //Growing must not move what is already stored.
    profile.record(-5, 0.9);
    REQUIRE(profile.best_for(-5) == 0.9);
    REQUIRE(profile.best_for(30) == 0.3);
    REQUIRE(profile.size() == 3);
}

TEST_CASE("JunctionBound starts unresolved", "[unit][span_profile]")
{
    JunctionBound bound;

    //An event with no junction on a side -- a V gene choice, or anything in a model where the
    //neighbour was never picked -- must read as unresolved, so iterate() skips the slot instead
    //of reading a profile that was never folded.
    REQUIRE_FALSE(bound.resolved());
    REQUIRE_FALSE(bound.folded());

    const SegmentSpan span = SegmentSpan::gap(static_cast<SeqTypeId>(V_gene_seq),
                                              static_cast<SeqTypeId>(D_gene_seq));
    bound.resolve(span, static_cast<SeqTypeId>(VD_ins_seq), 2, JunctionBound::Fold::Yes);

    REQUIRE(bound.resolved());
    REQUIRE(bound.folded());
    REQUIRE(bound.span() == span);
    REQUIRE(bound.proba_key() == static_cast<SeqTypeId>(VD_ins_seq));
    REQUIRE(bound.memory_layer() == 2);

    //A junction an event splits rather than measures is resolved but carries no folded profile.
    JunctionBound enclosing;
    enclosing.resolve(SegmentSpan::gap(static_cast<SeqTypeId>(V_gene_seq), static_cast<SeqTypeId>(J_gene_seq)),
                      static_cast<SeqTypeId>(VJ_ins_seq), 1, JunctionBound::Fold::No);
    REQUIRE(enclosing.resolved());
    REQUIRE_FALSE(enclosing.folded());
}
