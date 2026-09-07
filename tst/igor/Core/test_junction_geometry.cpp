/*
 * test_junction_geometry.cpp
 *
 *  Tests for JunctionGeometry::PendingModifierBounds (task S2).
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
 * PendingModifierBounds replaces the eight `*_min_del` / `*_max_del` scalars that Deletion
 * and Gene_choice each carry, and the two 60-line blocks that fill them (plan section 2.1).
 *
 * The load-bearing claim is therefore an *equivalence*: for every end, the interval this
 * reports must be the one the legacy scalars produce when fed through the legacy offset
 * arithmetic. `legacy_offset_delta()` below spells that arithmetic out, reading the legacy
 * accessors off the same event objects, so the comparison is against the old code's own
 * convention rather than against numbers typed into this file.
 *
 * See docs/ITERATE_GENERIC_REWRITE_PLAN.md sections 2.1 and 3.
 */

#include "test_utils.h"

#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/JunctionGeometry.h>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

using namespace IgorTestUtils;
using JunctionGeometry::PendingModifierBounds;

namespace {

SeqTypeId id_of(Seq_type seq_type)
{
    return static_cast<SeqTypeId>(seq_type);
}

/// Key an event exactly as Model_Parms::get_events_map() does, so the map under test is the
/// one the model really builds -- Dinucl_markov keyed with Undefined_side included.
void add_event(Events_map &events_map, const std::shared_ptr<Rec_Event> &event)
{
    const Seq_side map_side = (event->get_type() == Dinuclmarkov_t) ? Undefined_side : event->get_side();
    events_map.emplace(std::make_tuple(event->get_type(), event->get_seq_type(), map_side), event);
}

/**
 * The interval the legacy code computes for one end, from the legacy accessors.
 *
 * `Genechoice::initialize_event()` and `Deletion::initialize_event()` both cache
 *
 *     X_min_del = del->get_len_max();      X_max_del = del->get_len_min();
 *
 * -- so both scalars hold *negated* deletion counts -- and then reach the offset interval by
 * `off + X_max_del` .. `off + X_min_del` on a 3' end but `off - X_min_del` .. `off - X_max_del`
 * on a 5' one. That asymmetry is spelled out at each of the eight comparison sites; here it is
 * written once, as the reference the new bounds are checked against.
 *
 * The events built by make_deletion() use the range constructor, which pre-seeds both bounds,
 * so get_len_min() / get_len_max() are trustworthy here -- unlike the map-based constructor
 * the model reader uses (plan section 7.4).
 */
OffsetDelta legacy_offset_delta(const Deletion &deletion, Seq_side side)
{
    const int min_del = deletion.get_len_max();
    const int max_del = deletion.get_len_min();
    if (side == Three_prime) {
        return {max_del, min_del};
    }
    return {-min_del, -max_del};
}

SeqTypeRegistry make_registry(const std::vector<Seq_type_String> &ordering)
{
    SeqTypeRegistry registry;
    //As Model_Parms does: the six standard types keep their Seq_type enum ids whatever the
    //model's ordering, so a VJ model does not silently alias J_gene_seq onto D_gene_seq's id.
    registry.register_legacy_seq_types();
    registry.set_ordered_types(ordering);
    registry.freeze();
    return registry;
}

/// A VDJ model with a deletion on all four inner ends and an insertion in both junctions.
/// The D-5' range is deliberately palindromic so that end can travel in both directions.
struct VdjModel {
    std::shared_ptr<Gene_choice> v_choice = make_gene_choice(V_gene, {{"V1", "AAAAAAAAAA"}, {"V2", "AAAAAAAAAAAA"}}, 0);
    std::shared_ptr<Gene_choice> d_choice = make_gene_choice(D_gene, {{"D1", "CCCCC"}, {"D2", "CCCCCCCC"}}, 1);
    std::shared_ptr<Gene_choice> j_choice = make_gene_choice(J_gene, {{"J1", "GGGGGGGGG"}}, 2);
    std::shared_ptr<Deletion> del_v_3 = make_deletion(V_gene_seq, Three_prime, 0, 4, 3);
    std::shared_ptr<Deletion> del_d_5 = make_deletion(D_gene_seq, Five_prime, -2, 3, 4);
    std::shared_ptr<Deletion> del_d_3 = make_deletion(D_gene_seq, Three_prime, 0, 5, 5);
    std::shared_ptr<Deletion> del_j_5 = make_deletion(J_gene_seq, Five_prime, 0, 8, 6);
    std::shared_ptr<Insertion> ins_vd = make_insertion(VD_ins_seq, 0, 6, 7);
    std::shared_ptr<Insertion> ins_dj = make_insertion(DJ_ins_seq, 0, 7, 8);
    std::shared_ptr<Dinucl_markov> dinucl = std::make_shared<Dinucl_markov>(VD_ins_seq);

    SeqTypeRegistry registry = make_registry({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
    Events_map events_map;

    VdjModel()
    {
        dinucl->set_seq_type("VD_ins_seq");
        dinucl->set_seq_type_id(legacy_seq_type_registry().id("VD_ins_seq"));
        dinucl->update_event_name();
        for (const std::shared_ptr<Rec_Event> &event :
             std::vector<std::shared_ptr<Rec_Event>>{v_choice, d_choice, j_choice, del_v_3, del_d_5,
                                                     del_d_3, del_j_5, ins_vd, ins_dj, dinucl}) {
            add_event(events_map, event);
        }
    }

    PendingModifierBounds build(const std::unordered_set<Rec_Event_name> &processed = {}) const
    {
        PendingModifierBounds bounds;
        bounds.rebuild(registry, events_map, processed);
        return bounds;
    }
};

} // namespace

TEST_CASE("PendingModifierBounds: every end matches the legacy scalar arithmetic", "[junction_geometry][s2]")
{
    const VdjModel model;
    const PendingModifierBounds bounds = model.build();

    // The four ends a deletion bears on, each against the legacy formula for its own side.
    CHECK(bounds.offset_delta(id_of(V_gene_seq), Three_prime) == legacy_offset_delta(*model.del_v_3, Three_prime));
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == legacy_offset_delta(*model.del_d_5, Five_prime));
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Three_prime) == legacy_offset_delta(*model.del_d_3, Three_prime));
    CHECK(bounds.offset_delta(id_of(J_gene_seq), Five_prime) == legacy_offset_delta(*model.del_j_5, Five_prime));

    // The same four, as literals, so a sign flip on *both* sides of the comparison above
    // cannot pass unnoticed. A 3' end retreats (negative), a 5' end advances (positive), and
    // the palindromic D-5' range reaches two positions back the other way.
    CHECK(bounds.offset_delta(id_of(V_gene_seq), Three_prime) == OffsetDelta{-4, 0});
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == OffsetDelta{-2, 3});
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Three_prime) == OffsetDelta{-5, 0});
    CHECK(bounds.offset_delta(id_of(J_gene_seq), Five_prime) == OffsetDelta{0, 8});
}

TEST_CASE("PendingModifierBounds: an end no undecided event bears on cannot move", "[junction_geometry][s2]")
{
    const VdjModel model;

    SECTION("The outer ends of V and J have no deletion in the model")
    {
        // Indistinguishable from a pinned end by design -- both mean "will not move again".
        const PendingModifierBounds bounds = model.build();
        CHECK(bounds.offset_delta(id_of(V_gene_seq), Five_prime) == OffsetDelta{});
        CHECK(bounds.offset_delta(id_of(J_gene_seq), Three_prime) == OffsetDelta{});
    }

    SECTION("Insertion segments have no offsets of their own")
    {
        // What makes the generic B6 rule possible: an insertion's span is derived from where
        // its neighbours sit, so neither of its ends is a modifiable offset.
        const PendingModifierBounds bounds = model.build();
        for (const Seq_type junction : {VD_ins_seq, DJ_ins_seq}) {
            CHECK(bounds.offset_delta(id_of(junction), Five_prime) == OffsetDelta{});
            CHECK(bounds.offset_delta(id_of(junction), Three_prime) == OffsetDelta{});
        }
    }

    SECTION("A processed deletion no longer moves its end")
    {
        // Its realization was fixed further up the scenario tree. This is the branch the
        // legacy `processed_events.count(...) != 0 ? 0 : range` test takes.
        const PendingModifierBounds bounds = model.build({model.del_v_3->get_name()});
        CHECK(bounds.offset_delta(id_of(V_gene_seq), Three_prime) == OffsetDelta{});
        // Positive control: the other three ends are untouched by that decision.
        CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == OffsetDelta{-2, 3});
        CHECK(bounds.offset_delta(id_of(J_gene_seq), Five_prime) == OffsetDelta{0, 8});
    }

    SECTION("The event asking is not yet processed, so its own travel is included")
    {
        // Rec_Event::initialize_event() inserts the caller into processed_events *after* the
        // subclass has computed its bounds, so a deletion sees itself as still pending. The
        // legacy blocks behave the same way; pinned here because it is the kind of detail a
        // rewrite silently changes.
        const PendingModifierBounds bounds = model.build();
        REQUIRE(bounds.offset_delta(id_of(V_gene_seq), Three_prime) == OffsetDelta{-4, 0});
    }
}

TEST_CASE("PendingModifierBounds: length sums every event bearing on the segment", "[junction_geometry][s2]")
{
    const VdjModel model;
    const PendingModifierBounds bounds = model.build();

    SECTION("A gene segment: its template range, less what its deletions can remove")
    {
        // V: templates 10 and 12, one deletion of up to 4.
        CHECK(bounds.length(id_of(V_gene_seq)) == LengthContribution{6, 12});
    }

    SECTION("A D segment accumulates both of its deletions")
    {
        // Templates 5 and 8; D-5' in [-2, 3] contributes [-3, 2] (a palindrome lengthens it),
        // D-3' in [0, 5] contributes [-5, 0].
        //
        // The lower bound is 5 - 3 - 5 = -3: a *bound*, not a reachable length. The sum is
        // deliberately left unclamped -- clamping here would invent a semantics no consumer
        // has asked for, and the junction-length DP that will read this applies its own
        // floor. Pinned so that adding a clamp later is a visible decision.
        CHECK(bounds.length(id_of(D_gene_seq)) == LengthContribution{-3, 10});
    }

    SECTION("A junction: the insertion range, and nothing from the Dinucl_markov filling it")
    {
        // Dinucl_markov fills placeholders Insertion already allocated, so it adds no
        // nucleotides of its own -- the subclass whose len_min / len_max were never set.
        CHECK(bounds.length(id_of(VD_ins_seq)) == LengthContribution{0, 6});
        CHECK(bounds.length(id_of(DJ_ins_seq)) == LengthContribution{0, 7});
    }

    SECTION("A segment absent from the model contributes nothing")
    {
        // VJ_ins_seq is registered -- the six standard ids always are -- but no event in a
        // VDJ model touches it.
        CHECK(bounds.length(id_of(VJ_ins_seq)) == LengthContribution{});
    }

    SECTION("Processing the gene choice removes its template from the pending length")
    {
        const PendingModifierBounds pinned = model.build({model.v_choice->get_name()});
        CHECK(pinned.length(id_of(V_gene_seq)) == LengthContribution{-4, 0});
    }
}

TEST_CASE("PendingModifierBounds: a VJ model has no D ends at all", "[junction_geometry][s2]")
{
    auto v_choice = make_gene_choice(V_gene, {{"V1", "AAAAAAAAAA"}}, 0);
    auto j_choice = make_gene_choice(J_gene, {{"J1", "GGGGGGGGG"}}, 1);
    auto del_v_3 = make_deletion(V_gene_seq, Three_prime, 0, 4, 2);
    auto del_j_5 = make_deletion(J_gene_seq, Five_prime, 0, 8, 3);
    auto ins_vj = make_insertion(VJ_ins_seq, 0, 9, 4);

    const SeqTypeRegistry registry = make_registry({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    Events_map events_map;
    for (const std::shared_ptr<Rec_Event> &event :
         std::vector<std::shared_ptr<Rec_Event>>{v_choice, j_choice, del_v_3, del_j_5, ins_vj}) {
        add_event(events_map, event);
    }

    PendingModifierBounds bounds;
    bounds.rebuild(registry, events_map, {});

    // The ends that exist behave exactly as in the VDJ model.
    CHECK(bounds.offset_delta(id_of(V_gene_seq), Three_prime) == legacy_offset_delta(*del_v_3, Three_prime));
    CHECK(bounds.offset_delta(id_of(J_gene_seq), Five_prime) == legacy_offset_delta(*del_j_5, Five_prime));
    CHECK(bounds.length(id_of(VJ_ins_seq)) == LengthContribution{0, 9});

    // The D ends read as immobile rather than as an error -- which is what lets the same
    // generic body serve both topologies, and what the legacy `d_5_min_del = 0` fallback did.
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == OffsetDelta{});
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Three_prime) == OffsetDelta{});
    CHECK(bounds.length(id_of(D_gene_seq)) == LengthContribution{});
    // ...and the junctions a VDJ model would have are equally empty.
    CHECK(bounds.length(id_of(VD_ins_seq)) == LengthContribution{});
    CHECK(bounds.length(id_of(DJ_ins_seq)) == LengthContribution{});

    // Sized from the registry, so ids agree with every other id-keyed map in the traversal.
    CHECK(bounds.seq_type_count() == registry.total_count());
}

TEST_CASE("PendingModifierBounds: two modifiers on one end sum", "[junction_geometry][s2]")
{
    // No topology IGoR builds today puts two deletions on the same end -- Model_Parms keys
    // events by (type, seq_type, side), so a second one could not even enter the map. Tandem
    // D is the case that changes: a segment reachable by two paths accumulates the modifiers
    // of both. The key below is therefore synthetic; what is under test is this container's
    // own contract, that contributions compose by addition rather than by last-write.
    auto del_a = make_deletion(D_gene_seq, Five_prime, 0, 3, 0);
    auto del_b = make_deletion(D_gene_seq, Five_prime, 0, 5, 1);

    const SeqTypeRegistry registry = make_registry({"V_gene_seq", "D_gene_seq", "J_gene_seq"});
    Events_map events_map;
    events_map.emplace(std::make_tuple(Deletion_t, Seq_type_String("D_gene_seq"), Five_prime), del_a);
    events_map.emplace(std::make_tuple(Deletion_t, Seq_type_String("D2_synthetic"), Five_prime), del_b);
    REQUIRE(events_map.size() == 2);

    PendingModifierBounds bounds;
    bounds.rebuild(registry, events_map, {});
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == OffsetDelta{0, 8});
    CHECK(bounds.length(id_of(D_gene_seq)) == LengthContribution{-8, 0});

    // Pinning one of the two leaves the other's travel behind, rather than all or nothing.
    bounds.rebuild(registry, events_map, {del_a->get_name()});
    CHECK(bounds.offset_delta(id_of(D_gene_seq), Five_prime) == OffsetDelta{0, 5});
}

TEST_CASE("PendingModifierBounds: rebuild discards the previous accumulation", "[junction_geometry][s2]")
{
    // Rebuilt at every initialize_event(), so a stale term left behind would compound
    // silently -- and would look exactly like the `+=` this class relies on.
    const VdjModel model;
    PendingModifierBounds bounds;
    bounds.rebuild(model.registry, model.events_map, {});
    const OffsetDelta once = bounds.offset_delta(id_of(J_gene_seq), Five_prime);
    bounds.rebuild(model.registry, model.events_map, {});
    CHECK(bounds.offset_delta(id_of(J_gene_seq), Five_prime) == once);
    CHECK(bounds.length(id_of(D_gene_seq)) == LengthContribution{-3, 10});
}

TEST_CASE("PendingModifierBounds: queries outside the contract throw", "[junction_geometry][s2]")
{
    const VdjModel model;
    const PendingModifierBounds bounds = model.build();

    SECTION("Undefined_side is not a segment end")
    {
        // Dinucl_markov is keyed with Undefined_side in the events map, so the value is in
        // circulation; it is not an end anything can be anchored against.
        CHECK_THROWS_AS(bounds.offset_delta(id_of(V_gene_seq), Undefined_side), std::invalid_argument);
    }

    SECTION("An id the instance was not built for")
    {
        const auto beyond = static_cast<SeqTypeId>(model.registry.total_count());
        CHECK_THROWS_AS(bounds.offset_delta(beyond, Three_prime), std::out_of_range);
        CHECK_THROWS_AS(bounds.length(beyond), std::out_of_range);
        CHECK_THROWS_AS(bounds.offset_delta(kNoSeqType, Three_prime), std::out_of_range);
    }

    SECTION("A default-constructed instance holds no ids")
    {
        const PendingModifierBounds empty;
        CHECK(empty.seq_type_count() == 0);
        CHECK_THROWS_AS(empty.offset_delta(0, Three_prime), std::out_of_range);
    }
}

// ===========================================================================================
// S3 -- reachable() and check_overlap()
// ===========================================================================================

namespace {

using JunctionGeometry::check_overlap;
using JunctionGeometry::OffsetInterval;
using JunctionGeometry::Overlap;

/**
 * The four idioms the twelve legacy comparison sites are written in, each in its own
 * variables rather than in `lo` / `hi` -- writing them in the new vocabulary would make the
 * comparison below tautological.
 *
 * The scalars are the legacy ones: `X_min_del = get_len_max() = -min_del` and
 * `X_max_del = get_len_min() = -max_del`, i.e. both hold *negated* deletion counts. An
 * interval end is then reached by `off + X_*_del` on a 3' end and `off - X_*_del` on a 5' one.
 */
struct LegacyScalars {
    int min_del = 0;   ///< X_min_del: negated *minimum* deletion count
    int max_del = 0;   ///< X_max_del: negated *maximum* deletion count
};

/// Genechoice.cpp:265,269 / :440,444 -- the event places a 3' end that can still move, against
/// a 5' neighbour that can too.
Overlap legacy_gene_choice_own_three_prime(Seq_Offset own_off, LegacyScalars own,
                                           Seq_Offset other_off, LegacyScalars other)
{
    const Seq_Offset other_min_offset = other_off - other.min_del;
    const Seq_Offset other_max_offset = other_off - other.max_del;
    if ((own_off + own.max_del) >= other_max_offset) {
        return Overlap::Infeasible;
    }
    if ((own_off + own.min_del) < other_min_offset) {
        return Overlap::Safe;
    }
    return Overlap::Undetermined;
}

/// Genechoice.cpp:427,431 / :917,921 / :930,934 -- the event places a 5' end instead, so both
/// comparisons flip.
Overlap legacy_gene_choice_own_five_prime(Seq_Offset own_off, LegacyScalars own,
                                          Seq_Offset other_off, LegacyScalars other)
{
    const Seq_Offset other_min_offset = other_off + other.max_del;
    const Seq_Offset other_max_offset = other_off + other.min_del;
    if ((own_off - own.max_del) <= other_min_offset) {
        return Overlap::Infeasible;
    }
    if ((own_off - own.min_del) > other_max_offset) {
        return Overlap::Safe;
    }
    return Overlap::Undetermined;
}

/// Deletion.cpp:320,324 / :335,339 / :787,791 -- the deletion being drawn *is* the 3' end's
/// modifier, so its own offset is already final and enters as a bare value.
Overlap legacy_deletion_own_three_prime(Seq_Offset own_new_off, Seq_Offset other_off, LegacyScalars other)
{
    const Seq_Offset other_min_offset = other_off - other.min_del;
    const Seq_Offset other_max_offset = other_off - other.max_del;
    if (own_new_off >= other_max_offset) {
        return Overlap::Infeasible;
    }
    if (own_new_off < other_min_offset) {
        return Overlap::Safe;
    }
    return Overlap::Undetermined;
}

/// Deletion.cpp:565,569 / :1049,1053 / :1063,1067 -- the same, for a 5' end.
Overlap legacy_deletion_own_five_prime(Seq_Offset own_new_off, Seq_Offset other_off, LegacyScalars other)
{
    const Seq_Offset other_min_offset = other_off + other.max_del;
    const Seq_Offset other_max_offset = other_off + other.min_del;
    if (own_new_off <= other_min_offset) {
        return Overlap::Infeasible;
    }
    if (own_new_off > other_max_offset) {
        return Overlap::Safe;
    }
    return Overlap::Undetermined;
}

/// The legacy scalars for a deletion event, read off the event itself.
LegacyScalars scalars_of(const Deletion &deletion)
{
    return {deletion.get_len_max(), deletion.get_len_min()};
}

/// Offsets swept far enough either side of the fixed end to cross both verdict boundaries.
constexpr Seq_Offset kSweepLo = -25;
constexpr Seq_Offset kSweepHi = 25;

} // namespace

TEST_CASE("reachable(): the interval absorbs the 5'/3' sign flip", "[junction_geometry][s3]")
{
    const VdjModel model;
    const PendingModifierBounds bounds = model.build();

    // V-3' deletes 0..4: the end can only retreat.
    CHECK(bounds.reachable(id_of(V_gene_seq), Three_prime, 100) == OffsetInterval{96, 100});
    // J-5' deletes 0..8: the end can only advance.
    CHECK(bounds.reachable(id_of(J_gene_seq), Five_prime, 100) == OffsetInterval{100, 108});
    // D-5' deletes -2..3: a palindromic insertion moves it the other way, so the interval
    // spans both directions and lo/hi are *not* simply {current, current + something}.
    CHECK(bounds.reachable(id_of(D_gene_seq), Five_prime, 100) == OffsetInterval{98, 103});
    CHECK(bounds.reachable(id_of(D_gene_seq), Three_prime, 100) == OffsetInterval{95, 100});

    SECTION("An end nothing pending bears on degenerates to a point")
    {
        // The unification the generic predicate rests on: Gene_choice sees a proper interval
        // for its own end, Deletion sees a point, and both call the same check.
        const PendingModifierBounds pinned = model.build({model.del_v_3->get_name()});
        CHECK(pinned.reachable(id_of(V_gene_seq), Three_prime, 100) == OffsetInterval{100, 100});
        // An end with no modifier at all is a point for the same reason.
        CHECK(pinned.reachable(id_of(J_gene_seq), Three_prime, 100) == OffsetInterval{100, 100});
    }

    SECTION("The interval is relative, so it translates with the offset")
    {
        const OffsetInterval here = bounds.reachable(id_of(D_gene_seq), Five_prime, 100);
        const OffsetInterval there = bounds.reachable(id_of(D_gene_seq), Five_prime, 140);
        CHECK(there.lo - here.lo == 40);
        CHECK(there.hi - here.hi == 40);
    }
}

TEST_CASE("check_overlap(): the verdict boundaries", "[junction_geometry][s3]")
{
    // The constraint is left_3' + gap < right_5', strictly: two segments may not occupy the
    // same read position, and an empty segment between them is allowed. Both boundaries are
    // therefore checked at the exact transition and one position either side of it.
    const OffsetInterval left{10, 14};

    SECTION("Infeasible when even the best case overlaps")
    {
        // Best case is left.lo against right.hi. left.lo == 10, so right.hi must exceed 10.
        CHECK(check_overlap(left, OffsetInterval{9, 10}, 0) == Overlap::Infeasible);
        CHECK(check_overlap(left, OffsetInterval{9, 11}, 0) != Overlap::Infeasible);
        // Equality is an overlap, not a fit: the two ends would sit on the same position.
        CHECK(check_overlap(OffsetInterval{10, 10}, OffsetInterval{10, 10}, 0) == Overlap::Infeasible);
    }

    SECTION("Safe when even the worst case fits")
    {
        // Worst case is left.hi against right.lo. left.hi == 14, so right.lo must exceed 14.
        CHECK(check_overlap(left, OffsetInterval{15, 20}, 0) == Overlap::Safe);
        CHECK(check_overlap(left, OffsetInterval{14, 20}, 0) != Overlap::Safe);
    }

    SECTION("Undetermined in between")
    {
        CHECK(check_overlap(left, OffsetInterval{11, 20}, 0) == Overlap::Undetermined);
        CHECK(check_overlap(left, OffsetInterval{14, 20}, 0) == Overlap::Undetermined);
    }

    SECTION("The three verdicts partition the geometry")
    {
        // No pair of intervals falls through, and none satisfies two verdicts at once --
        // the callers rely on the else-branch of the Safe test meaning Undetermined.
        for (Seq_Offset lo = kSweepLo; lo <= kSweepHi; ++lo) {
            for (Seq_Offset span = 0; span != 6; ++span) {
                const OffsetInterval right{lo, lo + span};
                const Overlap verdict = check_overlap(left, right, 0);
                const bool infeasible = verdict == Overlap::Infeasible;
                const bool safe = verdict == Overlap::Safe;
                REQUIRE_FALSE((infeasible && safe));
                REQUIRE((infeasible || safe || verdict == Overlap::Undetermined));
            }
        }
    }
}

TEST_CASE("check_overlap(): the gap shifts both boundaries together", "[junction_geometry][s3]")
{
    // gap is 0 at every site today; it exists so the predicate stays correct once a tandem-D
    // ordering puts a segment with a non-zero minimum length between two checked ends.
    const OffsetInterval left{10, 14};
    const OffsetInterval right{15, 20};

    REQUIRE(check_overlap(left, right, 0) == Overlap::Safe);
    // Requiring one nucleotide in between costs exactly one position of clearance.
    CHECK(check_overlap(left, right, 1) == Overlap::Undetermined);
    // Infeasible only once the gap swallows the whole clearance: left.lo + gap >= right.hi,
    // i.e. 10 + gap >= 20.
    CHECK(check_overlap(left, right, 9) == Overlap::Undetermined);
    CHECK(check_overlap(left, right, 10) == Overlap::Infeasible);

    // Equivalently: a gap of g is the left interval translated right by g.
    for (int gap = 0; gap != 14; ++gap) {
        const OffsetInterval shifted{left.lo + gap, left.hi + gap};
        CHECK(check_overlap(left, right, gap) == check_overlap(shifted, right, 0));
    }
}

TEST_CASE("Every legacy comparison site reduces to the same predicate", "[junction_geometry][s3]")
{
    // The load-bearing test of S3. Each of the twelve sites is replayed in its own variables
    // by the lambdas above, and the generic predicate must agree at every offset -- including
    // both boundary positions, which the sweep crosses.
    const VdjModel model;
    const PendingModifierBounds pending = model.build();

    // The Deletion sites run with the drawn event already consumed, so its end is a point.
    const PendingModifierBounds v_3_drawn = model.build({model.del_v_3->get_name()});
    const PendingModifierBounds d_5_drawn = model.build({model.del_d_5->get_name()});
    const PendingModifierBounds d_3_drawn = model.build({model.del_d_3->get_name()});
    const PendingModifierBounds j_5_drawn = model.build({model.del_j_5->get_name()});

    constexpr Seq_Offset kFixed = 0;

    SECTION("Gene_choice, placing a 3' end (Genechoice.cpp:265,269 and :440,444)")
    {
        for (Seq_Offset own_off = kSweepLo; own_off <= kSweepHi; ++own_off) {
            // V-3' against D-5'
            REQUIRE(check_overlap(pending.reachable(id_of(V_gene_seq), Three_prime, own_off),
                                  pending.reachable(id_of(D_gene_seq), Five_prime, kFixed), 0)
                    == legacy_gene_choice_own_three_prime(own_off, scalars_of(*model.del_v_3),
                                                          kFixed, scalars_of(*model.del_d_5)));
            // D-3' against J-5'
            REQUIRE(check_overlap(pending.reachable(id_of(D_gene_seq), Three_prime, own_off),
                                  pending.reachable(id_of(J_gene_seq), Five_prime, kFixed), 0)
                    == legacy_gene_choice_own_three_prime(own_off, scalars_of(*model.del_d_3),
                                                          kFixed, scalars_of(*model.del_j_5)));
        }
    }

    SECTION("Gene_choice, placing a 5' end (Genechoice.cpp:427,431 and :917,921 and :930,934)")
    {
        for (Seq_Offset own_off = kSweepLo; own_off <= kSweepHi; ++own_off) {
            // D-5' against V-3'
            REQUIRE(check_overlap(pending.reachable(id_of(V_gene_seq), Three_prime, kFixed),
                                  pending.reachable(id_of(D_gene_seq), Five_prime, own_off), 0)
                    == legacy_gene_choice_own_five_prime(own_off, scalars_of(*model.del_d_5),
                                                         kFixed, scalars_of(*model.del_v_3)));
            // J-5' against V-3'
            REQUIRE(check_overlap(pending.reachable(id_of(V_gene_seq), Three_prime, kFixed),
                                  pending.reachable(id_of(J_gene_seq), Five_prime, own_off), 0)
                    == legacy_gene_choice_own_five_prime(own_off, scalars_of(*model.del_j_5),
                                                         kFixed, scalars_of(*model.del_v_3)));
            // J-5' against D-3'
            REQUIRE(check_overlap(pending.reachable(id_of(D_gene_seq), Three_prime, kFixed),
                                  pending.reachable(id_of(J_gene_seq), Five_prime, own_off), 0)
                    == legacy_gene_choice_own_five_prime(own_off, scalars_of(*model.del_j_5),
                                                         kFixed, scalars_of(*model.del_d_3)));
        }
    }

    SECTION("Deletion drawing a 3' end (Deletion.cpp:320,324 and :335,339 and :787,791)")
    {
        for (Seq_Offset new_off = kSweepLo; new_off <= kSweepHi; ++new_off) {
            // V-3' against D-5', and against J-5'
            REQUIRE(check_overlap(v_3_drawn.reachable(id_of(V_gene_seq), Three_prime, new_off),
                                  v_3_drawn.reachable(id_of(D_gene_seq), Five_prime, kFixed), 0)
                    == legacy_deletion_own_three_prime(new_off, kFixed, scalars_of(*model.del_d_5)));
            REQUIRE(check_overlap(v_3_drawn.reachable(id_of(V_gene_seq), Three_prime, new_off),
                                  v_3_drawn.reachable(id_of(J_gene_seq), Five_prime, kFixed), 0)
                    == legacy_deletion_own_three_prime(new_off, kFixed, scalars_of(*model.del_j_5)));
            // D-3' against J-5'
            REQUIRE(check_overlap(d_3_drawn.reachable(id_of(D_gene_seq), Three_prime, new_off),
                                  d_3_drawn.reachable(id_of(J_gene_seq), Five_prime, kFixed), 0)
                    == legacy_deletion_own_three_prime(new_off, kFixed, scalars_of(*model.del_j_5)));
        }
    }

    SECTION("Deletion drawing a 5' end (Deletion.cpp:565,569 and :1049,1053 and :1063,1067)")
    {
        for (Seq_Offset new_off = kSweepLo; new_off <= kSweepHi; ++new_off) {
            // D-5' against V-3'
            REQUIRE(check_overlap(d_5_drawn.reachable(id_of(V_gene_seq), Three_prime, kFixed),
                                  d_5_drawn.reachable(id_of(D_gene_seq), Five_prime, new_off), 0)
                    == legacy_deletion_own_five_prime(new_off, kFixed, scalars_of(*model.del_v_3)));
            // J-5' against V-3', and against D-3'
            REQUIRE(check_overlap(j_5_drawn.reachable(id_of(V_gene_seq), Three_prime, kFixed),
                                  j_5_drawn.reachable(id_of(J_gene_seq), Five_prime, new_off), 0)
                    == legacy_deletion_own_five_prime(new_off, kFixed, scalars_of(*model.del_v_3)));
            REQUIRE(check_overlap(j_5_drawn.reachable(id_of(D_gene_seq), Three_prime, kFixed),
                                  j_5_drawn.reachable(id_of(J_gene_seq), Five_prime, new_off), 0)
                    == legacy_deletion_own_five_prime(new_off, kFixed, scalars_of(*model.del_d_3)));
        }
    }

    SECTION("The sweep actually crosses every verdict")
    {
        // Without this the agreement above could be vacuous -- two functions returning
        // Undetermined everywhere agree perfectly.
        int infeasible = 0, safe = 0, undetermined = 0;
        for (Seq_Offset own_off = kSweepLo; own_off <= kSweepHi; ++own_off) {
            switch (check_overlap(pending.reachable(id_of(V_gene_seq), Three_prime, own_off),
                                  pending.reachable(id_of(D_gene_seq), Five_prime, kFixed), 0)) {
            case Overlap::Infeasible: ++infeasible; break;
            case Overlap::Safe: ++safe; break;
            case Overlap::Undetermined: ++undetermined; break;
            }
        }
        CHECK(infeasible > 0);
        CHECK(safe > 0);
        CHECK(undetermined > 0);
    }
}

namespace {

/**
 * An event whose capability queries answer with a reversed range.
 *
 * No real subclass can produce one -- every provider builds its pair already ordered -- which
 * is exactly why the guard needs a fake to be tested at all. Without this the guard would be
 * unreachable code that no mutation could reach either.
 */
class UnorderedDeltaEvent : public Rec_Event
{
public:
    explicit UnorderedDeltaEvent(bool break_offset) : break_offset_(break_offset)
    {
        this->name = break_offset ? "bad_offset_event" : "bad_length_event";
    }

    OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const override
    {
        return break_offset_ ? OffsetDelta{4, -4} : OffsetDelta{};
    }
    LengthContribution get_length_contribution(SeqTypeId) const override
    {
        return break_offset_ ? LengthContribution{} : LengthContribution{4, -4};
    }
    SeqConstructionRole get_seq_construction_role(SeqTypeId) const override { return SeqConstructionRole::None; }
    OffsetRole get_offset_role(SeqTypeId, Seq_side) const override { return OffsetRole::None; }

    // Everything below is inert: this event exists only to be asked the four queries above.
    void iterate(QuerySequenceContext &, const ModelContext &, ScenarioContext &, ExplorationContext &,
                 AccumulationContext &) override
    {
    }
    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const override
    {
        return {};
    }
    void write2txt(std::ofstream &) override {}
    void write2txt_legacy(std::ofstream &) override {}
    void write2txt_v2(std::ofstream &) override {}
    void initialize_event(
            std::unordered_set<Rec_Event_name> &, const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &,
            std::shared_ptr<Error_rate>, Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override
    {
    }
    void add_to_marginals(long double, Marginal_array_p &) const override {}
    std::shared_ptr<Rec_Event> copy() override { return nullptr; }
    bool has_effect_on(Seq_type) const override { return false; }
    void iterate_initialize_Len_proba(Seq_type, std::map<int, double> &, std::queue<std::shared_ptr<Rec_Event>> &,
                                      double &, const Marginal_array_p &, Index_map &, Seq_type_str_p_map &,
                                      int &) const override
    {
    }
    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &, const Marginal_array_p &,
                                    Index_map &) override
    {
    }

private:
    bool break_offset_;
};

} // namespace

TEST_CASE("An unordered range from a provider is rejected, not silently sorted", "[junction_geometry][s3]")
{
    // OffsetDelta and LengthContribution are contractually ordered intervals. reachable()
    // used to re-sort them, which made a provider bug invisible *and* made the re-sorting
    // itself untestable: no real subclass can produce an unordered pair, so nothing exercised
    // it. The guard belongs where the value enters, and it belongs as an error.
    const SeqTypeRegistry registry = make_registry({"V_gene_seq", "D_gene_seq", "J_gene_seq"});
    PendingModifierBounds bounds;

    SECTION("An unordered offset delta")
    {
        Events_map events_map;
        events_map.emplace(std::make_tuple(Deletion_t, Seq_type_String("V_gene_seq"), Three_prime),
                           std::make_shared<UnorderedDeltaEvent>(true));
        CHECK_THROWS_AS(bounds.rebuild(registry, events_map, {}), std::logic_error);
    }

    SECTION("An unordered length contribution")
    {
        Events_map events_map;
        events_map.emplace(std::make_tuple(Deletion_t, Seq_type_String("V_gene_seq"), Three_prime),
                           std::make_shared<UnorderedDeltaEvent>(false));
        CHECK_THROWS_AS(bounds.rebuild(registry, events_map, {}), std::logic_error);
    }

    SECTION("The real subclasses all honour it, on every end of every segment")
    {
        // The premise reachable() now relies on, checked against the providers themselves
        // rather than assumed. A palindromic deletion range is the case most likely to break
        // it, so the VDJ model's D-5' spans [-2, 3].
        const VdjModel model;
        CHECK_NOTHROW(model.build());
        for (const auto &[key, event] : model.events_map) {
            (void)key;
            for (SeqTypeId id = 0; id != static_cast<SeqTypeId>(model.registry.total_count()); ++id) {
                for (const Seq_side side : {Five_prime, Three_prime}) {
                    const OffsetDelta delta = event->get_offset_delta_bounds(id, side);
                    INFO(event->get_name() << " offset delta on id " << id << " side " << side);
                    REQUIRE(delta.min <= delta.max);
                }
                const LengthContribution length = event->get_length_contribution(id);
                INFO(event->get_name() << " length contribution on id " << id);
                REQUIRE(length.min <= length.max);
            }
        }
    }
}
