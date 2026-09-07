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
