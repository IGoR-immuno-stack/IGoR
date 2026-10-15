/**
 * @file test_traversal_plan.cpp
 * @brief SegmentIndex and TraversalPlan: who is placed and built when each event runs, the gaps
 * that need a downstream bound and the key each is held at, and the orders an engine refuses (T3).
 *
 * The plan reads nothing from the events but their A0 roles and R12's participation rule, so these
 * cases build real events and only change the order they run in.
 */

#include "Legacy/test_utils.h"

#include <igor/Core/SeqOffsetsMap.h>
#include <igor/Model/FlankCheck.h>
#include <igor/Model/SafetyMatrix.h>
#include <igor/Model/SegmentIndex.h>
#include <igor/Model/TraversalPlan.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using namespace igor::core;
using namespace igor::model;
using namespace igor::core::legacy;
using namespace igor::model::legacy;
using namespace IgorTestUtils;
using Catch::Matchers::ContainsSubstring;

namespace {

/// Move an event built on a legacy name onto another segment of `registry`.
template <class Event>
std::shared_ptr<Event> on(std::shared_ptr<Event> event, const SeqTypeRegistry &registry, const char *seq_type)
{
    event->set_seq_type(seq_type);
    event->set_seq_type_id(registry.id(seq_type));
    event->update_event_name();
    return event;
}

SegmentSpan gap(SeqTypeId left, SeqTypeId right) { return SegmentSpan::gap(left, right); }

bool contains(const std::vector<SeqTypeId> &ids, SeqTypeId id)
{
    return std::find(ids.begin(), ids.end(), id) != ids.end();
}

/// A VDJ model's eleven events.
struct Vdj {
    std::shared_ptr<Gene_choice> v = make_gene_choice(V_gene, {{"V1", "ACGTACGT"}}, 0);
    std::shared_ptr<Gene_choice> d = make_gene_choice(D_gene, {{"D1", "ACGT"}}, 1);
    std::shared_ptr<Gene_choice> j = make_gene_choice(J_gene, {{"J1", "ACGTACGT"}}, 2);
    std::shared_ptr<Deletion> v_3 = make_deletion(V_gene_seq, Three_prime, 0, 4, 3);
    std::shared_ptr<Deletion> d_5 = make_deletion(D_gene_seq, Five_prime, 0, 4, 4);
    std::shared_ptr<Deletion> d_3 = make_deletion(D_gene_seq, Three_prime, 0, 4, 5);
    std::shared_ptr<Deletion> j_5 = make_deletion(J_gene_seq, Five_prime, 0, 4, 6);
    std::shared_ptr<Insertion> vd_ins = make_insertion(VD_ins_seq, 0, 6, 7);
    std::shared_ptr<Dinucl_markov> vd_chain = make_dinucl_markov(VD_ins_seq, 8);
    std::shared_ptr<Insertion> dj_ins = make_insertion(DJ_ins_seq, 0, 6, 9);
    std::shared_ptr<Dinucl_markov> dj_chain = make_dinucl_markov(DJ_ins_seq, 10);

    /// The legacy priorities: V and J, D, the deletions, then each junction's insertion and chain.
    std::vector<const Rec_Event *> order() const
    {
        return {v.get(),   j.get(),      d.get(),        v_3.get(),    d_5.get(),     d_3.get(),
                j_5.get(), vd_ins.get(), vd_chain.get(), dj_ins.get(), dj_chain.get()};
    }
};

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

/// A V-D1-D2-J model's events, built on legacy names and moved onto the tandem segments.
struct Tandem {
    const SeqTypeRegistry &r = tandem_registry();
    std::shared_ptr<Gene_choice> v = make_gene_choice(V_gene, {{"V1", "ACGTACGT"}}, 0);
    std::shared_ptr<Gene_choice> d1 = on(make_gene_choice(D_gene, {{"D1", "ACGT"}}, 1), r, "D1_gene_seq");
    std::shared_ptr<Gene_choice> d2 = on(make_gene_choice(D_gene, {{"D2", "ACGT"}}, 2), r, "D2_gene_seq");
    std::shared_ptr<Gene_choice> j = make_gene_choice(J_gene, {{"J1", "ACGTACGT"}}, 3);
    std::shared_ptr<Deletion> d1_3 = on(make_deletion(D_gene_seq, Three_prime, 0, 4, 4), r, "D1_gene_seq");
    std::shared_ptr<Deletion> d2_5 = on(make_deletion(D_gene_seq, Five_prime, 0, 4, 5), r, "D2_gene_seq");
    std::shared_ptr<Insertion> vd1_ins = on(make_insertion(VD_ins_seq, 0, 6, 6), r, "VD1_ins_seq");
    std::shared_ptr<Dinucl_markov> vd1_chain = on(make_dinucl_markov(VD_ins_seq, 7), r, "VD1_ins_seq");
    std::shared_ptr<Insertion> d1d2_ins = on(make_insertion(VD_ins_seq, 0, 6, 8), r, "D1D2_ins_seq");
    std::shared_ptr<Dinucl_markov> d1d2_chain = on(make_dinucl_markov(VD_ins_seq, 9), r, "D1D2_ins_seq");
    std::shared_ptr<Insertion> dj_ins = make_insertion(DJ_ins_seq, 0, 6, 10);
    std::shared_ptr<Dinucl_markov> dj_chain = make_dinucl_markov(DJ_ins_seq, 11);

    SeqTypeId id(const char *seq_type) const { return r.id(seq_type); }

    std::vector<const Rec_Event *> order(bool d1_first = true) const
    {
        const Rec_Event *first = d1_first ? d1.get() : d2.get();
        const Rec_Event *second = d1_first ? d2.get() : d1.get();
        return {v.get(),       j.get(),         first,          second,          d1_3.get(),
                d2_5.get(),    vd1_ins.get(),   vd1_chain.get(), d1d2_ins.get(), d1d2_chain.get(),
                dj_ins.get(),  dj_chain.get()};
    }
};

} // namespace

TEST_CASE("SegmentIndex: who creates and who modifies each end and each content", "[traversal_plan]")
{
    const Vdj model;
    const std::vector<const Rec_Event *> events = model.order();
    const SegmentIndex index(vdj_seq_type_registry(), events);

    CHECK(index.end_creators(V_gene_seq, Three_prime) == SegmentIndex::Events{model.v.get()});
    CHECK(index.end_modifiers(V_gene_seq, Three_prime) == SegmentIndex::Events{model.v_3.get()});
    CHECK(index.end_modifiers(V_gene_seq, Five_prime).empty());
    CHECK(index.content_modifiers(D_gene_seq) == SegmentIndex::Events{model.d_5.get(), model.d_3.get()});
    // A junction: its ends placed by the insertion, its content built by the chain.
    CHECK(index.end_creators(VD_ins_seq, Five_prime) == SegmentIndex::Events{model.vd_ins.get()});
    CHECK(index.content_creators(VD_ins_seq) == SegmentIndex::Events{model.vd_chain.get()});

    SECTION("a factor of its own: templated or trimmed content, not content filled into a slot") {
        CHECK(index.has_own_factor(V_gene_seq));
        CHECK(index.has_own_factor(D_gene_seq));
        CHECK(index.has_own_factor(J_gene_seq));
        CHECK_FALSE(index.has_own_factor(VD_ins_seq));
        CHECK_FALSE(index.has_own_factor(DJ_ins_seq));
    }
}

TEST_CASE("TraversalPlan: what each event finds placed and built, in a VDJ order", "[traversal_plan]")
{
    const Vdj model;
    const TraversalPlan plan(vdj_seq_type_registry(), model.order());

    // V runs first: nothing to check against, nothing to bound.
    CHECK(plan.placement(*model.v).placed_left.empty());
    CHECK(plan.placement(*model.v).placed_right.empty());
    // J finds V.
    CHECK(plan.placement(*model.j).placed_left == std::vector<SeqTypeId>{V_gene_seq});
    // D finds both, and splits the gap V->J they leave.
    const EventPlacement &d = plan.placement(*model.d);
    CHECK(d.placed_left == std::vector<SeqTypeId>{V_gene_seq});
    CHECK(d.placed_right == std::vector<SeqTypeId>{J_gene_seq});
    CHECK(d.enclosing_gap() == gap(V_gene_seq, J_gene_seq));
    CHECK(d.gap_toward(Five_prime) == gap(V_gene_seq, D_gene_seq));
    CHECK(d.gap_toward(Three_prime) == gap(D_gene_seq, J_gene_seq));
    // V's 3' deletion is checked against D and J, nearest first, and widens V->D.
    const EventPlacement &v_3 = plan.placement(*model.v_3);
    CHECK(v_3.placed_right == std::vector<SeqTypeId>{D_gene_seq, J_gene_seq});
    CHECK(v_3.gap_toward(Three_prime) == gap(V_gene_seq, D_gene_seq));
    // The insertion sits in V->D, unbuilt; the chain fills the same gap.
    CHECK(plan.placement(*model.vd_ins).enclosing_gap() == gap(V_gene_seq, D_gene_seq));
    CHECK_FALSE(plan.placement(*model.vd_ins).target_built);
    CHECK(plan.placement(*model.vd_chain).enclosing_gap() == gap(V_gene_seq, D_gene_seq));

    SECTION("placed is not built: a slot its chain has not filled yet") {
        // D's 5' deletion moved between the VD insertion and its chain: it finds the slot placed,
        // so checks against it, but the gap it widens still runs from V, the nearest content built.
        auto order = model.order();
        order.erase(std::find(order.begin(), order.end(), model.d_5.get()));
        order.insert(std::find(order.begin(), order.end(), model.vd_chain.get()), model.d_5.get());
        const TraversalPlan reordered(vdj_seq_type_registry(), order);
        CHECK_NOTHROW(reordered.check_admissible());
        const EventPlacement &d_5 = reordered.placement(*model.d_5);
        CHECK(d_5.placed_left == std::vector<SeqTypeId>{VD_ins_seq, V_gene_seq});
        CHECK(d_5.built_left == V_gene_seq);
        CHECK(d_5.gap_toward(Five_prime) == gap(V_gene_seq, D_gene_seq));
    }

    SECTION("no gene choice or deletion finds an insertion's segment placed: they all run first") {
        for (const Rec_Event *event : {static_cast<const Rec_Event *>(model.v.get()), static_cast<const Rec_Event *>(model.d.get()),
                                       static_cast<const Rec_Event *>(model.j.get()), static_cast<const Rec_Event *>(model.v_3.get()),
                                       static_cast<const Rec_Event *>(model.d_5.get()), static_cast<const Rec_Event *>(model.d_3.get()),
                                       static_cast<const Rec_Event *>(model.j_5.get())}) {
            const EventPlacement &placed = plan.placement(*event);
            for (const SeqTypeId junction : {SeqTypeId{VD_ins_seq}, SeqTypeId{DJ_ins_seq}}) {
                CHECK_FALSE(contains(placed.placed_left, junction));
                CHECK_FALSE(contains(placed.placed_right, junction));
            }
        }
    }
}

TEST_CASE("TraversalPlan: a gap is held at the slot inside it, unless it encloses more", "[traversal_plan]")
{
    SECTION("VDJ: V->D and D->J at their insertions' ids, V->J at a key of its own") {
        const Vdj model;
        const SeqTypeRegistry &registry = vdj_seq_type_registry();
        const TraversalPlan plan(registry, model.order());
        CHECK(plan.gaps() == std::vector<SegmentSpan>{gap(V_gene_seq, J_gene_seq), gap(V_gene_seq, D_gene_seq),
                                                      gap(D_gene_seq, J_gene_seq)});
        CHECK(plan.gap_key(gap(V_gene_seq, D_gene_seq)) == VD_ins_seq);
        CHECK(plan.gap_key(gap(D_gene_seq, J_gene_seq)) == DJ_ins_seq);
        // VJ_ins_seq's id stays unused, so the product of the bound is unchanged.
        CHECK(plan.gap_key(gap(V_gene_seq, J_gene_seq)) == registry.total_count());
        CHECK(plan.bound_key_count() == registry.total_count() + 1);
        // Once the chain fills a junction, the two halves have nothing left to bound.
        CHECK_THROWS_AS(plan.gap_key(gap(V_gene_seq, VD_ins_seq)), std::out_of_range);
    }

    SECTION("VJ: V->J at the VJ insertion's id, nothing past the registry") {
        const SeqTypeRegistry &registry = vj_seq_type_registry();
        const auto v = make_gene_choice(V_gene, {{"V1", "ACGTACGT"}}, 0);
        const auto j = make_gene_choice(J_gene, {{"J1", "ACGTACGT"}}, 1);
        const auto v_3 = make_deletion(V_gene_seq, Three_prime, 0, 4, 2);
        const auto ins = make_insertion(VJ_ins_seq, 0, 6, 3);
        const auto chain = make_dinucl_markov(VJ_ins_seq, 4);
        const TraversalPlan plan(registry, std::vector<const Rec_Event *>{v.get(), j.get(), v_3.get(), ins.get(),
                                                                          chain.get()});
        CHECK(plan.gap_key(gap(V_gene_seq, J_gene_seq)) == VJ_ins_seq);
        CHECK(plan.bound_key_count() == registry.total_count());
        CHECK_NOTHROW(plan.check_admissible());
    }
}

TEST_CASE("TraversalPlan: tandem D", "[traversal_plan][tandem_d]")
{
    const Tandem model;
    const auto key = [&](const TraversalPlan &plan, const char *left, const char *right) {
        return plan.gap_key(gap(model.id(left), model.id(right)));
    };

    SECTION("D1 placed before D2") {
        const TraversalPlan plan(model.r, model.order());
        CHECK_NOTHROW(plan.check_admissible());

        // The gaps between neighbours are held at the insertion filling them...
        CHECK(key(plan, "V_gene_seq", "D1_gene_seq") == model.id("VD1_ins_seq"));
        CHECK(key(plan, "D1_gene_seq", "D2_gene_seq") == model.id("D1D2_ins_seq"));
        CHECK(key(plan, "D2_gene_seq", "J_gene_seq") == model.id("DJ_ins_seq"));
        // ...and the two this order opens around a D at keys of their own.
        const SeqTypeId v_to_j = key(plan, "V_gene_seq", "J_gene_seq");
        const SeqTypeId d1_to_j = key(plan, "D1_gene_seq", "J_gene_seq");
        CHECK(v_to_j >= model.r.total_count());
        CHECK(d1_to_j >= model.r.total_count());
        CHECK(v_to_j != d1_to_j);
        CHECK(plan.bound_key_count() == model.r.total_count() + 2);
        CHECK_THROWS_AS(key(plan, "V_gene_seq", "D2_gene_seq"), std::out_of_range); // never open

        // D2 is checked against D1 and V to its left, J to its right, and splits D1->J.
        const EventPlacement &d2 = plan.placement(*model.d2);
        CHECK(d2.placed_left == std::vector<SeqTypeId>{model.id("D1_gene_seq"), V_gene_seq});
        CHECK(d2.placed_right == std::vector<SeqTypeId>{J_gene_seq});
        CHECK(d2.enclosing_gap() == gap(model.id("D1_gene_seq"), J_gene_seq));
        // D1's 3' deletion is checked against D2 and J, and widens D1->D2.
        const EventPlacement &d1_3 = plan.placement(*model.d1_3);
        CHECK(d1_3.placed_right == std::vector<SeqTypeId>{model.id("D2_gene_seq"), J_gene_seq});
        CHECK(d1_3.gap_toward(Three_prime) == gap(model.id("D1_gene_seq"), model.id("D2_gene_seq")));
    }

    SECTION("D2 placed before D1: the other enclosing gap") {
        const TraversalPlan plan(model.r, model.order(/*d1_first=*/false));
        CHECK_NOTHROW(plan.check_admissible());
        CHECK(key(plan, "V_gene_seq", "D2_gene_seq") >= model.r.total_count());
        CHECK_THROWS_AS(key(plan, "D1_gene_seq", "J_gene_seq"), std::out_of_range);
    }
}

TEST_CASE("TraversalPlan: an order the events cannot run in is refused", "[traversal_plan]")
{
    const Vdj model;
    const SeqTypeRegistry &registry = vdj_seq_type_registry();
    const auto refused = [&](std::vector<const Rec_Event *> order) {
        return TraversalPlan(registry, order);
    };
    auto order = model.order();
    const auto move_before = [&](const Rec_Event *event, const Rec_Event *anchor) {
        order.erase(std::find(order.begin(), order.end(), event));
        order.insert(std::find(order.begin(), order.end(), anchor), event);
    };

    SECTION("the legacy order is admissible") { CHECK_NOTHROW(refused(order).check_admissible()); }
    SECTION("a deletion before the gene choice it trims") {
        move_before(model.v_3.get(), model.v.get());
        CHECK_THROWS_WITH(refused(order).check_admissible(),
                          ContainsSubstring("modifies the 3' end of V_gene_seq before any event creates it"));
    }
    SECTION("a chain before the insertion that places its slot") {
        move_before(model.vd_chain.get(), model.vd_ins.get());
        CHECK_THROWS_WITH(refused(order).check_admissible(),
                          ContainsSubstring("fills VD_ins_seq before both its ends are placed"));
    }
    SECTION("an insertion before the D it sits next to") {
        move_before(model.vd_ins.get(), model.d.get());
        CHECK_THROWS_WITH(refused(order).check_admissible(),
                          ContainsSubstring("places VD_ins_seq before the 5' end of D_gene_seq is placed"));
    }
    SECTION("two events creating the same segment") {
        const auto second_v = make_gene_choice(V_gene, {{"V2", "ACGT"}}, 11);
        order.push_back(second_v.get());
        CHECK_THROWS_WITH(refused(order).check_admissible(), ContainsSubstring("which an earlier event created"));
    }
    SECTION("a junction nothing fills") {
        order.erase(std::find(order.begin(), order.end(), model.dj_chain.get()));
        CHECK_THROWS_WITH(refused(order).check_admissible(),
                          ContainsSubstring("No event builds the content of DJ_ins_seq"));
    }
}

TEST_CASE("TraversalPlan: a deletion on a junction, after the chain filled it", "[traversal_plan]")
{
    // No model has one yet; the plan needs no notion of gene or junction to take it. The deletion
    // trims content the chain built, so the junction now has a factor of its own, and the gap
    // V->D can no longer be held at its id.
    const Vdj model;
    const SeqTypeRegistry &registry = vdj_seq_type_registry();
    const auto trims_junction = on(make_deletion(D_gene_seq, Three_prime, 0, 2, 11), registry, "VD_ins_seq");
    auto order = model.order();
    order.insert(std::find(order.begin(), order.end(), model.dj_ins.get()), trims_junction.get());

    const TraversalPlan plan(registry, order);
    CHECK_NOTHROW(plan.check_admissible());
    CHECK(plan.index().has_own_factor(VD_ins_seq));
    CHECK(plan.gap_key(gap(V_gene_seq, D_gene_seq)) >= registry.total_count());
    // It is checked against D, and widens the gap between the junction and D.
    const EventPlacement &placed = plan.placement(*trims_junction);
    CHECK(placed.placed_right.front() == D_gene_seq);
    CHECK(placed.gap_toward(Three_prime) == gap(VD_ins_seq, D_gene_seq));
    CHECK(plan.gap_key(gap(VD_ins_seq, D_gene_seq)) >= registry.total_count());

    SECTION("before the chain, it would trim content nobody built yet") {
        auto early = model.order();
        early.insert(std::find(early.begin(), early.end(), model.vd_chain.get()), trims_junction.get());
        CHECK_THROWS_WITH(TraversalPlan(registry, early).check_admissible(),
                          ContainsSubstring("modifies the content of VD_ins_seq before any event builds it"));
    }
}

TEST_CASE("TraversalPlan: one plan serves every copy of the model", "[traversal_plan]")
{
    // GenModel builds the plan once and shares it with every thread, whose events are copies.
    // copy() keeps an event's identifier, and the plan files placements by identifier.
    const Vdj model;
    const TraversalPlan plan(vdj_seq_type_registry(), model.order());

    const std::shared_ptr<Rec_Event> d_copy = model.d->copy();
    const std::shared_ptr<Rec_Event> v_3_copy = model.v_3->copy();
    REQUIRE(d_copy.get() != model.d.get());
    CHECK(plan.placement(*d_copy).placed_left == plan.placement(*model.d).placed_left);
    CHECK(plan.placement(*d_copy).enclosing_gap() == gap(V_gene_seq, J_gene_seq));
    CHECK(plan.placement(*v_3_copy).placed_right == std::vector<SeqTypeId>{D_gene_seq, J_gene_seq});

    SECTION("an event whose identifier is not in the order has no placement") {
        const auto stranger = make_gene_choice(V_gene, {{"V9", "ACGT"}}, 42);
        CHECK_THROWS_AS(plan.placement(*stranger), std::out_of_range);
    }
    SECTION("two events of one order may not share an identifier") {
        const auto twin = make_gene_choice(V_gene, {{"V2", "ACGT"}}, model.v->get_event_identifier());
        auto order = model.order();
        order.push_back(twin.get());
        CHECK_THROWS_WITH(TraversalPlan(vdj_seq_type_registry(), order), ContainsSubstring("shares its identifier"));
    }
}

TEST_CASE("resolve_flank_checks: the placed neighbours, 5' to 3', one safety layer per row",
          "[traversal_plan][tandem_d]")
{
    // Gene_choice and Deletion build their overlap checks with this one function. A tandem order
    // is the first to give a segment two placed neighbours on its 5' side, which is where the
    // order of the checks shows: Gene_choice takes the last 5' check as its nearest neighbour.
    const Tandem model;
    const TraversalPlan plan(model.r, model.order());
    const SeqTypeId d1 = model.id("D1_gene_seq");
    const SeqTypeId d2 = model.id("D2_gene_seq");
    SafetyMatrix safety(model.r);
    SeqOffsetsMap offsets(model.r);
    // The partners' facing ends, claimed by the events that placed them; V's 3' end twice, so
    // that its layer differs from D1's.
    offsets.request_layer(V_gene_seq, Three_prime);
    offsets.request_layer(V_gene_seq, Three_prime);
    offsets.request_layer(d1, Three_prime);
    offsets.request_layer(d2, Five_prime);
    offsets.request_layer(J_gene_seq, Five_prime);

    const auto partners = [](const std::vector<FlankCheck> &checks) {
        std::vector<SeqTypeId> ids;
        for (const FlankCheck &check : checks) {
            ids.push_back(check.partner_id);
        }
        return ids;
    };

    SECTION("a gene choice: both sides, 5' to 3'") {
        const std::vector<FlankCheck> checks =
                resolve_flank_checks(plan.placement(*model.d2), Undefined_side, safety, offsets);
        REQUIRE(partners(checks) == std::vector<SeqTypeId>{V_gene_seq, d1, J_gene_seq});
        CHECK_FALSE(checks[0].this_is_left);
        CHECK_FALSE(checks[1].this_is_left);
        CHECK(checks[2].this_is_left);
        CHECK(checks[1].partner_side == Three_prime);
        CHECK(checks[2].partner_side == Five_prime);
        CHECK(checks[1].safety_cell == safety.cell(d2, d1));
        CHECK(checks[0].partner_offset_layer == offsets.claimed_layer(V_gene_seq, Three_prime));
        CHECK(checks[1].partner_offset_layer == offsets.claimed_layer(d1, Three_prime));
        CHECK(checks[0].partner_offset_layer != checks[1].partner_offset_layer);
        // Three cells in three rows (V's, D1's and D2's): three claims.
        CHECK(checks[0].safety_layer == safety.claimed_layer(checks[0].safety_cell));
        CHECK(checks[2].safety_layer == safety.claimed_layer(checks[2].safety_cell));
    }

    SECTION("a 3' deletion: only the neighbours it trims toward, which share one row") {
        const std::vector<FlankCheck> checks =
                resolve_flank_checks(plan.placement(*model.d1_3), Three_prime, safety, offsets);
        REQUIRE(partners(checks) == std::vector<SeqTypeId>{d2, J_gene_seq});
        CHECK(checks[0].this_is_left);
        CHECK(checks[0].safety_cell.row == checks[1].safety_cell.row);
        CHECK(checks[0].safety_layer == checks[1].safety_layer);
        CHECK(checks[0].safety_layer == safety.claimed_layer(checks[0].safety_cell));
    }

    SECTION("a 5' deletion: only the neighbours to its 5' side, 5' to 3'") {
        const std::vector<FlankCheck> checks =
                resolve_flank_checks(plan.placement(*model.d2_5), Five_prime, safety, offsets);
        CHECK(partners(checks) == std::vector<SeqTypeId>{V_gene_seq, d1});
    }

    SECTION("the overlap verdict is oriented by the side the partner is on") {
        const std::vector<FlankCheck> checks =
                resolve_flank_checks(plan.placement(*model.d2), Undefined_side, safety, offsets);
        const JunctionGeometry::OffsetInterval low{2, 3};
        const JunctionGeometry::OffsetInterval high{10, 12};
        // D1 sits 5' of D2: D2's 5' end at 10..12 is clear of D1's 3' end at 2..3.
        CHECK(checks[1].overlap(high, low) == JunctionGeometry::Overlap::Safe);
        CHECK(checks[1].overlap(low, high) == JunctionGeometry::Overlap::Infeasible);
        // J sits 3' of D2: the other way round.
        CHECK(checks[2].overlap(low, high) == JunctionGeometry::Overlap::Safe);
    }
}
