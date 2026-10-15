/**
 * @file test_tandem_junctions.cpp
 * @brief Gene_choice and Deletion set up their checks and their junctions on a V-D1-D2-J
 * layout (plan T3).
 *
 * Both used to find the segments they are checked against through the legacy names V, D and J,
 * and the key of the junction they bound through the legacy junction seq types. A tandem model's
 * D1 and D2 are neither, and its D1->J gap has no legacy name at all, so initialize_event() could
 * not set either event up. Both now read their placement and their gap keys from the engine's
 * plan (TraversalPlan), which the harness builds over the order it runs the events in.
 *
 * The events are built on legacy names and moved onto the tandem segments, as
 * test_dinucl_markov_iterate.cpp does: on this branch an event cannot be built on a tandem name
 * (T1).
 */

#include "test_utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>

using namespace IgorTestUtils;

namespace {

const std::string kRead = "ACGTACGTACGTACGTACGTACGT"; // 24 nt

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

SeqTypeId id_of(const char *seq_type) { return tandem_registry().id(seq_type); }

/// Move an event built on a legacy name onto a segment of the tandem layout.
template <class Event>
std::shared_ptr<Event> on(std::shared_ptr<Event> event, const char *seq_type)
{
    event->set_seq_type(seq_type);
    event->set_seq_type_id(id_of(seq_type));
    event->update_event_name();
    return event;
}

std::string read_run(Seq_Offset five, Seq_Offset three)
{
    return kRead.substr(static_cast<std::size_t>(five), static_cast<std::size_t>(three - five + 1));
}

/// V at 0..7, D1 at 8..11 and J at 16..23, all placed; an insertion and its chain pending on each
/// junction it is given. Every marginal is 0.5, so a junction of length 0 is bounded at 0.5: the
/// chance of inserting nothing.
struct TandemModel {
    IterateTestState state = create_iterate_state(kRead, 1000, 32, tandem_registry());
    std::shared_ptr<Gene_choice> v_stub = make_gene_choice(V_gene, {{"V1", read_run(0, 7)}}, /*id=*/1);
    std::shared_ptr<Gene_choice> d1_stub =
            on(make_gene_choice(D_gene, {{"D1", read_run(8, 11)}}, /*id=*/2), "D1_gene_seq");
    std::shared_ptr<Gene_choice> j_stub = make_gene_choice(J_gene, {{"J1", read_run(16, 23)}}, /*id=*/3);

    TandemModel()
    {
        for (const auto &stub : {v_stub, d1_stub, j_stub}) {
            state.add_event(stub);
            state.mark_chosen(stub);
        }
        state.preset_segment(V_gene_seq, 0, 7, read_run(0, 7));
        state.preset_segment(id_of("D1_gene_seq"), 8, 11, read_run(8, 11));
        state.preset_segment(J_gene_seq, 16, 23, read_run(16, 23));
        for (std::size_t i = 0; i != 64; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }

    void add_pending_junction(const char *junction, int event_id)
    {
        state.add_downstream_event(on(make_insertion(VD_ins_seq, 0, 3, event_id), junction));
        state.add_downstream_event(on(make_dinucl_markov(VD_ins_seq, event_id + 1, Three_prime), junction));
    }

    SeqTypeId key(const char *left, const char *right) const
    {
        return state.plan().gap_key(SegmentSpan::gap(id_of(left), id_of(right)));
    }
};

}  // namespace

TEST_CASE("Tandem D: D2's gene choice is checked against D1 and bounds the gaps on both sides",
          "[tandem_d][gene_choice][iterate]")
{
    TandemModel model;
    model.add_pending_junction("D1D2_ins_seq", 10);
    model.add_pending_junction("DJ_ins_seq", 12);

    // Two placements: D2a fills 12..15, exactly between D1 and J; D2b at 10..13 overlaps D1.
    auto d2 = on(make_gene_choice(D_gene, {{"D2a", read_run(12, 15)}, {"D2b", read_run(10, 13)}}, /*id=*/0,
                                  /*fixed=*/false),
                 "D2_gene_seq");
    model.state.set_alignments(D_gene, {create_perfect_alignment("D2a", 12, 4), create_perfect_alignment("D2b", 10, 4)});

    const auto next = call_iterate_recording(d2, model.state);
    REQUIRE(next->call_count() == 1);

    // D1 is D2's nearest built neighbour on the left, J on the right: both gaps are bounded at
    // the insertion that fills them, and the D1->J gap D2 splits is set neutral at its own key.
    // Had D1 not been placed in the plan, the left gap would have run from V, which is never an
    // open gap in this order.
    const auto &bounds = next->calls.front().downstream_bounds;
    CHECK(model.key("D1_gene_seq", "D2_gene_seq") == id_of("D1D2_ins_seq"));
    CHECK(model.key("D2_gene_seq", "J_gene_seq") == id_of("DJ_ins_seq"));
    CHECK_THAT(bounds.at(model.key("D1_gene_seq", "D2_gene_seq")), Catch::Matchers::WithinRel(0.5, 1e-12));
    CHECK_THAT(bounds.at(model.key("D2_gene_seq", "J_gene_seq")), Catch::Matchers::WithinRel(0.5, 1e-12));
    CHECK(bounds.at(model.key("D1_gene_seq", "J_gene_seq")) == 1.0);
    CHECK_THROWS_AS(model.key("V_gene_seq", "D2_gene_seq"), std::out_of_range);
}

TEST_CASE("Tandem D: D1's 3' deletion is checked against D2 and widens the D1->D2 gap",
          "[tandem_d][deletion][iterate]")
{
    TandemModel model;
    // Two realizations, so its generated name differs from D1's stub, which has one: twin names
    // are T2's.
    auto d2_stub = on(make_gene_choice(D_gene, {{"D2", read_run(12, 15)}, {"D2x", "A"}}, /*id=*/4), "D2_gene_seq");
    model.state.add_event(d2_stub);
    model.state.mark_chosen(d2_stub);
    model.state.preset_segment(id_of("D2_gene_seq"), 12, 15, read_run(12, 15));
    model.state.preset_safety(id_of("D1_gene_seq"), id_of("D2_gene_seq"), false);
    model.add_pending_junction("D1D2_ins_seq", 10);

    // -1 is a one-nucleotide palindrome: D1's 3' end moves to 12, onto D2's 5' end.
    auto deletion = on(make_deletion(D_gene_seq, Three_prime, -1, 2, /*id=*/0), "D1_gene_seq");

    const auto next = call_iterate_recording(deletion, model.state);
    // 2, 1 and 0 leave a gap of 2, 1 and 0 for the pending insertion to fill; -1 collides with D2.
    // Checked against nothing, -1 would be handed off too.
    REQUIRE(next->call_count() == 3);
    const SeqTypeId d1_to_d2 = model.key("D1_gene_seq", "D2_gene_seq");
    for (const ScenarioSnapshot &call : next->calls) {
        CHECK(call.downstream_bounds.count(d1_to_d2) == 1);
    }
    // A gap of 0, the last realization handed off, is bounded at the chance of inserting nothing.
    CHECK_THAT(next->calls.back().downstream_bounds.at(d1_to_d2), Catch::Matchers::WithinRel(0.5, 1e-12));
}
