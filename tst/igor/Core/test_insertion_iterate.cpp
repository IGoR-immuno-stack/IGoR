/*
 * test_insertion_iterate.cpp
 *
 *  Characterization tests for Insertion::iterate (task 1a).
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
 * Written against the *unmodified* event, before B6 collapses its three seq_type branches
 * (plan step 1a). Follows docs/ITERATE_TEST_GUIDE.md: one TEST_CASE per pattern, one SECTION
 * per instance, every guard paired with a positive control.
 *
 * Insertion is the smallest iterate() body and the only one that never branches: its
 * realization is *derived* from where its neighbours sit, so a call produces exactly one
 * hand-off or none at all. Two rows of the matrix therefore change shape:
 *
 * - row 2 (realization enumeration) becomes "one hand-off, or none", and row 3's
 *   "realizations do not compound" has nothing to say, since there is only ever one;
 * - row 5 (neighbour dependence) has no chosen/unchosen distinction to make -- Insertion
 *   never asks whether a neighbour was chosen, it reads its offsets and assumes. The
 *   section pinning what happens when they are absent stands in for it.
 *
 * Everything else is present, including row 10, which the harness checks for free.
 */

#include "test_utils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <stdexcept>
#include <string>

using namespace IgorTestUtils;
using Catch::Approx;

namespace {

/**
 * A VD junction ready to iterate: V and D placed, the insertion event under test, and the
 * Dinucl_markov it needs to initialize at all.
 *
 * Offsets are inclusive read positions, so the junction length the event derives is
 * `off(D, 5') - off(V, 3') - 1` -- 3 for the defaults below.
 */
struct VdJunction {
    IterateTestState state = create_iterate_state("ACGTACGTACGTACGTACGTACGTACGT");
    std::shared_ptr<Insertion> insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/0);
    std::shared_ptr<Dinucl_markov> dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/1);

    explicit VdJunction(Seq_Offset v_three_prime = 10, Seq_Offset d_five_prime = 14)
    {
        state.preset_segment(V_gene_seq, 0, v_three_prime, segment_run(0, v_three_prime));
        state.preset_segment(D_gene_seq, d_five_prime, d_five_prime + 4,
                             segment_run(d_five_prime, d_five_prime + 4));
        state.add_downstream_event(dinucl);
        // A flat, non-zero marginal: every insertion length is equally probable, so a
        // discarded scenario is never explained by a zero in the model.
        for (std::size_t i = 0; i != 16; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }
};

/// The DJ arm: same shape as VdJunction, between D's 3' end and J's 5' end.
struct DjJunction {
    IterateTestState state = create_iterate_state("ACGTACGTACGTACGTACGTACGTACGT");
    std::shared_ptr<Insertion> insertion = make_insertion(DJ_ins_seq, 0, 6, /*event_id=*/0);

    explicit DjJunction(Seq_Offset d_three_prime = 10, Seq_Offset j_five_prime = 14)
    {
        state.preset_segment(D_gene_seq, d_three_prime - 4, d_three_prime,
                             segment_run(d_three_prime - 4, d_three_prime));
        state.preset_segment(J_gene_seq, j_five_prime, j_five_prime + 6,
                             segment_run(j_five_prime, j_five_prime + 6));
        state.add_downstream_event(make_dinucl_markov(DJ_ins_seq, /*event_id=*/1));
        for (std::size_t i = 0; i != 16; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }
};

/// The VJ arm: the whole V-to-J span, as a model with no D declares it.
struct VjJunction {
    IterateTestState state = create_iterate_state("ACGTACGTACGTACGTACGTACGTACGT", 1000, 32,
                                                  vj_seq_type_registry());
    std::shared_ptr<Insertion> insertion = make_insertion(VJ_ins_seq, 0, 6, /*event_id=*/0);

    explicit VjJunction(Seq_Offset v_three_prime = 10, Seq_Offset j_five_prime = 14)
    {
        state.preset_segment(V_gene_seq, 0, v_three_prime, segment_run(0, v_three_prime));
        state.preset_segment(J_gene_seq, j_five_prime, j_five_prime + 6,
                             segment_run(j_five_prime, j_five_prime + 6));
        state.add_downstream_event(make_dinucl_markov(VJ_ins_seq, /*event_id=*/1));
        for (std::size_t i = 0; i != 16; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }
};

} // namespace

TEST_CASE("Insertion: baseline write (G4)", "[insertion][iterate]")
{
    VdJunction fixture;
    const auto next = call_iterate_recording(fixture.insertion, fixture.state);

    REQUIRE(next->call_count() == 1);
    const auto &call = next->calls.front();

    // Where the junction sits: the three read positions strictly between the two neighbours.
    // The *segment* is not here, and its absence is the point -- under O12 (a') this event
    // decides the geometry and the Dinucl_markov that follows creates the sequence from it,
    // so no partially-constructed segment is ever handed on. What used to be asserted here
    // was three undetermined positions decoding to "...".
    REQUIRE(call.offsets.count(VD_ins_seq) == 1);
    CHECK(call.five_prime(VD_ins_seq) == 11);
    CHECK(call.three_prime(VD_ins_seq) == 13);
    CHECK_FALSE(has_constructed_sequence(fixture.state, VD_ins_seq));

    // It still writes no mismatch list, and no longer owes one: a list belongs to a segment
    // compared against the read, and this event constructs no segment (plan section 6.9, where
    // that defect is recorded as dissolved by the rescope rather than repaired).

    // The neighbours are handed on untouched.
    CHECK(call.three_prime(V_gene_seq) == 10);
    CHECK(call.five_prime(D_gene_seq) == 14);
}

TEST_CASE("Insertion: the length is derived from the neighbours (G9)", "[insertion][iterate]")
{
    // The pattern B6 generalises: an insertion has no realization of its own to enumerate.
    // It reads the two offsets around it and the arithmetic picks the realization.
    SECTION("Moving either neighbour moves the length")
    {
        VdJunction fixture(/*v_three_prime=*/10, /*d_five_prime=*/16);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        // Five positions now: 11..15.
        CHECK(next->calls.front().five_prime(VD_ins_seq) == 11);
        CHECK(next->calls.front().three_prime(VD_ins_seq) == 15);
    }

    SECTION("Adjacent neighbours give a placed-but-empty junction")
    {
        // Row 8 of the matrix, and the case B10's absent-segment semantics has to keep
        // distinct: the junction is placed and empty, which is not the same as never placed.
        // Emptiness is `3' == 5' - 1` -- a pair of offsets, not their absence, which is the
        // encoding that makes the distinction expressible at all.
        VdJunction fixture(/*v_three_prime=*/13, /*d_five_prime=*/14);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        REQUIRE(next->calls.front().offsets.count(VD_ins_seq) == 1);
        CHECK(next->calls.front().five_prime(VD_ins_seq) == 14);
        CHECK(next->calls.front().three_prime(VD_ins_seq) == 13);
    }

    SECTION("A neighbour that was never placed is an error, not an empty junction")
    {
        // Insertion assumes both its neighbours are already constructed; it has no
        // "not yet processed" branch. LayeredArray turns the assumption into a throw
        // rather than a read of uninitialized storage.
        IterateTestState state = create_iterate_state("ACGTACGTACGTACGTACGT");
        auto insertion = make_insertion(VD_ins_seq, 0, 6, 0);
        state.preset_segment(V_gene_seq, 0, 10, segment_run(0, 10));
        state.add_downstream_event(make_dinucl_markov(VD_ins_seq, 1));
        for (std::size_t i = 0; i != 16; ++i) {
            state.set_marginal(i, 0.5L);
        }
        CHECK_THROWS_AS(call_iterate(insertion, state), std::out_of_range);
    }
}

TEST_CASE("Insertion: exactly one hand-off, or none", "[insertion][iterate]")
{
    // The "non-branching" row: unlike Gene_choice, one call to iterate() can only ever
    // produce a single scenario, because the realization is determined rather than chosen.
    SECTION("One realization in range hands off once")
    {
        VdJunction fixture;
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 1);
    }

    SECTION("A length above the realization set is discarded")
    {
        // 0..6 insertions declared; this junction needs 9.
        VdJunction fixture(/*v_three_prime=*/0, /*d_five_prime=*/10);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 0);

        // No layer promise is owed on this path. The event claimed its downstream_proba_map
        // layer once at initialize_event() and leaves it unwritten here, which is sound
        // precisely because nothing downstream runs: the contract is "written before handing
        // off", and there is no hand-off. B6 must keep that property -- claiming per call
        // instead of once would turn this into the section 7.9 defect.
        CHECK(next->layer_violations.empty());
    }

    SECTION("The DJ and VJ arms discard out-of-range lengths the same way")
    {
        // Added after the coverage report showed the false side of `proba_contribution != 0`
        // was reached only on the VD arm: the three arms were covered for their hand-off but
        // not for their discard, so a rewrite could have lost two of the three guards
        // undetected.
        DjJunction dj(/*d_three_prime=*/0, /*j_five_prime=*/10);
        CHECK(call_iterate_recording(dj.insertion, dj.state)->call_count() == 0);

        VjJunction vj(/*v_three_prime=*/0, /*j_five_prime=*/10);
        CHECK(call_iterate_recording(vj.insertion, vj.state)->call_count() == 0);
    }

    SECTION("A negative length is discarded by the same lookup")
    {
        // Plan section 7.7: Insertion performs *no* overlap check. Neighbours that overlap
        // produce a negative junction length, and it is the realization lookup -- not a
        // geometric test -- that throws the scenario away. The distinction matters for B6:
        // the guard being replaced is a set-membership test, not an interval test.
        VdJunction fixture(/*v_three_prime=*/14, /*d_five_prime=*/14);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 0);
    }
}

TEST_CASE("Insertion: probability", "[insertion][iterate]")
{
    SECTION("The contribution is the marginal at base_index + realization_index")
    {
        VdJunction fixture;   // 3 insertions
        fixture.state.set_marginal(3, 0.25L);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().scenario_proba == Approx(0.25));
    }

    SECTION("The incoming probability multiplies through")
    {
        VdJunction fixture;
        fixture.state.set_marginal(3, 0.25L);
        fixture.state.set_scenario_proba(0.4);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().scenario_proba == Approx(0.1));
    }

    SECTION("base_index is honoured")
    {
        VdJunction fixture;
        fixture.state.set_base_index(/*event_id=*/0, /*base_index=*/8);
        fixture.state.set_marginal(3, 0.25L);    // would be read if base_index were ignored
        fixture.state.set_marginal(8 + 3, 0.75L);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().scenario_proba == Approx(0.75));
    }

    SECTION("A zero-probability length is discarded like an unreachable one")
    {
        // proba_contribution == 0 is the single gate for both, so a length the model declares
        // with probability zero takes the same path as one it cannot produce at all. Pinned
        // because a generic body that separates "not a realization" from "probability zero"
        // would change which scenarios reach the next event.
        VdJunction fixture;
        fixture.state.set_marginal(3, 0.0L);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 0);

        // Positive control: the same fixture with a non-zero entry at that one index.
        VdJunction control;
        control.state.set_marginal(3, 0.125L);
        const auto control_next = call_iterate_recording(control.insertion, control.state);
        CHECK(control_next->call_count() == 1);
    }
}

TEST_CASE("Insertion: downstream bound and memory layering", "[insertion][iterate]")
{
    SECTION("The junction's bound is the best probability achievable at that length")
    {
        // junction_length_best_proba_map is built by the reverse pass walking what remains of
        // the queue -- here the Dinucl_markov. With a flat 0.5 marginal the best a junction of
        // n nucleotides can reach is 0.5 for the insertion itself times 0.5 per filled
        // position, so a 3-nucleotide junction is 0.5^4.
        VdJunction fixture;
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        REQUIRE(next->calls.front().downstream_bounds.count(VD_ins_seq) == 1);
        CHECK(next->calls.front().downstream_bounds.at(VD_ins_seq) == Approx(0.0625));
    }

    SECTION("An empty junction is bounded by the insertion probability alone")
    {
        // Nothing to fill, so the dinucleotide factor drops out entirely. This is the row the
        // bound would silently lose if a rewrite folded the two factors together.
        VdJunction fixture(/*v_three_prime=*/13, /*d_five_prime=*/14);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().downstream_bounds.at(VD_ins_seq) == Approx(0.5));
    }

    SECTION("The write lands at the requested layer, leaving the layer below intact")
    {
        // Row 9. An upstream event's value at layer 0 must still be readable after the
        // insertion has written its own, because a downstream reader of `layer - 1` depends
        // on it.
        VdJunction fixture;
        fixture.state.exploration.downstream_proba_map.set(VD_ins_seq, 0.125, 0);

        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);

        CHECK(get_downstream_bound(fixture.state, VD_ins_seq, 0) == Approx(0.125));
        CHECK(get_downstream_bound(fixture.state, VD_ins_seq, 1)
              == Approx(next->calls.front().downstream_bounds.at(VD_ins_seq)));
    }
}

TEST_CASE("Insertion: pruning", "[insertion][iterate]")
{
    SECTION("A scenario below the threshold does not hand off")
    {
        VdJunction fixture;
        fixture.state.set_marginal(3, 1e-6L);
        fixture.state.set_pruning_threshold(1.0);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 0);
    }

    SECTION("Positive control: the same scenario without the threshold")
    {
        VdJunction fixture;
        fixture.state.set_marginal(3, 1e-6L);
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        CHECK(next->call_count() == 1);
    }
}

TEST_CASE("Insertion: the three junction branches are one body", "[insertion][iterate]")
{
    // B6 collapses the VD / DJ / VJ chain into a single generic body. These sections pin
    // that the three arms already behave identically once the neighbour pair is substituted
    // -- the anchor the collapse has to reproduce.

    SECTION("DJ: the same shape, between D's 3' end and J's 5' end")
    {
        DjJunction fixture;
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().five_prime(DJ_ins_seq) == 11);
        CHECK(next->calls.front().three_prime(DJ_ins_seq) == 13);
        CHECK(next->calls.front().downstream_bounds.count(DJ_ins_seq) == 1);
    }

    SECTION("VJ: the same shape again, spanning a model with no D")
    {
        VjJunction fixture;
        const auto next = call_iterate_recording(fixture.insertion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().five_prime(VJ_ins_seq) == 11);
        CHECK(next->calls.front().three_prime(VJ_ins_seq) == 13);
        CHECK(next->calls.front().downstream_bounds.count(VJ_ins_seq) == 1);
    }

    SECTION("The marginal index agrees across the arms, despite two derivations")
    {
        // VD and DJ compute new_index as `base_index + event_realizations.at(to_string(n)).index`
        // -- a string conversion and a hash lookup, in the hot loop, carrying its own FIXME --
        // while VJ uses the `realization_index` iterate_common() already resolved. The two must
        // agree, and B6 keeps only the second. add_to_marginals() is the cheapest place to
        // observe new_index, since it is the only thing that reads it.
        constexpr std::size_t kArraySize = 64;

        VdJunction vd;
        const auto vd_next = call_iterate_recording(vd.insertion, vd.state);
        REQUIRE(vd_next->call_count() == 1);
        Marginal_array_p vd_marginals(new long double[kArraySize]());
        vd.insertion->add_to_marginals(1.0L, vd_marginals);

        VjJunction vj;
        const auto vj_next = call_iterate_recording(vj.insertion, vj.state);
        REQUIRE(vj_next->call_count() == 1);
        Marginal_array_p vj_marginals(new long double[kArraySize]());
        vj.insertion->add_to_marginals(1.0L, vj_marginals);

        // Both junctions are three nucleotides long, so both must credit the same index.
        int vd_index = -1;
        int vj_index = -1;
        for (std::size_t i = 0; i != kArraySize; ++i) {
            if (vd_marginals[i] != 0.0L) {
                CHECK(vd_index == -1);   // exactly one index credited
                vd_index = static_cast<int>(i);
            }
            if (vj_marginals[i] != 0.0L) {
                CHECK(vj_index == -1);
                vj_index = static_cast<int>(i);
            }
        }
        CHECK(vd_index == 3);
        CHECK(vj_index == vd_index);
    }
}

TEST_CASE("Insertion: a junction with no segment on one side is rejected", "[insertion][iterate]")
{
    // B6 resolves the neighbours from the ordering at initialize_event() instead of naming
    // them per seq_type, which introduces a case the string comparisons could not have: an
    // insertion sitting at the end of the ordering, with nothing beyond it to bound the
    // junction. Caught at initialization rather than as a kNoSeqType subscript in the hot loop.
    static const SeqTypeRegistry truncated = [] {
        SeqTypeRegistry built;
        built.register_legacy_seq_types();
        built.set_ordered_types({"V_gene_seq", "VD_ins_seq"});   // nothing to the 3' side
        built.freeze();
        return built;
    }();

    IterateTestState state = create_iterate_state("ACGTACGTACGTACGT", 1000, 32, truncated);
    auto insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/0);
    state.preset_segment(V_gene_seq, 0, 10, segment_run(0, 10));
    state.add_downstream_event(make_dinucl_markov(VD_ins_seq, /*event_id=*/1));

    CHECK_THROWS_AS(call_iterate(insertion, state), std::runtime_error);
}

TEST_CASE("Insertion: the seq_type name and its registry id must agree", "[insertion][iterate]")
{
    // An event carries its seq_type twice -- the serialization name, and the registry id
    // resolved from it. Everything in the scenario is keyed by the id, so a disagreement does
    // not fail: it aliases this event's segment, offsets and bounds onto another seq_type's
    // keys and comes back as a wrong answer. That is the shape of the trap B2 hit with VJ.
    //
    // Checked once at initialize_event(), where both identities are in hand. Before B6 the
    // same site rejected only an *unrecognised* name, which caught less and did it by
    // enumerating the three legacy junctions.
    SECTION("A name that resolves to a different segment than the id")
    {
        IterateTestState state = create_iterate_state("ACGTACGTACGTACGT");
        auto insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/0);
        insertion->set_seq_type("V_gene_seq");   // id still points at VD_ins_seq
        CHECK_THROWS_AS(call_iterate(insertion, state), std::runtime_error);
    }

    SECTION("A name no registry knows")
    {
        IterateTestState state = create_iterate_state("ACGTACGTACGTACGT");
        auto insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/0);
        insertion->set_seq_type("not_a_seq_type");
        CHECK_THROWS_AS(call_iterate(insertion, state), std::runtime_error);
    }

    SECTION("An id never resolved against a registry")
    {
        // What an event holds before Model_Parms::finalize() runs. Reaching iterate() in that
        // state would subscript every scenario map with kNoSeqType.
        IterateTestState state = create_iterate_state("ACGTACGTACGTACGT");
        auto insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/0);
        insertion->set_seq_type_id(kNoSeqType);
        CHECK_THROWS_AS(call_iterate(insertion, state), std::runtime_error);
    }
}


// ===========================================================================================
// What used to be this file's [!shouldfail] block. Four confirmed defects stood here; R1+R3
// closed three of them and O12's rescope dissolved the fourth.
//
//   - no offsets for the junction, and none for an empty one -- repaired, and asserted by the
//     cases above rather than pinned separately: writing the offsets is not a property the
//     event has on the side, it is what the event now does.
//   - the segment written on an unrequested layer -- the write moved to Dinucl_markov, which
//     claims its layer, so the waiver in call_iterate_recording() is gone and every event is
//     now held to the ownership rule. The case below is what would notice it coming back.
//   - no mismatch list for the segment it creates -- **deleted rather than made to pass**.
//     Under (a') this event creates no sequence, so it needs no list for one. The reasoning
//     that made it a defect (a constructed segment implies a comparison against the read, and
//     "empty" is a different statement from "absent") now applies to Dinucl_markov, and is
//     recorded against it rather than reproduced here.
// ===========================================================================================

TEST_CASE("Insertion: everything it writes lands on a layer it owns", "[insertion][iterate][layers]")
{
    // The complement of the layer contract: a written layer must have been requested.
    // Requesting is what makes a layer *owned* -- a write with no request leaves claimed and
    // current out of step, so a downstream reader of `layer - 1` is trusting storage nobody
    // promised, and the contract itself is silent on that key because it only inspects keys
    // whose claimed mark the event raised.
    //
    // This event was the one standing exception, waived by name in call_iterate_recording()
    // so that the rule could be enforced everywhere else. The waiver is gone; this case is
    // what turns red if the exception returns.
    VdJunction fixture;
    const auto next = call_iterate_recording(fixture.insertion, fixture.state);
    REQUIRE(next->call_count() == 1);

    for (const OwnershipViolation &violation : next->ownership_violations) {
        UNSCOPED_INFO(violation.describe());
    }
    CHECK(next->ownership_violations.empty());
}
