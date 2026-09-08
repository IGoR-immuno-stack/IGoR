/*
 * test_dinucl_markov_iterate.cpp
 *
 *  Characterization tests for Dinucl_markov::iterate (task 2a).
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
 * Written against the *unmodified* event, before B7 replaces its hardcoded traversal specs and
 * two residual switches (plan step 2a). `Dinucl_markov::iterate` had **zero** coverage before
 * this file: nothing in the unit suite executed it at all (plan section 6.4).
 *
 * The event fills the placeholder nucleotides an Insertion allocated, copying them straight
 * from the read. Two rows of the matrix take a different shape:
 *
 * - like Insertion it never branches: one hand-off, or none if pruned. There is no realization
 *   to enumerate -- the "realization" is the read itself;
 * - it is the only event with a **reverse** traversal, and the only one that mutates a segment
 *   another event owns, in place, through a pointer.
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

/// The read every fixture below is aligned against, long enough for any window they take.
constexpr const char *kRead = "ACGTACGTACGTACGTACGTACGTACGT";

/**
 * A VD junction ready to be filled: V placed, the junction allocated as placeholders, and the
 * Insertion that allocated it registered so initialize_event() can size its index arrays.
 *
 * VD is the *forward* traversal: the anchor is V's 3' end, so the read window starts one
 * position past it and the seed nucleotide is V's last.
 */
struct VdFill {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Dinucl_markov> dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
    std::shared_ptr<Insertion> insertion = make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1);

    explicit VdFill(std::size_t junction_length = 3, Seq_Offset v_three_prime = 9,
                    const std::string &v_segment = "")
    {
        state.preset_segment(V_gene_seq, 0, v_three_prime,
                             v_segment.empty() ? segment_run(0, v_three_prime) : v_segment);
        state.preset_placeholders(VD_ins_seq, junction_length);
        state.add_event(insertion);
        //A flat dinucleotide marginal: every (previous, next) pair is equally probable, so a
        //probability that differs is the arithmetic and not the model.
        for (std::size_t i = 0; i != 32; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }
};

/**
 * A DJ junction: the *reverse* traversal. The anchor is J's 5' end, the read window ends one
 * position before it, and the seed nucleotide is J's first.
 */
struct DjFill {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Dinucl_markov> dinucl = make_dinucl_markov(DJ_ins_seq, /*event_id=*/0);
    std::shared_ptr<Insertion> insertion = make_insertion(DJ_ins_seq, 0, 6, /*event_id=*/1);

    explicit DjFill(std::size_t junction_length = 3, Seq_Offset j_five_prime = 14,
                    const std::string &j_segment = "")
    {
        state.preset_segment(J_gene_seq, j_five_prime, j_five_prime + 6,
                             j_segment.empty() ? segment_run(j_five_prime, j_five_prime + 6) : j_segment);
        state.preset_placeholders(DJ_ins_seq, junction_length);
        state.add_event(insertion);
        for (std::size_t i = 0; i != 32; ++i) {
            state.set_marginal(i, 0.5L);
        }
    }
};

} // namespace

TEST_CASE("Dinucl_markov: baseline fill, forward traversal (G9)", "[dinucl][iterate]")
{
    // V spans [0, 9] of "ACGTACGTACGT...", so its 3' end sits on read position 9 and the
    // junction takes positions 10, 11, 12 -- "GTA".
    VdFill fixture;
    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);

    REQUIRE(next->call_count() == 1);
    const auto &call = next->calls.front();
    REQUIRE(call.sequences.count(VD_ins_seq) == 1);
    CHECK(call.sequences.at(VD_ins_seq) == "GTA");

    // The anchor is handed on untouched, and the junction still has no offsets of its own.
    CHECK(call.three_prime(V_gene_seq) == 9);
    CHECK(call.offsets.count(VD_ins_seq) == 0);
}

TEST_CASE("Dinucl_markov: reverse traversal fills in the read's own orientation", "[dinucl][iterate]")
{
    // The DJ junction is filled backwards -- the read window is reversed, filled, and the
    // result reversed again -- because the Markov chain runs away from its anchor while the
    // segment is stored 5'->3'. What the next event sees must be the read window, unreversed.
    //
    // J spans [14, 20], so a 3-nucleotide junction occupies read positions 11, 12, 13.
    DjFill fixture;
    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);

    REQUIRE(next->call_count() == 1);
    CHECK(next->calls.front().sequences.at(DJ_ins_seq) == "TAC");

    // The double reversal is not a no-op the test could miss: reversing once would give "CAT".
    CHECK(next->calls.front().sequences.at(DJ_ins_seq) != "CAT");
}

TEST_CASE("Dinucl_markov: the read window is taken from the anchor's facing end", "[dinucl][iterate]")
{
    SECTION("Forward: it starts one position past the anchor's 3' end")
    {
        // Moving the anchor moves the window, and the filled junction with it.
        VdFill fixture(/*junction_length=*/3, /*v_three_prime=*/5);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().sequences.at(VD_ins_seq) == "GTA");   // positions 6, 7, 8
    }

    SECTION("Reverse: it ends one position before the anchor's 5' end")
    {
        // Chosen so the answer differs from the default fixture's: moving J from 14 to 12
        // moves the window from positions 11..13 to 9..11.
        DjFill fixture(/*junction_length=*/3, /*j_five_prime=*/12);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().sequences.at(DJ_ins_seq) == "CGT");   // positions 9, 10, 11
    }

    SECTION("A longer junction takes a longer window")
    {
        VdFill fixture(/*junction_length=*/5, /*v_three_prime=*/9);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().sequences.at(VD_ins_seq) == "GTACG");
    }
}

TEST_CASE("Dinucl_markov: probability chains the anchor's last nucleotide into the read", "[dinucl][iterate]")
{
    // The one piece of arithmetic worth pinning exactly. For a junction of n nucleotides the
    // contribution is a product of n conditional terms, indexed `base + previous * 4 + next`:
    //
    //   - the *first* term's `previous` is the anchor's own last nucleotide, not a read
    //     position -- this is the only place the anchor sequence is used at all;
    //   - every later term takes both nucleotides from the read.
    //
    // V spans [0, 9] of ACGTACGT..., so its last nucleotide is C (1) and the junction is
    // G(2) T(3) A(0). The three terms are therefore (1,2), (2,3) and (3,0).
    VdFill fixture;
    for (std::size_t i = 0; i != 32; ++i) {
        fixture.state.set_marginal(i, 0.0L);
    }
    fixture.state.set_marginal(1 * 4 + 2, 0.5L);
    fixture.state.set_marginal(2 * 4 + 3, 0.25L);
    fixture.state.set_marginal(3 * 4 + 0, 0.125L);

    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
    REQUIRE(next->call_count() == 1);
    CHECK(next->calls.front().scenario_proba == Approx(0.5 * 0.25 * 0.125));

    SECTION("The incoming probability multiplies through")
    {
        VdFill other;
        for (std::size_t i = 0; i != 32; ++i) {
            other.state.set_marginal(i, 0.0L);
        }
        other.state.set_marginal(1 * 4 + 2, 0.5L);
        other.state.set_marginal(2 * 4 + 3, 0.25L);
        other.state.set_marginal(3 * 4 + 0, 0.125L);
        other.state.set_scenario_proba(0.4);
        const auto other_next = call_iterate_recording(other.dinucl, other.state);
        REQUIRE(other_next->call_count() == 1);
        CHECK(other_next->calls.front().scenario_proba == Approx(0.4 * 0.5 * 0.25 * 0.125));
    }
}

TEST_CASE("Dinucl_markov: the seed really comes from the anchor", "[dinucl][iterate]")
{
    // A control for the section above: change *only* the anchor's last nucleotide and the
    // first term changes. Without this, a body that seeded from the read instead would still
    // pass every other assertion here, because the read's preceding position happens to be
    // adjacent to the window.
    // A V ending in A(0) makes the first pair (0,2); one ending in C(1) makes it (1,2).
    VdFill ends_in_c(/*junction_length=*/1, /*v_three_prime=*/9, "ACGTACGTAC");
    VdFill ends_in_a(/*junction_length=*/1, /*v_three_prime=*/9, "ACGTACGTAA");
    for (std::size_t i = 0; i != 32; ++i) {
        ends_in_c.state.set_marginal(i, 0.0L);
        ends_in_a.state.set_marginal(i, 0.0L);
    }
    ends_in_c.state.set_marginal(1 * 4 + 2, 0.5L);    // (C, G)
    ends_in_a.state.set_marginal(0 * 4 + 2, 0.25L);   // (A, G)

    const auto from_c = call_iterate_recording(ends_in_c.dinucl, ends_in_c.state);
    const auto from_a = call_iterate_recording(ends_in_a.dinucl, ends_in_a.state);
    REQUIRE(from_c->call_count() == 1);
    REQUIRE(from_a->call_count() == 1);
    CHECK(from_c->calls.front().scenario_proba == Approx(0.5));
    CHECK(from_a->calls.front().scenario_proba == Approx(0.25));

    // Both filled the same nucleotide from the read: only the probability differs.
    CHECK(from_c->calls.front().sequences.at(VD_ins_seq) == "G");
    CHECK(from_a->calls.front().sequences.at(VD_ins_seq) == "G");
}

TEST_CASE("Dinucl_markov: an empty junction is filled with nothing, not skipped", "[dinucl][iterate]")
{
    // Row 8. The event still runs, still writes its downstream bound and still hands off; the
    // probability contribution is the empty product. B7's skip-empty walk must keep this --
    // an empty junction is a state the model produces, not an error.
    VdFill fixture(/*junction_length=*/0);
    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);

    REQUIRE(next->call_count() == 1);
    REQUIRE(next->calls.front().sequences.count(VD_ins_seq) == 1);
    CHECK(next->calls.front().sequences.at(VD_ins_seq).empty());
    CHECK(next->calls.front().scenario_proba == Approx(1.0));
    CHECK(next->calls.front().downstream_bounds.count(VD_ins_seq) == 1);
}

TEST_CASE("Dinucl_markov: only placeholder positions are written", "[dinucl][iterate]")
{
    // The fill is guarded by `ins_seq.at(i) == -1`. A position already carrying a nucleotide
    // is left alone and contributes no probability term -- which is what makes the buffer
    // reusable across the scenarios that share it.
    VdFill fixture;
    // Overwrite the middle placeholder with a T(3) before the event runs.
    Int_Str *junction = const_cast<Int_Str *>(get_constructed_sequence(fixture.state, VD_ins_seq));
    REQUIRE(junction != nullptr);
    junction->at(1) = 3;

    for (std::size_t i = 0; i != 32; ++i) {
        fixture.state.set_marginal(i, 0.0L);
    }
    fixture.state.set_marginal(1 * 4 + 2, 0.5L);   // (anchor C, read G) -- position 0
    fixture.state.set_marginal(2 * 4 + 3, 0.25L);  // (G, T) -- would be position 1, if written
    fixture.state.set_marginal(3 * 4 + 0, 0.125L); // (T, A) -- position 2

    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
    REQUIRE(next->call_count() == 1);
    // Position 1 keeps the T it already had, which here happens to equal the read.
    CHECK(next->calls.front().sequences.at(VD_ins_seq) == "GTA");
    // ...and its term is absent from the product.
    CHECK(next->calls.front().scenario_proba == Approx(0.5 * 0.125));

    SECTION("The first position is guarded separately, and must be guarded too")
    {
        // The first position has its own copy of the guard, because its `previous` nucleotide
        // comes from the anchor rather than from the read. Overwriting position 1 above leaves
        // that copy untouched -- a mutation deleting it survived the whole suite until this
        // section existed.
        //
        // Position 0 is pre-set to A(0), which the read (G at position 10) disagrees with, so
        // an unguarded write is visible in the sequence and not only in the probability.
        VdFill first_filled;
        Int_Str *seq = const_cast<Int_Str *>(get_constructed_sequence(first_filled.state, VD_ins_seq));
        REQUIRE(seq != nullptr);
        seq->at(0) = 0;

        for (std::size_t i = 0; i != 32; ++i) {
            first_filled.state.set_marginal(i, 0.0L);
        }
        first_filled.state.set_marginal(1 * 4 + 2, 0.5L);    // (anchor C, read G) -- skipped
        first_filled.state.set_marginal(2 * 4 + 3, 0.25L);   // (G, T) -- position 1
        first_filled.state.set_marginal(3 * 4 + 0, 0.125L);  // (T, A) -- position 2

        const auto filled = call_iterate_recording(first_filled.dinucl, first_filled.state);
        REQUIRE(filled->call_count() == 1);
        CHECK(filled->calls.front().sequences.at(VD_ins_seq) == "ATA");
        CHECK(filled->calls.front().scenario_proba == Approx(0.25 * 0.125));

        // Later positions still take *both* nucleotides from the read, so the value sitting in
        // position 0 does not feed position 1's term.
        CHECK(filled->calls.front().scenario_proba != Approx(0.0));
    }
}

TEST_CASE("Dinucl_markov: downstream bound and memory layering", "[dinucl][iterate]")
{
    SECTION("The junction's bound is reset to 1")
    {
        // Dinucl_markov is the last event on its junction, so it has nothing left to bound:
        // it writes the neutral value rather than leaving the Insertion's estimate in place.
        VdFill fixture;
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().downstream_bounds.at(VD_ins_seq) == Approx(1.0));
    }

    SECTION("The write lands at the requested layer, leaving the layer below intact")
    {
        VdFill fixture;
        fixture.state.exploration.downstream_proba_map.set(VD_ins_seq, 0.125, 0);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(get_downstream_bound(fixture.state, VD_ins_seq, 0) == Approx(0.125));
        CHECK(get_downstream_bound(fixture.state, VD_ins_seq, 1) == Approx(1.0));
    }

}

TEST_CASE("Dinucl_markov: pruning", "[dinucl][iterate]")
{
    SECTION("A scenario below the threshold does not hand off")
    {
        VdFill fixture;
        fixture.state.set_scenario_proba(1e-9);
        fixture.state.set_pruning_threshold(1.0);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        CHECK(next->call_count() == 0);
    }

    SECTION("Positive control: the same scenario without the threshold")
    {
        VdFill fixture;
        fixture.state.set_scenario_proba(1e-9);
        const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
        CHECK(next->call_count() == 1);
    }
}

TEST_CASE("Dinucl_markov: the VJ arm behaves as the other two", "[dinucl][iterate]")
{
    // The third arm of the switch B7 collapses. It needs the VJ ordering, since a junction
    // event addresses seq_types that a VDJ registry orders differently.
    IterateTestState state = create_iterate_state(kRead, 1000, 32, vj_seq_type_registry());
    auto dinucl = make_dinucl_markov(VJ_ins_seq, /*event_id=*/0);
    state.preset_segment(V_gene_seq, 0, 9, segment_run(0, 9));
    state.preset_placeholders(VJ_ins_seq, 3);
    state.add_event(make_insertion(VJ_ins_seq, 0, 6, /*event_id=*/1));
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.5L);
    }

    const auto next = call_iterate_recording(dinucl, state);
    REQUIRE(next->call_count() == 1);
    // Forward traversal from V's 3' end, exactly as VD: read positions 10, 11, 12.
    CHECK(next->calls.front().sequences.at(VJ_ins_seq) == "GTA");
    CHECK(next->calls.front().downstream_bounds.at(VJ_ins_seq) == Approx(1.0));
    CHECK(next->calls.front().scenario_proba == Approx(0.5 * 0.5 * 0.5));
}

TEST_CASE("Dinucl_markov: an ambiguous read position is averaged, not indexed", "[dinucl][iterate]")
{
    // When either nucleotide of a pair is ambiguous the conditional probability cannot be read
    // straight from the marginal array, so the event uses dinuc_proba_matrix -- the mean over
    // the nucleotides the ambiguous code stands for -- and records -1 in the realization
    // indices instead of a marginal index.
    //
    // The matrix is built by update_event_internal_probas(), which GenModel calls and
    // initialize_event() does not, so the fixture builds it explicitly.
    IterateTestState state = create_iterate_state("ACGTACGTACNTACGTACGTACGT");
    auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
    state.preset_segment(V_gene_seq, 0, 9, segment_run(0, 9));
    state.preset_placeholders(VD_ins_seq, 1);
    state.add_event(make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1));

    // The anchor's last nucleotide is C(1); the read position under the junction is N. The
    // four entries the average is taken over are base + 1 * 4 + {A, C, G, T}.
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.0L);
    }
    state.set_marginal(1 * 4 + 0, 0.1L);
    state.set_marginal(1 * 4 + 1, 0.2L);
    state.set_marginal(1 * 4 + 2, 0.3L);
    state.set_marginal(1 * 4 + 3, 0.4L);

    std::unordered_map<Rec_Event_name, int> index_by_name{{dinucl->get_name(), 0}};
    dinucl->update_event_internal_probas(state.model.model_parameters, index_by_name);

    const auto next = call_iterate_recording(dinucl, state);
    REQUIRE(next->call_count() == 1);
    CHECK(next->calls.front().scenario_proba == Approx(0.25));

    // The ambiguous code is copied into the junction as it stands, not resolved to a base.
    CHECK(next->calls.front().sequences.at(VD_ins_seq) == "N");
}

TEST_CASE("Dinucl_markov: a seq_type with no traversal spec is rejected", "[dinucl][iterate]")
{
    // The specs come from a hardcoded switch over the three legacy junctions, so any other
    // seq_type yields an empty vector. B7 derives them from the registry instead, at which
    // point this becomes "a seq_type with no neighbours".
    IterateTestState state = create_iterate_state(kRead);
    auto dinucl = make_dinucl_markov(V_gene_seq, /*event_id=*/0);
    CHECK_THROWS_AS(call_iterate(dinucl, state), std::invalid_argument);
}

// ===========================================================================================
// The empty-anchor hazard (plan section 2.9, and now 7.12).
//
// Tagged [.] so it does NOT run by default. This is not squeamishness: the defect is an
// unguarded `previous_seq.back()` on an empty Int_Str -- a std::vector<int> whose data() is
// null -- so the failure mode is a segfault that takes the whole test binary down, which
// [!shouldfail] cannot catch and CI cannot survive. Run it deliberately:
//
//     ./build/bin/igor_tests "[dinucl][empty_anchor]"
//
// The assertion states what B7 owes: an anchor that carries no nucleotide must be rejected or
// skipped, never read. Once the skip-empty walk lands this becomes an ordinary test and the
// [.] tag comes off.
// ===========================================================================================
TEST_CASE("DEFECT (plan 7.12): a fully deleted anchor is read for a seed nucleotide",
          "[.][dinucl][empty_anchor][defect]")
{
    IterateTestState state = create_iterate_state(kRead);
    auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
    state.add_event(make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1));

    // A V deleted down to nothing: written, empty, and carrying the degenerate offsets that
    // say so. Deletion's 5' guard is a strict `>`, so removing exactly the whole segment is a
    // legal realization -- this state is reachable on the current corpus, not only under
    // tandem D.
    state.preset_segment(V_gene_seq, 0, -1, "");
    state.preset_placeholders(VD_ins_seq, 3);
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.5L);
    }

    CHECK_THROWS(call_iterate(dinucl, state));
}

TEST_CASE("DEFECT (plan 7.13): Dinucl_markov fills the Insertion's buffer instead of its own layer",
          "[dinucl][iterate][defect][!shouldfail]")
{
    // Every other event claims a memory layer for what it writes, so that a sibling scenario
    // restores the previous value on backtracking. Dinucl_markov does not: it takes the
    // pointer the Insertion stored in the constructed-sequence map and fills that buffer
    // through it, claiming nothing and writing the map not at all.
    //
    // It works today only because the two events behave as one -- neither branches, so there
    // is never a sibling to corrupt, and the Insertion re-assigns the buffer on its next call
    // anyway. That is a property of the current pair, not of the contract: an Insertion that
    // looped over lengths (an indel-aware error model) or a Dinucl_markov that branched over
    // ambiguous nucleotides would have siblings sharing one buffer, and the second would read
    // the first's nucleotides where it expects placeholders -- silently, since the
    // placeholder guard treats an already-written position as "someone filled this".
    //
    // Intended: the filled junction stands at this event's own layer, and the layer below
    // still holds the placeholders the Insertion wrote.
    VdFill fixture;
    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
    REQUIRE(next->call_count() == 1);

    const Int_Str *insertion_layer = get_constructed_sequence(fixture.state, VD_ins_seq, 0);
    REQUIRE(insertion_layer != nullptr);
    CHECK(int_str_to_nt(*insertion_layer) == "...");   // still unfilled, not 'N'

    const Int_Str *own_layer = get_constructed_sequence(fixture.state, VD_ins_seq, 1);
    REQUIRE(own_layer != nullptr);
    CHECK(int_str_to_nt(*own_layer) == "GTA");
}
