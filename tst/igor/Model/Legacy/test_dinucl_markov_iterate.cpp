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
        //The junction starts one position past the anchor, which is where the Insertion's
        //own derivation puts it. The segment itself is the event under test's to create.
        state.preset_junction(VD_ins_seq, v_three_prime + 1, junction_length);
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
        //Ends one position before the anchor: the mirror of VdFill, and the reason this
        //fixture exists at all.
        state.preset_junction(DJ_ins_seq, j_five_prime - static_cast<Seq_Offset>(junction_length),
                              junction_length);
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

    // The anchor is handed on untouched, and so is the junction's placement: this event reads
    // those offsets to size the segment and never writes one, which is what its
    // OffsetRole::None declares.
    CHECK(call.three_prime(V_gene_seq) == 9);
    CHECK(call.five_prime(VD_ins_seq) == 10);
    CHECK(call.three_prime(VD_ins_seq) == 12);
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

// The case that used to stand here -- "only placeholder positions are written" -- pinned the
// `ins_seq.at(i) == int_undefined` guard in iterate_common(), which let this event share one
// buffer with the Insertion across sibling scenarios. O12 (a') removed the sharing: the buffer
// is created here, per scenario, every position a placeholder by construction, which left the
// guard unreachable through iterate(). R11 deleted it, so iterate_common() now writes every
// position; the cases above that read the filled junction are what cover that.

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
    state.preset_junction(VJ_ins_seq, 10, 3);
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
    state.preset_junction(VD_ins_seq, 10, 1);
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

TEST_CASE("Dinucl_markov: a junction with nothing to seed from is rejected", "[dinucl][iterate]")
{
    // Before B7 this read "a seq_type outside the hardcoded switch"; now it is the general
    // statement, and there are three ways to fail it.
    SECTION("Nothing on the side the chain runs from")
    {
        // V_gene_seq is first in the ordering, so a chain seeding from its left has no anchor.
        IterateTestState state = create_iterate_state(kRead);
        auto dinucl = make_dinucl_markov(V_gene_seq, /*event_id=*/0, Three_prime);
        CHECK_THROWS_AS(call_iterate(dinucl, state), std::invalid_argument);
    }

    SECTION("No direction declared")
    {
        // event_side is what says which way the Markov chain runs. Without it the ordering
        // alone cannot choose between the two neighbours -- both are equally adjacent.
        IterateTestState state = create_iterate_state(kRead);
        auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
        dinucl->set_event_side(Undefined_side);
        CHECK_THROWS_AS(call_iterate(dinucl, state), std::invalid_argument);
    }

    SECTION("A seq_type the registry does not know")
    {
        IterateTestState state = create_iterate_state(kRead);
        auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
        dinucl->set_seq_type_id(kNoSeqType);
        CHECK_THROWS_AS(call_iterate(dinucl, state), std::invalid_argument);
    }
}

TEST_CASE("Dinucl_markov: the direction comes from the event, the anchor from the ordering",
          "[dinucl][iterate]")
{
    // The split B7 rests on. The registry says which segments are adjacent; event_side says
    // which of the two seeds the chain. Building a VD junction that seeds from its *right*
    // neighbour is not a model IGoR ships, but it is a model the code must now express --
    // and it proves the anchor is not being recovered from the seq_type name.
    IterateTestState state = create_iterate_state(kRead);
    auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0, /*chain_side=*/Five_prime);
    state.preset_segment(D_gene_seq, 14, 18, segment_run(14, 18));
    state.preset_junction(VD_ins_seq, 11, 3);
    state.add_event(make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1));
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.5L);
    }

    const auto next = call_iterate_recording(dinucl, state);
    REQUIRE(next->call_count() == 1);
    // Reverse traversal off D's 5' end at 14: read positions 11, 12, 13.
    CHECK(next->calls.front().sequences.at(VD_ins_seq) == "TAC");
}

// ===========================================================================================
// The empty-anchor hazard (plan section 2.9, then 7.12; repaired by R2).
//
// This was a [.]-hidden case for two weeks, and not out of squeamishness: the defect was an
// unguarded `previous_seq.back()` on an empty Int_Str -- a std::vector<int> whose data() is
// null -- so the failure mode was a segfault that took the whole test binary down, which
// [!shouldfail] cannot express and CI cannot survive. It could not even carry [dinucl], since
// Catch2 runs a hidden test whenever a filter names one of its tags.
//
// R2 makes the anchor's emptiness a throw, so the case is now an ordinary one and runs with
// the suite. What section 7.12 decided, and what these assert: throw, do not discard. A
// scenario whose anchor was fully deleted is geometrically legitimate, so discarding it would
// remove probability mass the model should account for and leave no trace; the user has to
// see it. And an anchor carrying no nucleotide must be rejected or skipped, *never read* --
// so a junction that needs no seed is not rejected at all.
//
// R2 is the throw half only. The occupancy-skipping walk that would seed from the next
// non-empty segment instead is deferred to B10 (sections 7.11 and 9).
// ===========================================================================================
TEST_CASE("Dinucl_markov: a fully deleted anchor is rejected, not read",
          "[dinucl][iterate][empty_anchor]")
{
    IterateTestState state = create_iterate_state(kRead);
    auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
    state.add_event(make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1));

    // A V deleted down to nothing: written, empty, and carrying the degenerate offsets that
    // say so. Deletion's 5' guard is a strict `>`, so removing exactly the whole segment is a
    // legal realization -- this state is reachable on the current corpus, not only under
    // tandem D.
    state.preset_segment(V_gene_seq, 0, -1, "");
    state.preset_junction(VD_ins_seq, 0, 3);
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.5L);
    }

    CHECK_THROWS(call_iterate(dinucl, state));
}

TEST_CASE("Dinucl_markov: an empty anchor is not read when the junction needs no seed",
          "[dinucl][iterate][empty_anchor]")
{
    // The other half of "rejected or skipped, never read", and the reason the guard asks
    // whether a seed is needed rather than whether the anchor is empty. A zero-length
    // junction chooses no nucleotide, so it never touches the anchor and stays perfectly
    // scoreable -- rejecting it would discard a legitimate scenario, which is the very thing
    // the throw exists to avoid doing silently.
    IterateTestState state = create_iterate_state(kRead);
    auto dinucl = make_dinucl_markov(VD_ins_seq, /*event_id=*/0);
    state.add_event(make_insertion(VD_ins_seq, 0, 6, /*event_id=*/1));

    state.preset_segment(V_gene_seq, 0, -1, "");
    state.preset_junction(VD_ins_seq, 0, 0);
    for (std::size_t i = 0; i != 32; ++i) {
        state.set_marginal(i, 0.5L);
    }

    // Recording rather than CHECK_NOTHROW: not throwing is the weaker claim, and it would
    // also hold if the event silently discarded the scenario. What R2 owes is that the
    // scenario is handed on.
    const auto next = call_iterate_recording(dinucl, state);
    REQUIRE(next->call_count() == 1);
    CHECK(next->calls.front().sequences.at(VD_ins_seq).empty());
}

TEST_CASE("Dinucl_markov: the junction it writes is one it created, at a layer it claimed",
          "[dinucl][iterate][layers]")
{
    // Plan section 7.13, repaired by R1 under O12's decision (a'). What stood here was a
    // [!shouldfail] case asking for the weaker of the two available repairs: the filled
    // junction at this event's own layer, the Insertion's placeholders still readable at the
    // layer below. (a') is stronger -- there are no placeholders to read, because the
    // partially-constructed segment never exists.
    //
    // The defect was that Dinucl_markov took the Int_Str * the Insertion had left in the
    // sequence map and filled that buffer through it, claiming nothing and writing the map not
    // at all. It worked only because the two events behave as one: neither branches, so there
    // was never a sibling to corrupt, and the Insertion re-assigned the buffer on its next call
    // anyway. That is a property of the current pair and not of the layer contract -- an
    // Insertion looping over lengths, or a Dinucl_markov branching over an ambiguous read
    // position, would give two siblings one buffer, and the second would read the first's
    // nucleotides where it expects placeholders. Silently: the placeholder guard reads an
    // already-written position as "someone filled this".
    VdFill fixture;
    // Nothing stands here before the event runs. The Insertion wrote offsets, not a segment.
    REQUIRE_FALSE(fixture.state.scenario.constructed_sequences.exists(VD_ins_seq));

    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
    REQUIRE(next->call_count() == 1);

    // Written at the layer this event claimed, which is the ownership statement the harness
    // checks for every other writer and had to waive for this pair.
    const int claimed = fixture.state.scenario.constructed_sequences.claimed_layer(VD_ins_seq);
    const Int_Str *own_layer = get_constructed_sequence(fixture.state, VD_ins_seq, claimed);
    REQUIRE(own_layer != nullptr);
    CHECK(int_str_to_nt(*own_layer) == "GTA");
}

// ===========================================================================================
// The completed-scenario invariant.
//
// int_undefined is a placeholder inside one scenario, never a value a consumer can interpret,
// so by the time a scenario reaches the error rate every position must be determined.
// Rec_Event::iterate_wrap_up asserts it at the leaf, under NDEBUG only -- the default build is
// RelWithDebInfo, so these test the predicate directly rather than the assertion.
// ===========================================================================================

TEST_CASE("first_unfilled_segment: reports the placeholder, not the ambiguity code", "[dinucl][invariant]")
{
    Seq_type_str_p_map sequences(legacy_seq_type_registry(), 4);

    SECTION("Nothing written at all")
    {
        CHECK(first_unfilled_segment(sequences) == kNoSeqType);
    }

    SECTION("Determined segments, ambiguity codes included")
    {
        // int_N is *determined*: the read says this position is ambiguous. It must not be
        // reported -- conflating it with a placeholder is exactly the confusion int_undefined
        // was moved to 15 to prevent.
        Int_Str gene = nt2int("ACGT");
        Int_Str ambiguous = nt2int("ANNT");
        claim_layer_zero(sequences, V_gene_seq);
        claim_layer_zero(sequences, D_gene_seq);
        sequences.set(V_gene_seq, &gene, 0);
        sequences.set(D_gene_seq, &ambiguous, 0);
        CHECK(first_unfilled_segment(sequences) == kNoSeqType);
    }

    SECTION("A junction still holding placeholders")
    {
        Int_Str gene = nt2int("ACGT");
        Int_Str junction(3, int_undefined);
        claim_layer_zero(sequences, V_gene_seq);
        claim_layer_zero(sequences, VD_ins_seq);
        sequences.set(V_gene_seq, &gene, 0);
        sequences.set(VD_ins_seq, &junction, 0);
        CHECK(first_unfilled_segment(sequences) == static_cast<SeqTypeId>(VD_ins_seq));
    }

    SECTION("A single unfilled position among filled ones is enough")
    {
        Int_Str partly = nt2int("ACGT");
        partly.at(2) = int_undefined;
        claim_layer_zero(sequences, DJ_ins_seq);
        sequences.set(DJ_ins_seq, &partly, 0);
        CHECK(first_unfilled_segment(sequences) == static_cast<SeqTypeId>(DJ_ins_seq));
    }

    SECTION("An empty segment and a null pointer are not unfilled")
    {
        // Both are legitimate states -- actively absent, and never written -- and neither
        // holds a placeholder.
        Int_Str empty;
        claim_layer_zero(sequences, VD_ins_seq);
        claim_layer_zero(sequences, DJ_ins_seq);
        sequences.set(VD_ins_seq, &empty, 0);
        sequences.set(DJ_ins_seq, nullptr, 0);
        CHECK(first_unfilled_segment(sequences) == kNoSeqType);
    }
}

TEST_CASE("Dinucl_markov leaves no unfilled position behind", "[dinucl][invariant]")
{
    // The invariant as the production pair actually maintains it. Under O12 (a') it is a
    // stronger statement than it was: there is no window in which the junction exists and
    // holds placeholders, so the event cannot hand one on even in principle. The precondition
    // below is what changed -- it used to assert that the placeholders were there.
    VdFill fixture;
    REQUIRE_FALSE(fixture.state.scenario.constructed_sequences.exists(VD_ins_seq));
    REQUIRE(first_unfilled_segment(fixture.state.scenario.constructed_sequences) == kNoSeqType);

    const auto next = call_iterate_recording(fixture.dinucl, fixture.state);
    REQUIRE(next->call_count() == 1);
    CHECK(first_unfilled_segment(fixture.state.scenario.constructed_sequences) == kNoSeqType);
}

TEST_CASE("first_unplaced_segment_end: the offsets half of the same invariant", "[dinucl][invariant]")
{
    // One invariant, two halves: by the time a scenario reaches a leaf every seq_type in the
    // model must have both a sequence and its offsets, not necessarily written by the same
    // event. This half could not be written before R3 -- while Insertion recorded no offsets
    // it would have fired on every scenario that has a junction, which is every scenario.

    SECTION("Nothing placed at all")
    {
        Seq_offsets_map offsets(vdj_seq_type_registry(), 4);
        const auto [id, side] = first_unplaced_segment_end(offsets);
        CHECK(id == static_cast<SeqTypeId>(V_gene_seq));
        CHECK(side == Five_prime);
    }

    SECTION("A complete VDJ scenario")
    {
        Seq_offsets_map offsets(vdj_seq_type_registry(), 4);
        for (const Seq_type type : {V_gene_seq, VD_ins_seq, D_gene_seq, DJ_ins_seq, J_gene_seq}) {
            claim_layer_zero(offsets, type, Five_prime);
            claim_layer_zero(offsets, type, Three_prime);
            offsets.set(type, Five_prime, 0, 0);
            offsets.set(type, Three_prime, 0, 0);
        }
        CHECK(first_unplaced_segment_end(offsets).first == kNoSeqType);
    }

    SECTION("One end placed and not the other")
    {
        // The key is per (SeqTypeId, Seq_side) and not per segment, because get_offset_role is
        // side-taking: the two ends can in principle be created by different events.
        Seq_offsets_map offsets(vdj_seq_type_registry(), 4);
        for (const Seq_type type : {V_gene_seq, VD_ins_seq, D_gene_seq, DJ_ins_seq, J_gene_seq}) {
            claim_layer_zero(offsets, type, Five_prime);
            claim_layer_zero(offsets, type, Three_prime);
            offsets.set(type, Five_prime, 0, 0);
            offsets.set(type, Three_prime, 0, 0);
        }
        Seq_offsets_map missing_three_prime(vdj_seq_type_registry(), 4);
        for (const Seq_type type : {V_gene_seq, VD_ins_seq, D_gene_seq, DJ_ins_seq, J_gene_seq}) {
            claim_layer_zero(missing_three_prime, type, Five_prime);
            missing_three_prime.set(type, Five_prime, 0, 0);
            if (type != DJ_ins_seq) {
                claim_layer_zero(missing_three_prime, type, Three_prime);
                missing_three_prime.set(type, Three_prime, 0, 0);
            }
        }
        const auto [id, side] = first_unplaced_segment_end(missing_three_prime);
        CHECK(id == static_cast<SeqTypeId>(DJ_ins_seq));
        CHECK(side == Three_prime);
    }

    SECTION("An empty junction is placed, not absent")
    {
        // `3' == 5' - 1` is the encoding, and the reason this predicate cannot be written as
        // "no width" -- an insertion that inserted nothing has placed both its ends.
        Seq_offsets_map offsets(vdj_seq_type_registry(), 4);
        for (const Seq_type type : {V_gene_seq, D_gene_seq, DJ_ins_seq, J_gene_seq}) {
            claim_layer_zero(offsets, type, Five_prime);
            claim_layer_zero(offsets, type, Three_prime);
            offsets.set(type, Five_prime, 0, 0);
            offsets.set(type, Three_prime, 0, 0);
        }
        claim_layer_zero(offsets, VD_ins_seq, Five_prime);
        claim_layer_zero(offsets, VD_ins_seq, Three_prime);
        offsets.set(VD_ins_seq, Five_prime, 11, 0);
        offsets.set(VD_ins_seq, Three_prime, 10, 0);
        CHECK(first_unplaced_segment_end(offsets).first == kNoSeqType);
    }

    SECTION("A VJ model is not held to the segments it does not have")
    {
        // The registry pins all six legacy names whatever the model is, so a VJ model carries
        // D_gene_seq and both flanking junctions as ids no event in it will ever place. The
        // sweep is over registry.ordering(), which is the model's actual segment layout; over
        // every registered id this would fire on every VJ scenario.
        Seq_offsets_map offsets(vj_seq_type_registry(), 4);
        for (const Seq_type type : {V_gene_seq, VJ_ins_seq, J_gene_seq}) {
            claim_layer_zero(offsets, type, Five_prime);
            claim_layer_zero(offsets, type, Three_prime);
            offsets.set(type, Five_prime, 0, 0);
            offsets.set(type, Three_prime, 0, 0);
        }
        CHECK(first_unplaced_segment_end(offsets).first == kNoSeqType);
    }
}

TEST_CASE("Dinucl_markov: a junction no Seq_type enum names still resolves", "[dinucl][iterate]")
{
    // The capability B7 exists for. A tandem-D model's D1D2 junction has no entry in the
    // Seq_type enum, so nothing about it can be recovered from a switch -- but the ordering
    // places it like any other segment and the event's own side says which way its chain runs.
    //
    // The event is told its neighbours, exactly as Model_Parms::finalize() tells them: the
    // registry is read in one place, and the event only combines that answer with its own
    // direction.
    SeqTypeRegistry tandem;
    tandem.register_legacy_seq_types();
    tandem.set_ordered_types({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "D1D2_ins_seq",
                              "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"});
    tandem.freeze();

    const SeqTypeId junction_id = tandem.id("D1D2_ins_seq");
    Dinucl_markov dinucl(VD_ins_seq);
    dinucl.set_seq_type("D1D2_ins_seq");
    dinucl.set_seq_type_id(junction_id);
    dinucl.set_adjacent_segments(tandem.left_neighbor(junction_id), tandem.right_neighbor(junction_id));
    dinucl.set_event_side(Three_prime);

    const DinuclTraversalSpec forward = dinucl.get_junction();
    CHECK(forward.target_id == junction_id);
    CHECK(forward.anchor_id == tandem.id("D_gene_seq"));
    CHECK(forward.anchor_side == Three_prime);
    // No Seq_type for either name, so the generation path's handles stay unset -- and are
    // flagged as such rather than left at a plausible-looking default.
    CHECK_FALSE(forward.legacy_enums_valid);

    SECTION("The same event seeded from the other side anchors on the other neighbour")
    {
        dinucl.set_event_side(Five_prime);
        CHECK(dinucl.get_junction().anchor_id == tandem.id("D2_gene_seq"));
    }

    SECTION("Deriving it twice gives the same answer")
    {
        // Nothing is stored, so there is no resolution step to run twice and no state to keep
        // in sync -- which is the point of deriving rather than resolving.
        CHECK(dinucl.get_junction().anchor_id == forward.anchor_id);
        CHECK(dinucl.get_junction().target_id == forward.target_id);
    }

    SECTION("Without adjacency there is no junction")
    {
        // What an event holds before Model_Parms::finalize() has run.
        Dinucl_markov unresolved(VD_ins_seq);
        unresolved.set_seq_type("D1D2_ins_seq");
        unresolved.set_seq_type_id(junction_id);
        unresolved.set_event_side(Three_prime);
        CHECK(unresolved.get_junction().anchor_id == kNoSeqType);
    }
}
