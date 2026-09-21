/*
 * test_deletion_iterate.cpp
 *
 *  Characterization tests for Deletion::iterate (task 4a).
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
 * Written against the *unmodified* event, before B5 collapses its four arms (plan step 4b).
 * Follows docs/ITERATE_TEST_GUIDE.md: one TEST_CASE per pattern, one SECTION per instance,
 * every guard paired with a positive control.
 *
 * `Deletion::iterate` was at **0 % line coverage** before this file: no unit test of any kind
 * executed it, and the only thing standing between a change in those 1004 lines and a wrong
 * answer was the regression corpus. It is also the largest of the four bodies, and the one
 * whose four arms -- V 3', D 5', D 3', J 5' -- differ in ways that are *not* all incidental:
 *
 * - **which end moves**: a 3' deletion retreats (`offset - k`) and truncates the tail, a 5'
 *   deletion advances (`offset + k`) and truncates the head. The mismatch list is trimmed to
 *   the surviving *contiguous* subrange accordingly -- a prefix for the 3' arms, a suffix for
 *   the 5' ones;
 * - **whether the template may vanish**: V and J forbid full deletion (`size() > value`), D
 *   allows it on both sides (`size() >= value`). That is a modelling decision, not an accident,
 *   and it is asserted as an asymmetry below;
 * - **how many prune checks there are**: V and J check twice (a `break`, then a `continue`),
 *   D once (a `continue`);
 * - **what guards the read boundary**: three arms guard, the J arm does not. See the last
 *   TEST_CASE.
 *
 * Every number below was read off a run of the unmodified event, never from a comment or from
 * what the code looks like it should do (guide section 8, rule 1).
 */

#include "test_utils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace IgorTestUtils;
using Catch::Approx;

namespace {

/// 28 nucleotides of ACGT. Periodic on purpose: a segment preset from it is one an aligner
/// could actually have produced, which the palindrome sections need -- a negative deletion
/// compares the reverse-complemented template against the read itself, so a template that
/// disagrees with the read would make those mismatch lists meaningless.
const std::string kRead = "ACGTACGTACGTACGTACGTACGTACGT";

/// The read's own nucleotides over the inclusive range [five, three]. Use this rather than
/// segment_run(), which always starts at 'A' and so agrees with the read only at offsets that
/// are multiples of four (guide trap 1, one layer down: the *length* is right either way, the
/// content is not).
std::string read_run(Seq_Offset five, Seq_Offset three)
{
    if (three < five) { return {}; }
    return kRead.substr(static_cast<std::size_t>(five), static_cast<std::size_t>(three - five + 1));
}

/**
 * The junction-length bound these fixtures produce at `distance`.
 *
 * Every event in the fold carries a flat 0.5 marginal and the Dinucl_markov contributes 0.5 per
 * inserted nucleotide, so a path that reaches `distance` by inserting `n` nucleotides is worth
 * `0.5^contributors * 0.5^n`. The fold keeps the best, which is the shortest insertion that can
 * reach the distance: `n = distance` when the deletions contribute nothing, and `n = 0`
 * whenever the distance is zero or negative (the deletions can always make up the difference).
 *
 * `contributors` is the number of realization-enumerating events the fold walks: the event
 * under test plus every length-affecting event *after* it in the queue. It differs per fixture,
 * which is why it is a parameter rather than a constant.
 */
double junction_bound_at(int distance, int contributors)
{
    return std::pow(0.5, static_cast<double>(std::max(distance, 0) + contributors));
}

/// Single_error_rate(rate)'s upper bound for a segment with `mismatches` errors over
/// `matches` genomic positions. Spelled out rather than hard-coded so a section says what it
/// means.
double err_bound(double rate, int mismatches, int matches)
{
    return std::pow(rate / 3.0, mismatches) * std::pow(1.0 - rate, matches);
}

// ---------------------------------------------------------------------------------------
// Fixtures. One per (segment, side) arm, plus the two a model without the neighbour needs.
//
// Each places the segment the deletion trims, the neighbour it checks against, and the events
// that feed the junction-length fold. Trap 4 of the guide is the reason for the last: without
// downstream events the junction map collapses to {0: 1.0} and every geometry that is not
// exactly adjacent is discarded by the junction guard, which turns any section into a test of
// that guard.
// ---------------------------------------------------------------------------------------

/// A V 3' deletion with D already chosen: the VD flank is the one it widens.
/// Fold contributors: this deletion, the D 5' deletion, the VD insertion -> 3.
struct VDel {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;
    std::shared_ptr<Gene_choice> d_stub = make_gene_choice(D_gene, {{"D1", "ACGTA"}}, /*id=*/1);
    static constexpr int kContributors = 3;

    VDel(Seq_Offset v_three = 10, Seq_Offset d_five = 14, int min_del = 0, int max_del = 4,
         Seq_Offset v_five = 0, const std::vector<std::size_t> &v_mis = {}, int ins_max = 20)
        : deletion(make_deletion(V_gene_seq, Three_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_safety(V_gene_seq, D_gene_seq, false);
        state.preset_segment(V_gene_seq, v_five, v_three, read_run(v_five, v_three), v_mis);
        state.add_event(d_stub);
        state.mark_chosen(d_stub);
        state.preset_segment(D_gene_seq, d_five, d_five + 4, read_run(d_five, d_five + 4));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 4, /*id=*/2));
        state.add_downstream_event(make_insertion(VD_ins_seq, 0, ins_max, /*id=*/3));
        state.add_downstream_event(make_dinucl_markov(VD_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A V 3' deletion with no neighbour chosen at all: no safety check, no junction, no bound.
struct VDelBare {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;

    VDelBare(Seq_Offset v_five, Seq_Offset v_three, int min_del, int max_del)
        : deletion(make_deletion(V_gene_seq, Three_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_segment(V_gene_seq, v_five, v_three, read_run(v_five, v_three));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A D 5' deletion with V already chosen.
/// Fold contributors: this deletion, the VD insertion -> 2. The V 3' deletion is registered but
/// sits *before* this event, so it bounds the safety check without entering the fold.
struct D5Del {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;
    std::shared_ptr<Gene_choice> v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACGTACG"}}, /*id=*/1);
    static constexpr int kContributors = 2;

    D5Del(Seq_Offset v_three = 10, Seq_Offset d_five = 14, Seq_Offset d_three = 18,
          int min_del = 0, int max_del = 4, const std::vector<std::size_t> &d_mis = {},
          int ins_max = 20)
        : deletion(make_deletion(D_gene_seq, Five_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_safety(V_gene_seq, D_gene_seq, false);
        state.preset_segment(D_gene_seq, d_five, d_three, read_run(d_five, d_three), d_mis);
        state.add_event(v_stub);
        state.mark_chosen(v_stub);
        state.preset_segment(V_gene_seq, 0, v_three, read_run(0, v_three));
        state.add_event(make_deletion(V_gene_seq, Three_prime, 0, 4, /*id=*/5));
        state.add_downstream_event(make_insertion(VD_ins_seq, 0, ins_max, /*id=*/3));
        state.add_downstream_event(make_dinucl_markov(VD_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A D 5' deletion with nothing chosen beside it: no safety check and no junction, which is
/// what a geometry the fold cannot represent needs in order to reach the arm's own guards.
struct D5DelBare {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;

    D5DelBare(Seq_Offset d_five, Seq_Offset d_three, int min_del, int max_del)
        : deletion(make_deletion(D_gene_seq, Five_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_segment(D_gene_seq, d_five, d_three, read_run(d_five, d_three));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A D 3' deletion with J already chosen. Fold contributors: this deletion, the DJ insertion.
struct D3Del {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;
    std::shared_ptr<Gene_choice> j_stub = make_gene_choice(J_gene, {{"J1", "ACGTACGT"}}, /*id=*/1);
    static constexpr int kContributors = 2;

    D3Del(Seq_Offset d_five = 10, Seq_Offset d_three = 14, Seq_Offset j_five = 18,
          int min_del = 0, int max_del = 4, const std::vector<std::size_t> &d_mis = {},
          int ins_max = 20)
        : deletion(make_deletion(D_gene_seq, Three_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_safety(D_gene_seq, J_gene_seq, false);
        state.preset_segment(D_gene_seq, d_five, d_three, read_run(d_five, d_three), d_mis);
        state.add_event(j_stub);
        state.mark_chosen(j_stub);
        state.preset_segment(J_gene_seq, j_five, 27, read_run(j_five, 27));
        state.add_event(make_deletion(J_gene_seq, Five_prime, 0, 4, /*id=*/5));
        state.add_downstream_event(make_insertion(DJ_ins_seq, 0, ins_max, /*id=*/3));
        state.add_downstream_event(make_dinucl_markov(DJ_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A D 3' deletion with nothing chosen beside it.
struct D3DelBare {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;

    D3DelBare(Seq_Offset d_five, Seq_Offset d_three, int min_del, int max_del)
        : deletion(make_deletion(D_gene_seq, Three_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_segment(D_gene_seq, d_five, d_three, read_run(d_five, d_three));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A J 5' deletion with D already chosen. Fold contributors: this deletion, the DJ insertion.
struct JDel {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;
    std::shared_ptr<Gene_choice> d_stub = make_gene_choice(D_gene, {{"D1", "ACGTA"}}, /*id=*/1);
    static constexpr int kContributors = 2;

    JDel(Seq_Offset d_three = 14, Seq_Offset j_five = 18, int min_del = 0, int max_del = 4,
         const std::vector<std::size_t> &j_mis = {}, int ins_max = 20)
        : deletion(make_deletion(J_gene_seq, Five_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_safety(D_gene_seq, J_gene_seq, false);
        state.preset_segment(J_gene_seq, j_five, 27, read_run(j_five, 27), j_mis);
        state.add_event(d_stub);
        state.mark_chosen(d_stub);
        state.preset_segment(D_gene_seq, d_three - 4, d_three, read_run(d_three - 4, d_three));
        state.add_event(make_deletion(D_gene_seq, Three_prime, 0, 4, /*id=*/5));
        state.add_downstream_event(make_insertion(DJ_ins_seq, 0, ins_max, /*id=*/3));
        state.add_downstream_event(make_dinucl_markov(DJ_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// A J 5' deletion with no neighbour chosen at all.
struct JDelBare {
    IterateTestState state = create_iterate_state(kRead);
    std::shared_ptr<Deletion> deletion;

    JDelBare(Seq_Offset j_five, Seq_Offset j_three, int min_del, int max_del)
        : deletion(make_deletion(J_gene_seq, Five_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_segment(J_gene_seq, j_five, j_three, read_run(j_five, j_three));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

/// The VJ arm: a model with no D, where a V deletion widens the V->J span.
/// Fold contributors: this deletion, the J 5' deletion, the VJ insertion -> 3.
struct VJDel {
    IterateTestState state = create_iterate_state(kRead, 1000, 32, vj_seq_type_registry());
    std::shared_ptr<Deletion> deletion;
    std::shared_ptr<Gene_choice> j_stub = make_gene_choice(J_gene, {{"J1", "ACGTACGT"}}, /*id=*/1);
    static constexpr int kContributors = 3;

    VJDel(Seq_Offset v_three = 10, Seq_Offset j_five = 14, int min_del = 0, int max_del = 4)
        : deletion(make_deletion(V_gene_seq, Three_prime, min_del, max_del, /*id=*/0))
    {
        state.preset_safety(V_gene_seq, J_gene_seq, false);
        state.preset_segment(V_gene_seq, 0, v_three, read_run(0, v_three));
        state.add_event(j_stub);
        state.mark_chosen(j_stub);
        state.preset_segment(J_gene_seq, j_five, 27, read_run(j_five, 27));
        state.add_downstream_event(make_deletion(J_gene_seq, Five_prime, 0, 4, /*id=*/2));
        state.add_downstream_event(make_insertion(VJ_ins_seq, 0, 20, /*id=*/3));
        state.add_downstream_event(make_dinucl_markov(VJ_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    }
};

} // namespace

// =======================================================================================
// Row 1 -- baseline write
// =======================================================================================

TEST_CASE("Deletion: baseline write (G4)", "[deletion][iterate]")
{
    // One realization's worth of state, in each of the four arms. The knob that distinguishes
    // them is *which end moves and in which direction*, which is the whole of B5's collapse.

    SECTION("V 3' -- the 3' end retreats and the tail is truncated")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/2, /*max_del=*/2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        const ScenarioSnapshot &call = next->calls.front();

        CHECK(call.five_prime(V_gene_seq) == 0);   // untouched
        CHECK(call.three_prime(V_gene_seq) == 8);  // 10 - 2
        CHECK(call.sequences.at(V_gene_seq) == read_run(0, 8));
        CHECK(call.mismatches.at(V_gene_seq).empty());

        // The neighbour is handed on exactly as it was found.
        CHECK(call.five_prime(D_gene_seq) == 14);
        CHECK(call.three_prime(D_gene_seq) == 18);
        CHECK(call.sequences.at(D_gene_seq) == read_run(14, 18));
    }

    SECTION("D 5' -- the 5' end advances and the head is truncated")
    {
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/2, /*max_del=*/2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        const ScenarioSnapshot &call = next->calls.front();

        CHECK(call.five_prime(D_gene_seq) == 16);  // 14 + 2
        CHECK(call.three_prime(D_gene_seq) == 18); // untouched
        CHECK(call.sequences.at(D_gene_seq) == read_run(16, 18));
        CHECK(call.mismatches.at(D_gene_seq).empty());
    }

    SECTION("D 3' -- the same segment, the other end")
    {
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/2, /*max_del=*/2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        const ScenarioSnapshot &call = next->calls.front();

        CHECK(call.five_prime(D_gene_seq) == 10);  // untouched
        CHECK(call.three_prime(D_gene_seq) == 12); // 14 - 2
        CHECK(call.sequences.at(D_gene_seq) == read_run(10, 12));
    }

    SECTION("J 5'")
    {
        JDel fixture(/*d_three=*/14, /*j_five=*/18, /*min_del=*/2, /*max_del=*/2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        const ScenarioSnapshot &call = next->calls.front();

        CHECK(call.five_prime(J_gene_seq) == 20);  // 18 + 2
        CHECK(call.three_prime(J_gene_seq) == 27); // untouched
        CHECK(call.sequences.at(J_gene_seq) == read_run(20, 27));
    }
}

// =======================================================================================
// Row 2 -- realization enumeration
// =======================================================================================

TEST_CASE("Deletion: realization enumeration (G4)", "[deletion][iterate]")
{
    SECTION("N realizations give N hand-offs, in decreasing deletion count")
    {
        // The order is not incidental: it is what makes the first prune check a `break`
        // rather than a `continue`, and what makes SpanProfile::best_for() a unit-stride
        // backwards scan (see SpanProfile.h). del_numb_compare sorts descending.
        VDel fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        for (std::size_t i = 0; i != 5; ++i) {
            const int deletions = 4 - static_cast<int>(i);
            INFO("hand-off " << i << ", " << deletions << " deletions");
            CHECK(next->calls[i].three_prime(V_gene_seq) == 10 - deletions);
        }
    }

    SECTION("Each hand-off carries its own state, not the last one written")
    {
        // The defect this row exists for: a body that writes through a shared member and
        // only publishes it once would hand every realization the same sequence.
        VDel fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        for (std::size_t i = 0; i != 5; ++i) {
            const int deletions = 4 - static_cast<int>(i);
            INFO("hand-off " << i);
            CHECK(next->calls[i].sequences.at(V_gene_seq) == read_run(0, 10 - deletions));
        }
    }

    SECTION("A one-realization event hands off exactly once")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/3, /*max_del=*/3);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        CHECK(next->call_count() == 1);
    }
}

// =======================================================================================
// Row 3 -- probability
// =======================================================================================

TEST_CASE("Deletion: probability (G4)", "[deletion][iterate]")
{
    SECTION("The incoming probability multiplies through")
    {
        VDel fixture;
        fixture.state.set_scenario_proba(0.25);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.scenario_proba == Approx(0.25 * 0.5));
        }
    }

    SECTION("Realizations do not compound, and base_index is honoured")
    {
        // Distinct marginals at 8..12 do both jobs at once: reading the wrong base index
        // gives 0.5 everywhere (the flat block the fixture writes at 0..63), and compounding
        // gives a decreasing product rather than the marginal itself.
        VDel fixture;
        fixture.state.set_base_index(0, 8);
        for (std::size_t i = 0; i != 5; ++i) {
            fixture.state.set_marginal(8 + i, 0.1L * static_cast<long double>(i + 1));
        }
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        // Enumerated in decreasing deletion count, so the first hand-off is realization
        // index 4 and reads marginal[8 + 4].
        const double expected[] = {0.5, 0.4, 0.3, 0.2, 0.1};
        for (std::size_t i = 0; i != 5; ++i) {
            INFO("hand-off " << i);
            CHECK(next->calls[i].scenario_proba == Approx(expected[i]));
        }
    }
}

// =======================================================================================
// Row 4 -- feasibility guards. One section per `continue`/`break`, each with a positive
// control differing in exactly one value.
// =======================================================================================

TEST_CASE("Deletion: the full-deletion guard is asymmetric between D and the flanking genes",
          "[deletion][iterate]")
{
    // A modelling decision, not a refactoring one: V and J test `size() > value` and so keep
    // at least one nucleotide of the template, D tests `size() >= value` on both sides and may
    // delete itself away entirely. Each section carries its own positive control -- the
    // realization one smaller, in the same fixture.

    SECTION("V forbids it")
    {
        // The V segment starts at read position 3, so a full deletion would still leave the
        // 3' offset inside the read: the guard is what rejects it, not the offset arithmetic.
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/7, /*max_del=*/8, /*v_five=*/3);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // 8 == size is rejected, 7 survives
        CHECK(next->calls.front().three_prime(V_gene_seq) == 3);
        CHECK(next->calls.front().sequences.at(V_gene_seq) == read_run(3, 3));
    }

    SECTION("J forbids it")
    {
        JDel fixture(/*d_three=*/16, /*j_five=*/18, /*min_del=*/9, /*max_del=*/10);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // 10 == size is rejected, 9 survives
        CHECK(next->calls.front().five_prime(J_gene_seq) == 27);
        CHECK(next->calls.front().sequences.at(J_gene_seq) == read_run(27, 27));
    }

    SECTION("D 5' allows it")
    {
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/4, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 2); // 5 == size survives
        CHECK(next->calls[0].sequences.at(D_gene_seq).empty());
        CHECK(next->calls[1].sequences.at(D_gene_seq) == read_run(18, 18));
    }

    SECTION("D 3' allows it")
    {
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/4, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 2);
        CHECK(next->calls[0].sequences.at(D_gene_seq).empty());
        CHECK(next->calls[1].sequences.at(D_gene_seq) == read_run(10, 10));
    }
}

TEST_CASE("Deletion: guards on the read boundary", "[deletion][iterate]")
{
    SECTION("D 5' -- the new 5' offset must stay inside the read")
    {
        // Deletion.cpp's `//THIS IS A TEMPORARY FIX //FIXME`. The D segment ends at the last
        // read position, so deleting all four of its nucleotides puts its 5' offset one past
        // the end.
        D5Del fixture(/*v_three=*/20, /*d_five=*/24, /*d_three=*/27, /*min_del=*/3, /*max_del=*/4);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // 4 puts 5' at 28, 3 puts it at 27
        CHECK(next->calls.front().five_prime(D_gene_seq) == 27);
    }

    SECTION("V negative -- the palindrome may not run off the end of the read")
    {
        VDelBare fixture(/*v_five=*/0, /*v_three=*/25, /*min_del=*/-3, /*max_del=*/-2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // -3 would reach offset 28, -2 reaches 27
        CHECK(next->calls.front().three_prime(V_gene_seq) == 27);
    }

    SECTION("V negative -- the palindrome may not be longer than the template")
    {
        VDelBare fixture(/*v_five=*/0, /*v_three=*/10, /*min_del=*/-12, /*max_del=*/-11);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // 12 > 11 nucleotides of template, 11 == 11 is allowed
        CHECK(next->calls.front().three_prime(V_gene_seq) == 21);
    }

    // The V arm's `if (v_3_new_offset < 0) continue;` has no section, and cannot have one:
    // it is unreachable. The guard above it already requires `value < size`, and with
    // `size == three_prime - five_prime + 1` and a 5' offset at or after the start of the read
    // that gives `value <= three_prime`, so the new offset is never negative. It becomes
    // reachable only for a segment whose length disagrees with its offsets, which no upstream
    // event produces (Gene_choice clips a V template that starts before the read). Recorded in
    // the plan for B5 rather than pinned here.
}

TEST_CASE("Deletion: overlap verdicts (G2/G3)", "[deletion][iterate]")
{
    // Three outcomes per flank, and the same three in every arm: the moving end is compared
    // against the interval the neighbour's own pending deletions can still reach.
    //
    //   past the far end of that interval -> Infeasible, discarded
    //   short of its near end             -> Safe, and the flag is raised
    //   inside it                         -> Undetermined, and the flag is lowered

    SECTION("V 3' vs D 5' -- all three, in one geometry")
    {
        // V's 3' end sweeps 16..20 while D's 5' end can reach anywhere in [14, 18].
        VDel fixture(/*v_three=*/20, /*d_five=*/14);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        // 20, 19 and 18 are at or past D's furthest reach: discarded.
        REQUIRE(next->call_count() == 2);
        CHECK(next->calls[0].three_prime(V_gene_seq) == 16);
        CHECK(next->calls[1].three_prime(V_gene_seq) == 17);
        // Both land inside the interval, so neither can be declared safe.
        CHECK(next->calls[0].safety.at({V_gene_seq, D_gene_seq}) == false);
        CHECK(next->calls[1].safety.at({V_gene_seq, D_gene_seq}) == false);
    }

    SECTION("V 3' vs D 5' -- Safe, when the V end clears the interval entirely")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.at({V_gene_seq, D_gene_seq}) == true);
        }
    }

    SECTION("D 5' vs V 3'")
    {
        // D's 5' end sweeps 5..10; V's 3' end can retreat to 6. The deletion range reaches one
        // past the template length on purpose: the *boundary* realization -- the one that lands
        // exactly on V's furthest retreat, and so separates `<=` from `<` -- sits at a negative
        // junction distance, and only a fold whose own deletions reach that far down has a
        // bound to offer it. A range one shorter leaves the boundary case discarded by the
        // junction guard instead, which is a section that no longer tests the comparison.
        D5Del fixture(/*v_three=*/10, /*d_five=*/5, /*d_three=*/9, /*min_del=*/0, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 4); // 6 is exactly at V's furthest retreat: discarded
        CHECK(next->calls[0].five_prime(D_gene_seq) == 10);
        CHECK(next->calls[3].five_prime(D_gene_seq) == 7);
        CHECK(next->calls[0].safety.at({V_gene_seq, D_gene_seq}) == false);
    }

    SECTION("D 3' vs J 5'")
    {
        // Same reasoning as the section above for the deletion range: the boundary realization
        // is the one worth having, and it needs the fold to reach a negative distance.
        D3Del fixture(/*d_five=*/10, /*d_three=*/24, /*j_five=*/20, /*min_del=*/0, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5); // 24 is exactly at J's furthest reach: discarded
        CHECK(next->calls[0].three_prime(D_gene_seq) == 19); // widest deletion first
        CHECK(next->calls[0].safety.at({D_gene_seq, J_gene_seq}) == true);  // clear of the interval
        CHECK(next->calls[4].three_prime(D_gene_seq) == 23);
        CHECK(next->calls[4].safety.at({D_gene_seq, J_gene_seq}) == false); // inside it
    }

    SECTION("J 5' vs D 3'")
    {
        JDel fixture(/*d_three=*/24, /*j_five=*/18, /*min_del=*/0, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 3); // 20 is exactly at D's furthest reach: discarded
        CHECK(next->calls[0].five_prime(J_gene_seq) == 23);
        CHECK(next->calls[2].five_prime(J_gene_seq) == 21);
        CHECK(next->calls[0].safety.at({D_gene_seq, J_gene_seq}) == false);
    }

    SECTION("Safe -- the other three arms, when the moving end clears the interval")
    {
        // The complement of the three sections above, one per arm. Without these, inverting a
        // Safe verdict into an Undetermined one is invisible: nothing else in the file reads
        // the flag these arms raise.
        {
            D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 5);
            for (const ScenarioSnapshot &call : next->calls) {
                CHECK(call.safety.at({V_gene_seq, D_gene_seq}) == true);
            }
        }
        {
            D3Del fixture;
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 5);
            for (const ScenarioSnapshot &call : next->calls) {
                CHECK(call.safety.at({D_gene_seq, J_gene_seq}) == true);
            }
        }
        {
            JDel fixture;
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 5);
            for (const ScenarioSnapshot &call : next->calls) {
                CHECK(call.safety.at({D_gene_seq, J_gene_seq}) == true);
            }
        }
    }

    SECTION("An already-safe flank is not re-checked, and the flag is carried forward")
    {
        // Exactly the geometry of the first section, with one value changed: the incoming
        // VD_safe flag. Nothing is discarded, because the comparison is skipped outright.
        VDel fixture(/*v_three=*/20, /*d_five=*/14);
        fixture.state.preset_safety(V_gene_seq, D_gene_seq, true);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.at({V_gene_seq, D_gene_seq}) == true);
        }
    }

    SECTION("...and likewise in the other three arms")
    {
        // Each of the four arms has its own copy of the short-circuit. Comparing against the
        // three sections of the overlap block above: same fixtures, incoming flag raised,
        // nothing discarded.
        {
            D5Del fixture(/*v_three=*/10, /*d_five=*/5, /*d_three=*/9, /*min_del=*/0, /*max_del=*/5);
            fixture.state.preset_safety(V_gene_seq, D_gene_seq, true);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            CHECK(next->call_count() == 5);
        }
        {
            D3Del fixture(/*d_five=*/10, /*d_three=*/24, /*j_five=*/20, /*min_del=*/0, /*max_del=*/5);
            fixture.state.preset_safety(D_gene_seq, J_gene_seq, true);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            CHECK(next->call_count() == 6);
        }
        {
            // Four rather than six: the two widest deletions leave a gap more negative than
            // this fold can represent, so the junction guard discards them whatever the
            // safety flag says. Three survive with the check on, four with it short-circuited.
            JDel fixture(/*d_three=*/24, /*j_five=*/18, /*min_del=*/0, /*max_del=*/5);
            fixture.state.preset_safety(D_gene_seq, J_gene_seq, true);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            CHECK(next->call_count() == 4);
        }
    }

    SECTION("The J arm checks both of its flanks too")
    {
        // The mirror of the V section below: J compares against D on one side and V on the
        // other, and reads two safety flags rather than one.
        JDel fixture;
        fixture.state.preset_safety(V_gene_seq, J_gene_seq, false);
        auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACG"}}, /*id=*/7);
        fixture.state.add_event(v_stub);
        fixture.state.mark_chosen(v_stub);
        fixture.state.preset_segment(V_gene_seq, 0, 6, read_run(0, 6));
        fixture.state.add_event(make_deletion(V_gene_seq, Three_prime, 0, 4, /*id=*/8));

        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.at({D_gene_seq, J_gene_seq}) == true);
            CHECK(call.safety.at({V_gene_seq, J_gene_seq}) == true);
        }
    }

    SECTION("...and the V flank can be the one that discards")
    {
        // One value changed from the section above: V reaches far enough right that J's 5' end
        // sits inside the interval its deletions can still cover.
        JDel fixture(/*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/5);
        fixture.state.preset_safety(V_gene_seq, J_gene_seq, false);
        auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACG"}}, /*id=*/7);
        fixture.state.add_event(v_stub);
        fixture.state.mark_chosen(v_stub);
        fixture.state.preset_segment(V_gene_seq, 0, 22, read_run(0, 22));
        fixture.state.add_event(make_deletion(V_gene_seq, Three_prime, 0, 4, /*id=*/8));

        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5); // the narrowest deletion lands on V's furthest reach

        // The widest deletion clears V's interval; the four narrower ones sit inside it.
        CHECK(next->calls[0].five_prime(J_gene_seq) == 23);
        CHECK(next->calls[0].safety.at({V_gene_seq, J_gene_seq}) == true);
        for (std::size_t i = 1; i != 5; ++i) {
            INFO("hand-off " << i);
            CHECK(next->calls[i].safety.at({V_gene_seq, J_gene_seq}) == false);
        }
    }

    SECTION("The V arm checks both of its flanks")
    {
        // D on the 3' side and J on the 5' side of the same V end. The VD comparison clears,
        // the VJ one does not, and it is the VJ verdict that decides.
        VDel fixture(/*v_three=*/14, /*d_five=*/16);
        fixture.state.preset_safety(V_gene_seq, J_gene_seq, false);
        auto j_stub = make_gene_choice(J_gene, {{"J1", "ACGTAC"}}, /*id=*/7);
        fixture.state.add_event(j_stub);
        fixture.state.mark_chosen(j_stub);
        fixture.state.preset_segment(J_gene_seq, 8, 13, read_run(8, 13));
        fixture.state.add_event(make_deletion(J_gene_seq, Five_prime, 0, 4, /*id=*/8));

        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 2);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.at({V_gene_seq, D_gene_seq}) == true);
            CHECK(call.safety.at({V_gene_seq, J_gene_seq}) == false);
        }
    }

    SECTION("An already-safe VJ flank short-circuits in both arms that read it")
    {
        // The V and J arms each carry a second short-circuit for the V<->J flank, distinct from
        // the one on their D flank. One value changed from the two sections above.
        {
            VDel fixture(/*v_three=*/14, /*d_five=*/16);
            fixture.state.preset_safety(V_gene_seq, J_gene_seq, true);
            auto j_stub = make_gene_choice(J_gene, {{"J1", "ACGTAC"}}, /*id=*/7);
            fixture.state.add_event(j_stub);
            fixture.state.mark_chosen(j_stub);
            fixture.state.preset_segment(J_gene_seq, 8, 13, read_run(8, 13));
            fixture.state.add_event(make_deletion(J_gene_seq, Five_prime, 0, 4, /*id=*/8));

            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 5); // nothing discarded: the flank is not compared
            for (const ScenarioSnapshot &call : next->calls) {
                CHECK(call.safety.at({V_gene_seq, J_gene_seq}) == true);
            }
        }
        {
            JDel fixture(/*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/5);
            fixture.state.preset_safety(V_gene_seq, J_gene_seq, true);
            auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACG"}}, /*id=*/7);
            fixture.state.add_event(v_stub);
            fixture.state.mark_chosen(v_stub);
            fixture.state.preset_segment(V_gene_seq, 0, 22, read_run(0, 22));
            fixture.state.add_event(make_deletion(V_gene_seq, Three_prime, 0, 4, /*id=*/8));

            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 6); // the one V discarded above is back
            for (const ScenarioSnapshot &call : next->calls) {
                CHECK(call.safety.at({V_gene_seq, J_gene_seq}) == true);
            }
        }
    }
}

// =======================================================================================
// Row 5 -- neighbour dependence
// =======================================================================================

TEST_CASE("Deletion: neighbour dependence (G5)", "[deletion][iterate]")
{
    SECTION("An unchosen neighbour is not checked and resolves no junction")
    {
        // Nothing to compare against and no span to bound: the deletion enumerates freely and
        // writes neither a safety flag nor a junction bound.
        VDelBare fixture(/*v_five=*/0, /*v_three=*/10, /*min_del=*/0, /*max_del=*/4);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.empty());
            CHECK(call.downstream_bounds.at(VD_ins_seq) == Approx(1.0));
        }
    }

    SECTION("Positive control -- the same geometry with D chosen")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.safety.count({V_gene_seq, D_gene_seq}) == 1);
            CHECK(call.downstream_bounds.at(VD_ins_seq) < 1.0);
        }
    }

    SECTION("With no D in the model a V deletion widens the V->J span instead")
    {
        // The `else if (j_chosen)` arm of initialize_event(). Same event, same side, a
        // different junction -- which is the whole reason B5 reads the junction off the
        // ordering rather than off the gene class.
        VJDel fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (std::size_t i = 0; i != 5; ++i) {
            const int deletions = 4 - static_cast<int>(i);
            INFO("hand-off " << i);
            CHECK(next->calls[i].three_prime(V_gene_seq) == 10 - deletions);
            CHECK(next->calls[i].safety.at({V_gene_seq, J_gene_seq}) == true);
            CHECK(next->calls[i].downstream_bounds.at(VJ_ins_seq)
                  == Approx(junction_bound_at(3 + deletions, VJDel::kContributors)));
            CHECK(next->calls[i].downstream_bounds.at(VD_ins_seq) == Approx(1.0));
        }
    }
}

// =======================================================================================
// Rows 1 and 7 of the per-event list -- the mismatch list
// =======================================================================================

TEST_CASE("Deletion: a positive deletion trims the mismatch list to a contiguous subrange",
          "[deletion][iterate]")
{
    // The list is ordered, so trimming is one walk from whichever end the deletion eats into.
    // The direction is the axis B5 generalises, and it is the opposite of what the arm's name
    // suggests: a 3' deletion keeps the *prefix*, a 5' deletion keeps the *suffix*.

    SECTION("A 3' deletion keeps the prefix at or below the new 3' offset")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{3, 7, 9});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        CHECK(next->calls[0].mismatches.at(V_gene_seq) == std::vector<std::size_t>{3});       // 3' at 6
        CHECK(next->calls[1].mismatches.at(V_gene_seq) == std::vector<std::size_t>{3, 7});    // 3' at 7
        CHECK(next->calls[2].mismatches.at(V_gene_seq) == std::vector<std::size_t>{3, 7});    // 3' at 8
        CHECK(next->calls[3].mismatches.at(V_gene_seq) == std::vector<std::size_t>{3, 7, 9}); // 3' at 9
        CHECK(next->calls[4].mismatches.at(V_gene_seq) == std::vector<std::size_t>{3, 7, 9}); // 3' at 10
    }

    SECTION("A 5' deletion keeps the suffix at or above the new 5' offset")
    {
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*d_mis=*/{14, 16, 18});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        CHECK(next->calls[0].mismatches.at(D_gene_seq) == std::vector<std::size_t>{18});         // 5' at 18
        CHECK(next->calls[1].mismatches.at(D_gene_seq) == std::vector<std::size_t>{18});         // 5' at 17
        CHECK(next->calls[2].mismatches.at(D_gene_seq) == std::vector<std::size_t>{16, 18});     // 5' at 16
        CHECK(next->calls[3].mismatches.at(D_gene_seq) == std::vector<std::size_t>{16, 18});     // 5' at 15
        CHECK(next->calls[4].mismatches.at(D_gene_seq) == std::vector<std::size_t>{14, 16, 18}); // 5' at 14
    }

    SECTION("D 3' keeps the prefix")
    {
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*d_mis=*/{10, 12, 14});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        CHECK(next->calls[0].mismatches.at(D_gene_seq) == std::vector<std::size_t>{10});
        CHECK(next->calls[4].mismatches.at(D_gene_seq) == std::vector<std::size_t>{10, 12, 14});
    }

    SECTION("J 5' keeps the suffix")
    {
        JDel fixture(/*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/4,
                     /*j_mis=*/{18, 20, 22});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        CHECK(next->calls[0].mismatches.at(J_gene_seq) == std::vector<std::size_t>{22});
        CHECK(next->calls[4].mismatches.at(J_gene_seq) == std::vector<std::size_t>{18, 20, 22});
    }

    SECTION("An empty incoming list stays empty")
    {
        VDel fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.mismatches.at(V_gene_seq).empty());
        }
    }

    SECTION("A deletion that outruns every mismatch clears the list")
    {
        // The `end_reached` walk, which each arm spells out separately.
        {
            VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/4, /*max_del=*/4, /*v_five=*/0,
                         /*v_mis=*/{8, 9, 10});
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1);
            CHECK(next->calls.front().mismatches.at(V_gene_seq).empty());
        }
        {
            D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/4,
                          /*max_del=*/4, /*d_mis=*/{14, 15, 16});
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1);
            CHECK(next->calls.front().mismatches.at(D_gene_seq).empty());
        }
        {
            D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/4,
                          /*max_del=*/4, /*d_mis=*/{12, 13, 14});
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1);
            CHECK(next->calls.front().mismatches.at(D_gene_seq).empty());
        }
        {
            JDel fixture(/*d_three=*/14, /*j_five=*/18, /*min_del=*/4, /*max_del=*/4,
                         /*j_mis=*/{18, 19, 20});
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1);
            CHECK(next->calls.front().mismatches.at(J_gene_seq).empty());
        }
    }
}

TEST_CASE("Deletion: a negative deletion appends the reverse-complemented template (palindrome)",
          "[deletion][iterate]")
{
    // A negative deletion is a palindromic insertion: the k nucleotides nearest the trimmed end
    // are reversed, complemented, and put back on the far side of it. The new nucleotides are
    // then compared against the read and any disagreement joins the mismatch list.

    SECTION("A 3' deletion appends at the 3' end")
    {
        VDelBare fixture(/*v_five=*/0, /*v_three=*/10, /*min_del=*/-3, /*max_del=*/-1);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 3);

        // Enumerated in decreasing deletion count, so -1 comes first.
        CHECK(next->calls[0].three_prime(V_gene_seq) == 11);
        CHECK(next->calls[0].sequences.at(V_gene_seq) == read_run(0, 10) + "C");
        CHECK(next->calls[2].three_prime(V_gene_seq) == 13);
        CHECK(next->calls[2].sequences.at(V_gene_seq) == read_run(0, 10) + "CGT");
    }

    SECTION("The appended nucleotides are scored against the read")
    {
        // read[11..13] is "TAC"; the palindrome is "CGT", so all three disagree.
        VDelBare fixture(/*v_five=*/0, /*v_three=*/10, /*min_del=*/-3, /*max_del=*/-3);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().mismatches.at(V_gene_seq) == std::vector<std::size_t>{11, 12, 13});
    }

    SECTION("D 3' appends at the 3' end as well")
    {
        // The fourth arm's copy of the palindrome branch, and the only one no other section
        // reaches. D's template starts at an odd read position here, so the reverse complement
        // disagrees with the read and the new positions join the mismatch list.
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/-3, /*max_del=*/0);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 4);

        CHECK(next->calls[0].three_prime(D_gene_seq) == 14); // no deletion
        CHECK(next->calls[1].three_prime(D_gene_seq) == 15);
        CHECK(next->calls[1].sequences.at(D_gene_seq) == read_run(10, 14) + "C");
        CHECK(next->calls[1].mismatches.at(D_gene_seq) == std::vector<std::size_t>{15});
        CHECK(next->calls[3].three_prime(D_gene_seq) == 17);
        CHECK(next->calls[3].mismatches.at(D_gene_seq) == std::vector<std::size_t>{15, 16, 17});
    }

    SECTION("A 5' deletion prepends at the 5' end")
    {
        // The range runs up to zero: a junction has to be able to close to nothing for the
        // widest palindrome's distance to be inside the fold's range at all.
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/-3, /*max_del=*/0);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 4);

        CHECK(next->calls[0].five_prime(D_gene_seq) == 14); // no deletion
        CHECK(next->calls[1].five_prime(D_gene_seq) == 13);
        CHECK(next->calls[1].sequences.at(D_gene_seq) == "C" + read_run(14, 18));
        CHECK(next->calls[3].five_prime(D_gene_seq) == 11);
        CHECK(next->calls[3].sequences.at(D_gene_seq) == "TAC" + read_run(14, 18));
    }

    SECTION("The 5' arms re-sort, so a new low mismatch lands before the ones already there")
    {
        // The prepended positions are all *below* everything in the incoming list, so without
        // the sort() they would sit at the end of it. D's template starts at an odd read
        // position here, which is what makes its palindrome disagree with the read at all.
        D5Del fixture(/*v_three=*/10, /*d_five=*/15, /*d_three=*/19, /*min_del=*/-3, /*max_del=*/0,
                      /*d_mis=*/{19});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 4);

        CHECK(next->calls[0].mismatches.at(D_gene_seq) == std::vector<std::size_t>{19});
        CHECK(next->calls[1].mismatches.at(D_gene_seq) == std::vector<std::size_t>{14, 19});
        CHECK(next->calls[3].mismatches.at(D_gene_seq) == std::vector<std::size_t>{12, 13, 14, 19});
    }

    SECTION("The J arm re-sorts likewise")
    {
        JDel fixture(/*d_three=*/15, /*j_five=*/19, /*min_del=*/-3, /*max_del=*/0,
                     /*j_mis=*/{25});
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 4);

        CHECK(next->calls[0].mismatches.at(J_gene_seq) == std::vector<std::size_t>{25});
        CHECK(next->calls[1].mismatches.at(J_gene_seq) == std::vector<std::size_t>{18, 25});
        CHECK(next->calls[3].mismatches.at(J_gene_seq) == std::vector<std::size_t>{16, 17, 18, 25});
    }

    SECTION("A 3' palindrome may not run off the end of the read, in the D arm either")
    {
        // The D 3' arm's own copy of the two guards the V arm has. No J chosen, so the fold has
        // nothing to say about the negative gaps these geometries leave and the guards under
        // test are what decides.
        {
            D3DelBare fixture(/*d_five=*/21, /*d_three=*/25, /*min_del=*/-3, /*max_del=*/-2);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1); // -3 would reach offset 28, -2 reaches 27
            CHECK(next->calls.front().three_prime(D_gene_seq) == 27);
        }
        {
            // A short template, so that the palindrome outgrows it well before it reaches the
            // end of the read -- otherwise the guard above catches the realization first and
            // this one is never reached.
            D3DelBare fixture(/*d_five=*/20, /*d_three=*/21, /*min_del=*/-3, /*max_del=*/-2);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1); // 3 > 2 nucleotides of template
            CHECK(next->calls.front().three_prime(D_gene_seq) == 23);
        }
    }

    SECTION("A 5' palindrome may not run off the front of the read")
    {
        // The behaviour, not the branch that produces it. The arm carries an explicit
        // `d_5_new_offset >= 0` test, but it is dead: the read-end guard above it compares the
        // *signed* offset against `int_sequence.size()`, so a negative offset wraps to a huge
        // unsigned value and is rejected there first. Recorded in the plan; what is pinned here
        // is that such a realization is discarded, whichever line does it.
        //
        // No V chosen: the geometry needs a negative junction length, which a fold holding only
        // this deletion and an insertion cannot represent, so the junction guard would
        // otherwise discard the realization before either branch ran.
        D5DelBare fixture(/*d_five=*/2, /*d_three=*/6, /*min_del=*/-3, /*max_del=*/-2);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 1); // -3 would put the 5' end at -1, -2 puts it at 0
        CHECK(next->calls.front().five_prime(D_gene_seq) == 0);
    }

    SECTION("A 5' palindrome may not be longer than the template either")
    {
        {
            // The range reaches up to zero so that every distance the arm asks for is one the
            // fold can answer; the only realization the junction does not account for is the
            // one the guard under test rejects.
            D5Del fixture(/*v_three=*/10, /*d_five=*/20, /*d_three=*/24, /*min_del=*/-6,
                          /*max_del=*/0);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 6); // -6 needs 6 nucleotides of a 5-nucleotide template
            CHECK(next->calls.back().five_prime(D_gene_seq) == 15);
        }
        {
            JDelBare fixture(/*j_five=*/20, /*j_three=*/24, /*min_del=*/-6, /*max_del=*/-5);
            const auto next = call_iterate_recording(fixture.deletion, fixture.state);
            REQUIRE(next->call_count() == 1);
            CHECK(next->calls.front().five_prime(J_gene_seq) == 15);
        }
    }

    SECTION("A palindrome that agrees with the read adds no mismatch")
    {
        // The positive control for the two sections above: the same event, one read position
        // over. ACGT is its own reverse complement, so an even-offset template reproduces the
        // read exactly.
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/-3, /*max_del=*/0);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 4);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.mismatches.at(D_gene_seq).empty());
        }
    }
}

// =======================================================================================
// Row 6 -- downstream bounds
// =======================================================================================

TEST_CASE("Deletion: the junction-length bound (G5/S4)", "[deletion][iterate]")
{
    SECTION("Each arm writes the bound for the flank it widens, keyed by that junction")
    {
        VDel fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        for (std::size_t i = 0; i != 5; ++i) {
            const int deletions = 4 - static_cast<int>(i);
            // The distance the arm looks up: the gap left between the two facing offsets.
            const int distance = 14 - (10 - deletions) - 1;
            INFO("hand-off " << i << ", distance " << distance);
            CHECK(next->calls[i].downstream_bounds.at(VD_ins_seq)
                  == Approx(junction_bound_at(distance, VDel::kContributors)));
            // ...and nothing is written for the junction on the other side.
            CHECK(next->calls[i].downstream_bounds.at(DJ_ins_seq) == Approx(1.0));
        }
    }

    SECTION("The D 5' arm measures against the V end it faces")
    {
        D5Del fixture;
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);
        for (std::size_t i = 0; i != 5; ++i) {
            const int deletions = 4 - static_cast<int>(i);
            const int distance = (14 + deletions) - 10 - 1;
            INFO("hand-off " << i);
            CHECK(next->calls[i].downstream_bounds.at(VD_ins_seq)
                  == Approx(junction_bound_at(distance, D5Del::kContributors)));
        }
    }

    SECTION("The DJ arms measure against the J and D ends respectively")
    {
        D3Del d3;
        const auto from_d = call_iterate_recording(d3.deletion, d3.state);
        REQUIRE(from_d->call_count() == 5);
        CHECK(from_d->calls[0].downstream_bounds.at(DJ_ins_seq)
              == Approx(junction_bound_at(18 - (14 - 4) - 1, D3Del::kContributors)));

        JDel j;
        const auto from_j = call_iterate_recording(j.deletion, j.state);
        REQUIRE(from_j->call_count() == 5);
        CHECK(from_j->calls[0].downstream_bounds.at(DJ_ins_seq)
              == Approx(junction_bound_at((18 + 4) - 14 - 1, JDel::kContributors)));
    }

    SECTION("A distance no completion can reach discards the realization")
    {
        // One knob against the section below: the insertion downstream can supply at most 4
        // nucleotides, so the two widest deletions ask for a gap nothing can fill.
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{}, /*ins_max=*/4);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 2);
        CHECK(next->calls[0].three_prime(V_gene_seq) == 9);
        CHECK(next->calls[1].three_prime(V_gene_seq) == 10);
    }

    SECTION("Positive control -- the same geometry with room in the junction")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{}, /*ins_max=*/20);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        CHECK(next->call_count() == 5);
    }

    SECTION("The other three arms discard on it too")
    {
        // Each arm has its own copy of the `best_for()` lookup and its own `continue`. One knob
        // per pair: the insertion range downstream.
        {
            D5Del narrow(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/0,
                         /*max_del=*/4, /*d_mis=*/{}, /*ins_max=*/4);
            D5Del wide(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/0,
                       /*max_del=*/4, /*d_mis=*/{}, /*ins_max=*/20);
            CHECK(call_iterate_recording(narrow.deletion, narrow.state)->call_count() == 2);
            CHECK(call_iterate_recording(wide.deletion, wide.state)->call_count() == 5);
        }
        {
            D3Del narrow(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/0,
                         /*max_del=*/4, /*d_mis=*/{}, /*ins_max=*/4);
            D3Del wide(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/0,
                       /*max_del=*/4, /*d_mis=*/{}, /*ins_max=*/20);
            CHECK(call_iterate_recording(narrow.deletion, narrow.state)->call_count() == 2);
            CHECK(call_iterate_recording(wide.deletion, wide.state)->call_count() == 5);
        }
        {
            JDel narrow(/*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/4,
                        /*j_mis=*/{}, /*ins_max=*/4);
            JDel wide(/*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*j_mis=*/{}, /*ins_max=*/20);
            CHECK(call_iterate_recording(narrow.deletion, narrow.state)->call_count() == 2);
            CHECK(call_iterate_recording(wide.deletion, wide.state)->call_count() == 5);
        }
    }
}

TEST_CASE("Deletion: what the junction fold composes", "[deletion][iterate][junction]")
{
    // The bound the sections above assert is a product over the events the fold walks, and one
    // of its factors is not a realization at all: Dinucl_markov contributes `p^L` for the
    // segment it fills, where L is the length the segment's *creator* published on that path.
    //
    // That factor's whole behavioural cover used to be a single assertion (plan 6.12). The
    // sections above widened it by accident -- every junction bound they assert carries it --
    // but they all place an Insertion in the fold, so the branch that decides what happens when
    // no creator ran was still untested. These two sections are that branch, and its control.
    //
    // The geometry is deliberately tight: with no Insertion the fold can only *shorten* the
    // junction, so the neighbours have to be adjacent for any distance to be reachable at all.

    const auto vd_bound_with = [](bool with_insertion) {
        IterateTestState state = create_iterate_state(kRead);
        auto deletion = make_deletion(V_gene_seq, Three_prime, 0, 0, /*id=*/0);
        state.preset_safety(V_gene_seq, D_gene_seq, false);
        state.preset_segment(V_gene_seq, 0, 13, read_run(0, 13));
        auto d_stub = make_gene_choice(D_gene, {{"D1", "ACGTA"}}, /*id=*/1);
        state.add_event(d_stub);
        state.mark_chosen(d_stub);
        state.preset_segment(D_gene_seq, 14, 18, read_run(14, 18));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 2, /*id=*/2));
        if (with_insertion) {
            state.add_downstream_event(make_insertion(VD_ins_seq, 0, 4, /*id=*/3));
        }
        state.add_downstream_event(make_dinucl_markov(VD_ins_seq, /*id=*/4));
        for (std::size_t i = 0; i != 64; ++i) {
            state.set_marginal(i, 0.5L);
        }
        const auto next = call_iterate_recording(deletion, state);
        REQUIRE(next->call_count() == 1);
        return next->calls.front().downstream_bounds.at(VD_ins_seq);
    };

    SECTION("A Dinucl_markov whose segment nothing created contributes 1, not p")
    {
        // Two events enumerate realizations here -- this deletion and the D 5' one -- and each
        // carries a flat 0.5. The Dinucl is in the fold (it is asked for a factor) but the
        // junction it fills has no creator on this path, so it must contribute nothing.
        CHECK(vd_bound_with(/*with_insertion=*/false) == Approx(0.25));
    }

    SECTION("Positive control -- with the creator in the fold the insertion's own factor appears")
    {
        // Same geometry, same distance of zero. The extra halving is the Insertion's marginal,
        // not the Dinucl's: at length zero `p^L` is still 1.
        CHECK(vd_bound_with(/*with_insertion=*/true) == Approx(0.125));
    }
}

TEST_CASE("Deletion: the segment's own error bound", "[deletion][iterate]")
{
    SECTION("V publishes the bound for what survives the deletion")
    {
        VDel fixture;
        fixture.state.set_error_rate(0.1);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        for (std::size_t i = 0; i != 5; ++i) {
            const int length = 7 + static_cast<int>(i); // 11 nucleotides less the deletion
            INFO("hand-off " << i);
            CHECK(next->calls[i].downstream_bounds.at(V_gene_seq)
                  == Approx(err_bound(0.1, 0, length)));
        }
    }

    SECTION("Mismatches that survive the deletion are charged")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{3, 7, 9});
        fixture.state.set_error_rate(0.1);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 5);

        CHECK(next->calls[0].downstream_bounds.at(V_gene_seq) == Approx(err_bound(0.1, 1, 6)));
        CHECK(next->calls[4].downstream_bounds.at(V_gene_seq) == Approx(err_bound(0.1, 3, 8)));
    }

    SECTION("D publishes it only once the opposite side has been processed")
    {
        // With both D deletions still to come the surviving mismatch count is not yet known,
        // so the arm declines to bound the segment and writes 1.0. That is a *weaker* bound,
        // not a wrong one, and the source marks the missing computation with a TODO -- pinned
        // as observed rather than as intended (guide section 11).
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*d_mis=*/{14, 16, 18});
        fixture.state.set_error_rate(0.1);
        fixture.state.add_event(make_deletion(D_gene_seq, Three_prime, 0, 4, /*id=*/6));
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.downstream_bounds.at(D_gene_seq) == Approx(1.0));
        }
    }

    SECTION("The D 3' arm declines in the same way")
    {
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*d_mis=*/{10, 12, 14});
        fixture.state.set_error_rate(0.1);
        fixture.state.add_event(make_deletion(D_gene_seq, Five_prime, 0, 4, /*id=*/6));
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        for (const ScenarioSnapshot &call : next->calls) {
            CHECK(call.downstream_bounds.at(D_gene_seq) == Approx(1.0));
        }
    }

    SECTION("Positive control -- with no D deletion on the other side it does publish")
    {
        // An absent opposite-side deletion reads as "already processed", so the same fixture
        // minus that one event bounds the segment for real.
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/0, /*max_del=*/4,
                      /*d_mis=*/{14, 16, 18});
        fixture.state.set_error_rate(0.1);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);

        REQUIRE(next->call_count() == 5);
        CHECK(next->calls[0].downstream_bounds.at(D_gene_seq) == Approx(err_bound(0.1, 1, 0)));
        CHECK(next->calls[4].downstream_bounds.at(D_gene_seq) == Approx(err_bound(0.1, 3, 2)));
    }
}

// =======================================================================================
// Row 7 -- pruning
// =======================================================================================

TEST_CASE("Deletion: pruning", "[deletion][iterate]")
{
    SECTION("Raising the threshold drops realizations from the wide end first")
    {
        // The bound falls as the deletion widens, because the junction it has to fill grows.
        VDel none;
        CHECK(call_iterate_recording(none.deletion, none.state)->call_count() == 5);

        VDel middle;
        middle.state.set_pruning_threshold(0.001);
        CHECK(call_iterate_recording(middle.deletion, middle.state)->call_count() == 3);

        VDel tight;
        tight.state.set_pruning_threshold(0.005);
        CHECK(call_iterate_recording(tight.deletion, tight.state)->call_count() == 1);
    }

    SECTION("What survives is the narrow end, with its own state intact")
    {
        VDel fixture;
        fixture.state.set_pruning_threshold(0.005);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().three_prime(V_gene_seq) == 10); // zero deletions
    }

    SECTION("The D arm prunes once where V and J prune twice")
    {
        // Same threshold, same 0.5 marginals, and a fold one contributor shorter -- so the D
        // arm's bound sits a factor of two higher and two realizations survive where the V arm
        // keeps one. The point of the section is that a *single* check is what stands here:
        // the D arms have no `break` stage at all.
        D5Del fixture;
        fixture.state.set_pruning_threshold(0.005);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        CHECK(next->call_count() == 2);
    }

    SECTION("The D 3' and J arms prune too")
    {
        // Same threshold and the same fold shape as the D 5' section: without these, disabling
        // either arm's check goes unnoticed.
        D3Del d3;
        d3.state.set_pruning_threshold(0.005);
        CHECK(call_iterate_recording(d3.deletion, d3.state)->call_count() == 2);

        JDel j;
        j.state.set_pruning_threshold(0.005);
        CHECK(call_iterate_recording(j.deletion, j.state)->call_count() == 2);

        D3Del d3_open;
        CHECK(call_iterate_recording(d3_open.deletion, d3_open.state)->call_count() == 5);
        JDel j_open;
        CHECK(call_iterate_recording(j_open.deletion, j_open.state)->call_count() == 5);
    }

    SECTION("The V arm's first check stops the enumeration outright")
    {
        // The `break`. It is reached here -- a segment whose surviving mismatches put its
        // error bound below the threshold -- but it cannot be told apart from a `continue` by
        // any observation: the first check's bound is the second's divided by two factors that
        // are both at most 1 (the realization's marginal, and the junction bound, which the
        // arm deliberately holds at 1.0 for this first pass). So the first check can only fire
        // where the second would have fired too, and both bounds fall monotonically along the
        // enumeration order. Recorded in the plan; the section pins the stopping, not the
        // mechanism.
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{3, 7, 9});
        fixture.state.set_error_rate(0.1);
        fixture.state.set_pruning_threshold(0.01);
        CHECK(call_iterate_recording(fixture.deletion, fixture.state)->call_count() == 0);
    }

    SECTION("Positive control -- the same fixture with pruning off")
    {
        VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/0, /*max_del=*/4, /*v_five=*/0,
                     /*v_mis=*/{3, 7, 9});
        fixture.state.set_error_rate(0.1);
        CHECK(call_iterate_recording(fixture.deletion, fixture.state)->call_count() == 5);
    }
}

// =======================================================================================
// Row 8 -- empty versus absent
// =======================================================================================

TEST_CASE("Deletion: a fully deleted D segment is written-and-empty, not absent",
          "[deletion][iterate]")
{
    SECTION("D 5' -- the segment exists with no content and degenerate offsets")
    {
        D5Del fixture(/*v_three=*/10, /*d_five=*/14, /*d_three=*/18, /*min_del=*/5, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);

        REQUIRE(next->calls.front().sequences.count(D_gene_seq) == 1);
        CHECK(next->calls.front().sequences.at(D_gene_seq).empty());
        // The empty-segment convention: three_prime == five_prime - 1.
        CHECK(next->calls.front().five_prime(D_gene_seq) == 19);
        CHECK(next->calls.front().three_prime(D_gene_seq) == 18);
        CHECK(has_constructed_sequence(fixture.state, D_gene_seq));
    }

    SECTION("D 3' -- the same, from the other end")
    {
        D3Del fixture(/*d_five=*/10, /*d_three=*/14, /*j_five=*/18, /*min_del=*/5, /*max_del=*/5);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);

        CHECK(next->calls.front().sequences.at(D_gene_seq).empty());
        CHECK(next->calls.front().five_prime(D_gene_seq) == 10);
        CHECK(next->calls.front().three_prime(D_gene_seq) == 9);
    }

    SECTION("A segment nothing touched is absent from the snapshot entirely")
    {
        VDelBare fixture(/*v_five=*/0, /*v_three=*/10, /*min_del=*/0, /*max_del=*/0);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().sequences.count(D_gene_seq) == 0);
        CHECK(next->calls.front().offsets.count(D_gene_seq) == 0);
    }
}

// =======================================================================================
// Row 9 -- memory layering
// =======================================================================================

TEST_CASE("Deletion: it reads the layer below and writes its own", "[deletion][iterate]")
{
    // Deletion is the event this row exists for: every read is `memory_layer_X - 1` and every
    // write is `memory_layer_X`, so a rewrite that reads its own layer sees the previous
    // realization's state instead of what the upstream event left.
    //
    // The deletion range excludes zero on purpose: with a zero deletion in it the last write
    // would reproduce the incoming value and the two layers would agree by accident.
    VDel fixture(/*v_three=*/10, /*d_five=*/14, /*min_del=*/2, /*max_del=*/4);
    const auto next = call_iterate_recording(fixture.deletion, fixture.state);
    REQUIRE(next->call_count() == 3);

    // Layer 0 still holds what preset_segment() wrote...
    CHECK(int_str_to_nt(*get_constructed_sequence(fixture.state, V_gene_seq, 0)) == read_run(0, 10));
    CHECK(get_seq_offset(fixture.state, V_gene_seq, Three_prime, 0) == 10);

    // ...and layer 1 holds the last realization the event handed off.
    CHECK(int_str_to_nt(*get_constructed_sequence(fixture.state, V_gene_seq, 1)) == read_run(0, 8));
    CHECK(get_seq_offset(fixture.state, V_gene_seq, Three_prime, 1) == 8);

    // The neighbour's offsets were read, never written, so they are still at layer 0.
    CHECK(get_seq_offset(fixture.state, D_gene_seq, Five_prime, 0) == 14);
}

// =======================================================================================
// The read-boundary asymmetry between the four arms
// =======================================================================================

TEST_CASE("Deletion: a model the generic body cannot read is rejected before it runs",
          "[deletion][iterate]")
{
    // Was "a D deletion with no side is rejected at the scenario node": the backstop inside
    // `case D_gene_seq`, and the only one of the body's two `default:` arms that was reachable.
    // B5 removed both -- the side stopped being a branch at all -- and put the question where
    // it belongs, in initialize_event(). The verdict is the same exception; what changed is
    // that it is raised once per model instead of once per scenario node, which is the shape
    // §2.3 asks of B5 for topology it cannot honour.
    //
    // What is *not* checked here is the model-level invariant §2.3's propagation proof rests
    // on -- "every pair adjacent in chosen order is checked by someone". No single event can
    // see it, and `Gene_choice`'s candidate partner list is still the legacy three by name, so
    // no ordering that violates it can be built yet. It stays phase C's, as §2.3 says.

    SECTION("A deletion that names neither end of its segment")
    {
        auto deletion = make_deletion(D_gene_seq, Undefined_side, /*min_del=*/0, /*max_del=*/4,
                                      /*id=*/0);
        IterateTestState state = create_iterate_state(kRead);
        state.preset_segment(D_gene_seq, 14, 18, read_run(14, 18));
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }

        CHECK_THROWS_AS(call_iterate_recording(deletion, state), std::invalid_argument);
    }

    SECTION("A deletion whose segment the ordering leaves out")
    {
        // A D deletion in a model with no D. D_gene_seq is registered -- the legacy six always
        // are -- but it is not in a VJ ordering, so it has no neighbours, no read end and no
        // pairs: every flag initialize_event() resolves would be answered by silence rather
        // than by the model. The four-arm body ran it anyway, on a segment nothing writes.
        auto deletion = make_deletion(D_gene_seq, Five_prime, /*min_del=*/0, /*max_del=*/4,
                                      /*id=*/0);
        IterateTestState state = create_iterate_state(kRead, 1000, 32, vj_seq_type_registry());
        for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }

        CHECK_THROWS_AS(call_iterate_recording(deletion, state), std::invalid_argument);
    }
}

TEST_CASE("Deletion: the J arm does not guard its palindrome against the start of the read",
          "[deletion][iterate]")
{
    // Three of the four arms bound the palindrome against the read before scoring it against
    // it -- V with `v_3_new_offset < sequence.size()`, D 5' with `d_5_new_offset >= 0`, D 3'
    // with `d_3_offset + 1 + i >= 0`. The J arm has neither, so a palindrome that reaches past
    // position 0 indexes int_sequence with a negative offset and the `.at()` throws.
    //
    // Pinned as observed, not as intended: *that* it should be guarded is clear from the other
    // three arms, but *how* is not -- the three of them disagree on whether to reject the
    // realization or to skip the out-of-range comparisons. Recorded in the plan for B5 to
    // settle; until then this section is what stops the behaviour changing unnoticed.

    SECTION("A J palindrome that reaches before position 0 throws")
    {
        JDelBare fixture(/*j_five=*/2, /*j_three=*/27, /*min_del=*/-3, /*max_del=*/-3);
        CHECK_THROWS_AS(call_iterate_recording(fixture.deletion, fixture.state),
                        std::out_of_range);
    }

    SECTION("Positive control -- two read positions further in, the same deletion is fine")
    {
        JDelBare fixture(/*j_five=*/4, /*j_three=*/27, /*min_del=*/-3, /*max_del=*/-3);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        REQUIRE(next->call_count() == 1);
        CHECK(next->calls.front().five_prime(J_gene_seq) == 1);
    }

    SECTION("The V arm, in the mirror-image geometry, discards instead")
    {
        VDelBare fixture(/*v_five=*/0, /*v_three=*/25, /*min_del=*/-3, /*max_del=*/-3);
        const auto next = call_iterate_recording(fixture.deletion, fixture.state);
        CHECK(next->call_count() == 0);
    }
}
