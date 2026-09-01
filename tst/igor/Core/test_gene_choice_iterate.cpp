/*
 * test_gene_choice_iterate.cpp
 *
 *  Characterization tests for Gene_choice::iterate().
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
 * These are **characterization** tests: every expected value was read off the current
 * implementation, not derived from what it ought to compute. Where a pinned value is known
 * to be wrong, the section name says so and points at the plan.
 *
 * The file is organised by the generic pattern each branch belongs to rather than by
 * V/D/J, because V/D/J is the axis B11 removes -- see
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md sections 6.1 and 4.
 */

#include "test_utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace IgorTestUtils;

namespace {

/// Dump a recorder's observations. Used while deriving expected values; kept because the
/// next person to add a section will want it too.
void dump(const std::shared_ptr<RecordingEvent> &rec, const char *label)
{
    UNSCOPED_INFO(label << ": " << rec->call_count() << " call(s)");
    for (std::size_t i = 0; i != rec->calls.size(); ++i) {
        const ScenarioSnapshot &s = rec->calls[i];
        UNSCOPED_INFO("  [" << i << "] proba=" << s.scenario_proba);
        for (const auto &[seq_type, off] : s.offsets) {
            UNSCOPED_INFO("      offsets[" << seq_type << "] = (" << off.first << ", " << off.second
                                           << ")");
        }
        for (const auto &[seq_type, seq] : s.sequences) {
            UNSCOPED_INFO("      seq[" << seq_type << "] = \"" << seq << "\"");
        }
        for (const auto &[safety, value] : s.safety) {
            UNSCOPED_INFO("      safety[" << safety << "] = " << value);
        }
        for (const auto &[seq_type, bound] : s.downstream_bounds) {
            UNSCOPED_INFO("      bound[" << seq_type << "] = " << bound);
        }
    }
}

} // namespace

TEST_CASE("Gene_choice::iterate baseline writes (G4)", "[gene_choice][iterate][baseline]")
{
    SECTION("V: one alignment writes offsets, sequence and mismatches")
    {
        const std::string gene = "ACGTACGTACGT";
        const std::string read = gene + "TTTTTT";

        auto state = create_iterate_state(read);
        auto v_event = make_gene_choice(V_gene, {{"V1", gene}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "V baseline");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        CHECK(s.five_prime(V_gene_seq) == 0);
        CHECK(s.three_prime(V_gene_seq) == static_cast<Seq_Offset>(gene.size()) - 1);
        CHECK(s.sequences.at(V_gene_seq) == gene);
        CHECK(s.mismatches.at(V_gene_seq).empty());
        CHECK(s.scenario_proba == 1.0);
    }

    SECTION("D: one alignment writes offsets from the alignment, not from the read")
    {
        const std::string d_gene = "TTTTAAAA";
        auto state = create_iterate_state("ACGTACGTTTTTAAAAGGGGCCCC");
        auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
        state.set_alignments(D_gene, {create_perfect_alignment("D1", 8, d_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(d_event, state);
        dump(rec, "D baseline");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        CHECK(s.five_prime(D_gene_seq) == 8);
        CHECK(s.three_prime(D_gene_seq) == 15);
        CHECK(s.sequences.at(D_gene_seq) == d_gene);
    }

    SECTION("J: one alignment writes offsets, sequence and mismatches")
    {
        const std::string j_gene = "GGGGCCCC";
        auto state = create_iterate_state("ACGTACGTACGTGGGGCCCC");
        auto j_event = make_gene_choice(J_gene, {{"J1", j_gene}}, 0, /*fixed=*/false);
        state.set_alignments(J_gene, {create_perfect_alignment("J1", 12, j_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(j_event, state);
        dump(rec, "J baseline");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        CHECK(s.five_prime(J_gene_seq) == 12);
        CHECK(s.three_prime(J_gene_seq) == 19);
        CHECK(s.sequences.at(J_gene_seq) == j_gene);
    }
}

TEST_CASE("Gene_choice::iterate template overhangs are trimmed away (G4, feeds B3)",
          "[gene_choice][iterate][baseline][flank]")
{
    // The two branches that clip a genomic template to the part the read can see. They are
    // mirror images and each handles exactly one side -- the side that faces outward:
    // V clips a 5' overhang (negative alignment offset), J clips a 3' overhang (alignment
    // running past the read end). Neither handles the other side.
    //
    // The clipped nucleotides are DISCARDED today. Task B3 of the parent plan turns them
    // into left_flank_seq / right_flank_seq instead, so these two sections are exactly the
    // ones that must change when B3 lands -- and until then they pin what is being lost.

    SECTION("V: a negative alignment offset drops the pre-alignment strip")
    {
        // Genechoice.cpp:253 -- gene_seq = value_str_int.substr(-offset), v_5_off = 0.
        const std::string v_gene = "ACGTACGTACGT"; // 12 nt
        auto state = create_iterate_state("TACGTACGTTTTTTT");
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", -3, v_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "V negative offset");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        // The first three template nucleotides ("ACG") are gone, not stored anywhere.
        CHECK(s.sequences.at(V_gene_seq) == "TACGTACGT");
        CHECK(s.sequences.at(V_gene_seq).size() == v_gene.size() - 3);
        // The 5' offset is clamped to the start of the read rather than going negative.
        CHECK(s.five_prime(V_gene_seq) == 0);
        CHECK(s.three_prime(V_gene_seq) == 8);
        // B3: the strip "ACG" becomes left_flank_seq here, and this section gains an
        // assertion on it rather than losing the ones above.
    }

    SECTION("J: an alignment running past the read end drops the tail")
    {
        // Genechoice.cpp:941 -- gene_seq = value_str_int.substr(0, sequence.size() - offset).
        const std::string j_gene = "GGGGCCCCTTTT"; // 12 nt
        const std::string read = "ACGTACGTACGTAAGGGGCC"; // 20 nt
        auto state = create_iterate_state(read);
        auto j_event = make_gene_choice(J_gene, {{"J1", j_gene}}, 0, /*fixed=*/false);
        state.set_alignments(J_gene, {create_perfect_alignment("J1", 14, j_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(j_event, state);
        dump(rec, "J past read end");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        // Only the 6 nucleotides that fit in the read survive; "CCTTTT" is discarded.
        CHECK(s.sequences.at(J_gene_seq) == "GGGGCC");
        CHECK(s.sequences.at(J_gene_seq).size() == read.size() - 14);
        CHECK(s.five_prime(J_gene_seq) == 14);
        // The 3' offset lands on the last read position, never past it.
        CHECK(s.three_prime(J_gene_seq) == static_cast<Seq_Offset>(read.size()) - 1);
        // B3: the tail "CCTTTT" becomes right_flank_seq here.
    }

    SECTION("V does not clip a 3' overhang, and J does not clip a 5' one")
    {
        // The asymmetry is deliberate in the sense that each gene only ever overhangs the
        // side it faces -- but nothing enforces it, and the generic B11 body will have one
        // clip path rather than two. Pinning it here so that unifying them is a visible
        // change rather than a silent one.
        const std::string v_gene = "ACGTACGTACGT";
        const std::string read = "ACGTACG"; // 7 nt: V runs 5 nt past the end
        auto state = create_iterate_state(read);
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "V past read end (not clipped)");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        // The full 12 nt template is stored even though the read holds only 7, and the 3'
        // offset points past the last read position.
        CHECK(s.sequences.at(V_gene_seq) == v_gene);
        CHECK(s.three_prime(V_gene_seq) == 11);
        CHECK(s.three_prime(V_gene_seq) > static_cast<Seq_Offset>(read.size()) - 1);
    }
}

TEST_CASE("Gene_choice::iterate overlap verdicts (G2/G3)", "[gene_choice][iterate][safety]")
{
    // V occupies read positions [0, 11]. Which verdict the V-D pair gets turns on where
    // D's 5' end sits and how far each end can still travel under its pending deletions:
    //
    //   Infeasible : v_3_off + v_3_max_del  >= d_5_off - d_5_max_del
    //   Safe       : v_3_off + v_3_min_del  <  d_5_off - d_5_min_del
    //   otherwise  : Undetermined
    //
    // *_max_del / *_min_del are Deletion::len_min / len_max, i.e. *negated* deletion
    // counts. That sign convention is what section 2.2 of the plan collapses into
    // reachable().
    //
    // Both deletion events and a VD insertion are registered downstream, so the
    // junction-length map spans a range of lengths rather than {0}. Without that, every
    // non-adjacent geometry is discarded by the junction guard before the safety verdict
    // can matter -- and the section would pass no matter what the safety check did.
    const std::string v_gene = "ACGTACGTACGT"; // 12 nt, 3' end at 11
    const std::string d_gene = "TTTT";
    const std::string read = "ACGTACGTACGTTTTTTTTTTTTTTTTTTTTT";

    struct Fixture {
        std::shared_ptr<Gene_choice> v_event;
        std::shared_ptr<Gene_choice> d_stub;
    };

    auto build = [&](IterateTestState &state, Seq_Offset d_five_prime, int d5_min_del,
                     int d5_max_del, int v3_min_del, int v3_max_del) {
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        auto d_stub = make_gene_choice(D_gene, {{"D1", d_gene}}, 1);

        state.add_event(d_stub);
        state.mark_chosen(d_stub);
        state.preset_segment(D_gene_seq, d_five_prime,
                             d_five_prime + static_cast<Seq_Offset>(d_gene.size()) - 1, d_gene);

        state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, v3_min_del, v3_max_del, 2));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, d5_min_del, d5_max_del, 3));

        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene.size())});
        state.set_marginal(0, 1.0L);
        return v_event;
    };

    SECTION("Infeasible: V 3' still overruns D 5' when both are maximally deleted")
    {
        // v_3_off + v_3_max_del = 11 - 4 = 7 ; d_5_off - d_5_max_del = 3 + 4 = 7 ; 7 >= 7.
        //
        // This section pins the *outcome*, not the overlap check specifically, and cannot
        // do better: the overlap early-out and the junction-length guard reject exactly the
        // same geometries. Writing L for the pre-deletion gap d_5_off - v_3_off - 1,
        //
        //   overlap rejects iff  L <= -max_del_V - max_del_D5 - 1
        //   guard   accepts iff  L >= min_ins - max_del_V - max_del_D5   (and L achievable)
        //
        // -- adjacent, mutually exclusive intervals when min_ins == 0, and the guard is
        // strictly stronger when min_ins > 0. Verified by mutation: deleting either check
        // on its own leaves every section here green. The overlap check is a pure early-out
        // in Gene_choice, worth keeping for the work it skips but carrying no semantics the
        // guard does not already enforce. See plan section 7.6.
        auto state = create_iterate_state(read);
        auto v_event = build(state, /*d_five_prime=*/3, 0, 4, 0, 4);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "VD infeasible");
        REQUIRE(rec->call_count() == 0);
    }

    SECTION("Positive control: one nucleotide more clearance and the same geometry survives")
    {
        // Only d_five_prime changes, 3 -> 4. This is the boundary on both checks at once:
        // the overlap comparison becomes 7 >= 8 (false) and L becomes -8, the smallest
        // achievable gap. Together with the section above it pins the boundary position,
        // though not which of the two checks enforces it.
        auto state = create_iterate_state(read);
        auto v_event = build(state, /*d_five_prime=*/4, 0, 4, 0, 4);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "VD feasible (boundary)");
        REQUIRE(rec->call_count() == 1);
    }

    SECTION("Safe: D 5' cannot travel left, so no deletion can create an overlap")
    {
        // D 5' deletions [0,4] move D's 5' end rightward only, so its leftmost reachable
        // position is 12 -- already clear of V's 3' end at 11.
        auto state = create_iterate_state(read);
        auto v_event = build(state, /*d_five_prime=*/12, 0, 4, 0, 4);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "VD safe");
        REQUIRE(rec->call_count() == 1);
        CHECK(rec->calls.at(0).safety.at(Event_safety::VD_safe) == true);
    }

    SECTION("Undetermined: palindromic D 5' insertions could still reach V")
    {
        // D 5' range [-4,4]: the palindromic end reaches offset 8, left of V's 3' end at
        // 11, so the verdict cannot be settled until the deletion is actually drawn.
        auto state = create_iterate_state(read);
        auto v_event = build(state, /*d_five_prime=*/12, -4, 4, 0, 4);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "VD undetermined");
        REQUIRE(rec->call_count() == 1);
        CHECK(rec->calls.at(0).safety.at(Event_safety::VD_safe) == false);
    }

    SECTION("Counterpart not chosen: no check is performed and the pair is marked unsafe")
    {
        // D exists in the model but has not been chosen, so its offset is unknown and
        // there is nothing to compare against.
        auto state = create_iterate_state(read);
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        state.add_event(make_gene_choice(D_gene, {{"D1", d_gene}}, 1));
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "VD not chosen");
        REQUIRE(rec->call_count() == 1);
        CHECK(rec->calls.at(0).safety.at(Event_safety::VD_safe) == false);
    }

    SECTION("No counterpart in the model at all: the pair is marked safe")
    {
        auto state = create_iterate_state(read);
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        REQUIRE(rec->call_count() == 1);
        CHECK(rec->calls.at(0).safety.at(Event_safety::VD_safe) == true);
        CHECK(rec->calls.at(0).safety.at(Event_safety::VJ_safe) == true);
    }
}

TEST_CASE("Gene_choice::iterate junction-length bound (G5)", "[gene_choice][iterate][junction]")
{
    // The guard at Genechoice.cpp:322: the gap between V's 3' end and the chosen D's 5'
    // end, measured *before* either pending deletion, must be a gap the downstream chain
    // can actually produce. The map is keyed by that pre-deletion gap, which is why
    // Deletion contributes -value_int and Insertion +value_int -- the identity enumerated
    // is  ins = gap + del_V + del_D5.
    //
    // Isolating this from the overlap early-out (section 7.6 of the plan) needs a gap that
    // is *too large*: the overlap check only ever fires on gaps that are too small, so a
    // positive out-of-range gap is rejected by the guard alone.
    const std::string v_gene = "ACGTACGTACGT"; // 12 nt, 3' end at 11
    const std::string d_gene = "TTTT";
    const std::string read = "ACGTACGTACGTTTTTTTTTTTTTTTTTTTTT";

    auto build = [&](IterateTestState &state, Seq_Offset d_five_prime) {
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        auto d_stub = make_gene_choice(D_gene, {{"D1", d_gene}}, 1);
        state.add_event(d_stub);
        state.mark_chosen(d_stub);
        state.preset_segment(D_gene_seq, d_five_prime,
                             d_five_prime + static_cast<Seq_Offset>(d_gene.size()) - 1, d_gene);
        state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 4, 2));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 4, 3));
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene.size())});
        state.set_marginal(0, 1.0L);
        return v_event;
    };

    SECTION("Achievable gap: the scenario survives and the VD bound is set")
    {
        // D at 12 leaves a gap of 0, reachable with zero deletions on both ends.
        auto state = create_iterate_state(read);
        auto v_event = build(state, 12);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "junction achievable");

        REQUIRE(rec->call_count() == 1);
        // The bound written for VD_ins_seq is the best probability the chain can reach at
        // this gap, not 1.0 by default.
        CHECK(rec->calls.at(0).downstream_bounds.count(VD_ins_seq) == 1);
    }

    SECTION("Unachievable gap: the scenario is discarded by the guard alone")
    {
        // D at 13 leaves a gap of 1. Deletions only ever widen the post-deletion gap, so a
        // pre-deletion gap of 1 needs at least one insertion, and no Insertion event is in
        // the chain. The overlap check passes here (7 >= 17 is false), so the guard is the
        // only thing rejecting this.
        auto state = create_iterate_state(read);
        auto v_event = build(state, 13);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "junction unachievable");

        REQUIRE(rec->call_count() == 0);
    }
}

TEST_CASE("Gene_choice::iterate endogenous-mismatch counting (G8)",
          "[gene_choice][iterate][endogenous]")
{
    // Mismatches that survive the maximum remaining deletion cannot be explained away, so
    // they set a floor on the error probability. This TEST_CASE covers the *counting* --
    // which mismatches fall inside the surviving core. The credited length is wrong in all
    // three branches and is covered by the [!shouldfail] cases below.
    const double kRate = 0.1;
    const std::string v_gene = "ACGTACGTACGT"; // 12 nt, 3' end at 11

    // V 3' deletions [0,4] => the surviving core is [0, 7]; a mismatch at 2 is endogenous,
    // one at 9 is not.
    auto bound_with_mismatch_at = [&](std::size_t position) {
        auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
        state.set_error_rate(kRate);
        auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
        state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 4, 2));
        state.set_alignments(
                V_gene, {create_alignment_with_mismatches("V1", 0, v_gene.size(), {position})});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        REQUIRE(rec->call_count() == 1);
        return rec->calls.at(0).downstream_bounds.at(V_gene_seq);
    };

    SECTION("A mismatch inside the core costs an error factor; one outside costs nothing")
    {
        // Asserting the *ratio* rather than either bound keeps this independent of the
        // credited-length defect: the two differ by exactly one endogenous mismatch, so
        // swapping an error-free position for an errored one multiplies by (r/3)/(1-r)
        // whatever the credited length happens to be.
        const double inside = bound_with_mismatch_at(2);
        const double outside = bound_with_mismatch_at(9);
        CHECK_THAT(inside / outside,
                   Catch::Matchers::WithinRel((kRate / 3.0) / (1.0 - kRate), 1e-9));
        CHECK(inside < outside);
    }
}

// ============================================================================
// Known defects
//
// Each of the following asserts the behaviour the code *should* have, and is tagged
// [!shouldfail] because it does not have it yet. Catch2 reports an expected failure as a
// pass, so the suite stays green -- and the moment the defect is fixed the case starts
// passing, which [!shouldfail] turns into a failure. That is the point: the tag has to be
// removed deliberately, so a fix cannot land unnoticed.
//
// Each case is section-free on purpose. [!shouldfail] is evaluated per test-case run, and
// Catch2 re-runs a case once per leaf section, so a case mixing passing and failing
// sections would report the passing ones as unexpected passes.
// ============================================================================

TEST_CASE("DEFECT (plan 7.1): V credits an error-free length larger than its surviving core",
          "[gene_choice][iterate][endogenous][defect][!shouldfail]")
{
    // The core that survives the maximum 3' deletion is [v_5_off, v_3_off + v_3_max_del]
    // = [0, 7], so at most 8 positions exist and, with one of them mismatched, at most 7
    // can be error-free. The code passes gene_seq.size() - v_3_max_del - endo, and
    // v_3_max_del is Deletion::len_min (negative), so this evaluates to 12 + 4 - 1 = 15.
    // The sign is inverted: it should be `+ v_3_max_del`.
    //
    // Consequence: the bound is too small by (1-r)^8, so the branch is pruned more
    // aggressively than the model justifies and scenarios that should contribute can be
    // discarded.
    const double kRate = 0.1;
    const std::string v_gene = "ACGTACGTACGT";

    auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
    state.set_error_rate(kRate);
    auto v_event = make_gene_choice(V_gene, {{"V1", v_gene}}, 0, /*fixed=*/false);
    state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 4, 2));
    state.set_alignments(V_gene,
                         {create_alignment_with_mismatches("V1", 0, v_gene.size(), {2, 9})});
    state.set_marginal(0, 1.0L);

    auto rec = call_iterate_recording(v_event, state);
    REQUIRE(rec->call_count() == 1);

    const double correct = std::pow(kRate / 3.0, 1) * std::pow(1.0 - kRate, 7);
    CHECK_THAT(rec->calls.at(0).downstream_bounds.at(V_gene_seq),
               Catch::Matchers::WithinRel(correct, 1e-9));
}

TEST_CASE("DEFECT (plan 7.1): J credits an error-free length larger than its surviving core",
          "[gene_choice][iterate][endogenous][defect][!shouldfail]")
{
    // Mirror of the V case. Core is [j_5_off - j_5_max_del, j_3_off] = [16, 19], 4
    // positions, one of them mismatched, so at most 3 can be error-free. The code credits
    // gene_seq.size() - j_5_max_del - endo = 8 + 4 - 1 = 11.
    const double kRate = 0.1;
    const std::string j_gene = "GGGGCCCC"; // 8 nt, aligned at 12 in a 20 nt read

    auto state = create_iterate_state("ACGTACGTACGTGGGGCCCC");
    state.set_error_rate(kRate);
    auto j_event = make_gene_choice(J_gene, {{"J1", j_gene}}, 0, /*fixed=*/false);
    state.add_downstream_event(make_deletion(J_gene_seq, Five_prime, 0, 4, 2));
    state.set_alignments(J_gene,
                         {create_alignment_with_mismatches("J1", 12, j_gene.size(), {13, 17})});
    state.set_marginal(0, 1.0L);

    auto rec = call_iterate_recording(j_event, state);
    REQUIRE(rec->call_count() == 1);

    const double correct = std::pow(kRate / 3.0, 1) * std::pow(1.0 - kRate, 3);
    CHECK_THAT(rec->calls.at(0).downstream_bounds.at(J_gene_seq),
               Catch::Matchers::WithinRel(correct, 1e-9));
}

TEST_CASE("DEFECT (plan 7.1): D credits one position fewer than its surviving core spans",
          "[gene_choice][iterate][endogenous][defect][!shouldfail]")
{
    // D is the only branch with deletions pending on both sides, so its core is a genuine
    // intersection: [d_5_off - d_5_max_del, d_3_off + d_3_max_del] = [10, 13] here. The
    // sign is right; what is wrong is that the credited length is the *difference* of the
    // two bounds, 13 - 10 = 3, for a span of 4 inclusive positions.
    //
    // A different defect from V and J, and the opposite direction: the bound comes out too
    // large, so it under-prunes. Harmless for correctness, but it is still an off-by-one.
    const double kRate = 0.1;
    const std::string d_gene = "TTTTTTTT"; // 8 nt, aligned at 8 => [8, 15]

    auto state = create_iterate_state("ACGTACGTTTTATTTTACGT");
    state.set_error_rate(kRate);
    auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
    state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 2, 2));
    state.add_downstream_event(make_deletion(D_gene_seq, Three_prime, 0, 2, 3));
    state.set_alignments(D_gene,
                         {create_alignment_with_mismatches("D1", 8, d_gene.size(), {9, 11})});
    state.set_marginal(0, 1.0L);

    auto rec = call_iterate_recording(d_event, state);
    REQUIRE(rec->call_count() == 1);

    const double correct = std::pow(kRate / 3.0, 1) * std::pow(1.0 - kRate, 3); // 4 - 1
    CHECK_THAT(rec->calls.at(0).downstream_bounds.at(D_gene_seq),
               Catch::Matchers::WithinRel(correct, 1e-9));
}

TEST_CASE("DEFECT: the no_d_align position map places D one nucleotide too far 5'",
          "[gene_choice][iterate][exhaustive][defect][!shouldfail]")
{
    // The two D realization paths disagree on what a junction length means.
    //
    //   alignment path, Genechoice.cpp:322 : L = d_5_off - v_3_off - 1
    //                                        i.e. d_5_off = v_3_off + L + 1
    //   position path,  Genechoice.cpp:556 : d_5_off = v_3_off + L
    //
    // and L comes from vj_length_d_position_proba, which is built from the *same*
    // vd_length_best_proba_map the alignment path's guard consults. So the position path is
    // one short: at L = 0, meaning "no insertions", it places D's 5' end on V's 3' end
    // rather than immediately after it, overlapping V's last nucleotide.
    //
    // Fixture: V over [0,11], J over [16,19], vj_len = 4, reachable only by
    // (vd_len, dj_len) = (0, 0). A zero-length VD junction should put D at 12.
    const std::string d_gene = "TTTT";
    auto state = create_iterate_state("ACGTACGTACGTTTAATTTT");

    auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
    auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACGTACGT"}}, 1);
    auto j_stub = make_gene_choice(J_gene, {{"J1", "TTTT"}}, 2);
    state.add_event(v_stub);
    state.add_event(j_stub);
    state.mark_chosen(v_stub);
    state.mark_chosen(j_stub);
    state.preset_segment(V_gene_seq, 0, 11, "ACGTACGTACGT");
    state.preset_segment(J_gene_seq, 16, 19, "TTTT");
    state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 2, 3));
    state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 2, 4));
    state.add_downstream_event(make_deletion(D_gene_seq, Three_prime, 0, 2, 5));
    state.add_downstream_event(make_deletion(J_gene_seq, Five_prime, 0, 2, 6));
    state.set_alignments(D_gene, {});
    state.set_marginal(0, 1.0L);

    auto rec = call_iterate_recording(d_event, state);
    REQUIRE(rec->call_count() == 1);

    // D should start immediately after V's 3' end, not on it.
    CHECK(rec->calls.at(0).five_prime(D_gene_seq) == 12);
}

TEST_CASE("Gene_choice::iterate pruning", "[gene_choice][iterate][pruning]")
{
    // should_prune(bound) is bound < seq_max_prob * proba_threshold_factor. Every other
    // section runs with the threshold at 0 so that nothing is dropped incidentally; this
    // one turns it on.
    const std::string v_gene_a = "ACGTACGTACGT";
    const std::string v_gene_b = "ACGTACGTACGT";

    SECTION("A realization below the threshold is skipped, one above is not")
    {
        auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
        auto v_event =
                make_gene_choice(V_gene, {{"V1", v_gene_a}, {"V2", v_gene_b}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene_a.size()),
                                      create_perfect_alignment("V2", 0, v_gene_b.size())});
        // With no mismatches and rate 0 every downstream bound is 1, so the upper bound is
        // just the realization probability.
        state.set_marginal(0, 0.9L);
        state.set_marginal(1, 0.1L);
        state.set_pruning_threshold(0.5);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "pruning");

        REQUIRE(rec->call_count() == 1);
        CHECK(rec->calls.at(0).scenario_proba == 0.9);
    }

    SECTION("Positive control: both survive with pruning off")
    {
        auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
        auto v_event =
                make_gene_choice(V_gene, {{"V1", v_gene_a}, {"V2", v_gene_b}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, v_gene_a.size()),
                                      create_perfect_alignment("V2", 0, v_gene_b.size())});
        state.set_marginal(0, 0.9L);
        state.set_marginal(1, 0.1L);

        auto rec = call_iterate_recording(v_event, state);
        REQUIRE(rec->call_count() == 2);
    }
}

TEST_CASE("Gene_choice::iterate exhaustive position fallback (G6)",
          "[gene_choice][iterate][exhaustive]")
{
    // When the aligner returns nothing for a D gene, Gene_choice falls back to trying the
    // template at every plausible position. The fallback lives inside `case D_gene`, so V
    // and J never reach it -- an accident of the switch rather than a stated policy, and
    // exactly what B11 must turn into an explicit per-event switch before the generic body
    // extends it to every gene class.

    SECTION("V with no alignments enumerates nothing")
    {
        // Contrast with the D section below: same empty-alignment situation, opposite
        // outcome, and the only thing that differs is the gene class.
        auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
        auto v_event = make_gene_choice(V_gene, {{"V1", "ACGTACGTACGT"}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        REQUIRE(rec->call_count() == 0);
    }

    SECTION("J with no alignments enumerates nothing")
    {
        auto state = create_iterate_state("ACGTACGTACGTTTTTTTT");
        auto j_event = make_gene_choice(J_gene, {{"J1", "GGGGCCCC"}}, 0, /*fixed=*/false);
        state.set_alignments(J_gene, {});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(j_event, state);
        REQUIRE(rec->call_count() == 0);
    }

    SECTION("D with no alignments and no chosen neighbours slides the template")
    {
        // The `else` branch at Genechoice.cpp:681. With neither V nor J chosen the window
        // runs from just inside the read to the last read position, one nucleotide at a
        // time, and every position is emitted because no junction guard applies.
        const std::string d_gene = "TTTT";
        const std::string read = "ACGTACGTACGTTTTTTTTT"; // 20 nt
        auto state = create_iterate_state(read);
        auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
        state.set_alignments(D_gene, {});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(d_event, state);
        dump(rec, "D sliding");

        // 14 positions, 5' offsets 2 through 15. The window starts at 2 -- "V cannot be
        // absent from the read, at least one nucleotide is present" plus one -- and runs
        // while the D 3' end is left of j_5_min_offset, which with no J chosen is the last
        // read position (19).
        REQUIRE(rec->call_count() == 14);
        CHECK(rec->calls.front().five_prime(D_gene_seq) == 2);
        CHECK(rec->calls.back().five_prime(D_gene_seq) == 15);
        CHECK(rec->calls.back().three_prime(D_gene_seq)
              == static_cast<Seq_Offset>(read.size()) - 2);
        // Positions are contiguous and advance by one.
        for (std::size_t i = 1; i < rec->calls.size(); ++i) {
            CHECK(rec->calls.at(i).five_prime(D_gene_seq)
                  == rec->calls.at(i - 1).five_prime(D_gene_seq) + 1);
        }
        // The template is placed whole at every position: 3' = 5' + len - 1.
        for (const ScenarioSnapshot &snapshot : rec->calls) {
            CHECK(snapshot.three_prime(D_gene_seq)
                  == snapshot.five_prime(D_gene_seq) + static_cast<Seq_Offset>(d_gene.size()) - 1);
            CHECK(snapshot.sequences.at(D_gene_seq) == d_gene);
        }
    }

    SECTION("D with no alignments and both neighbours chosen uses the position map")
    {
        // The other sub-branch, at Genechoice.cpp:546: instead of sliding, positions come
        // from vj_length_d_position_proba, built at :1444 by composing the VD and DJ length
        // maps with the D template length. This is the structure B11 has to generalise --
        // the identity is junction_len = d_gene.size() + vd_len + dj_len, with vd/dj taken
        // from VD_ins_seq and DJ_ins_seq by name.
        //
        //   V occupies [0, 11], J occupies [16, 19], so vj_len = 16 - 11 - 1 = 4.
        //   Deletions of at most 2 per end give VD and DJ length maps over [-4, 0], hence
        //   junction keys 4 + vd + dj over [-4, 4]. Key 4 is reachable only by (0, 0).
        const std::string d_gene = "TTTT";
        const std::string read = "ACGTACGTACGTTTAATTTT"; // 20 nt
        auto state = create_iterate_state(read);

        auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
        auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACGTACGT"}}, 1);
        auto j_stub = make_gene_choice(J_gene, {{"J1", "TTTT"}}, 2);
        state.add_event(v_stub);
        state.add_event(j_stub);
        state.mark_chosen(v_stub);
        state.mark_chosen(j_stub);
        state.preset_segment(V_gene_seq, 0, 11, "ACGTACGTACGT");
        state.preset_segment(J_gene_seq, 16, 19, "TTTT");

        state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 2, 3));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 2, 4));
        state.add_downstream_event(make_deletion(D_gene_seq, Three_prime, 0, 2, 5));
        state.add_downstream_event(make_deletion(J_gene_seq, Five_prime, 0, 2, 6));

        state.set_alignments(D_gene, {});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(d_event, state);
        dump(rec, "D position map");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        // The template is placed whole, 5' to 3'. Where it is placed is a defect and is
        // asserted separately below.
        CHECK(s.three_prime(D_gene_seq) == s.five_prime(D_gene_seq) + 3);
        CHECK(s.sequences.at(D_gene_seq) == d_gene);
    }

    SECTION("Mismatches are recomputed per position against the read")
    {
        // Same fixture. The template TTTT sits over read[11..14] = "TTTA", so exactly one
        // position mismatches -- and the mismatch list is rebuilt from the read at each
        // position rather than carried from an alignment, which is the tier-1 computation
        // section 2.10 of the plan wants lifted out of the scenario loop.
        const std::string d_gene = "TTTT";
        const std::string read = "ACGTACGTACGTTTAATTTT";
        auto state = create_iterate_state(read);

        auto d_event = make_gene_choice(D_gene, {{"D1", d_gene}}, 0, /*fixed=*/false);
        auto v_stub = make_gene_choice(V_gene, {{"V1", "ACGTACGTACGT"}}, 1);
        auto j_stub = make_gene_choice(J_gene, {{"J1", "TTTT"}}, 2);
        state.add_event(v_stub);
        state.add_event(j_stub);
        state.mark_chosen(v_stub);
        state.mark_chosen(j_stub);
        state.preset_segment(V_gene_seq, 0, 11, "ACGTACGTACGT");
        state.preset_segment(J_gene_seq, 16, 19, "TTTT");
        state.add_downstream_event(make_deletion(V_gene_seq, Three_prime, 0, 2, 3));
        state.add_downstream_event(make_deletion(D_gene_seq, Five_prime, 0, 2, 4));
        state.add_downstream_event(make_deletion(D_gene_seq, Three_prime, 0, 2, 5));
        state.add_downstream_event(make_deletion(J_gene_seq, Five_prime, 0, 2, 6));
        state.set_alignments(D_gene, {});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(d_event, state);
        dump(rec, "D position mismatches");

        REQUIRE(rec->call_count() == 1);
        // read[14] is 'A' where the template has 'T'; 11, 12 and 13 all match.
        CHECK(rec->calls.at(0).mismatches.at(D_gene_seq) == std::vector<std::size_t>{14});
    }
}
