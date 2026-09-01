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
