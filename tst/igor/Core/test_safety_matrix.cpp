/*
 * test_safety_matrix.cpp
 *
 *  Unit tests for the pairwise overlap-safety matrix (plan step S5, section 2.3).
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

#include <igor/Core/JunctionGeometry.h>
#include <igor/Core/SafetyMatrix.h>

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

namespace {

SeqTypeRegistry vdj_registry()
{
    SeqTypeRegistry registry;
    registry.register_legacy_seq_types();
    registry.set_ordered_types({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
    registry.freeze();
    return registry;
}

/// A VJ model: the ordering puts VJ_ins_seq at position 1, while its legacy id is 5 and
/// J_gene_seq's is 4. Ordering position and id disagree here, which is the whole reason the
/// matrix is keyed by the former.
SeqTypeRegistry vj_registry()
{
    SeqTypeRegistry registry;
    registry.register_legacy_seq_types();
    registry.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    registry.freeze();
    return registry;
}

} // namespace

TEST_CASE("SafetyMatrix names a pair by ordering position, not by seq_type id",
          "[safety_matrix]")
{
    SECTION("The 5'-most member is the row, whichever way round the pair is given")
    {
        const SeqTypeRegistry registry = vdj_registry();
        SafetyMatrix matrix(registry);

        const SafetyCell vd = matrix.cell(V_gene_seq, D_gene_seq);
        CHECK(vd.row == 0);      // V
        CHECK(vd.column == 2);   // D
        CHECK(matrix.cell(D_gene_seq, V_gene_seq) == vd);
    }

    SECTION("A VJ model's positions are not its ids")
    {
        // J_gene_seq is id 4 and VJ_ins_seq id 5, but in this ordering J sits *after* the
        // junction. Keying the bitmask by id would put J's column left of the insertion's and
        // make "the rest of the row" mean the wrong set of segments.
        const SeqTypeRegistry registry = vj_registry();
        SafetyMatrix matrix(registry);

        const SafetyCell vj = matrix.cell(V_gene_seq, J_gene_seq);
        CHECK(vj.row == 0);
        CHECK(vj.column == 2);
        CHECK(matrix.cell(V_gene_seq, VJ_ins_seq).column == 1);
    }

    SECTION("A segment the ordering leaves out is not half of any pair")
    {
        const SeqTypeRegistry registry = vj_registry();
        SafetyMatrix matrix(registry);

        CHECK(matrix.addresses(static_cast<SeqTypeId>(V_gene_seq)));
        CHECK_FALSE(matrix.addresses(static_cast<SeqTypeId>(D_gene_seq)));
        CHECK_THROWS_AS(matrix.cell(V_gene_seq, D_gene_seq), std::out_of_range);
    }

    SECTION("A segment does not overlap itself")
    {
        const SeqTypeRegistry registry = vdj_registry();
        SafetyMatrix matrix(registry);
        CHECK_THROWS_AS(matrix.cell(V_gene_seq, V_gene_seq), std::invalid_argument);
    }

    SECTION("An ordering wider than a row word is rejected, not truncated")
    {
        // The documented ceiling of the bitmask form. Tandem D needs 7 positions; going past
        // 32 is a modelling change that has to move the container, and silently aliasing
        // columns 32 and 0 would make two unrelated pairs the same flag.
        SeqTypeRegistry registry;
        std::vector<Seq_type_String> many;
        for (int i = 0; i != 33; ++i) {
            many.push_back("seg_" + std::to_string(i));
        }
        registry.set_ordered_types(many);
        registry.freeze();
        CHECK_THROWS_AS(SafetyMatrix(registry), std::length_error);
    }
}

TEST_CASE("SafetyMatrix: a row word is the whole state of that row", "[safety_matrix]")
{
    const SeqTypeRegistry registry = vdj_registry();
    SafetyMatrix matrix(registry);
    const SafetyCell vd = matrix.cell(V_gene_seq, D_gene_seq);
    const SafetyCell vj = matrix.cell(V_gene_seq, J_gene_seq);
    const SafetyCell dj = matrix.cell(D_gene_seq, J_gene_seq);

    SECTION("An untouched row has no verdict at all")
    {
        CHECK_FALSE(matrix.exists(vd));
        CHECK(matrix.claimed_layer(vd) == -1);
    }

    SECTION("Two cells of one row written at one layer both stand")
    {
        // The case that decides whether the two cells can share a layer, which is what lets a
        // Gene_choice with two neighbours on the same side claim once. A plain write would
        // make the second erase the first.
        matrix.set(vd, false, 0);
        matrix.set(vj, true, 0);
        CHECK(matrix.get(vd, 0) == false);
        CHECK(matrix.get(vj, 0) == true);
    }

    SECTION("A cell nobody writes at this depth keeps what the depth below left")
    {
        matrix.set(vd, true, 0);
        matrix.request_layer(vd);
        matrix.set(vj, false, 1);

        CHECK(matrix.get(vd, 1) == true);   // carried forward, not reset
        CHECK(matrix.get(vd, 0) == true);   // and the layer below is intact
        CHECK(matrix.get(vj, 1) == false);
    }

    SECTION("Rows are independent")
    {
        matrix.set(vd, true, 0);
        CHECK_FALSE(matrix.exists(dj));
    }

    SECTION("A claimed but unwritten layer is not readable")
    {
        // The same discipline the per-slot container enforced (section 7.9): requesting is a
        // promise to write, and a reader of `layer` must not be served a default because the
        // event that owns it handed off without writing.
        matrix.set(vd, true, 0);
        matrix.request_layer(vd);
        CHECK_THROWS_AS(matrix.get(vd, 1), std::out_of_range);
    }
}

TEST_CASE("SafetyMatrix: establishing a near pair marks the rest of the row",
          "[safety_matrix]")
{
    const SeqTypeRegistry registry = vdj_registry();
    SafetyMatrix matrix(registry);
    const SafetyCell vd = matrix.cell(V_gene_seq, D_gene_seq);
    const SafetyCell vj = matrix.cell(V_gene_seq, J_gene_seq);
    const SafetyCell dj = matrix.cell(D_gene_seq, J_gene_seq);

    SECTION("A safe verdict propagates 3', so one check answers for every further partner")
    {
        // Section 2.3's corollary, and the reason the work per event is O(1) per side rather
        // than one comparison against every segment further along.
        matrix.set(vd, true, 0);
        CHECK(matrix.get(vj, 0) == true);
    }

    SECTION("It does not propagate 5', nor into another row")
    {
        matrix.set(dj, true, 0);
        CHECK_FALSE(matrix.exists(vd));
        CHECK_FALSE(matrix.exists(vj));
    }

    SECTION("An unsafe verdict marks its own column and nothing else")
    {
        // "The deciding event must look" says nothing about pairs further away, so clearing
        // one column must not clear the corollary another check established.
        matrix.set(vd, true, 0);
        matrix.set(vd, false, 0);
        CHECK(matrix.get(vd, 0) == false);
        CHECK(matrix.get(vj, 0) == true);
    }
}

TEST_CASE("SafetyMatrix: propagation moves where a bad scenario dies, not whether",
          "[safety_matrix]")
{
    // Section 2.3's soundness argument, checked against the geometry predicate rather than
    // restated. Propagation can mark (A, C) safe when the completed scenario violates it. The
    // claim is that such a scenario violates (B, C) too -- and (B, C), being adjacent in
    // chosen order, is always checked for real.
    //
    // The offsets are pinned (no pending modifier), so check_overlap reduces to the raw
    // comparison a completed scenario is judged by, and `gap` is 0 as at every site today.
    const auto pinned = [](Seq_Offset at) { return JunctionGeometry::OffsetInterval{at, at}; };
    const auto separated = [&](Seq_Offset left_three, Seq_Offset right_five) {
        return JunctionGeometry::check_overlap(pinned(left_three), pinned(right_five), 0)
               != JunctionGeometry::Overlap::Infeasible;
    };

    const SeqTypeRegistry registry = vdj_registry();
    SafetyMatrix matrix(registry);
    const SafetyCell vd = matrix.cell(V_gene_seq, D_gene_seq);
    const SafetyCell vj = matrix.cell(V_gene_seq, J_gene_seq);

    SECTION("An empty middle segment is the tight case, and it holds")
    {
        // B10's degenerate convention: a segment that produced nothing has
        // three_prime == five_prime - 1, which is the one configuration where B.5' <= B.3' + 1
        // is an equality. If the argument survives anywhere it is here, and section 2.3 asks
        // for it by name rather than by comment.
        const Seq_Offset b_five = 10;
        const Seq_Offset b_three = b_five - 1; // B is empty

        for (Seq_Offset a_three = 0; a_three != b_five; ++a_three) {
            REQUIRE(separated(a_three, b_five)); // (A, B) is established safe
            matrix.set(vd, true, 0);
            REQUIRE(matrix.get(vj, 0) == true);  // ...so (A, C) reads safe by propagation

            for (Seq_Offset c_five = 0; c_five != 24; ++c_five) {
                if (separated(a_three, c_five)) {
                    continue; // the propagated verdict happens to be right; nothing to prove
                }
                INFO("A.3' = " << a_three << ", B = [" << b_five << ", " << b_three << "], C.5' = "
                               << c_five);
                CHECK_FALSE(separated(b_three, c_five)); // the (B, C) check discards it
            }
        }
    }

    SECTION("...and so does a middle segment with a body")
    {
        const Seq_Offset b_five = 10;
        for (Seq_Offset b_three = b_five; b_three != b_five + 8; ++b_three) {
            for (Seq_Offset a_three = 0; a_three != b_five; ++a_three) {
                REQUIRE(separated(a_three, b_five));
                for (Seq_Offset c_five = 0; c_five != 32; ++c_five) {
                    if (separated(a_three, c_five)) {
                        continue;
                    }
                    INFO("A.3' = " << a_three << ", B = [" << b_five << ", " << b_three
                                   << "], C.5' = " << c_five);
                    CHECK_FALSE(separated(b_three, c_five));
                }
            }
        }
    }
}
