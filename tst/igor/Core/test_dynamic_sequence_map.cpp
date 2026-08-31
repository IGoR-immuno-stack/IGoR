/**
 * @file test_dynamic_sequence_map.cpp
 * @brief Unit tests for DynamicSequenceMap and its ordered traversal.
 *
 * The traversal is what lets Insertion/Dinucl_markov/Deletion find their neighbouring
 * segments without naming a topology, so these tests are the executable form of the
 * B10 contract: "not yet processed" (layer -1) and "actively absent" (written, empty)
 * are both skipped, but they are different states.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/IntStr.h>
#include <igor/Core/SeqTypeRegistry.h>

#include <stdexcept>

namespace {

/// Standard VDJ ordering, frozen and ready to size a map from.
SeqTypeRegistry vdj_registry()
{
    SeqTypeRegistry r;
    r.set_ordered_types({"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
    r.freeze();
    return r;
}

/// Tandem-D ordering from the B7 worked example.
SeqTypeRegistry tandem_d_registry()
{
    SeqTypeRegistry r;
    r.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
                         "D1D2_ins_seq", "D2_gene_seq", "D2J_ins_seq", "J_gene_seq"});
    r.freeze();
    return r;
}

} // namespace

TEST_CASE("DynamicSequenceMap: requires a frozen registry", "[dynamic_sequence_map]")
{
    SeqTypeRegistry unfrozen;
    unfrozen.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    // The map sizes itself from the id space, so that space must not be able to grow
    // underneath it.
    CHECK_THROWS_AS(DynamicSequenceMap<Seq_Offset>(unfrozen), std::logic_error);

    unfrozen.freeze();
    CHECK_NOTHROW(DynamicSequenceMap<Seq_Offset>(unfrozen));
}

TEST_CASE("DynamicSequenceMap: sizes itself from the registry", "[dynamic_sequence_map]")
{
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Seq_Offset> m(r);
    CHECK(m.count() == r.total_count());
    CHECK(m.count() == 5);

    // An explicit smaller count is honoured (planned use: excluding flank types, B3).
    DynamicSequenceMap<Seq_Offset> narrow(r, 1, 3);
    CHECK(narrow.count() == 3);
}

TEST_CASE("DynamicSequenceMap: occupied() separates the three states", "[dynamic_sequence_map][b10]")
{
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Int_Str *> m(r);

    Int_Str present = {0, 1, 2};
    Int_Str absent;   // zero length: written, but nothing there

    const SeqTypeId v = r.id("V_gene_seq");
    const SeqTypeId vd = r.id("VD_ins_seq");
    const SeqTypeId d = r.id("D_gene_seq");

    m.set(v, &present, 0);
    m.set(vd, &absent, 0);
    // d is left untouched

    // not yet processed
    CHECK_FALSE(m.exists(d));
    CHECK_FALSE(m.occupied(d));
    // actively absent -- written, so it exists, but it is not occupied
    CHECK(m.exists(vd));
    CHECK_FALSE(m.occupied(vd));
    // present
    CHECK(m.exists(v));
    CHECK(m.occupied(v));

    // A null pointer counts as absent too.
    m.set(d, nullptr, 0);
    CHECK(m.exists(d));
    CHECK_FALSE(m.occupied(d));
}

TEST_CASE("DynamicSequenceMap: non-segment values are occupied once written",
          "[dynamic_sequence_map][b10]")
{
    // Offsets, probabilities and mismatch lists cannot be "absent": the written/unwritten
    // distinction is the whole story for them.
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Seq_Offset> m(r);
    const SeqTypeId v = r.id("V_gene_seq");

    CHECK_FALSE(m.occupied(v));
    m.set(v, 0, 0);            // zero is a perfectly good offset
    CHECK(m.occupied(v));
}

TEST_CASE("DynamicSequenceMap: traversal over a fully populated ordering",
          "[dynamic_sequence_map][traversal]")
{
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Int_Str *> m(r);

    Int_Str seg = {0};
    for (SeqTypeId id : r.ordering()) {
        m.set(id, &seg, 0);
    }

    const SeqTypeId d = r.id("D_gene_seq");
    CHECK(m.first_occupied_left(d) == r.id("VD_ins_seq"));
    CHECK(m.first_occupied_right(d) == r.id("DJ_ins_seq"));

    // The ends of the ordering report kNoSeqType rather than wrapping.
    CHECK(m.first_occupied_left(r.id("V_gene_seq")) == kNoSeqType);
    CHECK(m.first_occupied_right(r.id("J_gene_seq")) == kNoSeqType);
}

TEST_CASE("DynamicSequenceMap: traversal skips unwritten and absent segments alike",
          "[dynamic_sequence_map][traversal][b10]")
{
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Int_Str *> m(r);

    Int_Str present = {0};
    Int_Str absent;

    m.set(r.id("V_gene_seq"), &present, 0);
    m.set(r.id("VD_ins_seq"), &absent, 0);      // written but empty
    // D_gene_seq deliberately left unwritten
    m.set(r.id("J_gene_seq"), &present, 0);

    // From DJ_ins_seq, walking left passes over an unwritten D and an empty VD to reach V.
    CHECK(m.first_occupied_left(r.id("DJ_ins_seq")) == r.id("V_gene_seq"));
    CHECK(m.first_occupied_right(r.id("DJ_ins_seq")) == r.id("J_gene_seq"));
}

// The scenario the plan's B7 example is written around: a tandem-D model whose second D
// is absent must seed from the first D, with no conditional logic in the event.
TEST_CASE("DynamicSequenceMap: tandem-D falls back to D1 when D2 is absent",
          "[dynamic_sequence_map][traversal][tandem_d][b10]")
{
    const SeqTypeRegistry r = tandem_d_registry();
    DynamicSequenceMap<Int_Str *> m(r);

    Int_Str v_seq = {0, 0};
    Int_Str d1_seq = {2, 2};
    Int_Str empty;

    m.set(r.id("V_gene_seq"), &v_seq, 0);
    m.set(r.id("VD1_ins_seq"), &empty, 0);
    m.set(r.id("D1_gene_seq"), &d1_seq, 0);
    m.set(r.id("D1D2_ins_seq"), &empty, 0);
    m.set(r.id("D2_gene_seq"), &empty, 0);      // the mock/None D2: processed, not there

    SECTION("D2 absent: the D2J insertion seeds from D1")
    {
        CHECK(m.first_occupied_left(r.id("D2J_ins_seq")) == r.id("D1_gene_seq"));
    }

    SECTION("D2 present: the same call reaches D2, no code change")
    {
        Int_Str d2_seq = {4, 4};
        m.set(r.id("D2_gene_seq"), &d2_seq, 0);
        CHECK(m.first_occupied_left(r.id("D2J_ins_seq")) == r.id("D2_gene_seq"));
    }
}

TEST_CASE("DynamicSequenceMap: traversal follows the layer stack", "[dynamic_sequence_map][traversal]")
{
    const SeqTypeRegistry r = vdj_registry();
    DynamicSequenceMap<Int_Str *> m(r);

    Int_Str present = {0};
    Int_Str absent;

    m.set(r.id("V_gene_seq"), &present, 0);
    m.set(r.id("D_gene_seq"), &present, 0);

    // Backtracking must be visible to the traversal: pushing a layer on D and writing an
    // empty segment there makes D absent, and popping restores it.
    CHECK(m.first_occupied_left(r.id("DJ_ins_seq")) == r.id("D_gene_seq"));

    m.request_layer(r.id("D_gene_seq"));
    m.set(r.id("D_gene_seq"), &absent, 1);
    CHECK(m.first_occupied_left(r.id("DJ_ins_seq")) == r.id("V_gene_seq"));

    m.restore_layer(r.id("D_gene_seq"));
    CHECK(m.first_occupied_left(r.id("DJ_ins_seq")) == r.id("D_gene_seq"));
}

TEST_CASE("DynamicSequenceMap: a type outside the ordering has no neighbours",
          "[dynamic_sequence_map][traversal]")
{
    SeqTypeRegistry r;
    r.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    const SeqTypeId orphan = r.register_type("left_flank_seq");
    r.freeze();

    DynamicSequenceMap<Int_Str *> m(r);
    Int_Str seg = {0};
    for (SeqTypeId id : r.ordering()) {
        m.set(id, &seg, 0);
    }
    m.set(orphan, &seg, 0);

    // Registered but unordered: reachable by id, invisible to the traversal.
    CHECK(m.occupied(orphan));
    CHECK(m.first_occupied_left(orphan) == kNoSeqType);
    CHECK(m.first_occupied_right(orphan) == kNoSeqType);
}
