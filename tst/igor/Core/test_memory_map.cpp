/**
 * @file test_memory_map.cpp
 * @brief Tests for the Enum_fast_memory_map layered-storage container in Utils.h.
 *
 * The container backs the per-recursion-level scenario state used throughout
 * Rec_Event::iterate(): each recursion depth writes at its own "memory layer" so
 * that backtracking is a pointer move rather than a copy.
 *
 * Layers must be populated in order: layer n+1 is only writable once layer n has
 * been written. Both write paths uphold that:
 *
 *   request_memory_layer()  grows max_layer and advances memory_layer_ptr[key]
 *                           together, one step at a time. Used from
 *                           Rec_Event::initialize_event().
 *   set_value(k, v, layer)  asserts the ordering, then grows the buffer if the
 *                           layer lies past the current allocation.
 *
 * The growth in set_value() was missing: the constructor allocates only layer 0,
 * and a map driven purely by set_value() walking 0,1,2,... has no other growth
 * path, so the write landed past the end of value_ptr_arr.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/Utils.h>

TEST_CASE("Enum_fast_memory_map: set_value grows the buffer along a layer sequence",
          "[utils][memory_map]") {
    // The constructor allocates layer 0 only. Walking upwards through set_value()
    // alone must keep working past that point.
    Enum_fast_memory_map<Seq_type, int> map(6);

    for (int layer = 0; layer <= 3; ++layer) {
        REQUIRE_NOTHROW(map.set_value(V_gene_seq, 100 + layer, layer));
    }
    REQUIRE(map.at(V_gene_seq, 3) == 103);
}

TEST_CASE("Enum_fast_memory_map: growing for one key preserves the others",
          "[utils][memory_map]") {
    Enum_fast_memory_map<Seq_type, int> map(6);

    map.set_value(V_gene_seq, 10, 0);
    map.set_value(D_gene_seq, 20, 0);
    map.set_value(J_gene_seq, 30, 0);

    // Reallocating for V must copy every key's existing layers across.
    map.set_value(V_gene_seq, 11, 1);
    map.set_value(V_gene_seq, 12, 2);

    REQUIRE(map.at(V_gene_seq, 2) == 12);
    REQUIRE(map.at(D_gene_seq, 0) == 20);
    REQUIRE(map.at(J_gene_seq, 0) == 30);
}

TEST_CASE("Enum_fast_memory_map: set_value across every key at every layer",
          "[utils][memory_map]") {
    // Exercises the reallocation path for each key in turn, which is where an
    // undersized buffer shows up as heap corruption rather than a wrong value.
    Enum_fast_memory_map<Seq_type, int> map(6);
    const Seq_type keys[] = {V_gene_seq, VD_ins_seq, D_gene_seq,
                             DJ_ins_seq, J_gene_seq, VJ_ins_seq};

    for (int layer = 0; layer < 5; ++layer) {
        for (Seq_type k : keys) {
            REQUIRE_NOTHROW(map.set_value(k, layer * 100 + static_cast<int>(k), layer));
        }
    }
    for (Seq_type k : keys) {
        REQUIRE(map.at(k, 4) == 400 + static_cast<int>(k));
    }
}

TEST_CASE("Enum_fast_memory_map: request_memory_layer and set_value agree on ordering",
          "[utils][memory_map]") {
    // The route taken in production: initialize_event() reserves the layers via
    // request_memory_layer(), iterate() then writes them with set_value().
    Enum_fast_memory_map<Seq_type, int> map(6);

    map.request_memory_layer(V_gene_seq);  // layer 0, as Gene_choice would take
    const int gene_choice_layer = map.get_current_memory_layer(V_gene_seq);
    map.request_memory_layer(V_gene_seq);  // layer 1, as Deletion would take
    const int deletion_layer = map.get_current_memory_layer(V_gene_seq);

    REQUIRE(gene_choice_layer == 0);
    REQUIRE(deletion_layer == 1);

    REQUIRE_NOTHROW(map.set_value(V_gene_seq, 7, gene_choice_layer));
    REQUIRE_NOTHROW(map.set_value(V_gene_seq, 8, deletion_layer));
    REQUIRE(map.at(V_gene_seq, deletion_layer) == 8);
    REQUIRE(map.at(V_gene_seq, gene_choice_layer) == 7);
}

TEST_CASE("Enum_fast_memory_map: a fresh key is only writable at layer 0",
          "[utils][memory_map]") {
    // Pins the ordering invariant itself. Writing a never-written key at a layer
    // above 0 means the caller is using a layer index belonging to another map's
    // sequence -- the defect that motivated this file. The assertion catching it
    // is compiled out in release builds, so assert only on what holds in both.
    Enum_fast_memory_map<Seq_type, int> map(6);

    REQUIRE(map.get_current_memory_layer(V_gene_seq) == -1);
    REQUIRE_NOTHROW(map.set_value(V_gene_seq, 1, 0));
    REQUIRE(map.get_current_memory_layer(V_gene_seq) == 0);

    // A second key is still untouched, and likewise starts at layer 0.
    REQUIRE(map.get_current_memory_layer(D_gene_seq) == -1);
    REQUIRE_NOTHROW(map.set_value(D_gene_seq, 2, 0));
    REQUIRE(map.at(D_gene_seq, 0) == 2);
}
