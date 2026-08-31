/**
 * @file test_layered_array.cpp
 * @brief Unit tests for LayeredArray, the runtime-sized layered container.
 *
 * Covers the invariants that Enum_fast_memory_map documented only by behaviour:
 * the unwritten (-1) state, per-key independence of layers, growth, and the
 * separation of "read at a layer" from "rewind to a layer".
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <igor/Core/LayeredArray.h>
#include <igor/Core/LayeredMap.h>
#include <igor/Core/Utils.h>

#include <stdexcept>
#include <vector>

TEST_CASE("LayeredArray: unwritten keys are distinct from written-but-empty", "[layered_array]")
{
    LayeredArray<int> a(4);

    for (std::size_t k = 0; k != 4; ++k) {
        CHECK_FALSE(a.exists(k));
        CHECK(a.current_layer(k) == -1);
        CHECK_THROWS_AS(a.get(k), std::out_of_range);
    }

    // A key written with a "empty" sentinel value still exists -- this is the
    // distinction the ordered traversal depends on (actively absent vs not yet processed).
    a.set(2, 0, 0);
    CHECK(a.exists(2));
    CHECK(a.get(2) == 0);
    CHECK_FALSE(a.exists(1));
}

TEST_CASE("LayeredArray: layers are tracked per key independently", "[layered_array]")
{
    LayeredArray<int> a(3);
    a.set(0, 10, 0);
    a.set(1, 20, 0);
    a.set(2, 30, 0);

    a.request_layer(1);
    a.set(1, 21, 1);

    CHECK(a.current_layer(0) == 0);
    CHECK(a.current_layer(1) == 1);
    CHECK(a.current_layer(2) == 0);

    CHECK(a.get(0) == 10);
    CHECK(a.get(1) == 21);
    CHECK(a.get(2) == 30);

    // Backtracking key 1 leaves the others untouched and restores the old value.
    a.restore_layer(1);
    CHECK(a.get(1) == 20);
    CHECK(a.get(0) == 10);
    CHECK(a.get(2) == 30);
}

TEST_CASE("LayeredArray: a pushed layer does not disturb the one below", "[layered_array]")
{
    LayeredArray<int> a(1);
    a.set(0, 7, 0);
    a.request_layer(0);
    // The pushed layer's content is unspecified until written -- callers always set()
    // first -- but writing it must leave layer 0 intact.
    a.set(0, 9, 1);
    CHECK(a.get(0) == 9);
    a.restore_layer(0);
    CHECK(a.get(0) == 7);
}

TEST_CASE("LayeredArray: popping layer 0 returns a key to unwritten", "[layered_array]")
{
    LayeredArray<int> a(1);
    a.set(0, 5, 0);
    REQUIRE(a.exists(0));
    a.restore_layer(0);
    CHECK_FALSE(a.exists(0));
    CHECK_THROWS_AS(a.get(0), std::out_of_range);
    CHECK_THROWS_AS(a.restore_layer(0), std::out_of_range);
}

// The bug fixed on the AA-Pgen branch (Enum_fast_memory_map::set_value writing past the
// end of the buffer) reproduced against the new container: a map driven purely by set()
// walking 0,1,2,... must grow rather than corrupt the heap.
TEST_CASE("LayeredArray: set() beyond the allocated layers grows the storage", "[layered_array][growth]")
{
    LayeredArray<int> a(3, 1);   // one layer allocated up front
    REQUIRE(a.layer_capacity() == 1);

    for (std::size_t layer = 0; layer != 64; ++layer) {
        a.set(0, static_cast<int>(100 + layer), layer);
        a.set(1, static_cast<int>(200 + layer), layer);
        a.set(2, static_cast<int>(300 + layer), layer);
    }
    CHECK(a.layer_capacity() >= 64);

    // Every layer survived the reallocations.
    for (std::size_t layer = 0; layer != 64; ++layer) {
        CHECK(a.get(0, layer) == static_cast<int>(100 + layer));
        CHECK(a.get(1, layer) == static_cast<int>(200 + layer));
        CHECK(a.get(2, layer) == static_cast<int>(300 + layer));
    }
}

TEST_CASE("LayeredArray: request_layer beyond the allocated layers grows the storage", "[layered_array][growth]")
{
    LayeredArray<int> a(2, 1);
    a.set(0, 1, 0);
    for (int i = 0; i != 40; ++i) {
        a.request_layer(0);
    }
    CHECK(a.current_layer(0) == 40);
    a.set(0, 42, 40);
    CHECK(a.get(0) == 42);
    CHECK(a.get(0, 0) == 1);       // layer 0 survived every reallocation
    CHECK_FALSE(a.exists(1));      // untouched key unaffected by the reallocations
}

TEST_CASE("LayeredArray: layers must be filled bottom-up", "[layered_array]")
{
    LayeredArray<int> a(2);
    // Key 0 is unwritten (-1), so only layer 0 is legal.
    CHECK_THROWS_AS(a.set(0, 1, 1), std::out_of_range);
    a.set(0, 1, 0);
    CHECK_NOTHROW(a.set(0, 2, 1));
    CHECK_THROWS_AS(a.set(0, 3, 3), std::out_of_range);
}

TEST_CASE("LayeredArray: get(key, layer) does not rewind", "[layered_array]")
{
    LayeredArray<int> a(1);
    a.set(0, 1, 0);
    a.set(0, 2, 1);
    REQUIRE(a.current_layer(0) == 1);

    // Enum_fast_memory_map::at(key, layer) silently moved the current layer here.
    CHECK(a.get(0, 0) == 1);
    CHECK(a.current_layer(0) == 1);

    // The rewind is available, but only when asked for explicitly.
    a.set_current_layer(0, 0);
    CHECK(a.current_layer(0) == 0);
    CHECK(a.get(0) == 1);
    CHECK_THROWS_AS(a.set_current_layer(0, 5), std::out_of_range);
}

TEST_CASE("LayeredArray: bounds are checked on every key", "[layered_array]")
{
    LayeredArray<int> a(2);
    CHECK_THROWS_AS(a.set(2, 0, 0), std::out_of_range);
    CHECK_THROWS_AS(a.exists(2), std::out_of_range);
    CHECK_THROWS_AS(a.current_layer(9), std::out_of_range);
}

TEST_CASE("LayeredArray: row and layer views", "[layered_array]")
{
    LayeredArray<int> a(3, 2);
    a.init_first_layer(4);
    CHECK(a.exists(0));
    CHECK(a.exists(2));

    const std::span<const int> row = a.layer(0);
    REQUIRE(row.size() == 3);
    CHECK(row[0] == 4);
    CHECK(row[2] == 4);

    const std::span<const int> layers = a.current_layers();
    REQUIRE(layers.size() == 3);
    CHECK(layers[0] == 0);

    a.request_layer(1);
    CHECK(a.current_layers()[1] == 1);
    CHECK(a.current_layers()[0] == 0);
}

TEST_CASE("LayeredArray: copying is value semantics", "[layered_array]")
{
    // Enum_fast_memory_map owned raw new[] with no copy constructor: copying double-freed.
    LayeredArray<int> a(2);
    a.set(0, 1, 0);
    LayeredArray<int> b(a);
    b.set(0, 2, 0);
    CHECK(a.get(0) == 1);
    CHECK(b.get(0) == 2);
}

// ---------------------------------------------------------------------------
// Benchmarks against the container being replaced (plan decision D4, step 1).
//
// Two regimes, because they answer different questions:
//
//  - steady state: push/write/read/pop, balanced, so layer capacity stabilises and
//    no growth occurs. This is what the scenario traversal actually does, and it is
//    the number that matters for the B8 migration risk.
//  - growth: layers only ever pushed. Enum_fast_memory_map reallocates and copies the
//    whole buffer once per layer (quadratic); LayeredArray doubles (amortised linear).
// ---------------------------------------------------------------------------
TEST_CASE("LayeredArray vs Enum_fast_memory_map: steady state", "[!benchmark][layered_array]")
{
    constexpr std::size_t kKeys = 6;
    const std::size_t order[] = {0, 3, 1, 5, 2, 4};

    BENCHMARK("Enum_fast_memory_map: balanced push/write/read/pop")
    {
        Enum_fast_memory_map<int, Seq_Offset> m(kKeys);
        Seq_Offset sink = 0;
        for (std::size_t rep = 0; rep != 1000; ++rep) {
            for (std::size_t k : order) {
                m.request_memory_layer(static_cast<int>(k));
                const int l = m.get_current_memory_layer(static_cast<int>(k));
                m.set_value(static_cast<int>(k), static_cast<Seq_Offset>(k), l);
                sink += m.at(static_cast<int>(k));
                m.at(static_cast<int>(k), l - 1 >= 0 ? l - 1 : 0);   // pop == rewind
            }
        }
        return sink;
    };

    BENCHMARK("LayeredArray: balanced push/write/read/pop")
    {
        LayeredArray<Seq_Offset> a(kKeys, 8);
        Seq_Offset sink = 0;
        for (std::size_t rep = 0; rep != 1000; ++rep) {
            for (std::size_t k : order) {
                a.request_layer(k);
                a.set(k, static_cast<Seq_Offset>(k), static_cast<std::size_t>(a.current_layer(k)));
                sink += a.get(k);
                a.restore_layer(k);
            }
        }
        return sink;
    };
}

TEST_CASE("LayeredArray vs Enum_fast_memory_map: unbounded layer growth", "[!benchmark][layered_array]")
{
    constexpr std::size_t kKeys = 6;

    BENCHMARK("Enum_fast_memory_map: 2000 layers, realloc per layer")
    {
        Enum_fast_memory_map<int, Seq_Offset> m(kKeys);
        for (std::size_t d = 0; d != 2000; ++d) {
            for (std::size_t k = 0; k != kKeys; ++k) {
                m.request_memory_layer(static_cast<int>(k));
            }
        }
        return m.get_current_memory_layer(0);
    };

    BENCHMARK("LayeredArray: 2000 layers, doubling")
    {
        LayeredArray<Seq_Offset> a(kKeys, 1);
        for (std::size_t d = 0; d != 2000; ++d) {
            for (std::size_t k = 0; k != kKeys; ++k) {
                a.request_layer(k);
            }
        }
        return a.current_layer(0);
    };
}
