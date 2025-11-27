#include <catch2/catch_test_macros.hpp>
#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/SequenceTypes.h>
#include <memory>

TEST_CASE("DynamicSequenceMap Operations", "[DynamicSequenceMap]")
{
    // Register a dynamic type
    auto &registry = get_sequence_type_registry();
    auto dynamic_type_id = registry.register_type("Test_Dynamic_Type");

    // Instantiate map (should pick up the new size)
    DynamicSequenceMap<int> map;

    SECTION("Standard types access")
    {
        map.init_first_layer(-1);

        map.set_value(SequenceTypeRegistry::V_GENE_SEQ, 100, 0);
        REQUIRE(map.at(SequenceTypeRegistry::V_GENE_SEQ) == 100);

        map.set_value(SequenceTypeRegistry::J_GENE_SEQ, 200, 0);
        REQUIRE(map.at(SequenceTypeRegistry::J_GENE_SEQ) == 200);
    }

    SECTION("Dynamic type access")
    {
        map.init_first_layer(-1);

        // Use the dynamic ID as key
        map.set_value(dynamic_type_id, 300, 0);
        REQUIRE(map.at(dynamic_type_id) == 300);
    }

    SECTION("Memory layers")
    {
        map.init_first_layer(0);

        // Request new layer for dynamic type
        map.request_memory_layer(dynamic_type_id);
        map.set_value(dynamic_type_id, 400, 1);

        REQUIRE(map.at(dynamic_type_id, 1) == 400);
        REQUIRE(map.at(dynamic_type_id, 0) == 0); // Previous layer
    }

    SECTION("Bounds checking")
    {
        // Access out of bounds (beyond registered types)
        REQUIRE_THROWS_AS(map.at(registry.size() + 10), std::out_of_range);
    }
}
