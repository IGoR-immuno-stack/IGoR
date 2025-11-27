#include <catch2/catch_test_macros.hpp>
#include <igor/Core/SequenceTypes.h>

TEST_CASE("SequenceTypeRegistry Basic Operations", "[SequenceTypes]")
{
    SequenceTypeRegistry &registry = get_sequence_type_registry();

    SECTION("Predefined types exist")
    {
        REQUIRE(registry.get_type_id("V_gene_seq") == SequenceTypeRegistry::V_GENE_SEQ);
        REQUIRE(registry.get_type_id("VD_ins_seq") == SequenceTypeRegistry::VD_INS_SEQ);
        REQUIRE(registry.get_type_id("D_gene_seq") == SequenceTypeRegistry::D_GENE_SEQ);
        REQUIRE(registry.get_type_id("DJ_ins_seq") == SequenceTypeRegistry::DJ_INS_SEQ);
        REQUIRE(registry.get_type_id("J_gene_seq") == SequenceTypeRegistry::J_GENE_SEQ);
        REQUIRE(registry.get_type_id("VJ_ins_seq") == SequenceTypeRegistry::VJ_INS_SEQ);
    }

    SECTION("Register new type")
    {
        auto id = registry.register_type("New_Type");
        REQUIRE(id >= 6);
        REQUIRE(registry.get_type_id("New_Type") == id);
        REQUIRE(registry.get_type_name(id) == "New_Type");
    }

    SECTION("Register existing type returns same ID")
    {
        auto id1 = registry.register_type("Existing_Type");
        auto id2 = registry.register_type("Existing_Type");
        REQUIRE(id1 == id2);
    }

    SECTION("Tandem D helpers")
    {
        auto d1_id = registry.register_d_gene(1);
        REQUIRE(registry.get_type_name(d1_id) == "D1_gene_seq");

        auto d1d2_id = registry.register_d_insertion(1, 2);
        REQUIRE(registry.get_type_name(d1d2_id) == "D1D2_ins_seq");
    }

    SECTION("Error handling")
    {
        REQUIRE_THROWS_AS(registry.get_type_id("NonExistent"), std::runtime_error);
        REQUIRE_THROWS_AS(registry.get_type_name(9999), std::runtime_error);
    }
}
