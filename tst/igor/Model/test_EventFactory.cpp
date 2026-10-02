/*
 * test_EventFactory.cpp
 *
 * Tests for the event factory.
 *
 * The factory maps the type name a model document carries to a creator that builds the
 * concrete event FROM that document node. There is no blank object and no second
 * configuration phase: an event is complete as soon as it exists, so it cannot answer a
 * capability query with a value nobody chose. What each node must contain is the concrete
 * constructor's business, and the tests for that live in Core, in test_model_json.cpp.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Insertion.h>
#include <igor/Model/EventFactory.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <string>

namespace {

/// A minimal but complete node, of the shape Rec_Event::to_json() emits.
nlohmann::json deletion_node()
{
    return { { "type", "Deletion" },
             { "gene_class", "V_gene" },
             { "seq_type", "V_gene_seq" },
             { "side", "Three_prime" },
             { "priority", 5 },
             { "nickname", "v_3_del" },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "name", "0" }, { "value_int", 0 } },
                                       { { "index", 1 }, { "name", "1" }, { "value_int", 1 } } }) } };
}

nlohmann::json gene_choice_node()
{
    return { { "type", "GeneChoice" },
             { "gene_class", "D_gene" },
             { "seq_type", "D_gene_seq" },
             { "side", "Undefined_side" },
             { "priority", 6 },
             { "nickname", "d_gene" },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "name", "D1" }, { "value_str", "GGGG" } } }) } };
}

nlohmann::json insertion_node()
{
    return { { "type", "Insertion" },
             { "gene_class", "Undefined_gene" },
             { "seq_type", "VD_ins_seq" },
             { "side", "Undefined_side" },
             { "priority", 4 },
             { "nickname", "vd_ins" },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "name", "0" }, { "value_int", 0 } } }) } };
}

nlohmann::json dinucl_node()
{
    return { { "type", "DinucMarkov" },
             { "gene_class", "Undefined_gene" },
             { "seq_type", "VD_ins_seq" },
             { "side", "Three_prime" },
             { "priority", 3 },
             { "nickname", "vd_dinucl" },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "name", "A" }, { "value_str", "A" } },
                                       { { "index", 1 }, { "name", "C" }, { "value_str", "C" } },
                                       { { "index", 2 }, { "name", "G" }, { "value_str", "G" } },
                                       { { "index", 3 }, { "name", "T" }, { "value_str", "T" } } }) } };
}

}  // namespace

TEST_CASE("EventFactory: the four concrete types are registered", "[factory]")
{
    using namespace igor::model::event_factory;

    SECTION("by the name the document carries") {
        REQUIRE(is_registered("GeneChoice"));
        REQUIRE(is_registered("Deletion"));
        REQUIRE(is_registered("Insertion"));
        // The file says DinucMarkov, which is not how the enumerator is spelled.
        REQUIRE(is_registered("DinucMarkov"));
    }

    SECTION("and nothing else") {
        REQUIRE_FALSE(is_registered("Dinuclmarkov"));
        REQUIRE_FALSE(is_registered("Undefined"));
        REQUIRE_FALSE(is_registered(""));
        REQUIRE(registered_type_names().size() == 4);
    }
}

TEST_CASE("EventFactory: create builds the concrete class the node names", "[factory]")
{
    using namespace igor::model::event_factory;

    SECTION("Deletion") {
        const auto event = create(deletion_node());
        REQUIRE(std::dynamic_pointer_cast<Deletion>(event));
        REQUIRE(event->get_type() == Deletion_t);
    }

    SECTION("GeneChoice") {
        const auto event = create(gene_choice_node());
        REQUIRE(std::dynamic_pointer_cast<Gene_choice>(event));
        REQUIRE(event->get_class() == D_gene);
    }

    SECTION("Insertion") {
        const auto event = create(insertion_node());
        REQUIRE(std::dynamic_pointer_cast<Insertion>(event));
    }

    SECTION("DinucMarkov") {
        const auto event = create(dinucl_node());
        REQUIRE(std::dynamic_pointer_cast<Dinucl_markov>(event));
        // 4 states, so size() counts the 16 transitions while the tensor is 4 by 4.
        REQUIRE(event->size() == 16);
        REQUIRE(event->inherent_shape() == std::vector<std::size_t>{ 4, 4 });
    }
}

TEST_CASE("EventFactory: the event it returns is complete", "[factory]")
{
    using namespace igor::model::event_factory;

    const auto event = create(deletion_node());

    REQUIRE(event->get_seq_type() == "V_gene_seq");
    REQUIRE(event->get_side() == Three_prime);
    REQUIRE(event->get_priority() == 5);
    REQUIRE(event->get_nickname() == "v_3_del");
    REQUIRE(event->size() == 2);
    // Nothing had to be set afterwards, so the generated name is already the final one.
    REQUIRE(event->get_name() == "Deletion_V_gene_Three_prime_prio5_size2");
}

TEST_CASE("EventFactory: a node it cannot dispatch is rejected", "[factory]")
{
    using namespace igor::model::event_factory;

    SECTION("no type at all") {
        nlohmann::json node = deletion_node();
        node.erase("type");
        REQUIRE_THROWS_AS(create(node), std::runtime_error);
    }

    SECTION("a type that is not a string") {
        nlohmann::json node = deletion_node();
        node["type"] = 42;
        REQUIRE_THROWS_AS(create(node), std::runtime_error);
    }

    SECTION("an unknown type") {
        nlohmann::json node = deletion_node();
        node["type"] = "Hypermutation";
        REQUIRE_THROWS_AS(create(node), std::runtime_error);
    }

    SECTION("a node the concrete constructor refuses") {
        nlohmann::json node = deletion_node();
        node.erase("priority");
        REQUIRE_THROWS_AS(create(node), std::runtime_error);
    }
}

TEST_CASE("EventFactory: registering a null creator is refused", "[factory]")
{
    using namespace igor::model::event_factory;

    REQUIRE_THROWS_AS(register_creator("Whatever", nullptr), std::invalid_argument);
    REQUIRE_THROWS_AS(register_creator("", [](const nlohmann::json &) { return nullptr; }),
                      std::invalid_argument);
    REQUIRE_FALSE(is_registered("Whatever"));
}
