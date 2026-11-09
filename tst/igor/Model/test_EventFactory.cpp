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
#include <catch2/matchers/catch_matchers_string.hpp>

#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Insertion.h>
#include <igor/Model/EventFactory.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <string>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

namespace {

/// A minimal but complete node, of the shape Rec_Event::to_json() emits. The nickname is not in
/// it: it is the node's key in a document, and create()'s first argument.
nlohmann::json deletion_node()
{
    return { { "type", "Deletion" },
             { "priority", 5 },
             { "targets", nlohmann::json::array({ { { "segment", "V_gene_seq" }, { "side", "3'" } } }) },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "label", 0 } }, { { "index", 1 }, { "label", 1 } } }) } };
}

/// A gene choice names its library; its templates are the document's library's.
nlohmann::json gene_choice_node()
{
    return { { "type", "GeneChoice" },
             { "library", "D_gene" },
             { "priority", 6 },
             { "targets", nlohmann::json::array({ { { "segment", "D_gene_seq" } } }) },
             { "realizations", nlohmann::json::array({ { { "index", 0 }, { "label", "D1" } } }) } };
}

/// The document's "library" section the gene choice draws from.
const nlohmann::json kLibrary = { { "D_gene", { { "D1", { { "sequence", "GGGG" } } } } } };

nlohmann::json insertion_node()
{
    return { { "type", "Insertion" },
             { "priority", 4 },
             { "targets", nlohmann::json::array({ { { "segment", "VD_ins_seq" } } }) },
             { "realizations", nlohmann::json::array({ { { "index", 0 }, { "label", 0 } } }) } };
}

nlohmann::json dinucl_node()
{
    return { { "type", "DinucMarkov" },
             { "priority", 3 },
             { "targets", nlohmann::json::array({ { { "segment", "VD_ins_seq" }, { "side", "3'" } } }) },
             { "realizations",
               nlohmann::json::array({ { { "index", 0 }, { "label", "A" } },
                                       { { "index", 1 }, { "label", "C" } },
                                       { { "index", 2 }, { "label", "G" } },
                                       { { "index", 3 }, { "label", "T" } } }) } };
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
        const auto event = create("v_3_del", deletion_node(), kLibrary);
        REQUIRE(std::dynamic_pointer_cast<Deletion>(event));
        REQUIRE(event->get_type() == Deletion_t);
    }

    SECTION("GeneChoice, its template from the library") {
        const auto event = create("d_gene", gene_choice_node(), kLibrary);
        REQUIRE(std::dynamic_pointer_cast<Gene_choice>(event));
        REQUIRE(event->get_class() == D_gene);
        REQUIRE(event->realization_at(0).value_str == "GGGG");
    }

    SECTION("Insertion") {
        const auto event = create("vd_ins", insertion_node(), kLibrary);
        REQUIRE(std::dynamic_pointer_cast<Insertion>(event));
    }

    SECTION("DinucMarkov") {
        const auto event = create("vd_dinucl", dinucl_node(), kLibrary);
        REQUIRE(std::dynamic_pointer_cast<Dinucl_markov>(event));
        // 4 states, so size() counts the 16 transitions while the tensor is 4 by 4.
        REQUIRE(event->size() == 16);
        REQUIRE(event->inherent_shape() == std::vector<std::size_t>{ 4, 4 });
    }
}

TEST_CASE("EventFactory: the event it returns is complete", "[factory]")
{
    using namespace igor::model::event_factory;

    const auto event = create("v_3_del", deletion_node(), kLibrary);

    REQUIRE(event->get_seq_type() == "V_gene_seq");
    REQUIRE(event->get_side() == Three_prime);
    REQUIRE(event->get_priority() == 5);
    REQUIRE(event->get_nickname() == "v_3_del");
    REQUIRE(event->size() == 2);
    // Its name is its nickname (plan T2); its v1 name is complete with nothing set afterwards.
    REQUIRE(event->get_name() == "v_3_del");
    REQUIRE(event->get_v1_legacy_name() == "Deletion_V_gene_Three_prime_prio5_size2");
}

TEST_CASE("EventFactory: a node it cannot dispatch is rejected", "[factory]")
{
    using namespace igor::model::event_factory;

    SECTION("no type at all") {
        nlohmann::json node = deletion_node();
        node.erase("type");
        REQUIRE_THROWS_AS(create("v_3_del", node, kLibrary), std::runtime_error);
    }

    SECTION("a type that is not a string") {
        nlohmann::json node = deletion_node();
        node["type"] = 42;
        REQUIRE_THROWS_AS(create("v_3_del", node, kLibrary), std::runtime_error);
    }

    SECTION("an unknown type") {
        nlohmann::json node = deletion_node();
        node["type"] = "Hypermutation";
        REQUIRE_THROWS_AS(create("v_3_del", node, kLibrary), std::runtime_error);
    }

    SECTION("a node the concrete constructor refuses") {
        nlohmann::json node = deletion_node();
        node.erase("priority");
        REQUIRE_THROWS_AS(create("v_3_del", node, kLibrary), std::runtime_error);
    }

    SECTION("no nickname") {
        REQUIRE_THROWS_AS(create("", deletion_node(), kLibrary), std::runtime_error);
    }

    SECTION("a gene choice whose library is missing, or lacks its template") {
        REQUIRE_THROWS_AS(create("d_gene", gene_choice_node(), nlohmann::json()), std::runtime_error);
        nlohmann::json library = kLibrary;
        library["D_gene"].erase("D1");
        REQUIRE_THROWS_AS(create("d_gene", gene_choice_node(), library), std::runtime_error);
    }

    SECTION("a library of another shape, checked against its schema") {
        const nlohmann::json library = { { "D_gene", { { "D1", { { "seq", "GGGG" } } } } } };
        REQUIRE_THROWS_WITH(create("d_gene", gene_choice_node(), library),
                            Catch::Matchers::ContainsSubstring("/D_gene/D1: unknown key \"seq\""));
    }

    SECTION("a nickname inside the node, where the document's key is the nickname") {
        nlohmann::json node = deletion_node();
        node["nickname"] = "v_3_del";
        REQUIRE_THROWS_AS(create("v_3_del", node, kLibrary), std::runtime_error);
    }
}

TEST_CASE("EventFactory: registering a null creator is refused", "[factory]")
{
    using namespace igor::model::event_factory;

    REQUIRE_THROWS_AS(register_creator("Whatever", nullptr), std::invalid_argument);
    REQUIRE_THROWS_AS(register_creator("", [](const std::string &, const nlohmann::json &, const nlohmann::json &) {
                          return nullptr;
                      }),
                      std::invalid_argument);
    REQUIRE_FALSE(is_registered("Whatever"));
}
