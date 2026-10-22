/**
 * @file test_model_json.cpp
 * @brief Tests for the JSON model document (ModelJson.h): its layout, and the event nodes.
 *
 * What a document written from a v1 model looks like, section by section, against the layout of
 * SEGMENT_DECOMPOSITION_REVIEW.md §14; and that an event's node builds the event back. Reading
 * whole documents, and what a reader refuses in them, is test_model_formats.cpp.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/EventFactory.h>
#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <memory>
#include <set>
#include <string>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
using igor::model::event_factory::create;

static const std::string TEST_DATA_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/model_formats/";
static const std::string MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/models/";

namespace {

std::set<std::string> keys_of(const nlohmann::json &node)
{
    std::set<std::string> keys;
    for (const auto &entry : node.items())
        keys.insert(entry.key());
    return keys;
}

Model_Parms read_vdj()
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));
    return parms;
}

/// The event's own node, as Rec_Event::to_json() gives it: templates inline.
nlohmann::json node_of(const Model_Parms &parms, const std::string &nickname)
{
    return parms.get_event_pointer(nickname, true)->to_json();
}

}  // namespace

TEST_CASE("model json: a v1 model serializes into the planned layout", "[Core][json]")
{
    const Model_Parms parms = read_vdj();
    const nlohmann::json doc = model_parms_to_json(parms);

    SECTION("the document: its version and its sections") {
        REQUIRE(doc.at("schema_version").get<int>() == kModelJsonSchemaVersion);
        REQUIRE(kModelJsonSchemaVersion == 3);
        // A v1 model's events are not listed in nickname order, so the order goes with it as a hint.
        REQUIRE(keys_of(doc) == std::set<std::string>{ "schema_version", "structure", "library", "error_rate",
                                                       "engine_hints" });
        REQUIRE(keys_of(doc.at("structure")) == std::set<std::string>{ "layout", "events", "conditioning" });
    }

    SECTION("the layout is a Linear one, in file order") {
        REQUIRE(doc.at("structure").at("layout")
                == nlohmann::json{ { "type", "Linear" },
                                   { "segments", { "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq",
                                                   "J_gene_seq" } } });
    }

    SECTION("the events are keyed by nickname, and carry neither it nor their parents") {
        const nlohmann::json &events = doc.at("structure").at("events");
        REQUIRE(events.is_object());
        REQUIRE(events.size() == 11);
        for (const auto &[nickname, node] : events.items()) {
            INFO(nickname);
            const bool gene_choice = node.at("type") == "GeneChoice";
            std::set<std::string> expected{ "type", "priority", "targets", "realizations" };
            if (gene_choice)
                expected.insert("library");
            REQUIRE(keys_of(node) == expected);
        }
    }

    SECTION("a gene choice names its library; its template is in the library section") {
        const nlohmann::json &d = doc.at("structure").at("events").at("d_choice");
        REQUIRE(d.at("library") == "D_gene");
        REQUIRE(d.at("targets") == nlohmann::json::array({ { { "segment", "D_gene_seq" } } }));
        REQUIRE(d.at("realizations") == nlohmann::json::array({ { { "index", 0 }, { "label", "D_test" } } }));
        REQUIRE(doc.at("library").at("D_gene").at("D_test") == nlohmann::json{ { "sequence", "GGGG" } });
    }

    SECTION("a deletion targets an end, and its labels are its values") {
        const nlohmann::json &del = doc.at("structure").at("events").at("v_3_del");
        REQUIRE(del.at("targets")
                == nlohmann::json::array({ { { "segment", "V_gene_seq" }, { "side", "3'" } } }));
        REQUIRE(del.at("realizations") == nlohmann::json::array({ { { "index", 0 }, { "label", 0 } } }));
    }

    SECTION("a Markov chain keeps its direction, which the v1 line does not write") {
        const nlohmann::json &chain = doc.at("structure").at("events").at("vd_dinucl");
        REQUIRE(chain.at("targets")
                == nlohmann::json::array({ { { "segment", "VD_ins_seq" }, { "side", "3'" } } }));
        REQUIRE(chain.at("realizations").at(2) == nlohmann::json{ { "index", 2 }, { "label", "G" } });
    }

    SECTION("a model without edges has an empty conditioning") {
        REQUIRE(doc.at("structure").at("conditioning") == nlohmann::json::object());
    }

    SECTION("the error rate is serialized with its value") {
        REQUIRE(doc.at("error_rate") == nlohmann::json{ { "type", "SingleErrorRate" }, { "rate", 0.001 } });
    }

    SECTION("serialization is reproducible") {
        REQUIRE(model_parms_to_json(parms) == doc);
    }
}

TEST_CASE("model json: the conditioning is its own section, child to parents", "[Core][json][integration]")
{
    Model_Parms parms;
    const std::string shipped = MODELS_DIR + "human/tcr_beta/models/model_parms.txt";
    if (!std::ifstream(shipped).good())
        SKIP("models submodule not checked out: " + shipped);
    REQUIRE_NOTHROW(parms.read_model_parms(shipped));

    const nlohmann::json doc = model_parms_to_json(parms);
    const nlohmann::json &conditioning = doc.at("structure").at("conditioning");

    SECTION("every edge of the model, parents in nickname order") {
        REQUIRE(conditioning.at("d_gene") == nlohmann::json::array({ "j_choice", "v_choice" }));
        std::size_t edges = 0;
        for (const auto &event : parms.get_event_list()) {
            const auto parents = parms.get_parents(event);
            REQUIRE(conditioning.contains(event->get_nickname()) == not parents.empty());
            edges += parents.size();
        }
        std::size_t written = 0;
        for (const auto &[child, parents] : conditioning.items())
            written += parents.size();
        REQUIRE(written == edges);
    }

    SECTION("realizations come out sorted by index") {
        const nlohmann::json &v_choice = doc.at("structure").at("events").at("v_choice");
        REQUIRE(v_choice.at("realizations").size() == 89);
        int previous = -1;
        for (const auto &r : v_choice.at("realizations")) {
            const int index = r.at("index").get<int>();
            REQUIRE(index > previous);
            previous = index;
        }
        REQUIRE(doc.at("library").at("V_gene").size() == 89);
    }
}

TEST_CASE("model json: a legacy file serializes too, with the inferred order", "[Core][json]")
{
    const Model_Parms parms = read_vdj();
    const nlohmann::json doc = model_parms_to_json(parms);

    // A legacy file carries no order, so read_model_parms infers the standard VDJ one.
    REQUIRE(doc.at("structure").at("layout").at("segments")
            == nlohmann::json::array({ "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq" }));
}

TEST_CASE("model json: node to event to node is the identity", "[Core][json]")
{
    const Model_Parms parms = read_vdj();
    for (const auto &event : parms.get_event_list()) {
        INFO(event->get_nickname());
        const nlohmann::json node = event->to_json();
        const auto rebuilt = create(event->get_nickname(), node);
        REQUIRE(rebuilt->to_json() == node);
        REQUIRE(rebuilt->get_nickname() == event->get_nickname());
    }
}

TEST_CASE("model json: round trip on a shipped model", "[Core][json][integration]")
{
    const std::string shipped = MODELS_DIR + "human/tcr_beta/models/model_parms.txt";
    if (!std::ifstream(shipped).good())
        SKIP("models submodule not checked out: " + shipped);

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(shipped));

    REQUIRE(parms.get_event_list().size() == 11);
    for (const auto &event : parms.get_event_list()) {
        INFO(event->get_nickname());
        const nlohmann::json node = event->to_json();
        const auto rebuilt = create(event->get_nickname(), node);
        REQUIRE(rebuilt->to_json() == node);
        // The v1 name is what a v1 file's @Edges names the event by, so it has to survive the
        // round trip as well, not only the fields it is built from.
        REQUIRE(rebuilt->get_name() == event->get_name());
        REQUIRE(rebuilt->get_v1_legacy_name() == event->get_v1_legacy_name());
    }
}

TEST_CASE("model json: a rebuilt DinucMarkov gets the same v1 name as the text reader's",
          "[Core][json]")
{
    // A DinucMarkov's v1 name carries Undefined_side whatever its side, by decision: the
    // seq_type identifies the junction and the side is a direction only. Model_Parms keys it
    // that way in get_events_map(), and write2txt_legacy() writes it that way in the event line.
    // The v1 name is what a v1 file's @Edges section names the event by. Its name, the key of the
    // model's maps, is its nickname (plan T2).
    const Model_Parms parms = read_vdj();
    const nlohmann::json node = node_of(parms, "vd_dinucl");
    REQUIRE(node.at("targets").at(0).at("side") == "3'");

    const auto rebuilt = create("vd_dinucl", node);
    REQUIRE(rebuilt->get_side() == Three_prime);
    REQUIRE(rebuilt->get_name() == "vd_dinucl");
    REQUIRE(rebuilt->get_v1_legacy_name().find("Undefined_side") != std::string::npos);
    REQUIRE(rebuilt->get_v1_legacy_name()
            == parms.get_event_pointer("vd_dinucl", true)->get_v1_legacy_name());

    SECTION("the v1 name says Undefined_side, never the chain's direction") {
        REQUIRE(rebuilt->get_v1_legacy_name().find("Three_prime") == std::string::npos);
    }

    SECTION("the v1 name does not depend on the order of the setters") {
        const auto name = rebuilt->get_v1_legacy_name();

        rebuilt->set_event_side(Five_prime);
        REQUIRE(rebuilt->get_side() == Five_prime);
        REQUIRE(rebuilt->get_v1_legacy_name() == name);
    }

    SECTION("an event whose side IS its identity still carries it") {
        // The rule is specific to DinucMarkov. A deletion is identified by its side, so its
        // v1 name must keep it, and set_event_side() must move it.
        const auto deletion = create("v_3_del", node_of(parms, "v_3_del"));
        REQUIRE(deletion->get_v1_legacy_name().find("Three_prime") != std::string::npos);

        deletion->set_event_side(Five_prime);
        REQUIRE(deletion->get_v1_legacy_name().find("Five_prime") != std::string::npos);
        REQUIRE(deletion->get_v1_legacy_name().find("Three_prime") == std::string::npos);
    }
}

TEST_CASE("model json: a malformed event node is rejected", "[Core][json]")
{
    const Model_Parms parms = read_vdj();
    nlohmann::json node = node_of(parms, "v_3_del");

    SECTION("an unknown key, a typo for instance") {
        node["prioriti"] = 5;
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a missing key") {
        node.erase("priority");
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a type that does not match the constructor") {
        node["type"] = "Insertion";
        REQUIRE_THROWS_AS(std::make_shared<Deletion>("v_3_del", node), std::runtime_error);
    }

    SECTION("a gap in the realization indices") {
        node["realizations"].at(0)["index"] = 7;
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("two realizations under one label") {
        node["realizations"].push_back({ { "index", 1 }, { "label", 0 } });
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a label of the wrong kind") {
        node["realizations"].at(0)["label"] = "0";
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a realization in the v2-era spelling") {
        node["realizations"].at(0)["value_int"] = 0;
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("no target, two targets, or a target without a segment") {
        const nlohmann::json target = node.at("targets").at(0);
        node["targets"] = nlohmann::json::array();
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
        node["targets"] = nlohmann::json::array({ target, target });
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
        node["targets"] = nlohmann::json::array({ { { "side", "3'" } } });
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a side spelled otherwise than 5' or 3'") {
        node["targets"].at(0)["side"] = "Three_prime";
        REQUIRE_THROWS_AS(create("v_3_del", node), std::runtime_error);
    }

    SECTION("a library on an event that is not a gene choice") {
        for (const char *nickname : { "v_3_del", "vd_ins", "vd_dinucl" }) {
            nlohmann::json other = node_of(parms, nickname);
            other["library"] = "V_gene";
            REQUIRE_THROWS_AS(create(nickname, other), std::runtime_error);
        }
    }

    SECTION("a gene choice's realization without its template") {
        nlohmann::json gene = node_of(parms, "d_choice");
        gene["realizations"].at(0).erase("sequence");
        REQUIRE_THROWS_AS(create("d_choice", gene), std::runtime_error);
    }
}
