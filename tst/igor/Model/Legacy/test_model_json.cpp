/**
 * @file test_model_json.cpp
 * @brief Tests for the JSON serialization of a parsed model (ModelJson.h).
 *
 * Step 1 of the JSON format migration: Core still parses the text file, and this checks what
 * comes out the other side. The reading direction, through a factory and per-class
 * constructors taking a json node, is step 2; the round-trip test lands with it.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Insertion.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <memory>
#include <string>

namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy {}
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

static const std::string TEST_DATA_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/format_v2/";
static const std::string MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/models/";

namespace {

/// Every event node carries these, whatever its type.
void require_event_header(const nlohmann::json &node)
{
    for (const char *key : { "type", "gene_class", "seq_type", "side", "priority", "nickname",
                             "realizations", "parents" })
        REQUIRE(node.contains(key));
}

const nlohmann::json &event_by_nickname(const nlohmann::json &doc, const std::string &nickname)
{
    for (const auto &node : doc.at("events"))
        if (node.at("nickname").get<std::string>() == nickname)
            return node;
    throw std::runtime_error("no event with nickname \"" + nickname + "\" in the document");
}

/// Local dispatch, replaced by the event factory in the next step. Kept here so the
/// constructors can be tested before the factory exists.
std::shared_ptr<Rec_Event> event_from_json(const nlohmann::json &node)
{
    const auto type = node.at("type").get<std::string>();
    if (type == "GeneChoice")
        return std::make_shared<Gene_choice>(node);
    if (type == "Deletion")
        return std::make_shared<Deletion>(node);
    if (type == "Insertion")
        return std::make_shared<Insertion>(node);
    if (type == "DinucMarkov")
        return std::make_shared<Dinucl_markov>(node);
    throw std::runtime_error("unknown event type \"" + type + "\"");
}

/// An event does not know its parents: that edge list is model-level.
nlohmann::json without_parents(nlohmann::json node)
{
    node.erase("parents");
    return node;
}

}  // namespace

TEST_CASE("model json: a v2 model serializes with its segment order", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));

    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    SECTION("the document carries its schema version") {
        REQUIRE(doc.at("schema_version").get<int>() == igor::model::legacy::kModelJsonSchemaVersion);
    }

    SECTION("the segment order survives as an array, in file order") {
        REQUIRE(doc.at("seq_type_order")
                == nlohmann::json::array({ "V_gene_seq", "VD_ins_seq", "D_gene_seq",
                                           "DJ_ins_seq", "J_gene_seq" }));
    }

    SECTION("every event of the file is present, with its header") {
        REQUIRE(doc.at("events").size() == 11);
        for (const auto &node : doc.at("events"))
            require_event_header(node);
    }

    SECTION("a gene choice keeps its class, its seq_type and its sequence") {
        const auto &d = event_by_nickname(doc, "d_choice");
        REQUIRE(d.at("type").get<std::string>() == "GeneChoice");
        REQUIRE(d.at("gene_class").get<std::string>() == "D_gene");
        REQUIRE(d.at("seq_type").get<std::string>() == "D_gene_seq");
        REQUIRE(d.at("realizations").size() == 1);
        REQUIRE(d.at("realizations").at(0).at("value_str").get<std::string>() == "GGGG");
        REQUIRE_FALSE(d.at("realizations").at(0).contains("value_int"));
    }

    SECTION("a deletion keeps its side and its integer realizations") {
        const auto &del = event_by_nickname(doc, "v_3_del");
        REQUIRE(del.at("type").get<std::string>() == "Deletion");
        REQUIRE(del.at("seq_type").get<std::string>() == "V_gene_seq");
        REQUIRE(del.at("side").get<std::string>() == "Three_prime");
        REQUIRE(del.at("realizations").at(0).at("value_int").get<int>() == 0);
    }

    SECTION("the error rate is serialized with its value") {
        REQUIRE(doc.at("error_rate").at("type").get<std::string>() == "SingleErrorRate");
        REQUIRE(doc.at("error_rate").at("rate").get<double>() == 0.001);
    }

    SECTION("serialization is reproducible") {
        REQUIRE(igor::model::legacy::model_parms_to_json(parms) == doc);
    }
}

TEST_CASE("model json: a legacy file serializes too, with the inferred order", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    // A legacy file carries no order, so read_model_parms infers the standard VDJ one.
    REQUIRE(doc.at("seq_type_order")
            == nlohmann::json::array({ "V_gene_seq", "VD_ins_seq", "D_gene_seq",
                                       "DJ_ins_seq", "J_gene_seq" }));
}

TEST_CASE("model json: realizations come out sorted by index", "[Core][json][integration]")
{
    Model_Parms parms;
    const std::string shipped = MODELS_DIR + "human/tcr_beta/models/model_parms.txt";
    if (!std::ifstream(shipped).good())
        SKIP("models submodule not checked out: " + shipped);
    REQUIRE_NOTHROW(parms.read_model_parms(shipped));

    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);
    const auto &v_choice = event_by_nickname(doc, "v_choice");

    REQUIRE(v_choice.at("realizations").size() == 89);
    int previous = -1;
    for (const auto &r : v_choice.at("realizations")) {
        const int index = r.at("index").get<int>();
        REQUIRE(index > previous);
        previous = index;
    }

    // Parents are carried per event, so d_gene names both of its conditioning events.
    const auto &d_gene = event_by_nickname(doc, "d_gene");
    REQUIRE(d_gene.at("parents").size() == 2);
}

TEST_CASE("model json: json to object to json is the identity", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));
    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    for (const auto &node : doc.at("events")) {
        const auto rebuilt = event_from_json(node);
        REQUIRE(rebuilt->to_json() == without_parents(node));
    }
}

TEST_CASE("model json: round trip on a shipped model", "[Core][json][integration]")
{
    const std::string shipped = MODELS_DIR + "human/tcr_beta/models/model_parms.txt";
    if (!std::ifstream(shipped).good())
        SKIP("models submodule not checked out: " + shipped);

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(shipped));
    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    REQUIRE(doc.at("events").size() == 11);
    for (const auto &node : doc.at("events")) {
        const auto rebuilt = event_from_json(node);
        REQUIRE(rebuilt->to_json() == without_parents(node));
        // The generated name is what Model_marginals keys its index map by, so it has to
        // survive the round trip as well, not only the fields it is built from.
        REQUIRE(rebuilt->get_name() == parms.get_event_pointer(node.at("nickname")
                                                                       .get<std::string>(), true)
                                               ->get_name());
    }
}

TEST_CASE("model json: a rebuilt DinucMarkov gets the same name as the text reader's",
          "[Core][json]")
{
    // A DinucMarkov's generated name carries Undefined_side whatever its side, by decision:
    // the seq_type identifies the junction and the side is a direction only. Model_Parms
    // already keys it that way in get_events_map(), and write2txt_legacy() writes it that way
    // in the event line. This is what makes the name independent of the order of the setters,
    // and the name is a key — of Model_marginals::get_index_map(), of Model_Parms::edges, and
    // of the scenario and generation output columns.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));
    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    const nlohmann::json &node = event_by_nickname(doc, "vd_dinucl");
    REQUIRE(node.at("side").get<std::string>() == "Three_prime");

    const auto rebuilt = event_from_json(node);
    REQUIRE(rebuilt->get_side() == Three_prime);
    REQUIRE(rebuilt->get_name().find("Undefined_side") != std::string::npos);
    REQUIRE(rebuilt->get_name()
            == parms.get_event_pointer("vd_dinucl", true)->get_name());

    SECTION("the three name accessors agree on the side token") {
        // Only on the side. They disagree on the class token for a DinucMarkov, and that is a
        // separate matter: Dinucl_markov::update_event_name() derives it from ins_seq_type
        // ("VD_genes") while the two accessors read event_class, which the reader leaves
        // Undefined_gene. Worth fixing, not here.
        for (const auto &name : { rebuilt->get_name(), rebuilt->get_legacy_name(),
                                  rebuilt->get_v2_name() }) {
            INFO("name: " << name);
            REQUIRE(name.find("Undefined_side") != std::string::npos);
            REQUIRE(name.find("Three_prime") == std::string::npos);
        }
    }

    SECTION("the name no longer depends on the order of the setters") {
        const auto name = rebuilt->get_name();

        rebuilt->set_event_side(Five_prime);        // refreshes the name, as any setter does
        REQUIRE(rebuilt->get_side() == Five_prime);
        REQUIRE(rebuilt->get_name() == name);

        rebuilt->update_event_name();               // and calling it again changes nothing
        REQUIRE(rebuilt->get_name() == name);
    }

    SECTION("an event whose side IS its identity still carries it") {
        // The rule is specific to DinucMarkov. A deletion is identified by its side, so its
        // name must keep it, and set_event_side() must move it.
        const auto deletion = event_from_json(event_by_nickname(doc, "v_3_del"));
        REQUIRE(deletion->get_name().find("Three_prime") != std::string::npos);

        deletion->set_event_side(Five_prime);
        REQUIRE(deletion->get_name().find("Five_prime") != std::string::npos);
        REQUIRE(deletion->get_name().find("Three_prime") == std::string::npos);
    }
}

TEST_CASE("model json: a malformed event node is rejected", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));
    const nlohmann::json doc = igor::model::legacy::model_parms_to_json(parms);

    SECTION("an unknown key, a typo for instance") {
        nlohmann::json node = event_by_nickname(doc, "v_3_del");
        node["nickmame"] = "v_3_del";
        REQUIRE_THROWS_AS(event_from_json(node), std::runtime_error);
    }

    SECTION("a missing key") {
        nlohmann::json node = event_by_nickname(doc, "v_3_del");
        node.erase("priority");
        REQUIRE_THROWS_AS(event_from_json(node), std::runtime_error);
    }

    SECTION("a type that does not match the constructor") {
        nlohmann::json node = event_by_nickname(doc, "v_3_del");
        node["type"] = "Insertion";
        REQUIRE_THROWS_AS(std::make_shared<Deletion>(node), std::runtime_error);
    }

    SECTION("a gap in the realization indices") {
        nlohmann::json node = event_by_nickname(doc, "v_3_del");
        node["realizations"].at(0)["index"] = 7;
        REQUIRE_THROWS_AS(event_from_json(node), std::runtime_error);
    }

    SECTION("a gene class the seq_type contradicts") {
        nlohmann::json node = event_by_nickname(doc, "v_3_del");
        node["gene_class"] = "J_gene";
        REQUIRE_THROWS_AS(event_from_json(node), std::runtime_error);
    }
}
