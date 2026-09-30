/**
 * @file test_model_json.cpp
 * @brief Tests for the JSON serialization of a parsed model (ModelJson.h).
 *
 * Step 1 of the JSON format migration: Core still parses the text file, and this checks what
 * comes out the other side. The reading direction, through a factory and per-class
 * constructors taking a json node, is step 2; the round-trip test lands with it.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/ModelJson.h>
#include <igor/Core/Model_Parms.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <string>

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

}  // namespace

TEST_CASE("model json: a v2 model serializes with its segment order", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));

    const nlohmann::json doc = igor::model_parms_to_json(parms);

    SECTION("the document carries its schema version") {
        REQUIRE(doc.at("schema_version").get<int>() == igor::kModelJsonSchemaVersion);
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
        REQUIRE(igor::model_parms_to_json(parms) == doc);
    }
}

TEST_CASE("model json: a legacy file serializes too, with the inferred order", "[Core][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    const nlohmann::json doc = igor::model_parms_to_json(parms);

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

    const nlohmann::json doc = igor::model_parms_to_json(parms);
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
