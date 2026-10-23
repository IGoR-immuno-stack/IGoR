/*
 * test_ModelSchema.cpp
 *
 * The JSON Schemas of the model document (ModelSchema.h): that they are valid schemas, published
 * where their "$id" says, embedded as written, complete for the registered event types; and that
 * everything IGoR writes follows them.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <igor/Model/EventFactory.h>
#include <igor/Model/ModelSchema.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <nlohmann/json-schema.hpp>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
namespace schema = igor::model::schema;
using Catch::Matchers::ContainsSubstring;

namespace {

const std::string kSourceDir = std::string(IGOR_SOURCE_DIR) + "/src/igor/Model/schemas/0.3/";

std::string text_of(const std::string &path)
{
    std::ifstream in(path);
    return { std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>() };
}

/// A model read from \a parms_file, written whole: structure, library, parameters, error rate.
nlohmann::json complete_document(const std::string &parms_file, const std::string &marginals_file = "")
{
    Model_Parms parms;
    parms.read_model_parms(parms_file);
    Model_marginals marginals(parms);
    if (marginals_file.empty())
        marginals.uniform_initialize(parms);
    else
        marginals.read_marginals(marginals_file, parms);
    nlohmann::json doc = model_parms_to_json(parms);
    doc["parameters"] = model_marginals_to_json(parms, marginals);
    return doc;
}

}  // namespace

TEST_CASE("Model schemas: each file is a draft-07 schema, published where its $id says", "[model_format][schema]")
{
    // The meta-schema's formats: the library checks some, not the URI ones, which the $id check
    // below covers for this purpose.
    const auto formats = [](const std::string &format, const std::string &value) {
        if (format == "uri" or format == "uri-reference" or format == "regex")
            return;
        nlohmann::json_schema::default_string_format_check(format, value);
    };
    nlohmann::json_schema::json_validator meta(nlohmann::json_schema::draft7_schema_builtin, nullptr, formats);
    const auto files = schema::schema_files();
    REQUIRE(files.size() == 12);
    for (const auto &file : files) {
        INFO(file);
        const nlohmann::json content = nlohmann::json::parse(schema::schema_text(file));
        REQUIRE_NOTHROW(meta.validate(content));
        CHECK(content.at("$schema") == "http://json-schema.org/draft-07/schema#");
        CHECK(content.at("$id") == std::string(schema::kBaseUri) + file);
        CHECK(content.contains("description"));
    }
    CHECK(schema::document_schema_uri() == std::string(schema::kBaseUri) + "model.schema.json");
}

TEST_CASE("Model schemas: the embedded text is the source files'", "[model_format][schema]")
{
    std::set<std::string> on_disk;
    for (const auto &entry : std::filesystem::directory_iterator(kSourceDir))
        on_disk.insert(entry.path().filename().string());
    const auto files = schema::schema_files();
    REQUIRE(std::set<std::string>(files.begin(), files.end()) == on_disk);
    for (const auto &file : files) {
        INFO(file);
        CHECK(std::string(schema::schema_text(file)) == text_of(kSourceDir + file));
    }
    CHECK_THROWS_AS(schema::schema_text("absent.schema.json"), std::out_of_range);
}

TEST_CASE("Model schemas: the event union lists exactly the registered event types", "[model_format][schema]")
{
    const nlohmann::json union_schema = nlohmann::json::parse(schema::schema_text("event.schema.json"));
    const auto listed = union_schema.at("properties").at("type").at("enum").get<std::vector<std::string>>();
    const auto registered = igor::model::event_factory::registered_type_names();
    REQUIRE(std::set<std::string>(listed.begin(), listed.end())
            == std::set<std::string>(registered.begin(), registered.end()));

    std::set<std::string> dispatched;
    for (const auto &branch : union_schema.at("allOf")) {
        const std::string type = branch.at("if").at("properties").at("type").at("const");
        CHECK(branch.at("then").at("$ref") == "event-" + type + ".schema.json");
        dispatched.insert(type);
        const nlohmann::json own = nlohmann::json::parse(schema::schema_text("event-" + type + ".schema.json"));
        CHECK(own.at("properties").at("type").at("const") == type);
    }
    CHECK(dispatched == std::set<std::string>(listed.begin(), listed.end()));
}

TEST_CASE("Model schemas: what IGoR writes follows them", "[model_format][schema]")
{
    const std::string data = std::string(IGOR_SOURCE_DIR);

    SECTION("complete documents written from v1 models, with their parameters") {
        for (const auto &[parms, marginals] : std::vector<std::pair<std::string, std::string>>{
                     { data + "/tst/test_data/model_formats/test_legacy_model_parms.txt", "" },
                     { data + "/tst/test_data/model_formats/test_legacy_vdj_model_parms.txt", "" },
                     { data + "/scripts/tests/data/input/TRB_model_parms.txt",
                       data + "/scripts/tests/data/input/TRB_uniform_model_marginals.txt" },
                     { data + "/tst/test_data/inference/fixed_VJ_TRB_model_parms.txt",
                       data + "/tst/test_data/inference/fixed_VJ_TRB_marginals.txt" } }) {
            INFO(parms);
            const nlohmann::json doc = complete_document(parms, marginals);
            CHECK(doc.at("$schema") == schema::document_schema_uri());
            REQUIRE_NOTHROW(schema::validate(doc, "model.schema.json", "the document", { "structure", "parameters" }));
        }
    }

    SECTION("a shipped model") {
        const std::string shipped = data + "/models/human/tcr_beta/models/";
        if (!std::filesystem::exists(shipped + "model_parms.txt"))
            SKIP("models submodule not checked out");
        const nlohmann::json doc = complete_document(shipped + "model_parms.txt", shipped + "model_marginals.txt");
        REQUIRE_NOTHROW(schema::validate(doc, "model.schema.json", "the document", { "structure", "parameters" }));
    }

    SECTION("the JSON fixtures, and a fragment") {
        for (const char *name : { "vj_model.json", "vdj_model.json", "tandem_d_model.json" }) {
            INFO(name);
            std::ifstream in(data + "/tst/test_data/model_formats/" + name);
            REQUIRE_NOTHROW(schema::validate(parse_model_document(in), "model.schema.json", "the document"));
        }
        Model_Parms parms;
        parms.read_model_parms(data + "/tst/test_data/model_formats/vj_model.json");
        Model_marginals marginals(parms);
        marginals.uniform_initialize(parms);
        const std::string fragment = (std::filesystem::temp_directory_path() / "igor_test_schema_fragment.json").string();
        marginals.write_marginals_json(fragment, parms);
        // Read whole and closed before the removal, which Windows refuses on an open file.
        std::istringstream in(text_of(fragment));
        std::filesystem::remove(fragment);
        const nlohmann::json doc = parse_model_document(in);
        CHECK(doc.at("$schema") == schema::document_schema_uri());
        REQUIRE_NOTHROW(schema::validate(doc, "model.schema.json", "the fragment", { "parameters" }));
    }
}

TEST_CASE("Model schemas: a violation is named by its JSON pointer, every one of them", "[model_format][schema]")
{
    nlohmann::json doc = complete_document(std::string(IGOR_SOURCE_DIR) + "/scripts/tests/data/input/TRB_model_parms.txt");

    SECTION("each at its place") {
        doc["structure"]["events"]["v_3_del"]["targets"][0]["side"] = "Three_prime";
        doc["structure"]["events"]["vd_ins"]["realizations"][0]["label"] = -1;
        doc["error_rate"]["rate"] = 2;
        doc["engine_hints"]["legacy_dfs"]["order"] = nlohmann::json::array();
        REQUIRE_THROWS_WITH(schema::validate(doc, "model.schema.json", "the document"),
                            ContainsSubstring("/structure/events/v_3_del/targets/0/side: instance not found in "
                                              "required enum (\"Three_prime\")")
                                    && ContainsSubstring("/structure/events/vd_ins/realizations/0/label: instance is "
                                                         "below minimum")
                                    && ContainsSubstring("/error_rate/rate: instance exceeds maximum")
                                    && ContainsSubstring("/engine_hints/legacy_dfs: unknown key \"order\"")
                                    // Once each: neither an allOf's summary nor its case's prefix.
                                    && !ContainsSubstring("at least one subschema")
                                    && !ContainsSubstring("[combination"));
    }
    SECTION("ten of them, and how many more") {
        for (auto &value : doc["parameters"]["nodes"]["v_choice"]["values"])
            value = -1;
        REQUIRE_THROWS_WITH(schema::validate(doc, "model.schema.json", "the document"),
                            ContainsSubstring("/parameters/nodes/v_choice/values/9:")
                                    && !ContainsSubstring("/parameters/nodes/v_choice/values/10:")
                                    && ContainsSubstring("more"));
    }
    SECTION("a section a reader needs") {
        doc.erase("parameters");
        REQUIRE_THROWS_WITH(schema::validate(doc, "model.schema.json", "the document", { "parameters" }),
                            ContainsSubstring("/: required property 'parameters' not found"));
    }
    SECTION("an event node on its own, against its type's schema") {
        nlohmann::json node = doc["structure"]["events"]["v_3_del"];
        node.erase("priority");
        REQUIRE_THROWS_WITH(schema::validate(node, "event-Deletion.schema.json", "the node"),
                            ContainsSubstring("required property 'priority' not found"));
    }
}
