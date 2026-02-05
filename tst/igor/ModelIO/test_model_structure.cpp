/*
 * test_model_structure.cpp
 *
 *  Created on: Feb 5, 2026
 *      Author: IGoR Development Team
 *
 *  Unit tests for ModelIO structure reading and writing
 *
 *  Uses Catch2 v3 features:
 *  - TEST_CASE_METHOD for fixtures
 *  - SECTION for grouping related tests
 *  - Matchers for better error messages
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>

#include "ModelIOTestUtils.h"
#include <igor/ModelIO/ModelReader.h>
#include <igor/ModelIO/ModelWriter.h>

#include <fstream>

using namespace igor::modelio;
using namespace igor::modelio::test;
using namespace Catch::Matchers;

//==============================================================================
// Fixtures
//==============================================================================

/**
 * @brief Fixture for ModelIO structure tests with automatic cleanup
 *
 * Creates test model files in a temporary directory that is
 * automatically cleaned up when the fixture is destroyed.
 */
struct ModelIOStructureFixture
{
    TestDirectory test_dir{"igor_modelio_structure_tests"};

    std::string path(const std::string& name) { return test_dir.path(name); }

    /**
     * @brief Create a minimal JSON model file
     */
    std::string create_minimal_json()
    {
        auto filepath = path("minimal_model.json");
        ModelData model = create_minimal_model();
        structure::write_json(filepath, model);
        return filepath;
    }

    /**
     * @brief Create a JSON model file with realizations
     */
    std::string create_json_with_realizations(size_t count = 3)
    {
        auto filepath = path("model_with_realizations.json");
        ModelData model = create_model_with_realizations(count);
        structure::write_json(filepath, model);
        return filepath;
    }

    /**
     * @brief Create a complete VDJ model JSON file
     */
    std::string create_vdj_json()
    {
        auto filepath = path("vdj_model.json");
        ModelData model = create_vdj_model();
        structure::write_json(filepath, model);
        return filepath;
    }

    /**
     * @brief Create an invalid JSON file for error testing
     */
    std::string create_invalid_json()
    {
        auto filepath = path("invalid.json");
        std::ofstream file(filepath);
        file << "{ invalid json }";
        file.close();
        return filepath;
    }
};

//==============================================================================
// Basic Structure Tests
//==============================================================================
using namespace igor::modelio::test;
using namespace Catch::Matchers;

//==============================================================================
// Basic Structure Tests
//==============================================================================

TEST_CASE("ModelData: Basic structure creation", "[modelio][structure]")
{
    ModelData model = create_minimal_model();

    REQUIRE(model.format_version == FORMAT_VERSION);
    REQUIRE(model.metadata.species == "Homo sapiens");
    REQUIRE(model.metadata.chain == "TRB");
    REQUIRE_THAT(model.events, SizeIs(1));
    REQUIRE(model.events[0].id == "v_choice");
}

TEST_CASE("EventData: GeneChoice event creation", "[modelio][structure]")
{
    ModelData model = create_model_with_realizations(3);

    REQUIRE_THAT(model.events, SizeIs(1));
    const auto& event = model.events[0];

    REQUIRE(event.id == "v_choice");
    REQUIRE(event.type == "GeneChoice");
    REQUIRE_THAT(event.realizations, SizeIs(3));
    REQUIRE(event.realizations[0].name == "TRBV1*01");
    REQUIRE(event.realizations[1].name == "TRBV2*01");
    REQUIRE(event.realizations[2].name == "TRBV3*01");
}

//==============================================================================
// JSON I/O Tests
//==============================================================================

TEST_CASE_METHOD(ModelIOStructureFixture, "JSON: Write and read minimal model", "[modelio][structure][json]")
{
    SECTION("minimal model round-trip")
    {
        ModelData original = create_minimal_model();
        auto filepath = path("minimal.json");

        structure::write_json(filepath, original);
        REQUIRE(std::filesystem::exists(filepath));

        ModelData loaded = structure::read_json(filepath);

        REQUIRE(loaded.format_version == FORMAT_VERSION);
        REQUIRE(loaded.metadata.chain == "TRB");
        REQUIRE_THAT(loaded.events, SizeIs(1));
        REQUIRE(loaded.events[0].id == "v_choice");
        REQUIRE_THAT(loaded.edges, SizeIs(1));
        REQUIRE(loaded.edges[0].parent == "v_choice");
    }
}

TEST_CASE_METHOD(ModelIOStructureFixture, "JSON: Round-trip with realizations", "[modelio][structure][json]")
{
    SECTION("model with 3 realizations")
    {
        ModelData original = create_model_with_realizations(3);
        auto filepath = path("with_realizations.json");

        structure::write_json(filepath, original);
        ModelData loaded = structure::read_json(filepath);

        REQUIRE_THAT(loaded.events, SizeIs(1));
        REQUIRE_THAT(loaded.events[0].realizations, SizeIs(3));
        REQUIRE(loaded.events[0].realizations[0].name == "TRBV1*01");
        REQUIRE(loaded.events[0].realizations[1].name == "TRBV2*01");
        REQUIRE(loaded.events[0].realizations[2].name == "TRBV3*01");
    }

    SECTION("model with 10 realizations")
    {
        ModelData original = create_model_with_realizations(10);
        auto filepath = path("with_10_realizations.json");

        structure::write_json(filepath, original);
        ModelData loaded = structure::read_json(filepath);

        REQUIRE_THAT(loaded.events[0].realizations, SizeIs(10));
        REQUIRE(loaded.events[0].realizations[9].name == "TRBV10*01");
    }
}

TEST_CASE_METHOD(ModelIOStructureFixture, "JSON: Complete VDJ model", "[modelio][structure][json]")
{
    SECTION("full VDJ model round-trip")
    {
        ModelData original = create_vdj_model();
        auto filepath = path("vdj_complete.json");

        structure::write_json(filepath, original);
        ModelData loaded = structure::read_json(filepath);

        REQUIRE(loaded.metadata.model_type == "VDJ");
        REQUIRE_THAT(loaded.events, SizeIs(9)); // 3 gene choices + 4 deletions + 2 insertions
        REQUIRE_THAT(loaded.edges, SizeIs(9));
        REQUIRE_THAT(loaded.sequence_types.order, SizeIs(5));

        // Verify gene choice events have realizations
        REQUIRE_THAT(loaded.events[0].realizations, SizeIs(3)); // V genes
        REQUIRE_THAT(loaded.events[1].realizations, SizeIs(2)); // D genes
        REQUIRE_THAT(loaded.events[2].realizations, SizeIs(2)); // J genes
    }
}

TEST_CASE_METHOD(ModelIOStructureFixture, "JSON: Pretty printing", "[modelio][structure][json]")
{
    SECTION("pretty-printed output is readable")
    {
        ModelData model = create_minimal_model();
        auto filepath = path("pretty.json");

        structure::write_json(filepath, model, true);

        // Read file content
        std::ifstream file(filepath);
        std::string content((std::istreambuf_iterator<char>(file)),
                             std::istreambuf_iterator<char>());

        // Should contain newlines and proper indentation
        REQUIRE_THAT(content, ContainsSubstring("\n"));
        REQUIRE_THAT(content, ContainsSubstring("format_version"));
        REQUIRE_THAT(content, ContainsSubstring("  ")); // Indentation
    }
}

TEST_CASE_METHOD(ModelIOStructureFixture, "JSON: Sequence types preservation", "[modelio][structure][json]")
{
    SECTION("sequence type definitions round-trip")
    {
        ModelData model = create_vdj_model();
        auto filepath = path("seq_types.json");

        structure::write_json(filepath, model);
        ModelData loaded = structure::read_json(filepath);

        REQUIRE_THAT(loaded.sequence_types.order, SizeIs(5));
        REQUIRE(loaded.sequence_types.order[0] == "V_gene");
        REQUIRE(loaded.sequence_types.order[4] == "J_gene");

        REQUIRE(loaded.sequence_types.definitions.count("V_gene") == 1);
        REQUIRE(loaded.sequence_types.definitions["V_gene"].id == 0);

        REQUIRE(loaded.sequence_types.definitions.count("VD_ins") == 1);
        REQUIRE_THAT(loaded.sequence_types.definitions["VD_ins"].parents, SizeIs(1));
        REQUIRE(loaded.sequence_types.definitions["VD_ins"].parents[0] == "V_gene");
    }
}

//==============================================================================
// Error Handling Tests
//==============================================================================

TEST_CASE_METHOD(ModelIOStructureFixture, "Error handling: Invalid JSON", "[modelio][structure][error]")
{
    SECTION("malformed JSON throws exception")
    {
        auto filepath = create_invalid_json();
        REQUIRE_THROWS_AS(structure::read_json(filepath), std::runtime_error);
    }

    SECTION("error message contains useful information")
    {
        auto filepath = create_invalid_json();
        try {
            structure::read_json(filepath);
            FAIL("Expected std::runtime_error");
        } catch (const std::runtime_error& e) {
            std::string msg = e.what();
            REQUIRE_THAT(msg, ContainsSubstring("JSON parse error"));
        }
    }
}

TEST_CASE("Error handling: Non-existent file", "[modelio][structure][error]")
{
    SECTION("reading non-existent file throws exception")
    {
        REQUIRE_THROWS_AS(
            structure::read_json("/nonexistent/path/model.json"),
            std::runtime_error
        );
    }

    SECTION("error message is informative")
    {
        try {
            structure::read_json("/nonexistent/path/model.json");
            FAIL("Expected std::runtime_error");
        } catch (const std::runtime_error& e) {
            std::string msg = e.what();
            REQUIRE_THAT(msg, ContainsSubstring("Cannot open"));
        }
    }
}

//==============================================================================
// Comparison Helper Tests
//==============================================================================

TEST_CASE("Comparison: models_equal utility", "[modelio][structure][utility]")
{
    SECTION("identical models are equal")
    {
        ModelData model1 = create_minimal_model();
        ModelData model2 = create_minimal_model();
        REQUIRE(models_equal(model1, model2));
    }

    SECTION("different metadata makes models unequal")
    {
        ModelData model1 = create_minimal_model();
        ModelData model2 = create_minimal_model();
        model2.metadata.chain = "TRA";
        REQUIRE_FALSE(models_equal(model1, model2));
    }

    SECTION("different event count makes models unequal")
    {
        ModelData model1 = create_minimal_model();
        ModelData model2 = create_vdj_model();
        REQUIRE_FALSE(models_equal(model1, model2));
    }
}
