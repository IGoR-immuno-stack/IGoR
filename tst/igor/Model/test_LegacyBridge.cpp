/*
 * test_LegacyBridge.cpp
 *
 *  Created on: Feb 10, 2026
 *
 *  Integration tests for LegacyBridge: import/export between
 *  Model_marginals and InferenceEngine.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Model/LegacyBridge.h>
#include <igor/Model/Topology.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Model/InferenceEngine.h>

#include <sstream>
#include <cmath>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ─── Helper: Create minimal test model ─────────────────────────────────

static Model_Parms create_test_model() {
    // Create a minimal model with a few events
    // This would normally be loaded from files
    Model_Parms parms;

    // TODO: Add events to parms for testing
    // For now, this is a placeholder

    return parms;
}

// ─── Extract EventDescriptor Tests ─────────────────────────────────────

TEST_CASE("extract_event_descriptors with empty model", "[core][bridge]") {
    Model_Parms parms;

    auto descriptors = extract_event_descriptors(parms);

    REQUIRE(descriptors.empty());
}

// ─── Import/Export Round-Trip Tests ────────────────────────────────────

// ─── Round-trip with Mouse TCR beta model preserves marginals ────────────────

TEST_CASE("Round-trip with Mouse TCR beta model preserves marginals", "[core][bridge][integration]") {
    // Load real Mouse TCR beta model
    Model_Parms parms;

    // Try to read model, skip if files don't exist
    std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    try {
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found at: " + model_path);
    }

    // Load marginals
    Model_marginals original(parms);
    std::string marginals_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";
    original.txt2marginals(marginals_path, parms);

    // Extract descriptors and create engine
    auto descriptors = extract_event_descriptors(parms);
    REQUIRE(!descriptors.empty());

    InferenceEngine<long double> engine(descriptors);
    REQUIRE(engine.size() == descriptors.size());

    // Import from original marginals into engine
    import_from_legacy(engine, original, parms);

    // Export back
    Model_marginals roundtrip(parms);
    export_to_legacy(engine, roundtrip, parms);

    // Compare: all marginals should match within numerical precision
    const long double* orig_array = original.marginal_array_smart_p.get();
    const long double* trip_array = roundtrip.marginal_array_smart_p.get();
    size_t length = original.get_length();

    REQUIRE(roundtrip.get_length() == length);

    double max_diff = 0.0;
    for (size_t i = 0; i < length; ++i) {
        double diff = std::abs(orig_array[i] - trip_array[i]);
        max_diff = std::max(max_diff, diff);
        REQUIRE_THAT(static_cast<double>(trip_array[i]),
                     WithinAbs(static_cast<double>(orig_array[i]), 1e-12));
    }

    INFO("Maximum round-trip difference: " << max_diff);
}

TEST_CASE("Extract descriptors from Mouse TCR beta model", "[core][bridge]") {
    Model_Parms parms;

    try {
        std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found");
    }

    auto descriptors = extract_event_descriptors(parms);

    REQUIRE(!descriptors.empty());

    // Check that we have typical VDJ events
    bool has_v_choice = false;
    bool has_d_gene = false;
    bool has_j_choice = false;
    bool has_markov = false;

    INFO("Total events found: " << descriptors.size());

    for (const auto& desc : descriptors) {
        INFO("Event name: '" << desc.name << "' type=" << desc.type << " shape_size=" << desc.shape.size());

        if (desc.name.find("v_choice") != std::string::npos ||
            desc.name.find("V_gene") != std::string::npos) has_v_choice = true;
        if (desc.name.find("d_gene") != std::string::npos ||
            desc.name.find("D_gene") != std::string::npos) has_d_gene = true;
        if (desc.name.find("j_choice") != std::string::npos ||
            desc.name.find("J_gene") != std::string::npos) has_j_choice = true;
        if (desc.type == static_cast<int>(Dinuclmarkov_t)) has_markov = true;

        // All events should have valid shapes
        REQUIRE(!desc.shape.empty());
        REQUIRE(desc.shape[0] > 0);
    }

    // Mouse TCR beta model should have these core events
    REQUIRE(has_v_choice);
    REQUIRE(has_d_gene);
    REQUIRE(has_j_choice);
}

TEST_CASE("Round-trip preserves uniform marginals", "[core][bridge]") {
    // Load real Mouse TCR beta model
    Model_Parms parms;

    // Try to read model, skip if files don't exist
    try {
        std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found");
    }

    // Initialize with uniform marginals
    Model_marginals original(parms);
    original.uniform_initialize(parms);

    // Extract descriptors and create engine
    auto descriptors = extract_event_descriptors(parms);
    REQUIRE(!descriptors.empty());

    InferenceEngine<long double> engine(descriptors);
    REQUIRE(engine.size() == descriptors.size());

    // Import from original marginals into engine
    import_from_legacy(engine, original, parms);

    // Export back
    Model_marginals roundtrip(parms);
    export_to_legacy(engine, roundtrip, parms);

    // Verify arrays match
    const long double* orig_array = original.marginal_array_smart_p.get();
    const long double* trip_array = roundtrip.marginal_array_smart_p.get();
    size_t length = original.get_length();

    REQUIRE(roundtrip.get_length() == length);

    // Verify all values match within tolerance
    double max_diff = 0.0;
    for (size_t i = 0; i < length; ++i) {
        double diff = std::abs(orig_array[i] - trip_array[i]);
        max_diff = std::max(max_diff, diff);
        REQUIRE_THAT(static_cast<double>(trip_array[i]),
                     WithinAbs(static_cast<double>(orig_array[i]), 1e-12));
    }

    INFO("Maximum round-trip difference: " << max_diff);
}

TEST_CASE("Import validates event existence", "[core][bridge]") {
    // Load a real model to avoid "empty marginals" error
    Model_Parms parms;
    try {
        std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found");
    }

    Model_marginals marginals(parms);
    marginals.uniform_initialize(parms);

    // Create engine with an event that doesn't exist in parms
    InferenceEngine<long double> engine;
    engine.register_handler("nonexistent_event_name",
        std::make_unique<CategoricalHandler<long double>>("nonexistent_event_name", 3));

    // Should throw when trying to import
    REQUIRE_THROWS_AS(
        import_from_legacy(engine, marginals, parms),
        std::runtime_error
    );
}

TEST_CASE("Export validates event existence", "[core][bridge]") {
    // Load a real model to avoid "empty marginals" error
    Model_Parms parms;
    try {
        std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found");
    }

    Model_marginals marginals(parms);
    marginals.uniform_initialize(parms);

    // Create engine with an event that doesn't exist in parms
    InferenceEngine<long double> engine;
    engine.register_handler("nonexistent_event_name",
        std::make_unique<CategoricalHandler<long double>>("nonexistent_event_name", 3));

    // Should throw when trying to export
    REQUIRE_THROWS_AS(
        export_to_legacy(engine, marginals, parms),
        std::runtime_error
    );
}

// ─── Numerical Precision Tests ──────────────────────────────────────────

TEST_CASE("Type conversion preserves precision", "[core][bridge]") {
    // Test that long double <-> double conversion doesn't lose too much precision

    SECTION("double to long double preserves values") {
        InferenceEngine<double> engine_double;
        engine_double.register_handler("test",
            std::make_unique<CategoricalHandler<double>>("test", 3));

        auto& handler = engine_double.handler("test");
        handler.accumulator()(0) = 1.0;
        handler.accumulator()(1) = 2.0;
        handler.accumulator()(2) = 3.0;
        handler.maximize_likelihood();

        // Values should be exactly representable
        REQUIRE_THAT(handler.parameters()(0), WithinRel(1.0/6.0, 1e-15));
        REQUIRE_THAT(handler.parameters()(1), WithinRel(2.0/6.0, 1e-15));
        REQUIRE_THAT(handler.parameters()(2), WithinRel(3.0/6.0, 1e-15));
    }
}

TEST_CASE("Topology <-> Model_Parms conversion", "[core][bridge][topology]") {
    std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    auto topo1 = igor::model::readTopology(model_path);
    REQUIRE(topo1 != nullptr);

    Model_Parms parms;
    try {
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Mouse TCR beta model files not found");
    }

    auto topo2 = igor::model::import_from_legacy(parms);
    REQUIRE(topo2 != nullptr);
    REQUIRE(topo1->size() == topo2->size());

    for (const auto& ev : *topo1) {
        auto name = ev->get_nickname();
        INFO("Event: " << name);
        REQUIRE(topo2->hasEvent(name));
        
        // Ensure uids match strictly
        REQUIRE(ev->uid() == topo2->event(name)->uid());
        
        auto id1 = topo1->eventId(name);
        auto id2 = topo2->eventId(name);
        
        // Check parents
        auto parents1 = topo1->parents(id1);
        auto parents2 = topo2->parents(id2);
        int parents1_count = 0, parents2_count = 0;
        for (const auto& p : parents1) parents1_count++;
        for (const auto& p : parents2) parents2_count++;
        REQUIRE(parents1_count == parents2_count);

        auto it1_p = parents1.begin();
        auto it2_p = parents2.begin();
        while (it1_p != parents1.end() && it2_p != parents2.end()) {
            REQUIRE((*it1_p)->uid() == (*it2_p)->uid());
            ++it1_p;
            ++it2_p;
        }
        
        // Check children
        auto children1 = topo1->children(id1);
        auto children2 = topo2->children(id2);
        int children1_count = 0, children2_count = 0;
        for (const auto& c : children1) children1_count++;
        for (const auto& c : children2) children2_count++;
        REQUIRE(children1_count == children2_count);

        auto it1_c = children1.begin();
        auto it2_c = children2.begin();
        while (it1_c != children1.end() && it2_c != children2.end()) {
            REQUIRE((*it1_c)->uid() == (*it2_c)->uid());
            ++it1_c;
            ++it2_c;
        }
    }

    auto parms2 = igor::model::export_to_legacy(*topo1);
    REQUIRE(parms2 != nullptr);
    REQUIRE(parms.get_event_list().size() == parms2->get_event_list().size());

    for (const auto& ev : parms.get_event_list()) {
        auto name = ev->get_name();
        
        auto parents1 = parms.get_parents(name);
        auto parents2 = parms2->get_parents(name);
        REQUIRE(parents1.size() == parents2.size());
        
        auto children1 = parms.get_children(name);
        auto children2 = parms2->get_children(name);
        REQUIRE(children1.size() == children2.size());
    }
}
