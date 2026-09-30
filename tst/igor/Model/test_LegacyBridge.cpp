/*
 * test_LegacyBridge.cpp
 *
 *  Created on: Feb 10, 2026
 *  Updated on: Feb 21, 2026
 *
 *  Integration tests for LegacyBridge: import/export between
 *  Model_marginals and RecombinationModel, and Topology <-> Model_Parms.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Model/LegacyBridge.h>
#include <igor/Model/Topology.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Model_Parms.h>

#include <sstream>
#include <cmath>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// Topology <-> Model_Parms conversion

TEST_CASE("Topology <-> Model_Parms conversion", "[core][bridge][topology]") {
    std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    auto topo1 = igor::model::read_topology(model_path);
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
        REQUIRE(ev->uid() == topo2->event(name)->uid());

        auto id1 = topo1->eventId(name);
        auto id2 = topo2->eventId(name);

        auto parents1 = topo1->parents(id1);
        auto parents2 = topo2->parents(id2);
        int parents1_count = 0, parents2_count = 0;
        for (const auto& p : parents1) { (void)p; parents1_count++; }
        for (const auto& p : parents2) { (void)p; parents2_count++; }
        REQUIRE(parents1_count == parents2_count);

        auto children1 = topo1->children(id1);
        auto children2 = topo2->children(id2);
        int children1_count = 0, children2_count = 0;
        for (const auto& c : children1) { (void)c; children1_count++; }
        for (const auto& c : children2) { (void)c; children2_count++; }
        REQUIRE(children1_count == children2_count);
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

// RecombinationModel: read_parameters vs import_from_legacy

TEST_CASE("read_parameters matches import_from_legacy for Mouse TCR beta",
          "[core][bridge][RecombinationModel][integration]")
{
    const std::string parms_path     = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    const std::string marginals_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";

    igor::model::RecombinationModel<double> model_direct = [&]{
        try { return igor::model::recombination_model_from_files<double>(parms_path, marginals_path); }
        catch (...) { SKIP("Mouse TCR beta model files not found at: " + parms_path); throw; }
    }();

    Model_Parms parms;
    try { parms.read_model_parms(parms_path); }
    catch (...) { SKIP("Mouse TCR beta model_parms not found at: " + parms_path); }

    Model_marginals marginals(parms);
    try { marginals.txt2marginals(marginals_path, parms); }
    catch (...) { SKIP("Mouse TCR beta marginals not found at: " + marginals_path); }

    auto topology_legacy = igor::model::read_topology(parms_path);
    REQUIRE(topology_legacy);

    igor::model::RecombinationModel<double> model_legacy(
        std::make_unique<igor::model::Topology>(*topology_legacy));
    import_from_legacy(model_legacy, marginals);

    REQUIRE(model_direct.topology().size() > 0);
    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(model_direct.topology().size());
         ++uid)
    {
        const std::string name = model_direct.topology().event(uid)->get_nickname();
        INFO("Event: " << name << "  uid=" << uid);

        const auto& t_dir = model_direct.weight(uid);
        const auto& t_leg = model_legacy.weight(uid);

        REQUIRE(t_dir.size() == t_leg.size());

        for (std::size_t i = 0; i < t_dir.size(); ++i) {
            INFO("  index=" << i);
            REQUIRE_THAT(t_dir.data()[i], WithinAbs(t_leg.data()[i], 1e-6));
        }
    }
}
