/*
 * test_RecombinationModel.cpp
 *
 *  Created on: Feb 21, 2026
 *      Author: IGoR Agent
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Model/RecombinationModel.h>
#include <igor/Model/SamplingEngine.h>
#include <igor/Model/Topology.h>
#include <igor/Model/LegacyBridge.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>

#include <memory>
#include <numeric>

using namespace igor;
using namespace igor::model;
using Catch::Matchers::WithinAbs;

// ─── Helpers ─────────────────────────────────────────────────────────────────

static std::shared_ptr<Gene_choice> make_gene_choice(const std::string& nick, int n)
{
    std::unordered_map<std::string, Event_realization> r;
    for (int i = 0; i < n; ++i) {
        std::string s = "r_" + std::to_string(i);
        r.emplace(s, Event_realization(s, i, s, Int_Str(), i));
    }
    auto ev = std::make_shared<Gene_choice>(Undefined_gene, r);
    ev->set_nickname(nick);
    return ev;
}

// ─── Construction ────────────────────────────────────────────────────────────

TEST_CASE("RecombinationModel construction from Topology",
          "[Model][RecombinationModel]")
{
    auto topo = std::make_unique<Topology>();

    // Unconditional event with 4 realizations
    auto ev_a = make_gene_choice("ev_a", 4);
    index_type a = topo->addEvent(ev_a);

    // Conditional event with 3 realizations, depends on ev_a
    auto ev_b = make_gene_choice("ev_b", 3);
    index_type b = topo->addEvent(ev_b);
    topo->addEdge(a, b);

    RecombinationModel<double> model(std::move(topo));

    SECTION("size matches topology") {
        REQUIRE(model.size() == 2);
    }

    SECTION("unconditional tensor shape is [n_realizations]") {
        const auto& w = model.weight(a);
        REQUIRE(w.ndim() == 1);
        REQUIRE(w.shape()[0] == 4);
        REQUIRE(w.size() == 4);
    }

    SECTION("conditional tensor shape is [child_realizations, parent_realizations]") {
        const auto& w = model.weight(b);
        REQUIRE(w.ndim() == 2);
        REQUIRE(w.shape()[0] == 3);   // ev_b has 3 realizations
        REQUIRE(w.shape()[1] == 4);   // parent ev_a has 4 realizations
        REQUIRE(w.size() == 12);
    }

    SECTION("tensors are zero-initialised") {
        for (std::size_t i = 0; i < model.weight(a).size(); ++i)
            REQUIRE(model.weight(a).data()[i] == 0.0);
    }

    SECTION("access by name") {
        REQUIRE_NOTHROW(model.weight("ev_a"));
        REQUIRE_NOTHROW(model.weight("ev_b"));
        REQUIRE(model.weight("ev_a").size() == 4);
    }

    SECTION("access invalid UID throws") {
        REQUIRE_THROWS_AS(model.weight(99), std::out_of_range);
    }

    SECTION("access unknown name throws") {
        REQUIRE_THROWS(model.weight("nonexistent"));
    }

    SECTION("null topology throws") {
        REQUIRE_THROWS_AS(
            RecombinationModel<double>(std::unique_ptr<Topology>()),
            std::invalid_argument);
    }
}

// ─── Read parameters ─────────────────────────────────────────────────────────

TEST_CASE("RecombinationModel read_parameters from file",
          "[Model][RecombinationModel]")
{
    // Build the same topology that the Mouse TCR beta model uses
    const std::string parms_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    const std::string marginals_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";

    Model_Parms parms;
    try { parms.read_model_parms(parms_path); }
    catch (...) { SKIP("Mouse TCR beta model_parms not found"); }

    auto topology = import_from_legacy(parms);

    // Convert shared_ptr → unique_ptr for RecombinationModel
    RecombinationModel<double> model(
        std::make_unique<Topology>(std::move(*topology)));

    SECTION("read_parameters fills tensors correctly") {
        REQUIRE(read_parameters(marginals_path, model));

        // Every tensor should have at least one non-zero value
        for (index_type uid = 0;
             uid < static_cast<index_type>(model.size()); ++uid)
        {
            const auto& w = model.weight(uid);
            INFO("Event: " << model.topology().event(uid)->get_nickname()
                           << "  uid=" << uid);

            bool any_nonzero = false;
            for (std::size_t i = 0; i < w.size(); ++i) {
                if (w.data()[i] != 0.0) { any_nonzero = true; break; }
            }
            REQUIRE(any_nonzero);
        }
    }

    SECTION("read_parameters file not found returns false") {
        REQUIRE_FALSE(read_parameters("/no/such/path.txt", model));
    }
}

// ─── Consistency: SamplingEngine handlers borrow model tensors ───────────────

TEST_CASE("SamplingEngine handlers reference RecombinationModel tensors",
          "[Model][RecombinationModel]")
{
    const std::string parms_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    const std::string marginals_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";

    RecombinationModel<double> model = [&]{
        try { return recombination_model_from_files<double>(parms_path, marginals_path); }
        catch (...) { SKIP("Mouse TCR beta model files not found"); throw; }
    }();

    auto model_ptr = std::make_shared<const RecombinationModel<double>>(std::move(model));
    SamplingEngine<double> engine(model_ptr);

    // Every handler's weights() must point to the exact same data as the model tensor
    for (index_type uid = 0;
         uid < static_cast<index_type>(model_ptr->topology().size()); ++uid)
    {
        const auto& tensor = model_ptr->weight(uid);
        const auto& handler = engine.handler(uid);

        INFO("Event: " << model_ptr->topology().event(uid)->get_nickname());
        REQUIRE(handler.weights().data() == tensor.data());
        REQUIRE(handler.weights().size() == tensor.size());
    }
}

// ─── Iteration ───────────────────────────────────────────────────────────────

TEST_CASE("RecombinationModel iteration over weights",
          "[Model][RecombinationModel]")
{
    auto topo = std::make_unique<Topology>();
    topo->addEvent(make_gene_choice("g1", 5));
    topo->addEvent(make_gene_choice("g2", 3));

    RecombinationModel<double> model(std::move(topo));

    std::size_t count = 0;
    for (const auto& tensor : model) {
        REQUIRE(tensor.size() > 0);
        ++count;
    }
    REQUIRE(count == 2);
}

// ─── Ordered traversal ──────────────────────────────────────────────────────

TEST_CASE("RecombinationModel orderedWeights iterates in topological order",
          "[Model][RecombinationModel]")
{
    //  Build a diamond DAG:  a ──▶ b ──▶ d
    //                        a ──▶ c ──▶ d
    //
    //  A valid topological order must visit a before {b,c} and {b,c} before d.

    auto topo = std::make_unique<Topology>();
    auto ev_a = make_gene_choice("ev_a", 4);
    auto ev_b = make_gene_choice("ev_b", 3);
    auto ev_c = make_gene_choice("ev_c", 2);
    auto ev_d = make_gene_choice("ev_d", 5);

    index_type a = topo->addEvent(ev_a);
    index_type b = topo->addEvent(ev_b);
    index_type c = topo->addEvent(ev_c);
    index_type d = topo->addEvent(ev_d);

    topo->addEdge(a, b);
    topo->addEdge(a, c);
    topo->addEdge(b, d);
    topo->addEdge(c, d);

    RecombinationModel<double> model(std::move(topo));

    // Tag every tensor so we can identify it by its size
    // ev_a(4), ev_b(3×4=12), ev_c(2×4=8), ev_d(5×3×2=30)

    SECTION("orderedWeights visits all tensors") {
        auto nav = model.orderedWeights();
        REQUIRE(nav.size() == 4);

        std::size_t count = 0;
        for (const auto& tensor : nav) {
            REQUIRE(tensor.size() > 0);
            ++count;
        }
        REQUIRE(count == 4);
    }

    SECTION("topological order constraints are respected") {
        auto nav = model.orderedWeights();

        // Record the sizes encountered in order
        std::vector<std::size_t> sizes;
        for (const auto& tensor : nav)
            sizes.push_back(tensor.size());

        // Find positions
        auto pos = [&](std::size_t sz) {
            return std::distance(sizes.begin(),
                                 std::find(sizes.begin(), sizes.end(), sz));
        };

        auto pos_a = pos(model.weight(a).size());
        auto pos_b = pos(model.weight(b).size());
        auto pos_c = pos(model.weight(c).size());
        auto pos_d = pos(model.weight(d).size());

        // a must come before b, c, d
        REQUIRE(pos_a < pos_b);
        REQUIRE(pos_a < pos_c);
        REQUIRE(pos_a < pos_d);
        // b and c must come before d
        REQUIRE(pos_b < pos_d);
        REQUIRE(pos_c < pos_d);
    }

    SECTION("random access via operator[]") {
        auto nav = model.orderedWeights();
        REQUIRE(nav.size() >= 1);
        // First element should be ev_a (root)
        REQUIRE(nav[0].size() == model.weight(a).size());
    }
}

// ─── Topology access ─────────────────────────────────────────────────────────

TEST_CASE("RecombinationModel owns topology",
          "[Model][RecombinationModel]")
{
    auto topo = std::make_unique<Topology>();
    auto* raw = topo.get();
    topo->addEvent(make_gene_choice("ev", 2));

    RecombinationModel<double> model(std::move(topo));

    // The model now owns the topology — verify reference identity
    REQUIRE(&model.topology() == raw);
    REQUIRE(model.topology().size() == 1);
}

// ─── Mutability ──────────────────────────────────────────────────────────────

TEST_CASE("RecombinationModel tensors are mutable",
          "[Model][RecombinationModel]")
{
    auto topo = std::make_unique<Topology>();
    topo->addEvent(make_gene_choice("my_event", 3));

    RecombinationModel<double> model(std::move(topo));
    auto& w = model.weight("my_event");

    // Fill with a uniform distribution
    for (std::size_t i = 0; i < w.size(); ++i)
        w.data()[i] = 1.0 / static_cast<double>(w.size());

    // Verify
    double sum = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i)
        sum += w.data()[i];
    REQUIRE_THAT(sum, WithinAbs(1.0, 1e-12));
}

// ─── Factory from files ─────────────────────────────────────────────────────

TEST_CASE("recombination_model_from_files loads a model in one step",
          "[Model][RecombinationModel]")
{
    const std::string parms_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    const std::string marginals_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";

    RecombinationModel<double> model = [&] {
        try { return recombination_model_from_files(parms_path, marginals_path); }
        catch (...) { SKIP("Mouse TCR beta model files not found"); throw; }
    }();

    SECTION("topology is populated") {
        REQUIRE(model.size() > 0);
        REQUIRE(model.topology().size() == model.size());
    }

    SECTION("all tensors have non-zero values") {
        for (index_type uid = 0;
             uid < static_cast<index_type>(model.size()); ++uid)
        {
            const auto& w = model.weight(uid);
            INFO("Event: " << model.topology().event(uid)->get_nickname());

            bool any_nonzero = false;
            for (std::size_t i = 0; i < w.size(); ++i) {
                if (w.data()[i] != 0.0) { any_nonzero = true; break; }
            }
            REQUIRE(any_nonzero);
        }
    }

    SECTION("matches manual two-step construction") {
        // Build independently via the two-step approach
        auto topology2 = read_topology(parms_path);
        REQUIRE(topology2);
        RecombinationModel<double> model2(
            std::make_unique<Topology>(std::move(*topology2)));
        REQUIRE(read_parameters(marginals_path, model2));

        REQUIRE(model.size() == model2.size());
        for (index_type uid = 0;
             uid < static_cast<index_type>(model.size()); ++uid)
        {
            const auto& w1 = model.weight(uid);
            const auto& w2 = model2.weight(uid);
            REQUIRE(w1.size() == w2.size());
            for (std::size_t i = 0; i < w1.size(); ++i)
                REQUIRE(w1.data()[i] == w2.data()[i]);
        }
    }

    SECTION("bad parms file throws") {
        REQUIRE_THROWS_AS(
            recombination_model_from_files("/no/such/file.txt", marginals_path),
            std::runtime_error);
    }

    SECTION("bad marginals file throws") {
        REQUIRE_THROWS_AS(
            recombination_model_from_files(parms_path, "/no/such/file.txt"),
            std::runtime_error);
    }
}
