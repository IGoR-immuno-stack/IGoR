#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <igor/Generation/SamplingEngine.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Topology.h>
#include <igor/Generation/CategoricalSamplingHandler.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Rec_Event.h>

#include <memory>
#include <vector>
#include <string>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
using namespace igor::generation;

using namespace igor::model;

// Helper to create a dummy Event_realization map for GeneChoice
static std::unordered_map<std::string, Event_realization> create_realizations_se(int n) {
    std::unordered_map<std::string, Event_realization> realizations;
    for (int i = 0; i < n; ++i) {
        std::string name = "gene_" + std::to_string(i);
        realizations.emplace(name, Event_realization(name, i, name, Int_Str(), i));
    }
    return realizations;
}

TEST_CASE("SamplingEngine construction and factory builder", "[SamplingEngine]") {
    auto topology = std::make_unique<Topology>();

    // Create an event
    std::unordered_map<std::string, Event_realization> parent_realizations1 = create_realizations_se(3);
    auto event = std::make_shared<Gene_choice>(Undefined_gene, parent_realizations1);
    event->set_nickname("v_choice");

    topology->addEvent(event);

    auto model = std::make_shared<RecombinationModel<double>>(std::move(topology));
    // Fill with uniform probabilities so precomputeCDF succeeds
    auto& w = model->weight(0);
    for (std::size_t i = 0; i < w.size(); ++i) w.data()[i] = 1.0 / w.size();

    SamplingEngine<double> engine(model);

    SECTION("Handlers are automatically built") {
        REQUIRE_NOTHROW(engine.handler("v_choice"));
        auto& h_base = engine.handler("v_choice");
        REQUIRE(h_base.name() == std::string("v_choice"));
        REQUIRE(h_base.uid() == 0);
    }

    SECTION("Accessing unknown event throws") {
        REQUIRE_THROWS_AS(engine.handler("unknown_event"), std::out_of_range);
    }
}

TEST_CASE("SamplingEngine generation (Categorical)", "[SamplingEngine]") {
    // Parent -> Child topology
    auto topology = std::make_unique<Topology>();

    // Parent event (size 2)
    std::unordered_map<std::string, Event_realization> parent_realizations2 = create_realizations_se(2);
    auto parent_ev = std::make_shared<Gene_choice>(Undefined_gene, parent_realizations2);
    parent_ev->set_nickname("parent");
    igor::core::legacy::index_type parent_id = topology->addEvent(parent_ev);

    // Child event (size 3)
    std::unordered_map<std::string, Event_realization> child_realizations1 = create_realizations_se(3);
    auto child_ev = std::make_shared<Gene_choice>(Undefined_gene, child_realizations1);
    child_ev->set_nickname("child");
    igor::core::legacy::index_type child_id = topology->addEvent(child_ev);

    // Add edge
    topology->addEdge(parent_id, child_id);

    // Build RecombinationModel and fill tensors before constructing engine
    auto model = std::make_shared<RecombinationModel<double>>(std::move(topology));

    // Parent tensor (shape {2}): 100% prob for index 1
    {
        auto& w = model->weight(parent_id);
        w.data()[0] = 0.0;
        w.data()[1] = 1.0;
    }

    // Child tensor (shape {3, 2}): conditional on parent
    {
        auto& w = model->weight(child_id);
        std::fill(w.data(), w.data() + w.size(), 0.0);

        // If parent=0: 100% prob for child=0
        w.data()[0 * 2 + 0] = 1.0;

        // If parent=1: 100% prob for child=2
        w.data()[2 * 2 + 1] = 1.0;
    }

    // Construct engine — precomputes CDFs from model tensors
    SamplingEngine<double> engine(model);

    std::mt19937_64 rng(42);
    auto scenario = engine.run(rng);

    // Verification
    // Parent should be 1
    REQUIRE(scenario.events.size() > parent_id);
    REQUIRE(scenario.events[parent_id].indices.size() == 1);
    REQUIRE(scenario.events[parent_id].indices[0] == 1);

    // Child should be 2 (conditioned on parent=1)
    REQUIRE(scenario.events.size() > child_id);
    REQUIRE(scenario.events[child_id].indices.size() == 1);
    REQUIRE(scenario.events[child_id].indices[0] == 2);
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
