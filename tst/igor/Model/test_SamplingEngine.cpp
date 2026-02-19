#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <igor/Model/SamplingEngine.h>
#include <igor/Model/Topology.h>
#include <igor/Model/CategoricalSamplingHandler.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Rec_Event.h>

#include <memory>
#include <vector>
#include <string>

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

TEST_CASE("SamplingEngine construction and registration", "[SamplingEngine]") {
    auto topology = std::make_shared<Topology>();
    
    // Create an event
    auto event = std::make_shared<Gene_choice>(V_gene, GeneChoice_t, create_realizations_se(3));
    event->set_nickname("v_choice");
    
    topology->addEvent(event);
    
    SamplingEngine<double> engine(topology);
    
    SECTION("Register valid handler") {
        auto handler = std::make_unique<CategoricalSamplingHandler<double>>("v_choice", -1, std::vector<std::size_t>{3});
        handler->precomputeCDF();
        REQUIRE_NOTHROW(engine.registerHandler(std::move(handler)));
        REQUIRE(engine.handler("v_choice").name() == std::string("CategoricalSamplingHandler_") + typeid(double).name());
    }
    
    SECTION("Register handler for unknown event throws") {
        auto handler = std::make_unique<CategoricalSamplingHandler<double>>("unknown_event", -1, std::vector<std::size_t>{3});
        REQUIRE_THROWS_AS(engine.registerHandler(std::move(handler)), std::invalid_argument);
    }
    
    SECTION("Duplicate registration throws") {
        auto h1 = std::make_unique<CategoricalSamplingHandler<double>>("v_choice", -1, std::vector<std::size_t>{3});
        h1->precomputeCDF();
        engine.registerHandler(std::move(h1));
        
        auto h2 = std::make_unique<CategoricalSamplingHandler<double>>("v_choice", -1, std::vector<std::size_t>{3});
        REQUIRE_THROWS_AS(engine.registerHandler(std::move(h2)), std::logic_error);
    }
}

TEST_CASE("SamplingEngine generation (Categorical)", "[SamplingEngine]") {
    // Parent -> Child topology
    auto topology = std::make_shared<Topology>();
    
    // Parent event (size 2)
    auto parent = std::make_shared<Gene_choice>(V_gene, GeneChoice_t, create_realizations_se(2));
    parent->set_nickname("parent");
    igor::index_type parent_id = topology->addEvent(parent);
    
    // Child event (size 3)
    auto child = std::make_shared<Gene_choice>(J_gene, GeneChoice_t, create_realizations_se(3));
    child->set_nickname("child");
    igor::index_type child_id = topology->addEvent(child);
    
    // Add edge
    topology->addEdge(parent_id, child_id);
    
    SamplingEngine<double> engine(topology);
    
    // Create & Register Parent Handler (Unconditional)
    {
        auto h = std::make_unique<CategoricalSamplingHandler<double>>("parent", -1, std::vector<std::size_t>{2});
        // Set probas: 100% prob for index 1
        h->parameters().data()[0] = 0.0;
        h->parameters().data()[1] = 1.0;
        h->precomputeCDF();
        engine.registerHandler(std::move(h));
    }
    
    // Create & Register Child Handler (Conditional on parent)
    {
        // Shape: {3, 2} (3 realizations, 1 parent of size 2)
        auto h = std::make_unique<CategoricalSamplingHandler<double>>("child", -1, std::vector<std::size_t>{3, 2});
        
        // If parent=0: 100% prob for child=0
        // If parent=1: 100% prob for child=2
        auto& p = h->parameters(); // size 6
        // Flat index: child * n_parents + parent
        // child=0, parent=0 -> idx 0 (val 1.0)
        // child=1, parent=0 -> idx 2 (val 0.0)
        // child=2, parent=0 -> idx 4 (val 0.0)
        
        // child=0, parent=1 -> idx 1 (val 0.0)
        // child=1, parent=1 -> idx 3 (val 0.0)
        // child=2, parent=1 -> idx 5 (val 1.0)
        
        // Using correct indexing: 
        // offset = child * product(parent_dims) + parent_offset
        // here parent_dims=2.
        
        std::fill(p.data(), p.data() + 6, 0.0);
        
        // Case parent=0
        p.data()[0 * 2 + 0] = 1.0; 
        
        // Case parent=1
        p.data()[2 * 2 + 1] = 1.0;
        
        h->precomputeCDF();
        engine.registerHandler(std::move(h));
    }
    
    std::mt19937_64 rng(42);
    auto scenario = engine.generate(rng);
    
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
