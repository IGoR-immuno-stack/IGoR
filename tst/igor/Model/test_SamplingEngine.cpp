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

TEST_CASE("SamplingEngine construction and factory builder", "[SamplingEngine]") {
    auto topology = std::make_shared<Topology>();
    
    // Create an event
    std::unordered_map<std::string, Event_realization> parent_realizations1 = create_realizations_se(3);
    auto event = std::make_shared<Gene_choice>(Undefined_gene, parent_realizations1);
    event->set_nickname("v_choice");
    
    topology->addEvent(event);
    
    SamplingEngine<double> engine(topology);
    
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
    auto topology = std::make_shared<Topology>();
    
    // Parent event (size 2)
    std::unordered_map<std::string, Event_realization> parent_realizations2 = create_realizations_se(2);
    auto parent = std::make_shared<Gene_choice>(Undefined_gene, parent_realizations2);
    parent->set_nickname("parent");
    igor::index_type parent_id = topology->addEvent(parent);
    
    // Child event (size 3)
    std::unordered_map<std::string, Event_realization> child_realizations1 = create_realizations_se(3);
    auto child = std::make_shared<Gene_choice>(Undefined_gene, child_realizations1);
    child->set_nickname("child");
    igor::index_type child_id = topology->addEvent(child);
    
    // Add edge
    topology->addEdge(parent_id, child_id);
    
    SamplingEngine<double> engine(topology);
    
    // Fetch & setup Parent Handler (Unconditional)
    {
        auto& h_base = engine.handler("parent");
        auto* h = dynamic_cast<CategoricalSamplingHandler<double>*>(&h_base);
        REQUIRE(h != nullptr);
        
        // Set probas: 100% prob for index 1
        h->parameters().data()[0] = 0.0;
        h->parameters().data()[1] = 1.0;
        h->precomputeCDF();
    }
    
    // Fetch & setup Child Handler (Conditional on parent)
    {
        auto& h_base = engine.handler("child");
        auto* h = dynamic_cast<CategoricalSamplingHandler<double>*>(&h_base);
        REQUIRE(h != nullptr);
        
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
    }
    
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
