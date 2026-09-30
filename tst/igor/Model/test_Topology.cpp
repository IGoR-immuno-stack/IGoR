/*
 * test_Topology.cpp
 *
 *  Created on: Feb 19, 2026
 *      Author: IGoR Agent
 */

#include <catch2/catch_test_macros.hpp>
#include <igor/Model/Topology.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Genechoice.h>


using namespace igor;
using namespace igor::model;

TEST_CASE("Topology Construction", "[Model][Topology]") {
    Topology topology;
    auto e1 = std::make_shared<Gene_choice>(); e1->set_nickname("NickE1");
    auto e2 = std::make_shared<Gene_choice>(); e2->set_nickname("NickE2");

    REQUIRE(topology.size() == 0);

    // Names are generated based on priority, so use dynamic names for verification
    std::string n1 = e1->get_nickname();
    std::string n2 = e2->get_nickname();

    index_type id1 = topology.addEvent(e1);
    index_type id2 = topology.addEvent(e2);

    REQUIRE(id1 == 0);
    REQUIRE(id2 == 1);
    REQUIRE(topology.size() == 2);
    
    REQUIRE(e1->uid() == 0);
    REQUIRE(e2->uid() == 1);

    REQUIRE(topology.hasEvent(n1));
    REQUIRE(topology.hasEvent(n2));
    REQUIRE_FALSE(topology.hasEvent("NonExistent"));
    
    REQUIRE(topology.event(0) == e1);
    REQUIRE(topology.event(n2) == e2);
    REQUIRE(topology.eventId(n1) == 0);
}

TEST_CASE("Topology Edges & Adjacency", "[Model][Topology]") {
    Topology topology;
    auto root = std::make_shared<Gene_choice>(); root->set_priority(10);
    auto child1 = std::make_shared<Gene_choice>(); child1->set_priority(11);
    auto child2 = std::make_shared<Gene_choice>(); child2->set_priority(12);
    
    index_type r_id = topology.addEvent(root);
    index_type c1_id = topology.addEvent(child1);
    index_type c2_id = topology.addEvent(child2);

    topology.addEdge(r_id, c1_id);
    topology.addEdge(r_id, c2_id);

    REQUIRE(topology.hasEdge(r_id, c1_id));
    REQUIRE(topology.hasEdge(r_id, c2_id));
    REQUIRE_FALSE(topology.hasEdge(c1_id, r_id));

    // Test Adjacency Iteration (View)
    auto children = topology.children(r_id);
    REQUIRE(children.size() == 2);
    // Order in children list matches insertion order? Yes, push_back in addEdge.
    REQUIRE(children[0] == child1);
    REQUIRE(children[1] == child2);
    
    int count = 0;
    for(const auto& child : children) {
        REQUIRE((child == child1 || child == child2));
        count++;
    }
    REQUIRE(count == 2);

    auto parents1 = topology.parents(c1_id);
    REQUIRE(parents1.size() == 1);
    REQUIRE(*parents1.begin() == root);
}

TEST_CASE("Topology Modification", "[Model][Topology]") {
    Topology topology;
    auto p = std::make_shared<Gene_choice>(); p->set_priority(20);
    auto c = std::make_shared<Gene_choice>(); c->set_priority(21);
    
    index_type p_id = topology.addEvent(p);
    index_type c_id = topology.addEvent(c);
    
    topology.addEdge(p_id, c_id);
    REQUIRE(topology.hasEdge(p_id, c_id));

    SECTION("Invert Edge") {
        topology.invertEdge(p_id, c_id);
        REQUIRE_FALSE(topology.hasEdge(p_id, c_id));
        REQUIRE(topology.hasEdge(c_id, p_id));
        
        auto new_roots = topology.roots();
        REQUIRE(new_roots.size() == 1);
        REQUIRE(new_roots[0] == c_id);
    }
    
    SECTION("Remove Edge") {
        topology.removeEdge(p_id, c_id);
        REQUIRE_FALSE(topology.hasEdge(p_id, c_id));
        REQUIRE(topology.children(p_id).empty());
        REQUIRE(topology.parents(c_id).empty());
    }
}

TEST_CASE("Topology Algorithms", "[Model][Topology]") {
    Topology t;
    auto n1 = std::make_shared<Gene_choice>(); n1->set_priority(30);
    auto n2 = std::make_shared<Gene_choice>(); n2->set_priority(31);
    auto n3 = std::make_shared<Gene_choice>(); n3->set_priority(32);
    
    t.addEvent(n1); // 0
    t.addEvent(n2); // 1
    t.addEvent(n3); // 2
    
    t.addEdge(0, 1);
    t.addEdge(1, 2);
    
    SECTION("Ancestors") {
        auto ans = t.ancestors(2); // 2 -> 1 -> 0
        REQUIRE(ans.size() == 2);
        
        bool has0 = false, has1 = false;
        for(auto id : ans) {
            if(id == 0) has0 = true;
            if(id == 1) has1 = true;
        }
        REQUIRE(has0);
        REQUIRE(has1);
    }
    
    SECTION("Roots") {
        auto r = t.roots();
        REQUIRE(r.size() == 1);
        REQUIRE(r[0] == 0);
    }
}

TEST_CASE("Topology read_topology", "[Model][Topology]") {
    std::string model_path = std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    auto topo = igor::model::read_topology(model_path);
    REQUIRE(topo != nullptr);
    REQUIRE(topo->size() > 0);
    
    REQUIRE(topo->hasEvent("v_choice"));
    REQUIRE(topo->hasEvent("vd_ins"));
    REQUIRE(topo->hasEvent("d_gene"));
    REQUIRE(topo->hasEvent("dj_ins"));
    REQUIRE(topo->hasEvent("j_choice"));
    
    auto v_id = topo->eventId("v_choice");
    auto root_events = topo->roots();
    bool v_is_root = false;
    for (auto id : root_events) if (id == v_id) v_is_root = true;
    REQUIRE(v_is_root);
}
