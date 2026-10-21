/**
 * @file test_event_names.cpp
 * @brief What an event is called, and by whom (plan T2).
 *
 * The v1 name (get_v1_legacy_name()) is type, gene class, side, priority and size: what a v1
 * model file names an event by in its @Edges section, and the `legacy` output name. It tells two
 * events apart only where a v1 gene class does, so it is used at the v1 boundary only.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "ScratchDir.h"

#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Model_Parms.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace igor::core::legacy;
using namespace igor::model::legacy;

namespace {

const std::string kTestDataDir = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/format_v2/";
const std::string kHumanTrb = std::string(IGOR_MODELS_DIR) + "/human/tcr_beta/models/model_parms.txt";

std::string read_file(const std::string &path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/// Each event's children, by nickname: what a model's graph is, whatever its events are called.
using Graph = std::map<std::string, std::set<std::string>>;

Graph children_by_nickname(const Model_Parms &parms)
{
    Graph graph;
    for (const auto &event : parms.get_event_list()) {
        auto &children = graph[event->get_nickname()];
        for (const auto &child : parms.get_children(event)) {
            children.insert(child->get_nickname());
        }
    }
    return graph;
}

/// The same graph read from the other end, so that both adjacency lists are checked.
Graph parents_by_nickname(const Model_Parms &parms)
{
    Graph graph;
    for (const auto &event : parms.get_event_list()) {
        auto &parents = graph[event->get_nickname()];
        for (const auto &parent : parms.get_parents(event)) {
            parents.insert(parent->get_nickname());
        }
    }
    return graph;
}

Graph inverted(const Graph &children)
{
    Graph parents;
    for (const auto &[parent, its_children] : children) {
        parents[parent];
        for (const auto &child : its_children) {
            parents[child].insert(parent);
        }
    }
    return parents;
}

/// Event by event, so that a failure says which event lost or gained an edge.
void require_same(const Graph &actual, const Graph &expected, const char *relation)
{
    REQUIRE(actual.size() == expected.size());
    for (const auto &[nickname, events] : expected) {
        INFO(relation << " of " << nickname);
        REQUIRE(actual.count(nickname) == 1);
        CHECK(actual.at(nickname) == events);
    }
}

void require_graph(const Model_Parms &parms, const Graph &expected_children)
{
    require_same(children_by_nickname(parms), expected_children, "children");
    require_same(parents_by_nickname(parms), inverted(expected_children), "parents");
}

bool has_an_edge(const Graph &graph)
{
    return std::any_of(graph.begin(), graph.end(), [](const auto &entry) { return !entry.second.empty(); });
}

}  // namespace

TEST_CASE("A v1 file names each event of its edges by the event's v1 name", "[event_name][legacy]")
{
    const std::string source = GENERATE_COPY(kTestDataDir + "test_legacy_model_parms.txt", kHumanTrb);
    INFO("model " << source);
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(source));

    std::map<std::string, std::string> nickname_by_v1_name;
    for (const auto &event : parms.get_event_list()) {
        REQUIRE(nickname_by_v1_name.emplace(event->get_v1_legacy_name(), event->get_nickname()).second);
    }

    std::istringstream file(read_file(source));
    std::string line;
    while (std::getline(file, line) and line != "@Edges") {
    }
    Graph from_file;
    while (std::getline(file, line) and line[0] == '%') {
        const auto semicolon = line.find(';');
        const std::string parent = line.substr(1, semicolon - 1);
        const std::string child = line.substr(semicolon + 1);
        INFO("edge " << line);
        REQUIRE(nickname_by_v1_name.count(parent) == 1);
        REQUIRE(nickname_by_v1_name.count(child) == 1);
        from_file[nickname_by_v1_name.at(parent)].insert(nickname_by_v1_name.at(child));
    }
    REQUIRE(has_an_edge(from_file));
    for (const auto &event : parms.get_event_list()) {
        from_file[event->get_nickname()];
    }
    require_graph(parms, from_file);
}

TEST_CASE("A legacy model's edges survive a round trip", "[event_name][model_format]")
{
    const std::string source = GENERATE_COPY(kTestDataDir + "test_legacy_model_parms.txt", kHumanTrb);
    INFO("model " << source);
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(source));
    const Graph expected = children_by_nickname(parms);
    REQUIRE(has_an_edge(expected));

    ScratchDir scratch;
    const std::string written = scratch.file("written_parms.txt");
    SECTION("in v1") { parms.write_model_parms_legacy(written); }
    SECTION("in v2, whose edges name events by their v2 name") { parms.write_model_parms_v2(written); }

    Model_Parms read_back;
    REQUIRE_NOTHROW(read_back.read_model_parms(written));
    require_graph(read_back, expected);
}

TEST_CASE("A v1 edge token that names two events is refused", "[event_name][legacy]")
{
    using Catch::Matchers::ContainsSubstring;
    // Two V deletions alike but for their nickname: one v1 name, which the fixture's edge from
    // v_choice uses.
    std::string text = read_file(kTestDataDir + "test_legacy_model_parms.txt");
    const std::string line = "#Deletion;V_gene;Three_prime;5;v_3_del\n";
    const auto at = text.find(line);
    REQUIRE(at != std::string::npos);
    const auto block_end = text.find('#', at + 1);
    const std::string block = text.substr(at, block_end - at);
    std::string twin = block;
    twin.replace(0, line.size(), "#Deletion;V_gene;Three_prime;5;v_3_del_twin\n");
    text.insert(block_end, twin);

    ScratchDir scratch;
    Model_Parms parms;
    REQUIRE_THROWS_WITH(parms.read_model_parms(scratch.file("twin_parms.txt", text)),
                        ContainsSubstring("names two events"));
}

TEST_CASE("A model whose events the v1 format cannot name needs the extended format",
          "[event_name][model_format]")
{
    Model_Parms parms;
    const auto add_v_choice = [&] {
        auto v_choice = std::make_shared<Gene_choice>(V_gene);
        v_choice->add_realization("V1", "ACGT");
        v_choice->set_seq_type("V_gene_seq");
        v_choice->set_nickname("v_choice");
        parms.add_event(v_choice);
    };

    SECTION("a legacy model does not") {
        add_v_choice();
        CHECK_FALSE(parms.requires_extended_format());
    }
    SECTION("an event on a seq type its v1 line cannot name: a deletion on a junction") {
        add_v_choice();
        auto deletion = std::make_shared<Deletion>("VD_ins_seq", Five_prime, std::make_pair(0, 2));
        deletion->set_nickname("vd_5_del");
        parms.add_event(deletion);
        CHECK(parms.requires_extended_format());
    }
    SECTION("a gene choice whose v1 class names another segment") {
        auto j_choice = std::make_shared<Gene_choice>(V_gene);
        j_choice->add_realization("J1", "ACGT");
        j_choice->set_seq_type("J_gene_seq");
        j_choice->set_nickname("j_choice");
        parms.add_event(j_choice);
        CHECK(parms.requires_extended_format());
    }
    SECTION("two events sharing a v1 name") {
        for (const std::string nickname : {"d_5_del", "d_5_del_twin"}) {
            auto deletion = std::make_shared<Deletion>("D_gene_seq", Five_prime, std::make_pair(0, 2));
            deletion->set_nickname(nickname);
            parms.add_event(deletion);
        }
        CHECK(parms.requires_extended_format());
    }
}
