/**
 * @file test_model_parameters.cpp
 * @brief The parameters section of a JSON model document (ModelJson.h): the legacy marginals.
 *
 * A node per event, its axes labelled, its values nested one level per axis. What the writer lays
 * out, that values come back at double precision, and that the reader matches them by label: a
 * fragment fits a model whatever order that model gives its events, which the text file cannot do.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "ScratchDir.h"

#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
using Catch::Matchers::ContainsSubstring;

namespace {

/// VDJ with conditioning: d_3_del has two parents, d_gene and d_5_del.
const std::string kVdj = std::string(IGOR_SOURCE_DIR) + "/scripts/tests/data/input/TRB_model_parms.txt";
/// VJ with inferred, non-uniform marginals.
const std::string kVjParms = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/inference/fixed_VJ_TRB_model_parms.txt";
const std::string kVjMarginals = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/inference/fixed_VJ_TRB_marginals.txt";

Model_Parms read_parms(const std::string &path)
{
    Model_Parms parms;
    parms.read_model_parms(path);
    return parms;
}

/// Distinct values in every block, normalized: a value read at the wrong place shows.
Model_marginals distinct_marginals(const Model_Parms &parms)
{
    Model_marginals marginals(parms);
    for (std::size_t i = 0; i != marginals.get_length(); ++i)
        marginals.marginal_array_smart_p[i] = 1 + (i * 7919) % 101;
    const auto queue = parms.get_model_queue();
    marginals.normalize(marginals.get_inverse_offset_map(parms, queue), marginals.get_index_map(parms, queue), queue);
    return marginals;
}

nlohmann::json fragment(const Model_Parms &parms, const Model_marginals &marginals)
{
    nlohmann::json doc;
    doc["schema_version"] = kModelJsonSchemaVersion;
    doc["parameters"] = model_marginals_to_json(parms, marginals);
    return doc;
}

Model_marginals read_into(const Model_Parms &parms, const nlohmann::json &doc)
{
    Model_marginals marginals(parms);
    model_marginals_from_json(doc, parms, marginals);
    return marginals;
}

bool same_array(const Model_marginals &a, const Model_marginals &b)
{
    if (a.get_length() != b.get_length())
        return false;
    for (std::size_t i = 0; i != a.get_length(); ++i)
        if (a.marginal_array_smart_p[i] != b.marginal_array_smart_p[i])
            return false;
    return true;
}

/// A node's values flattened in row-major order.
void flatten(const nlohmann::json &values, std::vector<double> &out)
{
    if (values.is_array()) {
        for (const auto &v : values)
            flatten(v, out);
    } else {
        out.push_back(values.get<double>());
    }
}

/// Each value of a node keyed by its coordinates' labels, the parent axes sorted by event: what a
/// node says, whatever order it lists its axes and labels in.
std::map<std::string, double> by_labels(const nlohmann::json &node)
{
    const nlohmann::json &axes = node.at("axes");
    std::map<std::string, double> out;
    std::vector<std::size_t> position(axes.size(), 0);
    const auto walk = [&](const auto &self, const nlohmann::json &values, std::size_t depth) -> void {
        if (depth == axes.size()) {
            std::map<std::string, std::string> parents;
            std::string own;
            const std::string nickname = axes.back().at("event");
            for (std::size_t k = 0; k != axes.size(); ++k) {
                const std::string event = axes[k].at("event");
                const std::string label = axes[k].at("labels")[position[k]].dump();
                if (event == nickname)
                    own += label + ";";
                else
                    parents[event] = label;
            }
            std::string key;
            for (const auto &[event, label] : parents)
                key += event + "=" + label + ",";
            out[key + own] = values.get<double>();
            return;
        }
        for (std::size_t i = 0; i != values.size(); ++i) {
            position[depth] = i;
            self(self, values[i], depth + 1);
        }
    };
    walk(walk, node.at("values"), 0);
    return out;
}

/// The node with its axes put in \a order, and the labels of axis \a reversed (an index into the
/// new order) listed backwards; its values follow, so it says the same thing.
nlohmann::json rearranged(const nlohmann::json &node, const std::vector<std::size_t> &order, std::size_t reversed)
{
    const nlohmann::json &axes = node.at("axes");
    nlohmann::json out;
    out["axes"] = nlohmann::json::array();
    std::vector<std::size_t> sizes;
    for (std::size_t k = 0; k != order.size(); ++k) {
        nlohmann::json axis = axes[order[k]];
        if (k == reversed)
            std::reverse(axis["labels"].begin(), axis["labels"].end());
        sizes.push_back(axis.at("labels").size());
        out["axes"].push_back(std::move(axis));
    }
    std::vector<std::size_t> new_position(order.size()), old_position(order.size());
    const auto build = [&](const auto &self, std::size_t depth) -> nlohmann::json {
        if (depth == order.size()) {
            for (std::size_t k = 0; k != order.size(); ++k)
                old_position[order[k]] = (k == reversed) ? sizes[k] - 1 - new_position[k] : new_position[k];
            const nlohmann::json *value = &node.at("values");
            for (const std::size_t p : old_position)
                value = &(*value)[p];
            return *value;
        }
        nlohmann::json level = nlohmann::json::array();
        for (std::size_t i = 0; i != sizes[depth]; ++i) {
            new_position[depth] = i;
            level.push_back(self(self, depth + 1));
        }
        return level;
    };
    out["values"] = build(build, 0);
    return out;
}

std::string text_of(const std::string &path)
{
    std::ifstream in(path);
    return { std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>() };
}

}  // namespace

TEST_CASE("The parameters section lays each event out as the legacy array does", "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    const Model_marginals marginals = distinct_marginals(parms);
    const nlohmann::json section = model_marginals_to_json(parms, marginals);
    const auto queue = parms.get_model_queue();
    const auto index_map = marginals.get_index_map(parms, queue);

    REQUIRE(section.at("nodes").size() == parms.get_event_list().size());

    SECTION("the parents in the text file's order, then the own axis, labels in index order") {
        const nlohmann::json &axes = section.at("nodes").at("d_3_del").at("axes");
        REQUIRE(axes.size() == 3);
        CHECK(axes[0].at("event") == "d_gene");
        CHECK(axes[1].at("event") == "d_5_del");
        CHECK(axes[2].at("event") == "d_3_del");
        for (const auto &axis : axes) {
            const auto event = parms.get_event_pointer(axis.at("event").get<std::string>(), true);
            std::vector<std::pair<int, nlohmann::json>> by_index;
            for (const auto &entry : event->get_realizations_map())
                by_index.emplace_back(entry.second.index, realization_label(entry.second));
            std::sort(by_index.begin(), by_index.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
            nlohmann::json expected = nlohmann::json::array();
            for (const auto &[index, label] : by_index)
                expected.push_back(label);
            CHECK(axis.at("labels") == expected);
        }
        CHECK(axes[2].at("labels").front().is_number_integer());
    }

    SECTION("the values read straight through: row-major, they are the event's block") {
        for (const auto &event : parms.get_event_list()) {
            INFO(event->get_nickname());
            std::vector<double> flat;
            flatten(section.at("nodes").at(event->get_nickname()).at("values"), flat);
            const std::size_t base = index_map.at(event->get_name());
            for (std::size_t i = 0; i != flat.size(); ++i)
                REQUIRE(flat[i] == static_cast<double>(marginals.marginal_array_smart_p[base + i]));
        }
    }

    SECTION("a Markov chain has two own axes, from then to") {
        const nlohmann::json &chain = section.at("nodes").at("vd_dinucl");
        REQUIRE(chain.at("axes")
                == nlohmann::json::array({ { { "event", "vd_dinucl" }, { "labels", { "A", "C", "G", "T" } } },
                                           { { "event", "vd_dinucl" }, { "labels", { "A", "C", "G", "T" } } } }));
        const std::size_t base = index_map.at("vd_dinucl");
        CHECK(chain.at("values")[1][2].get<double>() == static_cast<double>(marginals.marginal_array_smart_p[base + 6]));
    }

    SECTION("the text file agrees on the parents' order") {
        ScratchDir scratch;
        const std::string text_file = scratch.file("marginals.txt");
        Model_marginals(marginals).write2txt(text_file, parms);
        CHECK_THAT(text_of(text_file), ContainsSubstring("@d_3_del\n$Dim[3,21,21]\n#[d_gene,0],[d_5_del,0]\n"));
    }
}

TEST_CASE("Marginals written as JSON read back at double precision", "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVjParms);
    Model_marginals original(parms);
    original.txt2marginals(kVjMarginals, parms);

    ScratchDir scratch;
    const std::string written = scratch.file("marginals.json");
    original.write_marginals_json(written, parms);
    REQUIRE(text_of(written).front() == '{');

    Model_marginals read_back(parms);
    read_back.read_marginals(written, parms);
    for (std::size_t i = 0; i != original.get_length(); ++i) {
        const long double a = original.marginal_array_smart_p[i];
        const long double b = read_back.marginal_array_smart_p[i];
        INFO("index " << i);
        REQUIRE(std::fabs(static_cast<double>(a - b)) <= 4e-16 * std::fabs(static_cast<double>(a)) + 1e-300);
    }

    SECTION("written again, the file is the same, text for text") {
        const std::string again = scratch.file("again.json");
        read_back.write_marginals_json(again, parms);
        REQUIRE(text_of(again) == text_of(written));
    }
    SECTION("and the text file still reads as text") {
        Model_marginals text(parms);
        REQUIRE_NOTHROW(text.read_marginals(kVjMarginals, parms));
        REQUIRE(same_array(text, original));
    }
}

TEST_CASE("A complete model document carries its parameters, read as a fragment is",
          "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    const Model_marginals marginals = distinct_marginals(parms);
    nlohmann::json complete = model_parms_to_json(parms);
    complete["parameters"] = model_marginals_to_json(parms, marginals);

    ScratchDir scratch;
    const std::string complete_file = scratch.file("model.json", complete.dump(2));
    const std::string fragment_file = scratch.file("marginals.json", fragment(parms, marginals).dump(2));

    Model_Parms from_complete;
    REQUIRE_NOTHROW(from_complete.read_model_parms(complete_file));
    CHECK(model_file_has_parameters(complete_file));
    CHECK(model_file_has_parameters(fragment_file));
    CHECK_FALSE(model_file_has_parameters(kVdj));
    CHECK_FALSE(model_file_has_parameters(scratch.file("structure.json", model_parms_to_json(parms).dump())));

    Model_marginals a(from_complete), b(from_complete);
    a.read_marginals(complete_file, from_complete);
    b.read_marginals(fragment_file, from_complete);
    REQUIRE(same_array(a, b));
}

TEST_CASE("Parameters are matched by label, whatever order the axes, the labels or the model's events come in",
          "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    const nlohmann::json doc = fragment(parms, distinct_marginals(parms));
    const Model_marginals as_written = read_into(parms, doc);

    SECTION("parent axes swapped, and labels listed backwards") {
        nlohmann::json edited = doc;
        auto &nodes = edited["parameters"]["nodes"];
        nodes["d_3_del"] = rearranged(nodes["d_3_del"], { 1, 0, 2 }, 2);
        nodes["j_choice"] = rearranged(nodes["j_choice"], { 0 }, 0);
        nodes["vd_dinucl"] = rearranged(nodes["vd_dinucl"], { 0, 1 }, 1);
        REQUIRE(edited != doc);
        REQUIRE(same_array(read_into(parms, edited), as_written));
    }

    SECTION("a model listing its events in another order, as one read from JSON without the hint") {
        nlohmann::json structure = model_parms_to_json(parms);
        structure.erase("engine_hints");
        Model_Parms reordered;
        model_parms_from_json(structure, reordered);
        const Model_marginals other = read_into(reordered, doc);

        // The positional text file could not have been moved across: its layout is not the same.
        ScratchDir scratch;
        Model_marginals(as_written).write2txt(scratch.file("a.txt"), parms);
        Model_marginals(other).write2txt(scratch.file("b.txt"), reordered);
        REQUIRE(text_of(scratch.file("a.txt")) != text_of(scratch.file("b.txt")));

        const nlohmann::json a = model_marginals_to_json(parms, as_written).at("nodes");
        const nlohmann::json b = model_marginals_to_json(reordered, other).at("nodes");
        for (const auto &[nickname, node] : a.items()) {
            INFO(nickname);
            REQUIRE(by_labels(node) == by_labels(b.at(nickname)));
        }
    }
}

TEST_CASE("Values that do not sum to one are normalized, as the text reader does", "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    nlohmann::json doc = fragment(parms, distinct_marginals(parms));
    auto &values = doc["parameters"]["nodes"]["vd_ins"]["values"];
    for (auto &value : values)
        value = 0;
    values[0] = 1;
    values[1] = 3;
    const Model_marginals read = read_into(parms, doc);
    const std::size_t base = read.get_index_map(parms).at("vd_ins");
    CHECK(read.marginal_array_smart_p[base] == 0.25L);
    CHECK(read.marginal_array_smart_p[base + 1] == 0.75L);
}

TEST_CASE("Parameters that do not fit the model are refused", "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    nlohmann::json doc = fragment(parms, distinct_marginals(parms));
    auto &nodes = doc["parameters"]["nodes"];
    const auto refused = [&](const std::string &part) {
        Model_marginals marginals(parms);
        REQUIRE_THROWS_WITH(model_marginals_from_json(doc, parms, marginals), ContainsSubstring(part));
    };

    SECTION("a missing node, or one for no event") {
        nlohmann::json kept = nodes["vd_ins"];
        nodes.erase("vd_ins");
        refused("no node for event \"vd_ins\"");
        nodes["vd_ins"] = kept;
        nodes["vd_insertion"] = kept;
        refused("\"vd_insertion\", which is no event");
    }
    SECTION("an unknown key in a node, or in an axis") {
        nodes["vd_ins"]["encoding"] = "decimal";
        refused("unknown key \"encoding\"");
        nodes["vd_ins"].erase("encoding");
        nodes["vd_ins"]["axes"][0]["index"] = 0;
        refused("unknown key \"index\"");
    }
    SECTION("the conditioning edited without the parameters") {
        nlohmann::json structure = model_parms_to_json(parms);
        structure["structure"]["conditioning"]["vd_ins"] = { "v_choice" };
        Model_Parms edited;
        model_parms_from_json(structure, edited);
        Model_marginals marginals(edited);
        REQUIRE_THROWS_WITH(model_marginals_from_json(doc, edited, marginals),
                            ContainsSubstring("\"vd_ins\"") && ContainsSubstring("1 parents"));
    }
    SECTION("an axis naming an event that is not a parent, or a parent twice") {
        nodes["d_3_del"]["axes"][1]["event"] = "v_choice";
        refused("\"v_choice\", which is not a parent");
        nodes["d_3_del"]["axes"][1]["event"] = "d_gene";
        refused("two axes name the parent \"d_gene\"");
    }
    SECTION("the own axis anywhere but last, or a parent where it comes") {
        const nlohmann::json kept = nodes["d_3_del"];
        nodes["d_3_del"] = rearranged(kept, { 2, 1, 0 }, 3);
        refused("own axis comes");
        // d_5_del has d_3_del's labels, so only the name can tell the axes apart.
        nodes["d_3_del"] = kept;
        nodes["d_3_del"]["axes"][2]["event"] = "d_5_del";
        refused("\"d_5_del\", where the event's own axis comes");
    }
    SECTION("a label the event does not have, one listed twice, or one missing") {
        auto &labels = nodes["v_3_del"]["axes"][1]["labels"];
        const nlohmann::json first = labels[0];
        labels[0] = 99;
        refused("no realization labelled 99");
        labels[0] = labels[1];
        refused("twice");
        labels.erase(labels.begin());
        refused("must list the 21 labels");
        labels.insert(labels.begin(), first);
        labels[0] = "-4";
        refused("no realization labelled \"-4\"");
    }
    SECTION("values of the wrong shape, negative, or not numbers") {
        nodes["j_choice"]["values"].erase(0);
        refused("not a list of 15");
        nodes = fragment(parms, distinct_marginals(parms))["parameters"]["nodes"];
        nodes["j_choice"]["values"][0] = -0.5;
        refused("not a finite, non-negative number");
        nodes["j_choice"]["values"][0] = "0.5";
        refused("not a finite, non-negative number");
    }
    SECTION("a document without parameters, or with another version") {
        nlohmann::json structure = model_parms_to_json(parms);
        Model_marginals marginals(parms);
        REQUIRE_THROWS_WITH(model_marginals_from_json(structure, parms, marginals),
                            ContainsSubstring("missing key \"parameters\""));
        doc["schema_version"] = "0.2";
        refused("schema version \"0.2\"");
    }
    SECTION("marginals sized for another model") {
        Model_marginals vj(read_parms(kVjParms));
        REQUIRE_THROWS_AS(model_marginals_from_json(doc, parms, vj), std::logic_error);
    }
}

TEST_CASE("A value too large for a double is refused", "[model_format][json][parameters]")
{
    const Model_Parms parms = read_parms(kVdj);
    nlohmann::json doc = fragment(parms, distinct_marginals(parms));
    ScratchDir scratch;

    SECTION("in the text: the parser refuses the number") {
        std::string text = doc.dump();
        const std::string values = "\"values\":[";
        const auto first = text.find(values, text.find("\"vd_ins\":{")) + values.size();
        text.replace(first, text.find_first_of(",]", first) - first, "1e999");
        Model_marginals marginals(parms);
        REQUIRE_THROWS_WITH(marginals.read_marginals(scratch.file("overflow.json", text), parms),
                            ContainsSubstring("not valid JSON") && ContainsSubstring("overflow"));
    }
    SECTION("in a document built in code: the reader refuses the infinity") {
        doc["parameters"]["nodes"]["vd_ins"]["values"][0] = std::numeric_limits<double>::infinity();
        Model_marginals marginals(parms);
        REQUIRE_THROWS_WITH(model_marginals_from_json(doc, parms, marginals),
                            ContainsSubstring("vd_ins") && ContainsSubstring("not a finite, non-negative number"));
    }
}
