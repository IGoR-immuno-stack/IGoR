/**
 * @file test_tandem_trb_model.cpp
 * @brief The V-D1-D2-J model of milestone 1 against the TRB model it is built from (capability
 * plan T5a).
 *
 * scripts/tandem_d/make_tandem_d_model.py writes tst/test_data/tandem_d/trb_tandem_d_model.json
 * from models/human/tcr_beta/models: TRB's D gene and its deletions twice, as D1 and D2, and its
 * two junctions as three. These cases read both models into IGoR. Each event of the tandem model
 * must be the TRB event it copies, on its own segment: the same node, the same templates, its
 * parents renamed within its D, and the same parameters once IGoR has normalized both.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

namespace {

const std::string kTrbDir = std::string(IGOR_SOURCE_DIR) + "/models/human/tcr_beta/models/";
const std::string kTandem = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/tandem_d/trb_tandem_d_model.json";

/// Each event of the tandem model, the TRB event it copies, and its segment: the script's COPIES.
struct Copy
{
    std::string nickname;
    std::string source;
    std::string segment;
};
const std::vector<Copy> kCopies = {
    { "v_choice", "v_choice", "V_gene_seq" },       { "v_3_del", "v_3_del", "V_gene_seq" },
    { "vd1_ins", "vd_ins", "VD1_ins_seq" },         { "vd1_dinucl", "vd_dinucl", "VD1_ins_seq" },
    { "d1_choice", "d_gene", "D1_gene_seq" },       { "d1_5_del", "d_5_del", "D1_gene_seq" },
    { "d1_3_del", "d_3_del", "D1_gene_seq" },       { "d1d2_ins", "vd_ins", "D1D2_ins_seq" },
    { "d1d2_dinucl", "vd_dinucl", "D1D2_ins_seq" }, { "d2_choice", "d_gene", "D2_gene_seq" },
    { "d2_5_del", "d_5_del", "D2_gene_seq" },       { "d2_3_del", "d_3_del", "D2_gene_seq" },
    { "d2j_ins", "dj_ins", "D2J_ins_seq" },         { "d2j_dinucl", "dj_dinucl", "D2J_ins_seq" },
    { "j_choice", "j_choice", "J_gene_seq" },       { "j_5_del", "j_5_del", "J_gene_seq" },
};

std::map<std::string, std::string> source_of()
{
    std::map<std::string, std::string> out;
    for (const Copy &copy : kCopies)
        out[copy.nickname] = copy.source;
    return out;
}

/// A model and its marginals, read as the CLI reads them.
struct Read
{
    Model_Parms parms;
    Model_marginals marginals;

    Read(const std::string &parms_file, const std::string &marginals_file)
    {
        parms.read_model_parms(parms_file);
        marginals = Model_marginals(parms);
        marginals.read_marginals(marginals_file, parms);
    }
};

const Read &trb()
{
    static const Read read(kTrbDir + "model_parms.txt", kTrbDir + "model_marginals.txt");
    return read;
}

/// The tandem model reads its parameters from its own document, as -set_custom_model does.
const Read &tandem()
{
    static const Read read(kTandem, kTandem);
    return read;
}

nlohmann::json event_node(const Model_Parms &parms, const std::string &nickname)
{
    return parms.get_event_pointer(nickname, true)->to_json();
}

/// Each value of a parameters node keyed by its coordinates' labels, its axes' events renamed by
/// \a rename and the parents sorted by event: what the node says, whatever order its axes come in.
std::map<std::string, double> by_labels(const nlohmann::json &node, const std::map<std::string, std::string> &rename)
{
    const nlohmann::json &axes = node.at("axes");
    const auto renamed = [&](const std::string &event) {
        const auto found = rename.find(event);
        return found == rename.end() ? event : found->second;
    };
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
                    parents[renamed(event)] = label;
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

} // namespace

TEST_CASE("The tandem TRB model reads, with its parameters", "[tandem_d][model_format][json]")
{
    if (!std::filesystem::exists(kTrbDir + "model_parms.txt"))
        SKIP("TRB model not found at " << kTrbDir << " (submodule may not be initialized)");

    const Model_Parms &parms = tandem().parms;
    CHECK(parms.get_seq_type_registry().get_ordered_types()
          == std::vector<std::string>{ "V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq", "D2_gene_seq",
                                       "D2J_ins_seq", "J_gene_seq" });
    std::set<std::string> nicknames;
    for (const auto &event : parms.get_event_list())
        nicknames.insert(event->get_nickname());
    std::set<std::string> copies;
    for (const Copy &copy : kCopies)
        copies.insert(copy.nickname);
    CHECK(nicknames == copies);
    CHECK(parms.requires_extended_format());
}

TEST_CASE("The tandem TRB model is the document IGoR writes for it", "[tandem_d][model_format][json]")
{
    if (!std::filesystem::exists(kTrbDir + "model_parms.txt"))
        SKIP("TRB model not found at " << kTrbDir << " (submodule may not be initialized)");

    // The parameters aside: IGoR writes them normalized, to 17 digits, where the script copies
    // TRB's 6. The cases below compare them.
    std::ifstream in(kTandem);
    nlohmann::json written = parse_model_document(in);
    REQUIRE(written.contains("parameters"));
    written.erase("parameters");
    CHECK(model_parms_to_json(tandem().parms) == written);
}

TEST_CASE("Each event of the tandem TRB model is the TRB event it copies", "[tandem_d][model_format][json]")
{
    if (!std::filesystem::exists(kTrbDir + "model_parms.txt"))
        SKIP("TRB model not found at " << kTrbDir << " (submodule may not be initialized)");

    for (const Copy &copy : kCopies) {
        CAPTURE(copy.nickname);
        nlohmann::json expected = event_node(trb().parms, copy.source);
        expected["targets"][0]["segment"] = copy.segment;
        CHECK(event_node(tandem().parms, copy.nickname) == expected);
    }

    SECTION("with TRB's templates, D1's and D2's once")
    {
        CHECK(model_parms_to_json(tandem().parms).at("library") == model_parms_to_json(trb().parms).at("library"));
    }
    SECTION("and TRB's error rate")
    {
        CHECK(model_parms_to_json(tandem().parms).at("error_rate")
              == model_parms_to_json(trb().parms).at("error_rate"));
    }
}

TEST_CASE("Each D of the tandem TRB model is conditioned as TRB's D is, within itself",
          "[tandem_d][model_format][json]")
{
    if (!std::filesystem::exists(kTrbDir + "model_parms.txt"))
        SKIP("TRB model not found at " << kTrbDir << " (submodule may not be initialized)");

    const std::map<std::string, std::set<std::string>> expected = {
        { "v_choice", {} },
        { "v_3_del", { "v_choice" } },
        { "j_choice", { "v_choice" } },
        { "j_5_del", { "j_choice" } },
        { "d1_choice", { "v_choice", "j_choice" } },
        { "d1_5_del", { "d1_choice" } },
        { "d1_3_del", { "d1_choice", "d1_5_del" } },
        { "d2_choice", { "v_choice", "j_choice" } },
        { "d2_5_del", { "d2_choice" } },
        { "d2_3_del", { "d2_choice", "d2_5_del" } },
        { "vd1_ins", {} },
        { "vd1_dinucl", {} },
        { "d1d2_ins", {} },
        { "d1d2_dinucl", {} },
        { "d2j_ins", {} },
        { "d2j_dinucl", {} },
    };
    const auto parents_of = [](const Model_Parms &parms, const std::string &nickname) {
        std::set<std::string> out;
        for (const auto &parent : parms.get_parents(parms.get_event_pointer(nickname, true)))
            out.insert(parent->get_nickname());
        return out;
    };
    const std::map<std::string, std::string> source = source_of();
    for (const Copy &copy : kCopies) {
        CAPTURE(copy.nickname);
        const std::set<std::string> parents = parents_of(tandem().parms, copy.nickname);
        CHECK(parents == expected.at(copy.nickname));

        // TRB's, once each parent is read back as the event it copies.
        std::set<std::string> as_sources;
        for (const std::string &parent : parents)
            as_sources.insert(source.at(parent));
        CHECK(as_sources == parents_of(trb().parms, copy.source));
    }
}

TEST_CASE("The parameters of the tandem TRB model are TRB's, event by event", "[tandem_d][model_format][json]")
{
    if (!std::filesystem::exists(kTrbDir + "model_parms.txt"))
        SKIP("TRB model not found at " << kTrbDir << " (submodule may not be initialized)");

    // Both files give TRB's values to 6 digits, each read as the nearest double and normalized
    // the same way: they match exactly.
    const nlohmann::json trb_nodes = model_marginals_to_json(trb().parms, trb().marginals).at("nodes");
    const nlohmann::json tandem_nodes = model_marginals_to_json(tandem().parms, tandem().marginals).at("nodes");
    const std::map<std::string, std::string> source = source_of();
    for (const Copy &copy : kCopies) {
        CAPTURE(copy.nickname);
        const auto copied = by_labels(tandem_nodes.at(copy.nickname), source);
        CHECK(copied.size() == by_labels(trb_nodes.at(copy.source), {}).size());
        CHECK(copied == by_labels(trb_nodes.at(copy.source), {}));
    }
}
