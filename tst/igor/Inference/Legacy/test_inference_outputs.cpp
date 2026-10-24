/**
 * @file test_inference_outputs.cpp
 * @brief The files an inference writes: a model the v1 format cannot hold, and its marginals, as
 * JSON documents named .json.
 *
 * No iteration is needed either: the initial and final files are written outside the EM loop, and
 * each iteration's through the same names. The last case runs one, on the same fixture.
 */

#include <catch2/catch_test_macros.hpp>

#include "ScratchDir.h"

#include <igor/Core/Legacy/AlignmentData.h>
#include <igor/Inference/Legacy/GenModel.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
using namespace igor::inference::legacy;

namespace {

Alignment_data aligned(const std::string &gene, int offset, std::size_t length)
{
    return Alignment_data(gene, offset, static_cast<std::size_t>(offset), offset + length - 1, length - 1, {}, {}, {},
                          100.0);
}

}  // namespace

TEST_CASE("An inference writes a JSON model's outputs, marginals included, as .json", "[model_format][json][genmodel]")
{
    Model_Parms parms;
    parms.read_model_parms(std::string(IGOR_SOURCE_DIR) + "/tst/test_data/model_formats/tandem_d_model.json");
    REQUIRE(parms.requires_extended_format());
    Model_marginals marginals(parms);
    marginals.uniform_initialize(parms);

    ScratchDir scratch;
    const std::string dir = scratch.file("");
    GenModel model(parms, marginals);
    // V_test, D1_test, D2_test and J_test of the fixture, end to end. D1 and D2 share the D alignments.
    const std::string read = std::string("ATCGATCG") + "GGGG" + "CCCC" + "GCTAGCTA";
    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;
    alignments[V_gene] = { aligned("V_test", 0, 8) };
    alignments[D_gene] = { aligned("D1_test", 8, 4), aligned("D2_test", 12, 4) };
    alignments[J_gene] = { aligned("J_test", 16, 8) };
    REQUIRE_FALSE(model.infer_model({ { 0, read, alignments } }, /*iterations=*/0, dir, /*fast_iter=*/true,
                                    /*likelihood_threshold=*/0.0, /*proba_threshold_factor=*/1e-4));

    for (const char *name : { "initial_model", "initial_marginals", "final_parms", "final_marginals" }) {
        INFO(name);
        CHECK(std::filesystem::exists(dir + name + ".json"));
        CHECK_FALSE(std::filesystem::exists(dir + name + ".txt"));
    }
    // Named .json and written as JSON: the name is GenModel's choice, the format the writer's.
    CHECK(model_file_has_parameters(dir + "initial_marginals.json"));
    CHECK(model_file_has_parameters(dir + "final_marginals.json"));

    SECTION("the final marginals read back with the final model") {
        Model_Parms final_parms;
        REQUIRE_NOTHROW(final_parms.read_model_parms(dir + "final_parms.json"));
        Model_marginals final_marginals(final_parms);
        REQUIRE_NOTHROW(final_marginals.read_marginals(dir + "final_marginals.json", final_parms));
        const Model_marginals &inferred = model.get_marginals();
        REQUIRE(final_marginals.get_length() == inferred.get_length());
        for (std::size_t i = 0; i != inferred.get_length(); ++i) {
            const double a = static_cast<double>(inferred.marginal_array_smart_p[i]);
            CHECK(std::fabs(static_cast<double>(final_marginals.marginal_array_smart_p[i]) - a) <= 4e-16 * std::fabs(a));
        }
    }
}

TEST_CASE("An iteration on the tandem fixture places each D at its own template's alignment", "[tandem_d][genmodel]")
{
    // Each D of the fixture holds one template, D1_test and D2_test, and both read the D
    // alignments: each passes over the other's. The read is the four templates end to end, so one
    // scenario explains it.
    Model_Parms parms;
    parms.read_model_parms(std::string(IGOR_SOURCE_DIR) + "/tst/test_data/model_formats/tandem_d_model.json");
    Model_marginals marginals(parms);
    marginals.uniform_initialize(parms);

    ScratchDir scratch;
    const std::string dir = scratch.file("");
    GenModel model(parms, marginals);
    const std::string read = std::string("ATCGATCG") + "GGGG" + "CCCC" + "GCTAGCTA";
    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;
    alignments[V_gene] = { aligned("V_test", 0, 8) };
    alignments[D_gene] = { aligned("D1_test", 8, 4), aligned("D2_test", 12, 4) };
    alignments[J_gene] = { aligned("J_test", 16, 8) };
    REQUIRE_FALSE(model.infer_model({ { 0, read, alignments } }, /*iterations=*/1, dir, /*fast_iter=*/true,
                                    /*likelihood_threshold=*/0.0, /*proba_threshold_factor=*/1e-4));

    // iteration_n;seq_processed;seq_index;nt_sequence;n_V_aligns;n_J_aligns;seq_likelihood;
    // seq_mean_n_errors;seq_n_scenarios;seq_best_scenario;time
    std::ifstream log(dir + "inference_logs.txt");
    std::string header, line;
    REQUIRE(std::getline(log, header));
    REQUIRE(std::getline(log, line));
    std::vector<std::string> fields;
    for (std::size_t start = 0, end; start <= line.size(); start = end + 1) {
        end = std::min(line.find(';', start), line.size());
        fields.push_back(line.substr(start, end - start));
    }
    REQUIRE(fields.size() == 11);
    CHECK(fields[8] == "1");
    CHECK(std::stod(fields[6]) > 0.0);
}
