/**
 * @file test_inference_outputs.cpp
 * @brief The files an inference writes: a model the v1 format cannot hold, and its marginals, as
 * JSON documents named .json.
 *
 * No iteration is needed either: the initial and final files are written outside the EM loop, and
 * each iteration's through the same names. An iteration on a tandem model waits for T7a, which makes
 * Single_error_rate read the mismatches of every segment rather than of the legacy ids.
 */

#include <catch2/catch_test_macros.hpp>

#include "ScratchDir.h"

#include <igor/Core/Legacy/AlignmentData.h>
#include <igor/Inference/Legacy/GenModel.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

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
