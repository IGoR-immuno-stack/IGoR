/**
 * @file test_tandem_guards.cpp
 * @brief The counters refuse a model they cannot read, and an inference refuses before it runs
 * (capability plan T7a).
 *
 * The four counters read V, D and J by their legacy names: on a tandem model the best scenarios and
 * the errors leave both Ds' mismatches out, the Pgen counter assembles V+J, and D coverage
 * terminates the process. Until CT2a of docs/COUNTER_REWRITE_PLAN.md moves them off the enum, they
 * refuse any ordering but V-D-J and V-J. GenModel::infer_model() asks the error model and every
 * counter before it starts. The error models' half is in
 * tst/igor/Model/Legacy/test_tandem_guards.cpp.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "ScratchDir.h"

#include <igor/Inference/Legacy/Bestscenarioscounter.h>
#include <igor/Inference/Legacy/Coverageerrcounter.h>
#include <igor/Inference/Legacy/Errorscounter.h>
#include <igor/Inference/Legacy/GenModel.h>
#include <igor/Inference/Legacy/Pgencounter.h>
#include <igor/Model/Legacy/Hypermutationglobalerrorrate.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

using namespace igor::core::legacy;
using namespace igor::model::legacy;
using namespace igor::inference::legacy;
using Catch::Matchers::ContainsSubstring;

namespace {

Model_Parms read_fixture(const std::string &name)
{
    Model_Parms parms;
    parms.read_model_parms(std::string(IGOR_SOURCE_DIR) + "/tst/test_data/model_formats/" + name);
    return parms;
}

/// The counter accepts the VDJ and VJ fixtures, and refuses the tandem one by name.
void check_refuses_tandem_only(const Counter &counter, const std::string &name)
{
    CHECK_NOTHROW(counter.check_admissible(read_fixture("vdj_model.json").get_seq_type_registry()));
    CHECK_NOTHROW(counter.check_admissible(read_fixture("vj_model.json").get_seq_type_registry()));
    CHECK_THROWS_WITH(counter.check_admissible(read_fixture("tandem_d_model.json").get_seq_type_registry()),
                      ContainsSubstring(name + " reads a V-D-J or a V-J model only")
                              && ContainsSubstring("where it cannot place VD1_ins_seq"));
}

}  // namespace

TEST_CASE("The best-scenarios counter refuses a tandem model", "[tandem_guard][tandem_d][counter]")
{
    check_refuses_tandem_only(Best_scenarios_counter(1), "BestScenarioCounter");
}

TEST_CASE("The coverage counter refuses a tandem model", "[tandem_guard][tandem_d][counter]")
{
    // On V too, whose coverage alone would not crash: it reads V and J at the enum's ids all the same.
    check_refuses_tandem_only(Coverage_err_counter(V_gene), "CoverageErrCounter");
}

TEST_CASE("The errors counter refuses a tandem model", "[tandem_guard][tandem_d][counter]")
{
    check_refuses_tandem_only(Errors_counter(1), "ErrorsCounter");
}

TEST_CASE("The Pgen counter refuses a tandem model", "[tandem_guard][tandem_d][counter]")
{
    check_refuses_tandem_only(Pgen_counter(), "PgenCounter");
}

TEST_CASE("An inference on a tandem model refuses a reader of V, D and J before it runs",
          "[tandem_guard][tandem_d][genmodel]")
{
    Model_Parms parms = read_fixture("tandem_d_model.json");
    Model_marginals marginals(parms);
    marginals.uniform_initialize(parms);
    ScratchDir scratch;
    const std::string dir = scratch.file("");

    SECTION("a counter, before it writes its file")
    {
        std::map<std::size_t, std::shared_ptr<Counter>> counters;
        counters.emplace(0, std::make_shared<Best_scenarios_counter>(1, dir, true));
        GenModel model(parms, marginals, counters);
        CHECK_THROWS_WITH(model.infer_model({}, /*iterations=*/1, dir, /*fast_iter=*/true, /*likelihood_threshold=*/0.0,
                                            /*proba_threshold_factor=*/1e-4),
                          ContainsSubstring("BestScenarioCounter reads a V-D-J or a V-J model only"));
        CHECK_FALSE(std::filesystem::exists(dir + "best_scenarios_counts.csv"));
        CHECK_FALSE(std::filesystem::exists(dir + "initial_model.json"));
    }
    SECTION("a hypermutation error model")
    {
        parms.set_error_ratep(std::make_shared<Hypermutation_global_errorrate>(3, VDJ_genes, VDJ_genes, 0.1));
        GenModel model(parms, marginals);
        CHECK_THROWS_WITH(model.infer_model({}, /*iterations=*/1, dir, /*fast_iter=*/true, /*likelihood_threshold=*/0.0,
                                            /*proba_threshold_factor=*/1e-4),
                          ContainsSubstring("HypermutationGlobalErrorRate reads a V-D-J or a V-J model only"));
        CHECK_FALSE(std::filesystem::exists(dir + "initial_model.json"));
    }
}
