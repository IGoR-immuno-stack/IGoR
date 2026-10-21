/**
 * @file test_best_scenarios_counter.cpp
 * @brief What the best-scenarios output names each event by (output.event_names, plan T2).
 */

#include <catch2/catch_test_macros.hpp>

#include "ScratchDir.h"

#include <igor/Inference/Legacy/Bestscenarioscounter.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>

#include <fstream>
#include <string>

using namespace igor::model::legacy;
using namespace igor::inference::legacy;

namespace {

std::string header_of(const ScratchDir &dir, Model_Parms &parms, Event_label labels)
{
    Model_marginals marginals(parms);
    Best_scenarios_counter counter(1, dir.file(""), true);
    counter.set_event_labels(labels);
    counter.initialize_counter(parms, marginals);

    std::ifstream in(dir.file("best_scenarios_counts.csv"));
    std::string header;
    std::getline(in, header);
    return header;
}

}  // namespace

TEST_CASE("The best-scenarios header names each event by nickname, or by its legacy name",
          "[counter][event_name]")
{
    ScratchDir dir;
    Model_Parms parms;
    parms.read_model_parms(std::string(IGOR_SOURCE_DIR) + "/tst/test_data/format_v2/test_legacy_model_parms.txt");

    SECTION("by nickname, the default") {
        CHECK(header_of(dir, parms, Event_label::Nickname)
              == "seq_index;scenario_rank;scenario_proba_cond_seq;"
                 "v_choice;j_choice;v_3_del;vj_ins;vj_dinucl;Mismatches");
    }
    SECTION("by the legacy name, which the regression goldens carry") {
        CHECK(header_of(dir, parms, Event_label::Legacy)
              == "seq_index;scenario_rank;scenario_proba_cond_seq;"
                 "GeneChoice_V_gene_Undefined_side_prio7_size1;GeneChoice_J_gene_Undefined_side_prio7_size1;"
                 "Deletion_V_gene_Three_prime_prio5_size2;Insertion_VJ_gene_Undefined_side_prio4_size2;"
                 "DinucMarkov_VJ_gene_Undefined_side_prio1_size16;Mismatches");
    }
}

TEST_CASE("The CLI's spellings of the event names", "[counter][event_name]")
{
    CHECK(str2EventLabel("nickname") == Event_label::Nickname);
    CHECK(str2EventLabel("legacy") == Event_label::Legacy);
    CHECK_THROWS_AS(str2EventLabel("Nickname"), std::invalid_argument);
}
