/*
 * test_SamplingHandlers.cpp
 *
 * Tests for CategoricalSamplingHandler and MarkovSamplingHandler.
 * Validates:
 *   (1) CDF precomputation produces monotone, [0,1]-bounded cumulative sums.
 *   (2) sample() output distribution matches the probability tensor
 *       (chi-squared / frequency check over N draws).
 *   (3) sample_sequence() drives the Markov chain correctly.
 *   (5) Scenario struct compiles and is index-addressable.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <igor/Generation/CategoricalSamplingHandler.h>
#include <igor/Generation/MarkovSamplingHandler.h>
#include <igor/Generation/SamplingHandler.h>
#include <igor/Generation/SamplingHandlerFactory.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Scenario.h>
#include <igor/Model/Topology.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy {}
namespace igor::generation::legacy {}
namespace igor::generation {}
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;
using namespace igor::generation::legacy;
using namespace igor::generation;

using namespace igor::model;
using namespace igor::math;
using namespace Catch;

// ─── Helpers ─────────────────────────────────────────────────────────────────

// Frequency-based test: draw N samples, check each bin's frequency is within
// tolerance of its expected probability.
static void check_distribution(
    auto& handler,
    std::mt19937_64& rng,
    const std::vector<double>& expected_probs,
    const std::vector<std::size_t>& parent_indices,
    std::size_t N = 20'000,
    double tol = 0.03)
{
    std::vector<std::size_t> counts(expected_probs.size(), 0);
    for (std::size_t i = 0; i < N; ++i) {
        auto idx = handler.sample(rng, parent_indices);
        REQUIRE(idx < expected_probs.size());
        ++counts[idx];
    }
    for (std::size_t k = 0; k < expected_probs.size(); ++k) {
        double freq = static_cast<double>(counts[k]) / N;
        INFO("bin " << k << ": expected " << expected_probs[k] << " got " << freq);
        REQUIRE(std::abs(freq - expected_probs[k]) < tol);
    }
}

// ─── Scenario ─────────────────────────────────────────────────────────────────

TEST_CASE("Scenario construction", "[Model][Sampling][Scenario]") {
    SampledScenario sc(3);
    REQUIRE(sc.events.size() == 3);

    sc.events[0] = { 0, {2} };
    sc.events[1] = { 1, {0, 3, 1} };
    sc.events[2] = { 2, {5} };

    REQUIRE(sc.index_of(0) == 2);
    REQUIRE(sc.index_of(2) == 5);
    REQUIRE(sc.events[1].indices.size() == 3);
}

// ─── CategoricalSamplingHandler ───────────────────────────────────────────────

TEST_CASE("CategoricalSamplingHandler: construction and CDF", "[Model][Sampling][Categorical]") {
    // 1D categorical — 4 realizations, no parents
    Tensor<double> weights({4});
    weights.data()[0] = 0.1;
    weights.data()[1] = 0.4;
    weights.data()[2] = 0.3;
    weights.data()[3] = 0.2;

    CategoricalSamplingHandler<double> sh("test_cat_1d", -1, weights);

    REQUIRE(sh.realizationCount() == 4);
    REQUIRE(sh.parentSliceCount() == 1);
    REQUIRE(sh.name() == "test_cat_1d");
    REQUIRE(sh.uid() == -1);  // unassigned

    sh.precomputeCDF();

    // CDF values for the single slice
    // Expected: [0.1, 0.5, 0.8, 1.0]
}

TEST_CASE("CategoricalSamplingHandler: sampling 1D", "[Model][Sampling][Categorical]") {
    std::vector<double> probs = {0.1, 0.4, 0.3, 0.2};
    Tensor<double> weights({4});
    for (std::size_t i = 0; i < 4; ++i) weights.data()[i] = probs[i];

    CategoricalSamplingHandler<double> sh("cat_1d", -1, weights);
    sh.precomputeCDF();

    std::mt19937_64 rng(42);
    check_distribution(sh, rng, probs, {}, 30'000, 0.025);
}

TEST_CASE("CategoricalSamplingHandler: sampling 2D (one parent)", "[Model][Sampling][Categorical]") {
    // shape: [2 parent states, 3 realizations]
    // Parent 0 → probs [0.6, 0.3, 0.1]
    // Parent 1 → probs [0.2, 0.5, 0.3]
    Tensor<double> weights({2, 3});
    double* d = weights.data();
    d[0] = 0.6; d[1] = 0.3; d[2] = 0.1;  // parent=0
    d[3] = 0.2; d[4] = 0.5; d[5] = 0.3;  // parent=1

    CategoricalSamplingHandler<double> sh("cat_2d", -1, weights);
    sh.precomputeCDF();

    std::mt19937_64 rng(123);

    SECTION("parent = 0") {
        check_distribution(sh, rng, {0.6, 0.3, 0.1}, {0}, 30'000, 0.02);
    }
    SECTION("parent = 1") {
        check_distribution(sh, rng, {0.2, 0.5, 0.3}, {1}, 30'000, 0.02);
    }
}

TEST_CASE("CategoricalSamplingHandler: uid protocol", "[Model][Sampling][Categorical]") {
    Tensor<double> weights({2});
    CategoricalSamplingHandler<double> sh("uid_cat", -1, weights);
    REQUIRE(sh.uid() == -1);
    sh.setUid(7);
    REQUIRE(sh.uid() == 7);
}

// ─── MarkovSamplingHandler ────────────────────────────────────────────────────

TEST_CASE("MarkovSamplingHandler: construction", "[Model][Sampling][Markov]") {
    // 4-state Markov, no parents
    Tensor<double> weights({4, 4});
    MarkovSamplingHandler<double> sh("test_markov", -1, weights);

    REQUIRE(sh.stateCount() == 4);
    REQUIRE(sh.parentSliceCount() == 1);
    REQUIRE(sh.name() == "test_markov");
}

TEST_CASE("MarkovSamplingHandler: first nucleotide (uniform marginal)", "[Model][Sampling][Markov]") {
    Tensor<double> weights({4, 4});
    double* d = weights.data();
    for (std::size_t i = 0; i < 16; ++i) d[i] = 0.25;

    MarkovSamplingHandler<double> sh("uniform_markov", -1, weights);
    sh.precomputeCDF();

    std::mt19937_64 rng(42);
    std::vector<double> expected(4, 0.25);
    check_distribution(sh, rng, expected, {}, 40'000, 0.02);
}

TEST_CASE("MarkovSamplingHandler: transition row sampling", "[Model][Sampling][Markov]") {
    // 3-state Markov, no parents
    // Row 0: [0.7, 0.2, 0.1]
    // Row 1: [0.1, 0.8, 0.1]
    // Row 2: [0.3, 0.3, 0.4]
    Tensor<double> weights({3, 3});
    double* d = weights.data();
    d[0] = 0.7; d[1] = 0.2; d[2] = 0.1;  // row 0
    d[3] = 0.1; d[4] = 0.8; d[5] = 0.1;  // row 1
    d[6] = 0.3; d[7] = 0.3; d[8] = 0.4;  // row 2

    MarkovSamplingHandler<double> sh("biased_markov", -1, weights);
    sh.precomputeCDF();

    std::mt19937_64 rng(99);

    SECTION("from state 0") {
        check_distribution(sh, rng, {0.7, 0.2, 0.1}, {0}, 30'000, 0.025);
    }
    SECTION("from state 1") {
        check_distribution(sh, rng, {0.1, 0.8, 0.1}, {1}, 30'000, 0.025);
    }
    SECTION("from state 2") {
        check_distribution(sh, rng, {0.3, 0.3, 0.4}, {2}, 30'000, 0.025);
    }
}

TEST_CASE("MarkovSamplingHandler: conditioned on a parent", "[Model][Sampling][Markov]") {
    // Shape [2 parent values, 4 from-states, 4 to-states]: the shape a Dinucl_markov event gets
    // as soon as it has a parent. Every (parent, from) row is a one-hot on one to-state, so a
    // draw is deterministic and names the row that was actually read.
    //
    // The table is chosen so that the two rows a swapped index would confuse never share a
    // value: precomputeCDF() stores row ps * 4 + from, and the reader used to look up row
    // from * 2 + ps, which maps (0,1)->(0,2), (0,2)->(1,0), (0,3)->(1,2), (1,0)->(0,1),
    // (1,1)->(0,3), (1,2)->(1,1). Each of those pairs differs below.
    const std::size_t expected_to[2][4] = { { 0, 3, 2, 1 },     // parent 0
                                            { 1, 0, 3, 2 } };   // parent 1

    Tensor<double> weights({2, 4, 4});
    std::fill(weights.begin(), weights.end(), 0.0);
    for (std::size_t ps = 0; ps < 2; ++ps)
        for (std::size_t from = 0; from < 4; ++from)
            weights.data()[ps * 16 + from * 4 + expected_to[ps][from]] = 1.0;

    MarkovSamplingHandler<double> sh("cond_markov", -1, weights);
    REQUIRE(sh.stateCount() == 4);
    REQUIRE(sh.parentSliceCount() == 2);
    sh.precomputeCDF();

    std::mt19937_64 rng(7);

    SECTION("sample() reads the row of its (parent, from) pair") {
        for (std::size_t ps = 0; ps < 2; ++ps) {
            for (std::size_t from = 0; from < 4; ++from) {
                INFO("parent " << ps << ", from " << from);
                // The from-state comes first, the parents after it: the order sampleSequence() builds.
                REQUIRE(sh.sample(rng, { from, ps }) == expected_to[ps][from]);
            }
        }
    }

    SECTION("sampleSequence() carries the parent through the chain") {
        // Under parent 1 the table alternates 0 -> 1 -> 0 -> 1; under parent 0 it stays on 0.
        REQUIRE(sh.sampleSequence(rng, /*first_state=*/0, /*n_steps=*/3, { 1 })
                == std::vector<std::size_t>{ 0, 1, 0, 1 });
        REQUIRE(sh.sampleSequence(rng, /*first_state=*/0, /*n_steps=*/3, { 0 })
                == std::vector<std::size_t>{ 0, 0, 0, 0 });
    }

    SECTION("a missing parent index is refused") {
        REQUIRE_THROWS_AS(sh.sample(rng, { 0 }), std::invalid_argument);
    }
}

TEST_CASE("MarkovSamplingHandler: sample_sequence", "[Model][Sampling][Markov]") {
    // 2-state Markov: deterministic alternating
    // Row 0 → always go to 1: [0.0, 1.0]
    // Row 1 → always go to 0: [1.0, 0.0]
    Tensor<double> weights({2, 2});
    double* d = weights.data();
    d[0] = 0.0; d[1] = 1.0;
    d[2] = 1.0; d[3] = 0.0;

    MarkovSamplingHandler<double> sh("det_markov", -1, weights);
    sh.precomputeCDF();

    std::mt19937_64 rng(0);
    auto chain = sh.sampleSequence(rng, /*first_state=*/0, /*n_steps=*/5);

    // Expected: 0 → 1 → 0 → 1 → 0 → 1
    REQUIRE(chain.size() == 6);
    REQUIRE(chain[0] == 0);
    REQUIRE(chain[1] == 1);
    REQUIRE(chain[2] == 0);
    REQUIRE(chain[3] == 1);
    REQUIRE(chain[4] == 0);
    REQUIRE(chain[5] == 1);
}

TEST_CASE("MarkovSamplingHandler: uid protocol", "[Model][Sampling][Markov]") {
    Tensor<double> weights({2, 2});
    MarkovSamplingHandler<double> sh("uid_mrk", -1, weights);
    REQUIRE(sh.uid() == -1);
    sh.setUid(3);
    REQUIRE(sh.uid() == 3);
}

// ─── SamplingHandlerFactory ───────────────────────────────────────────────────

TEST_CASE("SamplingHandlerFactory: build from RecombinationModel", "[Model][Sampling][Factory]") {
    // 1. Create a minimal topology: Root1, Root2 -> Child -> GrandChild
    auto topology = std::make_unique<Topology>();

    std::unordered_map<std::string, Event_realization> root1_realizations;
    root1_realizations.emplace("R1_1", Event_realization("R1_1", 0, "R1_1", Int_Str(), 0));
    root1_realizations.emplace("R1_2", Event_realization("R1_2", 1, "R1_2", Int_Str(), 1));
    root1_realizations.emplace("R1_3", Event_realization("R1_3", 2, "R1_3", Int_Str(), 2));
    auto root1 = std::make_shared<Gene_choice>(Undefined_gene, root1_realizations);
    root1->set_nickname("Root1");
    root1->set_priority(10);

    std::unordered_map<std::string, Event_realization> root2_realizations;
    root2_realizations.emplace("R2_1", Event_realization("R2_1", 0, "R2_1", Int_Str(), 0));
    root2_realizations.emplace("R2_2", Event_realization("R2_2", 1, "R2_2", Int_Str(), 1));
    root2_realizations.emplace("R2_3", Event_realization("R2_3", 2, "R2_3", Int_Str(), 2));
    root2_realizations.emplace("R2_4", Event_realization("R2_4", 3, "R2_4", Int_Str(), 3));
    auto root2 = std::make_shared<Gene_choice>(Undefined_gene, root2_realizations);
    root2->set_nickname("Root2");
    root2->set_priority(5);

    std::unordered_map<std::string, Event_realization> child_realizations;
    child_realizations.emplace("D1", Event_realization("D1", 1, "", Int_Str(), 0));
    child_realizations.emplace("D2", Event_realization("D2", 2, "", Int_Str(), 1));
    // A deletion and a dinucleotide chain are identified by their target segment now, not by a
    // gene class, and neither can be default-constructed: the segment goes in the constructor.
    auto child = std::make_shared<Deletion>(V_gene_seq, Three_prime, child_realizations);
    child->set_nickname("Child");

    auto grandchild = std::make_shared<Dinucl_markov>(VD_ins_seq);
    grandchild->set_nickname("GrandChild");

    igor::core::legacy::index_type root1_id = topology->addEvent(root1);
    igor::core::legacy::index_type root2_id = topology->addEvent(root2);
    igor::core::legacy::index_type child_id = topology->addEvent(child);
    igor::core::legacy::index_type grandchild_id = topology->addEvent(grandchild);

    topology->addEdge(root1_id, child_id);
    topology->addEdge(root2_id, child_id);
    topology->addEdge(child_id, grandchild_id);
    topology->addEdge(root1_id, grandchild_id);

    // 2. Test topological ordering based on priority.
    //    Core's Event_comparator orders by DESCENDING priority, so Root1 (priority 10)
    //    must be visited before Root2 (priority 5).
    std::vector<igor::core::legacy::index_type> order = topology->topologicalOrder();
    REQUIRE(order.size() == 4);
    REQUIRE(order[0] == root1_id);
    REQUIRE(order[1] == root2_id);
    REQUIRE(order[2] == child_id);
    REQUIRE(order[3] == grandchild_id);

    // 3. Build RecombinationModel from the topology and then build handlers
    RecombinationModel<double> model(std::move(topology));
    auto handlers = sampling_handler_factory::build<double>(model);

    REQUIRE(handlers.size() == 4);

    // 4. Verify handlers

    // Root1 handler: Categorical, shape {3}
    auto* h_root1 = static_cast<CategoricalSamplingHandler<double>*>(handlers[root1_id].get());
    REQUIRE(h_root1 != nullptr);
    REQUIRE(h_root1->name() == std::string("Root1"));
    REQUIRE(h_root1->uid() == root1_id);
    REQUIRE(h_root1->realizationCount() == 3);
    REQUIRE(h_root1->parentSliceCount() == 1);

    // Root2 handler: Categorical, shape {4}
    auto* h_root2 = static_cast<CategoricalSamplingHandler<double>*>(handlers[root2_id].get());
    REQUIRE(h_root2 != nullptr);
    REQUIRE(h_root2->name() == std::string("Root2"));
    REQUIRE(h_root2->uid() == root2_id);
    REQUIRE(h_root2->realizationCount() == 4);
    REQUIRE(h_root2->parentSliceCount() == 1);

    // Child handler: Categorical, shape {2, 3, 4} (2 self, 3 from Root1, 4 from Root2)
    auto* h_child = static_cast<CategoricalSamplingHandler<double>*>(handlers[child_id].get());
    REQUIRE(h_child != nullptr);
    REQUIRE(h_child->uid() == child_id);
    REQUIRE(h_child->realizationCount() == 2);
    REQUIRE(h_child->parentSliceCount() == 12); // 3 (Root1) * 4 (Root2) = 12

    // Grandchild handler: Markov, shape {4, 4, 2, 3}
    // inherent_shape for Dinucl_markov is {4, 4}, parents are Child (inherent {2}) and Root1 (inherent {3})
    auto ptr = handlers[grandchild_id].get();
    REQUIRE(ptr != nullptr);

    auto* h_grandchild = static_cast<MarkovSamplingHandler<double>*>(ptr);

    REQUIRE(h_grandchild != nullptr);
    REQUIRE(h_grandchild->uid() == grandchild_id);
    REQUIRE(h_grandchild->stateCount() == 4);
    REQUIRE(h_grandchild->parentSliceCount() == 6); // 2 (Child) * 3 (Root1) = 6
    REQUIRE(h_grandchild->weights().size() == 4 * 4 * 6); // 4x4 matrix for each of the 6 parent combinations
}

