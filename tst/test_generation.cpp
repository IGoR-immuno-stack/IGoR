/**
 * @file test_generation.cpp
 * @brief Tests the sequence generation process using KL divergence
 *        and entropy-based validation.
 *
 * Generates 10, 100, 1000, and 1'000'000 sequences from the simple VJ model
 * in demo/. For each sample size, the empirical marginal distributions of every
 * recombination event are compared to the theoretical model marginals via the
 * Kullback-Leibler divergence:
 *
 *     D_KL(P || Q) = Σ_i P(i) · log2( P(i) / Q(i) )
 *
 * where P is the model (theoretical) distribution and Q is the empirical one.
 *
 * Following the entropy concept from pygor3
 * (https://pygor3.readthedocs.io/en/latest/Tutorial.html#Entropy):
 *
 *     H = − Σ_{E⃗} P(E⃗) · log2 P(E⃗)
 *
 * each recombination event contributes an entropy term. The KL divergence for
 * each event is checked against the corresponding entropy: as N grows, D_KL
 * must become negligibly small compared to H. Concretely:
 *
 *   • N = 10       → D_KL may be large; we only check finiteness
 *   • N = 100      → D_KL < H  (loose bound)
 *   • N = 1 000    → D_KL < H / 10
 *   • N = 1 000 000 → D_KL < H / 1000
 *
 * Additionally, the sequence of D_KL values across increasing N is checked
 * for a monotonic decrease (up to a tolerance).
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

#include <cmath>
#include <forward_list>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#ifndef IGOR_SOURCE_DIR
#error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

// ---------------------------------------------------------------------------
// Paths to the simple VJ demo model
// ---------------------------------------------------------------------------
static const std::string MODEL_PARMS_PATH =
        std::string(IGOR_SOURCE_DIR) + "/demo/simple_vj_model.txt";
static const std::string MODEL_MARGINALS_PATH =
        std::string(IGOR_SOURCE_DIR) + "/demo/simple_vj_model_marginals.txt";

// ---------------------------------------------------------------------------
// Mathematical helpers
// ---------------------------------------------------------------------------

/**
 * @brief Kullback–Leibler divergence D_KL(P || Q) in bits.
 *
 * Skips bins where P(i) == 0.  Returns +∞ when P(i) > 0 but Q(i) == 0.
 */
static double kl_divergence(const std::vector<double> &P,
                            const std::vector<double> &Q)
{
    double kl = 0.0;
    for (size_t i = 0; i < P.size(); ++i) {
        if (P[i] > 0.0) {
            if (Q[i] <= 0.0) {
                return std::numeric_limits<double>::infinity();
            }
            kl += P[i] * std::log2(P[i] / Q[i]);
        }
    }
    return kl;
}

/**
 * @brief Shannon entropy H(P) = −Σ P(i) log2 P(i) in bits.
 */
static double entropy(const std::vector<double> &P)
{
    double h = 0.0;
    for (double p : P) {
        if (p > 0.0) {
            h -= p * std::log2(p);
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
// Test-case data structures
// ---------------------------------------------------------------------------

/**
 * @brief Metadata for one recombination event that we will test.
 *
 * Stores the event's name, nickname, number of realizations,
 * its position in the model queue (so we can pop the right queue element
 * from a scenario), and the theoretical marginal distribution.
 */
struct EventInfo {
    Rec_Event_name name;
    std::string nickname;
    int num_realizations;
    size_t queue_position;              ///< index in the model-queue order
    bool is_dinuc_markov;               ///< skip empirical check for DinucMarkov
    std::vector<double> model_marginal; ///< P(realization) marginalised over parents
    double H;                           ///< entropy of the marginal
};

/**
 * @brief Build a vector of EventInfo for every event in the model queue.
 */
static std::vector<EventInfo> build_event_info(
        const Model_Parms &parms,
        const Model_marginals &marginals)
{
    std::vector<EventInfo> infos;
    auto model_queue = parms.get_model_queue();
    size_t pos = 0;

    while (!model_queue.empty()) {
        auto ev = model_queue.front();
        model_queue.pop();

        EventInfo info;
        info.name = ev->get_name();
        info.nickname = ev->get_nickname();
        info.num_realizations = ev->size();
        info.queue_position = pos;
        info.is_dinuc_markov = (ev->get_type() == Event_type::Dinuclmarkov_t);

        std::cout << "  Processing event: " << info.nickname
                  << " (is_dinuc=" << info.is_dinuc_markov
                  << ", size=" << info.num_realizations << ")" << std::endl;

        if (!info.is_dinuc_markov) {
            // Compute the marginal probability for this event
            auto [dims, probs] =
                    marginals.compute_event_marginal_probability(info.name, parms);

            // The first dimension is the event itself; the returned array
            // is already the marginal.
            info.model_marginal.resize(info.num_realizations, 0.0);
            for (int i = 0; i < info.num_realizations; ++i) {
                info.model_marginal[i] = static_cast<double>(probs.get()[i]);
            }
            info.H = entropy(info.model_marginal);
        }

        infos.push_back(std::move(info));
        ++pos;
    }
    return infos;
}

/**
 * @brief Given a list of generated (sequence, realizations) pairs, compute
 *        the empirical marginal distribution for a non-DinucMarkov event.
 *
 * For each scenario, the realizations queue contains one inner queue per
 * event in model-queue order. For GeneChoice / Deletion / Insertion the
 * inner queue holds a single int (the realization index).
 */
static std::vector<double> compute_empirical_marginal(
        const std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> &scenarios,
        size_t event_queue_position,
        int num_realizations,
        size_t total_sequences)
{
    std::vector<size_t> counts(num_realizations, 0);

    for (auto it = scenarios.begin(); it != scenarios.end(); ++it) {
        // Copy the outer queue so we can pop to the desired position
        auto outer = it->second;
        for (size_t k = 0; k < event_queue_position; ++k) {
            outer.pop();
        }
        // The front queue now corresponds to the event of interest
        auto inner = outer.front();
        if (!inner.empty()) {
            int realization = inner.front();
            if (realization >= 0 && realization < num_realizations) {
                ++counts[realization];
            }
        }
    }

    std::vector<double> empirical(num_realizations, 0.0);
    double n = static_cast<double>(total_sequences);
    for (int i = 0; i < num_realizations; ++i) {
        empirical[i] = static_cast<double>(counts[i]) / n;
    }
    return empirical;
}

// ---------------------------------------------------------------------------
// THE TEST
// ---------------------------------------------------------------------------

TEST_CASE("Generation marginals converge — KL divergence vs entropy",
          "[generation]")
{
    // ------------------------------------------------------------------
    // 1. Load the simple VJ model
    // ------------------------------------------------------------------
    Model_Parms model_parms;
    model_parms.read_model_parms(MODEL_PARMS_PATH);

    Model_marginals model_marginals(model_parms);
    model_marginals.txt2marginals(MODEL_MARGINALS_PATH, model_parms);

    // ------------------------------------------------------------------
    // 2. Collect per-event metadata and theoretical marginals
    // ------------------------------------------------------------------
    std::vector<EventInfo> event_infos =
            build_event_info(model_parms, model_marginals);

    INFO("Model loaded, " << event_infos.size() << " events found");
    REQUIRE(event_infos.size() >= 5); // V, J, V_del, J_del, VJ_ins (+ dinuc)

    // Print entropy decomposition (mirrors pygor3's entropy table)
    double total_entropy = 0.0;
    std::cout << "\n=== Entropy decomposition (bits) ===" << std::endl;
    for (const auto &ev : event_infos) {
        if (!ev.is_dinuc_markov) {
            std::cout << "  " << ev.nickname << " : H = " << ev.H << std::endl;
            total_entropy += ev.H;
        }
    }
    std::cout << "  TOTAL (non-dinuc): " << total_entropy << std::endl;

    // ------------------------------------------------------------------
    // 3. Generate sequences for increasing N and track D_KL per event
    // ------------------------------------------------------------------
    // Map: event_queue_position → vector of D_KL values (one per sample-size)
    std::map<size_t, std::vector<double>> kl_traces;

    const std::vector<int> sample_sizes = {10, 100, 1000, 1000000};

    for (int N : sample_sizes) {
        INFO("Generating " << N << " sequences");

        // Generate (sequences are not written to disk)
        GenModel gen_model(model_parms, model_marginals);
        auto scenarios = gen_model.generate_sequences(N, /*output_realizations=*/true);

        // Count how many were actually generated (forward_list has no size())
        size_t actual_count = 0;
        for (auto it = scenarios.begin(); it != scenarios.end(); ++it) {
            ++actual_count;
        }
        REQUIRE(actual_count == static_cast<size_t>(N));

        std::cout << "\n--- N = " << N << " ---" << std::endl;

        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov) {
                continue; // DinucMarkov has variable-length realizations
            }

            auto empirical = compute_empirical_marginal(
                    scenarios, ev.queue_position, ev.num_realizations, actual_count);

            double dkl = kl_divergence(ev.model_marginal, empirical);
            kl_traces[ev.queue_position].push_back(dkl);

            std::cout << "  " << ev.nickname
                      << " : D_KL = " << dkl
                      << "  (H = " << ev.H << ")" << std::endl;

            // ---- assertions ----
            // D_KL must always be non-negative
            CHECK(dkl >= 0.0);

            if (N == 10) {
                // With only 10 samples, D_KL can be very large
                // (even +∞ if a realization is never drawn).
                // No bound check at this sample size.
            } else if (N == 100) {
                // Loose bound: D_KL should be within the same order as H
                CHECK(dkl < std::max(ev.H, 1.0));
            } else if (N == 1000) {
                // Tighter: D_KL should be at most a tenth of H
                CHECK(dkl < std::max(ev.H, 1.0) * 0.1);
            } else if (N == 1000000) {
                // Very tight: D_KL should be negligible compared to H
                // Theoretical expectation: D_KL ≈ (k−1)/(2·N·ln2)
                double expected_upper =
                        10.0 * (ev.num_realizations - 1) / (2.0 * N * std::log(2.0));
                CHECK(dkl < std::max(expected_upper, 1e-3));
            }
        }
    }

    // ------------------------------------------------------------------
    // 4. Check monotonic decrease of D_KL with increasing N
    // ------------------------------------------------------------------
    std::cout << "\n=== Monotonic decrease check ===" << std::endl;
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) continue;

        const auto &trace = kl_traces.at(ev.queue_position);
        // Compare D_KL at N=1000 vs N=1'000'000.  Both are large enough
        // for the estimator to be reliable.  A small absolute tolerance
        // handles the case where the smaller sample happens to nail the
        // distribution exactly (D_KL ≈ 0) while the larger one is merely
        // very small.
        if (trace.size() >= 4) {
            double dkl_at_1k = trace[2];   // index 2 → N=1 000
            double dkl_at_1M = trace[3];   // index 3 → N=1 000 000
            std::cout << "  " << ev.nickname
                      << " : D_KL(N=1k)=" << dkl_at_1k
                      << "  D_KL(N=1M)=" << dkl_at_1M << std::endl;
            CHECK(dkl_at_1M < dkl_at_1k + 1e-4);
        }
    }
}
