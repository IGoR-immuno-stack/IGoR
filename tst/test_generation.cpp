/**
 * @file test_generation.cpp
 * @brief Tests the sequence generation process using KL divergence
 *        and entropy-based validation.
 *
 * Uses the modern Topology + SamplingEngine architecture to generate
 * 10, 100, 1000, and 1'000'000 scenarios from the human TCR alpha
 * (VJ) and human TCR beta (VDJ) models. For each sample size, the
 * empirical marginal distributions of every recombination event are
 * compared to the theoretical model marginals via the Kullback-Leibler
 * divergence:
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
 *   • N = 1 000    → D_KL < H / 10,   uncovered mass < 10%
 *   • N = 1 000 000 → D_KL < H / 1000, uncovered mass < 0.01%
 *
 * Bins where P(i) > 0 but Q(i) == 0 (unsampled realizations) are skipped
 * in the D_KL sum; their total probability mass is tracked separately as
 * "uncovered mass" and is expected to shrink with increasing N.
 *
 * Additionally, the sequence of D_KL values across increasing N is checked
 * for a monotonic decrease (up to a tolerance).
 *
 * For the TCR beta model, the computed entropy decomposition is cross-
 * validated against the reference values from pygor3's tutorial:
 * https://pygor3.readthedocs.io/en/latest/Tutorial.html#Entropy
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/LegacyBridge.h>
#include <igor/Model/SamplingEngine.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Topology.h>
#include <igor/Model/Scenario.h>

#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

#include <igor/Model/Topology.h>

#include "entropy_test_helpers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef IGOR_SOURCE_DIR
#  error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

// Model base directory
static const std::string MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/models";

/**
 * @brief Build EventInfo metadata by iterating the Topology in
 *        topological order.
 *
 * Uses the Topology's UIDs as position keys (stored in
 * EventInfo::queue_position) so that events map directly to
 * SampledScenario::events[uid].
 *
 * Theoretical marginals are still computed from the legacy
 * Model_marginals, which is the only source that provides the
 * marginalization over parent dimensions.
 */
static std::vector<EventInfo>
build_event_info_from_topology(const igor::model::Topology &topology,
                               const Model_Parms &parms,
                               const Model_marginals &marginals)
{
    std::vector<EventInfo> infos;
    auto index_map = marginals.get_index_map(parms);

    for (auto uid : topology.topologicalOrder()) {
        auto ev = topology.event(uid);

        EventInfo info;
        info.name = ev->get_name();
        info.nickname = ev->get_nickname();
        info.num_realizations = ev->size();
        info.queue_position = static_cast<size_t>(uid); // UID as key
        info.is_dinuc_markov = (ev->get_type() == Event_type::Dinuclmarkov_t);
        info.gene_class = ev->get_class();
        info.dinuc_T.fill(0.0);
        info.dinuc_entropy_rate = 0.0;

        std::cout << "  Processing event: " << info.nickname
                  << " (uid=" << uid
                  << ", is_dinuc=" << info.is_dinuc_markov
                  << ", size=" << info.num_realizations << ")"
                  << std::endl;

        if (info.is_dinuc_markov) {
            int base_idx = index_map.at(info.name);
            for (int k = 0; k < 16; ++k) {
                info.dinuc_T[k] = static_cast<double>(
                    marginals.marginal_array_smart_p[base_idx + k]);
            }
            info.dinuc_entropy_rate = markov_entropy_rate(info.dinuc_T);

            info.model_marginal.resize(16, 0.0);
            for (int k = 0; k < 16; ++k) {
                info.model_marginal[k] = info.dinuc_T[k];
            }
            info.H = info.dinuc_entropy_rate;
        } else {
            auto [dims, probs] =
                marginals.compute_event_marginal_probability(info.name, parms);

            info.model_marginal.resize(info.num_realizations, 0.0);
            for (int i = 0; i < info.num_realizations; ++i) {
                info.model_marginal[i] = static_cast<double>(probs.get()[i]);
            }
            info.H = entropy(info.model_marginal);
        }

        infos.push_back(std::move(info));
    }
    return infos;
}

/**
 * @brief Compute empirical marginal distributions from SampledScenario output.
 *
 * For each non-DinucMarkov event, counts the realization index
 * (SampledEvent::indices[0]) across all scenarios and normalises
 * to a probability vector.
 *
 * Events are keyed by their Topology UID (stored in
 * EventInfo::queue_position).
 */
static std::map<size_t, std::vector<double>>
compute_all_empirical_marginals(
    const std::vector<igor::model::SampledScenario> &scenarios,
    const std::vector<EventInfo> &event_infos,
    size_t total_sequences)
{
    // Prepare per-event count arrays (non-DinucMarkov only)
    std::map<size_t, std::vector<size_t>> counts;
    for (const auto &ev : event_infos) {
        if (!ev.is_dinuc_markov) {
            counts[ev.queue_position].assign(ev.num_realizations, 0);
        }
    }

    // Single pass over all generated scenarios
    for (const auto &scenario : scenarios) {
        for (auto &[uid, cvec] : counts) {
            if (uid < scenario.events.size()
                && !scenario.events[uid].indices.empty()) {
                auto realization =
                    static_cast<int>(scenario.events[uid].indices[0]);
                int num_real = static_cast<int>(cvec.size());
                if (realization < 0 || realization >= num_real) {
                    throw std::out_of_range(
                        "Realization index " + std::to_string(realization)
                        + " out of bounds [0, "
                        + std::to_string(num_real) + ")");
                }
                ++cvec[realization];
            }
        }
    }

    // Convert counts to probabilities
    std::map<size_t, std::vector<double>> result;
    double n = static_cast<double>(total_sequences);
    for (auto &[uid, cvec] : counts) {
        std::vector<double> empirical(cvec.size(), 0.0);
        for (size_t i = 0; i < cvec.size(); ++i) {
            empirical[i] = static_cast<double>(cvec[i]) / n;
        }
        result[uid] = std::move(empirical);
    }
    return result;
}

/**
 * @brief Compute empirical marginal distributions from SampledScenario output.
 *
 * For each non-DinucMarkov event, counts the sampled realization index
 * (SampledEvent::indices[0]) and normalises to a probability distribution.
 *
 * @param scenarios  The generated scenarios (one per sample).
 * @param topology   The model topology (maps UID → Rec_Event).
 * @param event_infos Reference event info (from build_event_info).
 * @param nick_to_uid  Mapping from event nickname to topology UID.
 * @return Map from topology UID to empirical probability vector.
 */
static std::map<igor::index_type, std::vector<double>>
compute_all_empirical_marginals(
    const std::vector<igor::model::SampledScenario>& scenarios,
    const std::vector<EventInfo>& event_infos,
    const std::map<std::string, igor::index_type>& nick_to_uid)
{
    // Prepare per-event count arrays (non-DinucMarkov only)
    std::map<igor::index_type, std::vector<std::size_t>> counts;
    for (const auto& ev : event_infos) {
        if (!ev.is_dinuc_markov) {
            auto it = nick_to_uid.find(ev.nickname);
            if (it != nick_to_uid.end()) {
                counts[it->second].assign(
                    static_cast<std::size_t>(ev.num_realizations), 0);
            }
        }
    }

    // Single pass over all scenarios
    for (const auto& scenario : scenarios) {
        for (auto& [uid, cvec] : counts) {
            std::size_t val = scenario.index_of(uid);
            if (val >= cvec.size()) {
                throw std::out_of_range(
                    "Realization index " + std::to_string(val) +
                    " out of bounds [0, " + std::to_string(cvec.size()) + ")");
            }
            ++cvec[val];
        }
    }

    // Convert counts to probabilities
    std::map<igor::index_type, std::vector<double>> result;
    double n = static_cast<double>(scenarios.size());
    for (auto& [uid, cvec] : counts) {
        std::vector<double> empirical(cvec.size(), 0.0);
        for (std::size_t i = 0; i < cvec.size(); ++i) {
            empirical[i] = static_cast<double>(cvec[i]) / n;
        }
        result[uid] = std::move(empirical);
    }
    return result;
}

// ---------------------------------------------------------------------------
// THE TEST
// ---------------------------------------------------------------------------

TEST_CASE("Generation marginals converge - KL divergence vs entropy", "[generation]")
{
    using namespace igor::model;

    // ------------------------------------------------------------------
    // 1. Select model to test (TCR alpha or TCR beta)
    // ------------------------------------------------------------------
    std::string model_parms_path;
    std::string model_marginals_path;
    std::string model_label;
    int min_events = 0;
    // Reference entropy values from pygor3 for cross-validation.
    // Keys are event nicknames; values in bits.
    std::map<std::string, double> pygor3_reference_entropy;

    SECTION("human TCR alpha (VJ)")
    {
        model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        model_label = "human/tcr_alpha";
        min_events = 5; // V, J, V_del, J_del, VJ_ins (+ dinuc)
    }

    SECTION("human TCR beta (VDJ)")
    {
        model_parms_path = MODELS_DIR + "/human/tcr_beta/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_beta/models/model_marginals.txt";
        model_label = "human/tcr_beta";
        min_events = 9; // V, J, D, 4 dels, 2 ins (+ 2 dinuc)
        // Reference entropies from pygor3 tutorial:
        // https://pygor3.readthedocs.io/en/latest/Tutorial.html#Entropy
        // Insertion rows are the combined H(ℓ) + E[ℓ]·h values.
        pygor3_reference_entropy = {
            { "v_choice", 5.252905 }, { "d_gene", 1.141779 },  { "j_choice", 3.609102 },
            { "vd_ins", 14.894931 },  { "dj_ins", 14.981991 }, { "v_3_del", 3.147511 },
            { "d_3_del", 2.778230 },  { "d_5_del", 3.634137 }, { "j_5_del", 3.356340 },
        };
    }

    // ------------------------------------------------------------------
    // 2. Load legacy model, build Topology + SamplingEngine
    // ------------------------------------------------------------------
    INFO("Testing model: " << model_label);
    Model_Parms model_parms;
    model_parms.read_model_parms(model_parms_path);

    Model_marginals model_marginals(model_parms);
    model_marginals.txt2marginals(model_marginals_path, model_parms);

    // Build the modern Topology graph from the legacy model
    auto topology = import_from_legacy(model_parms);
    REQUIRE(topology != nullptr);
    REQUIRE(topology->size() > 0);

    // Build the RecombinationModel and import marginals, then wrap in SamplingEngine
    igor::model::RecombinationModel<double> model_obj(
        std::make_unique<igor::model::Topology>(*topology));
    import_from_legacy(model_obj, model_marginals);
    auto model_ptr = std::make_shared<const igor::model::RecombinationModel<double>>(
        std::move(model_obj));
    SamplingEngine<double> engine(model_ptr);

    // ------------------------------------------------------------------
    // 2b. Collect per-event metadata (using Topology UIDs)
    // ------------------------------------------------------------------
    std::vector<EventInfo> event_infos =
        build_event_info_from_topology(*topology, model_parms, model_marginals);

    INFO("Model loaded, " << event_infos.size() << " events found");
    REQUIRE(event_infos.size() >= static_cast<size_t>(min_events));

    // ------------------------------------------------------------------
    // Compute combined insertion+dinucleotide entropy.
    // Following pygor3's get_df_entropy_decomposition(), insertion events
    // like "vj_ins" are paired with their DinucMarkov partner "vj_dinucl"
    // (sharing the same Gene_class). The combined entropy is:
    //     H_total = H(ℓ) + H_ss + h · [E[ℓ] − (1 − P(0))]
    // where H(ℓ) is the insertion-length entropy, H_ss is the entropy of
    // the Markov chain's stationary distribution, E[ℓ] the expected length,
    // h the Markov chain entropy rate, and P(0) the probability of zero
    // insertions.
    // ------------------------------------------------------------------
    struct InsDinucPair
    {
        const EventInfo *ins_event = nullptr;
        const EventInfo *dinuc_event = nullptr;
        double combined_H = 0.0;
    };
    std::map<Gene_class, InsDinucPair> ins_dinuc_pairs;

    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) {
            ins_dinuc_pairs[ev.gene_class].dinuc_event = &ev;
        } else if (ev.nickname.find("_ins") != std::string::npos) {
            ins_dinuc_pairs[ev.gene_class].ins_event = &ev;
        }
    }

    for (auto &[gc, pair] : ins_dinuc_pairs) {
        if (pair.ins_event && pair.dinuc_event) {
            // Get insertion length values from the event's realizations
            auto ev_ptr = topology->event(
                static_cast<igor::index_type>(pair.ins_event->queue_position));
            auto realizations = ev_ptr->get_realizations_map();
            std::vector<int> ins_lengths(pair.ins_event->num_realizations, 0);
            for (const auto &[key, real] : realizations) {
                ins_lengths[real.index] = real.value_int;
            }

            pair.combined_H =
                    insertion_dinuc_entropy(pair.ins_event->model_marginal, ins_lengths, pair.dinuc_event->dinuc_T);
        }
    }

    // Print entropy decomposition (mirrors pygor3's entropy table)
    double total_entropy = 0.0;
    std::cout << "\n=== Entropy decomposition (bits) ===" << std::endl;
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) {
            // Already accounted for via the ins+dinuc pair
            continue;
        }

        // Check if this is an insertion event that has a DinucMarkov partner
        auto it = ins_dinuc_pairs.find(ev.gene_class);
        if (it != ins_dinuc_pairs.end() && it->second.ins_event == &ev && it->second.dinuc_event != nullptr) {
            double combined = it->second.combined_H;
            std::cout << "  " << ev.nickname << " + " << it->second.dinuc_event->nickname << " : H = " << combined
                      << "  (H_len=" << ev.H << ", H_dinuc=" << (combined - ev.H)
                      << ", h=" << it->second.dinuc_event->dinuc_entropy_rate << ")" << std::endl;
            total_entropy += combined;
        } else {
            std::cout << "  " << ev.nickname << " : H = " << ev.H << std::endl;
            total_entropy += ev.H;
        }
    }
    // Also print standalone DinucMarkov entropy rates for reference
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) {
            std::cout << "  (" << ev.nickname << " entropy rate h = " << ev.dinuc_entropy_rate << " bits/nt)"
                      << std::endl;
        }
    }
    std::cout << "  TOTAL: " << total_entropy << std::endl;

    // ------------------------------------------------------------------
    // Cross-validate against pygor3 reference values (if available).
    // The reference comes from mdl_hb.get_df_entropy_decomposition()
    // in https://pygor3.readthedocs.io/en/latest/Tutorial.html#Entropy
    //
    // IMPORTANT: pygor3 decomposes the total scenario entropy as
    //     H(scenario) = Σ_i H(X_i | Pa(X_i)),
    // i.e. the conditional entropy of each event given its Bayesian-
    // network parents. Our test computes the *marginal* entropy H(X_i)
    // via compute_event_marginal_probability(). By the information-
    // processing inequality, H(X) >= H(X|Pa), so for events WITH
    // parents our value is an upper bound on pygor3's.
    //
    // For combined insertion+dinuc entropies, a ~0.05 bit residual
    // exists because pygor3's get_P_stationary_state_dinucl() uses
    // np.linalg.eig(T) (right eigenvectors) instead of the transpose
    // np.linalg.eig(T.T) (left eigenvectors = true stationary dist).
    // For a row-stochastic T, the right eigenvector with eigenvalue 1
    // is always [1,1,1,1], yielding a uniform π=[0.25,0.25,0.25,0.25]
    // regardless of T.  This causes pygor3 to overestimate both H_ss
    // (always 2.0 bits = log₂4) and the entropy rate h.  Our power-
    // iteration approach computes the correct stationary distribution,
    // so our values are more accurate than pygor3's references.
    // ------------------------------------------------------------------
    if (!pygor3_reference_entropy.empty()) {
        std::cout << "\n=== Cross-validation against pygor3 (" << model_label << ") ===" << std::endl;
        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov)
                continue;

            // Determine the reported entropy for this event
            double reported_H;
            auto it_pair = ins_dinuc_pairs.find(ev.gene_class);
            bool is_combined_ins = it_pair != ins_dinuc_pairs.end() && it_pair->second.ins_event == &ev
                    && it_pair->second.dinuc_event != nullptr;
            if (is_combined_ins) {
                reported_H = it_pair->second.combined_H;
            } else {
                reported_H = ev.H;
            }

            auto it_ref = pygor3_reference_entropy.find(ev.nickname);
            if (it_ref != pygor3_reference_entropy.end()) {
                double ref = it_ref->second;
                double diff = std::abs(reported_H - ref);
                bool has_parents = !topology->parentsIds(
                    static_cast<igor::index_type>(ev.queue_position)).empty();

                std::cout << "  " << ev.nickname << ": computed=" << reported_H << "  pygor3=" << ref
                          << "  diff=" << diff << (has_parents ? "  (marginal; parents present)" : "") << std::endl;

                if (is_combined_ins) {
                    // Combined insertion + dinuc: pygor3 uses a uniform
                    // stationary distribution (bug in eig computation)
                    // which inflates its reference by ~0.05 bits.
                    // Our value is more accurate; allow 0.1 bit tolerance.
                    CHECK(diff < 0.1);
                } else if (has_parents) {
                    // Marginal entropy >= conditional entropy.
                    // The gap is bounded by I(X; Pa), the mutual
                    // information, typically < 1–2 bits for gene/del events.
                    CHECK(reported_H >= ref - 0.01);
                    CHECK(diff < 1.5);
                } else {
                    // Root event (no parents): marginal = conditional.
                    CHECK(diff < 0.01);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // 3. Generate scenarios for increasing N and track D_KL per event
    // ------------------------------------------------------------------

    // Update ev.H for insertion events to the full combined entropy
    // (insertion-length + dinucleotide Markov) so that D_KL bounds and
    // printouts reference the total information content of the pair.
    for (auto &ev : event_infos) {
        if (ev.is_dinuc_markov)
            continue;
        auto it = ins_dinuc_pairs.find(ev.gene_class);
        if (it != ins_dinuc_pairs.end() && it->second.ins_event
            && it->second.ins_event->queue_position == ev.queue_position && it->second.dinuc_event) {
            ev.H = it->second.combined_H;
        }
    }

    // Map: event UID → vector of (D_KL, uncovered_mass) per sample-size
    std::map<size_t, std::vector<std::pair<double, double>>> kl_traces;
    // Map: event UID → vector of empirical entropy per sample-size
    std::map<size_t, std::vector<double>> entropy_traces;
    // FastGenerator empirical marginals at N=1M, saved for cross-validation
    // against SamplingEngine (section 6).
    std::map<size_t, std::vector<double>> fg_empiricals_1M;

    const std::vector<int> sample_sizes = { 10, 100, 1000, 1000000 };

    // Seed the RNG for reproducibility
    std::mt19937_64 rng(42);

    for (int N : sample_sizes) {
        INFO("Generating " << N << " scenarios");

        // Generate N scenarios using SamplingEngine
        std::vector<SampledScenario> scenarios;
        scenarios.reserve(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) {
            scenarios.push_back(engine.run(rng));
        }

        REQUIRE(scenarios.size() == static_cast<size_t>(N));

        std::cout << "\n--- N = " << N << " ---" << std::endl;

        // Compute all empirical marginals in a single pass
        auto all_empiricals = compute_all_empirical_marginals(
            scenarios, event_infos, scenarios.size());

        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov) {
                continue; // DinucMarkov has variable-length realizations
            }

            const auto &empirical = all_empiricals.at(ev.queue_position);

            double uncovered = 0.0;
            double dkl = kl_divergence(ev.model_marginal, empirical, &uncovered);
            kl_traces[ev.queue_position].emplace_back(dkl, uncovered);

            // Empirical entropy of the generated marginal
            double H_emp = entropy(empirical);
            entropy_traces[ev.queue_position].push_back(H_emp);

            std::cout << "  " << ev.nickname << " : D_KL = " << dkl << "  H_emp = " << H_emp << "  H_theo = " << ev.H
                      << "  uncovered = " << uncovered << std::endl;

            // ---- assertions ----
            // D_KL must always be finite (partial KL can be negative when
            // uncovered bins are skipped, so we only check finiteness)
            CHECK(std::isfinite(dkl));

            if (N == 10) {
                // With only 10 samples, D_KL can be very large
                // and many bins will be uncovered.
                // No bound check at this sample size.
            } else if (N == 100) {
                // Loose bound: D_KL should be within the same order as H
                CHECK(dkl < (std::max)(ev.H, 1.0));
            } else if (N == 1000) {
                // Tighter: D_KL should be at most a tenth of H
                CHECK(dkl < (std::max)(ev.H, 1.0) * 0.1);
                // Most of the distribution should be covered
                CHECK(uncovered < 0.1);
            } else if (N == 1000000) {
                // Very tight: D_KL should be negligible compared to H.
                // Asymptotic expectation: D_KL ≈ (k−1)/(2·N·ln2).
                // The estimator variance follows χ²(k−1), so we use a
                // 50× safety factor to avoid flaky CI failures across
                // many events and two model sections.
                double expected_upper = 50.0 * (ev.num_realizations - 1) / (2.0 * N * std::log(2.0));
                CHECK(dkl < (std::max)(expected_upper, 1e-3));
                // At 1M samples, essentially all bins should be covered
                // (rare alleles with P~1e-6 may still be missed)
                CHECK(uncovered < 1e-3);
            }
        }

        // Save empiricals for cross-validation (overwritten each iteration;
        // the last iteration is N=1M which is the one we need).
        fg_empiricals_1M = all_empiricals;
    }

    // ------------------------------------------------------------------
    // 4. Empirical entropy convergence to the theoretical entropy
    // ------------------------------------------------------------------
    // The empirical entropy H(Q) of the generated marginal should
    // converge to the theoretical entropy H(P) of the model.  The
    // finite-sample estimator has a known negative bias:
    //     E[H_emp] ≈ H − (k−1) / (2·N·ln2)
    // so we check that the difference |H_emp − H| shrinks with N.
    // ------------------------------------------------------------------
    std::cout << "\n=== Empirical entropy vs theoretical ==" << std::endl;
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov)
            continue;

        const auto &htrace = entropy_traces.at(ev.queue_position);
        if (htrace.size() >= 4) {
            double H_emp_1M = htrace[3]; // index 3 → N=1 000 000
            // H_emp is the empirical entropy of the marginal distribution
            // for this event only (e.g. insertion-length distribution),
            // so it converges to H(model_marginal), not to the full
            // combined entropy ev.H (which includes dinuc Markov).
            double H_theo = entropy(ev.model_marginal);
            // The finite-sample entropy estimator has:
            //   bias  ≈ (k−1) / (2·N·ln2)
            //   stdev ≈ sqrt(Var[−log2 P]) / sqrt(N)
            // We bound stdev by sqrt(H²) / sqrt(N) = H/sqrt(N) as a
            // conservative upper estimate (Var ≤ E²).
            // Use 5σ + 50× bias for a CI-safe tolerance.
            double bias = 50.0 * (ev.num_realizations - 1) / (2.0 * 1e6 * std::log(2.0));
            double stdev_bound = 5.0 * H_theo / std::sqrt(1e6);
            double tolerance = bias + stdev_bound;
            double diff = std::abs(H_emp_1M - H_theo);

            std::cout << "  " << ev.nickname << " : H_theo=" << H_theo << "  H_emp(1M)=" << H_emp_1M
                      << "  diff=" << diff << "  tol=" << tolerance << std::endl;

            // The empirical entropy should be close to theoretical;
            // allow bias + 5σ variance, with a minimum floor.
            CHECK(diff < (std::max)(tolerance, 1e-3));
        }
    }

    // ------------------------------------------------------------------
    // 5. Check monotonic decrease of D_KL with increasing N
    // ------------------------------------------------------------------
    std::cout << "\n=== Monotonic decrease check ===" << std::endl;
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov)
            continue;

        const auto &trace = kl_traces.at(ev.queue_position);
        // Compare D_KL at N=1000 vs N=1'000'000.  Both are large enough
        // for the estimator to be reliable.  A small absolute tolerance
        // handles the case where the smaller sample happens to nail the
        // distribution exactly (D_KL ≈ 0) while the larger one is merely
        // very small.
        //
        // The partial KL (which skips uncovered bins) can go negative
        // when many bins are missed at smaller N, making the comparison
        // meaningless.  We only check monotonicity when uncovered mass
        // at N=1k is small enough for the partial KL to be a faithful
        // estimate of the true KL.
        if (trace.size() >= 4) {
            auto [dkl_at_1k, uncov_1k] = trace[2]; // index 2 → N=1 000
            auto [dkl_at_1M, uncov_1M] = trace[3]; // index 3 → N=1 000 000
            std::cout << "  " << ev.nickname << " : D_KL(N=1k)=" << dkl_at_1k << "  D_KL(N=1M)=" << dkl_at_1M
                      << "  uncov(1k)=" << uncov_1k << "  uncov(1M)=" << uncov_1M << std::endl;
            if (uncov_1k < 0.01) {
                // Both estimates are reliable; D_KL should decrease
                CHECK(dkl_at_1M < dkl_at_1k + 1e-4);
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. SamplingEngine – modern architecture cross-validation
    // ------------------------------------------------------------------
    // Validate that the modern SamplingEngine produces marginals
    // consistent with the theoretical model AND with the FastGenerator.
    // ------------------------------------------------------------------
    {
        using namespace igor::model;

        std::cout << "\n=== SamplingEngine cross-validation (" << model_label
                  << ") ===" << std::endl;

        // Build the modern pipeline
        auto se_model = recombination_model_from_files<double>(
            model_parms_path, model_marginals_path);
        REQUIRE(se_model.topology().size() == event_infos.size());

        auto se_model_ptr = std::make_shared<const RecombinationModel<double>>(
            std::move(se_model));
        SamplingEngine<double> engine(se_model_ptr);

        // Build nickname → UID map for cross-referencing
        const auto& se_topology = se_model_ptr->topology();
        std::map<std::string, igor::index_type> nick_to_uid;
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(se_topology.size()); ++uid) {
            nick_to_uid[se_topology.event(uid)->get_nickname()] = uid;
        }

        // Verify every legacy event has a matching topology node
        for (const auto &ev : event_infos) {
            INFO("Missing topology node for event: " << ev.nickname);
            CHECK(nick_to_uid.count(ev.nickname) == 1);
        }

        // Generate 1M scenarios
        constexpr std::size_t SE_N = 1'000'000;
        std::mt19937_64 rng(42);

        std::cout << "  Generating " << SE_N
                  << " scenarios with SamplingEngine ..." << std::flush;
        std::vector<SampledScenario> scenarios;
        scenarios.reserve(SE_N);
        for (std::size_t i = 0; i < SE_N; ++i) {
            scenarios.push_back(engine.run(rng));
        }
        std::cout << " done." << std::endl;

        auto se_empiricals = compute_all_empirical_marginals(
            scenarios, event_infos, nick_to_uid);

        // D_KL checks and cross-validation against FastGenerator
        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov)
                continue;

            auto uid_it = nick_to_uid.find(ev.nickname);
            REQUIRE(uid_it != nick_to_uid.end());
            igor::index_type uid = uid_it->second;

            auto emp_it = se_empiricals.find(uid);
            REQUIRE(emp_it != se_empiricals.end());
            const auto &se_emp = emp_it->second;

            // --- Check vs theoretical model marginal ---
            double uncovered = 0.0;
            double dkl = kl_divergence(ev.model_marginal, se_emp, &uncovered);
            double expected_upper = 50.0 * (ev.num_realizations - 1)
                                    / (2.0 * static_cast<double>(SE_N) * std::log(2.0));

            INFO("SamplingEngine N=1M  Event: " << ev.nickname
                     << "  D_KL=" << dkl << "  uncovered=" << uncovered
                     << "  upper=" << expected_upper);
            CHECK(dkl < (std::max)(expected_upper, 1e-3));
            CHECK(uncovered < 1e-3);

            // --- Cross-validate vs FastGenerator empirical at N=1M ---
            // Both generators sample from the same model, so at 1M samples
            // the two empirical marginals should be nearly identical.
            // D_KL between two empirical samples of size N has
            //   E[D_KL] ≈ (k−1) / (N·ln2)  (doubled since both are estimates).
            const auto &fg_emp = fg_empiricals_1M.at(ev.queue_position);
            double dkl_cross = kl_divergence(fg_emp, se_emp);
            double cross_upper = 100.0 * (ev.num_realizations - 1)
                                 / (2.0 * static_cast<double>(SE_N) * std::log(2.0));

            std::cout << "  " << ev.nickname
                      << " : D_KL(model||SE)=" << dkl
                      << "  D_KL(FG||SE)=" << dkl_cross
                      << "  uncov=" << uncovered << std::endl;

            INFO("Cross-validation D_KL(FG||SE)=" << dkl_cross);
            CHECK(dkl_cross < (std::max)(cross_upper, 1e-3));
        }
    }
}

TEST_CASE("SamplingEngine marginals converge - modern architecture",
          "[generation][SamplingEngine]")
{
    using namespace igor::model;

    // ------------------------------------------------------------------
    // 1. Select model (same two-SECTION structure as the legacy test)
    // ------------------------------------------------------------------
    std::string model_parms_path;
    std::string model_marginals_path;
    std::string model_label;
    int min_events = 0;

    SECTION("human TCR alpha (VJ)")
    {
        model_parms_path     = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        model_label          = "TCR alpha";
        min_events           = 4;
    }
    SECTION("human TCR beta (VDJ)")
    {
        model_parms_path     = MODELS_DIR + "/human/tcr_beta/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_beta/models/model_marginals.txt";
        model_label          = "TCR beta";
        min_events           = 8;
    }

    std::cout << "\n=== SamplingEngine generation test: " << model_label
              << " ===" << std::endl;

    // ------------------------------------------------------------------
    // 2. Legacy pipeline – for reference marginals
    // ------------------------------------------------------------------
    Model_Parms parms;
    parms.read_model_parms(model_parms_path);
    Model_marginals marginals(parms);
    marginals.txt2marginals(model_marginals_path, parms);
    auto event_infos = build_event_info(parms, marginals);
    REQUIRE(event_infos.size() >= static_cast<std::size_t>(min_events));

    // ------------------------------------------------------------------
    // 3. Modern pipeline – RecombinationModel + SamplingEngine
    // ------------------------------------------------------------------
    auto model_opt = recombination_model_from_files<double>(
        model_parms_path, model_marginals_path);
    REQUIRE(model_opt.topology().size() == event_infos.size());

    auto model_ptr = std::make_shared<const RecombinationModel<double>>(
        std::move(model_opt));
    SamplingEngine<double> engine(model_ptr);

    // Build nickname → UID map for cross-referencing
    const auto& topology = model_ptr->topology();
    std::map<std::string, igor::index_type> nick_to_uid;
    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(topology.size()); ++uid) {
        auto ev = topology.event(uid);
        nick_to_uid[ev->get_nickname()] = uid;
    }

    // Also build nickname → EventInfo* for quick lookup
    std::map<std::string, const EventInfo*> info_by_nick;
    for (const auto& ev : event_infos) {
        info_by_nick[ev.nickname] = &ev;
    }

    // Verify every legacy event has a matching topology node
    for (const auto& ev : event_infos) {
        INFO("Missing topology node for event: " << ev.nickname);
        CHECK(nick_to_uid.count(ev.nickname) == 1);
    }

    // ------------------------------------------------------------------
    // 4. Generate scenarios with SamplingEngine
    // ------------------------------------------------------------------
    // Sample sizes: {1k, 100k, 1M}
    // Much faster than the legacy test (no queue overhead) so we can
    // afford bigger runs and still get cleaner statistics.
    constexpr std::array<std::size_t, 3> sample_sizes = {1'000, 100'000, 1'000'000};

    // Track D_KL evolution per event across sample sizes
    // key = nickname, value = vector of (D_KL, uncovered) pairs
    std::map<std::string, std::vector<std::pair<double, double>>> traces;

    std::mt19937_64 rng(42);

    for (std::size_t N : sample_sizes) {
        std::cout << "  Generating " << N << " scenarios ..." << std::flush;

        std::vector<SampledScenario> scenarios;
        scenarios.reserve(N);
        for (std::size_t i = 0; i < N; ++i) {
            scenarios.push_back(engine.run(rng));
        }
        std::cout << " done." << std::endl;

        auto empiricals = compute_all_empirical_marginals(
            scenarios, event_infos, nick_to_uid);

        // ----------------------------------------------------------
        // 5. Compare empirical vs theoretical for each event
        // ----------------------------------------------------------
        for (const auto& ev : event_infos) {
            if (ev.is_dinuc_markov) continue;

            auto uid_it = nick_to_uid.find(ev.nickname);
            REQUIRE(uid_it != nick_to_uid.end());
            igor::index_type uid = uid_it->second;

            auto emp_it = empiricals.find(uid);
            REQUIRE(emp_it != empiricals.end());
            const auto& empirical = emp_it->second;

            double uncovered = 0.0;
            double dkl = kl_divergence(ev.model_marginal, empirical, &uncovered);

            traces[ev.nickname].emplace_back(dkl, uncovered);

            // Asymptotic expected upper bound:
            //   E[D_KL] ≈ (k-1) / (2·N·ln2)
            // Use a generous safety factor (50×) for statistical fluctuations.
            double expected_upper = 50.0 * static_cast<double>(ev.num_realizations - 1)
                                    / (2.0 * static_cast<double>(N) * std::log(2.0));

            // At N=1k the bound may be loose; at 100k+ it should be tight.
            INFO("N=" << N << "  Event: " << ev.nickname
                      << "  D_KL=" << dkl << "  uncovered=" << uncovered
                      << "  upper=" << expected_upper);

            if (N >= 100'000) {
                CHECK(dkl < std::max(expected_upper, 1e-3));
                CHECK(uncovered < 0.01);
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. Monotonic decrease of D_KL with increasing N
    // ------------------------------------------------------------------
    std::cout << "  D_KL convergence summary:" << std::endl;
    for (const auto& ev : event_infos) {
        if (ev.is_dinuc_markov) continue;

        auto& trace = traces[ev.nickname];
        if (trace.size() >= 3) {
            auto [dkl_1k,  uncov_1k]  = trace[0]; // N=1k
            auto [dkl_100k, uncov_100k] = trace[1]; // N=100k
            auto [dkl_1M,  uncov_1M]  = trace[2]; // N=1M

            std::cout << "    " << ev.nickname
                      << " : D_KL(1k)=" << dkl_1k
                      << "  D_KL(100k)=" << dkl_100k
                      << "  D_KL(1M)=" << dkl_1M
                      << "  uncov(1k)=" << uncov_1k
                      << "  uncov(1M)=" << uncov_1M << std::endl;

            // With 100k→1M, D_KL should decrease when estimate is reliable
            if (uncov_100k < 0.01) {
                CHECK(dkl_1M < dkl_100k + 1e-4);
            }
        }
    }
}
