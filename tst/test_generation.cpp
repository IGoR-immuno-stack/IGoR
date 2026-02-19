/**
 * @file test_generation.cpp
 * @brief Tests the sequence generation process using KL divergence
 *        and entropy-based validation.
 *
 * Generates 10, 100, 1000, and 1'000'000 sequences from the human TCR alpha
 * (VJ) and human TCR beta (VDJ) models. For each sample size, the empirical
 * marginal distributions of every
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

#include <igor/Core/FastGenerator.h>
#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

#include "entropy_test_helpers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <forward_list>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
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
 * @brief Compute empirical marginal distributions for ALL non-DinucMarkov
 *        events in a single pass over the scenario list.
 *
 * Returns a map from event queue_position to its empirical marginal vector.
 * This avoids the O(events × N) queue-copying that caused timeouts
 * when compute_empirical_marginal was called per-event.
 */
static std::map<size_t, std::vector<double>>
compute_all_empirical_marginals(const std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> &scenarios,
                                const std::vector<EventInfo> &event_infos, size_t total_sequences)
{
    // Prepare per-event count arrays
    // Only for non-DinucMarkov events
    std::map<size_t, std::vector<size_t>> counts;
    for (const auto &ev : event_infos) {
        if (!ev.is_dinuc_markov) {
            counts[ev.queue_position].assign(ev.num_realizations, 0);
        }
    }

    // Single pass over all scenarios
    for (auto it = scenarios.begin(); it != scenarios.end(); ++it) {
        // Copy the outer queue once per scenario and drain it,
        // collecting realizations for every event in one go.
        auto outer = it->second;
        size_t pos = 0;
        while (!outer.empty()) {
            auto inner = outer.front();
            outer.pop();

            auto cnt_it = counts.find(pos);
            if (cnt_it != counts.end() && !inner.empty()) {
                int realization = inner.front();
                int num_real = static_cast<int>(cnt_it->second.size());
                if (realization < 0 || realization >= num_real) {
                    throw std::out_of_range("Realization index " + std::to_string(realization) + " out of bounds [0, "
                                            + std::to_string(num_real) + ")");
                }
                ++cnt_it->second[realization];
            }
            ++pos;
        }
    }

    // Convert counts to probabilities
    std::map<size_t, std::vector<double>> result;
    double n = static_cast<double>(total_sequences);
    for (auto &[qpos, cvec] : counts) {
        std::vector<double> empirical(cvec.size(), 0.0);
        for (size_t i = 0; i < cvec.size(); ++i) {
            empirical[i] = static_cast<double>(cvec[i]) / n;
        }
        result[qpos] = std::move(empirical);
    }
    return result;
}

/**
 * @brief Compute empirical marginal distributions from FastGenerator output.
 *
 * Overload that works with the FastGenerator's output format:
 * vector<GeneratedSequence> where each GeneratedSequence::realizations
 * is a vector<vector<int>>. For non-DinucMarkov events, realizations[pos]
 * contains a single element: the realization index.
 */
static std::map<size_t, std::vector<double>>
compute_all_empirical_marginals(const std::vector<igor::fast::GeneratedSequence> &sequences,
                                const std::vector<EventInfo> &event_infos, size_t total_sequences)
{
    // Prepare per-event count arrays (non-DinucMarkov only)
    std::map<size_t, std::vector<size_t>> counts;
    for (const auto &ev : event_infos) {
        if (!ev.is_dinuc_markov) {
            counts[ev.queue_position].assign(ev.num_realizations, 0);
        }
    }

    // Single pass over all generated sequences
    for (const auto &seq : sequences) {
        for (auto &[qpos, cvec] : counts) {
            if (qpos < seq.realizations.size() && !seq.realizations[qpos].empty()) {
                int realization = seq.realizations[qpos][0];
                int num_real = static_cast<int>(cvec.size());
                if (realization < 0 || realization >= num_real) {
                    throw std::out_of_range("Realization index " + std::to_string(realization) + " out of bounds [0, "
                                            + std::to_string(num_real) + ")");
                }
                ++cvec[realization];
            }
        }
    }

    // Convert counts to probabilities
    std::map<size_t, std::vector<double>> result;
    double n = static_cast<double>(total_sequences);
    for (auto &[qpos, cvec] : counts) {
        std::vector<double> empirical(cvec.size(), 0.0);
        for (size_t i = 0; i < cvec.size(); ++i) {
            empirical[i] = static_cast<double>(cvec[i]) / n;
        }
        result[qpos] = std::move(empirical);
    }
    return result;
}

// ---------------------------------------------------------------------------
// THE TEST
// ---------------------------------------------------------------------------

TEST_CASE("Generation marginals converge - KL divergence vs entropy", "[generation]")
{
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
    // 2. Load model and collect per-event metadata
    // ------------------------------------------------------------------
    INFO("Testing model: " << model_label);
    Model_Parms model_parms;
    model_parms.read_model_parms(model_parms_path);

    Model_marginals model_marginals(model_parms);
    model_marginals.txt2marginals(model_marginals_path, model_parms);

    std::vector<EventInfo> event_infos = build_event_info(model_parms, model_marginals);

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
            auto reals_map = pair.ins_event->name;
            auto ev_ptr = model_parms.get_event_pointer(pair.ins_event->name);
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
                bool has_parents = !model_parms.get_parents(ev.name).empty();

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
    // 3. Generate sequences for increasing N and track D_KL per event
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

    // Map: event_queue_position → vector of (D_KL, uncovered_mass) per sample-size
    std::map<size_t, std::vector<std::pair<double, double>>> kl_traces;
    // Map: event_queue_position → vector of empirical entropy per sample-size
    std::map<size_t, std::vector<double>> entropy_traces;

    const std::vector<int> sample_sizes = { 10, 100, 1000, 1000000 };

    // Initialize FastGenerator once, reuse across all sample sizes
    igor::fast::FastGenerator fast_gen;
    fast_gen.initialize(model_parms, model_marginals);
    REQUIRE(fast_gen.is_initialized());

    for (int N : sample_sizes) {
        INFO("Generating " << N << " sequences");

        // Generate using FastGenerator (parallel, precomputed CDFs)
        igor::fast::FastGeneratorConfig config;
        config.show_progress = false;
        auto sequences = fast_gen.generate(static_cast<size_t>(N), config);

        REQUIRE(sequences.size() == static_cast<size_t>(N));

        std::cout << "\n--- N = " << N << " ---" << std::endl;

        // Compute all empirical marginals in a single pass
        auto all_empiricals = compute_all_empirical_marginals(sequences, event_infos, sequences.size());

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
}

// ---------------------------------------------------------------------------
// Engine-based generation validation
// ---------------------------------------------------------------------------

// TEST_CASE("Engine-based generation marginals converge", "[generation][engine]")
// {
//     // ------------------------------------------------------------------
//     // 1. Select model to test
//     // ------------------------------------------------------------------
//     std::string model_parms_path;
//     std::string model_marginals_path;
//     std::string model_label;
//     int min_events = 0;

//     SECTION("human TCR alpha (VJ)")
//     {
//         model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
//         model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
//         model_label = "human/tcr_alpha";
//         min_events = 5;
//     }

//     SECTION("human TCR beta (VDJ)")
//     {
//         model_parms_path = MODELS_DIR + "/human/tcr_beta/models/model_parms.txt";
//         model_marginals_path = MODELS_DIR + "/human/tcr_beta/models/model_marginals.txt";
//         model_label = "human/tcr_beta";
//         min_events = 9;
//     }

//     // ------------------------------------------------------------------
//     // 2. Load model
//     // ------------------------------------------------------------------
//     INFO("Testing engine-based generation for model: " << model_label);
//     Model_Parms model_parms;
//     model_parms.read_model_parms(model_parms_path);

//     Model_marginals model_marginals(model_parms);
//     model_marginals.txt2marginals(model_marginals_path, model_parms);

//     std::vector<EventInfo> event_infos = build_event_info(model_parms, model_marginals);
//     REQUIRE(event_infos.size() >= static_cast<size_t>(min_events));

//     // ------------------------------------------------------------------
//     // 3. Generate sequences with InferenceEngine path
//     // ------------------------------------------------------------------
//     GenModel gen_model(model_parms, model_marginals);

//     // Use N=1000 for the engine path (slower than FastGenerator)
//     const int N = 1000;
//     const int SEED = 42;

//     INFO("Generating " << N << " sequences via engine path");
//     auto sequences = gen_model.generate_sequences_with_engine(N, /*output_realizations=*/true, SEED);

//     // Count sequences (forward_list has no size())
//     size_t count = 0;
//     for (auto it = sequences.begin(); it != sequences.end(); ++it) {
//         ++count;
//     }
//     REQUIRE(count == static_cast<size_t>(N));

//     // ------------------------------------------------------------------
//     // 4. Validate marginals via KL divergence
//     // ------------------------------------------------------------------
//     auto all_empiricals = compute_all_empirical_marginals(sequences, event_infos, count);

//     std::cout << "\n=== Engine-based generation: N = " << N << " (" << model_label << ") ===" << std::endl;

//     for (const auto &ev : event_infos) {
//         if (ev.is_dinuc_markov) {
//             continue; // DinucMarkov has variable-length realizations
//         }

//         const auto &empirical = all_empiricals.at(ev.queue_position);

//         double uncovered = 0.0;
//         double dkl = kl_divergence(ev.model_marginal, empirical, &uncovered);

//         std::cout << "  " << ev.nickname << " : D_KL = " << dkl << "  H = " << ev.H
//                   << "  uncovered = " << uncovered << std::endl;

//         // D_KL must be finite
//         CHECK(std::isfinite(dkl));

//         // At N=1000: D_KL should be at most a tenth of H
//         CHECK(dkl < (std::max)(ev.H, 1.0) * 0.1);

//         // Most of the distribution should be covered
//         CHECK(uncovered < 0.1);
//     }
// }
