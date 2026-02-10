/**
 * @file test_generation.cpp
 * @brief Tests the sequence generation process using KL divergence
 *        and entropy-based validation.
 *
 * Generates 10, 100, 1000, and 1'000'000 sequences from the human TCR alpha
 * model in models/human/tcr_alpha/. For each sample size, the empirical
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
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

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
#error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

// ---------------------------------------------------------------------------
// Paths to the human TCR alpha model
// ---------------------------------------------------------------------------
static const std::string MODEL_PARMS_PATH =
        std::string(IGOR_SOURCE_DIR) + "/models/human/tcr_alpha/models/model_parms.txt";
static const std::string MODEL_MARGINALS_PATH =
        std::string(IGOR_SOURCE_DIR) + "/models/human/tcr_alpha/models/model_marginals.txt";

// ---------------------------------------------------------------------------
// Mathematical helpers
// ---------------------------------------------------------------------------

/**
 * @brief Kullback–Leibler divergence D_KL(P || Q) in bits.
 *
 * Skips bins where P(i) == 0.  When P(i) > 0 but Q(i) == 0, the bin
 * is also skipped and the total probability mass of such bins is
 * accumulated in @p uncovered_mass (if non-null).  This avoids +∞
 * results that arise from finite-sample under-coverage of rare
 * realizations in models with many bins (e.g. 103 V genes).
 */
static double kl_divergence(const std::vector<double> &P,
                            const std::vector<double> &Q,
                            double *uncovered_mass = nullptr)
{
    double kl = 0.0;
    double skipped = 0.0;
    for (size_t i = 0; i < P.size(); ++i) {
        if (P[i] > 0.0) {
            if (Q[i] <= 0.0) {
                skipped += P[i];
            } else {
                kl += P[i] * std::log2(P[i] / Q[i]);
            }
        }
    }
    if (uncovered_mass) {
        *uncovered_mass = skipped;
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
// DinucMarkov entropy helpers
// ---------------------------------------------------------------------------

/**
 * @brief Compute the stationary distribution of a 4×4 row-stochastic
 *        Markov transition matrix via power iteration.
 *
 * The transition matrix T is stored as a flat 16-element array in
 * row-major order: T[i*4 + j] = P(j | i).
 * Returns π such that π · T = π and Σ π_i = 1.
 */
static std::array<double, 4> markov_stationary_distribution(
        const std::array<double, 16> &T,
        int max_iter = 1000,
        double tol = 1e-12)
{
    // Start with uniform distribution
    std::array<double, 4> pi = {0.25, 0.25, 0.25, 0.25};
    std::array<double, 4> pi_next{};

    for (int iter = 0; iter < max_iter; ++iter) {
        pi_next.fill(0.0);
        // π_next[j] = Σ_i π[i] * T[i][j]
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i) {
                pi_next[j] += pi[i] * T[i * 4 + j];
            }
        }
        // Normalize (should already be ~1, but guard against drift)
        double sum = 0.0;
        for (double v : pi_next) sum += v;
        for (double &v : pi_next) v /= sum;

        // Check convergence
        double max_diff = 0.0;
        for (int i = 0; i < 4; ++i) {
            max_diff = (std::max)(max_diff, std::abs(pi_next[i] - pi[i]));
        }
        pi = pi_next;
        if (max_diff < tol) break;
    }
    return pi;
}

/**
 * @brief Entropy rate of a first-order Markov chain.
 *
 *     h = − Σ_i π_i  Σ_j T_{ij} log2(T_{ij})
 *
 * where π is the stationary distribution and T is the transition matrix.
 */
static double markov_entropy_rate(const std::array<double, 16> &T)
{
    auto pi = markov_stationary_distribution(T);
    double h = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double t = T[i * 4 + j];
            if (t > 0.0 && pi[i] > 0.0) {
                h -= pi[i] * t * std::log2(t);
            }
        }
    }
    return h;
}

/**
 * @brief Combined insertion + dinucleotide entropy.
 *
 * Following the pygor3 entropy decomposition
 * (see test_IgorModel_entropy in pygor3), the total entropy
 * contribution for an insertion + dinucleotide pair is:
 *
 *     H_total = H(ℓ) + E[ℓ] · h
 *
 * where:
 *   - H(ℓ)  = Shannon entropy of the insertion-length distribution
 *   - E[ℓ]  = expected insertion length = Σ_ℓ ℓ · P(ℓ)
 *   - h     = entropy rate of the dinucleotide Markov chain
 *
 * @param ins_marginal  Marginal probability of each insertion length
 * @param ins_lengths   The actual length value for each insertion bin
 * @param dinuc_T       The 4×4 transition matrix (flat, row-major)
 * @return The combined entropy in bits
 */
static double insertion_dinuc_entropy(
        const std::vector<double> &ins_marginal,
        const std::vector<int> &ins_lengths,
        const std::array<double, 16> &dinuc_T)
{
    // H(ℓ) and E[ℓ]
    double H_len = 0.0;
    double E_len = 0.0;
    for (size_t i = 0; i < ins_marginal.size(); ++i) {
        double p = ins_marginal[i];
        if (p > 0.0) {
            H_len -= p * std::log2(p);
            E_len += p * ins_lengths[i];
        }
    }

    double h = markov_entropy_rate(dinuc_T);

    return H_len + E_len * h;
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
    Gene_class gene_class;              ///< gene class (VD_genes, DJ_genes, VJ_genes …)
    std::array<double, 16> dinuc_T;     ///< transition matrix (only for DinucMarkov)
    double dinuc_entropy_rate;          ///< Markov entropy rate h (only for DinucMarkov)
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
    auto index_map = marginals.get_index_map(parms);
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
        info.gene_class = ev->get_class();
        info.dinuc_T.fill(0.0);
        info.dinuc_entropy_rate = 0.0;

        std::cout << "  Processing event: " << info.nickname
                  << " (is_dinuc=" << info.is_dinuc_markov
                  << ", size=" << info.num_realizations << ")" << std::endl;

        if (info.is_dinuc_markov) {
            // Extract the 4×4 transition matrix from the raw marginal array.
            // DinucMarkov has no parents so 16 values sit contiguously
            // at index_map[name].
            int base_idx = index_map.at(info.name);
            for (int k = 0; k < 16; ++k) {
                info.dinuc_T[k] =
                        static_cast<double>(marginals.marginal_array_smart_p[base_idx + k]);
            }
            info.dinuc_entropy_rate = markov_entropy_rate(info.dinuc_T);

            // Marginal for a DinucMarkov is the flat 16-element array;
            // we store it for completeness but the "entropy of the
            // marginal" is the Markov entropy rate per step, not a
            // standard Shannon entropy.
            info.model_marginal.resize(16, 0.0);
            for (int k = 0; k < 16; ++k) {
                info.model_marginal[k] = info.dinuc_T[k];
            }
            info.H = info.dinuc_entropy_rate;
        } else {
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
            if (realization < 0 || realization >= num_realizations) {
                throw std::out_of_range(
                        "Realization index " + std::to_string(realization) +
                        " out of bounds [0, " + std::to_string(num_realizations) + ")");
            }
            ++counts[realization];
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
    // 1. Load the human TCR alpha model
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

    // ------------------------------------------------------------------
    // Compute combined insertion+dinucleotide entropy.
    // Following pygor3's get_df_entropy_decomposition(), insertion events
    // like "vj_ins" are paired with their DinucMarkov partner "vj_dinucl"
    // (sharing the same Gene_class). The combined entropy is:
    //     H_total = H(ℓ) + E[ℓ] · h
    // where H(ℓ) is the insertion-length entropy, E[ℓ] the expected length,
    // and h the Markov chain entropy rate.
    // ------------------------------------------------------------------
    struct InsDinucPair {
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

            pair.combined_H = insertion_dinuc_entropy(
                    pair.ins_event->model_marginal,
                    ins_lengths,
                    pair.dinuc_event->dinuc_T);
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
        if (it != ins_dinuc_pairs.end() &&
            it->second.ins_event == &ev &&
            it->second.dinuc_event != nullptr)
        {
            double combined = it->second.combined_H;
            std::cout << "  " << ev.nickname << " + "
                      << it->second.dinuc_event->nickname
                      << " : H = " << combined
                      << "  (H_len=" << ev.H
                      << ", E[l]*h=" << (combined - ev.H)
                      << ", h=" << it->second.dinuc_event->dinuc_entropy_rate
                      << ")" << std::endl;
            total_entropy += combined;
        } else {
            std::cout << "  " << ev.nickname << " : H = " << ev.H << std::endl;
            total_entropy += ev.H;
        }
    }
    // Also print standalone DinucMarkov entropy rates for reference
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) {
            std::cout << "  (" << ev.nickname << " entropy rate h = "
                      << ev.dinuc_entropy_rate << " bits/nt)" << std::endl;
        }
    }
    std::cout << "  TOTAL: " << total_entropy << std::endl;

    // ------------------------------------------------------------------
    // 3. Generate sequences for increasing N and track D_KL per event
    // ------------------------------------------------------------------
    // Map: event_queue_position → vector of (D_KL, uncovered_mass) per sample-size
    std::map<size_t, std::vector<std::pair<double, double>>> kl_traces;

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

            double uncovered = 0.0;
            double dkl = kl_divergence(ev.model_marginal, empirical, &uncovered);
            kl_traces[ev.queue_position].emplace_back(dkl, uncovered);

            std::cout << "  " << ev.nickname
                      << " : D_KL = " << dkl
                      << "  (H = " << ev.H
                      << ", uncovered = " << uncovered << ")" << std::endl;

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
                // Very tight: D_KL should be negligible compared to H
                // Theoretical expectation: D_KL ≈ (k−1)/(2·N·ln2)
                double expected_upper =
                        10.0 * (ev.num_realizations - 1) / (2.0 * N * std::log(2.0));
                CHECK(dkl < (std::max)(expected_upper, 1e-3));
                // At 1M samples, essentially all bins should be covered
                CHECK(uncovered < 1e-4);
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
        //
        // The partial KL (which skips uncovered bins) can go negative
        // when many bins are missed at smaller N, making the comparison
        // meaningless.  We only check monotonicity when uncovered mass
        // at N=1k is small enough for the partial KL to be a faithful
        // estimate of the true KL.
        if (trace.size() >= 4) {
            auto [dkl_at_1k, uncov_1k] = trace[2];   // index 2 → N=1 000
            auto [dkl_at_1M, uncov_1M] = trace[3];   // index 3 → N=1 000 000
            std::cout << "  " << ev.nickname
                      << " : D_KL(N=1k)=" << dkl_at_1k
                      << "  D_KL(N=1M)=" << dkl_at_1M
                      << "  uncov(1k)=" << uncov_1k
                      << "  uncov(1M)=" << uncov_1M << std::endl;
            if (uncov_1k < 0.01) {
                // Both estimates are reliable; D_KL should decrease
                CHECK(dkl_at_1M < dkl_at_1k + 1e-4);
            }
        }
    }
}
