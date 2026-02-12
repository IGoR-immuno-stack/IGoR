/**
 * @file entropy_test_helpers.h
 * @brief Common helper functions for IGoR tests
 * 
 * Provides mathematical utilities for KL divergence, entropy calculations,
 * and event metadata handling shared across generation and inference tests.
 */

#pragma once

#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

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
static inline double kl_divergence(const std::vector<double> &P,
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
static inline double entropy(const std::vector<double> &P)
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
static inline std::array<double, 4> markov_stationary_distribution(
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
 * @brief Cross-entropy between two first-order Markov chain transition matrices.
 *
 *     H(P, Q) = − Σ_i π_P[i]  Σ_j P_{ij} log2(Q_{ij})
 *
 * where π_P is the stationary distribution of P (the "true" distribution),
 * P is the true transition matrix, and Q is the predicted/inferred matrix.
 *
 * Uses epsilon (1e-100) for numerical stability if Q[i,j] = 0 when P[i,j] > 0.
 */
static inline double markov_cross_entropy(
        const std::array<double, 16> &P,
        const std::array<double, 16> &Q,
        const std::array<double, 4> &pi_P)
{
    constexpr double EPSILON_PROB = 1e-100;
    double h = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (pi_P[i] <= 0.0) continue;
        for (int j = 0; j < 4; ++j) {
            double p = P[i * 4 + j];
            double q = Q[i * 4 + j];
            if (p > 0.0) {
                // Use epsilon for numerical stability if q=0
                if (q <= 0.0) {
                    q = EPSILON_PROB;
                }
                h -= pi_P[i] * p * std::log2(q);
            }
        }
    }
    return h;
}

/**
 * @brief Entropy rate of a first-order Markov chain.
 *
 *     h(P) = H(P, P) = − Σ_i π_i  Σ_j P_{ij} log2(P_{ij})
 *
 * where π is the stationary distribution and P is the transition matrix.
 */
static inline double markov_entropy_rate(const std::array<double, 16> &P)
{
    auto pi = markov_stationary_distribution(P);
    return markov_cross_entropy(P, P, pi);
}

/**
 * @brief KL divergence between two first-order Markov chain transition matrices.
 *
 *     D_KL(P || Q) = H(P, Q) - h(P)
 *                  = Σ_i π_P[i] Σ_j P_{ij} log2(P_{ij} / Q_{ij})
 *
 * where π_P is the stationary distribution of P (the "true" distribution).
 */
static inline double markov_kl_divergence(
        const std::array<double, 16> &P,
        const std::array<double, 16> &Q)
{
    auto pi_P = markov_stationary_distribution(P);
    double cross_ent = markov_cross_entropy(P, Q, pi_P);
    double ent = markov_cross_entropy(P, P, pi_P);
    return cross_ent - ent;
}

/**
 * @brief Combined insertion + dinucleotide entropy.
 *
 * Following pygor3's get_df_Insertion_entropy_contribution() and
 * get_conditional_entropy_dinucl_function_l_ins(), the total entropy
 * contribution for an insertion + DinucMarkov pair is:
 *
 *     H_total = H(ℓ) + Σ_ℓ P(ℓ) · f(ℓ)
 *
 * where f(ℓ) is the entropy of the dinucleotide Markov chain of
 * length ℓ, which (under the stationary-start approximation) gives:
 *
 *     f(ℓ=0) = H_ss           (pygor3 convention; strictly should be 0)
 *     f(ℓ≥1) = H_ss + (ℓ−1)·h
 *
 * with:
 *   - H_ss = entropy of the Markov chain's stationary distribution π
 *   - h    = entropy rate = −Σ_i π_i Σ_j T_ij log₂ T_ij
 *
 * This simplifies to:
 *
 *     H_total = H(ℓ) + H_ss + h · [E[ℓ] − (1 − P(0))]
 *
 * @param ins_marginal  Marginal probability of each insertion length
 * @param ins_lengths   The actual length value for each insertion bin
 * @param dinuc_T       The 4×4 transition matrix (flat, row-major)
 * @return The combined entropy in bits
 */
static inline double insertion_dinuc_entropy(
        const std::vector<double> &ins_marginal,
        const std::vector<int> &ins_lengths,
        const std::array<double, 16> &dinuc_T)
{
    // H(ℓ) = Shannon entropy of the insertion-length distribution
    // E[ℓ] = expected insertion length
    // P0   = P(ℓ=0)
    double H_len = 0.0;
    double E_len = 0.0;
    double P0 = 0.0;
    for (size_t i = 0; i < ins_marginal.size(); ++i) {
        double p = ins_marginal[i];
        if (p > 0.0) {
            H_len -= p * std::log2(p);
            E_len += p * ins_lengths[i];
            if (ins_lengths[i] == 0) P0 += p;
        }
    }

    // Markov chain entropy rate and stationary distribution entropy
    auto pi = markov_stationary_distribution(dinuc_T);
    double h = 0.0;
    double H_ss = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (pi[i] > 0.0) {
            H_ss -= pi[i] * std::log2(pi[i]);
        }
        for (int j = 0; j < 4; ++j) {
            double t = dinuc_T[i * 4 + j];
            if (t > 0.0 && pi[i] > 0.0) {
                h -= pi[i] * t * std::log2(t);
            }
        }
    }

    // pygor3 formula:
    //   H_total = H(ℓ) + H_ss + h · [E[ℓ] − (1 − P(0))]
    return H_len + H_ss + h * (E_len - (1.0 - P0));
}

/**
 * @brief Cross-entropy between two insertion + dinucleotide Markov models.
 *
 * Computes H(P, Q) where P is the true model and Q is the predicted model.
 * Follows the same decomposition as insertion_dinuc_entropy:
 *
 *     H(P, Q) = H(P_len, Q_len) + H_ss(P, Q) + h_cross(P, Q) · [E_P[ℓ] − (1 − P_P(0))]
 *
 * where:
 *   - H(P_len, Q_len) = cross-entropy of insertion length distributions
 *   - H_ss(P, Q) = cross-entropy of stationary distributions
 *   - h_cross(P, Q) = cross-entropy rate of Markov chains
 *
 * @return The cross-entropy in bits
 */
static inline double insertion_dinuc_cross_entropy(
        const std::vector<double> &P_ins_marginal,
        const std::vector<double> &Q_ins_marginal,
        const std::vector<int> &ins_lengths,
        const std::array<double, 16> &P_dinuc_T,
        const std::array<double, 16> &Q_dinuc_T)
{
    // Use tiny epsilon for numerical stability when Q has zeros
    // This gives finite but very large penalty: log₂(1e-100) ≈ -332 bits
    constexpr double EPSILON_PROB = 1e-100;
    
    // H(P_len, Q_len) = cross-entropy of insertion-length distributions
    // E_P[ℓ] = expected insertion length under P
    // P_P(0) = P(ℓ=0) under P
    double H_len_cross = 0.0;
    double E_len_P = 0.0;
    double P0_P = 0.0;
    for (size_t i = 0; i < P_ins_marginal.size(); ++i) {
        double p = P_ins_marginal[i];
        double q = Q_ins_marginal[i];
        if (p > 0.0) {
            // Use epsilon for numerical stability if q=0
            if (q <= 0.0) {
                q = EPSILON_PROB;
            }
            H_len_cross -= p * std::log2(q);
            E_len_P += p * ins_lengths[i];
            if (ins_lengths[i] == 0) P0_P += p;
        }
    }

    // Stationary distribution cross-entropy: H(π_P, π_Q)
    auto pi_P = markov_stationary_distribution(P_dinuc_T);
    auto pi_Q = markov_stationary_distribution(Q_dinuc_T);
    
    double H_ss_cross = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (pi_P[i] > 0.0) {
            double q = pi_Q[i];
            // Use epsilon for numerical stability if q=0
            if (q <= 0.0) {
                q = EPSILON_PROB;
            }
            H_ss_cross -= pi_P[i] * std::log2(q);
        }
    }

    // Markov chain cross-entropy rate: h_cross(P, Q)
    double h_cross = markov_cross_entropy(P_dinuc_T, Q_dinuc_T, pi_P);

    // Combined cross-entropy
    return H_len_cross + H_ss_cross + h_cross * (E_len_P - (1.0 - P0_P));
}

// ---------------------------------------------------------------------------
// Event metadata structure
// ---------------------------------------------------------------------------

/**
 * @brief Metadata for one recombination event.
 *
 * Stores the event's name, nickname, number of realizations,
 * its position in the model queue, and the theoretical marginal distribution.
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
static inline std::vector<EventInfo> build_event_info(
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

// ---------------------------------------------------------------------------
// Combined insertion + dinucleotide entropy
// ---------------------------------------------------------------------------

/**
 * @brief Pairing of insertion event with its DinucMarkov partner
 */
struct InsDinucPair {
    const EventInfo *ins_event = nullptr;
    const EventInfo *dinuc_event = nullptr;
    double combined_H = 0.0;
};

/**
 * @brief Compute combined insertion+dinucleotide entropy and update EventInfo.H
 * 
 * Following pygor3's get_df_entropy_decomposition(), insertion events
 * like "vj_ins" are paired with their DinucMarkov partner "vj_dinucl"
 * (sharing the same Gene_class). The combined entropy H_total is computed
 * using insertion_dinuc_entropy() and the EventInfo.H for insertion events
 * is updated to reflect the full combined entropy.
 * 
 * @param events Vector of EventInfo (modified in-place)
 * @param parms Model parameters (to extract insertion length realizations)
 * @param print_decomposition Whether to print entropy decomposition table
 * @return Map of Gene_class to InsDinucPair for reference
 */
static inline std::map<Gene_class, InsDinucPair> 
compute_combined_insertion_dinuc_entropy(
        std::vector<EventInfo> &events,
        const Model_Parms &parms,
        bool print_decomposition = true)
{
    std::map<Gene_class, InsDinucPair> ins_dinuc_pairs;

    // 1. Identify insertion and DinucMarkov event pairs by Gene_class
    for (const auto &ev : events) {
        if (ev.is_dinuc_markov) {
            ins_dinuc_pairs[ev.gene_class].dinuc_event = &ev;
        } else if (ev.nickname.find("_ins") != std::string::npos) {
            ins_dinuc_pairs[ev.gene_class].ins_event = &ev;
        }
    }

    // 2. Compute combined entropy for each pair
    for (auto &[gc, pair] : ins_dinuc_pairs) {
        if (pair.ins_event && pair.dinuc_event) {
            // Get insertion length values from the event's realizations
            auto ev_ptr = parms.get_event_pointer(pair.ins_event->name);
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

    // 3. Update EventInfo.H for insertion events to the full combined entropy
    for (auto &ev : events) {
        if (ev.is_dinuc_markov) continue;
        auto it = ins_dinuc_pairs.find(ev.gene_class);
        if (it != ins_dinuc_pairs.end() &&
            it->second.ins_event &&
            it->second.ins_event->queue_position == ev.queue_position &&
            it->second.dinuc_event)
        {
            ev.H = it->second.combined_H;
        }
    }

    // 4. Print entropy decomposition (mirrors pygor3's entropy table)
    if (print_decomposition) {
        double total_entropy = 0.0;
        std::cout << "\n=== Entropy decomposition (bits) ===" << std::endl;
        for (const auto &ev : events) {
            if (ev.is_dinuc_markov) {
                // Already accounted for via the ins+dinuc pair
                continue;
            }

            // Check if this is an insertion event with a DinucMarkov partner
            auto it = ins_dinuc_pairs.find(ev.gene_class);
            if (it != ins_dinuc_pairs.end() &&
                it->second.ins_event == &ev &&
                it->second.dinuc_event != nullptr)
            {
                double combined = it->second.combined_H;
                std::cout << "  " << ev.nickname << " + "
                          << it->second.dinuc_event->nickname
                          << " : H = " << combined
                          << "  (H_len=" << entropy(ev.model_marginal)
                          << ", H_dinuc=" << (combined - entropy(ev.model_marginal))
                          << ", h=" << it->second.dinuc_event->dinuc_entropy_rate
                          << ")" << std::endl;
                total_entropy += combined;
            } else {
                std::cout << "  " << ev.nickname << " : H = " << ev.H << std::endl;
                total_entropy += ev.H;
            }
        }
        // Also print standalone DinucMarkov entropy rates for reference
        for (const auto &ev : events) {
            if (ev.is_dinuc_markov) {
                std::cout << "  (" << ev.nickname << " entropy rate h = "
                          << ev.dinuc_entropy_rate << " bits/nt)" << std::endl;
            }
        }
        std::cout << "  TOTAL: " << total_entropy << std::endl;
    }

    return ins_dinuc_pairs;
}
