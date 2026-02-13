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
static inline double kl_divergence(const std::vector<double> &P, const std::vector<double> &Q,
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
static inline std::array<double, 4> markov_stationary_distribution(const std::array<double, 16> &T, int max_iter = 1000,
                                                                   double tol = 1e-12)
{
    // Start with uniform distribution
    std::array<double, 4> pi = { 0.25, 0.25, 0.25, 0.25 };
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
        for (double v : pi_next)
            sum += v;
        for (double &v : pi_next)
            v /= sum;

        // Check convergence
        double max_diff = 0.0;
        for (int i = 0; i < 4; ++i) {
            max_diff = (std::max)(max_diff, std::abs(pi_next[i] - pi[i]));
        }
        pi = pi_next;
        if (max_diff < tol)
            break;
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
static inline double markov_entropy_rate(const std::array<double, 16> &T)
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
static inline double insertion_dinuc_entropy(const std::vector<double> &ins_marginal,
                                             const std::vector<int> &ins_lengths, const std::array<double, 16> &dinuc_T)
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
            if (ins_lengths[i] == 0)
                P0 += p;
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

// ---------------------------------------------------------------------------
// Event metadata structure
// ---------------------------------------------------------------------------

/**
 * @brief Metadata for one recombination event.
 *
 * Stores the event's name, nickname, number of realizations,
 * its position in the model queue, and the theoretical marginal distribution.
 */
struct EventInfo
{
    Rec_Event_name name;
    std::string nickname;
    int num_realizations;
    size_t queue_position; ///< index in the model-queue order
    bool is_dinuc_markov; ///< skip empirical check for DinucMarkov
    std::vector<double> model_marginal; ///< P(realization) marginalised over parents
    double H; ///< entropy of the marginal
    Gene_class gene_class; ///< gene class (VD_genes, DJ_genes, VJ_genes …)
    std::array<double, 16> dinuc_T; ///< transition matrix (only for DinucMarkov)
    double dinuc_entropy_rate; ///< Markov entropy rate h (only for DinucMarkov)
};

/**
 * @brief Build a vector of EventInfo for every event in the model queue.
 */
static inline std::vector<EventInfo> build_event_info(const Model_Parms &parms, const Model_marginals &marginals)
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
        info.gene_class = static_cast<Gene_class>(ev->get_class());
        info.dinuc_T.fill(0.0);
        info.dinuc_entropy_rate = 0.0;

        std::cout << "  Processing event: " << info.nickname << " (is_dinuc=" << info.is_dinuc_markov
                  << ", size=" << info.num_realizations << ")" << std::endl;

        if (info.is_dinuc_markov) {
            // Extract the 4×4 transition matrix from the raw marginal array.
            // DinucMarkov has no parents so 16 values sit contiguously
            // at index_map[name].
            int base_idx = index_map.at(info.name);
            for (int k = 0; k < 16; ++k) {
                info.dinuc_T[k] = static_cast<double>(marginals.marginal_array_smart_p[base_idx + k]);
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
            auto [dims, probs] = marginals.compute_event_marginal_probability(info.name, parms);

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
