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
#include <igor/Core/FastGenerator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <forward_list>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

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
 * @brief Pair of insertion and dinucleotide Markov events for combined entropy.
 */
struct InsDinucPair {
    const EventInfo *ins_event = nullptr;
    const EventInfo *dinuc_event = nullptr;
    double combined_H = 0.0;
};

// ---------------------------------------------------------------------------
// Empirical marginals computation
// ---------------------------------------------------------------------------

/**
 * @brief Compute empirical marginal distributions for ALL non-DinucMarkov
 *        events in a single pass over the scenario list.
 *
 * Returns a map from event queue_position to its empirical marginal vector.
 * This avoids the O(events × N) queue-copying that caused timeouts
 * when compute_empirical_marginal was called per-event.
 */
static inline std::map<size_t, std::vector<double>> compute_all_empirical_marginals(
        const std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> &scenarios,
        const std::vector<EventInfo> &event_infos,
        size_t total_sequences)
{
    // Prepare per-event count arrays (only for non-DinucMarkov events)
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
                    throw std::out_of_range(
                            "Realization index " + std::to_string(realization) +
                            " out of bounds [0, " + std::to_string(num_real) + ")");
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
static inline std::map<size_t, std::vector<double>> compute_all_empirical_marginals(
        const std::vector<igor::fast::GeneratedSequence> &sequences,
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

    // Single pass over all generated sequences
    for (const auto &seq : sequences) {
        for (auto &[qpos, cvec] : counts) {
            if (qpos < seq.realizations.size() && !seq.realizations[qpos].empty()) {
                int realization = seq.realizations[qpos][0];
                int num_real = static_cast<int>(cvec.size());
                if (realization < 0 || realization >= num_real) {
                    throw std::out_of_range(
                            "Realization index " + std::to_string(realization) +
                            " out of bounds [0, " + std::to_string(num_real) + ")");
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

/**
 * @brief Results of comparing two distributions for table display.
 * 
 * Generic structure used by both test_generation (theoretical vs empirical)
 * and test_inference (ground truth vs inferred).
 */
struct ComparisonRow {
    std::string event_nickname;
    double kl_divergence = 0.0;
    double H_reference = 0.0;      ///< H_truth or H_theoretical
    double H_compared = 0.0;       ///< H_inferred or H_empirical
    double uncovered_mass = 0.0;
    bool passes = false;
    
    // Sampling baseline (for inference tests)
    double kl_sampling_baseline = -1.0;  ///< D_KL(truth || empirical), -1 if not computed
    
    // For insertion+dinuc pairs: decomposition into components
    bool is_insertion_dinuc_pair = false;
    double kl_combined = 0.0;
    double H_combined_reference = 0.0;
    double H_combined_compared = 0.0;
    double kl_length = 0.0;
    double H_length_reference = 0.0;
    double H_length_compared = 0.0;
    double uncovered_length = 0.0;
    bool passes_length = false;
    double kl_length_sampling_baseline = -1.0;  ///< Sampling baseline for length distribution
    double kl_combined_sampling_baseline = -1.0; ///< Sampling baseline for combined (if available)
    double kl_dinuc = 0.0;
    double h_dinuc_reference = 0.0;
    double h_dinuc_compared = 0.0;
    bool passes_combined = false;
    bool passes_dinuc = false;
};

/**
 * @brief Build comparison rows from event info and compared marginals.
 * 
 * Handles both inference (truth vs inferred) and generation (model vs empirical) cases.
 * For insertion+dinuc pairs, shows decomposition if compute_combined_kl is true.
 * 
 * @param event_infos List of events with theoretical marginals
 * @param compared_marginals Map from queue_position to compared marginal distribution
 * @param ins_dinuc_pairs Map of insertion+dinuc pairs with combined entropy
 * @param compared_parms Optional Model_Parms for computing combined D_KL
 * @param compared_model Optional Model_marginals for computing combined D_KL
 * @param compute_combined_kl If true, compute combined/dinuc D_KL (requires compared_model)
 * @param kl_threshold_func Function to compute pass threshold from entropy
 * @param sampling_baseline_kl Optional map of sampling baseline D_KL values (for inference tests)
 */
static inline std::vector<ComparisonRow> build_comparison_rows(
        const std::vector<EventInfo> &event_infos,
        const std::map<size_t, std::vector<double>> &compared_marginals,
        const std::map<Gene_class, InsDinucPair> &ins_dinuc_pairs,
        const Model_Parms *compared_parms = nullptr,
        const Model_marginals *compared_model = nullptr,
        bool compute_combined_kl = false,
        std::function<bool(double, double, int)> kl_threshold_func = nullptr,
        const std::map<size_t, double> *sampling_baseline_kl = nullptr)
{
    std::vector<ComparisonRow> rows;
    
    // Get index map if we need to compute combined D_KL
    std::unordered_map<std::string, int> compared_index_map;
    if (compute_combined_kl && compared_model && compared_parms) {
        compared_index_map = compared_model->get_index_map(*compared_parms);
    }
    
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) {
            continue; // Skip - shown as part of insertion pairs
        }
        
        const auto &compared = compared_marginals.at(ev.queue_position);
        
        // Compute D_KL and uncovered mass for marginal distribution
        double uncovered = 0.0;
        double dkl = kl_divergence(ev.model_marginal, compared, &uncovered);
        double H_compared = entropy(compared);
        
        // Check if this is part of an insertion+dinuc pair
        auto it_pair = ins_dinuc_pairs.find(ev.gene_class);
        bool is_ins_event = it_pair != ins_dinuc_pairs.end() &&
                            it_pair->second.ins_event &&
                            it_pair->second.ins_event->queue_position == ev.queue_position &&
                            it_pair->second.dinuc_event;
        
        ComparisonRow row;
        row.event_nickname = ev.nickname;
        
        // Populate sampling baseline if provided
        if (sampling_baseline_kl && sampling_baseline_kl->count(ev.queue_position) > 0) {
            row.kl_sampling_baseline = sampling_baseline_kl->at(ev.queue_position);
        }
        
        if (is_ins_event) {
            // Insertion+dinuc pair: show decomposition
            const auto &pair = it_pair->second;
            row.is_insertion_dinuc_pair = true;
            
            // Combined row
            row.H_combined_reference = pair.combined_H;
            row.passes_combined = true;
            
            // Length component (always available)
            row.kl_length = dkl;
            row.H_length_reference = entropy(ev.model_marginal);
            row.H_length_compared = H_compared;
            row.uncovered_length = uncovered;
            row.kl_length_sampling_baseline = row.kl_sampling_baseline; // Same baseline for length dist
            
            // Dinuc component (theoretical)
            row.h_dinuc_reference = pair.dinuc_event->dinuc_entropy_rate;
            row.passes_dinuc = true;
            
            // Compute combined/dinuc D_KL if requested
            if (compute_combined_kl && compared_model && compared_parms) {
                const auto *dinuc_ev = pair.dinuc_event;
                
                // Get compared dinuc transition matrix
                int dinuc_base_idx = compared_index_map.at(dinuc_ev->name);
                std::array<double, 16> compared_dinuc_T;
                for (int k = 0; k < 16; ++k) {
                    compared_dinuc_T[k] = static_cast<double>(
                            compared_model->marginal_array_smart_p[dinuc_base_idx + k]);
                }
                
                // Get insertion lengths
                auto ev_ptr = compared_parms->get_event_pointer(ev.name);
                auto realizations = ev_ptr->get_realizations_map();
                std::vector<int> ins_lengths(ev.num_realizations, 0);
                for (const auto &[key, real] : realizations) {
                    ins_lengths[real.index] = real.value_int;
                }
                
                // Combined entropy for compared model
                row.H_combined_compared = insertion_dinuc_entropy(
                        compared, ins_lengths, compared_dinuc_T);
                
                // Combined D_KL
                double combined_cross_entropy = insertion_dinuc_cross_entropy(
                        ev.model_marginal, compared,
                        ins_lengths, dinuc_ev->dinuc_T, compared_dinuc_T);
                row.kl_combined = combined_cross_entropy - pair.combined_H;
                
                // Dinuc D_KL
                row.kl_dinuc = markov_kl_divergence(dinuc_ev->dinuc_T, compared_dinuc_T);
                row.h_dinuc_compared = markov_entropy_rate(compared_dinuc_T);
            } else {
                // Generation test: can't compute combined/dinuc D_KL
                row.H_combined_compared = 0.0;
                row.kl_combined = 0.0;
                row.kl_dinuc = 0.0;
                row.h_dinuc_compared = 0.0;
            }
            
            // Apply threshold function if provided
            if (kl_threshold_func) {
                row.passes_length = kl_threshold_func(row.kl_length, row.H_length_reference, ev.num_realizations);
                if (compute_combined_kl) {
                    row.passes_combined = kl_threshold_func(row.kl_combined, row.H_combined_reference, -1);
                    row.passes_dinuc = kl_threshold_func(row.kl_dinuc, row.h_dinuc_reference, -1);
                }
            } else {
                row.passes_length = true;
            }
        } else {
            // Regular event
            row.kl_divergence = dkl;
            row.H_reference = ev.H;
            row.H_compared = H_compared;
            row.uncovered_mass = uncovered;
            
            // Apply threshold function if provided
            if (kl_threshold_func) {
                row.passes = kl_threshold_func(dkl, ev.H, ev.num_realizations);
            } else {
                row.passes = true;
            }
        }
        
        rows.push_back(row);
    }
    
    return rows;
}

/**
 * @brief Print a formatted table of comparison results.
 * 
 * @param rows Vector of comparison rows to display
 * @param reference_label Label for reference column (e.g. "H_truth" or "H_model")
 * @param compared_label Label for compared column (e.g. "H_infer" or "H_emp")
 * @param show_sampling_baseline If true, show D_KL(sampling) column for inference tests
 */
static inline void print_comparison_table(
        const std::vector<ComparisonRow> &rows,
        const std::string &reference_label = "H_ref",
        const std::string &compared_label = "H_cmp",
        bool show_sampling_baseline = false)
{
    // Check if any row has sampling baseline
    bool has_baseline = show_sampling_baseline && std::any_of(rows.begin(), rows.end(),
            [](const ComparisonRow& r) { return r.kl_sampling_baseline >= 0.0; });
    
    // Print header
    std::cout << "\nEvent                  | D_KL(R||C) ";
    if (has_baseline) {
        std::cout << "| D_KL(samp) ";
    }
    std::cout << "| " << std::setw(7) << reference_label << " | " 
              << std::setw(7) << compared_label << " | Uncovered | Pass" << std::endl;
    
    std::cout << "---------------------- | ---------- ";
    if (has_baseline) {
        std::cout << "| ---------- ";
    }
    std::cout << "| ------- | ------- | --------- | ----" << std::endl;
    
    for (const auto& row : rows) {
        if (row.is_insertion_dinuc_pair) {
            // Main row (empty metrics)
            std::cout << std::left << std::setw(22) << row.event_nickname << " | "
                      << std::fixed 
                      << "          " << " | ";
            if (has_baseline) {
                std::cout << "          " << " | ";
            }
            std::cout << "       " << " | "
                      << "       " << " | "
                      << "         " << " | "
                      << " " << std::endl;
            // Combined row
            std::cout << "  ├─ combined           | "
                      << std::setprecision(4) << std::setw(10);
            if (row.kl_combined > 0.0) {
                std::cout << row.kl_combined;
            } else {
                std::cout << "-";
            }
            std::cout << " | ";
            if (has_baseline) {
                std::cout << std::setw(10);
                if (row.kl_combined_sampling_baseline >= 0.0) {
                    std::cout << row.kl_combined_sampling_baseline;
                } else {
                    std::cout << "-";
                }
                std::cout << " | ";
            }
            std::cout << std::setw(7) << row.H_combined_reference << " | " << std::setw(7);
            if (row.H_combined_compared > 0.0) {
                std::cout << row.H_combined_compared;
            } else {
                std::cout << "-";
            }
            std::cout << " | "
                      << "         " << " | "
                      << (row.passes_combined ? "✓" : "✗") << std::endl;
            // Length component
            std::cout << "  ├─ ins length         | "
                      << std::setw(10) << row.kl_length << " | ";
            if (has_baseline) {
                std::cout << std::setw(10);
                if (row.kl_length_sampling_baseline >= 0.0) {
                    std::cout << row.kl_length_sampling_baseline;
                } else {
                    std::cout << "-";
                }
                std::cout << " | ";
            }
            std::cout << std::setw(7) << row.H_length_reference << " | "
                      << std::setw(7) << row.H_length_compared << " | "
                      << std::setw(9) << row.uncovered_length << " | "
                      << (row.passes_length ? "✓" : "✗") << std::endl;
            // Dinuc component
            std::cout << "  └─ dinuc Markov (h)   | "
                      << std::setw(10);
            if (row.kl_dinuc > 0.0) {
                std::cout << row.kl_dinuc;
            } else {
                std::cout << "-";
            }
            std::cout << " | ";
            if (has_baseline) {
                std::cout << std::setw(10) << "-" << " | ";  // No sampling baseline for dinuc
            }
            std::cout << std::setw(7) << row.h_dinuc_reference << " | " << std::setw(7);
            if (row.h_dinuc_compared > 0.0) {
                std::cout << row.h_dinuc_compared;
            } else {
                std::cout << "-";
            }
            std::cout << " | "
                      << "         " << " | "
                      << (row.passes_dinuc ? "✓" : "✗") << std::endl;
        } else {
            // Regular event row
            std::cout << std::left << std::setw(22) << row.event_nickname << " | "
                      << std::fixed << std::setprecision(4) << std::setw(10) << row.kl_divergence << " | ";
            if (has_baseline) {
                std::cout << std::setw(10);
                if (row.kl_sampling_baseline >= 0.0) {
                    std::cout << row.kl_sampling_baseline;
                } else {
                    std::cout << "-";
                }
                std::cout << " | ";
            }
            std::cout << std::setw(7) << row.H_reference << " | "
                      << std::setw(7) << row.H_compared << " | "
                      << std::setw(9) << row.uncovered_mass << " | "
                      << (row.passes ? "✓" : "✗") << std::endl;
        }
    }
}

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
        info.gene_class = ev->get_class();
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

// ---------------------------------------------------------------------------
// Combined insertion + dinucleotide entropy
// ---------------------------------------------------------------------------

/**
 * @brief Pairing of insertion event with its DinucMarkov partner
 */
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
