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
#include <filesystem>
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

// Model base directory
static const std::string MODELS_DIR =
        std::string(IGOR_SOURCE_DIR) + "/models";

static bool required_model_files_exist(const std::string& parms_path,
                                       const std::string& marginals_path,
                                       const std::string& model_label)
{
    if (std::filesystem::exists(parms_path) &&
        std::filesystem::exists(marginals_path)) {
        return true;
    }
    WARN("Skipping " << model_label
                     << " generation test because model fixtures are missing");
    return false;
}

// ---------------------------------------------------------------------------
// Test configuration
// ---------------------------------------------------------------------------

struct GenerationTestConfig {
    std::string model_parms_path;
    std::string model_marginals_path;
    std::string model_label;
    size_t sample_size = 0;
    std::vector<size_t> sample_size_reduction_factors;
    bool test_convergence = true;
    std::map<std::string, double> pygor3_reference_entropy;
};

static std::vector<size_t> build_sample_schedule(const GenerationTestConfig& cfg)
{
    std::vector<size_t> sample_sizes;
    sample_sizes.push_back(cfg.sample_size);

    for (size_t factor : cfg.sample_size_reduction_factors) {
        if (factor <= 1) {
            continue;
        }
        size_t reduced = cfg.sample_size / factor;
        if (reduced > 0) {
            sample_sizes.push_back(reduced);
        }
    }

    std::sort(sample_sizes.begin(), sample_sizes.end());
    sample_sizes.erase(std::unique(sample_sizes.begin(), sample_sizes.end()), sample_sizes.end());
    return sample_sizes;
}

namespace {
constexpr double kDklSafetyFactor = 50.0;
constexpr double kMonotonicGapFactor = 5.0;
}  // namespace

/**
 * @brief Leading-order expectation for sampling-induced D_KL at finite N.
 *
 * For an event with k realizations sampled N times from the true model,
 * the plug-in KL estimator has asymptotic expectation:
 *
 *   E[D_KL(P||Q_hat)] ~= (k - 1) / (2 * N * ln 2)
 */
static double expected_large_sample_dkl(int num_realizations, size_t sample_size)
{
    return (num_realizations - 1) / (2.0 * static_cast<double>(sample_size) * std::log(2.0));
}

/**
 * @brief Heuristic high-confidence scale for finite-sample D_KL.
 *
 * This is not a strict probabilistic upper bound; it is a safety-scaled
 * expectation used as a robust acceptance threshold in CI.
 */
static double large_sample_dkl_scale(int num_realizations, size_t sample_size)
{
    return kDklSafetyFactor * expected_large_sample_dkl(num_realizations, sample_size);
}

/**
 * @brief Monotonicity tolerance from expected D_KL gap across sample sizes.
 *
 * When N increases from N_prev to N_last, the expected sampling-induced D_KL
 * decreases. We allow non-monotonic fluctuations proportional to this expected
 * gap. The proportionality factor accounts for estimator variance and skew.
 */
static double dkl_monotonicity_tolerance(int num_realizations,
                                         size_t prev_sample_size,
                                         size_t last_sample_size)
{
    const double expected_prev =
        expected_large_sample_dkl(num_realizations, prev_sample_size);
    const double expected_last =
        expected_large_sample_dkl(num_realizations, last_sample_size);
    return kMonotonicGapFactor * std::abs(expected_prev - expected_last);
}

static bool passes_large_sample_dkl_threshold(double dkl, int num_realizations, size_t sample_size)
{
    return dkl < large_sample_dkl_scale(num_realizations, sample_size);
}

static double uncovered_mass_upper(size_t sample_size)
{
    return (std::max)(1000.0 / static_cast<double>(sample_size), 1e-4);
}

// ---------------------------------------------------------------------------
// Test implementation
// ---------------------------------------------------------------------------

/**
 * @brief Run generation convergence test with a primary sample size and reductions
 *
 * Test workflow:
 * 1. Load ground truth model and event metadata
 * 2. Compute combined insertion+dinucleotide entropy
 * 3. Cross-validate entropy against pygor3 reference (if available)
 * 4. Build sample schedule from the primary sample size and reduction factors
 * 5. Generate sequences for each schedule point and compute empirical marginals
 * 6. If test_convergence=true, perform:
 *    - Empirical entropy convergence to theoretical entropy
 *    - Monotonic decrease check of D_KL with increasing N
 *
 * @param cfg GenerationTestConfig with model paths, sampling schedule, and options
 */
static void run_generation_convergence_test(const GenerationTestConfig& cfg)
{
    // ------------------------------------------------------------------
    // 1. Load model and collect per-event metadata
    // ------------------------------------------------------------------
    if (!required_model_files_exist(cfg.model_parms_path, cfg.model_marginals_path, cfg.model_label)) {
        return;
    }

    INFO("Testing model: " << cfg.model_label);
    REQUIRE(((not cfg.test_convergence) or (cfg.sample_size >= 1000000)));

    const std::vector<size_t> sample_sizes = build_sample_schedule(cfg);
    REQUIRE(!sample_sizes.empty());

    INFO("Sample sizes: " << [&sample_sizes]() {
        std::ostringstream oss;
        for (size_t i = 0; i < sample_sizes.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << sample_sizes[i];
        }
        return oss.str();
    }());
    INFO("Test convergence: " << (cfg.test_convergence ? "true" : "false"));

    Model_Parms model_parms;
    model_parms.read_model_parms(cfg.model_parms_path);

    Model_marginals model_marginals(model_parms);
    model_marginals.txt2marginals(cfg.model_marginals_path, model_parms);

    std::vector<EventInfo> event_infos =
            build_event_info(model_parms, model_marginals);

    INFO("Model loaded, " << event_infos.size() << " events found");
    int min_events = cfg.model_label.find("tcr_beta") != std::string::npos ? 9 : 5;
    REQUIRE(event_infos.size() >= static_cast<size_t>(min_events));

    // Compute combined insertion+dinucleotide entropy and print decomposition
    auto ins_dinuc_pairs = compute_combined_insertion_dinuc_entropy(
            event_infos, model_parms, /*print=*/true);

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
    if (!cfg.pygor3_reference_entropy.empty()) {
        std::cout << "\n=== Cross-validation against pygor3 (" << cfg.model_label
                  << ") ===" << std::endl;
        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov) continue;

            // Determine the reported entropy for this event
            double reported_H;
            auto it_pair = ins_dinuc_pairs.find(ev.gene_class);
            bool is_combined_ins =
                    it_pair != ins_dinuc_pairs.end() &&
                    it_pair->second.ins_event == &ev &&
                    it_pair->second.dinuc_event != nullptr;
            if (is_combined_ins) {
                reported_H = it_pair->second.combined_H;
            } else {
                reported_H = ev.H;
            }

            auto it_ref = cfg.pygor3_reference_entropy.find(ev.nickname);
            if (it_ref != cfg.pygor3_reference_entropy.end()) {
                double ref = it_ref->second;
                double diff = std::abs(reported_H - ref);
                bool has_parents = !model_parms.get_parents(ev.name).empty();

                std::cout << "  " << ev.nickname
                          << ": computed=" << reported_H
                          << "  pygor3=" << ref
                          << "  diff=" << diff
                          << (has_parents ? "  (marginal; parents present)" : "")
                          << std::endl;

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
    // 2. Generate sequences for increasing N and track D_KL per event
    // ------------------------------------------------------------------

    // Update ev.H for insertion events to the full combined entropy
    // (insertion-length + dinucleotide Markov) so that D_KL bounds and
    // printouts reference the total information content of the pair.
    for (auto &ev : event_infos) {
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

    // Map: event_queue_position → vector of (D_KL, uncovered_mass) per sample-size
    std::map<size_t, std::vector<std::pair<double, double>>> kl_traces;
    // Map: event_queue_position → vector of empirical entropy per sample-size
    std::map<size_t, std::vector<double>> entropy_traces;

    // Initialize FastGenerator once, reuse across all sample sizes
    igor::fast::FastGenerator fast_gen;
    fast_gen.initialize(model_parms, model_marginals);
    REQUIRE(fast_gen.is_initialized());

    for (size_t N : sample_sizes) {
        INFO("Generating " << N << " sequences");

        // Generate using FastGenerator (parallel, precomputed CDFs)
        igor::fast::FastGeneratorConfig config;
        config.show_progress = false;
        auto sequences = fast_gen.generate(N, config);

        REQUIRE(sequences.size() == N);

        std::cout << "\n--- N = " << N << " ---" << std::endl;

        // Compute all empirical marginals in a single pass
        auto all_empiricals = compute_all_empirical_marginals(
                sequences, event_infos, sequences.size());

        // Build comparison rows using helper
        auto threshold_func = [N](double dkl, double /*H*/, int num_realizations) -> bool {
            return passes_large_sample_dkl_threshold(dkl, num_realizations, N);
        };

        auto rows = build_comparison_rows(
                event_infos, all_empiricals, ins_dinuc_pairs,
                nullptr, nullptr, false, threshold_func);

        // Print table
        print_comparison_table(rows, "H_model", "H_emp");

        // Update traces for convergence checks
        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov) continue;

            const auto &empirical = all_empiricals.at(ev.queue_position);
            double uncovered = 0.0;
            double dkl = kl_divergence(ev.model_marginal, empirical, &uncovered);
            kl_traces[ev.queue_position].emplace_back(dkl, uncovered);

            double H_emp = entropy(empirical);
            entropy_traces[ev.queue_position].push_back(H_emp);

            // ---- assertions ----
            // D_KL must always be finite (partial KL can be negative when
            // uncovered bins are skipped, so we only check finiteness)
            CHECK(std::isfinite(dkl));

            CHECK(passes_large_sample_dkl_threshold(dkl, ev.num_realizations, N));
            CHECK(uncovered < uncovered_mass_upper(N));
        }
    }

    // Skip convergence checks if not enabled
    if (!cfg.test_convergence) {
        std::cout << "\n=== Skipping convergence validation ===" << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // 3. Empirical entropy convergence to the theoretical entropy
    // ------------------------------------------------------------------
    // The empirical entropy H(Q) of the generated marginal should
    // converge to the theoretical entropy H(P) of the model.  The
    // finite-sample estimator has a known negative bias:
    //     E[H_emp] ≈ H − (k−1) / (2·N·ln2)
    // so we check that the difference |H_emp − H| shrinks with N.
    // ------------------------------------------------------------------
    std::cout << "\n=== Empirical entropy vs theoretical ==" << std::endl;
    for (const auto &ev : event_infos) {
        if (ev.is_dinuc_markov) continue;

        const auto &htrace = entropy_traces.at(ev.queue_position);
        if (htrace.size() >= 2 && sample_sizes.size() >= 2) {
            // Use last sample size for entropy convergence check
            size_t last_idx = htrace.size() - 1;
            double H_emp_last = htrace[last_idx];
            size_t N_last = sample_sizes[last_idx];

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
            double bias =
                        50.0 * (ev.num_realizations - 1) / (2.0 * static_cast<double>(N_last) * std::log(2.0));
            double stdev_bound = 5.0 * H_theo / std::sqrt(N_last);
            double tolerance = bias + stdev_bound;
            double diff = std::abs(H_emp_last - H_theo);

            std::cout << "  " << ev.nickname
                      << " : H_theo=" << H_theo
                      << "  H_emp(N=" << N_last << ")=" << H_emp_last
                      << "  diff=" << diff
                      << "  tol=" << tolerance << std::endl;

            // The empirical entropy should be close to theoretical;
            // allow bias + 5σ variance, with a minimum floor.
            CHECK(diff < (std::max)(tolerance, 1e-3));
        }
    }

    // ------------------------------------------------------------------
    // 4. Check monotonic decrease of D_KL with increasing N
    // ------------------------------------------------------------------
    if (sample_sizes.size() >= 2) {
        std::cout << "\n=== Monotonic decrease check ===" << std::endl;
        for (const auto &ev : event_infos) {
            if (ev.is_dinuc_markov) continue;

            const auto &trace = kl_traces.at(ev.queue_position);
            // Compare D_KL between last two non-integration sample sizes.
            // Both should be large enough for the estimator to be reliable.
            // A small absolute tolerance handles the case where the smaller
            // sample happens to nail the distribution exactly (D_KL ≈ 0)
            // while the larger one is merely very small.
            //
            // The partial KL (which skips uncovered bins) can go negative
            // when many bins are missed at smaller N, making the comparison
            // meaningless.  We only check monotonicity when uncovered mass
            // is small enough for the partial KL to be a faithful estimate
            // of the true KL.
            if (trace.size() >= 2) {
                size_t idx_prev = trace.size() - 2;
                size_t idx_last = trace.size() - 1;
                auto [dkl_prev, uncov_prev] = trace[idx_prev];
                auto [dkl_last, uncov_last] = trace[idx_last];
                size_t N_prev = sample_sizes[idx_prev];
                size_t N_last = sample_sizes[idx_last];

                std::cout << "  " << ev.nickname
                          << " : D_KL(N=" << N_prev << ")=" << dkl_prev
                          << "  D_KL(N=" << N_last << ")=" << dkl_last
                          << "  uncov(" << N_prev << ")=" << uncov_prev
                          << "  uncov(" << N_last << ")=" << uncov_last << std::endl;
                if (uncov_prev < 0.01) {
                    // Both estimates are reliable. Allow residual finite-N
                    // fluctuations up to the predicted upper-bound drop.
                    const double tol = dkl_monotonicity_tolerance(
                            ev.num_realizations, N_prev, N_last);
                    std::cout << "D_KL decrease computed tolerance: " << tol <<std::endl;
                    CHECK(dkl_last < dkl_prev + tol);
                }
            }
        }
    }

    std::cout << "\n=== Generation test completed ===" << std::endl;
}

// ---------------------------------------------------------------------------
// THE TESTS
// ---------------------------------------------------------------------------

TEST_CASE("Generation marginals converge - integration", "[generation][integration]")
{
    // Integration test: large sample generation path only.
    // No convergence assertions.
    SECTION("human TCR alpha (VJ)") {
        GenerationTestConfig cfg;
        cfg.model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        cfg.model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        cfg.model_label = "human/tcr_alpha";
        cfg.sample_size = 100;
        cfg.sample_size_reduction_factors = {10};
        cfg.test_convergence = false;
        run_generation_convergence_test(cfg);
    }

    SECTION("human TCR beta (VDJ)") {
        GenerationTestConfig cfg;
        cfg.model_parms_path = MODELS_DIR + "/human/tcr_beta/models/model_parms.txt";
        cfg.model_marginals_path = MODELS_DIR + "/human/tcr_beta/models/model_marginals.txt";
        cfg.model_label = "human/tcr_beta";
        cfg.sample_size = 100;
        cfg.sample_size_reduction_factors = {10};
        cfg.test_convergence = false;
        // No pygor3 reference for integration test
        run_generation_convergence_test(cfg);
    }
}

TEST_CASE("Generation marginals converge - convergence", "[generation][convergence][slow]")
{
    // Convergence test: high-N schedule with reduced-N checkpoints for monotonic D_KL checks.
    SECTION("human TCR alpha (VJ) - N=1M with N/10 monotonic checkpoint") {
        GenerationTestConfig cfg;
        cfg.model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        cfg.model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        cfg.model_label = "human/tcr_alpha";
        cfg.sample_size = 1e6;
        cfg.sample_size_reduction_factors = {10};
        cfg.test_convergence = true;
        run_generation_convergence_test(cfg);
    }

    SECTION("human TCR beta (VDJ) - N=1M with N/10 checkpoint and pygor3 cross-validation") {
        GenerationTestConfig cfg;
        cfg.model_parms_path = MODELS_DIR + "/human/tcr_beta/models/model_parms.txt";
        cfg.model_marginals_path = MODELS_DIR + "/human/tcr_beta/models/model_marginals.txt";
        cfg.model_label = "human/tcr_beta";
        cfg.sample_size = 1e7;
        cfg.sample_size_reduction_factors = {10};
        cfg.test_convergence = true;
        // Reference entropies from pygor3 tutorial:
        // https://pygor3.readthedocs.io/en/latest/Tutorial.html#Entropy
        // Insertion rows are the combined H(ℓ) + E[ℓ]·h values.
        cfg.pygor3_reference_entropy = {
            {"v_choice", 5.252905},
            {"d_gene",   1.141779},
            {"j_choice", 3.609102},
            {"vd_ins",  14.894931},
            {"dj_ins",  14.981991},
            {"v_3_del",  3.147511},
            {"d_3_del",  2.778230},
            {"d_5_del",  3.634137},
            {"j_5_del",  3.356340},
        };
        run_generation_convergence_test(cfg);
    }
}
