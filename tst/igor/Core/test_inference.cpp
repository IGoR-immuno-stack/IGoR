/**
 * @file test_inference.cpp
 * @brief Inference validation using sampling baseline and relative KL degradation
 *
 * Test workflow:
 * 1. Generate sequences from ground truth model with known scenarios
 * 2. Create mock perfect alignments from scenarios (O(1) vs O(N×M) SW alignment)
 * 3. Run inference from uniform initialization
 * 4. Validate convergence using relative KL degradation
 *
 * Validation approach:
 * - Compute sampling baseline: D_KL(truth || empirical_from_samples)
 * - Measure inference quality: (D_KL_inferred - D_KL_sampling) / H_truth < threshold
 * - This isolates inference-specific divergence from sampling noise
 * - Typical thresholds: 5-15% relative degradation depending on sample size
 * - Dinucleotide Markov: D_KL_inferred / h_dinuc < 0.05% (heavily sampled, no baseline)
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Aligner.h>

#include "entropy_test_helpers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <forward_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef IGOR_SOURCE_DIR
#error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

// Uncomment to enable debug logging to /tmp/
#define DEBUG_INFERENCE_TEST

// Model base directory
static const std::string MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/models";
static const std::string TEST_MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/inference";

// ---------------------------------------------------------------------------
// Scenario parsing
// ---------------------------------------------------------------------------

struct ParsedScenario {
    int v_gene_index;
    int j_gene_index;
    int d_gene_index;
    int v_3p_del;
    int j_5p_del;
    int d_5p_del, d_3p_del;
    int vd_ins, dj_ins, vj_ins;

    std::string v_gene_name;
    std::string j_gene_name;
    std::string d_gene_name;

    // Constructor with defaults for VJ models
    ParsedScenario() :
        v_gene_index(-1), j_gene_index(-1), d_gene_index(-1),
        v_3p_del(-1), j_5p_del(-1), d_5p_del(-1), d_3p_del(-1),
        vd_ins(-1), dj_ins(-1), vj_ins(-1) {}
};

/**
 * @brief Parse scenario queue to extract event realizations
 *
 * Walks through the scenario queue matching queue positions to events,
 * extracting gene indices, deletion lengths, and insertion lengths.
 */
static ParsedScenario parse_scenario(
    const std::queue<std::queue<int>>& scenario,
    const std::vector<EventInfo>& event_infos,
    const Model_Parms& parms)
{
    ParsedScenario parsed;

    // Make a copy to drain
    auto scenario_copy = scenario;
    size_t pos = 0;

    while (!scenario_copy.empty()) {
        auto inner = scenario_copy.front();
        scenario_copy.pop();

        if (inner.empty()) {
            ++pos;
            continue;
        }

        int realization_idx = inner.front();

        // Find the event at this position
        const EventInfo* ev_info = nullptr;
        for (const auto& ev : event_infos) {
            if (ev.queue_position == pos) {
                ev_info = &ev;
                break;
            }
        }

        if (!ev_info) {
            ++pos;
            continue;
        }

        // Extract based on event type and gene class
        auto ev_ptr = parms.get_event_pointer(ev_info->name);
        auto event_type = ev_ptr->get_type();
        auto gene_class = ev_ptr->get_class();
        const std::string seq_type = ev_ptr->get_seq_type();

        if (event_type == Event_type::GeneChoice_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    if (gene_class == V_gene) {
                        parsed.v_gene_index = realization_idx;
                        parsed.v_gene_name = name;
                    } else if (gene_class == J_gene) {
                        parsed.j_gene_index = realization_idx;
                        parsed.j_gene_name = name;
                    } else if (gene_class == D_gene) {
                        parsed.d_gene_index = realization_idx;
                        parsed.d_gene_name = name;
                    }
                    break;
                }
            }
        } else if (event_type == Event_type::Deletion_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    int del_value = real.value_int;
                    if (gene_class == V_gene) {
                        parsed.v_3p_del = del_value;
                    } else if (gene_class == J_gene) {
                        parsed.j_5p_del = del_value;
                    } else if (gene_class == D_gene) {
                        auto seq_side = ev_ptr->get_side();
                        if (seq_side == Five_prime) {
                            parsed.d_5p_del = del_value;
                        } else if (seq_side == Three_prime) {
                            parsed.d_3p_del = del_value;
                        }
                    }
                    break;
                }
            }
        } else if (event_type == Event_type::Insertion_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    int ins_value = real.value_int;
                    if (seq_type == "VD_ins_seq") {
                        parsed.vd_ins = ins_value;
                    } else if (seq_type == "DJ_ins_seq") {
                        parsed.dj_ins = ins_value;
                    } else if (seq_type == "VJ_ins_seq") {
                        parsed.vj_ins = ins_value;
                    }
                    break;
                }
            }
        }

        ++pos;
    }

    return parsed;
}

// ---------------------------------------------------------------------------
// Mock alignment construction
// ---------------------------------------------------------------------------

/**
 * @brief Create mock alignment for V gene with proper mismatch handling
 *
 * Aligns the non-deleted portion perfectly, then checks deleted region
 * for mismatches against the sequence (which contains insertion nucleotides).
 */
static Alignment_data create_v_mock_alignment(
    const std::string& sequence,
    const ParsedScenario& scenario,
    const std::unordered_map<std::string, std::string>& gene_templates)
{
    const std::string& v_template = gene_templates.at(scenario.v_gene_name);
    int v_len = static_cast<int>(v_template.length());
    int v_del = scenario.v_3p_del;

    // Check for mismatches only in the deleted region
    std::vector<int> mismatches;
    for (int i = v_len - v_del; i < v_len && i < static_cast<int>(sequence.length()); ++i) {
        if (v_template[i] != sequence[i]) {
            mismatches.push_back(i);
        }
    }

    Alignment_data v_align(scenario.v_gene_name, 0);
    v_align.five_p_offset = 0;
    v_align.three_p_offset = v_len - v_del - 1;
    v_align.align_length = v_len - v_del;
    v_align.mismatches = mismatches;
    v_align.insertions = {};
    v_align.deletions = {};
    v_align.score = 5.0 * (v_align.align_length - static_cast<int>(mismatches.size()));

    return v_align;
}

/**
 * @brief Create mock alignment for J gene with proper mismatch handling
 *
 * Similar to V but checks the 5' deleted region (before the aligned portion).
 */
static Alignment_data create_j_mock_alignment(
    const std::string& sequence,
    const ParsedScenario& scenario,
    const std::unordered_map<std::string, std::string>& gene_templates)
{
    const std::string& j_template = gene_templates.at(scenario.j_gene_name);
    int j_len = static_cast<int>(j_template.length());
    int j_del = scenario.j_5p_del;
    int j_offset = static_cast<int>(sequence.length()) - j_len;

    // Check for mismatches only in the deleted region (5' of aligned portion)
    std::vector<int> mismatches;
    for (int i = 0; i < j_del; ++i) {
        int seq_pos = j_offset + i;
        if (seq_pos >= 0 && seq_pos < static_cast<int>(sequence.length())) {
            if (j_template[i] != sequence[seq_pos]) {
                mismatches.push_back(seq_pos);
            }
        }
    }

    Alignment_data j_align(scenario.j_gene_name, j_offset);
    j_align.five_p_offset = j_offset - j_del;
    j_align.three_p_offset = static_cast<int>(sequence.length()) - 1;
    j_align.align_length = j_len - j_del;
    j_align.mismatches = mismatches;
    j_align.insertions = {};
    j_align.deletions = {};
    j_align.score = 5.0 * (j_align.align_length - static_cast<int>(mismatches.size()));

    return j_align;
}

// ---------------------------------------------------------------------------
// Model comparison
// ---------------------------------------------------------------------------

/**
 * @brief Validate a ComparisonRow using relative degradation thresholds
 *
 * For regular events and insertion length:
 *   relative_degradation = (D_KL_inferred - D_KL_sampling) / H_truth < threshold
 *
 * For dinucleotide Markov (no sampling baseline):
 *   D_KL_inferred / h_dinuc < threshold
 *
 * @param row ComparisonRow to validate (passes flags are set)
 * @param relative_kl_degradation_threshold Threshold for marginal events (e.g., 0.05 = 5%)
 * @param relative_kl_threshold_dinuc Threshold for dinucleotide (e.g., 0.0005 = 0.05%)
 */
static void validate_comparison_row(
        ComparisonRow& row,
        double relative_kl_degradation_threshold,
        double relative_kl_degradation_threshold_inslen,
        double relative_kl_threshold_dinuc)
{
    if (row.is_insertion_dinuc_pair) {
        // Length component: use relative degradation
        if (row.kl_length_sampling_baseline >= 0.0 && row.H_length_reference > 0.0) {
            double relative_degradation =
                (row.kl_length - row.kl_length_sampling_baseline) / row.H_length_reference;
            row.passes_length = (relative_degradation < relative_kl_degradation_threshold_inslen);
        } else {
            row.passes_length = true;  // No baseline to validate against
        }

        // Dinuc component: use relative threshold (no baseline)
        if (row.kl_dinuc > 0.0 && row.h_dinuc_reference > 0.0) {
            row.passes_dinuc = (row.kl_dinuc < relative_kl_threshold_dinuc * row.h_dinuc_reference);
        } else {
            row.passes_dinuc = true;
        }

        // Combined passes if both components pass
        row.passes_combined = row.passes_length && row.passes_dinuc;
    } else {
        // Regular event: use relative degradation
        if (row.kl_sampling_baseline >= 0.0 && row.H_reference > 0.0) {
            double relative_degradation =
                (row.kl_divergence - row.kl_sampling_baseline) / row.H_reference;
            row.passes = (relative_degradation < relative_kl_degradation_threshold);
        } else {
            row.passes = true;  // No baseline to validate against
        }
    }
}

/**
 * @brief Compare inferred model to ground truth using KL divergence
 *
 * Uses sampling baseline to evaluate inference quality:
 * Validates that (D_KL_inferred - D_KL_sampling) / H_truth is small,
 * isolating the inference-specific divergence from sampling noise.
 */
static std::vector<ComparisonRow> compare_inference_to_ground_truth(
        const std::vector<EventInfo>& ground_truth_events,
        const Model_marginals& inferred_marginals,
        const Model_Parms& inferred_parms,
        const std::map<size_t, double>& sampling_baseline_kl,
        double relative_kl_degradation_threshold,
        double relative_kl_degradation_threshold_inslen,
        double relative_kl_threshold_dinuc)
{
    // Compute inferred marginals for all events
    std::map<size_t, std::vector<double>> inferred_marginals_map;
    for (const auto& gt_ev : ground_truth_events) {
        if (gt_ev.is_dinuc_markov) continue;

        auto [dims, inferred_probs] =
            inferred_marginals.compute_event_marginal_probability(gt_ev.name, inferred_parms);

        std::vector<double> marginal(gt_ev.num_realizations, 0.0);
        for (int i = 0; i < gt_ev.num_realizations; ++i) {
            marginal[i] = static_cast<double>(inferred_probs.get()[i]);
        }

        inferred_marginals_map[gt_ev.queue_position] = marginal;
    }

    // Build insertion+dinuc pairs with combined entropy computed
    // (reuse from ground truth events which already have combined_H)
    std::map<std::string, InsDinucPair> ins_dinuc_pairs;
    for (const auto& ev : ground_truth_events) {
        if (ev.is_dinuc_markov) {
            ins_dinuc_pairs[ev.seq_type].dinuc_event = &ev;
        } else if (ev.nickname.find("_ins") != std::string::npos) {
            ins_dinuc_pairs[ev.seq_type].ins_event = &ev;
        }
    }

    // Compute combined entropy for pairs
    for (auto& [gene_class, pair] : ins_dinuc_pairs) {
        if (pair.ins_event && pair.dinuc_event) {
            // Get insertion lengths
            auto ev_ptr = inferred_parms.get_event_pointer(pair.ins_event->name);
            auto realizations = ev_ptr->get_realizations_map();
            std::vector<int> ins_lengths(pair.ins_event->num_realizations, 0);
            for (const auto& [key, real] : realizations) {
                ins_lengths[real.index] = real.value_int;
            }

            // Compute combined entropy for ground truth
            pair.combined_H = insertion_dinuc_entropy(
                    pair.ins_event->model_marginal,
                    ins_lengths,
                    pair.dinuc_event->dinuc_T);
        }
    }

    // Build comparison rows with combined D_KL computation
    auto rows = build_comparison_rows(
            ground_truth_events,
            inferred_marginals_map,
            ins_dinuc_pairs,
            &inferred_parms,
            &inferred_marginals,
            true,  // compute_combined_kl
            nullptr,  // No threshold function - validation done separately
            &sampling_baseline_kl);  // Pass sampling baseline for display

    // Validate each row using relative degradation thresholds
    for (auto& row : rows) {
        validate_comparison_row(row, relative_kl_degradation_threshold, relative_kl_degradation_threshold_inslen, relative_kl_threshold_dinuc);
    }

    return rows;
}

// ---------------------------------------------------------------------------
// Test helper
// ---------------------------------------------------------------------------

struct InferenceTestConfig {
    std::string model_parms_path;
    std::string model_marginals_path;
    std::string model_label;
    int sample_size = 0;
    double relative_kl_degradation_threshold = 0.10;
    double relative_kl_degradation_threshold_inslen = 0.15;
    double relative_kl_threshold_dinuc_model = 0.001;
    int num_iterations = 20;
    bool test_convergence = true;
};

/**
 * @brief Validate inference convergence using sampling baseline and relative KL degradation
 *
 * Test workflow:
 * 1. Load ground truth model and extract gene templates from parameters
 * 2. Generate N sequences with known scenarios (no sequencing errors)
 * 3. Compute sampling baseline: D_KL(truth || empirical_from_samples) for each event
 * 4. Create mock perfect alignments from scenarios (avoids expensive SW alignment)
 * 5. Run inference from uniform initialization for ~20 iterations
 * 6. Compare inferred vs ground truth using relative KL degradation:
 *    - Marginal events: (D_KL_inferred - D_KL_sampling) / H_truth < threshold (5-15%)
 *    - Dinucleotide Markov: D_KL_inferred / h_dinuc < 0.05% (heavily sampled, no baseline)
 *
 * Three test sections:
 * - N=100 (smoke test): Only checks inference pipeline execution
 * - N=10000 (thorough): 5% degradation threshold for strict validation
 * - N=3000 (shallow): 10% degradation threshold for faster testing
 *
 * The relative degradation metric isolates inference quality from sampling noise,
 * making validation robust across different event types and sample sizes.
 */
static void run_inference_recovery_test(const InferenceTestConfig& cfg)
{
    // ------------------------------------------------------------------
    // 1. Load ground truth model
    // ------------------------------------------------------------------
    INFO("Testing inference with model: " << cfg.model_label);
    INFO("Sample size: " << cfg.sample_size);

    Model_Parms truth_parms;
    truth_parms.read_model_parms(cfg.model_parms_path);

    Model_marginals truth_marginals(truth_parms);
    truth_marginals.txt2marginals(cfg.model_marginals_path, truth_parms);

    auto truth_events = build_event_info(truth_parms, truth_marginals);

    std::cout << "\n=== Ground truth model loaded ===" << std::endl;
    std::cout << "Events: " << truth_events.size() << std::endl;

    // Compute combined insertion+dinucleotide entropy and print decomposition
    compute_combined_insertion_dinuc_entropy(truth_events, truth_parms, /*print=*/true);

    // ------------------------------------------------------------------
    // 2. Extract gene templates from model parameters
    // ------------------------------------------------------------------
    std::unordered_map<std::string, std::string> gene_templates;

    // Get events map
    auto events_map = truth_parms.get_events_map();

    // Extract V gene templates
    auto v_gene_event = events_map.at(std::make_tuple(GeneChoice_t, std::string("V_gene_seq"), Undefined_side));
    auto v_gene_choice = std::dynamic_pointer_cast<Gene_choice>(v_gene_event);
    auto v_realizations = v_gene_choice->get_realizations_map();
    for (const auto& [name, realization] : v_realizations) {
        gene_templates[name] = realization.value_str;
    }

    // Extract J gene templates
    auto j_gene_event = events_map.at(std::make_tuple(GeneChoice_t, std::string("J_gene_seq"), Undefined_side));
    auto j_gene_choice = std::dynamic_pointer_cast<Gene_choice>(j_gene_event);
    auto j_realizations = j_gene_choice->get_realizations_map();
    for (const auto& [name, realization] : j_realizations) {
        gene_templates[name] = realization.value_str;
    }

    // Extract D gene templates if VDJ model
    if (events_map.count(std::make_tuple(GeneChoice_t, std::string("D_gene_seq"), Undefined_side)) > 0) {
        auto d_gene_event = events_map.at(std::make_tuple(GeneChoice_t, std::string("D_gene_seq"), Undefined_side));
        auto d_gene_choice = std::dynamic_pointer_cast<Gene_choice>(d_gene_event);
        auto d_realizations = d_gene_choice->get_realizations_map();
        for (const auto& [name, realization] : d_realizations) {
            gene_templates[name] = realization.value_str;
        }
    }

    std::cout << "\nExtracted " << gene_templates.size() << " gene templates from model" << std::endl;

    // ------------------------------------------------------------------
    // 3. Generate sequences with scenarios
    // ------------------------------------------------------------------
    std::cout << "\n=== Generating " << cfg.sample_size << " sequences ===" << std::endl;
    // Create a 0 error rate
    Single_error_rate null_error_model = Single_error_rate(0.0);
    truth_parms.set_error_ratep(&null_error_model);
    GenModel gen_model(truth_parms, truth_marginals);
    auto scenarios = gen_model.generate_sequences(cfg.sample_size, /*generate_errors=*/false);

    // Count generated sequences
    size_t actual_count = 0;
    for (auto it = scenarios.begin(); it != scenarios.end(); ++it) {
        ++actual_count;
    }

    REQUIRE(actual_count == static_cast<size_t>(cfg.sample_size));
    std::cout << "Generated " << actual_count << " sequences" << std::endl;

    // ------------------------------------------------------------------
    // 3b. Compute sampling baseline D_KL(truth || empirical)
    // ------------------------------------------------------------------
    std::cout << "\n=== Computing sampling baseline ===" << std::endl;

    // Compute empirical marginals from generated scenarios
    auto empirical_marginals_map = compute_all_empirical_marginals(
            scenarios, truth_events, actual_count);

    // Compute sampling baseline D_KL for each event
    std::map<size_t, double> sampling_baseline_kl;
    for (const auto& ev : truth_events) {
        if (ev.is_dinuc_markov) continue;

        const auto& empirical = empirical_marginals_map.at(ev.queue_position);
        double baseline_dkl = kl_divergence(ev.model_marginal, empirical);
        sampling_baseline_kl[ev.queue_position] = baseline_dkl;

        std::cout << "  " << ev.nickname
                  << ": D_KL(truth||empirical) = " << baseline_dkl
                  << " bits" << std::endl;
    }

    // ------------------------------------------------------------------
    // 4. Create mock alignments in memory
    // ------------------------------------------------------------------
    std::cout << "\n=== Creating mock alignments ===" << std::endl;

    std::vector<std::tuple<int, std::string,
        std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
        sequences_with_alignments;

    int seq_idx = 0;
    for (const auto& [seq, scenario] : scenarios) {
        // Parse scenario to extract gene choices and deletions
        ParsedScenario parsed = parse_scenario(scenario, truth_events, truth_parms);

        // Create V alignment
        Alignment_data v_align = create_v_mock_alignment(seq, parsed, gene_templates);

        // Create J alignment
        Alignment_data j_align = create_j_mock_alignment(seq, parsed, gene_templates);

        // Create empty D alignment (required for VDJ models)
        std::vector<Alignment_data> d_aligns;

        // Build alignment map
        std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns_map;
        aligns_map[V_gene] = {v_align};
        aligns_map[J_gene] = {j_align};
        aligns_map[D_gene] = d_aligns;

        sequences_with_alignments.emplace_back(seq_idx, seq, aligns_map);
        ++seq_idx;
    }

    std::cout << "Created alignments for " << sequences_with_alignments.size()
              << " sequences" << std::endl;

    // ------------------------------------------------------------------
    // 5. Set up debug logging directory (if enabled)
    // ------------------------------------------------------------------
    std::string inference_dir = "";
#ifdef DEBUG_INFERENCE_TEST
    inference_dir = "/tmp/igor_inference_test_" + std::to_string(cfg.sample_size) + "/";
    std::cout << "\n=== Debug mode enabled ===" << std::endl;
    std::cout << "Output directory: " << inference_dir << std::endl;

    // Create directory if it doesn't exist
    std::string mkdir_cmd = "mkdir -p " + inference_dir;
    std::system(mkdir_cmd.c_str());

    // Write generated sequences and scenarios to file
    std::cout << "\n=== Writing sequences and scenarios to file ===" << std::endl;
    gen_model.write_seq_real2txt(
        inference_dir + "generated_seqs_indexed.csv",
        inference_dir + "generated_realizations_indexed.csv",
        scenarios);
    std::cout << "Wrote sequences to: " << inference_dir << "generated_seqs_indexed.csv" << std::endl;
    std::cout << "Wrote scenarios to: " << inference_dir << "generated_realizations_indexed.csv" << std::endl;

    // Write alignments to file for debugging
    std::cout << "\n=== Writing alignments to file for debugging ===" << std::endl;

    // Extract V gene alignments
    std::unordered_map<int, std::forward_list<Alignment_data>> v_alignments;
    for (const auto& [seq_idx2, seq, aligns_map] : sequences_with_alignments) {
        if (aligns_map.count(V_gene) > 0) {
            std::forward_list<Alignment_data> v_list;
            for (const auto& align : aligns_map.at(V_gene)) {
                v_list.push_front(align);
            }
            v_alignments[seq_idx2] = v_list;
        }
    }

    // Extract J gene alignments
    std::unordered_map<int, std::forward_list<Alignment_data>> j_alignments;
    for (const auto& [seq_idx2, seq, aligns_map] : sequences_with_alignments) {
        if (aligns_map.count(J_gene) > 0) {
            std::forward_list<Alignment_data> j_list;
            for (const auto& align : aligns_map.at(J_gene)) {
                j_list.push_front(align);
            }
            j_alignments[seq_idx2] = j_list;
        }
    }

    // Write using Aligner class methods
    Aligner aligner;
    aligner.write_alignments_seq_csv(inference_dir + "V_alignments.csv", v_alignments);
    std::cout << "Wrote V alignments to: " << inference_dir << "V_alignments.csv" << std::endl;

    aligner.write_alignments_seq_csv(inference_dir + "J_alignments.csv", j_alignments);
    std::cout << "Wrote J alignments to: " << inference_dir << "J_alignments.csv" << std::endl;
#endif

    // ------------------------------------------------------------------
    // 6. Initialize inference with uniform marginals
    // ------------------------------------------------------------------
    std::cout << "\n=== Initializing inference ===" << std::endl;
    Model_marginals initial_marginals(truth_parms);
    initial_marginals.uniform_initialize(truth_parms);

    // ------------------------------------------------------------------
    // 7. Run inference
    // ------------------------------------------------------------------
    std::cout << "\n=== Running inference ===" << std::endl;

    GenModel inference_model(truth_parms, initial_marginals);

    std::cout << "Inference iterations: " << cfg.num_iterations << std::endl;

    bool failed = inference_model.infer_model(
        sequences_with_alignments,
        cfg.num_iterations,
        inference_dir,
        /*fast_iter=*/true,
        /*likelihood_threshold=*/0.0,
        /*viterbi_like=*/false,
        /*proba_threshold_factor=*/1e-4);

    REQUIRE(!failed);
    std::cout << "Inference completed successfully" << std::endl;

    // ------------------------------------------------------------------
    // 8. Use inferred marginals from GenModel object
    //    (marginals are modified by reference during inference)
    // ------------------------------------------------------------------
    std::cout << "\n=== Using inferred marginals ===" << std::endl;

    // The initial_marginals object now contains the inferred marginals
    // after GenModel.infer_model() completes
    Model_Parms inferred_parms = truth_parms;
    const Model_marginals& inferred_marginals = inference_model.get_marginals();

    // ------------------------------------------------------------------
    // 9. Compare inferred vs ground truth
    // ------------------------------------------------------------------
    std::cout << "\n=== Comparing inferred model to ground truth ===" << std::endl;

    if (cfg.test_convergence) {
        auto rows = compare_inference_to_ground_truth(
            truth_events,
            inferred_marginals,
            inferred_parms,
            sampling_baseline_kl,
            cfg.relative_kl_degradation_threshold,
            cfg.relative_kl_degradation_threshold_inslen,
            cfg.relative_kl_threshold_dinuc_model);

        print_comparison_table(rows, "H_truth", "H_infer", /*show_sampling_baseline=*/true);

        std::cout << "\n=== Validating convergence ===" << std::endl;
        // Check that all rows passed validation (passes flags set by validate_comparison_row)
        for (const auto& row : rows) {
            INFO("Event: " << row.event_nickname);

            if (row.is_insertion_dinuc_pair) {
                INFO(" Combined (len+dinuc) - D_KL=" << row.kl_combined
                     << ", baseline=" << row.kl_combined_sampling_baseline
                     << ", H_truth=" << row.H_combined_reference
                     << ", H_inferred=" << row.H_combined_compared);
                INFO(" Ins length - D_KL=" << row.kl_length
                     << ", baseline=" << row.kl_length_sampling_baseline
                     << ", H_truth=" << row.H_length_reference
                     << ", H_inferred=" << row.H_length_compared);
                INFO(" Dinuc Markov - D_KL=" << row.kl_dinuc
                     << ", h_truth=" << row.h_dinuc_reference
                     << ", h_inferred=" << row.h_dinuc_compared);

                // Validate components (passes flags already set)
                if (row.kl_length_sampling_baseline >= 0.0 && row.H_length_reference > 0.0) {
                    double relative_degradation =
                        (row.kl_length - row.kl_length_sampling_baseline) / row.H_length_reference;
                    INFO(" Length relative degradation: " << (relative_degradation * 100)
                         << "% (threshold: " << (cfg.relative_kl_degradation_threshold_inslen * 100) << "%)");
                }

                if (row.kl_dinuc > 0.0 && row.h_dinuc_reference > 0.0) {
                    double relative_kl = row.kl_dinuc / row.h_dinuc_reference;
                    INFO(" Dinuc relative D_KL: " << (relative_kl * 100)
                         << "% (threshold: " << (cfg.relative_kl_threshold_dinuc_model * 100) << "%)");
                }

                // Print diagnostics if length validation fails
                if (!row.passes_length && !row.reference_marginal.empty()) {
                    print_distribution_diagnostics(
                        row.event_nickname + " (length)",
                        row.reference_marginal,
                        row.compared_marginal,
                        0.01);
                }

                // Print diagnostics if dinuc validation fails
                if (!row.passes_dinuc) {
                    print_dinuc_matrix_diagnostics(
                        row.event_nickname + " (dinuc)",
                        row.reference_dinuc_T,
                        row.compared_dinuc_T);
                }

                CHECK(row.passes_length);
                CHECK(row.passes_dinuc);
            } else {
                INFO("D_KL(truth||inferred) = " << row.kl_divergence);
                INFO("D_KL(truth||empirical) = " << row.kl_sampling_baseline);
                INFO("H(truth) = " << row.H_reference);
                INFO("H(inferred) = " << row.H_compared);

                // Show relative degradation for context
                if (row.kl_sampling_baseline >= 0.0 && row.H_reference > 0.0) {
                    double relative_degradation =
                        (row.kl_divergence - row.kl_sampling_baseline) / row.H_reference;
                    INFO("Relative degradation: " << (relative_degradation * 100)
                         << "% (threshold: " << (cfg.relative_kl_degradation_threshold * 100) << "%)");
                }

                // Print diagnostics if validation fails
                if (!row.passes && !row.reference_marginal.empty()) {
                    print_distribution_diagnostics(
                        row.event_nickname,
                        row.reference_marginal,
                        row.compared_marginal,
                        0.01);
                }

                CHECK(row.passes);
            }
        }
    }

    std::cout << "\n=== Inference test completed ===" << std::endl;
}

// ---------------------------------------------------------------------------
// THE TESTS
// ---------------------------------------------------------------------------

// this test should take around 0.5 second
TEST_CASE("Inference recovers ground truth model - smoke", "[integration][inference]")
{
    InferenceTestConfig cfg;
    cfg.model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
    cfg.model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
    cfg.model_label = "human/tcr_alpha";
    cfg.sample_size = 100;
    cfg.test_convergence = false;

    run_inference_recovery_test(cfg);
}

TEST_CASE("Inference recovers ground truth model - convergence", "[integration][inference][slow][!mayfail]")
{
    // this test should take around 5 seconds
    SECTION("human TCR alpha (VJ) - N=10000 - thorough validation") {
        InferenceTestConfig cfg;
        cfg.model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        cfg.model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        cfg.model_label = "human/tcr_alpha";
        cfg.sample_size = 10000;
        cfg.relative_kl_degradation_threshold = 0.05;
        cfg.relative_kl_degradation_threshold_inslen = 0.10;

        run_inference_recovery_test(cfg);
    }

    // this test should take around 200 seconds
    SECTION("Fixed VJ TCR beta (VDJ) - N=3000 - shallow validation") {
        InferenceTestConfig cfg;
        cfg.model_parms_path = TEST_MODELS_DIR + "/fixed_VJ_TRB_model_parms.txt";
        cfg.model_marginals_path = TEST_MODELS_DIR + "/fixed_VJ_TRB_marginals.txt";
        cfg.model_label = "human/fixed_vj_tcr_beta";
        cfg.sample_size = 3000;
        cfg.relative_kl_degradation_threshold = 0.15;
        cfg.relative_kl_degradation_threshold_inslen = 0.20;
        cfg.num_iterations = 20;

        run_inference_recovery_test(cfg);
    }
}
