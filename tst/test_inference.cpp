/**
 * @file test_inference.cpp
 * @brief Tests the inference machinery using KL divergence validation
 *
 * Generates sequences from a known ground truth model, creates mock alignments
 * from the known scenarios, runs inference, and validates that the inferred
 * model parameters match the ground truth using KL divergence and entropy.
 *
 * The mock alignments simulate perfect SW alignment without the O(N×M) cost:
 * - Non-deleted regions align perfectly
 * - Deleted regions are compared against sequence to find mismatches
 * - IGoR's Gene_choice::iterate correctly filters deletable mismatches
 *
 * KL divergence thresholds:
 *   • N = 1000:  D_KL < H / 5   (fast smoke test)
 *   • N = 10000: D_KL < H / 20  (thorough validation)
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
        
        if (event_type == Event_type::GeneChoice_t) {
            // Get gene name from realization
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
                    if (gene_class == VD_genes) {
                        parsed.vd_ins = ins_value;
                    } else if (gene_class == DJ_genes) {
                        parsed.dj_ins = ins_value;
                    } else if (gene_class == VJ_genes) {
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
    int v_len = v_template.length();
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
    int j_len = j_template.length();
    int j_del = scenario.j_5p_del;
    int j_offset = sequence.length() - j_len;
    
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
    j_align.three_p_offset = sequence.length() - 1;
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
 * @brief Compare inferred model to ground truth using KL divergence
 */
static std::vector<ComparisonRow> compare_inference_to_ground_truth(
        const std::vector<EventInfo>& ground_truth_events,
        const Model_marginals& inferred_marginals,
        const Model_Parms& inferred_parms,
        double kl_threshold_factor = 10.0)
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
    std::map<Gene_class, InsDinucPair> ins_dinuc_pairs;
    for (const auto& ev : ground_truth_events) {
        if (ev.is_dinuc_markov) {
            ins_dinuc_pairs[ev.gene_class].dinuc_event = &ev;
        } else if (ev.nickname.find("_ins") != std::string::npos) {
            ins_dinuc_pairs[ev.gene_class].ins_event = &ev;
        }
    }
    
    // Compute combined entropy for pairs
    for (auto& [gene_class, pair] : ins_dinuc_pairs) {
        if (pair.ins_event && pair.dinuc_event) {
            // Get insertion lengths
            auto ev_ptr = inferred_parms.get_event_pointer(pair.ins_event->name);
            auto realizations = ev_ptr->get_realizations_map();
            std::vector<int> ins_lengths(pair.ins_event->num_realizations, 0);
            for (const auto &[key, real] : realizations) {
                ins_lengths[real.index] = real.value_int;
            }
            
            // Compute combined entropy for ground truth
            pair.combined_H = insertion_dinuc_entropy(
                    pair.ins_event->model_marginal, 
                    ins_lengths, 
                    pair.dinuc_event->dinuc_T);
        }
    }
    
    // Threshold function for inference
    auto threshold_func = [kl_threshold_factor](double dkl, double H, int) -> bool {
        double threshold = (std::max)(H / kl_threshold_factor, 0.01);
        return dkl < threshold;
    };
    
    // Use helper to build comparison rows with combined D_KL computation
    return build_comparison_rows(
            ground_truth_events,
            inferred_marginals_map,
            ins_dinuc_pairs,
            &inferred_parms,
            &inferred_marginals,
            true,  // compute_combined_kl
            threshold_func);
}

// ---------------------------------------------------------------------------
// THE TEST
// ---------------------------------------------------------------------------

TEST_CASE("Inference recovers ground truth model", "[inference]")
{
    std::string model_parms_path;
    std::string model_marginals_path;
    std::string model_label;
    int sample_size = 0;
    double kl_threshold_factor = 5.0;
    
    SECTION("human TCR alpha (VJ) - N=1000 - smoke test") {
        model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        model_label = "human/tcr_alpha";
        sample_size = 1000;
        kl_threshold_factor = 5.0;  // Loose threshold for fast test
    }
    
    SECTION("human TCR alpha (VJ) - N=10000 - thorough validation") {
        model_parms_path = MODELS_DIR + "/human/tcr_alpha/models/model_parms.txt";
        model_marginals_path = MODELS_DIR + "/human/tcr_alpha/models/model_marginals.txt";
        model_label = "human/tcr_alpha";
        sample_size = 10000;
        kl_threshold_factor = 20.0;  // Stricter threshold
    }

    SECTION("Fixed VJ TCR beta (VDJ) - N=10000 - thorough validation") {
        model_parms_path = TEST_MODELS_DIR + "/fixed_VJ_TRB_model_parms.txt";
        model_marginals_path = TEST_MODELS_DIR + "/fixed_VJ_TRB_marginals.txt";
        model_label = "human/fixed_vj_tcr_beta";
        sample_size = 10000;
        kl_threshold_factor = 20.0;  // Stricter threshold
    }
    
    // ------------------------------------------------------------------
    // 1. Load ground truth model
    // ------------------------------------------------------------------
    INFO("Testing inference with model: " << model_label);
    INFO("Sample size: " << sample_size);
    
    Model_Parms truth_parms;
    truth_parms.read_model_parms(model_parms_path);
    
    Model_marginals truth_marginals(truth_parms);
    truth_marginals.txt2marginals(model_marginals_path, truth_parms);
    
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
    auto v_gene_event = events_map.at(std::make_tuple(GeneChoice_t, V_gene, Undefined_side));
    auto v_gene_choice = std::dynamic_pointer_cast<Gene_choice>(v_gene_event);
    auto v_realizations = v_gene_choice->get_realizations_map();
    for (const auto& [name, realization] : v_realizations) {
        gene_templates[name] = realization.value_str;
    }
    
    // Extract J gene templates
    auto j_gene_event = events_map.at(std::make_tuple(GeneChoice_t, J_gene, Undefined_side));
    auto j_gene_choice = std::dynamic_pointer_cast<Gene_choice>(j_gene_event);
    auto j_realizations = j_gene_choice->get_realizations_map();
    for (const auto& [name, realization] : j_realizations) {
        gene_templates[name] = realization.value_str;
    }
    
    // Extract D gene templates if VDJ model
    if (events_map.count(std::make_tuple(GeneChoice_t, D_gene, Undefined_side)) > 0) {
        auto d_gene_event = events_map.at(std::make_tuple(GeneChoice_t, D_gene, Undefined_side));
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
    std::cout << "\n=== Generating " << sample_size << " sequences ===" << std::endl;
    GenModel gen_model(truth_parms, truth_marginals);
    auto scenarios = gen_model.generate_sequences(sample_size, /*generate_errors=*/false);
    
    // Count generated sequences
    size_t actual_count = 0;
    for (auto it = scenarios.begin(); it != scenarios.end(); ++it) {
        ++actual_count;
    }
    REQUIRE(actual_count == static_cast<size_t>(sample_size));
    std::cout << "Generated " << actual_count << " sequences" << std::endl;
    
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
        
        // Build alignment map
        std::unordered_map<Gene_class, std::vector<Alignment_data>> aligns_map;
        aligns_map[V_gene] = {v_align};
        aligns_map[J_gene] = {j_align};
        
        sequences_with_alignments.emplace_back(seq_idx, seq, aligns_map);
        ++seq_idx;
    }
    
    std::cout << "Created alignments for " << sequences_with_alignments.size() 
              << " sequences" << std::endl;
    
    // ------------------------------------------------------------------
    // 5. Set up debug logging directory (if enabled)
    // ------------------------------------------------------------------
    std::string inference_dir = "";  // Disable file logging by default
    #ifdef DEBUG_INFERENCE_TEST
    inference_dir = "/tmp/igor_inference_test_" + std::to_string(sample_size) + "/";
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
    for (const auto& [seq_idx, seq, aligns_map] : sequences_with_alignments) {
        if (aligns_map.count(V_gene) > 0) {
            std::forward_list<Alignment_data> v_list;
            for (const auto& align : aligns_map.at(V_gene)) {
                v_list.push_front(align);
            }
            v_alignments[seq_idx] = v_list;
        }
    }
    
    // Extract J gene alignments
    std::unordered_map<int, std::forward_list<Alignment_data>> j_alignments;
    for (const auto& [seq_idx, seq, aligns_map] : sequences_with_alignments) {
        if (aligns_map.count(J_gene) > 0) {
            std::forward_list<Alignment_data> j_list;
            for (const auto& align : aligns_map.at(J_gene)) {
                j_list.push_front(align);
            }
            j_alignments[seq_idx] = j_list;
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
    
    int num_iterations = 20;
    std::cout << "Inference iterations: " << num_iterations << std::endl;
    
    bool failed = inference_model.infer_model(
        sequences_with_alignments,
        num_iterations,
        inference_dir,
        /*fast_iter=*/true,
        /*likelihood_threshold=*/1e-60,
        /*viterbi_like=*/false,
        /*proba_threshold_factor=*/1e-3);
    
    REQUIRE(!failed);
    std::cout << "Inference completed successfully" << std::endl;
    
    // ------------------------------------------------------------------
    // 8. Use inferred marginals from GenModel object
    //    (marginals are modified by reference during inference)
    // ------------------------------------------------------------------
    std::cout << "\n=== Using inferred marginals ===" << std::endl;
    
    // The initial_marginals object now contains the inferred marginals
    // after GenModel.infer_model() completes
    Model_Parms inferred_parms = truth_parms;  // Structure is the same
    const Model_marginals& inferred_marginals = inference_model.get_marginals();
    
    // ------------------------------------------------------------------
    // 9. Compare inferred vs ground truth
    // ------------------------------------------------------------------
    std::cout << "\n=== Comparing inferred model to ground truth ===" << std::endl;
    
    auto rows = compare_inference_to_ground_truth(
        truth_events,
        inferred_marginals,
        inferred_parms,
        kl_threshold_factor);
    
    print_comparison_table(rows, "H_truth", "H_infer");
    
    std::cout << "\n=== Inference test completed ===" << std::endl;
    
    // Check assertions
    for (const auto& row : rows) {
        INFO("Event: " << row.event_nickname);
        
        if (row.is_insertion_dinuc_pair) {
            INFO("  Combined (len+dinuc) - D_KL=" << row.kl_combined
                 << ", H_truth=" << row.H_combined_reference
                 << ", H_inferred=" << row.H_combined_compared);
            INFO("  Ins length - D_KL=" << row.kl_length 
                 << ", H_truth=" << row.H_length_reference 
                 << ", H_inferred=" << row.H_length_compared);
            INFO("  Dinuc Markov - D_KL=" << row.kl_dinuc
                 << ", h_truth=" << row.h_dinuc_reference
                 << ", h_inferred=" << row.h_dinuc_compared);
            
            CHECK(row.passes_combined);
            CHECK(row.passes_length);
            CHECK(row.passes_dinuc);
        } else {
            INFO("D_KL(truth||inferred) = " << row.kl_divergence);
            INFO("H(truth) = " << row.H_reference);
            INFO("H(inferred) = " << row.H_compared);
            CHECK(row.passes);
        }
    }
    
    std::cout << "\n=== Inference test completed ===" << std::endl;
}
