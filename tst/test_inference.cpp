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
// #define DEBUG_INFERENCE_TEST

// Model base directory
static const std::string MODELS_DIR = std::string(IGOR_SOURCE_DIR) + "/models";

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

struct InferenceComparison {
    std::string event_nickname;
    double kl_divergence_forward;   // D_KL(truth || inferred) - penalizes missing support
    double kl_divergence_reverse;   // D_KL(inferred || truth) - penalizes impossible events
    double entropy_truth;
    double entropy_inferred;
    double uncovered_mass;
    bool passes_threshold;
    
    // For insertion+dinuc pairs, store decomposition
    bool is_insertion_dinuc_pair = false;
    double kl_combined = 0.0;       // D_KL(truth || inferred) for combined (length + dinuc)
    double kl_dinuc = 0.0;          // D_KL(truth || inferred) for dinuc Markov chain alone
    bool passes_combined_threshold = true;
    bool passes_dinuc_threshold = true;
    double H_combined_truth = 0.0;
    double H_combined_inferred = 0.0;
    double H_len_truth = 0.0;
    double H_len_inferred = 0.0;
    double h_dinuc_truth = 0.0;     // entropy rate
    double h_dinuc_inferred = 0.0;  // entropy rate
    std::string dinuc_nickname;
};

/**
 * @brief Compare inferred model to ground truth using KL divergence
 */
static std::vector<InferenceComparison> compare_inference_to_ground_truth(
        const std::vector<EventInfo>& ground_truth_events,
        const Model_marginals& inferred_marginals,
        const Model_Parms& inferred_parms,
        double kl_threshold_factor = 10.0)
{
    std::vector<InferenceComparison> comparisons;
    auto inferred_index_map = inferred_marginals.get_index_map(inferred_parms);
    
    // Build insertion+dinuc pairs for inferred model (same structure as truth)
    std::map<Gene_class, InsDinucPair> inferred_ins_dinuc_pairs;
    for (const auto& gt_ev : ground_truth_events) {
        if (gt_ev.is_dinuc_markov) {
            inferred_ins_dinuc_pairs[gt_ev.gene_class].dinuc_event = &gt_ev;
        } else if (gt_ev.nickname.find("_ins") != std::string::npos) {
            inferred_ins_dinuc_pairs[gt_ev.gene_class].ins_event = &gt_ev;
        }
    }
    
    for (const auto& gt_ev : ground_truth_events) {
        InferenceComparison cmp;
        cmp.event_nickname = gt_ev.nickname;
        cmp.entropy_truth = gt_ev.H;
        
        if (gt_ev.is_dinuc_markov) {
            // Skip DinucMarkov events - they're shown as part of insertion pairs
            continue;
        }
        
        // Non-DinucMarkov event: standard marginal comparison
        auto [dims, inferred_probs] = 
            inferred_marginals.compute_event_marginal_probability(gt_ev.name, inferred_parms);
        
        std::vector<double> inferred_marginal(gt_ev.num_realizations, 0.0);
        for (int i = 0; i < gt_ev.num_realizations; ++i) {
            inferred_marginal[i] = static_cast<double>(inferred_probs.get()[i]);
        }
        
        // Check if this is an insertion event with a DinucMarkov partner
        auto it = inferred_ins_dinuc_pairs.find(gt_ev.gene_class);
        if (it != inferred_ins_dinuc_pairs.end() &&
            it->second.ins_event == &gt_ev &&
            it->second.dinuc_event != nullptr)
        {
            // Compute combined insertion+dinuc entropy for inferred model
            const auto* dinuc_ev = it->second.dinuc_event;
            
            // Get inferred dinuc transition matrix
            int dinuc_base_idx = inferred_index_map.at(dinuc_ev->name);
            std::array<double, 16> inferred_dinuc_T;
            for (int k = 0; k < 16; ++k) {
                inferred_dinuc_T[k] = static_cast<double>(
                        inferred_marginals.marginal_array_smart_p[dinuc_base_idx + k]);
            }
            
            // Get insertion length values
            auto ev_ptr = inferred_parms.get_event_pointer(gt_ev.name);
            auto realizations = ev_ptr->get_realizations_map();
            std::vector<int> ins_lengths(gt_ev.num_realizations, 0);
            for (const auto &[key, real] : realizations) {
                ins_lengths[real.index] = real.value_int;
            }
            
            // Compute combined entropy for both models
            double H_len_truth = entropy(gt_ev.model_marginal);
            double H_len_inferred = entropy(inferred_marginal);
            
            double combined_H_truth = insertion_dinuc_entropy(
                    gt_ev.model_marginal, ins_lengths, dinuc_ev->dinuc_T);
            double combined_H_inferred = insertion_dinuc_entropy(
                    inferred_marginal, ins_lengths, inferred_dinuc_T);
            
            // Compute combined cross-entropy and D_KL(Truth || Inferred)
            double combined_cross_entropy = insertion_dinuc_cross_entropy(
                    gt_ev.model_marginal, inferred_marginal,
                    ins_lengths, dinuc_ev->dinuc_T, inferred_dinuc_T);
            double kl_combined = combined_cross_entropy - combined_H_truth;
            
            // Store decomposition
            cmp.is_insertion_dinuc_pair = true;
            cmp.H_combined_truth = combined_H_truth;
            cmp.H_combined_inferred = combined_H_inferred;
            cmp.H_len_truth = H_len_truth;
            cmp.H_len_inferred = H_len_inferred;
            cmp.h_dinuc_truth = dinuc_ev->dinuc_entropy_rate;
            cmp.h_dinuc_inferred = markov_entropy_rate(inferred_dinuc_T);
            cmp.dinuc_nickname = dinuc_ev->nickname;
            
            // Compute D_KL(Truth || Inferred) for individual components
            cmp.kl_combined = kl_combined;
            cmp.kl_dinuc = markov_kl_divergence(dinuc_ev->dinuc_T, inferred_dinuc_T);
            
            // Check thresholds
            double combined_threshold = (std::max)(combined_H_truth / kl_threshold_factor, 0.01);
            cmp.passes_combined_threshold = (kl_combined < combined_threshold);
            
            double dinuc_threshold = (std::max)(dinuc_ev->dinuc_entropy_rate / kl_threshold_factor, 0.01);
            cmp.passes_dinuc_threshold = (cmp.kl_dinuc < dinuc_threshold);
            
            // Display combined entropy
            cmp.entropy_truth = combined_H_truth;
            cmp.entropy_inferred = combined_H_inferred;
        } else {
            cmp.entropy_inferred = entropy(inferred_marginal);
        }
        
        // D_KL(truth || inferred) on marginal distribution
        cmp.kl_divergence_forward = kl_divergence(
            gt_ev.model_marginal, inferred_marginal, &cmp.uncovered_mass);
        
        // D_KL(inferred || truth) - useful for symmetric view
        double dummy;
        cmp.kl_divergence_reverse = kl_divergence(
            inferred_marginal, gt_ev.model_marginal, &dummy);
        
        // Check threshold (for insertion events, uses length entropy not combined)
        double threshold = (std::max)(entropy(gt_ev.model_marginal) / kl_threshold_factor, 0.01);
        cmp.passes_threshold = (cmp.kl_divergence_forward < threshold);
        
        comparisons.push_back(cmp);
    }
    
    return comparisons;
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
    // 2. Load gene templates
    // ------------------------------------------------------------------
    std::unordered_map<std::string, std::string> gene_templates;
    
    // Get V gene templates
    auto v_genomic = read_genomic_fasta(MODELS_DIR + "/human/tcr_alpha/ref_genome/genomicVs.fasta");
    for (const auto& [name, seq] : v_genomic) {
        gene_templates[name] = seq;
    }
    
    // Get J gene templates
    auto j_genomic = read_genomic_fasta(MODELS_DIR + "/human/tcr_alpha/ref_genome/genomicJs.fasta");
    for (const auto& [name, seq] : j_genomic) {
        gene_templates[name] = seq;
    }
    
    std::cout << "\nLoaded " << gene_templates.size() << " gene templates" << std::endl;
    
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
    #endif
    
    #ifdef DEBUG_INFERENCE_TEST
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
    
    auto comparisons = compare_inference_to_ground_truth(
        truth_events,
        inferred_marginals,
        inferred_parms,
        kl_threshold_factor);
    
    std::cout << "\nEvent                  | D_KL(T||I) | H_truth | H_infer | Uncovered | Pass" << std::endl;
    std::cout << "---------------------- | ---------- | ------- | ------- | --------- | ----" << std::endl;
    
    for (const auto& cmp : comparisons) {

        
        // If this is an insertion+dinuc pair, show decomposition
        if (cmp.is_insertion_dinuc_pair) {
                        // Main row
            std::cout << std::left << std::setw(22) << cmp.event_nickname << " | "
                    << std::fixed 
                    << "          " << " | "
                    << "       " << " | "
                    << "       " << " | "
                    << "         " << " | "
                    << (" ") << std::endl;
            std::cout << "  ├─ combined           | "
                      << std::setw(10) << cmp.kl_combined << " | "
                      << std::setw(7) << cmp.H_combined_truth << " | "
                      << std::setw(7) << cmp.H_combined_inferred << " | "
                      << "         " << " | "
                      << (cmp.passes_combined_threshold ? "✓" : "✗") << std::endl;
            std::cout << "  ├─ ins length         | "
                      << std::setw(10) << cmp.kl_divergence_forward << " | "
                      << std::setw(7) << cmp.H_len_truth << " | "
                      << std::setw(7) << cmp.H_len_inferred << " | "
                      << std::setw(9) << cmp.uncovered_mass << " | "
                      << (cmp.passes_threshold ? "✓" : "✗") << std::endl;
            std::cout << "  └─ dinuc Markov (h)   | "
                      << std::setw(10) << cmp.kl_dinuc << " | "
                      << std::setw(7) << cmp.h_dinuc_truth << " | "
                      << std::setw(7) << cmp.h_dinuc_inferred << " | "
                      << "         " << " | "
                      << (cmp.passes_dinuc_threshold ? "✓" : "✗") << std::endl;
        } else {
            // Main row
            std::cout << std::left << std::setw(22) << cmp.event_nickname << " | "
                    << std::fixed << std::setprecision(4) << std::setw(10) << cmp.kl_divergence_forward << " | "
                    << std::setw(7) << cmp.entropy_truth << " | "
                    << std::setw(7) << cmp.entropy_inferred << " | "
                    << std::setw(9) << cmp.uncovered_mass << " | "
                    << (cmp.passes_threshold ? "✓" : "✗") << std::endl;
        }
        
        INFO("Event: " << cmp.event_nickname);
        INFO("D_KL(truth||inferred) = " << cmp.kl_divergence_forward);
        INFO("H(truth) = " << cmp.entropy_truth);
        INFO("H(inferred) = " << cmp.entropy_inferred);
        if (cmp.is_insertion_dinuc_pair) {
            INFO("  Combined (len+dinuc) - D_KL=" << cmp.kl_combined
                 << ", H_truth=" << cmp.H_combined_truth
                 << ", H_inferred=" << cmp.H_combined_inferred
                 << ", threshold=" << cmp.H_combined_truth / kl_threshold_factor);
            INFO("  Ins length - D_KL=" << cmp.kl_divergence_forward 
                 << ", H_truth=" << cmp.H_len_truth 
                 << ", H_inferred=" << cmp.H_len_inferred
                 << ", threshold=" << cmp.H_len_truth / kl_threshold_factor);
            INFO("  Dinuc Markov - D_KL=" << cmp.kl_dinuc
                 << ", h_truth=" << cmp.h_dinuc_truth
                 << ", h_inferred=" << cmp.h_dinuc_inferred
                 << ", threshold=" << cmp.h_dinuc_truth / kl_threshold_factor);
        } else {
            INFO("Threshold = " << cmp.entropy_truth / kl_threshold_factor);
        }
        
        CHECK(cmp.passes_threshold);
        if (cmp.is_insertion_dinuc_pair) {
            CHECK(cmp.passes_combined_threshold);
            CHECK(cmp.passes_dinuc_threshold);
        }
    }
    
    std::cout << "\n=== Inference test completed ===" << std::endl;
}
