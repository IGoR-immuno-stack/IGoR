/**
 * @file test_platform_diagnostic.cpp
 * @brief Minimal diagnostic test to trace numerical differences across platforms.
 *
 * Purpose: identify whether per-sequence likelihoods differ between macOS and
 * Ubuntu (small FP precision differences in the iterate() accumulation), or
 * whether only the output ordering differs (unordered_map non-determinism).
 *
 * Usage:
 *   OMP_NUM_THREADS=1 ./igor_tests "[platform-diagnostic]" 2>/dev/null \
 *       | grep -v "^Filters:\|^Randomness" \
 *       | tee /tmp/igor_platform_diag.txt
 *
 * On the other platform:
 *   OMP_NUM_THREADS=1 ./igor_tests "[platform-diagnostic]" 2>/dev/null \
 *       | grep -v "^Filters:\|^Randomness" \
 *       > /tmp/igor_platform_diag_other.txt
 *   diff /tmp/igor_platform_diag.txt /tmp/igor_platform_diag_other.txt
 *
 * Output is sorted by (event_name, realization_name) and seq_index,
 * so it is independent of unordered_map iteration order.
 *
 * Inference parameters match the igor CLI defaults (main.cpp lines 201-202):
 *   likelihood_threshold = 1e-60, proba_threshold_factor = 1e-5, fast_iter = true
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Core/Aligner.h>
#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#ifndef IGOR_SOURCE_DIR
#error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

static const std::string REG_INPUT_DIR =
        std::string(IGOR_SOURCE_DIR) + "/scripts/tests/data/input";
static const std::string REG_ALIGNS_DIR =
        std::string(IGOR_SOURCE_DIR) + "/scripts/tests/data/reference/aligns";

// ---------------------------------------------------------------------------
// Parse inference_logs.txt (after column-2 removal by test script)
// Returns map: seq_index -> seq_likelihood (for iteration 0)
// ---------------------------------------------------------------------------
static std::map<int, double> parse_inference_logs(const std::string &path)
{
    std::map<int, double> result;
    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // header
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string iter_s, seq_proc_s, seq_idx_s, seq_s, nv_s, nj_s, lkl_s;
        std::getline(ss, iter_s, ';');
        std::getline(ss, seq_proc_s, ';');
        std::getline(ss, seq_idx_s, ';');
        std::getline(ss, seq_s, ';');
        std::getline(ss, nv_s, ';');
        std::getline(ss, nj_s, ';');
        std::getline(ss, lkl_s, ';');
        if (iter_s != "0") continue;
        int seq_idx = std::stoi(seq_idx_s);
        double lkl  = std::stod(lkl_s);
        result[seq_idx] = lkl;
    }
    return result;
}

TEST_CASE("Platform diagnostic - per-sequence likelihoods and marginals",
          "[platform-diagnostic]")
{
    // Force single-threaded (set before GenModel runs)
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif

    // ------------------------------------------------------------------
    // 1. Load model
    // ------------------------------------------------------------------
    Model_Parms parms;
    parms.read_model_parms(REG_INPUT_DIR + "/TRB_model_parms.txt");
    Model_marginals marginals(parms);
    marginals.txt2marginals(REG_INPUT_DIR + "/TRB_uniform_model_marginals.txt", parms);

    // ------------------------------------------------------------------
    // 2. Load alignments for all 300 default sequences
    // ------------------------------------------------------------------
    auto seqs = read_indexed_csv(REG_ALIGNS_DIR + "/default_indexed_sequences.csv");

    auto alignments = read_alignments_seq_csv_score_range(
            REG_ALIGNS_DIR + "/default_V_alignments.csv", V_gene, 55.0, false, seqs);
    alignments = read_alignments_seq_csv_score_range(
            REG_ALIGNS_DIR + "/default_D_alignments.csv", D_gene, 35.0, false, seqs, alignments);
    alignments = read_alignments_seq_csv_score_range(
            REG_ALIGNS_DIR + "/default_J_alignments.csv", J_gene, 10.0, false, seqs, alignments);

    auto sequences_vec = map2vect(alignments);

    // ------------------------------------------------------------------
    // 3. Run 1 EM iteration, write to temp dir
    // ------------------------------------------------------------------
    const std::string out_dir = "/tmp/igor_platform_diag/";
    std::filesystem::create_directories(out_dir);

    GenModel gen_model(parms, marginals);
    // Match defaults from igor CLI (main.cpp lines 201-202)
    bool failed = gen_model.infer_model(
            sequences_vec,
            /*iterations=*/1,
            out_dir,
            /*fast_iter=*/true,
            /*likelihood_threshold=*/1e-60,
            /*viterbi_like=*/false,
            /*proba_threshold_factor=*/1e-5);

    REQUIRE_FALSE(failed);

    // ------------------------------------------------------------------
    // 4. Parse per-sequence likelihoods from inference_logs.txt
    //    Print sorted by seq_index
    // ------------------------------------------------------------------
    auto likelihoods = parse_inference_logs(out_dir + "inference_logs.txt");

    std::cout << "\n=== PER-SEQUENCE LIKELIHOODS (iteration 0, sorted by seq_index) ===\n";
    std::cout << std::setprecision(15);
    for (const auto &[idx, lkl] : likelihoods) {
        std::cout << "seq_idx=" << idx << " likelihood=" << lkl << "\n";
    }

    // ------------------------------------------------------------------
    // 5. Print per-event marginals after iteration 1, sorted by
    //    (event_name, realization_name) for platform-independent comparison
    // ------------------------------------------------------------------
    const Model_marginals &inferred = gen_model.get_marginals();
    auto index_map = inferred.get_index_map(parms);

    // Collect: (event_name, realization_name, marginal_value)
    std::vector<std::tuple<std::string, std::string, long double>> marginal_rows;

    for (const auto &event_ptr : parms.get_event_list()) {
        const std::string &ev_name = event_ptr->get_name();

        // Only consider events without parent dependencies for simplicity
        // (v_choice, j_choice, insertion length, dinucl have simple structure)
        // For events with parents (d_gene depends on j_choice), the marginals
        // are multi-dimensional — skip for this minimal diagnostic.
        //
        // v_choice and j_choice in TRB have no parents → 1-D marginals.

        auto it = index_map.find(ev_name);
        if (it == index_map.end()) continue;
        int base = it->second;

        for (const auto &[real_name, real] : event_ptr->get_realizations_map()) {
            long double val = inferred.marginal_array_smart_p[base + real.index];
            marginal_rows.emplace_back(ev_name, real_name, val);
        }
    }

    // Sort by (event_name, realization_name)
    std::sort(marginal_rows.begin(), marginal_rows.end());

    std::cout << "\n=== PER-EVENT MARGINALS after 1 iteration (sorted, simple events only) ===\n";
    std::cout << "event_name;realization_name;marginal_value\n";
    for (const auto &[ev, real, val] : marginal_rows) {
        std::cout << ev << ";" << real << ";" << std::setprecision(15) << val << "\n";
    }
}
