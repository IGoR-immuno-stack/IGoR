/**
 * @title Parallel EM Validation on Real Data
 * @brief Integration test running side-by-side legacy and parallel EM on the Mouse TCR beta demo dataset.
 *
 * This test loads the full Mouse TCR beta demo model and sequences, runs 2 EM iterations,
 * and asserts numerical equivalence between the legacy Model_marginals and the new InferenceEngine.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <filesystem>
#include <algorithm> // for min
#include <memory>

#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/Singleerrorrate.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Aligner.h>
#include <igor/Core/Utils.h>

#include <igor/Model/InferenceEngine.h>
#include <igor/Core/LegacyBridge.h>

#include "ParallelRecursion.h" // Contains infer_model_parallel_validation

using namespace std;
namespace fs = std::filesystem;

// Helper to construct full path to demo data
#ifndef IGOR_DATA_DIR
    // Fallback relative path assuming running from build/bin
    #define IGOR_DATA_DIR "../../demo"
#endif

TEST_CASE("Phase 2: Parallel EM Validation on Demo Dataset", "[phase2][integration][slow]") {

    // Paths
    string data_dir = IGOR_DATA_DIR;
    // Check if directory exists
    if (!fs::exists(data_dir)) {
        if (fs::exists("../demo")) data_dir = "../demo";
        else if (fs::exists("../../demo")) data_dir = "../../demo";
    }

    // Verify paths exist
    REQUIRE(fs::exists(data_dir + "/genomicVs_with_primers.fasta"));

    // --- 1. Model Setup (Load from test data directory) ---

    // Use the test data directory which has matching model/genomic/sequences
    string test_data_dir = data_dir + "/../scripts/tests/data/input";
    REQUIRE(fs::exists(test_data_dir + "/TRB_model_parms.txt"));

    // Load model parameters
    Model_Parms parms;
    string model_path = test_data_dir + "/TRB_model_parms.txt";

    try {
        parms.read_model_parms(model_path);
    } catch (...) {
        SKIP("Test model files not found at: " + model_path);
    }

    // Load marginals
    Model_marginals legacy_marginals(parms);
    string marginals_path = test_data_dir + "/TRB_uniform_model_marginals.txt";

    try {
        legacy_marginals.txt2marginals(marginals_path, parms);
    } catch (...) {
        SKIP("Test marginals file not found at: " + marginals_path);
    }

    // Load genomic templates from test data directory (matches the model)
    vector<pair<string, string>> v_genomic = read_genomic_fasta(test_data_dir + "/genomicVs_with_primers.fasta");
    vector<pair<string, string>> d_genomic = read_genomic_fasta(test_data_dir + "/genomicDs.fasta");
    vector<pair<string, string>> j_genomic = read_genomic_fasta(test_data_dir + "/genomicJs_all_curated.fasta");

    // --- 2. Load and Align Sequences ---

    // Substitution matrix (nuc44 with ambiguous bases - from demo)
    // A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N
    double nuc44_vect[] = {
        5,   -14, -14, -14, -14, 2,   -14, 2,   2,   -14, -14, 1,   1,   1,   0,
        -14, 5,   -14, -14, -14, 2,   2,   -14, -14, 2,   1,   -14, 1,   1,   0,
        -14, -14, 5,   -14, 2,   -14, 2,   -14, 2,   -14, 1,   1,   -14, 1,   0,
        -14, -14, -14, 5,   2,   -14, -14, 2,   -14, 2,   1,   1,   1,   -14, 0,
        -14, -14, 2,   2,   1.5, -14, -12, -12, -12, -12, 1,   1,   -13, -13, 0,
        2,   2,   -14, -14, -14, 1.5, -12, -12, -12, -12, -13, -13, 1,   1,   0,
        -14, 2,   2,   -14, -12, -12, 1.5, -14, -12, -12, 1,   -13, -13, 1,   0,
        2,   -14, -14, 2,   -12, -12, -14, 1.5, -12, -12, -13, 1,   1,   -13, 0,
        2,   -14, 2,   -14, -12, -12, -12, -12, 1.5, -14, -13, 1,   -13, 1,   0,
        -14, 2,   -14, 2,   -12, -12, -12, -12, -14, 1.5, 1,   -13, 1,   -13, 0,
        -14, 1,   1,   1,   1,   -13, 1,   -13, -13, 1,   0.5, -12, -12, -12, 0,
        1,   -14, 1,   1,   1,   -13, -13, 1,   1,   -13, -12, 0.5, -12, -12, 0,
        1,   1,   -14, 1,   -13, 1,   -13, 1,   -13, 1,   -12, -12, 0.5, -12, 0,
        1,   1,   1,   -14, -13, 1,   1,   -13, 1,   -13, -12, -12, -12, 0.5, 0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
    };
    Matrix<double> nuc44_sub_matrix(15, 15, nuc44_vect);

    // Load sequences path - use test data sequences that match the model
    string seq_path = test_data_dir + "/murugan_naive1_noncoding_demo_seqs.txt";
    REQUIRE(fs::exists(seq_path));

    auto indexed_seqlist = read_txt(seq_path);

    // Aligner Setup
    Aligner v_aligner(nuc44_sub_matrix, 50, V_gene);
    v_aligner.set_genomic_sequences(v_genomic);

    Aligner d_aligner(nuc44_sub_matrix, 50, D_gene);
    d_aligner.set_genomic_sequences(d_genomic);

    Aligner j_aligner(nuc44_sub_matrix, 50, J_gene);
    j_aligner.set_genomic_sequences(j_genomic);

    // Manual Alignment Loop to build sequences_vec
    vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> sequences_vec;
    sequences_vec.reserve(indexed_seqlist.size());

    // Use 5 sequences for meaningful test (changed from 1)
    size_t count_limit = std::min(indexed_seqlist.size(), (size_t)5);

    for(size_t i=0; i<count_limit; ++i) {
        int seq_id = indexed_seqlist[i].first;
        string seq_str = indexed_seqlist[i].second;

        unordered_map<Gene_class, vector<Alignment_data>> aligns_map;

        // Single sequence alignment needs vector wrapping for align_seqs
        vector<pair<const int, const string>> one_seq;
        one_seq.emplace_back(seq_id, seq_str);

        // V (use moderate threshold to ensure we get alignments)
        auto v_map_res = v_aligner.align_seqs(one_seq, 30.0, false);  // Lower threshold, keep all alignments
        if(v_map_res.count(seq_id)) {
            vector<Alignment_data> vec(v_map_res.at(seq_id).begin(), v_map_res.at(seq_id).end());
            aligns_map.emplace(V_gene, vec);
        }

        // D (use moderate threshold with best_only)
        auto d_map_res = d_aligner.align_seqs(one_seq, 20.0, true);
        if(d_map_res.count(seq_id)) {
             vector<Alignment_data> vec(d_map_res.at(seq_id).begin(), d_map_res.at(seq_id).end());
             aligns_map.emplace(D_gene, vec);
        }

        // J (use moderate threshold with best_only to get at least some alignments)
        auto j_map_res = j_aligner.align_seqs(one_seq, 30.0, true);
        if(j_map_res.count(seq_id)) {
             vector<Alignment_data> vec(j_map_res.at(seq_id).begin(), j_map_res.at(seq_id).end());
             aligns_map.emplace(J_gene, vec);
        }

        sequences_vec.emplace_back(seq_id, seq_str, aligns_map);

        std::cout << "Sequence " << i << " (ID=" << seq_id << "): ";
        std::cout << "V=" << (aligns_map.count(V_gene) ? std::to_string(aligns_map.at(V_gene).size()) : "0") << " ";
        std::cout << "D=" << (aligns_map.count(D_gene) ? std::to_string(aligns_map.at(D_gene).size()) : "0") << " ";
        std::cout << "J=" << (aligns_map.count(J_gene) ? std::to_string(aligns_map.at(J_gene).size()) : "0") << " alignments" << std::endl;
    }

    // --- 3. Initialize Inference Engine ---

    // Build engine from parms (legacy bridge helper)
    auto descriptors = igor::core::extract_event_descriptors(parms);
    igor::model::InferenceEngine<long double> engine_instance(descriptors);

    // Import initial uniform values
    igor::core::import_from_legacy(engine_instance, legacy_marginals, parms);

    // --- 4. Run Parallel Validation ---

    SECTION("Two EM iterations produce identical results") {
        std::cout << "=== Starting parallel validation with " << sequences_vec.size() << " sequences ===" << std::endl;
        // Run 2 iterations
        infer_model_parallel_validation(
            sequences_vec,
            parms,
            legacy_marginals,
            engine_instance,
            2,        // n_iterations
            INFINITY  // error threshold
        );
    }
}
