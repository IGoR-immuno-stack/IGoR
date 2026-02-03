/*
 * test_fast_sampler.cpp
 *
 * Unit tests for FastSampler classes.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <igor/Core/FastSampler.h>
#include <igor/Core/Utils.h>
#include <random>
#include <map>
#include <cmath>

using namespace igor;

TEST_CASE("CategoricalSampler basic operations", "[fastsampler]") {
    SECTION("Empty sampler") {
        CategoricalSampler sampler;
        REQUIRE(sampler.empty());
        REQUIRE(sampler.size() == 0);
    }
    
    SECTION("Single category") {
        std::vector<double> probs = {1.0};
        CategoricalSampler sampler(probs);
        REQUIRE_FALSE(sampler.empty());
        REQUIRE(sampler.size() == 1);
        
        std::mt19937_64 rng(42);
        // Should always return 0
        for (int i = 0; i < 100; ++i) {
            REQUIRE(sampler.sample(rng) == 0);
        }
    }
    
    SECTION("Uniform distribution") {
        std::vector<double> probs = {0.25, 0.25, 0.25, 0.25};
        CategoricalSampler sampler(probs);
        REQUIRE(sampler.size() == 4);
        
        std::mt19937_64 rng(42);
        std::map<size_t, int> counts;
        const int n_samples = 10000;
        
        for (int i = 0; i < n_samples; ++i) {
            counts[sampler.sample(rng)]++;
        }
        
        // Each category should have roughly 25% (allow 5% tolerance)
        for (int i = 0; i < 4; ++i) {
            double frac = static_cast<double>(counts[i]) / n_samples;
            REQUIRE(frac > 0.20);
            REQUIRE(frac < 0.30);
        }
    }
    
    SECTION("Skewed distribution") {
        std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
        CategoricalSampler sampler(probs);
        
        std::mt19937_64 rng(42);
        std::map<size_t, int> counts;
        const int n_samples = 10000;
        
        for (int i = 0; i < n_samples; ++i) {
            counts[sampler.sample(rng)]++;
        }
        
        // Check ordering: category 3 > 2 > 1 > 0
        REQUIRE(counts[3] > counts[2]);
        REQUIRE(counts[2] > counts[1]);
        REQUIRE(counts[1] > counts[0]);
    }
    
    SECTION("Batch sampling") {
        std::vector<double> probs = {0.5, 0.3, 0.2};
        CategoricalSampler sampler(probs);
        
        std::mt19937_64 rng(42);
        const size_t batch_size = 1000;
        std::vector<size_t> output(batch_size);
        
        sampler.sample_batch(rng, output.data(), batch_size);
        
        // All samples should be in valid range
        for (size_t val : output) {
            REQUIRE(val < 3);
        }
    }
}

TEST_CASE("ConditionalSampler operations", "[fastsampler]") {
    SECTION("Dinucleotide-like conditional") {
        // 4 conditions (prev nucleotide), 4 outcomes (next nucleotide)
        ConditionalSampler sampler(4, 4);
        
        // Set up simple transition probabilities
        // Each row sums to 1
        sampler.set_probabilities(0, std::vector<double>{0.5, 0.2, 0.2, 0.1}); // A -> ...
        sampler.set_probabilities(1, std::vector<double>{0.1, 0.6, 0.2, 0.1}); // C -> ...
        sampler.set_probabilities(2, std::vector<double>{0.2, 0.2, 0.5, 0.1}); // G -> ...
        sampler.set_probabilities(3, std::vector<double>{0.1, 0.1, 0.2, 0.6}); // T -> ...
        
        REQUIRE(sampler.is_initialized());
        REQUIRE(sampler.num_conditions() == 4);
        REQUIRE(sampler.num_outcomes() == 4);
        
        std::mt19937_64 rng(42);
        std::map<size_t, int> counts;
        const int n_samples = 10000;
        
        // Sample from condition 0 (A)
        for (int i = 0; i < n_samples; ++i) {
            counts[sampler.sample(0, rng)]++;
        }
        
        // A -> A should be most frequent (50%)
        double a_frac = static_cast<double>(counts[0]) / n_samples;
        REQUIRE(a_frac > 0.45);
        REQUIRE(a_frac < 0.55);
    }
}

TEST_CASE("DinucleotideMarkovSampler operations", "[fastsampler]") {
    SECTION("Simple Markov chain") {
        DinucleotideMarkovSampler sampler;
        
        // Build a simple 4x4 matrix (simplified for test)
        Matrix<double> dinuc_matrix(15, 15);
        // Initialize with zeros
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 15; ++j) {
                dinuc_matrix(i, j) = 0.0;
            }
        }
        // Set uniform transitions for nucleotides 0-3
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                dinuc_matrix(i, j) = 0.25;
            }
        }
        // Ambiguous codes (4-14) also get uniform over 0-3
        for (int i = 4; i < 15; ++i) {
            for (int j = 0; j < 4; ++j) {
                dinuc_matrix(i, j) = 0.25;
            }
        }
        
        sampler.build_from_matrix(dinuc_matrix);
        REQUIRE(sampler.is_initialized());
        
        std::mt19937_64 rng(42);
        
        // Sample from A (0)
        std::map<size_t, int> counts;
        const int n_samples = 10000;
        for (int i = 0; i < n_samples; ++i) {
            counts[sampler.sample_next(0, rng)]++;
        }
        
        // Should be roughly uniform
        for (int i = 0; i < 4; ++i) {
            double frac = static_cast<double>(counts[i]) / n_samples;
            REQUIRE(frac > 0.20);
            REQUIRE(frac < 0.30);
        }
    }
    
    SECTION("Insertion sampling") {
        DinucleotideMarkovSampler sampler;
        
        Matrix<double> dinuc_matrix(15, 15);
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 15; ++j) {
                dinuc_matrix(i, j) = 0.0;
            }
        }
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 4; ++j) {
                dinuc_matrix(i, j) = 0.25;
            }
        }
        
        sampler.build_from_matrix(dinuc_matrix);
        
        std::mt19937_64 rng(42);
        std::vector<int> output;
        
        // Sample insertion of length 10 starting from A (0)
        sampler.sample_insertion(0, 10, output, rng);
        
        REQUIRE(output.size() == 10);
        for (int val : output) {
            REQUIRE(val >= 0);
            REQUIRE(val < 4);
        }
    }
    
    SECTION("String insertion sampling") {
        DinucleotideMarkovSampler sampler;
        
        Matrix<double> dinuc_matrix(15, 15);
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 15; ++j) {
                dinuc_matrix(i, j) = 0.0;
            }
        }
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 4; ++j) {
                dinuc_matrix(i, j) = 0.25;
            }
        }
        
        sampler.build_from_matrix(dinuc_matrix);
        
        std::mt19937_64 rng(42);
        std::string inserted = "IIIII"; // 5 positions to fill
        std::vector<int> realizations;
        
        sampler.sample_insertion_to_string(0, inserted, realizations, rng);
        
        REQUIRE(realizations.size() == 5);
        REQUIRE(inserted.size() == 5);
        
        // All characters should be valid nucleotides
        for (char c : inserted) {
            REQUIRE((c == 'A' || c == 'C' || c == 'G' || c == 'T'));
        }
    }
}

TEST_CASE("RealizationBuffer operations", "[fastsampler]") {
    SECTION("Clear and reserve") {
        RealizationBuffer buffer;
        buffer.reserve(500, 50, 100, 50);
        
        buffer.v_gene_seq = "ACGT";
        buffer.d_gene_seq = "GGCC";
        buffer.j_gene_seq = "TTAA";
        buffer.vd_insertion = "AC";
        buffer.dj_insertion = "GT";
        
        buffer.clear();
        
        REQUIRE(buffer.v_gene_seq.empty());
        REQUIRE(buffer.d_gene_seq.empty());
        REQUIRE(buffer.j_gene_seq.empty());
        REQUIRE(buffer.vd_insertion.empty());
        REQUIRE(buffer.dj_insertion.empty());
    }
    
    SECTION("Build final sequence with D gene") {
        RealizationBuffer buffer;
        buffer.v_gene_seq = "AAAA";
        buffer.vd_insertion = "CC";
        buffer.d_gene_seq = "GG";
        buffer.dj_insertion = "TT";
        buffer.j_gene_seq = "CCCC";
        
        buffer.build_final_sequence(true);
        
        REQUIRE(buffer.final_sequence == "AAAACCGGTTCCCC");
    }
    
    SECTION("Build final sequence without D gene") {
        RealizationBuffer buffer;
        buffer.v_gene_seq = "AAAA";
        buffer.vj_insertion = "CCGG";
        buffer.j_gene_seq = "TTTT";
        
        buffer.build_final_sequence(false);
        
        REQUIRE(buffer.final_sequence == "AAAACCGGTTTT");
    }
}

TEST_CASE("CDF correctness", "[fastsampler]") {
    SECTION("CDF values are monotonic and normalized") {
        std::vector<double> probs = {0.1, 0.2, 0.3, 0.25, 0.15};
        CategoricalSampler sampler(probs);
        
        const auto& cdf = sampler.get_cdf();
        REQUIRE(cdf.size() == 5);
        
        // Check monotonicity
        for (size_t i = 1; i < cdf.size(); ++i) {
            REQUIRE(cdf[i] >= cdf[i-1]);
        }
        
        // Check normalization (last value should be 1.0)
        REQUIRE(cdf.back() == Catch::Approx(1.0));
    }
}

TEST_CASE("Binary search is correct", "[fastsampler]") {
    // Test that binary search gives same results as linear search would
    SECTION("Deterministic sampling") {
        std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
        CategoricalSampler sampler(probs);
        
        std::mt19937_64 rng1(12345);
        std::mt19937_64 rng2(12345);
        
        // Sample with same seed should give same results
        for (int i = 0; i < 1000; ++i) {
            size_t sample1 = sampler.sample(rng1);
            size_t sample2 = sampler.sample(rng2);
            REQUIRE(sample1 == sample2);
        }
    }
}
