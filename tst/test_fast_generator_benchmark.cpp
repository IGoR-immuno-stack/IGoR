/*
 * test_fast_generator_benchmark.cpp
 *
 * Benchmark tests comparing fast vs original generation.
 * These tests are marked [benchmark] and can be run separately.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <igor/Core/FastSampler.h>
#include <igor/Core/Utils.h>
#include <random>
#include <chrono>
#include <iostream>

using namespace igor;

TEST_CASE("CategoricalSampler performance", "[benchmark][fastsampler]") {
    // Compare binary search vs linear search
    std::vector<double> probs(100);
    for (int i = 0; i < 100; ++i) {
        probs[i] = 1.0 / 100.0;
    }
    
    CategoricalSampler sampler(probs);
    std::mt19937_64 rng(42);
    
    BENCHMARK("Binary search sampling (100 categories, 10000 samples)") {
        size_t sum = 0;
        for (int i = 0; i < 10000; ++i) {
            sum += sampler.sample(rng);
        }
        return sum;
    };
    
    BENCHMARK("Binary search sampling (100 categories, 100000 samples)") {
        size_t sum = 0;
        for (int i = 0; i < 100000; ++i) {
            sum += sampler.sample(rng);
        }
        return sum;
    };
}

TEST_CASE("DinucleotideMarkovSampler performance", "[benchmark][fastsampler]") {
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
    
    BENCHMARK("Markov insertion sampling (length 20, 10000 times)") {
        std::vector<int> output;
        for (int i = 0; i < 10000; ++i) {
            sampler.sample_insertion(0, 20, output, rng);
        }
        return output.size();
    };
    
    BENCHMARK("Markov string insertion (length 20, 10000 times)") {
        std::string inserted(20, 'I');
        std::vector<int> realizations;
        for (int i = 0; i < 10000; ++i) {
            std::fill(inserted.begin(), inserted.end(), 'I');
            sampler.sample_insertion_to_string(0, inserted, realizations, rng);
        }
        return realizations.size();
    };
}

TEST_CASE("RealizationBuffer performance", "[benchmark][fastsampler]") {
    BENCHMARK("Buffer clear and build sequence (10000 times)") {
        RealizationBuffer buffer;
        buffer.reserve(500, 50, 100, 50);
        
        for (int i = 0; i < 10000; ++i) {
            buffer.clear();
            buffer.v_gene_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
            buffer.d_gene_seq = "GGCCGGCCGGCC";
            buffer.j_gene_seq = "TTAATTAATTAATTAATTAATTAATTAA";
            buffer.vd_insertion = "ACGTACGT";
            buffer.dj_insertion = "GTACGTAC";
            buffer.build_final_sequence(true);
        }
        return buffer.final_sequence.size();
    };
}

// Manual timing test for detailed performance analysis
TEST_CASE("Detailed timing analysis", "[benchmark][fastsampler]") {
    const int NUM_SAMPLES = 1000000;
    
    // Test CategoricalSampler with different sizes
    SECTION("CategoricalSampler scaling") {
        std::vector<int> sizes = {4, 16, 64, 256};
        
        for (int size : sizes) {
            std::vector<double> probs(size, 1.0 / size);
            CategoricalSampler sampler(probs);
            std::mt19937_64 rng(42);
            
            auto start = std::chrono::high_resolution_clock::now();
            size_t sum = 0;
            for (int i = 0; i < NUM_SAMPLES; ++i) {
                sum += sampler.sample(rng);
            }
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            double samples_per_sec = NUM_SAMPLES * 1e6 / duration.count();
            
            std::cout << "CategoricalSampler(" << size << " categories): " 
                      << std::fixed << std::setprecision(0) << samples_per_sec 
                      << " samples/sec" << std::endl;
            
            REQUIRE(sum > 0);  // Prevent optimization
        }
    }
    
    SECTION("DinucleotideMarkovSampler insertion lengths") {
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
        
        std::vector<int> lengths = {5, 10, 20, 50};
        std::mt19937_64 rng(42);
        std::vector<int> output;
        
        for (int length : lengths) {
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < NUM_SAMPLES / 10; ++i) {
                sampler.sample_insertion(0, length, output, rng);
            }
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            double inserts_per_sec = (NUM_SAMPLES / 10) * 1e6 / duration.count();
            
            std::cout << "Markov insertion(length=" << length << "): " 
                      << std::fixed << std::setprecision(0) << inserts_per_sec 
                      << " insertions/sec" << std::endl;
            
            REQUIRE(output.size() == static_cast<size_t>(length));
        }
    }
}
