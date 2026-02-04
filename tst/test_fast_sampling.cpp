/*
 * test_fast_sampling.cpp
 *
 *  Created on: Feb 4, 2026
 *      Unit tests for fast sampling classes
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <igor/Core/FastSampling.h>
#include <igor/Core/FastGenerator.h>
#include <cmath>
#include <map>
#include <numeric>

using namespace igor::fast;
using Catch::Matchers::WithinAbs;

//==============================================================================
// CategoricalSampler Tests
//==============================================================================

TEST_CASE("CategoricalSampler initialization", "[CategoricalSampler]") {
    CategoricalSampler sampler;
    std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
    sampler.initialize(probs, false);

    REQUIRE(sampler.is_initialized());
    REQUIRE(sampler.size() == 4);
}

TEST_CASE("CategoricalSampler with alias method", "[CategoricalSampler]") {
    CategoricalSampler sampler;
    std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
    sampler.initialize(probs, true);

    REQUIRE(sampler.is_initialized());
    REQUIRE(sampler.size() == 4);
}

TEST_CASE("CategoricalSampler sampling distribution", "[CategoricalSampler]") {
    // Test that sampling produces the expected distribution
    CategoricalSampler sampler;
    std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
    sampler.initialize(probs, false);

    std::mt19937_64 rng(42);
    std::map<size_t, int> counts;
    const int n_samples = 100000;

    for (int i = 0; i < n_samples; ++i) {
        size_t idx = sampler.sample(rng);
        counts[idx]++;
    }

    // Check that empirical distribution matches expected
    for (size_t i = 0; i < probs.size(); ++i) {
        double empirical = static_cast<double>(counts[i]) / n_samples;
        REQUIRE_THAT(empirical, WithinAbs(probs[i], 0.01));
    }
}

TEST_CASE("CategoricalSampler alias sampling distribution", "[CategoricalSampler]") {
    // Test alias method produces correct distribution
    CategoricalSampler sampler;
    std::vector<double> probs = {0.1, 0.2, 0.3, 0.4};
    sampler.initialize(probs, true);  // Use alias

    std::mt19937_64 rng(42);
    std::map<size_t, int> counts;
    const int n_samples = 100000;

    for (int i = 0; i < n_samples; ++i) {
        size_t idx = sampler.sample(rng);
        counts[idx]++;
    }

    // Check that empirical distribution matches expected
    for (size_t i = 0; i < probs.size(); ++i) {
        double empirical = static_cast<double>(counts[i]) / n_samples;
        REQUIRE_THAT(empirical, WithinAbs(probs[i], 0.01));
    }
}

TEST_CASE("CategoricalSampler normalizes unnormalized", "[CategoricalSampler]") {
    CategoricalSampler sampler;
    std::vector<double> probs = {1.0, 2.0, 3.0, 4.0};  // Sum = 10
    sampler.initialize(probs, false);

    // Should normalize to 0.1, 0.2, 0.3, 0.4
    REQUIRE_THAT(sampler.probability(0), WithinAbs(0.1, 1e-10));
    REQUIRE_THAT(sampler.probability(1), WithinAbs(0.2, 1e-10));
    REQUIRE_THAT(sampler.probability(2), WithinAbs(0.3, 1e-10));
    REQUIRE_THAT(sampler.probability(3), WithinAbs(0.4, 1e-10));
}

TEST_CASE("CategoricalSampler handles all zeros", "[CategoricalSampler]") {
    CategoricalSampler sampler;
    std::vector<double> probs = {0.0, 0.0, 0.0, 0.0};
    sampler.initialize(probs, false);

    // Should create uniform distribution
    REQUIRE_THAT(sampler.probability(0), WithinAbs(0.25, 1e-10));
    REQUIRE_THAT(sampler.probability(1), WithinAbs(0.25, 1e-10));
    REQUIRE_THAT(sampler.probability(2), WithinAbs(0.25, 1e-10));
    REQUIRE_THAT(sampler.probability(3), WithinAbs(0.25, 1e-10));
}

//==============================================================================
// ConditionalSampler Tests
//==============================================================================

TEST_CASE("ConditionalSampler 2D initialization", "[ConditionalSampler]") {
    ConditionalSampler sampler;
    std::vector<std::vector<double>> probs = {
        {0.1, 0.9},
        {0.5, 0.5},
        {0.9, 0.1}
    };
    sampler.initialize(probs, false);

    REQUIRE(sampler.is_initialized());
    REQUIRE(sampler.num_conditions() == 3);
}

TEST_CASE("ConditionalSampler conditional sampling", "[ConditionalSampler]") {
    ConditionalSampler sampler;
    std::vector<std::vector<double>> probs = {
        {0.1, 0.9},  // condition 0: heavily biased to outcome 1
        {0.9, 0.1}   // condition 1: heavily biased to outcome 0
    };
    sampler.initialize(probs, false);

    std::mt19937_64 rng(42);
    const int n_samples = 10000;

    // Sample from condition 0
    int count_1_given_0 = 0;
    for (int i = 0; i < n_samples; ++i) {
        if (sampler.sample(0, rng) == 1) {
            count_1_given_0++;
        }
    }
    double p_1_given_0 = static_cast<double>(count_1_given_0) / n_samples;
    REQUIRE_THAT(p_1_given_0, WithinAbs(0.9, 0.02));

    // Sample from condition 1
    int count_0_given_1 = 0;
    for (int i = 0; i < n_samples; ++i) {
        if (sampler.sample(1, rng) == 0) {
            count_0_given_1++;
        }
    }
    double p_0_given_1 = static_cast<double>(count_0_given_1) / n_samples;
    REQUIRE_THAT(p_0_given_1, WithinAbs(0.9, 0.02));
}

//==============================================================================
// DinucleotideMarkovSampler Tests
//==============================================================================

TEST_CASE("DinucleotideMarkovSampler initialization", "[DinucleotideMarkovSampler]") {
    DinucleotideMarkovSampler sampler;

    // Create a simple transition matrix (uniform)
    std::vector<double> probs(16, 0.25);  // 4x4 matrix, uniform
    sampler.initialize(probs.data(), 4);

    REQUIRE(sampler.is_initialized());
}

TEST_CASE("DinucleotideMarkovSampler sample next deterministic", "[DinucleotideMarkovSampler]") {
    DinucleotideMarkovSampler sampler;

    // Create transition matrix where each nucleotide always produces itself
    std::vector<double> probs(16, 0.0);
    probs[0*4 + 0] = 1.0;  // P(A|A) = 1
    probs[1*4 + 1] = 1.0;  // P(C|C) = 1
    probs[2*4 + 2] = 1.0;  // P(G|G) = 1
    probs[3*4 + 3] = 1.0;  // P(T|T) = 1
    sampler.initialize(probs.data(), 4);

    std::mt19937_64 rng(42);

    // Given A, should always produce A
    for (int i = 0; i < 100; ++i) {
        REQUIRE(sampler.sample_next(0, rng) == 0);
    }

    // Given T, should always produce T
    for (int i = 0; i < 100; ++i) {
        REQUIRE(sampler.sample_next(3, rng) == 3);
    }
}

TEST_CASE("DinucleotideMarkovSampler generate sequence", "[DinucleotideMarkovSampler]") {
    DinucleotideMarkovSampler sampler;

    // Uniform transitions
    std::vector<double> probs(16, 0.25);
    sampler.initialize(probs.data(), 4);

    std::mt19937_64 rng(42);
    std::string seq = sampler.generate_sequence_str(0, 10, rng);

    REQUIRE(seq.size() == 10);
    // Each character should be a valid nucleotide
    for (char c : seq) {
        REQUIRE((c == 'A' || c == 'C' || c == 'G' || c == 'T'));
    }
}

//==============================================================================
// RealizationBuffer Tests
//==============================================================================

TEST_CASE("RealizationBuffer clear", "[RealizationBuffer]") {
    RealizationBuffer buffer(1000);

    buffer.v_gene_seq_ = "ACGT";
    buffer.final_seq_ = "ACGTACGT";

    buffer.clear();

    REQUIRE(buffer.v_gene_seq_.empty());
    REQUIRE(buffer.final_seq_.empty());
}

//==============================================================================
// BufferPool Tests
//==============================================================================

TEST_CASE("BufferPool acquire and release", "[BufferPool]") {
    BufferPool pool(2, 1000);

    REQUIRE(pool.available() == 2);

    auto buf1 = pool.acquire();
    REQUIRE(pool.available() == 1);

    auto buf2 = pool.acquire();
    REQUIRE(pool.available() == 0);

    // Acquiring when empty should still work (creates new buffer)
    auto buf3 = pool.acquire();
    REQUIRE(pool.available() == 0);

    pool.release(std::move(buf1));
    REQUIRE(pool.available() == 1);

    pool.release(std::move(buf2));
    REQUIRE(pool.available() == 2);
}

//==============================================================================
// Utility Function Tests
//==============================================================================

TEST_CASE("get_optimal_thread_count returns positive", "[Utility]") {
    size_t threads = get_optimal_thread_count();
    REQUIRE(threads > 0);
}

TEST_CASE("draw_random_seed generates different seeds", "[Utility]") {
    uint64_t seed1 = draw_random_seed();
    uint64_t seed2 = draw_random_seed();

    // Seeds should be different (with very high probability)
    REQUIRE(seed1 != seed2);
}

//==============================================================================
// Performance Comparison Test
//==============================================================================

TEST_CASE("Binary search vs alias method produce same distribution", "[Performance]") {
    // Compare sampling performance between binary search and alias method
    std::vector<double> probs(100);
    for (size_t i = 0; i < 100; ++i) {
        probs[i] = static_cast<double>(i + 1);  // Increasing probabilities
    }

    CategoricalSampler binary_sampler;
    binary_sampler.initialize(probs, false);  // Binary search

    CategoricalSampler alias_sampler;
    alias_sampler.initialize(probs, true);  // Alias method

    std::mt19937_64 rng(42);
    const int n_samples = 1000000;

    // Both should produce similar distributions
    std::map<size_t, int> binary_counts, alias_counts;

    for (int i = 0; i < n_samples; ++i) {
        binary_counts[binary_sampler.sample(rng)]++;
        alias_counts[alias_sampler.sample(rng)]++;
    }

    // Verify distributions are similar
    double total_prob = std::accumulate(probs.begin(), probs.end(), 0.0);
    for (size_t i = 0; i < probs.size(); ++i) {
        double expected = probs[i] / total_prob;
        double binary_empirical = static_cast<double>(binary_counts[i]) / n_samples;
        double alias_empirical = static_cast<double>(alias_counts[i]) / n_samples;

        REQUIRE_THAT(binary_empirical, WithinAbs(expected, 0.01));
        REQUIRE_THAT(alias_empirical, WithinAbs(expected, 0.01));
    }
}
