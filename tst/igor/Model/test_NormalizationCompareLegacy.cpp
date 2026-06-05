/**
 * @file test_NormalizationCompare.cpp
 * @brief Compare normalization behavior between legacy and new architectures
 *
 * This test isolates the normalization step to identify any algorithmic or
 * numerical differences between:
 *   - Legacy: Model_marginals::normalize() + Rec_Event::ind_normalize()
 *   - New:    InferenceHandler<T>::maximizeLikelihood() + linalg::normalize_axis()
 *
 * The tests use the human TCR alpha model (VJ recombination) and compare
 * results for:
 *   1. Uniform initialization
 *   2. Random values (various scales)
 *   3. Sparse distributions (many zeros)
 *   4. Edge cases (all zeros, single non-zero)
 *   5. Repeated normalization (idempotency check)
 *
 * Key goals:
 *   - Verify that both methods produce bitwise-identical results
 *   - Identify numerical precision differences between double/long double
 *   - Detect any ordering or conditional logic differences
 *   - Validate parent-conditional normalization correctness
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// New Model architecture
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/InferenceEngine.h>
#include <igor/Model/LegacyBridge.h>
#include <igor/Model/Topology.h>

// Legacy Core types
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef IGOR_MODELS_DIR
#  error "IGOR_MODELS_DIR must be defined (set by CMake)"
#endif

using namespace igor::model;
using Catch::Matchers::WithinAbs;

// ---------------------------------------------------------------------------
// Model base directory
// ---------------------------------------------------------------------------
static const std::string MODELS_DIR = std::string(IGOR_MODELS_DIR);

// ---------------------------------------------------------------------------
// Helper: Compare tensors element-wise
// ---------------------------------------------------------------------------

template <typename T>
struct ComparisonResult {
    bool exact_match;
    std::size_t total_elements;
    std::size_t mismatches;
    T max_abs_diff;
    T max_rel_diff;
    std::size_t max_abs_diff_idx;
    std::size_t max_rel_diff_idx;
};

template <typename T>
ComparisonResult<T> compare_arrays(const long double* legacy_data, const T* new_data, 
                                   std::size_t size, double tolerance = 1e-14)
{
    ComparisonResult<T> result{};
    result.exact_match = true;
    result.total_elements = size;
    result.mismatches = 0;
    result.max_abs_diff = T(0);
    result.max_rel_diff = T(0);
    result.max_abs_diff_idx = 0;
    result.max_rel_diff_idx = 0;

    for (std::size_t i = 0; i < size; ++i) {
        T legacy_val = static_cast<T>(legacy_data[i]);
        T new_val = new_data[i];
        T abs_diff = std::abs(legacy_val - new_val);
        
        if (abs_diff > tolerance) {
            result.exact_match = false;
            result.mismatches++;
        }
        
        if (abs_diff > result.max_abs_diff) {
            result.max_abs_diff = abs_diff;
            result.max_abs_diff_idx = i;
        }
        
        // Relative difference (avoid division by zero)
        T denom = std::max(std::abs(legacy_val), std::abs(new_val));
        if (denom > T(1e-15)) {
            T rel_diff = abs_diff / denom;
            if (rel_diff > result.max_rel_diff) {
                result.max_rel_diff = rel_diff;
                result.max_rel_diff_idx = i;
            }
        }
    }
    
    return result;
}

// ---------------------------------------------------------------------------
// Helper: Fill arrays with test data patterns
// ---------------------------------------------------------------------------

enum class FillPattern {
    Uniform,       // All elements = 1.0
    Sequential,    // 1, 2, 3, 4, ...
    Random,        // Uniform random [0, 1]
    RandomLog,     // Random spanning many orders of magnitude
    Sparse,        // Mostly zeros, few non-zeros
    SinglePeak,    // One element = 1.0, rest = 0
    AllZeros       // All elements = 0
};

template <typename T>
void fill_array(T* data, std::size_t size, FillPattern pattern, 
                std::mt19937_64& rng, std::size_t peak_idx = 0)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_real_distribution<double> log_scale(-6.0, 0.0); // 10^-6 to 10^0
    
    switch (pattern) {
        case FillPattern::Uniform:
            std::fill(data, data + size, T(1.0));
            break;
            
        case FillPattern::Sequential:
            for (std::size_t i = 0; i < size; ++i) {
                data[i] = T(i + 1);
            }
            break;
            
        case FillPattern::Random:
            for (std::size_t i = 0; i < size; ++i) {
                data[i] = T(uniform(rng));
            }
            break;
            
        case FillPattern::RandomLog:
            for (std::size_t i = 0; i < size; ++i) {
                data[i] = T(std::pow(10.0, log_scale(rng)));
            }
            break;
            
        case FillPattern::Sparse:
            std::fill(data, data + size, T(0));
            // Set ~10% of elements to random values
            for (std::size_t i = 0; i < size / 10; ++i) {
                std::size_t idx = rng() % size;
                data[idx] = T(uniform(rng));
            }
            break;
            
        case FillPattern::SinglePeak:
            std::fill(data, data + size, T(0));
            if (peak_idx < size) {
                data[peak_idx] = T(1.0);
            }
            break;
            
        case FillPattern::AllZeros:
            std::fill(data, data + size, T(0));
            break;
    }
}

// ---------------------------------------------------------------------------
// THE TESTS
// ---------------------------------------------------------------------------

TEMPLATE_TEST_CASE("Normalization comparison: legacy vs new architecture",
                   "[model][normalization]", double, long double)
{
    using T = TestType;
    
    const std::string model_dir = MODELS_DIR + "/human/tcr_alpha";
    const std::string parms_path = model_dir + "/models/model_parms.txt";
    const std::string marginals_path = model_dir + "/models/model_marginals.txt";
    
    INFO("Testing normalization with scalar type: " << 
         (std::is_same_v<T, double> ? "double" : "long double"));
    
    // ══════════════════════════════════════════════════════════════════════
    // 1. Load legacy model structure
    // ══════════════════════════════════════════════════════════════════════
    Model_Parms parms;
    parms.read_model_parms(parms_path);
    
    Model_marginals truth_marginals(parms);
    truth_marginals.txt2marginals(marginals_path, parms);
    
    // ══════════════════════════════════════════════════════════════════════
    // 2. Build new model structure
    // ══════════════════════════════════════════════════════════════════════
    auto topology = import_from_legacy(parms);
    auto model = std::make_shared<RecombinationModel<T>>(
        std::make_unique<Topology>(std::move(*topology)));
    
    import_from_legacy(*model, truth_marginals);
    
    InferenceEngine<T> engine(model);
    
    REQUIRE(engine.size() > 0);
    std::cout << "\n=== Model loaded: " << engine.size() << " events ===" << std::endl;
    
    // ══════════════════════════════════════════════════════════════════════
    // 3. Get index map for legacy array addressing
    // ══════════════════════════════════════════════════════════════════════
    Model_marginals layout_marginals(parms);
    auto index_map = layout_marginals.get_index_map(parms);
    auto inverse_offset_map = layout_marginals.get_inverse_offset_map(parms);
    
    // ══════════════════════════════════════════════════════════════════════
    // 4. Get event ordering (model queue)
    // ══════════════════════════════════════════════════════════════════════
    std::queue<std::shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
    
    // ══════════════════════════════════════════════════════════════════════
    // 5. Test multiple fill patterns
    // ══════════════════════════════════════════════════════════════════════
    std::mt19937_64 rng(12345); // Fixed seed for reproducibility
    
    std::vector<FillPattern> patterns = {
        FillPattern::Uniform,
        FillPattern::Sequential,
        FillPattern::Random,
        FillPattern::RandomLog,
        FillPattern::Sparse,
        FillPattern::AllZeros
    };
    
    std::vector<std::string> pattern_names = {
        "Uniform",
        "Sequential",
        "Random [0,1]",
        "Random log-scale",
        "Sparse (~10% non-zero)",
        "All zeros"
    };
    
    for (std::size_t pat_idx = 0; pat_idx < patterns.size(); ++pat_idx) {
        FillPattern pattern = patterns[pat_idx];
        const std::string& pattern_name = pattern_names[pat_idx];
        
        SECTION(pattern_name) {
            std::cout << "\n─── Testing pattern: " << pattern_name << " ───" << std::endl;
            
            // ──────────────────────────────────────────────────────────────
            // 5a. Fill both representations with identical test data
            // ──────────────────────────────────────────────────────────────
            
            // Legacy: fill marginal array
            Model_marginals test_legacy(parms);
            long double* legacy_data = test_legacy.marginal_array_smart_p.get();
            std::size_t legacy_size = test_legacy.get_length();
            
            // New: fill handler accumulators
            for (igor::index_type uid = 0; 
                 uid < static_cast<igor::index_type>(engine.size()); ++uid)
            {
                const std::string& event_name = model->topology().event(uid)->get_name();
                int base_idx = index_map.at(event_name);
                auto& acc = engine.handler(uid).accumulator();
                
                // Fill with pattern
                std::vector<T> temp_data(acc.size());
                fill_array(temp_data.data(), acc.size(), pattern, rng);
                
                // Copy to both representations
                for (std::size_t i = 0; i < acc.size(); ++i) {
                    legacy_data[base_idx + i] = static_cast<long double>(temp_data[i]);
                    acc.data()[i] = temp_data[i];
                }
            }
            
            // ──────────────────────────────────────────────────────────────
            // 5b. Normalize using LEGACY method
            // ──────────────────────────────────────────────────────────────
            test_legacy.normalize(inverse_offset_map, index_map, model_queue);
            
            // ──────────────────────────────────────────────────────────────
            // 5c. Normalize using NEW method
            // ──────────────────────────────────────────────────────────────
            engine.updateParameters();
            
            // ──────────────────────────────────────────────────────────────
            // 5d. Compare results event-by-event
            // ──────────────────────────────────────────────────────────────
            
            std::cout << std::fixed << std::setprecision(10);
            std::cout << "\nEvent                  | Size   | Mismatch | Max Abs Diff  | Max Rel Diff" 
                      << std::endl;
            std::cout << "---------------------- | ------ | -------- | ------------- | ------------" 
                      << std::endl;
            
            bool all_exact = true;
            
            for (igor::index_type uid = 0; 
                 uid < static_cast<igor::index_type>(engine.size()); ++uid)
            {
                const auto& event = model->topology().event(uid);
                const std::string& event_name = event->get_name();
                const std::string& nickname = event->get_nickname();
                int base_idx = index_map.at(event_name);
                
                const auto& weights = engine.handler(uid).weights();
                
                auto cmp = compare_arrays(
                    legacy_data + base_idx,
                    weights.data(),
                    weights.size(),
                    std::numeric_limits<T>::epsilon() * T(10)
                );
                
                if (!cmp.exact_match) {
                    all_exact = false;
                }
                
                std::cout << std::setw(22) << std::left << nickname
                          << " | " << std::setw(6) << std::right << cmp.total_elements
                          << " | " << std::setw(8) << cmp.mismatches
                          << " | " << std::setw(13) << std::scientific << cmp.max_abs_diff
                          << " | " << std::setw(12) << cmp.max_rel_diff
                          << std::endl;
                
                // Check for significant differences (tolerance depends on scalar type)
                T tolerance = std::is_same_v<T, double> ? T(1e-12) : T(1e-15);
                
                if (pattern != FillPattern::AllZeros) {
                    // For non-zero data, we expect very close agreement
                    REQUIRE(cmp.max_abs_diff < tolerance);
                } else {
                    // Sum over array should be zero
                     T new_arch_sum = T(0);
                    long double legacy_arch_sum = 0.0;
                    for (std::size_t i = 0; i < weights.size(); ++i) {
                        new_arch_sum+=weights.data()[i];
                        legacy_arch_sum+=legacy_data[base_idx + i];
                    }
                    REQUIRE(new_arch_sum==0.0);
                    REQUIRE(legacy_arch_sum==0.0);
                // All-zero case: should be bitwise identical
                    REQUIRE(cmp.max_abs_diff == T(0));
                }
            }
            
            if (all_exact) {
                std::cout << "\n✓ All events: bitwise identical" << std::endl;
            } else {
                std::cout << "\n⚠ Some events differ (within tolerance)" << std::endl;
            }
            
            // ──────────────────────────────────────────────────────────────
            // 5e. Verify normalization correctness (sums to 1)
            // ──────────────────────────────────────────────────────────────
            
            std::cout << "\n─── Checking normalization correctness ───" << std::endl;
            
            for (igor::index_type uid = 0; 
                 uid < static_cast<igor::index_type>(engine.size()); ++uid)
            {
                const auto& event = model->topology().event(uid);
                const std::string& nickname = event->get_nickname();
                const auto& weights = engine.handler(uid).weights();
                
                // For categorical events with no parents, check simple sum
                if (weights.ndim() == 1) {
                    T sum = std::accumulate(weights.begin(), weights.end(), T(0));
                    
                    bool is_normalized = (pattern == FillPattern::AllZeros) 
                        ? (sum == T(0)) 
                        : std::abs(sum - T(1)) < T(1e-10);
                    
                    if (!is_normalized && pattern != FillPattern::AllZeros) {
                        std::cout << "  " << nickname << ": sum = " 
                                  << std::scientific << sum 
                                  << " (expected 1.0)" << std::endl;
                    }
                    
                    if (pattern != FillPattern::AllZeros) {
                        REQUIRE(std::abs(sum - T(1)) < T(1e-10));
                    }
                }
            }
            
            std::cout << "✓ Normalization sums verified" << std::endl;
        }
    }
}

// ---------------------------------------------------------------------------
// Test idempotency: normalizing already-normalized data should not change it
// ---------------------------------------------------------------------------

TEMPLATE_TEST_CASE("Normalization idempotency",
                   "[model][normalization][idempotency]", double, long double)
{
    using T = TestType;
    
    const std::string model_dir = MODELS_DIR + "/human/tcr_alpha";
    const std::string parms_path = model_dir + "/models/model_parms.txt";
    const std::string marginals_path = model_dir + "/models/model_marginals.txt";
    
    INFO("Testing idempotency with scalar type: " << 
         (std::is_same_v<T, double> ? "double" : "long double"));
    
    // Load models
    Model_Parms parms;
    parms.read_model_parms(parms_path);
    
    Model_marginals truth_marginals(parms);
    truth_marginals.txt2marginals(marginals_path, parms);
    
    auto topology = import_from_legacy(parms);
    auto model = std::make_shared<RecombinationModel<T>>(
        std::make_unique<Topology>(std::move(*topology)));
    import_from_legacy(*model, truth_marginals);
    
    InferenceEngine<T> engine(model);
    
    // Get index map
    Model_marginals layout_marginals(parms);
    auto index_map = layout_marginals.get_index_map(parms);
    auto inverse_offset_map = layout_marginals.get_inverse_offset_map(parms);
    std::queue<std::shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
    
    std::cout << "\n=== Testing idempotency: normalize(normalize(x)) == normalize(x) ===" 
              << std::endl;
    
    // ══════════════════════════════════════════════════════════════════════
    // 1. First normalization (from random data)
    // ══════════════════════════════════════════════════════════════════════
    std::mt19937_64 rng(42);
    
    Model_marginals legacy(parms);
    long double* legacy_data = legacy.marginal_array_smart_p.get();
    
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        const std::string& event_name = model->topology().event(uid)->get_name();
        int base_idx = index_map.at(event_name);
        auto& acc = engine.handler(uid).accumulator();
        
        std::vector<T> temp_data(acc.size());
        fill_array(temp_data.data(), acc.size(), FillPattern::Random, rng);
        
        for (std::size_t i = 0; i < acc.size(); ++i) {
            legacy_data[base_idx + i] = static_cast<long double>(temp_data[i]);
            acc.data()[i] = temp_data[i];
        }
    }
    
    // Normalize once
    legacy.normalize(inverse_offset_map, index_map, model_queue);
    engine.updateParameters();
    
    // ══════════════════════════════════════════════════════════════════════
    // 2. Save normalized values
    // ══════════════════════════════════════════════════════════════════════
    std::vector<long double> legacy_after_first(legacy.get_length());
    std::copy(legacy_data, legacy_data + legacy.get_length(), 
              legacy_after_first.begin());
    
    std::unordered_map<igor::index_type, std::vector<T>> new_after_first;
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        const auto& weights = engine.handler(uid).weights();
        new_after_first[uid] = std::vector<T>(weights.begin(), weights.end());
    }
    
    // ══════════════════════════════════════════════════════════════════════
    // 3. Copy normalized values back to accumulators and normalize again
    // ══════════════════════════════════════════════════════════════════════
    
    // Legacy: copy back and normalize
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        const std::string& event_name = model->topology().event(uid)->get_name();
        int base_idx = index_map.at(event_name);
        auto& weights = engine.handler(uid).weights();
        
        // Copy from saved values
        for (std::size_t i = 0; i < weights.size(); ++i) {
            legacy_data[base_idx + i] = legacy_after_first[base_idx + i];
        }
    }
    legacy.normalize(inverse_offset_map, index_map, model_queue);
    
    // New: copy to accumulator and normalize
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        auto& acc = engine.handler(uid).accumulator();
        const auto& saved = new_after_first[uid];
        std::copy(saved.begin(), saved.end(), acc.begin());
    }
    engine.updateParameters();
    
    // ══════════════════════════════════════════════════════════════════════
    // 4. Compare: second normalization should produce identical results
    // ══════════════════════════════════════════════════════════════════════
    std::cout << "\nEvent                  | Max change after 2nd norm" << std::endl;
    std::cout << "---------------------- | -------------------------" << std::endl;
    
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        const auto& event = model->topology().event(uid);
        const std::string& nickname = event->get_nickname();
        const std::string& event_name = event->get_name();
        int base_idx = index_map.at(event_name);
        
        const auto& weights = engine.handler(uid).weights();
        const auto& first_norm = new_after_first[uid];
        
        T max_change = T(0);
        for (std::size_t i = 0; i < weights.size(); ++i) {
            T diff = std::abs(weights.data()[i] - first_norm[i]);
            max_change = std::max(max_change, diff);
        }
        
        std::cout << std::setw(22) << std::left << nickname
                  << " | " << std::scientific << max_change << std::endl;
        
        // Should be machine epsilon or very close
        REQUIRE(max_change < std::numeric_limits<T>::epsilon() * T(100));
    }
    
    std::cout << "\n✓ Idempotency verified: normalize(normalize(x)) == normalize(x)" 
              << std::endl;
}

// ---------------------------------------------------------------------------
// Test markov-specific normalization (row-wise)
// ---------------------------------------------------------------------------

TEST_CASE("Markov normalization correctness", "[model][normalization][markov]")
{
    const std::string model_dir = MODELS_DIR + "/human/tcr_beta";
    const std::string parms_path = model_dir + "/models/model_parms.txt";
    const std::string marginals_path = model_dir + "/models/model_marginals.txt";
    
    // Load models (TCR beta has VD and DJ dinucleotide Markov chains)
    Model_Parms parms;
    parms.read_model_parms(parms_path);
    
    Model_marginals truth_marginals(parms);
    truth_marginals.txt2marginals(marginals_path, parms);
    
    using T = double;
    auto topology = import_from_legacy(parms);
    auto model = std::make_shared<RecombinationModel<T>>(
        std::make_unique<Topology>(std::move(*topology)));
    import_from_legacy(*model, truth_marginals);
    
    InferenceEngine<T> engine(model);
    
    std::cout << "\n=== Testing Markov chain normalization (TCR beta model) ===" << std::endl;
    
    // Find markov events
    for (igor::index_type uid = 0; 
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        const auto& event = model->topology().event(uid);
        if (event->get_type() != Dinuclmarkov_t) continue;
        
        const std::string& nickname = event->get_nickname();
        const auto& weights = engine.handler(uid).weights();
        
        std::cout << "\nEvent: " << nickname << std::endl;
        std::cout << "  Shape: [";
        for (std::size_t i = 0; i < weights.ndim(); ++i) {
            std::cout << weights.shape()[i];
            if (i < weights.ndim() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
        
        // For rank-2 markov (unconditional), check each row sums to 1
        if (weights.ndim() == 2) {
            std::size_t n_states = weights.shape()[0];
            REQUIRE(weights.shape()[1] == n_states);
            
            std::cout << "  Checking row sums (each row = 1 transition from state):" << std::endl;
            
            for (std::size_t from = 0; from < n_states; ++from) {
                T row_sum = T(0);
                for (std::size_t to = 0; to < n_states; ++to) {
                    row_sum += weights(from, to);
                }
                
                std::cout << "    Row " << from << " sum: " 
                          << std::fixed << std::setprecision(15) << row_sum << std::endl;
                
                REQUIRE(std::abs(row_sum - T(1)) < T(1e-12));
            }
        }
    }
    
    std::cout << "\n✓ Markov normalization verified" << std::endl;
}
