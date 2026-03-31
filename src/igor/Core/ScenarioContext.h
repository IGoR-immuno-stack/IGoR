#pragma once

#include <igor/Core/Utils.h>
#include <igor/Core/IntStr.h>

/**
 * @brief Encapsulates per-path mutable state during scenario exploration
 * 
 * ScenarioContext holds state that changes as we recursively explore
 * possible recombination scenarios. Critically, this keeps the hot-path
 * unified: probability computation and sequence building are handled
 * together (accessed in tight loops).
 * 
 * Contains:
 * - Running probability of current scenario path
 * - Partially constructed sequences (gene/junction segments)
 * - Offsets of sequence fragments (5'/3' positions)
 * - Mismatches accumulated for error rate calculation
 * 
 * Design principle: Zero-copy - store references/pointers to minimize
 * overhead in hot path (called millions of times per sequence).
 */
struct ScenarioContext {
    // Scenario probability (multiplied at each recursion level)
    double& scenario_proba;
    
    // Constructed sequences (built up during recursion)
    // Maps Seq_type → string pointer with memory layer
    Seq_type_str_p_map& constructed_sequences;
    
    // Offsets of each sequence fragment (5'/3' positions)
    Seq_offsets_map& seq_offsets;
    
    // Mismatch positions for error rate calculation
    Mismatch_vectors_map& mismatches_lists;
    
    /**
     * @brief Constructor - binds references to existing state
     * 
     * All parameters are references to existing state that will be
     * modified during scenario exploration.
     */
    ScenarioContext(
        double& scenario_proba_,
        Seq_type_str_p_map& constructed_sequences_,
        Seq_offsets_map& seq_offsets_,
        Mismatch_vectors_map& mismatches_lists_
    ) : scenario_proba(scenario_proba_),
        constructed_sequences(constructed_sequences_),
        seq_offsets(seq_offsets_),
        mismatches_lists(mismatches_lists_)
    {}
    
    // Prevent copying (would create dangling references)
    ScenarioContext(const ScenarioContext&) = delete;
    ScenarioContext& operator=(const ScenarioContext&) = delete;
    
    // Allow moving (just rebinds references)
    ScenarioContext(ScenarioContext&&) = default;
    ScenarioContext& operator=(ScenarioContext&&) = default;
};
