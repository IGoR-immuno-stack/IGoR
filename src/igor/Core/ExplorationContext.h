#pragma once

#include <igor/Core/Utils.h>

/**
 * @brief Encapsulates tree exploration policy and pruning
 * 
 * ExplorationContext holds parameters that guide how we explore the
 * scenario tree:
 * - Downstream probability bounds (for pruning)
 * - Maximum probability scenario found (for adaptive pruning)
 * - Pruning threshold factor
 * - Parent realization indirect tracking (index_map for marginal indexing)
 * 
 * Note: index_map tracks which parent realizations were chosen during
 * exploration (needed for computing indices into marginal arrays). It's
 * modified during recursion as we explore different paths.
 * 
 * Separating this from accumulation enables testing different
 * exploration strategies (aggressive vs conservative pruning)
 * without changing result collection logic.
 */
struct ExplorationContext {
    // Downstream probability bounds for pruning
    // Maps event → max probability of all downstream paths
    Downstream_scenario_proba_bound_map& downstream_proba_map;
    
    // Maximum probability scenario seen for this sequence
    // Used for adaptive pruning threshold
    double& seq_max_prob_scenario;
    
    // Probability threshold factor (scenarios below this are pruned)
    // Scenarios with prob < seq_max_prob_scenario * factor are skipped
    const double proba_threshold_factor;
    
    // Parent event state tracking (which realization chosen at each level)
    // Updated during recursion for marginal array indexing
    Index_map& index_map;
    
    // Next event to call in iteration order (traversal utility)
    // Array indexed by event_index to find which event comes next
    std::shared_ptr<Next_event_ptr>& next_event_ptr_arr;
    
    // Safety flags for gene overlap validation (exploration policy)
    // Used to decide which branches are safe to explore (VD_safe, DJ_safe, VJ_safe)
    // Guides exploration: determines if overlap checks are needed until reaching leaf
    // Once at leaf node, this doesn't describe the scenario - only guided exploration
    Safety_bool_map& safety_set;
    
    /**
     * @brief Constructor
     */
    ExplorationContext(
        Downstream_scenario_proba_bound_map& downstream_proba_map_,
        double& seq_max_prob_scenario_,
        double proba_threshold_factor_,
        Index_map& index_map_,
        std::shared_ptr<Next_event_ptr>& next_event_ptr_arr_,
        Safety_bool_map& safety_set_
    ) : downstream_proba_map(downstream_proba_map_),
        seq_max_prob_scenario(seq_max_prob_scenario_),
        proba_threshold_factor(proba_threshold_factor_),
        index_map(index_map_),
        next_event_ptr_arr(next_event_ptr_arr_),
        safety_set(safety_set_)
    {}
    
    // Prevent copying
    ExplorationContext(const ExplorationContext&) = delete;
    ExplorationContext& operator=(const ExplorationContext&) = delete;
    
    // Allow moving
    ExplorationContext(ExplorationContext&&) = default;
    ExplorationContext& operator=(ExplorationContext&&) = default;
    
    /**
     * @brief Check if scenario probability is above pruning threshold
     */
    inline bool should_prune(double scenario_upper_bound_proba) const {
        return scenario_upper_bound_proba < 
               (seq_max_prob_scenario * proba_threshold_factor);
    }
    
    /**
     * @brief Update maximum probability if current scenario is better
     */
    inline void update_max_prob(double scenario_prob) {
        if (scenario_prob > seq_max_prob_scenario) {
            seq_max_prob_scenario = scenario_prob;
        }
    }
    
    /**
     * @brief Compute upper bound probability for pruning
     * 
     * Multiplies base scenario probability by all downstream probability
     * bounds to get the maximum possible probability for this branch.
     * Used for pruning: if even the best case is below threshold, skip branch.
     * 
     * @param base_proba Base scenario probability before downstream events
     * @param downstream_layers Memory layers to query for downstream bounds
     * @return Upper bound probability (base * all downstream bounds)
     */
    inline double compute_upper_bound(
        double base_proba,
        int* downstream_layers
    ) const {
        double upper_bound = base_proba;
        downstream_proba_map.multiply_all(upper_bound, downstream_layers);
        return upper_bound;
    }
    
    // ========================================================================
    // Gene Overlap Safety Tracking
    // ========================================================================
    // Semantic operations for safety_set (gene overlap validation).
    // These methods guide exploration by marking which gene arrangements
    // are safe to explore (no invalid overlaps).
    // 
    // PERFORMANCE: All inline - zero overhead.
    // ========================================================================
    
    /**
     * @brief Check if gene overlap is safe
     * 
     * Used by Gene_choice to verify V-D or D-J don't overlap in invalid ways.
     * Guides exploration: determines if we can safely explore this branch.
     * 
     * @param safety_type Safety check type (VD_safe, DJ_safe, VJ_safe)
     * @param memory_layer Memory layer to query
     * @return true if overlap is safe, false otherwise
     * 
     * PERFORMANCE: Inline, single map lookup
     */
    inline bool is_overlap_safe(
        Event_safety safety_type,
        size_t memory_layer
    ) const {
        return safety_set.at(safety_type, memory_layer);
    }
    
    /**
     * @brief Mark gene overlap as safe or unsafe
     * 
     * Sets safety flag for future overlap checks during exploration.
     * 
     * @param safety_type Safety check type (VD_safe, DJ_safe, VJ_safe)
     * @param is_safe Whether this arrangement is safe
     * @param memory_layer Memory layer for storage
     * 
     * PERFORMANCE: Inline, single map write
     */
    inline void set_overlap_safety(
        Event_safety safety_type,
        bool is_safe,
        size_t memory_layer
    ) {
        safety_set.set_value(safety_type, is_safe, memory_layer);
    }
};
