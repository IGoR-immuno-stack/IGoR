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
    
    /**
     * @brief Constructor
     */
    ExplorationContext(
        Downstream_scenario_proba_bound_map& downstream_proba_map_,
        double& seq_max_prob_scenario_,
        double proba_threshold_factor_,
        Index_map& index_map_,
        std::shared_ptr<Next_event_ptr>& next_event_ptr_arr_
    ) : downstream_proba_map(downstream_proba_map_),
        seq_max_prob_scenario(seq_max_prob_scenario_),
        proba_threshold_factor(proba_threshold_factor_),
        index_map(index_map_),
        next_event_ptr_arr(next_event_ptr_arr_)
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
};
