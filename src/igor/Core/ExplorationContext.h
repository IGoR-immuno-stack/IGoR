#pragma once

#include <igor/Core/SafetyMatrix.h>
#include <igor/Core/BoundTightness.h>
#include <igor/Core/Utils.h>
#include <span>

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

    // Which pairs of segments are already known not to overlap (exploration policy)
    // Guides exploration: determines if overlap checks are needed until reaching leaf
    // Once at leaf node, this doesn't describe the scenario - only guided exploration
    SafetyMatrix& safety_set;

    /**
     * @brief NT-floor mismatch positions per Seq_type (conservative pruning bound)
     *
     * floor[seg] holds the positions where the genomic nucleotide is incompatible
     * with EVERY nucleotide the query allows at that position, i.e.
     * !comp_nt_int(gene_nt, query.int_sequence[p]).
     *
     * Invariant: floor[seg] is a subset of mismatches_lists.get(seg) (the upper-bound
     * track in ScenarioContext). For exact NT queries the two tracks are identical,
     * so this is a no-op for standard inference.
     *
     * Consumed by the downstream probability bound: pruning must use the count that
     * cannot be reduced by any choice of query branch, otherwise scenarios that do
     * contribute to the result can be discarded.
     */
    Pruning_mismatch_floor_map& pruning_mismatch_floor;

    /**
     * @brief Constructor
     */
    ExplorationContext(
        Downstream_scenario_proba_bound_map& downstream_proba_map_,
        double& seq_max_prob_scenario_,
        double proba_threshold_factor_,
        Index_map& index_map_,
        std::shared_ptr<Next_event_ptr>& next_event_ptr_arr_,
        SafetyMatrix& safety_set_,
        Pruning_mismatch_floor_map& pruning_mismatch_floor_
    ) : downstream_proba_map(downstream_proba_map_),
        seq_max_prob_scenario(seq_max_prob_scenario_),
        proba_threshold_factor(proba_threshold_factor_),
        index_map(index_map_),
        next_event_ptr_arr(next_event_ptr_arr_),
        safety_set(safety_set_),
        pruning_mismatch_floor(pruning_mismatch_floor_)
    {}

    // Prevent copying
    ExplorationContext(const ExplorationContext&) = delete;
    ExplorationContext& operator=(const ExplorationContext&) = delete;

    // Allow moving
    ExplorationContext(ExplorationContext&&) = default;
    ExplorationContext& operator=(ExplorationContext&&) = delete;

    /**
     * @brief The threshold comparison itself, with nothing observing it
     *
     * Used where the value being tested is a *realized* probability rather than a bound -- the
     * leaf, where the comparison accepts or rejects a completed scenario instead of pruning a
     * subtree. Keeping the two callable separately is what lets the instrumentation count
     * subtree prunings without the leaf test inflating the count.
     */
    inline bool is_below_threshold(double proba) const {
        return proba < (seq_max_prob_scenario * proba_threshold_factor);
    }

    /**
     * @brief Check if scenario probability is above pruning threshold
     *
     * The outcome is offered to the bound instrumentation, which attributes it to the node whose
     * child event took the decision: a node all of whose children were tested and rejected is one
     * a tighter bound would have deleted, and a node whose children were never tested at all died
     * on geometry rather than on probability. `note_prune` is an empty inline unless
     * `IGOR_BOUND_INSTRUMENTATION` is defined.
     */
    inline bool should_prune(double scenario_upper_bound_proba) const {
        const bool prune = is_below_threshold(scenario_upper_bound_proba);
        BoundTightness::note_prune(prune);
        return prune;
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
     * @param downstream_layers Snapshot of per-key memory layers to read the bounds from
     * @return Upper bound probability (base * all downstream bounds)
     */
    inline double compute_upper_bound(
        double base_proba,
        std::span<const int> downstream_layers
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
    // The pair is named by a SafetyCell, resolved once in initialize_event() from the
    // model's 5'->3' ordering. It used to be named by the Event_safety enum, which could
    // spell only the three pairs a single-D model has -- the tandem-D ceiling S5 lifted.
    //
    // PERFORMANCE: All inline - zero overhead.
    // ========================================================================

    /**
     * @brief Is this pair of segments already known not to overlap?
     *
     * @param cell         the pair, resolved at initialization
     * @param memory_layer Memory layer to query -- the caller's own layer minus one, i.e.
     *                     what the enclosing depth left for this pair
     * @return true if no realization still pending can make them overlap
     *
     * PERFORMANCE: Inline, one array read and one bit test
     */
    inline bool is_overlap_safe(
        SafetyCell cell,
        size_t memory_layer
    ) const {
        return safety_set.get(cell, memory_layer);
    }

    /**
     * @brief Record whether this pair is established as non-overlapping
     *
     * A `true` verdict also marks every pair between the same left segment and anything
     * further 3' -- section 2.3's corollary, which moves where such a scenario is
     * discarded without changing whether it is. See SafetyMatrix.
     *
     * @param cell         the pair, resolved at initialization
     * @param is_safe      Whether this arrangement is established safe
     * @param memory_layer Memory layer for storage
     *
     * PERFORMANCE: Inline, one read-modify-write of the row word
     */
    inline void set_overlap_safety(
        SafetyCell cell,
        bool is_safe,
        size_t memory_layer
    ) {
        safety_set.set(cell, is_safe, memory_layer);
    }
};
