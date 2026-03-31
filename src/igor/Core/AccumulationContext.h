#pragma once

#include <igor/Core/Utils.h>
#include <map>
#include <memory>

// Forward declarations
class Counter;
class Error_rate;

/**
 * @brief Encapsulates result collection at leaf nodes
 * 
 * AccumulationContext holds state for collecting results when we
 * reach complete scenarios (leaf nodes in the scenario tree):
 * - Accumulated posterior marginals
 * - Statistics counters
 * - Error rate accumulation
 * 
 * Note: error_rate is here (not ModelContext) because Error_rate
 * implements accumulation logic that modifies internal state during
 * iterate_wrap_up() calls at leaf nodes. This is a write operation,
 * making it incompatible with ModelContext's immutability.
 * 
 * Separating this from exploration enables:
 * - Testing accumulation logic without exploration policy
 * - Mocking counters for unit tests
 * - Clear separation of "where to explore" vs "what to collect"
 */
struct AccumulationContext {
    // Accumulated posterior marginals (written to at leaf nodes)
    Marginal_array_p& updated_marginals;
    
    // Statistics counters (called at leaf nodes)
    std::map<size_t, std::shared_ptr<Counter>>& counters;
    
    // Error rate accumulation (modified at leaf nodes)
    std::shared_ptr<Error_rate>& error_rate;
    
    /**
     * @brief Constructor
     */
    AccumulationContext(
        Marginal_array_p& updated_marginals_,
        std::map<size_t, std::shared_ptr<Counter>>& counters_,
        std::shared_ptr<Error_rate>& error_rate_
    ) : updated_marginals(updated_marginals_),
        counters(counters_),
        error_rate(error_rate_)
    {}
    
    // Prevent copying
    AccumulationContext(const AccumulationContext&) = delete;
    AccumulationContext& operator=(const AccumulationContext&) = delete;
    
    // Allow moving
    AccumulationContext(AccumulationContext&&) = default;
    AccumulationContext& operator=(AccumulationContext&&) = default;
};
