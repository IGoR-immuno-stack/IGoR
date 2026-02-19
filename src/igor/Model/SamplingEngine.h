#pragma once

#include <igor/Model/SamplingHandler.h>
#include <igor/Model/Topology.h>
#include <igor/Model/Scenario.h>

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <random>

namespace igor::model {

template <typename T = double>
class SamplingEngine {
public:
    using HandlerPtr = std::unique_ptr<SamplingHandler<T>>;

    // ─── Construction ──────────────────────────────────────────────────────────

    /**
     * @brief Constructed with a shared Topology.
     * The topology defines the graph structure and node indexing (UIDs).
     * Automatically creates empty SamplingHandlers for all events in the topology.
     */
    explicit SamplingEngine(std::shared_ptr<const Topology> topology);

    // ─── Configuration ─────────────────────────────────────────────────────────

    /**
     * @brief registerHandler
     * Adds or replaces a sampling handler for a specific event node.
     * Takes ownership. 
     * Verifies that event exists in Topology and assigns its UID to the handler.
     * Throws if handler name doesn't match any event in Topology.
     */
    void registerHandler(HandlerPtr handler);
    
    // ─── Handler Access (Mutable) ──────────────────────────────────────────────
    
    // Access handler by name for filling parameters (e.g. from LegacyBridge)
    SamplingHandler<T>& handler(const std::string& name);
    SamplingHandler<T>& handler(index_type uid);


    // ─── Generation ────────────────────────────────────────────────────────────

    /**
     * @brief generate
     * Produces a single Scenario realization.
     * Iterates through the topology in topological order, sampling each event
     * given its parents' realized values in the current scenario.
     */
    SampledScenario generate(std::mt19937_64& rng) const;
    
    // ─── Accessors ─────────────────────────────────────────────────────────────

    const SamplingHandler<T>& handler(const std::string& name) const;
    const SamplingHandler<T>& handler(index_type uid) const;

private:
    std::shared_ptr<const Topology> m_topology;
    
    // Handlers stored by UID for fast access during generation loop.
    // Size matches topology size. Nullptr if no handler registered for that node.
    std::vector<HandlerPtr> m_handlers;
    
    // Cached topological order for generation loop (optimization)
    std::vector<index_type> m_execution_order;
};

} // namespace igor::model

#include <igor/Model/SamplingEngine.tpp>
