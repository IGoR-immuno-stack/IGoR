#pragma once

#include "igor/Model/Navigator.h"
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
class SamplingEngine
{
public:
    using HandlerPtr = std::unique_ptr<SamplingHandler<T>>;
    using Adjacency_t = Navigator<SamplingHandler<T>, HandlerPtr>;
    using OrderedList = Navigator<SamplingHandler<T>, HandlerPtr>;

    /**
     * @brief Constructed with a shared Topology.
     * The topology defines the graph structure and node indexing (UIDs).
     * Automatically creates empty SamplingHandlers for all events in the topology.
     */
    explicit SamplingEngine(std::shared_ptr<const Topology> topology);

    SamplingHandler<T>& handler(const std::string& name);
    SamplingHandler<T>& handler(index_type uid);

    /**
     * @brief generate
     * Produces a single Scenario realization.
     * Iterates through the topology in topological order, sampling each event
     * given its parents' realized values in the current scenario.
     */
    SampledScenario run(std::mt19937_64& rng) const;

    const SamplingHandler<T>& handler(const std::string& name) const;
    const SamplingHandler<T>& handler(index_type uid) const;

    auto begin(void) const { return m_handlers.begin(); }
    auto end(void) const { return m_handlers.end(); }
    auto begin(void) { return m_handlers.begin(); }
    auto end(void) { return m_handlers.end(); }

    const Topology& topology(void) const { return *m_topology; }

    /**
     * @brief Navigator over all handlers in strictly topological order.
     */
    OrderedList orderedHandlers(void) const;

    /**
     * @brief Navigator over the parents of a specific event.
     */
    Adjacency_t parents(index_type uid) const;

    /**
     * @brief Navigator over the children of a specific event.
     */
    Adjacency_t children(index_type uid) const;

private:
    std::shared_ptr<const Topology> m_topology;

    // Handlers stored by UID for fast access during generation loop.
    // Size matches topology size. Nullptr if no handler registered for that node.
    std::vector<HandlerPtr> m_handlers;

    // Cached topological order for generation loop (optimization)
    std::vector<index_type> m_execution_order;
};

/**
 * @brief Read marginal probabilities from a Model_marginals text file into a SamplingEngine.
 *
 * The topology of @p engine must already be consistent with the model encoded in the
 * file (i.e., same events and graph structure).  The function can be called multiple
 * times on the same engine to hot-swap parameters while keeping the topology intact.
 *
 * @param filename Path to the marginals text file (same format as Model_marginals::txt2marginals).
 * @param engine   SamplingEngine whose handlers will be overwritten with the file values.
 * @return true on success; false if the file could not be opened/parsed.
 * @throws std::runtime_error if the file is consistent but an event or bound check fails.
 */
template <typename T = double>
bool readParameters(const std::string& filename, SamplingEngine<T>& engine);

} // namespace igor::model

#include <igor/Model/SamplingEngine.tpp>
