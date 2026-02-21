#pragma once

#include "igor/Model/Navigator.h"
#include <igor/Model/SamplingHandler.h>
#include <igor/Model/Scenario.h>

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <random>

namespace igor::model {

template <typename T> class RecombinationModel; // Forward declaration

template <typename T = double>
class SamplingEngine
{
public:
    using HandlerPtr = std::unique_ptr<SamplingHandler<T>>;
    using Adjacency_t = Navigator<SamplingHandler<T>, HandlerPtr>;
    using OrderedList = Navigator<SamplingHandler<T>, HandlerPtr>;

    /**
     * @brief Constructed with a shared RecombinationModel.
     *
     * The model provides both the topology (graph structure) and the
     * probability tensors. Each handler borrows a const reference to its
     * tensor in the model and precomputes its CDF tables immediately.
     *
     * The model must outlive the engine (shared ownership via shared_ptr).
     */
    explicit SamplingEngine(std::shared_ptr<const RecombinationModel<T>> model);

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

    /**
     * @brief Access the underlying RecombinationModel.
     */
    const RecombinationModel<T>& model(void) const { return *m_model; }

private:
    std::shared_ptr<const RecombinationModel<T>> m_model;

    // Handlers stored by UID for fast access during generation loop.
    std::vector<HandlerPtr> m_handlers;

    // Cached topological order for generation loop (optimization)
    std::vector<index_type> m_execution_order;
};

} // namespace igor::model

#include <igor/Model/SamplingEngine.tpp>
