/*
 * InferenceEngine.h
 *
 *  Created on: Feb 10, 2026
 *
 *  Central orchestrator for EM inference on marginal parameters.
 *  Owns a collection of InferenceHandler<T> instances (one per event)
 *  and coordinates reset, accumulation, and M-step across all of them.
 *
 *  Each handler borrows a mutable reference to its weight tensor from the
 *  shared RecombinationModel, so the M-step writes normalised values
 *  directly back into the model.
 *
 *  Template parameter T is the scalar storage type (double or long double).
 */

#pragma once

#include <igor/Model/InferenceHandler.h>
#include <igor/Model/Navigator.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/InferenceHandlerFactory.h>

#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

namespace igor::model {

template <typename T>
class InferenceEngine
{
public:
    using scalar_type = T;
    using HandlerPtr = std::unique_ptr<InferenceHandler<T>>;
    using Adjacency_t = Navigator<InferenceHandler<T>, HandlerPtr>;
    using OrderedList = Navigator<InferenceHandler<T>, HandlerPtr>;

    // ─── Construction ──────────────────────────────────────────────────

    /// Build from a shared RecombinationModel — one handler per topology node.
    /// Handlers borrow mutable references to the model's weight tensors.
    explicit InferenceEngine(std::shared_ptr<RecombinationModel<T>> model);

    // ─── Handler Access ────────────────────────────────────────────────

    /// Access handler by event name (throws if not found)
          InferenceHandler<T>& handler(const std::string& name);
    const InferenceHandler<T>& handler(const std::string& name) const;

    /// Access handler by topology uid
          InferenceHandler<T>& handler(igor::index_type uid);
    const InferenceHandler<T>& handler(igor::index_type uid) const;

    /// Check if a handler exists for the given event name
    bool hasHandler(const std::string& name) const;

    /// Number of registered handlers
    std::size_t size(void) const;

    // ─── Model Access ──────────────────────────────────────────────────

    /// The underlying RecombinationModel
    const RecombinationModel<T>& model(void) const;
          RecombinationModel<T>& model(void);

    // ─── EM Operations ─────────────────────────────────────────────────

    /**
     * @brief Run one full EM iteration.
     *
     * Executes the three phases in order:
     *   1. resetAccumulators()
     *   2. eStep(*this)          — caller-provided accumulation logic
     *   3. updateParameters()    — M-step: normalise accumulators → weights
     *
     * @param eStep  Callable with signature `void(InferenceEngine<T>&)`
     *               that accumulates counts into each handler's accumulator.
     */
    template <typename Func>
    void run(Func&& eStep);

    /// Reset all accumulators to zero (E-step preparation)
    void resetAccumulators(void);

    /// Update all parameters from accumulators (M-step)
    void updateParameters(void);

    /// Merge another engine's accumulators into this one (for parallel EM).
    /// Both engines must have handlers for events with the same names.
    void combineAccumulators(const InferenceEngine<T>& other);

    // ─── Iteration ─────────────────────────────────────────────────────

    auto begin(void) const { return m_handlers.begin(); }
    auto end(void)   const { return m_handlers.end(); }
    auto begin(void)       { return m_handlers.begin(); }
    auto end(void)         { return m_handlers.end(); }

    /// Navigator over all handlers in strictly topological order.
    OrderedList orderedHandlers(void) const;

    /// Navigator over the parents of a specific event.
    Adjacency_t parents(igor::index_type uid) const;

    /// Navigator over the children of a specific event.
    Adjacency_t children(igor::index_type uid) const;

    /// Iterate over all handlers in topological order (const)
    template <typename Func>
    void forEachHandler(Func&& func) const;

    /// Iterate over all handlers in topological order (mutable)
    template <typename Func>
    void forEachHandler(Func&& func);

private:
    std::shared_ptr<RecombinationModel<T>> m_model;
    std::vector<HandlerPtr>                m_handlers;        // indexed by topology uid
    std::vector<igor::index_type>          m_execution_order;  // topological order
};

// ─── Type Aliases ──────────────────────────────────────────────────────

/// Legacy-compatible engine using long double (matches Marginal_array_p)
using LegacyEngine = InferenceEngine<long double>;

/// Default engine using double (better performance, sufficient precision)
using Engine = InferenceEngine<double>;

} // namespace igor::model

#include <igor/Model/InferenceEngine.tpp>
