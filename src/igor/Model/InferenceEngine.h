/*
 * InferenceEngine.h
 *
 *  Created on: Feb 10, 2026
 *
 *  Central orchestrator for EM inference on marginal parameters.
 *  Owns a collection of MarginalHandler<T> instances (one per event)
 *  and coordinates reset, accumulation, and M-step across all of them.
 *
 *  Template parameter T is the scalar storage type (double or long double).
 */

#pragma once

#include <igor/Model/MarginalHandler.h>
#include <igor/Model/CategoricalHandler.h>
#include <igor/Model/MarkovHandler.h>

#include <unordered_map>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

namespace igor::model {

template <typename T = double>
class InferenceEngine {
public:
    using scalar_type = T;
    using handler_ptr = std::unique_ptr<MarginalHandler<T>>;

    // ─── Construction ──────────────────────────────────────────────────

    /// Default constructor (empty engine, register handlers manually)
    InferenceEngine() = default;

    /// Construct from a list of EventDescriptors
    explicit InferenceEngine(const std::vector<EventDescriptor>& descriptors);

    // ─── Handler Registration ──────────────────────────────────────────

    /// Register a handler for a named event. Takes ownership.
    void register_handler(std::string name, handler_ptr handler);

    // ─── Handler Access ────────────────────────────────────────────────

    /// Access handler by event name (throws if not found)
    MarginalHandler<T>& handler(const std::string& name);
    const MarginalHandler<T>& handler(const std::string& name) const;

    /// Check if a handler exists for the given event name
    bool has_handler(const std::string& name) const;

    /// Number of registered handlers
    std::size_t size() const { return event_order_.size(); }

    // ─── EM Operations ─────────────────────────────────────────────────

    /// Reset all accumulators to zero (E-step preparation)
    void reset_accumulators();

    /// Update all parameters from accumulators (M-step)
    void update_parameters();

    /// Merge another engine's accumulators into this one (for parallel EM).
    /// Both engines must have the same handlers registered in the same order.
    void combine_accumulators(const InferenceEngine<T>& other);


    // ─── I/O ───────────────────────────────────────────────────────────

    /// Write all parameters to stream (in event order)
    void write_marginals(std::ostream& out) const;

    /// Read all parameters from stream (in event order)
    void read_marginals(std::istream& in);

    // ─── Iteration ─────────────────────────────────────────────────────

    /// Ordered event names (registration order, typically topological)
    const std::vector<std::string>& event_names() const { return event_order_; }

    /// Iterate over all handlers (const)
    template <typename Func>
    void for_each_handler(Func&& func) const;

    /// Iterate over all handlers (mutable)
    template <typename Func>
    void for_each_handler(Func&& func);

private:
    std::unordered_map<std::string, handler_ptr> handlers_;
    std::vector<std::string> event_order_;
};

// ─── Type Aliases ──────────────────────────────────────────────────────

/// Legacy-compatible engine using long double (matches Marginal_array_p)
using LegacyEngine = InferenceEngine<long double>;

/// Default engine using double (better performance, sufficient precision)
using Engine = InferenceEngine<double>;

} // namespace igor::model

#include <igor/Model/InferenceEngine.tpp>
