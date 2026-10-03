// RecombinationModel.h ---

#pragma once

#include <igor/Model/Navigator.h>
#include <igor/Model/Topology.h>
#include <igor/Math/Tensor.h>
#include <igor/Core/Legacy/Typedef.h>

#include <memory>
#include <string>
#include <vector>

namespace igor::model {

// ─── RecombinationModel<T> ───────────────────────────────────────────────────
//
// Pairs a Topology (graph structure) with one probability tensor per node.
//
// This is the central data structure for the VDJ recombination model:
//   - Topology provides the directed acyclic graph of recombination events,
//     event metadata, and topological ordering.
//   - The weight vector stores the conditional probability table (CPT) for
//     each event, indexed by the topology UID.
//
// Tensor shapes are computed from the Topology at construction time:
//   [event_dim0, ..., parent1_dim, parent2_dim, ...]
//
// RecombinationModel itself can be shared between a SamplingEngine and
// an InferenceEngine via std::shared_ptr<RecombinationModel<T>>.

template <typename T = double>
class RecombinationModel
{
public:
    using OrderedList = Navigator<math::Tensor<T>, math::Tensor<T>>;

    /// Build from a Topology — creates zero-initialised tensors with correct shapes.
    explicit RecombinationModel(std::unique_ptr<Topology> topology);

    // ── Topology access ──────────────────────────────────────────────────────

    Topology&       topology(void)       { return *m_topology; }
    const Topology& topology(void) const { return *m_topology; }

    // ── Weight access by UID ─────────────────────────────────────────────────

    math::Tensor<T>&       weight(igor::core::legacy::index_type uid);
    const math::Tensor<T>& weight(igor::core::legacy::index_type uid) const;

    // ── Weight access by event nickname ──────────────────────────────────────

    math::Tensor<T>&       weight(const std::string& name);
    const math::Tensor<T>& weight(const std::string& name) const;

    // ── Ordered traversal ─────────────────────────────────────────────────────

    /// Navigator over all weights in strictly topological order.
    OrderedList orderedWeights(void) const;

    // ── Size / iteration ─────────────────────────────────────────────────────

    std::size_t size(void) const { return m_weights.size(); }

    auto begin(void)       { return m_weights.begin(); }
    auto end(void)         { return m_weights.end(); }
    auto begin(void) const { return m_weights.begin(); }
    auto end(void)   const { return m_weights.end(); }

private:
    std::unique_ptr<Topology> m_topology;
    std::vector<math::Tensor<T>>     m_weights;         // indexed by topology UID
    std::vector<igor::core::legacy::index_type>    m_execution_order;  // cached topological order
};

// ─── Free functions ──────────────────────────────────────────────────────────

/// Build a fully-loaded RecombinationModel in one step from files.
/// Core reads both files, LegacyBridge builds the Topology and fills the tensors from the
/// flat marginal array. There is deliberately no second reader of model_marginals here:
/// Core's renormalizes after parsing, and its values are the ones inference uses.
template <typename T = double>
RecombinationModel<T> recombination_model_from_files(
    const std::string& file_model_parms,
    const std::string& file_model_marginals);

} // namespace igor::model

#include <igor/Model/RecombinationModel.tpp>

//
// RecombinationModel.h ends here
