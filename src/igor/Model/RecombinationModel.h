// RecombinationModel.h ---

#pragma once

#include <igor/Model/Navigator.h>
#include <igor/Model/Topology.h>
#include <igor/Math/Tensor.h>
#include <igor/Core/Typedef.h>

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
    explicit RecombinationModel(std::unique_ptr<const Topology> topology);

    // ── Topology access ──────────────────────────────────────────────────────

    const Topology& topology(void) const { return *m_topology; }

    // ── Weight access by UID ─────────────────────────────────────────────────

    math::Tensor<T>&       weight(igor::index_type uid);
    const math::Tensor<T>& weight(igor::index_type uid) const;

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
    std::unique_ptr<const Topology> m_topology;
    std::vector<math::Tensor<T>>     m_weights;         // indexed by topology UID
    std::vector<igor::index_type>    m_execution_order;  // cached topological order
};

// ─── Free functions ──────────────────────────────────────────────────────────

/// Load probability tensors from a model_marginals text file.
/// Returns false if the file cannot be opened.
template <typename T>
bool read_parameters(const std::string& filename, RecombinationModel<T>& model);

/// Build a fully-loaded RecombinationModel in one step from files.
/// Reads the topology from model_parms, constructs the model, then
/// loads the probability tensors from model_marginals.
template <typename T = double>
RecombinationModel<T> recombination_model_from_files(
    const std::string& file_model_parms,
    const std::string& file_model_marginals);

} // namespace igor::model

#include <igor/Model/RecombinationModel.tpp>

//
// RecombinationModel.h ends here
