// RecombinationModel.tpp ---

#pragma once

#include <igor/Model/RecombinationModel.h>
// Included after RecombinationModel.h on purpose: LegacyBridge's template body needs the
// complete class, and this .tpp is itself included at the end of that header.
#include <igor/Model/LegacyBridge.h>

#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>

namespace igor::model {

// ─── Constructor ─────────────────────────────────────────────────────────────

template <typename T>
RecombinationModel<T>::RecombinationModel(std::unique_ptr<Topology> topology)
    : m_topology(std::move(topology))
{
    if (!m_topology) {
        throw std::invalid_argument("RecombinationModel: topology must not be null");
    }

    const auto n = m_topology->size();
    m_weights.reserve(n);

    for (const auto& event : *m_topology) {
        // Shape = [parent1_dims..., parent2_dims..., event_dims...]
        // Row-major: parents first, child last
        std::vector<std::size_t> shape;
        for (const auto& parent : m_topology->parents(event->uid())) {
            auto parent_shape = parent->inherent_shape();
            shape.insert(shape.end(), parent_shape.begin(), parent_shape.end());
        }
        auto event_shape = event->inherent_shape();
        shape.insert(shape.end(), event_shape.begin(), event_shape.end());

        m_weights.emplace_back(std::move(shape));   // zero-initialised Tensor<T>
    }

    m_execution_order = m_topology->topologicalOrder();
}

// ─── Ordered traversal ──────────────────────────────────────────────────────

template <typename T>
auto RecombinationModel<T>::orderedWeights() const -> OrderedList
{
    return OrderedList(m_weights, m_execution_order);
}

// ─── Weight access by UID ────────────────────────────────────────────────────

template <typename T>
math::Tensor<T>& RecombinationModel<T>::weight(igor::index_type uid)
{
    if (uid < 0 || uid >= static_cast<igor::index_type>(m_weights.size())) {
        throw std::out_of_range(
            "RecombinationModel::weight: invalid UID " + std::to_string(uid));
    }
    return m_weights[uid];
}

template <typename T>
const math::Tensor<T>& RecombinationModel<T>::weight(igor::index_type uid) const
{
    if (uid < 0 || uid >= static_cast<igor::index_type>(m_weights.size())) {
        throw std::out_of_range(
            "RecombinationModel::weight: invalid UID " + std::to_string(uid));
    }
    return m_weights[uid];
}

// ─── Weight access by name ───────────────────────────────────────────────────

template <typename T>
math::Tensor<T>& RecombinationModel<T>::weight(const std::string& name)
{
    return weight(m_topology->eventId(name));
}

template <typename T>
const math::Tensor<T>& RecombinationModel<T>::weight(const std::string& name) const
{
    return weight(m_topology->eventId(name));
}

// read_parameters() lived here: a second reader of model_marginals, kept until step 3 of the
// JSON migration. It is gone because it disagreed with Core's, and a test caught it. Core's
// txt2marginals() renormalizes after reading, on purpose, "to deal with the problem of float
// precision output from the text file" (Model_marginals.cpp), so its values are the ones
// inference actually uses. This reader copied the rounded text as-is: on mouse TCR beta, one
// v_choice entry came out 0.0440805 where Core has 0.04408051162843897. Reading marginals is
// Core's job; the tensors are filled from its flat array by LegacyBridge.

// ─── recombination_model_from_files ───────────────────────────────────────────

template <typename T>
RecombinationModel<T> recombination_model_from_files(
    const std::string& file_model_parms,
    const std::string& file_model_marginals)
{
    // 1. Core reads. It is the only tokenizer of the text format, so the v2 sections come for
    //    free: @Version, @Seq_type_order, and the six-field event lines that carry the
    //    seq_type. Both calls throw std::runtime_error when a file is missing or malformed.
    Model_Parms parms;
    parms.read_model_parms(file_model_parms);

    Model_marginals marginals(parms);
    marginals.txt2marginals(file_model_marginals, parms);

    // 2. The bridge builds: events cloned, edges translated by nickname, segment order carried.
    auto topology = import_from_legacy(parms);

    // 3. Tensors are shaped from the topology, then filled from the flat array.
    RecombinationModel<T> model(std::make_unique<Topology>(std::move(*topology)));
    import_from_legacy(model, marginals);

    return model;
}

} // namespace igor::model

//
// RecombinationModel.tpp ends here
