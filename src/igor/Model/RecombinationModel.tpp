// RecombinationModel.tpp ---

#pragma once

#include <igor/Model/RecombinationModel.h>

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

// ─── read_parameters ─────────────────────────────────────────────────────────
//
// File format (model_marginals text):
//   @nickname        – start of an event block
//   $Dim[d1,d2,...]  – ignored (dimension annotation)
//   #[parent,idx]    – ignored (conditioning context header)
//   %v1,v2,v3,...    – probability values
//
// All %-lines for a given event, concatenated in file order, form the exact
// flat memory layout of the tensor:
//   [slice_0_probs..., slice_1_probs..., ...]

template <typename T>
bool read_parameters(const std::string& filename, RecombinationModel<T>& model)
{
    std::ifstream infile(filename);
    if (!infile) return false;

    // ── Pass 1: collect flat value vectors per event nickname ─────────────
    std::unordered_map<std::string, std::vector<T>> values_by_event;
    std::string current_event;
    std::string line;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        switch (line[0]) {
        case '@':
            current_event = line.substr(1);
            values_by_event.emplace(current_event, std::vector<T>{});
            break;

        case '%': {
            if (current_event.empty()) break;
            auto& vec = values_by_event[current_event];
            std::size_t pos = 1;
            while (pos < line.size()) {
                std::size_t comma = line.find(',', pos);
                std::string token = (comma == std::string::npos)
                    ? line.substr(pos)
                    : line.substr(pos, comma - pos);
                if (!token.empty())
                    vec.push_back(static_cast<T>(std::stod(token)));
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
            break;
        }

        default: break;
        }
    }

    // ── Pass 2: fill tensors ─────────────────────────────────────────────
    const auto& topo = model.topology();

    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(topo.size()); ++uid)
    {
        const std::string& nickname = topo.event(uid)->get_nickname();
        auto& tensor = model.weight(uid);

        auto it = values_by_event.find(nickname);
        if (it == values_by_event.end()) {
            throw std::runtime_error(
                "read_parameters: event '" + nickname +
                "' not found in file '" + filename + "'");
        }

        const auto& vals = it->second;
        if (vals.size() != tensor.size()) {
            throw std::runtime_error(
                "read_parameters: size mismatch for '" + nickname +
                "': file has " + std::to_string(vals.size()) +
                " values, model expects " + std::to_string(tensor.size()));
        }

        std::copy(vals.begin(), vals.end(), tensor.data());
    }

    return true;
}

// ─── recombination_model_from_files ───────────────────────────────────────────

template <typename T>
RecombinationModel<T> recombination_model_from_files(
    const std::string& file_model_parms,
    const std::string& file_model_marginals)
{
    // 1. Build topology from model_parms
    auto topology = read_topology(file_model_parms);
    if (!topology) {
        throw std::runtime_error(
            "recombination_model_from_files: failed to read topology from '"
            + file_model_parms + "'");
    }

    // 2. Construct model (move shared_ptr content into unique_ptr)
    RecombinationModel<T> model(
        std::make_unique<Topology>(std::move(*topology)));

    // 3. Load marginals
    if (!read_parameters(file_model_marginals, model)) {
        throw std::runtime_error(
            "recombination_model_from_files: failed to read parameters from '"
            + file_model_marginals + "'");
    }

    return model;
}

} // namespace igor::model

//
// RecombinationModel.tpp ends here
