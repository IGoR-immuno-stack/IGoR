#pragma once

#include <igor/Model/SamplingEngine.h>
#include <igor/Model/SamplingHandlerFactory.h>

#include <stdexcept>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>


namespace igor::model {

template <typename T>
SamplingEngine<T>::SamplingEngine(std::shared_ptr<const Topology> topology)
    : m_topology(std::move(topology))
{
    if (!m_topology) {
        throw std::invalid_argument("SamplingEngine: topology must not be null");
    }

    // Use the factory to automatically build the handlers with the correct shapes
    m_handlers = sampling_handler_factory::build<T>(*m_topology);
    m_execution_order = m_topology->topologicalOrder();
}

template <typename T>
const SamplingHandler<T>& SamplingEngine<T>::handler(const std::string& name) const
{
    if (!m_topology->hasEvent(name)) {
        throw std::out_of_range("SamplingEngine: event not found: " + name);
    }
    index_type uid = m_topology->eventId(name);
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for: " + name);
    }
    return *m_handlers[uid];
}

template <typename T>
const SamplingHandler<T>& SamplingEngine<T>::handler(index_type uid) const
{
    if (uid >= static_cast<index_type>(m_handlers.size()) || uid < 0) {
        throw std::out_of_range("SamplingEngine: invalid handler UID: " + std::to_string(uid));
    }
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for UID: " + std::to_string(uid));
    }
    return *m_handlers[uid];
}

template <typename T>
SamplingHandler<T>& SamplingEngine<T>::handler(const std::string& name)
{
    if (!m_topology->hasEvent(name)) {
        throw std::out_of_range("SamplingEngine: event not found: " + name);
    }
    index_type uid = m_topology->eventId(name);
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for: " + name);
    }
    return *m_handlers[uid];
}

template <typename T>
SamplingHandler<T>& SamplingEngine<T>::handler(index_type uid)
{
    if (uid >= static_cast<index_type>(m_handlers.size()) || uid < 0) {
        throw std::out_of_range("SamplingEngine: invalid handler UID: " + std::to_string(uid));
    }
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for UID: " + std::to_string(uid));
    }
    return *m_handlers[uid];
}

template <typename T>
SampledScenario SamplingEngine<T>::run(std::mt19937_64& rng) const
{
    // 1. Initialize Scenario (empty)
    SampledScenario result(m_topology->size());

    // 2. Iterate in topological order
    for (const auto& handler_ptr : this->orderedHandlers()) {

        // Skip nodes without handlers (though they shouldn't exist in orderedHandlers if topological order is valid, but checking here)
        if (!handler_ptr) {
            throw std::logic_error("SamplingEngine: missing handler during generation");
        }

        const auto& current_handler = *handler_ptr;
        index_type current_uid = current_handler.uid();

        // 3. Resolve parent realizations
        auto parents_nav = this->parents(current_uid);

        std::vector<std::size_t> parent_realizations;
        parent_realizations.reserve(parents_nav.size());

        for (const auto& parent_ptr : parents_nav) {

            if (!parent_ptr) {
                 throw std::logic_error("SamplingEngine: parent handler missing");
            }
            index_type parent_uid = parent_ptr->uid();

            if (parent_uid >= static_cast<index_type>(result.events.size()) || result.events[parent_uid].indices.empty()) {
                throw std::logic_error("SamplingEngine: parent not realized");
            }
            parent_realizations.push_back(result.index_of(parent_uid));
        }

        // 4. Sample
        std::size_t val = current_handler.sample(rng, parent_realizations);
        result.events[current_uid] = {current_uid, {val}};
    }

    return result;
}

template <typename T>
auto SamplingEngine<T>::orderedHandlers() const -> OrderedList
{
    return OrderedList(m_handlers, m_execution_order);
}

template <typename T>
auto SamplingEngine<T>::parents(index_type uid) const -> Adjacency_t
{
    return Adjacency_t(m_handlers, m_topology->parentsIds(uid));
}

template <typename T>
auto SamplingEngine<T>::children(index_type uid) const -> Adjacency_t
{
    return Adjacency_t(m_handlers, m_topology->childrenIds(uid));
}

// ─── read_parameters ─────────────────────────────────────────────────────────
//
// File format (model_marginals text):
//   @nickname        – start of an event block
//   $Dim[d1,d2,...]  – ignored (dimension annotation)
//   #[parent,idx]    – ignored (conditioning context header)
//   %v1,v2,v3,...    – probability values
//
// All %‑lines for a given event, concatenated in file order, form the exact
// flat memory layout of the SamplingHandler's parameter tensor:
//   [slice_0_probs..., slice_1_probs..., ...]   (m_slice_count × m_realization_count)
// where slice ordering matches the parent-conditioning order in the topology.

template <typename T>
bool read_parameters(const std::string& filename, SamplingEngine<T>& engine)
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
            // Parse comma-separated doubles starting after the '%'
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

        default: break;  // '$', '#', or blank annotation lines
        }
    }

    // ── Pass 2: fill handlers ─────────────────────────────────────────────
    for (auto& handler_ptr : engine) {
        if (!handler_ptr) continue;
        SamplingHandler<T>& h = *handler_ptr;

        auto it = values_by_event.find(h.name());
        if (it == values_by_event.end()) {
            throw std::runtime_error(
                "read_parameters: event '" + h.name() + "' not found in file '" + filename + "'");
        }

        const auto& vals = it->second;
        if (vals.size() != h.rawDataSize()) {
            throw std::runtime_error(
                "read_parameters: size mismatch for '" + h.name() +
                "': file has " + std::to_string(vals.size()) +
                " values, handler expects " + std::to_string(h.rawDataSize()));
        }

        std::copy(vals.begin(), vals.end(), h.rawData());
        h.precomputeCDF();
    }

    return true;
}

} // namespace igor::model
