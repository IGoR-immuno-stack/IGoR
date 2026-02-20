#pragma once

#include "igor/Model/SamplingEngine.h"
#include <stdexcept>
#include <vector>
#include <string>

#include <igor/Model/SamplingHandlerFactory.h>
#include <igor/Model/MarkovSamplingHandler.h>

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

template <typename T>
bool readParameters(const std::string& filename, SamplingEngine<T>& engine)
{
    const Topology& topo = engine.topology();

    // 1. Reconstruct a Model_Parms from the current topology (needed by Model_marginals)
    auto parms = export_to_legacy(topo);

    // 2. Allocate and read the marginals file
    Model_marginals marginals(*parms);
    try {
        marginals.txt2marginals(filename, *parms);
    } catch (const std::exception&) {
        return false;
    }

    // 3. Build the flat-array index map (keyed by full event name)
    auto index_map            = marginals.get_index_map(*parms);
    const long double* source = marginals.marginal_array_smart_p.get();

    // 4. Write values into each handler and rebuild its CDF
    for (index_type uid = 0;
         uid < static_cast<index_type>(topo.size());
         ++uid)
    {
        auto event = topo.event(uid);

        // index_map uses the full event name, not the nickname
        std::string full_name = event->get_name();

        auto index_it = index_map.find(full_name);
        if (index_it == index_map.end()) {
            throw std::runtime_error(
                "readParameters: no index found for event '" + full_name + "'");
        }
        const int base_index = index_it->second;

        SamplingHandler<T>& h = engine.handler(uid);

        // Resolve the concrete parameter tensor via downcast
        math::Tensor<T>* params = nullptr;
        if (auto* cat = dynamic_cast<CategoricalSamplingHandler<T>*>(&h)) {
            params = &cat->parameters();
        } else if (auto* mrk = dynamic_cast<MarkovSamplingHandler<T>*>(&h)) {
            params = &mrk->parameters();
        } else {
            throw std::runtime_error(
                "readParameters: unknown handler type for event '" + full_name + "'");
        }

        // Bounds check
        if (static_cast<std::size_t>(base_index) + params->size() > marginals.get_length()) {
            throw std::runtime_error(
                "readParameters: boundary check failed for event '" + full_name + "'");
        }

        // Copy long double → T
        for (std::size_t i = 0; i < params->size(); ++i) {
            params->data()[i] = static_cast<T>(source[base_index + i]);
        }

        // Invalidate and rebuild the cached CDF slices
        h.precomputeCDF();
    }

    return true;
}

} // namespace igor::model
