#pragma once

#include <igor/Model/SamplingEngine.h>
#include <igor/Model/SamplingHandlerFactory.h>
#include <igor/Model/RecombinationModel.h>

#include <stdexcept>
#include <string>
#include <vector>


namespace igor::model {

template <typename T>
SamplingEngine<T>::SamplingEngine(std::shared_ptr<const RecombinationModel<T>> model)
    : m_model(std::move(model))
{
    if (!m_model) {
        throw std::invalid_argument("SamplingEngine: model must not be null");
    }

    // Build handlers — each one borrows a const reference to its tensor in the model
    m_handlers = sampling_handler_factory::build<T>(*m_model);
    m_execution_order = m_model->topology().topologicalOrder();

    // Precompute CDFs immediately — the model's tensors are already filled
    for (auto& h : m_handlers) {
        if (h) h->precomputeCDF();
    }
}

template <typename T>
const SamplingHandler<T>& SamplingEngine<T>::handler(const std::string& name) const
{
    const auto& topo = m_model->topology();
    if (!topo.hasEvent(name)) {
        throw std::out_of_range("SamplingEngine: event not found: " + name);
    }
    index_type uid = topo.eventId(name);
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
    const auto& topo = m_model->topology();
    if (!topo.hasEvent(name)) {
        throw std::out_of_range("SamplingEngine: event not found: " + name);
    }
    index_type uid = topo.eventId(name);
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
    SampledScenario result(m_model->topology().size());

    // 2. Iterate in topological order
    for (const auto& handler_ptr : this->orderedHandlers()) {

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
    return Adjacency_t(m_handlers, m_model->topology().parentsIds(uid));
}

template <typename T>
auto SamplingEngine<T>::children(index_type uid) const -> Adjacency_t
{
    return Adjacency_t(m_handlers, m_model->topology().childrenIds(uid));
}

} // namespace igor::model
