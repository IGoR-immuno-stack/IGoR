#pragma once

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
    
    // Initialize handler storage
    std::size_t size = m_topology->size();
    m_handlers.resize(size);
    
    // Cache execution order
    m_execution_order = m_topology->topologicalOrder();
    
    // Instantiate factory (it doesn't need topology in constructor)
    // The factory uses static registrations.

    // Automatically instantiate empty handlers for all events in topology
    for (index_type uid = 0; uid < static_cast<index_type>(size); ++uid) {
        // Use Topology::eventName logic (helper available inside loop)
        // Hmm, factory takes nickname. 
        // We can get nickname from Topology::eventName(uid)
        std::string nickname = m_topology->eventName(uid);
        
        try {
            // factory functions are in namespace scope
            auto handler = sampling_handler_factory::create<T>(Event_type::Undefined_t, nullptr, std::vector<std::size_t>{}); 
            // TODO: Actually instantiating the right handler needs correct type/event/shape extraction
            // For now just fix the compilation issue by calling the namespace function correctly
            if (handler) {
                handler->setUid(uid);
                m_handlers[uid] = std::move(handler);
            }
        } catch (const std::exception& e) {
             // If creation fails (e.g. event structure issues), we log or rethrow?
             // Since this is the core engine init, we should probably rethrow or ensure robust handling.
             // Given it's auto-init, failing here means the Topology is unusable for Sampling.
             throw; 
        }
    }
}

template <typename T>
void SamplingEngine<T>::registerHandler(HandlerPtr handler) {
    if (!handler) {
        throw std::invalid_argument("SamplingEngine: attempt to register null handler");
    }
    
    std::string name = handler->name();
    
    // Validate that event exists in structure
    if (!m_topology->hasEvent(name)) {
        throw std::invalid_argument("SamplingEngine: unknown event in topology: " + name);
    }
    
    index_type uid = m_topology->eventId(name);
    
    // Replaces potentially auto-created handler
    
    // Handler UID must match Topology UID
    handler->setUid(uid);
    m_handlers[uid] = std::move(handler);
}

template <typename T>
const SamplingHandler<T>& SamplingEngine<T>::handler(const std::string& name) const {
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
const SamplingHandler<T>& SamplingEngine<T>::handler(index_type uid) const {
    if (uid >= static_cast<index_type>(m_handlers.size()) || uid < 0) {
        throw std::out_of_range("SamplingEngine: invalid handler UID: " + std::to_string(uid));
    }
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for UID: " + std::to_string(uid));
    }
    return *m_handlers[uid];
}

template <typename T>
SamplingHandler<T>& SamplingEngine<T>::handler(const std::string& name) {
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
SamplingHandler<T>& SamplingEngine<T>::handler(index_type uid) {
    if (uid >= static_cast<index_type>(m_handlers.size()) || uid < 0) {
        throw std::out_of_range("SamplingEngine: invalid handler UID: " + std::to_string(uid));
    }
    if (!m_handlers[uid]) {
        throw std::logic_error("SamplingEngine: no handler registered for UID: " + std::to_string(uid));
    }
    return *m_handlers[uid];
}

template <typename T>
SampledScenario SamplingEngine<T>::generate(std::mt19937_64& rng) const {
    // 1. Initialize Scenario (empty)
    SampledScenario result(m_topology->size());
    
    // 2. Iterate in topological order
    for (index_type current_uid : m_execution_order) {
        
        // Skip nodes without handlers
        if (!m_handlers[current_uid]) {
            throw std::logic_error("SamplingEngine: missing handler during generation for UID: " 
                + std::to_string(current_uid));
        }
        
        const auto& current_handler = *m_handlers[current_uid];
        std::string current_name = current_handler.name();
        
        // 3. Resolve parent realizations
        // Use Topology to get parent IDs
        const std::vector<index_type>& parent_ids = m_topology->parentsIds(current_uid);
        
        std::vector<std::size_t> parent_realizations;
        parent_realizations.reserve(parent_ids.size());
        
        for (index_type parent_uid : parent_ids) {
            
            // Look up parent handler to get its name (nickname)
            if (!m_handlers[parent_uid]) {
                 throw std::logic_error("SamplingEngine: parent handler missing for UID: " + std::to_string(parent_uid));
            }
            std::string parent_key = m_handlers[parent_uid]->name();
            
            if (result.index_of(parent_uid) >= result.events.size()) {
                throw std::logic_error("SamplingEngine: parent not realized");
            }
            parent_realizations.push_back(result.index_of(parent_uid));
        }
        
        // 4. Sample
        if (auto* markov = dynamic_cast<const MarkovSamplingHandler<T>*>(&current_handler)) {
             // For Markov, similar logic as before
             std::size_t val = current_handler.sample(rng, parent_realizations);
             result.events[current_uid] = {current_uid, {val}};
        } else {
             // Default (Categorical)
             std::size_t val = current_handler.sample(rng, parent_realizations);
             result.events[current_uid] = {current_uid, {val}};
        }
    }
    
    return result;
}

} // namespace igor::model
