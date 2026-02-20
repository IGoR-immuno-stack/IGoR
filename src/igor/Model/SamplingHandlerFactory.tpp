#pragma once

#include <stdexcept>
#include <igor/Model/Topology.h>

namespace igor::model::sampling_handler_factory {

// Instantiated in header, nothing to declare here in detail

template <typename T>
void register_creator(Event_type type, Creator<T> func) 
{
    if (!func) {
        throw std::invalid_argument("SamplingHandlerFactory: Null creator");
    }
    detail::creators<T>[type] = func;
}

template <typename T> 
HandlerPtr<T> create(Event_type type, EventPtr event, const std::vector<std::size_t>& shape) 
{
    auto it = detail::creators<T>.find(type);
    if (it == detail::creators<T>.end()) {
        throw std::runtime_error("SamplingHandlerFactory: No creator registered");
    }
    return it->second(event, shape);
}

template <typename T>
std::vector<HandlerPtr<T>> build(const igor::model::Topology& topology) 
{
    std::vector<HandlerPtr<T>> handlers(topology.size());

    for (igor::index_type uid = 0; uid < static_cast<igor::index_type>(topology.size()); ++uid) {
        EventPtr event = topology.event(uid);
        
        // 1. Compute Shape: [current_shape..., parent1_shape..., parent2_shape...]
        std::vector<std::size_t> shape = event->inherent_shape();
        
        // Iterate over the parent events (as defined by Topology) to append conditional dimensions
        for (igor::index_type parent_uid : topology.parentsIds(uid)) {
            EventPtr parent = topology.event(parent_uid);
            std::vector<std::size_t> parent_shape = parent->inherent_shape();
            shape.insert(shape.end(), parent_shape.begin(), parent_shape.end());
        }

        // 2. Delegate to the core factory function
        auto handler = create<T>(event->get_type(), event, shape);

        // 3. Assign the Topology UID to the allocated handler
        handler->setUid(uid);

        // 4. Store in the exact same topology order
        handlers[uid] = std::move(handler);
    }

    return handlers;
}

}