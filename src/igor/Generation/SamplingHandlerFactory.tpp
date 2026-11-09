#pragma once

#include <stdexcept>
#include <igor/Model/RecombinationModel.h>

namespace igor::generation::sampling_handler_factory {

template <typename T>
void register_creator(core::EventType type, Creator<T> func)
{
    if (!func) {
        throw std::invalid_argument("SamplingHandlerFactory: Null creator");
    }
    detail::get_creators<T>()[type] = func;
}

template <typename T>
HandlerPtr<T> create(core::EventType type, EventPtr event, const math::Tensor<T>& weights)
{
    auto it = detail::get_creators<T>().find(type);
    if (it == detail::get_creators<T>().end()) {
        throw std::runtime_error("SamplingHandlerFactory: No creator registered");
    }
    return it->second(event, weights);
}

template <typename T>
std::vector<HandlerPtr<T>> build(const igor::generation::RecombinationModel<T>& model)
{
    const auto& topology = model.topology();
    std::vector<HandlerPtr<T>> handlers(topology.size());

    for (igor::core::index_type uid = 0; uid < static_cast<igor::core::index_type>(topology.size()); ++uid) {
        EventPtr event = topology.event(uid);
        const auto& weights = model.weight(uid);

        auto handler = create<T>(event->get_type(), event, weights);
        handlers[uid] = std::move(handler);
    }

    return handlers;
}

}