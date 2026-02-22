#pragma once

#include <stdexcept>
#include <algorithm>
#include <cmath>

#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Topology.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>

namespace igor::model {

// ─── Convert Topology <-> Model_Parms ──────────────────────────────────────

inline std::shared_ptr<Topology> import_from_legacy(const Model_Parms& legacy_model) {
    auto topo = std::make_shared<Topology>();

    // Copy events into new Topology
    auto events = legacy_model.get_event_list();
    for (const auto& ev : events) {
        topo->addEvent(ev->copy());
    }

    // Translate edges
    for (const auto& ev : events) {
        auto parents = legacy_model.get_parents(ev->get_name());
        for (const auto& parent : parents) {
            auto parent_id = topo->eventId(parent->get_nickname());
            auto child_id = topo->eventId(ev->get_nickname());
            if (parent_id >= 0 && child_id >= 0) {
                if (!topo->hasEdge(parent_id, child_id)) {
                    topo->addEdge(parent_id, child_id);
                }
            } else {
                throw std::runtime_error("LegacyBridge: Failed to map edge from " + parent->get_nickname() + " to " + ev->get_nickname());
            }
        }
    }
    return topo;
}

inline std::shared_ptr<Model_Parms> export_to_legacy(const Topology& topology) {
    std::list<std::shared_ptr<Rec_Event>> event_list;
    for (const auto& ev : topology) {
        event_list.push_back(ev->copy());
    }
    auto parms = std::make_shared<Model_Parms>(event_list);

    for (const auto& ev : topology) {
        index_type child_id = topology.eventId(ev->get_nickname());
        auto parents = topology.parents(child_id);
        for (const auto& parent : parents) {
            parms->add_edge(parent->get_name(), ev->get_name());
        }
    }
    return parms;
}

// ─── Import Model_marginals into RecombinationModel ──────────────────────────

template <typename T>
void import_from_legacy(RecombinationModel<T>& model,
                        const Model_marginals& marginals)
{
    const auto& topology = model.topology();

    // Reconstruct a Model_Parms from the topology so we can use
    // Model_marginals::get_index_map(), which requires it.
    auto parms = export_to_legacy(topology);

    auto index_map            = marginals.get_index_map(*parms);
    const long double* source = marginals.marginal_array_smart_p.get();

    for (index_type uid = 0; uid < static_cast<index_type>(topology.size()); ++uid) {
        // index_map is keyed by full event name, not nickname
        const std::string full_name = topology.event(uid)->get_name();

        auto it = index_map.find(full_name);
        if (it == index_map.end()) {
            throw std::runtime_error(
                "import_from_legacy (RecombinationModel): no index for event '" + full_name + "'");
        }
        const int base_index = it->second;

        auto& tensor = model.weight(uid);

        if (static_cast<std::size_t>(base_index) + tensor.size() > marginals.get_length()) {
            throw std::runtime_error(
                "import_from_legacy (RecombinationModel): boundary check failed for '" + full_name + "'");
        }

        T* dest = tensor.data();
        for (std::size_t i = 0; i < tensor.size(); ++i)
            dest[i] = static_cast<T>(source[base_index + i]);
    }
}

} // namespace igor::model
