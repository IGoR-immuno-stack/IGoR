#pragma once

#include <stdexcept>
#include <algorithm>
#include <cmath>

#include <igor/Model/CategoricalSamplingHandler.h>
#include <igor/Model/MarkovSamplingHandler.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Topology.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>

namespace igor::model {

// ─── Helper: Extract EventDescriptor from Rec_Event ────────────────────

inline EventDescriptor make_event_descriptor(std::shared_ptr<const Rec_Event> event,
                                             const Model_Parms& parms) {
    EventDescriptor desc;
    // Use nickname for compatibility with Model_marginals text format
    // The marginals file uses nicknames like "v_choice" not full names
    desc.name = event->get_nickname();
    desc.type = static_cast<int>(event->get_type());
    desc.gene_class = static_cast<int>(event->get_class());
    desc.side = static_cast<int>(event->get_side());

    // Build shape vector: shape interpretation depends on event type
    //
    // For CATEGORICAL events (GeneChoice, Insertion, Deletion):
    //   - Unconditional: shape = {n_realizations}           → 1D tensor
    //   - Conditional:   shape = {n_realizations, p1, p2...} → multi-D tensor
    //
    // For MARKOV events (Dinuclmarkov):
    //   - Unconditional: shape = {n_states, n_states}           → 2D tensor (square matrix)
    //   - Conditional:   shape = {n_states, n_states, p1, p2...} → multi-D tensor
    //
    // The handlers check ndim() to distinguish:
    //   - CategoricalHandler: ndim==1 for unconditional, ndim>1 for conditional
    //   - MarkovHandler:      ndim==2 for unconditional, ndim>2 for conditional

    // For Markov: second dimension is also n_states (square transition matrix)
    if (event->get_type() == Dinuclmarkov_t) {
        // Markov event: square matrix (n_states * n_states)
        int n_states = static_cast<int>(std::sqrt(static_cast<double>(event->size())));
        desc.shape.push_back(n_states); // dim 0: from state
        desc.shape.push_back(n_states); // dim 1: to state
    } else {
        // Categorical event
        desc.shape.push_back(event->size());
    }

    // Add parent dimensions for parent-conditioned events
    // Use event name to get parents (avoids const issues)
    auto parents = parms.get_parents(event->get_name());
    for (const auto& parent : parents) {
        desc.shape.push_back(parent->size());
    }

    return desc;
}

// ─── Extract EventDescriptor Vector ─────────────────────────────────────

inline std::vector<EventDescriptor> extract_event_descriptors(const Model_Parms& parms) {
    std::vector<EventDescriptor> descriptors;

    // Get events in topological order
    auto event_queue = parms.get_model_queue();

    while (!event_queue.empty()) {
        auto event = event_queue.front();
        event_queue.pop();

        descriptors.push_back(make_event_descriptor(event, parms));
    }

    return descriptors;
}

// ─── Import from Legacy Model_marginals ─────────────────────────────────

template <typename T>
void import_from_legacy(InferenceEngine<T>& engine,
                       const Model_marginals& marginals,
                       const Model_Parms& parms) {

    // Get index map to locate events in the marginal array
    auto index_map = marginals.get_index_map(parms);

    // Access the raw marginal array
    const long double* marginal_array = marginals.marginal_array_smart_p.get();

    // Iterate through each handler in the engine
    engine.for_each_handler([&](const std::string& name, MarginalHandler<T>& handler) {

        // Find the event in Model_Parms by nickname (not full name)
        auto event = parms.get_event_pointer(name, true);  // true = search by nickname
        if (!event) {
            throw std::runtime_error("Event not found in Model_Parms: " + name);
        }

        // index_map is keyed by FULL event name, not nickname
        std::string full_name = event->get_name();

        // Get base index for this event
        auto index_it = index_map.find(full_name);
        if (index_it == index_map.end()) {
            throw std::runtime_error("No index found for event: " + full_name);
        }

        int base_index = index_it->second;

        // Copy values from legacy array to handler parameters
        auto& params = handler.parameters();

        // Verify array bounds
        if (static_cast<size_t>(base_index) + params.size() > marginals.get_length()) {
             throw std::runtime_error("LegacyBridge import: Boundary check failed for event " + full_name);
        }

        for (std::size_t i = 0; i < params.size(); ++i) {
            // Convert long double to T
            params.data()[i] = static_cast<T>(marginal_array[base_index + i]);
        }
    });
}

// ─── Export to Legacy Model_marginals ───────────────────────────────────

template <typename T>
void export_to_legacy(const InferenceEngine<T>& engine,
                     Model_marginals& marginals,
                     const Model_Parms& parms) {

    // Get index map to locate events in the marginal array
    auto index_map = marginals.get_index_map(parms);

    // Access the raw marginal array
    long double* marginal_array = marginals.marginal_array_smart_p.get();

    // Iterate through each handler in the engine
    engine.for_each_handler([&](const std::string& name, const MarginalHandler<T>& handler) {

        // Find the event in Model_Parms by nickname (not full name)
        auto event = parms.get_event_pointer(name, true);  // true = search by nickname
        if (!event) {
            throw std::runtime_error("Event not found in Model_Parms: " + name);
        }

        // index_map is keyed by FULL event name, not nickname
        std::string full_name = event->get_name();

        // Get base index for this event
        auto index_it = index_map.find(full_name);
        if (index_it == index_map.end()) {
            throw std::runtime_error("No index found for event: " + full_name);
        }

        int base_index = index_it->second;

        // Copy values from handler parameters to legacy array
        const auto& params = handler.parameters();

        // Verify array bounds
        if (static_cast<size_t>(base_index) + params.size() > marginals.get_length()) {
             throw std::runtime_error("LegacyBridge export: Boundary check failed for event " + full_name);
        }

        for (std::size_t i = 0; i < params.size(); ++i) {
            // Convert T to long double
            marginal_array[base_index + i] = static_cast<long double>(params.data()[i]);
        }
    });
}

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
