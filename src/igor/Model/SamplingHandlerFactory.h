// SamplingHandlerFactory.h ---

#pragma once

#include <igor/Model/Export.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Model/SamplingHandler.h>
#include <igor/Core/Utils.h>

#include <memory>
#include <functional>
#include <string>

#include <initializer_list>

namespace igor::model::sampling_handler_factory {

using EventPtr = std::shared_ptr<Rec_Event>;

template <typename T>
using HandlerPtr = std::unique_ptr<SamplingHandler<T>>;

template <typename T>
using Creator = std::function<HandlerPtr<T>(EventPtr, const std::vector<std::size_t>&)>;

namespace detail {
template <typename T>
std::unordered_map<Event_type, Creator<T>> creators;
}

/**
 * @brief Register a creator function for an event type
 * @param type The Event_type enum value
 * @param func Function that creates a default-constructed event
 *
 * The factory uses two-phase initialization:
 * 1. Factory creates default-constructed event
 * 2. Caller populates event using setters/add_realization methods
 */
template <typename T>
void register_creator(Event_type type, Creator<T> func);

/**
 * @brief Create a default-constructed event of the specified type
 * @param type The Event_type to create
 * @return Shared pointer to default-constructed Rec_Event
 * @throws std::runtime_error if type not registered
 *
 * After creation, use event-specific methods to populate:
 * - Gene_choice: set_genomic_templates(), add_realization(name, seq)
 * - Deletion: add_realization(int)
 * - Insertion: add_realization(int)
 * - Dinucl_markov: Uses default initialization with Gene_class
 */
template <typename T> 
HandlerPtr<T> create(Event_type type, EventPtr event, const std::vector<std::size_t>& shape);

/**
 * @brief Check if an event type has a registered creator
 * @param type The Event_type to check
 * @return true if registered, false otherwise
 */
MODEL_EXPORT bool is_registered(Event_type type);

} // namespace igor::model::sampling_handler_factory

namespace igor::model {
    class Topology; // Forward declaration
}

namespace igor::model::sampling_handler_factory {

/**
 * @brief Build a list of sampling handlers based on a Topology
 * @param topology The Topology to build handlers for
 * @return A vector of default-constructed sampling handlers
 */
template <typename T>
std::vector<HandlerPtr<T>> build(const igor::model::Topology& topology);

/**
 * @brief Template for automatic event registration using static initialization
 * @tparam Type The Event_type enum value
 * @tparam EventClass The concrete event class (e.g., Gene_choice, Deletion)
 *
 * Usage:
 *   static Registrar<GeneChoice_t, Gene_choice> gene_choice_registrar;
 *
 * This will automatically register the creator function when the static
 * object is constructed during program initialization.
 */
template<typename T, typename HandlerClass>
struct Registrar {
    Registrar(std::initializer_list<Event_type> types) {
        for (auto type : types) {
            register_creator<T>(type, [](EventPtr event, const std::vector<std::size_t>& shape) { 
                return std::make_unique<HandlerClass>(event->get_nickname(), event->uid(), shape); 
            });
        }
    }
};

}

#include <igor/Model/SamplingHandlerFactory.tpp>

//
// SamplingHandlerFactory.h ends here