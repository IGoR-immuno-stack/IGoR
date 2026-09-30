// EventFactory.h ---

#pragma once

#include <igor/Model/Export.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Model/SamplingHandler.h>
#include <igor/Core/Utils.h>

#include <memory>
#include <functional>

namespace igor::model::event_factory {

using EventPtr = std::shared_ptr<Rec_Event>;
using EventCreator = std::function<EventPtr()>;

/**
 * @brief Register a creator function for an event type
 * @param type The Event_type enum value
 * @param func Function that creates a default-constructed event
 *
 * The factory uses two-phase initialization:
 * 1. Factory creates default-constructed event
 * 2. Caller populates event using setters/add_realization methods
 */
MODEL_EXPORT void register_creator(Event_type type, EventCreator func);

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
MODEL_EXPORT EventPtr create(Event_type type);

/**
 * @brief Check if an event type has a registered creator
 * @param type The Event_type to check
 * @return true if registered, false otherwise
 */
MODEL_EXPORT bool is_registered(Event_type type);

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
template<Event_type Type, typename EventClass>
struct Registrar {
    Registrar() {
        register_creator(Type, []() { return std::make_shared<EventClass>(); });
    }
};

}

//
// EventFactory.h ends here