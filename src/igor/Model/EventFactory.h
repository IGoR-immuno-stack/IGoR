// EventFactory.h ---

#pragma once

#include <igor/Model/Export.h>
#include <igor/Core/Rec_Event.h>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace igor::model::event_factory {

using EventPtr = std::shared_ptr<Rec_Event>;

/// Builds one event from its serialized node. The creator receives the descriptor, so the
/// event is complete as soon as it exists: no default construction, no window during which
/// the object answers its capability queries with values nobody chose.
using EventCreator = std::function<EventPtr(const nlohmann::json &)>;

/**
 * @brief Register a creator under the type name the document carries.
 *
 * Keyed by the string, not by Event_type, because the string is what the file has and what
 * create() dispatches on. The enum stays what runtime code switches on; it is not the
 * serialization identity.
 */
MODEL_EXPORT void register_creator(const std::string &type_name, EventCreator func);

/**
 * @brief Build the event described by one node of the model document.
 *
 * Reads node["type"], looks the name up and hands the whole node to the creator, so the
 * caller has no dispatch of its own to write.
 *
 * @throws std::runtime_error if the node has no type, or the type is not registered. What
 *         the node's own fields have to satisfy is the concrete constructor's business.
 */
MODEL_EXPORT EventPtr create(const nlohmann::json &node);

/// Whether a creator is registered for that type name.
MODEL_EXPORT bool is_registered(const std::string &type_name);

/// Every registered type name, for diagnostics and tests.
MODEL_EXPORT std::vector<std::string> registered_type_names();

/**
 * @brief Static registration of one concrete event class.
 *
 * Usage, in the translation unit that knows the class:
 *   static Registrar<Deletion> deletion_registrar{"Deletion"};
 *
 * The class only has to offer a constructor taking a json node. Adding an event type means
 * adding a class and one line here, and no existing code changes.
 */
template <typename EventClass>
struct Registrar {
    explicit Registrar(const std::string &type_name)
    {
        register_creator(type_name, [](const nlohmann::json &node) -> EventPtr {
            return std::make_shared<EventClass>(node);
        });
    }
};

}  // namespace igor::model::event_factory

//
// EventFactory.h ends here
