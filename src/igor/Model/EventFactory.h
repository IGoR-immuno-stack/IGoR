// EventFactory.h ---

#pragma once

#include <igor/Model/Export.h>
#include <igor/Model/Legacy/Rec_Event.h>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace igor::model::event_factory {

using EventPtr = std::shared_ptr<legacy::Rec_Event>;

/// Builds one event from its nickname, its serialized node and the document's library. The creator
/// receives the descriptor, so the event is complete as soon as it exists: no default
/// construction, no window during which the object answers its capability queries with values
/// nobody chose.
using EventCreator =
        std::function<EventPtr(const std::string &nickname, const nlohmann::json &node, const nlohmann::json &library)>;

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
 * The nickname is the node's key in the document's "events" object, so it comes beside the
 * node rather than in it. \a library is the document's "library" section, null if it has none: a
 * gene choice takes its templates from it. Reads node["type"], looks the name up and hands all
 * three to the creator, so the caller has no dispatch of its own to write.
 *
 * @throws std::runtime_error if the nickname is empty, the node has no type, or the type is not
 *         registered. What the node's own fields have to satisfy is the concrete constructor's
 *         business: it checks the node against its type's schema (ModelSchema.h).
 */
MODEL_EXPORT EventPtr create(const std::string &nickname, const nlohmann::json &node, const nlohmann::json &library);

/// Whether a creator is registered for that type name.
MODEL_EXPORT bool is_registered(const std::string &type_name);

/// Every registered type name, for diagnostics and tests.
MODEL_EXPORT std::vector<std::string> registered_type_names();

/**
 * @brief Whether an event class draws templates from the document's library, its JSON constructor
 * then taking the library after the node. Opted into by specialization, next to the registration:
 * nlohmann::json converts to nearly anything, so asking whether the three-argument constructor
 * exists would find a legacy one that takes a pair of ints.
 */
template <typename EventClass>
inline constexpr bool takes_library = false;

/**
 * @brief Static registration of one concrete event class.
 *
 * Usage, in the translation unit that knows the class:
 *   static Registrar<Deletion> deletion_registrar{"Deletion"};
 *
 * The class only has to offer a constructor taking a nickname and a json node, and the library
 * too if it draws templates from it (takes_library), and a schema,
 * schemas/<version>/event-<type>.schema.json, which the event union lists. Adding an event type
 * means adding a class, a schema and one line here, and no existing code changes.
 */
template <typename EventClass>
struct Registrar {
    explicit Registrar(const std::string &type_name)
    {
        register_creator(type_name, [](const std::string &nickname, const nlohmann::json &node,
                                       const nlohmann::json &library) -> EventPtr {
            if constexpr (takes_library<EventClass>)
                return std::make_shared<EventClass>(nickname, node, library);
            else
                return std::make_shared<EventClass>(nickname, node);
        });
    }
};

}  // namespace igor::model::event_factory

//
// EventFactory.h ends here
