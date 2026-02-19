// EventFactory.cpp ---

#include <igor/Model/EventFactory.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <unordered_map>
#include <stdexcept>
#include <string>

namespace igor::model::event_factory {

namespace detail {

// Registry mapping Event_type to event creator functions
std::unordered_map<Event_type, EventCreator> creators;

// Helper to get event type name for error messages
std::string event_type_name(Event_type type) {
    switch (type) {
        case GeneChoice_t: return "GeneChoice";
        case Deletion_t: return "Deletion";
        case Insertion_t: return "Insertion";
        case Dinuclmarkov_t: return "Dinuclmarkov";
        default: return "Unknown(" + std::to_string(static_cast<int>(type)) + ")";
    }
}

}

void register_creator(Event_type type, EventCreator func)
{
    if (!func) {
        throw std::invalid_argument(
            "EventFactory::register_creator: Cannot register null creator for " +
            detail::event_type_name(type));
    }
    detail::creators[type] = func;
}

EventPtr create(Event_type type)
{
    auto it = detail::creators.find(type);
    if (it == detail::creators.end()) {
        throw std::runtime_error(
            "EventFactory: No creator registered for type: " +
            detail::event_type_name(type));
    }

    return it->second();
}

bool is_registered(Event_type type)
{
    return detail::creators.find(type) != detail::creators.end();
}

// Register Core events within this translation unit to handle dependency cycle
// (This creates a Model->Core dependency in this file, but avoids Core->Model headers in Core files)
namespace {
static Registrar<GeneChoice_t, Gene_choice> gene_choice_registrar;
static Registrar<Deletion_t, Deletion> deletion_registrar;
static Registrar<Insertion_t, Insertion> insertion_registrar;
static Registrar<Dinuclmarkov_t, Dinucl_markov> dinucl_markov_registrar;
}

} // namespace igor::model::event_factory

//
// EventFactory.cpp ends here