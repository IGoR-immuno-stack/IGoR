// EventFactory.cpp ---

#include <igor/Model/EventFactory.h>

#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Insertion.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

namespace igor::model::event_factory {

namespace detail {

/// Function-local static rather than a namespace-scope object: the registrars below run
/// during static initialization, and this makes their order irrelevant.
std::map<std::string, EventCreator> &get_creators()
{
    static std::map<std::string, EventCreator> creators;
    return creators;
}

std::string known_type_names()
{
    std::string out;
    for (const auto &entry : get_creators())
        out += (out.empty() ? "" : ", ") + entry.first;
    return out;
}

}  // namespace detail

void register_creator(const std::string &type_name, EventCreator func)
{
    if (!func)
        throw std::invalid_argument(
                "EventFactory::register_creator: null creator for \"" + type_name + "\"");
    if (type_name.empty())
        throw std::invalid_argument("EventFactory::register_creator: empty type name");
    detail::get_creators()[type_name] = std::move(func);
}

EventPtr create(const nlohmann::json &node)
{
    const auto type = node.find("type");
    if (type == node.end())
        throw std::runtime_error("EventFactory: event node has no \"type\"");
    if (!type->is_string())
        throw std::runtime_error("EventFactory: event node \"type\" is not a string");

    const auto type_name = type->get<std::string>();
    const auto creator = detail::get_creators().find(type_name);
    if (creator == detail::get_creators().end())
        throw std::runtime_error("EventFactory: no creator registered for \"" + type_name
                                 + "\"; known types are " + detail::known_type_names());

    return creator->second(node);
}

bool is_registered(const std::string &type_name)
{
    return detail::get_creators().find(type_name) != detail::get_creators().end();
}

std::vector<std::string> registered_type_names()
{
    std::vector<std::string> out;
    out.reserve(detail::get_creators().size());
    for (const auto &entry : detail::get_creators())
        out.push_back(entry.first);
    return out;
}

// Registered here, in the one translation unit that may know the concrete Core classes: Core
// must not name Model. The strings are the ones the model file carries, which is why the last
// one is "DinucMarkov" and not "Dinuclmarkov".
namespace {
const Registrar<legacy::Gene_choice> gene_choice_registrar{ "GeneChoice" };
const Registrar<legacy::Deletion> deletion_registrar{ "Deletion" };
const Registrar<legacy::Insertion> insertion_registrar{ "Insertion" };
const Registrar<legacy::Dinucl_markov> dinucl_markov_registrar{ "DinucMarkov" };
}  // namespace

}  // namespace igor::model::event_factory

//
// EventFactory.cpp ends here
