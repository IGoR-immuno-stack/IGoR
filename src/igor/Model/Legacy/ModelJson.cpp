#include <igor/Model/Legacy/ModelJson.h>

#include <igor/Model/EventFactory.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Singleerrorrate.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

namespace {

/// The error rate hierarchy has no serialization of its own yet: write2txt() takes an
/// ofstream, so it cannot be captured into a string, and only Single_error_rate exposes its
/// parameter. Every model IGoR ships uses that one; the hypermutation variants live in
/// supplementary_models and carry N-mer coefficients this schema does not describe yet.
/// They raise rather than serialize a document that would rebuild a different model.
nlohmann::json error_rate_to_json(const std::shared_ptr<Error_rate> &error_rate)
{
    if (auto single = std::dynamic_pointer_cast<Single_error_rate>(error_rate)) {
        nlohmann::json out;
        out["type"] = single->type();
        out["rate"] = single->get_model_rate();
        return out;
    }
    throw std::runtime_error(
            "model_parms_to_json: error rate of type \"" + error_rate->type()
            + "\" cannot be serialized yet; only SingleErrorRate is described by schema version "
            + std::to_string(kModelJsonSchemaVersion));
}

/// The keys a model document has; an event node's are checked by its constructor.
constexpr std::initializer_list<const char *> kModelKeys = {"schema_version", "seq_type_order", "events",
                                                            "error_rate"};

const nlohmann::json &require_model_key(const nlohmann::json &doc, const char *key)
{
    const auto it = doc.find(key);
    if (it == doc.end())
        throw std::runtime_error(std::string("model json: missing key \"") + key + "\"");
    return *it;
}

/// The reader's half of error_rate_to_json(), for the one error rate the schema describes.
std::shared_ptr<Error_rate> error_rate_from_json(const nlohmann::json &node)
{
    const auto type = require_model_key(node, "type").get<std::string>();
    if (type != Single_error_rate(0.0).type()) {
        throw std::runtime_error("model json: error rate of type \"" + type
                                 + "\" is not described by schema version "
                                 + std::to_string(kModelJsonSchemaVersion)
                                 + "; the hypermutation models' JSON form is step HJ of "
                                   "docs/ERROR_MODEL_REWRITE_PLAN.md");
    }
    for (const auto &entry : node.items()) {
        if (entry.key() != "type" and entry.key() != "rate")
            throw std::runtime_error("model json: unknown key \"" + entry.key() + "\" in the error rate");
    }
    return std::make_shared<Single_error_rate>(require_model_key(node, "rate").get<double>());
}

}  // namespace

nlohmann::json model_parms_to_json(const Model_Parms &parms)
{
    nlohmann::json doc;
    doc["schema_version"] = kModelJsonSchemaVersion;

    // The segment order is model data and is not derivable from the event graph: the graph
    // says which event conditions which, not which segment sits left of which.
    doc["seq_type_order"] = parms.get_seq_type_registry().get_ordered_types();

    nlohmann::json events = nlohmann::json::array();
    for (const auto &event : parms.get_event_list()) {
        nlohmann::json node = event->to_json();

        nlohmann::json parents = nlohmann::json::array();
        for (const auto &parent : parms.get_parents(event))
            parents.push_back(parent->get_nickname());
        node["parents"] = std::move(parents);

        events.push_back(std::move(node));
    }
    doc["events"] = std::move(events);

    if (auto error_rate = const_cast<Model_Parms &>(parms).get_err_rate_p())
        doc["error_rate"] = error_rate_to_json(error_rate);

    return doc;
}


void model_parms_from_json(const nlohmann::json &doc, Model_Parms &parms)
{
    if (not parms.get_event_list().empty())
        throw std::logic_error("model_parms_from_json: the model already has events");
    for (const auto &entry : doc.items()) {
        const bool known = std::any_of(kModelKeys.begin(), kModelKeys.end(),
                                       [&](const char *key) { return entry.key() == key; });
        if (not known)
            throw std::runtime_error("model json: unknown key \"" + entry.key() + "\"");
    }
    const int version = require_model_key(doc, "schema_version").get<int>();
    if (version != kModelJsonSchemaVersion)
        throw std::runtime_error("model json: schema version " + std::to_string(version) + ", where "
                                 + std::to_string(kModelJsonSchemaVersion) + " is the one this reader knows");

    parms.set_seq_type_order(require_model_key(doc, "seq_type_order").get<std::vector<Seq_type_String>>());

    // The events first, then the edges: a parent may be listed after its child.
    const nlohmann::json &events = require_model_key(doc, "events");
    std::unordered_map<std::string, std::shared_ptr<Rec_Event>> by_nickname;
    for (const auto &node : events) {
        const std::shared_ptr<Rec_Event> event = igor::model::event_factory::create(node);
        parms.add_event(event);
        by_nickname.emplace(event->get_nickname(), event);
    }
    for (const auto &node : events) {
        const auto parents = node.find("parents");
        if (parents == node.end())
            continue;
        const auto &child = by_nickname.at(node.at("nickname").get<std::string>());
        for (const auto &parent_nickname : *parents) {
            const auto parent = by_nickname.find(parent_nickname.get<std::string>());
            if (parent == by_nickname.end())
                throw std::runtime_error("model json: event \"" + child->get_nickname()
                                         + "\" names an unknown parent \"" + parent_nickname.get<std::string>()
                                         + "\"");
            parms.add_edge(parent->second, child);
        }
    }

    parms.set_error_ratep(error_rate_from_json(require_model_key(doc, "error_rate")));
    parms.finalize();
}

} // namespace igor::model::legacy
