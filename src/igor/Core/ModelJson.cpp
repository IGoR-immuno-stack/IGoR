#include <igor/Core/ModelJson.h>

#include <igor/Core/Model_Parms.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Singleerrorrate.h>

#include <nlohmann/json.hpp>

#include <stdexcept>

namespace igor {

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

}  // namespace igor
