#include <igor/Model/Legacy/ModelJson.h>

#include <igor/Model/EventFactory.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Singleerrorrate.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <istream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
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
            + kModelJsonSchemaVersion);
}

/// The keys of the document, of its structure section and of a layout.
constexpr std::initializer_list<const char *> kModelKeys = { "schema_version", "structure", "library",
                                                             "error_rate", "engine_hints" };
constexpr std::initializer_list<const char *> kStructureKeys = { "layout", "events", "conditioning" };
constexpr std::initializer_list<const char *> kLayoutKeys = { "type", "segments" };

const nlohmann::json &require_model_key(const nlohmann::json &node, const char *key)
{
    const auto it = node.find(key);
    if (it == node.end())
        throw std::runtime_error(std::string("model json: missing key \"") + key + "\"");
    return *it;
}

void reject_unknown_model_keys(const nlohmann::json &node, std::initializer_list<const char *> allowed,
                               const std::string &where)
{
    if (not node.is_object())
        throw std::runtime_error("model json: " + where + " is not an object");
    for (const auto &entry : node.items()) {
        const bool known = std::any_of(allowed.begin(), allowed.end(),
                                       [&](const char *key) { return entry.key() == key; });
        if (not known)
            throw std::runtime_error("model json: unknown key \"" + entry.key() + "\" in " + where);
    }
}

/// The reader's half of error_rate_to_json(), for the one error rate the schema describes.
std::shared_ptr<Error_rate> error_rate_from_json(const nlohmann::json &node)
{
    const auto type = require_model_key(node, "type").get<std::string>();
    if (type != Single_error_rate(0.0).type()) {
        throw std::runtime_error("model json: error rate of type \"" + type
                                 + "\" is not described by schema version "
                                 + kModelJsonSchemaVersion
                                 + "; the hypermutation models' JSON form is step HJ of "
                                   "docs/ERROR_MODEL_REWRITE_PLAN.md");
    }
    reject_unknown_model_keys(node, { "type", "rate" }, "the error rate");
    return std::make_shared<Single_error_rate>(require_model_key(node, "rate").get<double>());
}

/// The Linear layout's segments. A DAG of segments is synthesis step 5.
std::vector<Seq_type_String> layout_segments(const nlohmann::json &structure)
{
    const nlohmann::json &layout = require_model_key(structure, "layout");
    reject_unknown_model_keys(layout, kLayoutKeys, "the layout");
    const auto type = require_model_key(layout, "type").get<std::string>();
    if (type != "Linear")
        throw std::runtime_error("model json: a layout of type \"" + type + "\"; this reader knows \"Linear\" only");

    const auto segments = require_model_key(layout, "segments").get<std::vector<Seq_type_String>>();
    if (segments.empty())
        throw std::runtime_error("model json: the layout has no segment");
    std::set<Seq_type_String> seen;
    for (const auto &segment : segments) {
        if (segment.empty() or not seen.insert(segment).second)
            throw std::runtime_error("model json: the layout lists \"" + segment + "\" "
                                     + (segment.empty() ? "as a segment" : "twice"));
    }
    return segments;
}

/// The library section: per library, per label, the template. Checked whole, used or not.
void check_library(const nlohmann::json &doc)
{
    const auto library = doc.find("library");
    if (library == doc.end())
        return;
    if (not library->is_object())
        throw std::runtime_error("model json: \"library\" is not an object");
    for (const auto &[name, templates] : library->items()) {
        if (not templates.is_object())
            throw std::runtime_error("model json: library \"" + name + "\" is not an object keyed by label");
        for (const auto &[label, entry] : templates.items()) {
            reject_unknown_model_keys(entry, { "sequence" }, "template \"" + label + "\" of library \"" + name + "\"");
            if (not require_model_key(entry, "sequence").is_string())
                throw std::runtime_error("model json: template \"" + label + "\" of library \"" + name
                                         + "\" has no sequence string");
        }
    }
}

/// The node with each realization's template put back from the library, as a gene choice's
/// constructor takes it. A node without a library is returned as it is.
nlohmann::json with_templates(const nlohmann::json &doc, const nlohmann::json &node)
{
    const auto library = node.find("library");
    if (library == node.end() or not library->is_string())
        return node;
    const auto realizations = node.find("realizations");
    if (realizations == node.end() or not realizations->is_array())
        return node;    // the constructor names what is missing

    const std::string &name = library->get_ref<const std::string &>();
    const nlohmann::json &section = require_model_key(doc, "library");
    const auto templates = section.find(name);
    if (templates == section.end())
        throw std::runtime_error("model json: the library section has no \"" + name + "\"");

    nlohmann::json resolved = node;
    for (nlohmann::json &realization : resolved.at("realizations")) {
        if (realization.contains("sequence"))
            throw std::runtime_error("model json: a template belongs in the library section, not on a realization");
        const auto label = realization.find("label");
        if (label == realization.end() or not label->is_string())
            continue;    // the constructor names what is wrong
        const auto entry = templates->find(label->get_ref<const std::string &>());
        if (entry == templates->end())
            throw std::runtime_error("model json: library \"" + name + "\" has no template \""
                                     + label->get<std::string>() + "\"");
        realization["sequence"] = entry->at("sequence");
    }
    return resolved;
}

/**
 * The order the events are built in: the legacy engine's hint if the document has one, nickname
 * order otherwise.
 *
 * The legacy engine queues events by priority, then by the event list, and the marginals file
 * is laid out by that queue. A v1 file's list is its lines; a document has no order (§14.2), so
 * model_parms_to_json() records a list that is not in nickname order as a hint: the legacy order
 * the migration's importer is to record (ENGINE_MIGRATION_PLAN Q4). Until the marginals' axes
 * are labelled (Q3), a Topology lays out its tensors by uid, so it follows the hint too.
 */
std::vector<std::string> event_order(const nlohmann::json &doc, const std::set<std::string> &nicknames)
{
    const auto hints = doc.find("engine_hints");
    if (hints == doc.end())
        return { nicknames.begin(), nicknames.end() };
    reject_unknown_model_keys(*hints, { "legacy_dfs" }, "\"engine_hints\"");
    const nlohmann::json &legacy = require_model_key(*hints, "legacy_dfs");
    reject_unknown_model_keys(legacy, { "event_order" }, "the legacy_dfs hints");

    const auto order = require_model_key(legacy, "event_order").get<std::vector<std::string>>();
    if (std::set<std::string>(order.begin(), order.end()) != nicknames or order.size() != nicknames.size())
        throw std::runtime_error("model json: the legacy_dfs event order does not list every event once");
    return order;
}

/// The (parent, child) pairs of the conditioning section, children and their parents in
/// nickname order: the order a list of parents is written in is not content.
std::vector<std::pair<std::string, std::string>> conditioning_edges(const nlohmann::json &structure,
                                                                    const std::set<std::string> &nicknames)
{
    std::vector<std::pair<std::string, std::string>> edges;
    const auto conditioning = structure.find("conditioning");
    if (conditioning == structure.end())
        return edges;
    if (not conditioning->is_object())
        throw std::runtime_error("model json: \"conditioning\" is not an object keyed by nickname");

    std::map<std::string, std::vector<std::string>> parents_of;
    for (const auto &[child, parents] : conditioning->items()) {
        if (nicknames.count(child) == 0)
            throw std::runtime_error("model json: the conditioning names an unknown event \"" + child + "\"");
        if (not parents.is_array())
            throw std::runtime_error("model json: the parents of \"" + child + "\" are not a list");
        std::vector<std::string> &sorted = parents_of[child];
        for (const auto &parent : parents) {
            const auto nickname = parent.get<std::string>();
            if (nicknames.count(nickname) == 0)
                throw std::runtime_error("model json: event \"" + child + "\" names an unknown parent \""
                                         + nickname + "\"");
            if (nickname == child)
                throw std::runtime_error("model json: event \"" + child + "\" names itself as a parent");
            sorted.push_back(nickname);
        }
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::runtime_error("model json: event \"" + child + "\" names a parent twice");
    }

    // The graph must be acyclic. Model_Parms::add_edge() checks it too, Topology::addEdge() only
    // asserts, so the check is here, once, for both readers.
    std::map<std::string, int> in_degree;
    std::map<std::string, std::vector<std::string>> children_of;
    for (const auto &nickname : nicknames)
        in_degree[nickname] = 0;
    for (const auto &[child, parents] : parents_of) {
        in_degree[child] = static_cast<int>(parents.size());
        for (const auto &parent : parents)
            children_of[parent].push_back(child);
    }
    std::vector<std::string> ready;
    for (const auto &[nickname, degree] : in_degree)
        if (degree == 0)
            ready.push_back(nickname);
    std::size_t reached = 0;
    while (not ready.empty()) {
        const std::string nickname = ready.back();
        ready.pop_back();
        ++reached;
        for (const auto &child : children_of[nickname])
            if (--in_degree[child] == 0)
                ready.push_back(child);
    }
    if (reached != nicknames.size()) {
        std::string cycle;
        for (const auto &[nickname, degree] : in_degree)
            if (degree > 0)
                cycle += (cycle.empty() ? "\"" : ", \"") + nickname + "\"";
        throw std::runtime_error("model json: the conditioning has a cycle among " + cycle);
    }

    for (const auto &[child, parents] : parents_of)
        for (const auto &parent : parents)
            edges.emplace_back(parent, child);
    return edges;
}

}  // namespace

nlohmann::json parse_model_document(std::istream &in)
{
    // One frame per object being parsed: the keys seen so far, and the key it sits under, for
    // the message. Objects inside an array are named by the array's key.
    struct Frame {
        std::set<std::string> keys;
        std::string last_key;
        std::string name;
    };
    std::vector<Frame> open;
    const nlohmann::json::parser_callback_t on_event = [&](int, nlohmann::json::parse_event_t event,
                                                           nlohmann::json &parsed) {
        switch (event) {
        case nlohmann::json::parse_event_t::object_start:
            open.push_back({ {}, {}, open.empty() ? "the document" : "\"" + open.back().last_key + "\"" });
            break;
        case nlohmann::json::parse_event_t::object_end:
            open.pop_back();
            break;
        case nlohmann::json::parse_event_t::key: {
            Frame &frame = open.back();
            frame.last_key = parsed.get<std::string>();
            if (not frame.keys.insert(frame.last_key).second)
                throw std::runtime_error("model json: key \"" + frame.last_key + "\" appears twice in "
                                         + frame.name);
            break;
        }
        default:
            break;
        }
        return true;
    };
    return nlohmann::json::parse(in, on_event, /*allow_exceptions=*/true, /*ignore_comments=*/true);
}

DocumentStructure read_document_structure(const nlohmann::json &doc)
{
    reject_unknown_model_keys(doc, kModelKeys, "the document");
    // Only the current version: a pre-release 0.x is not read once the next one exists.
    const nlohmann::json &version = require_model_key(doc, "schema_version");
    if (version != kModelJsonSchemaVersion)
        throw std::runtime_error("model json: schema version " + version.dump() + ", where \""
                                 + kModelJsonSchemaVersion + "\" is the one this reader knows");
    check_library(doc);

    const nlohmann::json &structure = require_model_key(doc, "structure");
    reject_unknown_model_keys(structure, kStructureKeys, "\"structure\"");

    DocumentStructure out;
    out.segments = layout_segments(structure);

    // Built in nickname order, or the legacy engine's, whatever order the parser keeps the
    // members in.
    const nlohmann::json &events = require_model_key(structure, "events");
    if (not events.is_object())
        throw std::runtime_error("model json: \"events\" is not an object keyed by nickname");
    std::set<std::string> nicknames;
    for (const auto &entry : events.items())
        nicknames.insert(entry.key());
    for (const auto &nickname : event_order(doc, nicknames)) {
        try {
            out.events.push_back(igor::model::event_factory::create(nickname, with_templates(doc, events.at(nickname))));
        } catch (const std::exception &e) {
            throw std::runtime_error("model json: event \"" + nickname + "\": " + e.what());
        }
    }

    out.edges = conditioning_edges(structure, nicknames);
    return out;
}

nlohmann::json model_parms_to_json(const Model_Parms &parms)
{
    nlohmann::json doc;
    doc["schema_version"] = kModelJsonSchemaVersion;

    // The segment order is model data and is not derivable from the event graph: the graph
    // says which event conditions which, not which segment sits left of which.
    const auto &segments = parms.get_seq_type_registry().get_ordered_types();
    if (segments.empty())
        throw std::runtime_error("model_parms_to_json: the model has no segment order; finalize() it first");
    nlohmann::json structure;
    structure["layout"] = { { "type", "Linear" }, { "segments", segments } };

    nlohmann::json events = nlohmann::json::object();
    nlohmann::json conditioning = nlohmann::json::object();
    nlohmann::json library = nlohmann::json::object();
    std::vector<std::string> listed;
    for (const auto &event : parms.get_event_list()) {
        listed.push_back(event->get_nickname());
        nlohmann::json node = event->to_json();

        // The templates move to the library, under the label that names them there.
        const auto library_name = node.find("library");
        if (library_name != node.end()) {
            nlohmann::json &templates = library[library_name->get<std::string>()];
            for (nlohmann::json &realization : node.at("realizations")) {
                const auto label = realization.at("label").get<std::string>();
                nlohmann::json entry = { { "sequence", realization.at("sequence") } };
                if (templates.contains(label) and templates.at(label) != entry)
                    throw std::runtime_error("model_parms_to_json: library \"" + library_name->get<std::string>()
                                             + "\" would hold two templates labelled \"" + label + "\"");
                templates[label] = std::move(entry);
                realization.erase("sequence");
            }
        }
        events[event->get_nickname()] = std::move(node);

        std::vector<std::string> parents;
        for (const auto &parent : parms.get_parents(event))
            parents.push_back(parent->get_nickname());
        std::sort(parents.begin(), parents.end());
        if (not parents.empty())
            conditioning[event->get_nickname()] = std::move(parents);
    }
    structure["events"] = std::move(events);
    structure["conditioning"] = std::move(conditioning);
    doc["structure"] = std::move(structure);
    if (not library.empty())
        doc["library"] = std::move(library);
    if (not std::is_sorted(listed.begin(), listed.end()))
        doc["engine_hints"]["legacy_dfs"]["event_order"] = std::move(listed);

    if (auto error_rate = const_cast<Model_Parms &>(parms).get_err_rate_p())
        doc["error_rate"] = error_rate_to_json(error_rate);

    return doc;
}

void model_parms_from_json(const nlohmann::json &doc, Model_Parms &parms)
{
    if (not parms.get_event_list().empty())
        throw std::logic_error("model_parms_from_json: the model already has events");

    const DocumentStructure structure = read_document_structure(doc);
    parms.set_seq_type_order(structure.segments);
    for (const auto &event : structure.events)
        parms.add_event(event);
    for (const auto &[parent, child] : structure.edges)
        parms.add_edge(parent, child);

    const auto error_rate = doc.find("error_rate");
    if (error_rate != doc.end())
        parms.set_error_ratep(error_rate_from_json(*error_rate));
    parms.finalize();
}

} // namespace igor::model::legacy
