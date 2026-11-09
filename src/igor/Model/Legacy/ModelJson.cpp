#include <igor/Model/Legacy/ModelJson.h>

#include <igor/Model/EventFactory.h>
#include <igor/Model/ModelSchema.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Singleerrorrate.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
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

/**
 * The checks every reader makes first: the version, then the schema, with the sections the reader
 * needs. The version comes first, and in its own words: another version's document would break
 * this version's schema in ways that say nothing of the real cause.
 */
void check_document(const nlohmann::json &doc, std::initializer_list<const char *> required)
{
    const auto version = doc.is_object() ? doc.find("schema_version") : doc.end();
    if (version != doc.end() and *version != kModelJsonSchemaVersion)
        throw std::runtime_error("model json: schema version " + version->dump() + ", where \""
                                 + kModelJsonSchemaVersion + "\" is the one this reader knows");
    igor::model::schema::validate(doc, "model.schema.json", "the document", required);
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
    if (hints == doc.end() or not hints->contains("legacy_dfs"))
        return { nicknames.begin(), nicknames.end() };
    const auto order = hints->at("legacy_dfs").at("event_order").get<std::vector<std::string>>();
    if (std::set<std::string>(order.begin(), order.end()) != nicknames)
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

    std::map<std::string, std::vector<std::string>> parents_of;
    for (const auto &[child, parents] : conditioning->items()) {
        if (nicknames.count(child) == 0)
            throw std::runtime_error("model json: the conditioning names an unknown event \"" + child + "\"");
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

/// An event's realizations as (index, label), in index order.
std::vector<std::pair<int, nlohmann::json>> realizations_by_index(const Rec_Event &event)
{
    std::vector<std::pair<int, nlohmann::json>> out;
    for (const auto &entry : event.get_realizations_map())
        out.emplace_back(entry.second.index, realization_label(entry.second));
    std::sort(out.begin(), out.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    return out;
}

/// One axis of a parameters node, as the legacy array lays it out: its event, its stride in the
/// event's block, and the realization index at each position of the node's values along it.
struct LegacyAxis {
    const Rec_Event *event;
    std::size_t stride;
    std::vector<int> index_at;
};

using InverseOffsets =
        std::unordered_map<Rec_Event_name, std::list<std::pair<std::shared_ptr<const Rec_Event>, int>>>;

/// An event's parents with their strides, outermost first: the order of the text file's $Dim line.
std::list<std::pair<std::shared_ptr<const Rec_Event>, int>> parents_by_stride(const InverseOffsets &inverse_offsets,
                                                                            const Rec_Event &event)
{
    const auto it = inverse_offsets.find(event.get_name());
    if (it == inverse_offsets.end())
        return {};
    auto parents = it->second;
    parents.sort(offset_comp());
    return parents;
}

/// The strides of an event's own axes: one axis, or a Markov chain's two, from then to, the from
/// axis stepping over a row of the four nucleotides.
std::vector<std::size_t> own_strides(const Rec_Event &event)
{
    if (event.get_type() == Dinuclmarkov_t)
        return { event.get_realizations_map().size(), 1 };
    return { 1 };
}

nlohmann::json nested_values(const long double *block, const std::vector<LegacyAxis> &axes, std::size_t depth,
                             std::size_t offset)
{
    nlohmann::json out = nlohmann::json::array();
    for (const int index : axes[depth].index_at) {
        const std::size_t at = offset + axes[depth].stride * static_cast<std::size_t>(index);
        if (depth + 1 == axes.size())
            out.push_back(static_cast<double>(block[at]));
        else
            out.push_back(nested_values(block, axes, depth + 1, at));
    }
    return out;
}

/// A node's axes, checked against its event: its parents in any order, then its own axes. Each
/// axis's labels must be its event's, each once, in any order. An axis without labels follows its
/// event's index order, as the model reading it holds it.
std::vector<LegacyAxis> read_axes(const nlohmann::json &node, const Rec_Event &event,
                                  const std::list<std::pair<std::shared_ptr<const Rec_Event>, int>> &parents)
{
    const nlohmann::json &axes = node.at("axes");
    const std::vector<std::size_t> own = own_strides(event);
    if (axes.size() != parents.size() + own.size())
        throw std::runtime_error("the axes must be the event's " + std::to_string(parents.size())
                                 + " parents, then its " + std::to_string(own.size()) + " own");

    std::vector<LegacyAxis> out;
    std::set<std::string> named_parents;
    for (std::size_t k = 0; k != axes.size(); ++k) {
        const std::string where = "axis " + std::to_string(k);
        const auto name = axes[k].at("event").get<std::string>();
        LegacyAxis axis{ &event, 0, {} };
        if (k >= parents.size()) {
            if (name != event.get_nickname())
                throw std::runtime_error(where + " names \"" + name + "\", where the event's own axis comes");
            axis.stride = own[k - parents.size()];
        } else {
            if (name == event.get_nickname())
                throw std::runtime_error(where + " names the event itself, whose own axis comes after its parents'");
            const auto parent = std::find_if(parents.begin(), parents.end(),
                                             [&](const auto &entry) { return entry.first->get_nickname() == name; });
            if (parent == parents.end())
                throw std::runtime_error(where + " names \"" + name + "\", which is not a parent of the event");
            if (not named_parents.insert(name).second)
                throw std::runtime_error("two axes name the parent \"" + name + "\"");
            axis.event = parent->first.get();
            axis.stride = static_cast<std::size_t>(parent->second);
        }

        const auto domain = realizations_by_index(*axis.event);
        const auto labels = axes[k].find("labels");
        if (labels == axes[k].end()) {
            for (const auto &entry : domain)
                axis.index_at.push_back(entry.first);
            out.push_back(std::move(axis));
            continue;
        }
        std::map<std::string, int> index_of;
        for (const auto &[index, label] : domain)
            index_of.emplace(label.dump(), index);
        if (labels->size() != index_of.size())
            throw std::runtime_error(where + " must list the " + std::to_string(index_of.size()) + " labels of \""
                                     + name + "\"");
        std::set<int> placed;
        for (const auto &label : *labels) {
            const auto found = index_of.find(label.dump());
            if (found == index_of.end())
                throw std::runtime_error(where + ": \"" + name + "\" has no realization labelled " + label.dump());
            if (not placed.insert(found->second).second)
                throw std::runtime_error(where + " lists " + label.dump() + " twice");
            axis.index_at.push_back(found->second);
        }
        out.push_back(std::move(axis));
    }
    return out;
}

/// The node's values into the event's block, each at the realizations its labels name. The
/// schema has them non-negative numbers; a document built in code may still hold an infinity.
void place_values(const nlohmann::json &values, long double *block, const std::vector<LegacyAxis> &axes,
                  std::size_t depth, std::size_t offset)
{
    const LegacyAxis &axis = axes[depth];
    if (not values.is_array() or values.size() != axis.index_at.size())
        throw std::runtime_error("the values along axis " + std::to_string(depth) + " are not a list of "
                                 + std::to_string(axis.index_at.size()));
    for (std::size_t position = 0; position != axis.index_at.size(); ++position) {
        const std::size_t at = offset + axis.stride * static_cast<std::size_t>(axis.index_at[position]);
        if (depth + 1 != axes.size()) {
            place_values(values[position], block, axes, depth + 1, at);
            continue;
        }
        const nlohmann::json &value = values[position];
        if (not value.is_number() or not std::isfinite(value.get<double>()))
            throw std::runtime_error("value " + value.dump() + " is not a finite number");
        block[at] = value.get<double>();
    }
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

bool starts_like_json_document(std::istream &in)
{
    in >> std::ws;
    const auto first = in.peek();
    return first == '{' or first == '/';
}

bool model_file_has_parameters(const std::string &filename)
{
    std::ifstream in(filename);
    if (not in or not starts_like_json_document(in))
        return false;
    return parse_model_document(in).contains("parameters");
}

DocumentStructure read_document_structure(const nlohmann::json &doc)
{
    check_document(doc, { "structure" });
    const nlohmann::json &structure = doc.at("structure");

    DocumentStructure out;
    out.segments = structure.at("layout").at("segments").get<std::vector<Seq_type_String>>();

    // Built in nickname order, or the legacy engine's, whatever order the parser keeps the
    // members in. A gene choice takes its templates from the library.
    static const nlohmann::json no_library;
    const nlohmann::json &library = doc.contains("library") ? doc.at("library") : no_library;
    const nlohmann::json &events = structure.at("events");
    std::set<std::string> nicknames;
    for (const auto &entry : events.items())
        nicknames.insert(entry.key());
    for (const auto &nickname : event_order(doc, nicknames)) {
        try {
            out.events.push_back(igor::model::event_factory::create(nickname, events.at(nickname), library));
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
    doc["$schema"] = igor::model::schema::document_schema_uri();
    doc["schema_version"] = kModelJsonSchemaVersion;

    // The segment order is model data and is not derivable from the event graph: the graph
    // says which event conditions which, not which segment sits left of which.
    const auto &segments = parms.get_seq_type_registry().get_ordered_types();
    if (segments.empty())
        throw std::runtime_error("model_parms_to_json: the model has no segment order; finalize() it first");
    nlohmann::json structure;
    structure["layout"]["type"] = "Linear";
    structure["layout"]["segments"] = segments;

    nlohmann::json events = nlohmann::json::object();
    nlohmann::json conditioning = nlohmann::json::object();
    nlohmann::json library = nlohmann::json::object();
    std::vector<std::string> listed;
    for (const auto &event : parms.get_event_list()) {
        listed.push_back(event->get_nickname());
        nlohmann::json node = event->to_json();

        // A templated event's templates go to the library it names, under their labels.
        const auto library_name = node.find("library");
        if (library_name != node.end()) {
            nlohmann::json &templates = library[library_name->get<std::string>()];
            for (const auto &entry : event->get_realizations_map()) {
                const Event_realization &realization = entry.second;
                nlohmann::json value;
                value["sequence"] = realization.value_str;
                if (templates.contains(realization.name) and templates.at(realization.name) != value)
                    throw std::runtime_error("model_parms_to_json: library \"" + library_name->get<std::string>()
                                             + "\" would hold two templates labelled \"" + realization.name + "\"");
                templates[realization.name] = std::move(value);
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

    // The schema allows SingleErrorRate only: the hypermutation models' JSON form is step HJ of
    // docs/ERROR_MODEL_REWRITE_PLAN.md.
    const auto error_rate = doc.find("error_rate");
    if (error_rate != doc.end())
        parms.set_error_ratep(std::make_shared<Single_error_rate>(error_rate->at("rate").get<double>()));
    parms.finalize();
}

nlohmann::json model_marginals_to_json(const Model_Parms &parms, const Model_marginals &marginals)
{
    const auto queue = parms.get_model_queue();
    const auto index_map = marginals.get_index_map(parms, queue);
    const InverseOffsets inverse_offsets = marginals.get_inverse_offset_map(parms, queue);

    nlohmann::json nodes = nlohmann::json::object();
    for (const auto &event : parms.get_event_list()) {
        std::vector<LegacyAxis> axes;
        nlohmann::json axes_json = nlohmann::json::array();
        const auto add_axis = [&](const Rec_Event &axis_event, std::size_t stride) {
            LegacyAxis axis{ &axis_event, stride, {} };
            nlohmann::json labels = nlohmann::json::array();
            for (auto &[index, label] : realizations_by_index(axis_event)) {
                axis.index_at.push_back(index);
                labels.push_back(std::move(label));
            }
            nlohmann::json entry;
            entry["event"] = axis_event.get_nickname();
            entry["labels"] = std::move(labels);
            axes_json.push_back(std::move(entry));
            axes.push_back(std::move(axis));
        };
        for (const auto &[parent, stride] : parents_by_stride(inverse_offsets, *event))
            add_axis(*parent, static_cast<std::size_t>(stride));
        for (const std::size_t stride : own_strides(*event))
            add_axis(*event, stride);

        const long double *block = marginals.marginal_array_smart_p.get() + index_map.at(event->get_name());
        nlohmann::json node;
        node["axes"] = std::move(axes_json);
        node["values"] = nested_values(block, axes, 0, 0);
        nodes[event->get_nickname()] = std::move(node);
    }
    nlohmann::json section;
    section["nodes"] = std::move(nodes);
    return section;
}

void model_marginals_from_json(const nlohmann::json &doc, const Model_Parms &parms, Model_marginals &marginals)
{
    check_document(doc, { "parameters" });
    const nlohmann::json &nodes = doc.at("parameters").at("nodes");

    const auto queue = parms.get_model_queue();
    const auto index_map = marginals.get_index_map(parms, queue);
    const InverseOffsets inverse_offsets = marginals.get_inverse_offset_map(parms, queue);

    // The blocks of every event, end to end, are the whole array.
    std::size_t size = 0;
    for (const auto &event : parms.get_event_list()) {
        std::size_t block = event->size();
        for (const auto &entry : parents_by_stride(inverse_offsets, *event))
            block *= entry.first->size();
        size += block;
    }
    if (marginals.get_length() != size)
        throw std::logic_error("model_marginals_from_json: the marginals are not sized for this model");

    std::set<std::string> nicknames;
    for (const auto &event : parms.get_event_list()) {
        const std::string &nickname = event->get_nickname();
        nicknames.insert(nickname);
        const auto node = nodes.find(nickname);
        if (node == nodes.end())
            throw std::runtime_error("model json: the parameters have no node for event \"" + nickname + "\"");
        try {
            const auto axes = read_axes(*node, *event, parents_by_stride(inverse_offsets, *event));
            place_values(node->at("values"), marginals.marginal_array_smart_p.get() + index_map.at(event->get_name()),
                         axes, 0, 0);
        } catch (const std::exception &e) {
            throw std::runtime_error("model json: the parameters of \"" + nickname + "\": " + e.what());
        }
    }
    for (const auto &entry : nodes.items()) {
        if (nicknames.count(entry.key()) == 0)
            throw std::runtime_error("model json: the parameters have a node for \"" + entry.key()
                                     + "\", which is no event");
    }

    // As txt2marginals() does: values written to a finite precision need not sum to one exactly.
    // Values that already do, to a double's precision, as IGoR writes them, are kept as written:
    // normalizing would move each by a few units in its last place wherever long double is double
    // (arm64, MSVC), and the document written again would no longer be the same text. The bound is
    // far above that rounding (about 1e-16 per term summed), far below a text file's 6 digits.
    constexpr long double kAlreadyNormalized = 1e-12L;
    long double *array = marginals.marginal_array_smart_p.get();
    const std::vector<long double> as_written(array, array + size);
    marginals.normalize(inverse_offsets, index_map, queue);
    if (std::equal(as_written.begin(), as_written.end(), array, [](long double written, long double normalized) {
            return std::fabs(normalized - written) <= kAlreadyNormalized * std::fabs(written);
        }))
        std::copy(as_written.begin(), as_written.end(), array);
}

} // namespace igor::model::legacy
