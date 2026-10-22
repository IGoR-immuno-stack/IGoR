#pragma once

/**
 * \file ModelJson.h
 * \brief The JSON model document: written from a Model_Parms, read into one or into a Topology.
 *
 * The layout is the one SEGMENT_DECOMPOSITION_REVIEW.md §14 and ARCHITECTURE_SYNTHESIS.md
 * (principle 6) plan, restricted to what the legacy model holds:
 *
 *     {
 *       "$schema": "https://raw.githubusercontent.com/.../schemas/0.3/model.schema.json",
 *       "schema_version": "0.3",
 *       "structure": {
 *         "layout":       { "type": "Linear", "segments": ["V_gene_seq", "VD_ins_seq", ...] },
 *         "events":       { "v_choice": { "type": "GeneChoice", "library": "V_gene", "priority": 7,
 *                                         "targets": [ { "segment": "V_gene_seq" } ],
 *                                         "realizations": [ { "index": 0, "label": "TRBV1*01" }, ... ] },
 *                           "v_3_del":  { "type": "Deletion", "priority": 5,
 *                                         "targets": [ { "segment": "V_gene_seq", "side": "3'" } ],
 *                                         "realizations": [ { "index": 0, "label": -4 }, ... ] }, ... },
 *         "conditioning": { "v_3_del": ["v_choice"], ... }
 *       },
 *       "library":    { "V_gene": { "TRBV1*01": { "sequence": "GATACT..." }, ... }, ... },
 *       "parameters": { "nodes": { "v_3_del": { "axes": [ { "event": "v_choice", "labels": ["TRBV1*01", ...] },
 *                                                          { "event": "v_3_del", "labels": [-4, ...] } ],
 *                                               "values": [ [0.01, ...], ... ] }, ... } },
 *       "error_rate": { "type": "SingleErrorRate", "rate": 0.000846 },
 *       "engine_hints": { "legacy_dfs": { "event_order": ["v_choice", "d_gene", ...] } }
 *     }
 *
 * - Every reader checks the document against its JSON Schema first (ModelSchema.h): keys, types,
 *   enumerations, ranges. What follows is checked in code, as no schema can relate one part of a
 *   document to another. "$schema" names the schema, so that an editor validates and completes.
 * - The events are an object keyed by nickname, the identity. The order of its members is never
 *   read (§14.2), and a key repeated in one object is refused at parse (parse_model_document()).
 *   The readers build the events in byte-wise nickname order, unless the document carries the
 *   legacy engine's event order: the writer records it for a model whose event list is in
 *   another order, a v1 model's for one, so that its queue and its marginals file still fit
 *   (ENGINE_MIGRATION_PLAN Q4).
 * - The conditioning graph is its own section, child to parents, apart from the events.
 * - A realization is its label, its identity, and its index, its coordinate in the marginal
 *   array. Indices run from 0 without gaps: the legacy marginals file is read by position.
 * - Templates sit in the library, keyed by label, so that D1 and D2 share one. A gene choice's
 *   node names its library, and its constructor takes the templates from there.
 * - The parameters, the marginals of the legacy model, are optional. A document may hold them
 *   with the rest of the model, or alone with its version: a fragment, which a marginals file is.
 *   Each node lists its axes, its event's parents then its own, each axis with the labels of its
 *   coordinates, so values are matched to realizations by label (REVIEW §15.4), whatever order
 *   the model reading them gives its events. Values are doubles. An axis may leave its labels out,
 *   to prototype: it then follows the index order of the model reading it, and only its size is
 *   checked. IGoR always writes them.
 *
 * Comments are allowed, and dropped when IGoR writes a model.
 *
 * Not yet, and refused when present: the `provenance` and `accumulators` sections, the
 * parameters' structure fingerprint and compact axes (synthesis step 4b), several targets per
 * event, a layout other than Linear (step 5), a parameterization other than the event type's, and
 * the authoring form's shorthands and `description` fields (§14.4-14.5). `error_rate` is the one
 * key outside the plan: it stands for the error model until that is a node of the structure with
 * its value in `parameters` (ERROR_MODEL_REWRITE_PLAN.md, E4).
 */

#include <igor/Core/Types.h>
#include <igor/Model/Export.h>

#include <nlohmann/json_fwd.hpp>

#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace igor::model::legacy {

class Model_marginals;
class Model_Parms;
class Rec_Event;

/// Schema version written into every document, and the only version the reader accepts.
///
/// A 0.x version is a pre-release one. Each replaces the one before it, and no reader keeps
/// support for an older 0.x: none leaves the repository, so none has a user to keep, as the v2
/// text format had none. A released version, from 1.0, stays readable by every later reader.
inline constexpr const char *kModelJsonSchemaVersion = "0.3";

/**
 * \brief Parse a model document.
 *
 * Comments, // and C-style, are allowed anywhere whitespace is, so that a model or a test
 * fixture can say what it is (REVIEW §14.5). IGoR does not keep them: a model it writes has none.
 * \throws std::runtime_error naming the key and its object if a key appears twice in one
 *         object: nlohmann::json would keep the last one without a word.
 * \throws nlohmann::json::parse_error if the text is not JSON.
 */
MODEL_EXPORT nlohmann::json parse_model_document(std::istream &in);

/// Whether a model or marginals file is a JSON document: one opens with its object or with a
/// comment, a v1 file with an @ section. Skips the leading whitespace.
MODEL_EXPORT bool starts_like_json_document(std::istream &in);

/// Whether \a filename is a JSON document that carries a "parameters" section.
MODEL_EXPORT bool model_file_has_parameters(const std::string &filename);

/// The structure of a document, checked: what model_parms_from_json() and topology_from_json()
/// both build from.
struct DocumentStructure {
    std::vector<core::SeqTypeString> segments;           ///< the Linear layout, 5' to 3'
    std::vector<std::shared_ptr<Rec_Event>> events;      ///< built by the event factory, in build order
    std::vector<std::pair<std::string, std::string>> edges;  ///< (parent, child) nicknames
};

/**
 * \brief Check a document and build its events.
 *
 * Each event goes through the event factory, a gene choice with its templates taken from the
 * library, in nickname order or the legacy engine's. The conditioning is checked to name events
 * only, each parent once, without a cycle.
 * \throws std::runtime_error on another schema version, an unknown or missing key, or any of
 *         the above; an event's error names its nickname.
 */
MODEL_EXPORT DocumentStructure read_document_structure(const nlohmann::json &doc);

/**
 * \brief Serialize a model into a document.
 * \param parms a finalized model, so its seq_type registry holds its segment order.
 * \throws std::runtime_error if the model carries an error rate this function cannot serialize
 *         yet, or two templates under one label in one library. Failing loudly is deliberate: a
 *         document silently missing either would look complete and rebuild a different model.
 */
MODEL_EXPORT nlohmann::json model_parms_to_json(const Model_Parms &parms);

/**
 * \brief Fill an empty model from a document.
 *
 * read_document_structure(), then the segment order, the events, the edges and the error rate
 * into \a parms, and Model_Parms::finalize().
 * \throws std::runtime_error as read_document_structure() does, or on an error rate other than
 *         SingleErrorRate: the hypermutation models have no JSON form yet (step HJ of
 *         docs/ERROR_MODEL_REWRITE_PLAN.md).
 * \throws std::logic_error if \a parms already has events.
 */
MODEL_EXPORT void model_parms_from_json(const nlohmann::json &doc, Model_Parms &parms);

/**
 * \brief The "parameters" section of \a marginals, laid out by \a parms.
 *
 * One node per event, keyed by nickname. Its axes come in the legacy layout's order, the parents
 * in the order of the marginals file's $Dim line, then its own axis, or a Markov chain's two (from,
 * then to). Each lists its event's labels in index order, so the values, nested one level per axis,
 * read straight through on the model that wrote them. A value is the long double rounded to a
 * double: 17 significant digits, where the text file has 6.
 */
MODEL_EXPORT nlohmann::json model_marginals_to_json(const Model_Parms &parms, const Model_marginals &marginals);

/**
 * \brief Fill \a marginals from the "parameters" section of a document, then normalize them as
 * Model_marginals::txt2marginals() does.
 *
 * \a doc may be a complete model document or a fragment, its version and "parameters" only.
 * Every event of \a parms needs its node, and no node may name another event. A node's axes are
 * its event's parents, in any order, then its own axis or axes. An axis's labels are its event's,
 * each once, in any order: a value goes to the realization its label names. Nothing therefore
 * depends on the order \a parms gives its events. An axis without labels follows the index order
 * of \a parms, by position.
 * \throws std::runtime_error naming the event on any mismatch, or on a value that is negative or
 *         not a finite number.
 * \throws std::logic_error if \a marginals is not sized for \a parms.
 */
MODEL_EXPORT void model_marginals_from_json(const nlohmann::json &doc, const Model_Parms &parms,
                                            Model_marginals &marginals);

} // namespace igor::model::legacy
