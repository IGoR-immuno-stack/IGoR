#pragma once

/**
 * \file ModelJson.h
 * \brief JSON serialization of a parsed model.
 *
 * Step 1 of the move to a JSON model format. Core keeps its text parser: this header does
 * not read anything, it serializes a Model_Parms that has already been read. That way there
 * is exactly one tokenizer for the text format, and the JSON schema is defined from the
 * objects rather than from a second parser that could drift from the first.
 *
 * The schema differs from the text format in four places, all deliberate:
 *   - `seq_type_order` is a real array instead of a bracket-delimited line;
 *   - edges are carried per event as a list of parent nicknames, so there is no separate
 *     `@Edges` section and no dependence on the generated event names;
 *   - realizations always come out sorted by index, so the output is reproducible;
 *   - only a gene choice carries `gene_class`, the gene library it draws from. On the other
 *     events the text format's class field is a v1 token their seq type already gives.
 *
 * model_parms_from_json() reads it back into a Model_Parms, through the event factory and the
 * per-class constructors. It is the model file format for any model the v1 format cannot hold,
 * a tandem D's for instance (plan T2), and replaced the v2 text format.
 */

#include <igor/Model/Export.h>

#include <nlohmann/json_fwd.hpp>

namespace igor::model::legacy {

class Model_Parms;

/// Schema version written into every document, and the only version the future reader
/// will accept without a migration.
constexpr int kModelJsonSchemaVersion = 1;

/**
 * \brief Serialize a model into the JSON schema.
 * \param parms a model that has already been read, so its seq_type registry is frozen.
 * \throws std::runtime_error if the model carries an error rate this function cannot
 *         serialize yet. Failing loudly is deliberate: a document silently missing its
 *         error rate would look complete and rebuild a different model.
 */
MODEL_EXPORT nlohmann::json model_parms_to_json(const Model_Parms &parms);

/**
 * \brief Fill an empty model from a document model_parms_to_json() wrote.
 *
 * The segment order, the events through the event factory, the edges from each event's
 * parents, the error rate; then Model_Parms::finalize(). A key the schema does not have is
 * refused, as it is in an event node.
 * \throws std::runtime_error on another schema version, an unknown key, a parent no event is
 *         named, or an error rate other than SingleErrorRate: the hypermutation models have no
 *         JSON form yet (step HJ of docs/ERROR_MODEL_REWRITE_PLAN.md).
 * \throws std::logic_error if \a parms already has events.
 */
MODEL_EXPORT void model_parms_from_json(const nlohmann::json &doc, Model_Parms &parms);

} // namespace igor::model::legacy