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
 * The schema differs from the text format in three places, all deliberate:
 *   - `seq_type_order` is a real array instead of a bracket-delimited line;
 *   - edges are carried per event as a list of parent nicknames, so there is no separate
 *     `@Edges` section and no dependence on the generated event names;
 *   - realizations always come out sorted by index, so the output is reproducible.
 *
 * Reading this schema back, through a factory and per-class constructors, is step 2.
 */

#include <igor/Model/Export.h>

#include <nlohmann/json_fwd.hpp>

class Model_Parms;

namespace igor {

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

}  // namespace igor
