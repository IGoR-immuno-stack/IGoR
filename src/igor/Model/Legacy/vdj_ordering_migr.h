#pragma once

#include <igor/Core/Legacy/SeqTypeRegistry.h>

#include <igor/Model/Export.h>

#include <string>

/**
 * Migration code (T7a): deleted, with its source, once nothing reads V, D and J by the `Seq_type`
 * enum. Its callers go at T7b for the two hypermutation error models and at CT2a of
 * docs/COUNTER_REWRITE_PLAN.md for the four counters.
 */

namespace igor::model::legacy {
namespace migration {

/**
 * \brief Refuse a model whose ordering is neither V-D-J nor V-J, by the legacy segment names.
 *
 * For the readers that still find V, D and J at the ids the `Seq_type` enum pins and assemble a
 * scenario in that order: the two hypermutation error models and the four counters. On any other
 * ordering, a tandem D's V-D1-D2-J or a VDJ model whose segments are named otherwise, they would
 * answer without an error, and wrongly. Each calls this from its `check_admissible()`, which
 * GenModel::infer_model() runs before anything else.
 *
 * \param reader what refuses, named in the message
 * \throws std::invalid_argument naming the ordering and the first segment of it the reader cannot
 *         place, or the segment it lacks
 */
MODEL_EXPORT void require_vdj_or_vj_ordering(const igor::core::SeqTypeRegistry &registry, const std::string &reader);

} // namespace migration
} // namespace igor::model::legacy
