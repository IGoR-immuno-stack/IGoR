#pragma once

/**
 * \file LegacyFormat.h
 * \brief The legacy (v1) model file's gene classes, translated to and from seq types.
 *
 * A v1 file names a segment by a Gene_class_legacy: its event lines say `#Deletion;V_gene;...`
 * or `#Insertion;VD_genes;...`, and its `@Edges` section lists generated event names that carry
 * the same token. A model in memory names a segment by its seq type. These two functions are
 * that translation, in both directions.
 *
 * Only the six legacy seq types have a gene class. Any other, a tandem D's included, gets
 * Undefined_gene; Model_Parms::write_model_parms() writes such a model as JSON.
 */

#include <igor/Core/Legacy/Utils.h>

#include <igor/Model/Export.h>

namespace igor::model::legacy {
using namespace igor::core::legacy;

/// The seq type a v1 gene class names: V_gene to V_gene_seq, VD_genes to VD_ins_seq, and so
/// on. "Undefined_seq" for VDJ_genes and Undefined_gene, which name no single segment.
MODEL_EXPORT Seq_type_String legacy_gene_class_to_seq_type(Gene_class_legacy gene_class);

/**
 * \brief The v1 gene class an event of this type carries for this seq type.
 *
 * The inverse of legacy_gene_class_to_seq_type() on the six legacy seq types. The event type
 * picks the half: a gene choice or a deletion acts on a gene segment and gets V, D or J; an
 * insertion or a Markov chain fills a junction and gets VD, DJ or VJ. Anything else, an
 * insertion on a gene segment included, gets Undefined_gene.
 */
MODEL_EXPORT Gene_class_legacy legacy_gene_class_of(Event_type event_type, const Seq_type_String &seq_type);

} // namespace igor::model::legacy
