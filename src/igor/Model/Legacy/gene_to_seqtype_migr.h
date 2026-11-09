#pragma once

#include <igor/Core/Legacy/Utils.h>

#include <igor/Model/Export.h>

/**
 * Migration code (#59, CAP B1): the legacy `Seq_type` of a gene class's gene or N region segment.
 * Deleted, with its source, at CAP T6b, once nothing calls the two `try_*` families:
 * FastGenerator after CAP T4b, and the test helper make_gene_choice().
 */

namespace igor::model::legacy {
using namespace igor::core::legacy;
namespace migration {

MODEL_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class_legacy gene, Seq_type &seq_type);
MODEL_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type);

MODEL_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class_legacy gene_pair, Seq_type &seq_type);
MODEL_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type);

} // namespace migration
} // namespace igor::model::legacy