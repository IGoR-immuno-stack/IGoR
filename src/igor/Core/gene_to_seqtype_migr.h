#pragma once

#include <igor/Core/Utils.h>

#include <igorCoreExport.h>

namespace igor {
namespace migration {

CORE_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class_legacy gene, Seq_type &seq_type);
CORE_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type);

CORE_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class_legacy gene_pair, Seq_type &seq_type);
CORE_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type);

} // namespace migration
} // namespace igor
