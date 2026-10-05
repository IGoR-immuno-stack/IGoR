#pragma once

#include <igor/Core/Legacy/Utils.h>

#include <igor/Model/Export.h>

namespace igor::model::legacy {
using namespace igor::core::legacy;
namespace migration {

MODEL_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class_legacy gene, Seq_type &seq_type);
MODEL_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type);

MODEL_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class_legacy gene_pair, Seq_type &seq_type);
MODEL_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type);

} // namespace migration
} // namespace igor::model::legacy