#include <igor/Model/Legacy/gene_to_seqtype_migr.h>

using namespace std;

namespace igor {
namespace migration {

bool try_gene_class_to_gene_seq_type(Gene_class_legacy gene, Seq_type &seq_type)
{
  switch (gene) {
  case V_gene_legacy:
    seq_type = V_gene_seq;
    return true;
  case D_gene_legacy:
    seq_type = D_gene_seq;
    return true;
  case J_gene_legacy:
    seq_type = J_gene_seq;
    return true;
  default:
    return false;
  }
}

bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type)
{
  switch (gene) {
  case V_gene:
    seq_type = V_gene_seq;
    return true;
  case D_gene:
    seq_type = D_gene_seq;
    return true;
  case J_gene:
    seq_type = J_gene_seq;
    return true;
  default:
    return false;
  }
}

bool try_insertion_gene_class_to_seq_type(Gene_class_legacy gene_pair, Seq_type &seq_type)
{
  switch (gene_pair) {
  case VD_genes:
    seq_type = VD_ins_seq;
    return true;
  case DJ_genes:
    seq_type = DJ_ins_seq;
    return true;
  case VJ_genes:
    seq_type = VJ_ins_seq;
    return true;
  default:
    return false;
  }
}

bool try_insertion_gene_class_to_seq_type(Gene_class /*gene_pair*/, Seq_type & /*seq_type*/)
{
  // Gene_class (new slim enum) has no junction/insertion types; always fails.
  return false;
}

} // namespace migration
} // namespace igor
