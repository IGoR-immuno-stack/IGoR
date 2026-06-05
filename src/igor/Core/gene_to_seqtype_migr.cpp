#include <igor/Core/gene_to_seqtype_migr.h>
#include <igor/Core/EventUtils.h> // for try_get_event

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

bool try_insertion_seq_type_to_gene_class(Seq_type seq_type, Gene_class_legacy &gene_pair)
{
  switch (seq_type) {
  case VD_ins_seq:
    gene_pair = VD_genes;
    return true;
  case DJ_ins_seq:
    gene_pair = DJ_genes;
    return true;
  case VJ_ins_seq:
    gene_pair = VJ_genes;
    return true;
  default:
    return false;
  }
}

bool try_event_key_to_seq_key(
    Event_type event_type,
    Gene_class_legacy gene_class,
    Seq_side seq_side,
    std::tuple<Event_type, Seq_type, Seq_side> &seq_key)
{
  Seq_type seq_type = V_gene_seq;
  bool mapped = false;

  switch (event_type) {
  case GeneChoice_t:
  case Deletion_t:
    mapped = try_gene_class_to_gene_seq_type(gene_class, seq_type);
    break;
  case Insertion_t:
  case Dinuclmarkov_t:
    mapped = try_insertion_gene_class_to_seq_type(gene_class, seq_type);
    break;
  default:
    mapped = false;
    break;
  }

  if (!mapped) {
    return false;
  }

  seq_key = std::make_tuple(event_type, seq_type, seq_side);
  return true;
}

const LegacyEventsMap &empty_legacy_events_map()
{
    static const LegacyEventsMap kEmptyLegacyEventsMap;
    return kEmptyLegacyEventsMap;
}

} // namespace migration
} // namespace igor
