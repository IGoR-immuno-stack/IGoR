#include <igor/Core/gene_to_seqtype_migr.h>
#include <igor/Core/EventUtils.h> // for try_get_event

using namespace std;

namespace igor {
namespace migration {

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

bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type)
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

bool try_insertion_seq_type_to_gene_class(Seq_type seq_type, Gene_class &gene_pair)
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
    Gene_class gene_class,
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

SeqEventsMap build_seq_events_map(const LegacyEventsMap &events_map)
{
    SeqEventsMap seq_events_map;
    seq_events_map.reserve(events_map.size());
    for (const auto &entry : events_map) {
        tuple<Event_type, Seq_type, Seq_side> seq_key;
        if (!try_event_key_to_seq_key(get<0>(entry.first), get<1>(entry.first), get<2>(entry.first),
                                                   seq_key)) {
            continue;
        }
        seq_events_map.emplace(seq_key, entry.second);
    }
    return seq_events_map;
}

GeneChoiceStatus check_gene_choice_seq_type(
        Gene_class gene, const SeqEventsMap &events_map, const unordered_set<Rec_Event_name> &processed_events)
{
    GeneChoiceStatus status{ false, false, nullptr };
    Seq_type gene_seq = V_gene_seq;
    if (!try_gene_class_to_gene_seq_type(gene, gene_seq)) {
        return status;
    }

    shared_ptr<Rec_Event> event_ptr;
    if (!EventUtils::try_get_event(events_map, GeneChoice_t, gene_seq, Undefined_side, event_ptr)) {
        return status;
    }

    status.exists = true;
    status.chosen = processed_events.count(event_ptr->get_name()) != 0;
    status.event_ptr = event_ptr;
    return status;
}

} // namespace migration
} // namespace igor
