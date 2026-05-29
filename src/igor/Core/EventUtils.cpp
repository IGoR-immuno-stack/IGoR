#include <igor/Core/EventUtils.h>
#include <igor/Core/Rec_Event.h>

namespace EventUtils {

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

bool has_insertion_seq_type(
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Seq_type seq_type)
{
  Gene_class gene_pair = Undefined_gene;
  if (!try_insertion_seq_type_to_gene_class(seq_type, gene_pair)) {
    return false;
  }

  auto key = std::make_tuple(Insertion_t, gene_pair, Undefined_side);
  return events_map.find(key) != events_map.end();
}

bool has_insertion_seq_type(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Seq_type seq_type)
{
  auto key = std::make_tuple(Insertion_t, seq_type, Undefined_side);
  return events_map.find(key) != events_map.end();
}

bool try_get_event(
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Event_type event_type,
    Gene_class gene_class,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr)
{
  auto key = std::make_tuple(event_type, gene_class, seq_side);
  auto event_it = events_map.find(key);
  if (event_it == events_map.end()) {
    event_ptr = nullptr;
    return false;
  }

  event_ptr = event_it->second;
  return true;
}

bool try_get_event(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Event_type event_type,
    Seq_type seq_type,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr)
{
  auto key = std::make_tuple(event_type, seq_type, seq_side);
  auto event_it = events_map.find(key);
  if (event_it == events_map.end()) {
    event_ptr = nullptr;
    return false;
  }

  event_ptr = event_it->second;
  return true;
}

Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins) {

  Int_Str scenario_resulting_sequence;

  if (has_v) {
    scenario_resulting_sequence += (*constructed_sequences[V_gene_seq]);
  }

  if (has_d) {
    if (has_vd_ins) {
      scenario_resulting_sequence += (*constructed_sequences[VD_ins_seq]);
    }
    scenario_resulting_sequence += (*constructed_sequences[D_gene_seq]);
    if (has_dj_ins) {
      scenario_resulting_sequence += (*constructed_sequences[DJ_ins_seq]);
    }
  } else {
    if (has_vj_ins) {
      scenario_resulting_sequence += (*constructed_sequences[VJ_ins_seq]);
    }
  }

  if (has_j) {
    scenario_resulting_sequence += (*constructed_sequences[J_gene_seq]);
  }

  return scenario_resulting_sequence;
}

int get_insertion_len_max(
    Gene_class gene_pair,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map) {
  auto key = std::make_tuple(Insertion_t, gene_pair, Undefined_side);
  auto event_it = events_map.find(key);
  if (event_it != events_map.end()) {
    return event_it->second->get_len_max();
  }
  return 0;
}

int get_insertion_len_max(
    Seq_type seq_type,
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map) {
  auto key = std::make_tuple(Insertion_t, seq_type, Undefined_side);
  auto event_it = events_map.find(key);
  if (event_it != events_map.end()) {
    return event_it->second->get_len_max();
  }
  return 0;
}

void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets) {

  for (auto iter = offset_vector.begin(); iter != offset_vector.end(); ++iter) {
    int event_identitfier = (*iter).first->get_event_identifier();
    index_map.request_memory_layer(event_identitfier);
    memory_and_offsets.emplace_front(
        event_identitfier,
        index_map.get_current_memory_layer(event_identitfier), (*iter).second);
  }
}
} // namespace EventUtils
