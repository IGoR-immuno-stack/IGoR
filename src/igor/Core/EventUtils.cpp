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

GeneChoiceStatus check_gene_choice(
    Gene_class gene,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events) {

  GeneChoiceStatus status;
  status.exists = false;
  status.chosen = false;
  status.event_ptr = nullptr;

  auto key = std::make_tuple(GeneChoice_t, gene, Undefined_side);
  if (events_map.count(key) != 0) {
    status.exists = true;
    status.event_ptr = events_map.at(key);
    if (processed_events.count(status.event_ptr->get_name()) != 0) {
      status.chosen = true;
    }
  }
  return status;
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
  if (events_map.count(key) != 0) {
    return events_map.at(key)->get_len_max();
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
