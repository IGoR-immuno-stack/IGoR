#include <igor/Core/EventUtils.h>
#include <igor/Core/Rec_Event.h>

namespace EventUtils {

Seq_type_String seq_type_to_string(Seq_type seq_type)
{
    switch (seq_type) {
    case V_gene_seq:  return "V_gene_seq";
    case VD_ins_seq:  return "VD_ins_seq";
    case D_gene_seq:  return "D_gene_seq";
    case DJ_ins_seq:  return "DJ_ins_seq";
    case J_gene_seq:  return "J_gene_seq";
    case VJ_ins_seq:  return "VJ_ins_seq";
    default:          return "Undefined_seq";
    }
}

bool has_insertion_seq_type(
    const Events_map &events_map,
    Seq_type seq_type)
{
    auto key = std::make_tuple(Insertion_t, seq_type_to_string(seq_type), Undefined_side);
    return events_map.count(key) != 0;
}

bool try_get_event(
    const Events_map &events_map,
    Event_type event_type,
    Seq_type seq_type,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr)
{
    auto key = std::make_tuple(event_type, seq_type_to_string(seq_type), seq_side);
    auto it = events_map.find(key);
    if (it != events_map.end()) {
        event_ptr = it->second;
        return true;
    }
    return false;
}

igor::migration::GeneChoiceStatus check_gene_choice(
  const Seq_type_String &gene_seq_type,
    const Events_map &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events) {

  igor::migration::GeneChoiceStatus status;
  status.exists = false;
  status.chosen = false;
  status.event_ptr = nullptr;

  auto key = std::make_tuple(GeneChoice_t, gene_seq_type, Undefined_side);
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

Int_Str build_scenario_sequence(
    const SeqTypeRegistry &registry,
  const std::unordered_map<Seq_type_String, const Int_Str *> &constructed_sequences) {

  Int_Str scenario_resulting_sequence;
  for (const auto &seq_type : registry.get_ordered_types()) {
    auto it = constructed_sequences.find(seq_type);
    if (it != constructed_sequences.end() && it->second != nullptr) {
      scenario_resulting_sequence += *(it->second);
    }
  }
  return scenario_resulting_sequence;
}

int get_insertion_len_max(
  const Seq_type_String &ins_seq_type,
    const Events_map &events_map) {
  auto key = std::make_tuple(Insertion_t, ins_seq_type, Undefined_side);
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
