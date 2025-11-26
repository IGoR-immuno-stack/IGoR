#pragma once

#include <forward_list>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <igor/Core/IntStr.h>

namespace EventUtils {

struct GeneChoiceStatus {
  bool exists;
  bool chosen;
  std::shared_ptr<const Rec_Event> event_ptr;
};

GeneChoiceStatus check_gene_choice(
    Gene_class gene,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events);

Int_Str build_scenario_sequence(const Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets);

int get_insertion_len_max(
    Gene_class gene_pair,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map);
} // namespace EventUtils
