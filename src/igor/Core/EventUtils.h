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
#include <igorCoreExport.h>

namespace EventUtils {

struct GeneChoiceStatus {
  bool exists;
  bool chosen;
  std::shared_ptr<const Rec_Event> event_ptr;
};

// gene_seq_type: seq_type of the GeneChoice event to look up (e.g. "V_gene_seq")
CORE_EXPORT GeneChoiceStatus check_gene_choice(
    const std::string &gene_seq_type,
    const Events_map &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events);

CORE_EXPORT Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

CORE_EXPORT void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets);

// ins_seq_type: seq_type of the Insertion event to look up (e.g. "VD_ins_seq")
CORE_EXPORT int get_insertion_len_max(
    const std::string &ins_seq_type,
    const Events_map &events_map);
} // namespace EventUtils
