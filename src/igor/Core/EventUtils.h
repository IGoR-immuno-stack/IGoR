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

CORE_EXPORT bool has_insertion_seq_type(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Seq_type seq_type);

CORE_EXPORT bool try_get_event(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Event_type event_type,
    Seq_type seq_type,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr);



CORE_EXPORT Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

CORE_EXPORT void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets);

CORE_EXPORT int get_insertion_len_max(
    Seq_type seq_type,
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map);
} // namespace EventUtils
