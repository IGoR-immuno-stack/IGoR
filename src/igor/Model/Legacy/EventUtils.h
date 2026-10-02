#pragma once

#include <forward_list>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Core/Legacy/SeqTypeRegistry.h>
#include <igor/Core/Legacy/Utils.h>
#include <igor/Model/Legacy/EventTypedefs.h>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <igor/Core/Legacy/IntStr.h>
#include <igor/Model/Export.h>

namespace EventUtils {

MODEL_EXPORT bool has_insertion_seq_type(
    const Events_map &events_map,
    Seq_type seq_type);

MODEL_EXPORT bool try_get_event(
    const Events_map &events_map,
    Event_type event_type,
    Seq_type seq_type,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr);

/// Convert a Seq_type enum value to its canonical string name used in Events_map keys.
MODEL_EXPORT Seq_type_String seq_type_to_string(Seq_type seq_type);

/// Whether a model has a GeneChoice on a segment, and whether it has been realized yet.
struct GeneChoiceStatus {
  bool exists;
  bool chosen;
  std::shared_ptr<const Rec_Event> event_ptr;
};

// gene_seq_type: seq_type of the GeneChoice event to look up (e.g. "V_gene_seq")
MODEL_EXPORT GeneChoiceStatus check_gene_choice(
    const Seq_type_String &gene_seq_type,
    const Events_map &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events);

/// Legacy (boolean-flag) overload — assembles in hardcoded VDJ or VJ order.
MODEL_EXPORT Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

/// Registry-based overload — assembles in the order defined by the SeqTypeRegistry.
/// @param registry   Ordered seq_type names (e.g. from Model_Parms::get_seq_type_registry()).
/// @param constructed_sequences  Map from seq_type string to the corresponding
///                               sequence fragment (nullptr entries are skipped).
MODEL_EXPORT Int_Str build_scenario_sequence(
    const SeqTypeRegistry &registry,
    const std::unordered_map<Seq_type_String, const Int_Str *> &constructed_sequences);

MODEL_EXPORT void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets);

// ins_seq_type: seq_type of the Insertion event to look up (e.g. "VD_ins_seq")
MODEL_EXPORT int get_insertion_len_max(
    const Seq_type_String &ins_seq_type,
    const Events_map &events_map);
} // namespace EventUtils
