#pragma once

#include <forward_list>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/SeqTypeRegistry.h>
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

CORE_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type);

CORE_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type);

CORE_EXPORT bool try_insertion_seq_type_to_gene_class(Seq_type seq_type, Gene_class &gene_pair);

CORE_EXPORT bool try_event_key_to_seq_key(
    Event_type event_type,
    Gene_class gene_class,
    Seq_side seq_side,
    std::tuple<Event_type, Seq_type, Seq_side> &seq_key);

CORE_EXPORT bool has_insertion_seq_type(
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Seq_type seq_type);

CORE_EXPORT bool has_insertion_seq_type(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Seq_type seq_type);

CORE_EXPORT bool try_get_event(
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Event_type event_type,
    Gene_class gene_class,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr);

CORE_EXPORT bool try_get_event(
    const std::unordered_map<std::tuple<Event_type, Seq_type, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    Event_type event_type,
    Seq_type seq_type,
    Seq_side seq_side,
    std::shared_ptr<Rec_Event> &event_ptr);

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

/// Legacy (boolean-flag) overload — assembles in hardcoded VDJ or VJ order.
CORE_EXPORT Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

/// Registry-based overload — assembles in the order defined by the SeqTypeRegistry.
/// @param registry   Ordered seq_type names (e.g. from Model_Parms::get_seq_type_registry()).
/// @param constructed_sequences  Map from seq_type string to the corresponding
///                               sequence fragment (nullptr entries are skipped).
CORE_EXPORT Int_Str build_scenario_sequence(
    const SeqTypeRegistry &registry,
    const std::unordered_map<std::string, const Int_Str *> &constructed_sequences);

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
