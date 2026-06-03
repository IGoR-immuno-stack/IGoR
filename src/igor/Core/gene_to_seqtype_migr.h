#pragma once

#include <forward_list>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>

#include <igorCoreExport.h>

namespace igor {
namespace migration {

using LegacyEventsMap = std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>;

CORE_EXPORT bool try_gene_class_to_gene_seq_type(Gene_class gene, Seq_type &seq_type);

CORE_EXPORT bool try_insertion_gene_class_to_seq_type(Gene_class gene_pair, Seq_type &seq_type);

CORE_EXPORT bool try_insertion_seq_type_to_gene_class(Seq_type seq_type, Gene_class &gene_pair);

CORE_EXPORT bool try_event_key_to_seq_key(
    Event_type event_type,
    Gene_class gene_class,
    Seq_side seq_side,
    std::tuple<Event_type, Seq_type, Seq_side> &seq_key);

struct GeneChoiceStatus {
  bool exists;
  bool chosen;
  std::shared_ptr<const Rec_Event> event_ptr;
};

CORE_EXPORT const LegacyEventsMap &empty_legacy_events_map();

} // namespace migration
} // namespace igor
