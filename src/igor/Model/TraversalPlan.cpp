/*
 * TraversalPlan.cpp
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <igor/Model/TraversalPlan.h>

#include <igor/Model/Legacy/Rec_Event.h>

#include <limits>
#include <stdexcept>

namespace igor::model {

using core::kNoSeqType;
using core::SegmentSpan;
using core::SeqSide;
using core::SeqTypeId;

namespace {

constexpr std::size_t kNotInOrdering = std::numeric_limits<std::size_t>::max();

std::size_t end_index(SeqTypeId segment, SeqSide side)
{
    return 2 * static_cast<std::size_t>(segment) + (side == core::Three_prime ? 1 : 0);
}

std::string event_label(const legacy::Rec_Event &event)
{
    return event.get_nickname().empty() ? event.get_name() : event.get_nickname();
}

std::string side_label(SeqSide side) { return side == core::Five_prime ? "5'" : "3'"; }

/// Ordering positions by id, kNotInOrdering for an id outside the ordering.
std::vector<std::size_t> positions(const core::SeqTypeRegistry &registry)
{
    std::vector<std::size_t> position(registry.total_count(), kNotInOrdering);
    const auto &ordering = registry.ordering();
    for (std::size_t p = 0; p != ordering.size(); ++p) {
        position[ordering[p]] = p;
    }
    return position;
}

bool in_ordering(const std::vector<std::size_t> &position, SeqTypeId segment)
{
    return segment != kNoSeqType and static_cast<std::size_t>(segment) < position.size()
           and position[segment] != kNotInOrdering;
}

/// Mark what exists before the first event runs.
void apply(const InitialState &initial, std::vector<bool> &end_created, std::vector<bool> &built)
{
    for (const SeqTypeId segment : initial.placed) {
        end_created.at(end_index(segment, core::Five_prime)) = true;
        end_created.at(end_index(segment, core::Three_prime)) = true;
    }
    for (const SeqTypeId segment : initial.built) {
        end_created.at(end_index(segment, core::Five_prime)) = true;
        end_created.at(end_index(segment, core::Three_prime)) = true;
        built.at(segment) = true;
    }
}

} // namespace

std::optional<SegmentSpan> EventPlacement::enclosing_gap() const
{
    if (built_left == kNoSeqType or built_right == kNoSeqType) {
        return std::nullopt;
    }
    return SegmentSpan::gap(built_left, built_right);
}

std::optional<SegmentSpan> EventPlacement::gap_toward(SeqSide side) const
{
    if (not target_built) {
        return enclosing_gap();
    }
    if (side == core::Five_prime) {
        return built_left == kNoSeqType ? std::nullopt : std::optional(SegmentSpan::gap(built_left, target));
    }
    return built_right == kNoSeqType ? std::nullopt : std::optional(SegmentSpan::gap(target, built_right));
}

TraversalPlan::TraversalPlan(const core::SeqTypeRegistry &registry,
                             std::span<const legacy::Rec_Event *const> order, InitialState initial)
    : registry_(registry), order_(order.begin(), order.end()), initial_(std::move(initial)), index_(registry, order)
{
    const auto &ordering = registry.ordering();
    const std::vector<std::size_t> position = positions(registry);
    std::vector<bool> end_created(2 * registry.total_count(), false);
    std::vector<bool> built(registry.total_count(), false);
    apply(initial_, end_created, built);

    //The stretches between consecutive built segments, 5' to 3'.
    const auto open_gaps = [&]() {
        std::vector<SegmentSpan> gaps;
        SeqTypeId previous = kNoSeqType;
        for (const SeqTypeId segment : ordering) {
            if (not built[segment]) {
                continue;
            }
            if (previous != kNoSeqType) {
                gaps.push_back(SegmentSpan::gap(previous, segment));
            }
            previous = segment;
        }
        return gaps;
    };

    //A gap needs a bound while something in it is still to be decided: a segment inside it not
    //built yet, or an event from `from` on that takes part in it (R12's rule). Two built
    //neighbours with nothing between and nothing pending have nothing left to bound.
    const auto needs_bound = [&](const SegmentSpan &gap, std::size_t from) {
        for (std::size_t p = position[gap.left.id] + 1; p < position[gap.right.id]; ++p) {
            if (not built[ordering[p]]) {
                return true;
            }
        }
        for (std::size_t j = from; j < order_.size(); ++j) {
            if (order_[j]->participates_in_span(gap, registry)) {
                return true;
            }
        }
        return false;
    };

    const auto record = [&](std::size_t from) {
        for (const SegmentSpan &gap : open_gaps()) {
            const auto key = std::make_pair(gap.left.id, gap.right.id);
            if (gap_keys_.count(key) == 0 and needs_bound(gap, from)) {
                gap_keys_.emplace(key, kNoSeqType);
                gaps_.push_back(gap);
            }
        }
    };

    for (std::size_t i = 0; i != order_.size(); ++i) {
        const legacy::Rec_Event &event = *order_[i];
        EventPlacement placed;
        placed.target = event.get_seq_type_id();
        if (in_ordering(position, placed.target)) {
            const std::size_t here = position[placed.target];
            for (std::size_t p = here; p-- > 0;) {
                const SeqTypeId segment = ordering[p];
                if (end_created[end_index(segment, core::Three_prime)]) {
                    placed.placed_left.push_back(segment);
                }
                if (built[segment] and placed.built_left == kNoSeqType) {
                    placed.built_left = segment;
                }
            }
            for (std::size_t p = here + 1; p < ordering.size(); ++p) {
                const SeqTypeId segment = ordering[p];
                if (end_created[end_index(segment, core::Five_prime)]) {
                    placed.placed_right.push_back(segment);
                }
                if (built[segment] and placed.built_right == kNoSeqType) {
                    placed.built_right = segment;
                }
            }
        }

        record(i); //the gaps open while it runs, itself among what is still to run

        for (const SeqTypeId segment : ordering) {
            for (const SeqSide side : {core::Five_prime, core::Three_prime}) {
                if (event.get_offset_role(segment, side) == legacy::OffsetRole::Creates) {
                    end_created[end_index(segment, side)] = true;
                }
            }
            const legacy::SeqConstructionRole content = event.get_seq_construction_role(segment);
            if (content == legacy::SeqConstructionRole::Creates or content == legacy::SeqConstructionRole::Fills) {
                built[segment] = true;
            }
        }
        placed.target_built = in_ordering(position, placed.target) and built[placed.target];

        record(i + 1); //the gaps it opens: a placed gene splits the one it sat in into two

        //By identifier, not by address: a copy of the model has other events with the same ones.
        if (not placements_.emplace(event.get_event_identifier(), std::move(placed)).second) {
            throw std::invalid_argument("TraversalPlan: event " + event_label(event) + " shares its identifier "
                                        + std::to_string(event.get_event_identifier())
                                        + " with an earlier event of the order");
        }
    }

    //Keys, now that every gap is known. See the class comment for why a gap may take the id of
    //the one segment inside it.
    std::size_t next_key = registry.total_count();
    for (const SegmentSpan &gap : gaps_) {
        const std::size_t first_inside = position[gap.left.id] + 1;
        const std::size_t past_inside = position[gap.right.id];
        SeqTypeId key = kNoSeqType;
        if (past_inside == first_inside + 1 and not index_.has_own_factor(ordering[first_inside])) {
            key = ordering[first_inside];
        } else {
            if (next_key >= static_cast<std::size_t>(kNoSeqType)) {
                throw std::length_error("TraversalPlan: too many gap keys for the SeqTypeId range");
            }
            key = static_cast<SeqTypeId>(next_key++);
        }
        gap_keys_.at({gap.left.id, gap.right.id}) = key;
    }
    bound_key_count_ = next_key;
}

const EventPlacement &TraversalPlan::placement(const legacy::Rec_Event &event) const
{
    const auto it = placements_.find(event.get_event_identifier());
    if (it == placements_.end()) {
        throw std::out_of_range("TraversalPlan::placement(): event " + event_label(event) + " (identifier "
                                + std::to_string(event.get_event_identifier())
                                + ") is not in the order this plan was built for");
    }
    return it->second;
}

SeqTypeId TraversalPlan::gap_key(SegmentSpan gap) const
{
    const auto it = gap_keys_.find({gap.left.id, gap.right.id});
    if (it == gap_keys_.end()) {
        const auto name = [this](SeqTypeId id) {
            return id < registry_.total_count() ? registry_.name(id) : std::to_string(id);
        };
        throw std::out_of_range("TraversalPlan::gap_key(): " + name(gap.left.id) + " -> " + name(gap.right.id)
                                + " is never an open gap with something left to bound in this order");
    }
    return it->second;
}

void TraversalPlan::check_admissible() const
{
    const auto &ordering = registry_.ordering();
    const std::vector<std::size_t> position = positions(registry_);
    std::vector<bool> end_created(2 * registry_.total_count(), false);
    std::vector<bool> built(registry_.total_count(), false);
    apply(initial_, end_created, built);

    const auto refuse = [](const legacy::Rec_Event &event, const std::string &why) {
        throw std::invalid_argument("Event " + event_label(event) + " cannot run in this order: " + why);
    };

    for (const legacy::Rec_Event *event_p : order_) {
        const legacy::Rec_Event &event = *event_p;
        const SeqTypeId target = event.get_seq_type_id();
        if (not in_ordering(position, target)) {
            refuse(event, "its segment " + (target == kNoSeqType ? std::string("is unresolved") : registry_.name(target))
                                  + " is not in the model's 5'->3' ordering");
        }

        for (const SeqTypeId segment : ordering) {
            const std::string &name = registry_.name(segment);
            bool creates_an_end = false;
            for (const SeqSide side : {core::Five_prime, core::Three_prime}) {
                const legacy::OffsetRole role = event.get_offset_role(segment, side);
                if (role == legacy::OffsetRole::Creates) {
                    creates_an_end = true;
                    if (end_created[end_index(segment, side)]) {
                        refuse(event, "it creates the " + side_label(side) + " end of " + name
                                              + ", which an earlier event created");
                    }
                } else if (role == legacy::OffsetRole::Modifies and not end_created[end_index(segment, side)]) {
                    refuse(event, "it modifies the " + side_label(side) + " end of " + name
                                          + " before any event creates it");
                }
            }

            const legacy::SeqConstructionRole content = event.get_seq_construction_role(segment);
            const bool creates_content =
                    content == legacy::SeqConstructionRole::Creates or content == legacy::SeqConstructionRole::Fills;
            if (creates_content and built[segment]) {
                refuse(event, "it creates the content of " + name + ", which an earlier event created");
            }
            if (content == legacy::SeqConstructionRole::Modifies and not built[segment]) {
                refuse(event, "it modifies the content of " + name + " before any event builds it");
            }
            if (creates_content and not creates_an_end
                and not(end_created[end_index(segment, core::Five_prime)]
                        and end_created[end_index(segment, core::Three_prime)])) {
                refuse(event, "it fills " + name + " before both its ends are placed");
            }
            if (creates_an_end and not creates_content) {
                //A slot: its ends are derived from its neighbours' facing ends, so those come first
                //(B10 decides what "neighbour" means once a segment can be absent)...
                const SeqTypeId left = registry_.left_neighbor(segment);
                const SeqTypeId right = registry_.right_neighbor(segment);
                if (left == kNoSeqType or not end_created[end_index(left, core::Three_prime)]) {
                    refuse(event, "it places " + name + " before the 3' end of "
                                          + (left == kNoSeqType ? std::string("a 5' neighbour") : registry_.name(left))
                                          + " is placed");
                }
                if (right == kNoSeqType or not end_created[end_index(right, core::Five_prime)]) {
                    refuse(event, "it places " + name + " before the 5' end of "
                                          + (right == kNoSeqType ? std::string("a 3' neighbour") : registry_.name(right))
                                          + " is placed");
                }
                //...and it has to sit in a gap, between built content on both sides.
                bool built_before = false;
                for (std::size_t p = 0; p < position[segment]; ++p) {
                    built_before = built_before or built[ordering[p]];
                }
                bool built_after = false;
                for (std::size_t p = position[segment] + 1; p < ordering.size(); ++p) {
                    built_after = built_after or built[ordering[p]];
                }
                if (not built_before or not built_after) {
                    refuse(event, "it places " + name + " outside a gap: nothing is built on its "
                                          + (built_before ? "3'" : "5'") + " side yet");
                }
            }
        }

        for (const SeqTypeId segment : ordering) {
            for (const SeqSide side : {core::Five_prime, core::Three_prime}) {
                if (event.get_offset_role(segment, side) == legacy::OffsetRole::Creates) {
                    end_created[end_index(segment, side)] = true;
                }
            }
            const legacy::SeqConstructionRole content = event.get_seq_construction_role(segment);
            if (content == legacy::SeqConstructionRole::Creates or content == legacy::SeqConstructionRole::Fills) {
                built[segment] = true;
            }
        }
    }

    for (const SeqTypeId segment : ordering) {
        const std::string &name = registry_.name(segment);
        for (const SeqSide side : {core::Five_prime, core::Three_prime}) {
            if (not end_created[end_index(segment, side)]) {
                throw std::invalid_argument("No event creates the " + side_label(side) + " end of " + name);
            }
        }
        if (not built[segment]) {
            throw std::invalid_argument("No event builds the content of " + name);
        }
    }
}

} // namespace igor::model
