/*
 * SegmentIndex.h
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

#pragma once

#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Types.h>
#include <igor/Model/Legacy/Rec_Event.h>

#include <cstddef>
#include <span>
#include <vector>

namespace igor::model {

/**
 * \brief Which events create and which modify each segment of a model: its two ends, and its
 * content.
 *
 * Read from the events' A0 roles (`get_offset_role()`, `get_seq_construction_role()`), over the
 * segments of the registry's ordering. It does not depend on the order the events run in, so every
 * engine reads the same index; what does depend on it is TraversalPlan's. This is the synthesis'
 * "SegmentIndex (creators and modifiers per segment and end)" (ARCHITECTURE_SYNTHESIS §2).
 *
 * Nothing here knows a gene from a junction. A gene choice creates a segment's ends and its
 * content at once; an insertion creates the ends of a slot, and the Markov chain after it the
 * content; a deletion modifies an end and the content. Those are the facts the index holds.
 */
class SegmentIndex
{
public:
    using Events = std::vector<const legacy::Rec_Event *>;

    SegmentIndex(const core::SeqTypeRegistry &registry, std::span<const legacy::Rec_Event *const> events)
        : end_creators_(2 * registry.total_count()),
          end_modifiers_(2 * registry.total_count()),
          content_creators_(registry.total_count()),
          content_modifiers_(registry.total_count())
    {
        for (const legacy::Rec_Event *event : events) {
            for (const core::SeqTypeId segment : registry.ordering()) {
                for (const core::SeqSide side : {core::Five_prime, core::Three_prime}) {
                    switch (event->get_offset_role(segment, side)) {
                    case legacy::OffsetRole::Creates: end_creators_[end(segment, side)].push_back(event); break;
                    case legacy::OffsetRole::Modifies: end_modifiers_[end(segment, side)].push_back(event); break;
                    case legacy::OffsetRole::None: break;
                    }
                }
                switch (event->get_seq_construction_role(segment)) {
                case legacy::SeqConstructionRole::Creates:
                case legacy::SeqConstructionRole::Fills: content_creators_[segment].push_back(event); break;
                case legacy::SeqConstructionRole::Modifies: content_modifiers_[segment].push_back(event); break;
                case legacy::SeqConstructionRole::None: break;
                }
            }
        }
    }

    const Events &end_creators(core::SeqTypeId segment, core::SeqSide side) const
    {
        return end_creators_.at(end(segment, side));
    }
    const Events &end_modifiers(core::SeqTypeId segment, core::SeqSide side) const
    {
        return end_modifiers_.at(end(segment, side));
    }
    const Events &content_creators(core::SeqTypeId segment) const { return content_creators_.at(segment); }
    const Events &content_modifiers(core::SeqTypeId segment) const { return content_modifiers_.at(segment); }

    /**
     * Does the segment carry a factor of its own in the downstream bound, apart from the gaps on
     * either side of it?
     *
     * It does when its content can still change against the read: an event modifies it (a
     * deletion trims it), or the event creating its content also places it (a gene choice aligns a
     * template, which can mismatch the read). Content filled into a slot another event placed (the
     * Markov chain's nucleotides, read off the read itself) has no factor of its own: its bound is
     * the bound of the gap it fills. TraversalPlan's key layout reads this, and nothing else does.
     */
    bool has_own_factor(core::SeqTypeId segment) const
    {
        if (not content_modifiers(segment).empty()) {
            return true;
        }
        for (const legacy::Rec_Event *creator : content_creators(segment)) {
            for (const core::SeqSide side : {core::Five_prime, core::Three_prime}) {
                if (creator->get_offset_role(segment, side) == legacy::OffsetRole::Creates) {
                    return true;
                }
            }
        }
        return false;
    }

private:
    static std::size_t end(core::SeqTypeId segment, core::SeqSide side)
    {
        return 2 * static_cast<std::size_t>(segment) + (side == core::Three_prime ? 1 : 0);
    }

    std::vector<Events> end_creators_;
    std::vector<Events> end_modifiers_;
    std::vector<Events> content_creators_;
    std::vector<Events> content_modifiers_;
};

} // namespace igor::model
