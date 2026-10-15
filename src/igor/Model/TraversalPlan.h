/*
 * TraversalPlan.h
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

#include <igor/Core/SegmentSpan.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Types.h>
#include <igor/Model/SegmentIndex.h>

#include <igor/Model/Export.h>

#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace igor::model {

namespace legacy {
class Rec_Event;
}

/**
 * \brief What one event can lean on when it runs, in the order an engine runs the events.
 *
 * Two different questions, which coincide in every model shipped today and need not:
 *   - **placed**: segments whose end facing this event has been created. They are what a placement
 *     or a deletion is checked against for overlap.
 *   - **built**: segments whose content has been created. A *gap* is the stretch between two built
 *     segments with nothing built in between; its downstream bound is one factor, refined by every
 *     event that changes it. An insertion places a slot inside a gap without building it, so the
 *     gap stays one until the Markov chain fills the slot.
 */
struct MODEL_EXPORT EventPlacement {
    core::SeqTypeId target = core::kNoSeqType;
    std::vector<core::SeqTypeId> placed_left;  ///< 3' end created before this event, nearest first
    std::vector<core::SeqTypeId> placed_right; ///< 5' end created before this event, nearest first
    core::SeqTypeId built_left = core::kNoSeqType;  ///< nearest built content to the 5' side
    core::SeqTypeId built_right = core::kNoSeqType; ///< nearest built content to the 3' side
    bool target_built = false; ///< the target's content is built once this event has run

    /// The gap the target sits in, between the nearest built content on each side.
    std::optional<core::SegmentSpan> enclosing_gap() const;

    /// The gap from the target to the nearest built content on `side`, once the target is built.
    /// While it is not, the gap it sits in, which is the one it can change.
    std::optional<core::SegmentSpan> gap_toward(core::SeqSide side) const;
};

/**
 * \brief Segments that exist before the first event of an order runs.
 *
 * Empty for a whole model. An engine running part of one starts from such a state -- a cluster of
 * Phase D's, with its neighbours' boundary given -- and so does a test fixture that writes a
 * segment instead of running the event that would create it.
 */
struct MODEL_EXPORT InitialState {
    std::vector<core::SeqTypeId> placed; ///< both ends created, the content not built (a slot)
    std::vector<core::SeqTypeId> built;  ///< both ends created and the content built
};

/**
 * \brief The structure an engine reads off its own order of events: who is placed and built when
 * each event runs, which gaps need a downstream bound, and the key each one is held at.
 *
 * Built once per model and order (tier 0a of SEGMENT_DECOMPOSITION_REVIEW §13), from the order
 * alone: it needs nothing from the events but their A0 roles and R12's participation rule. The
 * events are told their placement before initialize_event(), as Model_Parms::finalize() tells
 * them their neighbours; they no longer search the model for them. Owned by GenModel's set-up
 * until LegacyDfsEngine exists.
 *
 * Placements are filed by event identifier, which Rec_Event::copy() keeps, so one plan serves
 * every copy of the model it was built from: GenModel builds it once, and the events of every
 * thread, in every EM iteration, read the same plan. It is never modified after construction, so
 * concurrent reads are safe. It keeps pointers to the events it was built from, which must outlive
 * it; only check_admissible() and index() use them.
 *
 * ### Keys
 *
 * A segment's own factor is held at its SeqTypeId, as before. A gap is held at the id of the one
 * segment strictly inside it when there is exactly one and it has no factor of its own (the slot
 * an insertion places, which then bounds the same gap at that key); otherwise at a key of its own,
 * past the registry's ids, in the order the walk first meets it. The first rule is a layout rule,
 * not a meaning: the bound is a product taken in key order, and it is what keeps every shipped
 * model's product, and so its pruning, unchanged.
 */
class MODEL_EXPORT TraversalPlan
{
public:
    TraversalPlan(const core::SeqTypeRegistry &registry, std::span<const legacy::Rec_Event *const> order,
                  InitialState initial = {});

    const SegmentIndex &index() const { return index_; }

    /// The placement of `event`, or of the event with its identifier in the model the plan was
    /// built from. \throws std::out_of_range for an identifier not in this plan's order.
    const EventPlacement &placement(const legacy::Rec_Event &event) const;

    /// \throws std::out_of_range for a span that is never an open gap needing a bound in this
    ///         order: the event asking and the plan disagree about the model's structure.
    core::SeqTypeId gap_key(core::SegmentSpan gap) const;

    /// Every gap that needs a bound, in the order the walk first meets it.
    const std::vector<core::SegmentSpan> &gaps() const { return gaps_; }

    /// Size of the downstream bound map: the registry's ids and the gaps' own keys.
    std::size_t bound_key_count() const { return bound_key_count_; }

    /**
     * \brief Refuse an order in which the events cannot be coherent.
     *
     * The checks an engine owns, in place of the ones events used to make on themselves:
     *   - every event's target is a segment of the ordering;
     *   - every end and every content is created by exactly one event;
     *   - an end or a content is modified only after it is created;
     *   - content filled into a slot comes after the slot's ends are created;
     *   - a slot comes after its neighbours' facing ends, which its own are derived from, and
     *     sits inside a gap, between built content on both sides.
     * \throws std::invalid_argument naming the event and the reason.
     */
    void check_admissible() const;

private:
    const core::SeqTypeRegistry &registry_;
    std::vector<const legacy::Rec_Event *> order_;
    InitialState initial_;
    SegmentIndex index_;
    std::unordered_map<int, EventPlacement> placements_; ///< by event identifier
    std::vector<core::SegmentSpan> gaps_;
    std::map<std::pair<core::SeqTypeId, core::SeqTypeId>, core::SeqTypeId> gap_keys_;
    std::size_t bound_key_count_ = 0;
};

} // namespace igor::model
