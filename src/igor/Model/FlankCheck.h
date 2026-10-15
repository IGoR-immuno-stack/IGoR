/*
 * FlankCheck.h
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

#include <igor/Core/SeqOffsetsMap.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Types.h>
#include <igor/Model/JunctionGeometry.h>
#include <igor/Model/SafetyMatrix.h>

#include <igor/Model/Export.h>

#include <vector>

namespace igor::model {

struct EventPlacement;

/**
 * \brief One neighbouring segment end that a placement or a deletion is checked against for
 * overlap.
 *
 * Resolved once, by resolve_flank_checks() in initialize_event(); iterate() reads it.
 * `Gene_choice` checks the neighbours on both sides of the segment it places, `Deletion` those on
 * the side it trims. What either reads per scenario stays with it.
 */
struct FlankCheck {
    core::SeqTypeId partner_id = core::kNoSeqType;
    bool this_is_left = true;                       ///< the partner sits 3' of the checked segment
    core::SeqSide partner_side = core::Five_prime;  ///< the partner end facing the checked segment
    SafetyCell safety_cell;                         ///< the pair, by ordering position (S5)
    int safety_layer = -1;                          ///< one per *row*, so same-row checks share it
    int partner_offset_layer = -1;                  ///< the layer the partner's facing end is read at

    /// Can the checked segment's facing end, which can still reach `mine`, avoid the partner's,
    /// which can still reach `theirs`? The one predicate of section 2.2, oriented by the side the
    /// partner is on.
    JunctionGeometry::Overlap overlap(JunctionGeometry::OffsetInterval mine,
                                      JunctionGeometry::OffsetInterval theirs) const
    {
        return this_is_left ? JunctionGeometry::check_overlap(mine, theirs, 0)
                            : JunctionGeometry::check_overlap(theirs, mine, 0);
    }
};

/**
 * The checks of `placed.target` against the neighbours whose facing end is placed when its event
 * runs, on `side` of it, or on both sides for `Undefined_side`.
 *
 * - **Order:** 5' to 3'. In a VDJ model, a D gene choice gets the former
 *   `memory_layer_safety_1` / `_2` pair in that order, which keeps the claimed layers where they
 *   were.
 * - **Safety layers:** one claim per *row*, not per check. Two neighbours on the same side name
 *   two cells of one row, which share a word and therefore a layer. Claiming twice would leave the
 *   second layer unwritten, and an unwritten layer is unreadable by design (section 7.9). A
 *   downstream deletion checking the same pair reads this layer as the one before its own.
 * - **Partner offsets:** each check records the layer its partner's facing end was claimed at.
 */
MODEL_EXPORT std::vector<FlankCheck> resolve_flank_checks(const EventPlacement &placed, core::SeqSide side,
                                                          SafetyMatrix &safety_set,
                                                          const core::SeqOffsetsMap &seq_offsets);

} // namespace igor::model
