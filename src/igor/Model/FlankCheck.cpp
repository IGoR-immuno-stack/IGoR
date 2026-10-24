/*
 * FlankCheck.cpp
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

#include <igor/Model/FlankCheck.h>

#include <igor/Model/TraversalPlan.h>

namespace igor::model {

using core::SeqSide;
using core::SeqTypeId;

std::vector<FlankCheck> resolve_flank_checks(const EventPlacement &placed, SeqSide side, SafetyMatrix &safety_set,
                                             const core::SeqOffsetsMap &seq_offsets)
{
    std::vector<FlankCheck> checks;
    const auto check_against = [&](SeqTypeId partner, bool this_is_left) {
        FlankCheck check;
        check.partner_id = partner;
        check.this_is_left = this_is_left;
        check.partner_side = this_is_left ? core::Five_prime : core::Three_prime;
        check.safety_cell = safety_set.cell(placed.target, partner);
        checks.push_back(check);
    };
    //The nearest placed neighbour on each side, and no other (ITER section 2.3): a pair further
    //away follows from the near ones, and a safe verdict written here propagates along the row
    //(SafetyMatrix::set()). The plan lists each side nearest first. The checks run 5' to 3'.
    if (side != core::Three_prime and not placed.placed_left.empty()) {
        check_against(placed.placed_left.front(), /*this_is_left=*/false);
    }
    if (side != core::Five_prime and not placed.placed_right.empty()) {
        check_against(placed.placed_right.front(), /*this_is_left=*/true);
    }

    //One claim per check: the 5' neighbour's pair sits in that neighbour's row, the 3' one's in
    //the target's, so the two never share a row.
    for (FlankCheck &check : checks) {
        safety_set.request_layer(check.safety_cell);
        check.safety_layer = safety_set.claimed_layer(check.safety_cell);
        check.partner_offset_layer = seq_offsets.claimed_layer(check.partner_id, check.partner_side);
    }
    return checks;
}

} // namespace igor::model
