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
    //The plan lists each side nearest first: the 5' side is walked backwards, so the checks run
    //5' to 3'.
    if (side != core::Three_prime) {
        for (auto partner = placed.placed_left.rbegin(); partner != placed.placed_left.rend(); ++partner) {
            check_against(*partner, /*this_is_left=*/false);
        }
    }
    if (side != core::Five_prime) {
        for (const SeqTypeId partner : placed.placed_right) {
            check_against(partner, /*this_is_left=*/true);
        }
    }

    for (std::size_t i = 0; i != checks.size(); ++i) {
        FlankCheck &check = checks[i];
        std::size_t earlier = 0;
        while (earlier != i and checks[earlier].safety_cell.row != check.safety_cell.row) {
            ++earlier;
        }
        if (earlier != i) {
            check.safety_layer = checks[earlier].safety_layer;
        } else {
            safety_set.request_layer(check.safety_cell);
            check.safety_layer = safety_set.claimed_layer(check.safety_cell);
        }
        check.partner_offset_layer = seq_offsets.claimed_layer(check.partner_id, check.partner_side);
    }
    return checks;
}

} // namespace igor::model
