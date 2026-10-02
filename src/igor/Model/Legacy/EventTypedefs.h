/*
 * EventTypedefs.h
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

// The typedefs that name Rec_Event. They used to live in Core's Utils.h behind a forward
// declaration; since step 1b of doc/LAYER_REFACTORING_PROPOSAL.md Core names no Model type,
// so they live here, next to the class they point to, in a header small enough to be included
// by the contexts without pulling Rec_Event.h in.

#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>

#include <igor/Core/Legacy/CoreEnums.h>
#include <igor/Core/Legacy/StdTypedefs.h>
#include <igor/Core/Legacy/Utils.h>


namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

class Rec_Event;

//Typedef used for getting the next event ptr
//typedef std::shared_ptr<Rec_Event> Next_event_ptr; //Does not work for some reason
typedef Rec_Event *Next_event_ptr;

// v2.0 events map: keyed by (Event_type, seq_type string, Seq_side) so that
// multiple events of the same type and side but different seq_types (e.g. two
// D-gene deletions on D1 vs D2 in a tandem-D model) can coexist unambiguously.
typedef std::unordered_map<std::tuple<Event_type, Seq_type_String, Seq_side>,
                           std::shared_ptr<Rec_Event>>
        Events_map;

struct inverse_offset_comparator
{
    bool operator()(const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_1,
                    const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_2)
    {
        return inv_offset_1.second < inv_offset_2.second;
    }
};

} // namespace igor::model::legacy
