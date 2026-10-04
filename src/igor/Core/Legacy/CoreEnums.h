/*
 * CoreEnums.h
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

#include <igor/Core/Types.h>

// CoreEnums.h was merged into igor/Core/Types.h (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).
// This stub keeps the legacy include path and the legacy names for the code that has not
// been promoted yet; it goes when its last consumer switches.
namespace igor::core::legacy {

using Event_type = igor::core::EventType;
using igor::core::GeneChoice_t;
using igor::core::Deletion_t;
using igor::core::Insertion_t;
using igor::core::Dinuclmarkov_t;
using igor::core::Undefined_t;

using Seq_side = igor::core::SeqSide;
using igor::core::Five_prime;
using igor::core::Three_prime;
using igor::core::Undefined_side;

using Seq_type = igor::core::SeqType;
using igor::core::V_gene_seq;
using igor::core::VD_ins_seq;
using igor::core::D_gene_seq;
using igor::core::DJ_ins_seq;
using igor::core::J_gene_seq;
using igor::core::VJ_ins_seq;

// The textual forms of these enums moved with them; legacy code calls them unqualified.
using igor::core::to_string;
using igor::core::operator<<;
using igor::core::operator+;

} // namespace igor::core::legacy
