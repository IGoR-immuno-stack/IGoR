/*
 * IntStr.h
 *
 *  Created on: Jul 21, 2016
 *      Author: Quentin Marcou
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

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <igor/Core/IntStr.h>

// IntStr.h was promoted out of Legacy/ (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).
// This stub keeps the legacy include path and the legacy names for the code that has not
// been promoted yet; it goes when its last consumer switches.
namespace igor::core::legacy {

using Int_Str = igor::core::IntStr;

} // namespace igor::core::legacy
