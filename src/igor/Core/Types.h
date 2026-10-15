/*
 * Types.h
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

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>

#include <igor/Core/Export.h>

/**
 * \file Types.h
 * \brief The dependency-free vocabulary of IGoR: the enums describing events and sequence
 *        segments, and the small typedefs every layer names.
 *
 * Merged from the legacy CoreEnums.h, StdTypedefs.h and Typedef.h (step 1c of
 * doc/LAYER_REFACTORING_PROPOSAL.md). Nothing here depends on any other IGoR header. The
 * Gene_class enums stay in Legacy/Utils.h for now because they carry export annotations and a
 * cluster of conversion functions.
 */

namespace igor::core {

enum EventType { GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t, Undefined_t };

/// Which end of a constructed sequence segment an offset or a deletion refers to.
enum SeqSide { Five_prime = 0, Three_prime = 1, Undefined_side = 2 };

/**
 * The six sequence types of the legacy VDJ topology. Their values are the SeqTypeIds the
 * legacy registry assigns, which Model_Parms::read_model_parms() relies on, so enum-keyed code
 * still addresses the right slots.
 */
enum SeqType { V_gene_seq = 0, VD_ins_seq = 1, D_gene_seq = 2, DJ_ins_seq = 3, J_gene_seq = 4, VJ_ins_seq = 5 };

/// Name of a sequence type, as written in model files.
using SeqTypeString = std::string;

/// Type used as key for unordered maps, since an event cannot be instantiated.
using EventName = std::string;

/// Array of long doubles holding the marginal values.
using MarginalArrayPtr = std::unique_ptr<long double[]>;

/// Offset of an aligned sequence in the sequence_offsets maps. Characterizes the beginning
/// and the end of a sequence piece on the data sequence.
using SeqOffset = int;

using CodonTable = std::unordered_map<std::string, std::string>;

/// Index of a node, an event or a tensor dimension.
using index_type = std::int64_t;

/// Textual forms of the vocabulary enums, as written in model files and messages.
CORE_EXPORT SeqTypeString to_string(const SeqType);
CORE_EXPORT std::string to_string(const SeqSide);
CORE_EXPORT std::ostream &operator<<(std::ostream &, SeqSide);
CORE_EXPORT std::string operator+(const std::string &, SeqSide);
CORE_EXPORT std::string operator+(const std::string &, EventType);

} // namespace igor::core
