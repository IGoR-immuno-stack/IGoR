/*
 * StdTypedefs.h
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

/**
 * \file StdTypedefs.h
 * \brief IGoR typedefs whose definitions involve only standard-library types.
 *
 * Split out of Utils.h so that low-level headers can name these types without pulling in
 * Utils.h wholesale -- SeqTypeRegistry.h needs Seq_type_String, and including Utils.h for
 * it creates a cycle (Utils.h -> DynamicSequenceMap.h -> SeqTypeRegistry.h).
 *
 * The dividing line is deliberate: everything here depends only on <string>, <memory> and
 * friends. Typedefs that name IGoR's own classes (Int_Str_ptr, Next_event_ptr, the layered
 * map instantiations) stay with those classes, since they cannot be used without them
 * anyway. Utils.h includes this header, so existing consumers are unaffected.
 */

#include <memory>
#include <string>
#include <unordered_map>

/// Name of a sequence type (e.g. "V_gene_seq", "D1_gene_seq"). The serialization identity
/// of a sequence type; SeqTypeRegistry maps it to the SeqTypeId used at runtime.
using Seq_type_String = std::string;

/// Type used as key for unordered maps, since Rec_Event cannot be instantiated.
typedef std::string Rec_Event_name;

/// Array of long doubles holding the marginal values.
typedef std::unique_ptr<long double[]> Marginal_array_p;

/// Offset of an aligned sequence in the sequence_offsets maps. Characterizes the beginning
/// and the end of a sequence piece on the data sequence.
typedef int Seq_Offset;

typedef std::unordered_map<std::string, std::string> UMCodonTable;
