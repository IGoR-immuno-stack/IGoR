/*
 * SeqOffsetsMap.h
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

#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/CoreEnums.h>
#include <igor/Core/StdTypedefs.h>

/**
 * \class Seq_offsets_map SeqOffsetsMap.h
 * \brief The 5' and 3' offsets of each constructed sequence segment.
 *
 * Replaces Enum_fast_memory_dual_key_map<Seq_type, Seq_side, Seq_Offset>. The dual-key map
 * addressed a single array as `key1 + range_key1 * key2`, with range_key1 fixed at the six
 * Seq_type enum values -- which is exactly what a model with more sequence types cannot
 * use. Storage is now **two independent single-key maps**, each runtime-sized from the
 * registry, so Seq_side stops being a dimension of the address space and becomes the
 * identity of the map you are addressing.
 *
 * Both are public: code that knows its side at compile time should say
 * `five_prime.get(id)` and skip the dispatch. The dual-key methods below are kept because
 * they let the ~100 existing call sites move over without being rewritten, and because
 * ScenarioContext's wrappers genuinely take the side as a parameter. As Deletion,
 * Gene_choice and Insertion are rewritten (B5/B6/B11) their sites become direct, and what
 * remains of this API is the handful that dispatch for real.
 *
 * Undefined_side is not a valid offset key -- an offset is one end of a segment or the
 * other -- and was never used as one even though the old map allocated a row for it.
 */
class Seq_offsets_map
{
public:
    explicit Seq_offsets_map(const SeqTypeRegistry &registry, std::size_t initial_layers = 1)
        : five_prime(registry, initial_layers), three_prime(registry, initial_layers)
    { }

    /// The map for one end. \throws std::out_of_range for Undefined_side.
    DynamicSequenceMap<Seq_Offset> &side(Seq_side seq_side)
    {
        return const_cast<DynamicSequenceMap<Seq_Offset> &>(
            static_cast<const Seq_offsets_map *>(this)->side(seq_side));
    }

    const DynamicSequenceMap<Seq_Offset> &side(Seq_side seq_side) const
    {
        switch (seq_side) {
        case Five_prime:
            return five_prime;
        case Three_prime:
            return three_prime;
        default:
            throw std::out_of_range("Seq_offsets_map: an offset must be Five_prime or Three_prime");
        }
    }

    Seq_Offset get(SeqTypeId type_id, Seq_side seq_side) const { return side(seq_side).get(type_id); }

    Seq_Offset get(SeqTypeId type_id, Seq_side seq_side, std::size_t layer) const
    {
        return side(seq_side).get(type_id, layer);
    }

    void set(SeqTypeId type_id, Seq_side seq_side, Seq_Offset offset, std::size_t layer)
    {
        side(seq_side).set(type_id, offset, layer);
    }

    void set_current(SeqTypeId type_id, Seq_side seq_side, Seq_Offset offset)
    {
        side(seq_side).set_current(type_id, offset);
    }

    bool exists(SeqTypeId type_id, Seq_side seq_side) const { return side(seq_side).exists(type_id); }

    int claimed_layer(SeqTypeId type_id, Seq_side seq_side) const
    {
        return side(seq_side).claimed_layer(type_id);
    }

    void request_layer(SeqTypeId type_id, Seq_side seq_side) { side(seq_side).request_layer(type_id); }

    void restore_layer(SeqTypeId type_id, Seq_side seq_side) { side(seq_side).restore_layer(type_id); }

    DynamicSequenceMap<Seq_Offset> five_prime;
    DynamicSequenceMap<Seq_Offset> three_prime;
};
