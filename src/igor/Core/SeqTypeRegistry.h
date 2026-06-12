/*
 * SeqTypeRegistry.h
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

#include <igor/Core/Utils.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * \class SeqTypeRegistry SeqTypeRegistry.h
 * \brief Registry for the ordered sequence types in a constructed sequence (v2.0 format).
 *
 * Stores the left-to-right (5' to 3') ordering of sequence segment identifiers
 * (e.g. V_gene_seq, VD_ins_seq, D_gene_seq, ...) and provides neighbor lookup
 * for linked-list-style traversal used by the DinucMarkov model.
 *
 * For legacy (v1.x) model files the ordering is inferred automatically by
 * Model_Parms::read_model_parms() and stored here so the rest of the code
 * has a single authoritative source.
 */
class SeqTypeRegistry
{
public:
    SeqTypeRegistry() = default;

    /**
     * Replace the registry contents with a new ordered list of seq_type names.
     * The order must reflect the 5'→3' direction of the constructed sequence.
     */
    void set_ordered_types(const std::vector<Seq_type_String> &types)
    {
        ordered_seq_types = types;
        type_to_index.clear();
        for (size_t i = 0; i < types.size(); ++i) {
            type_to_index[types[i]] = i;
        }
    }

    /** Return the ordered list of seq_type names (5'→3'). */
    const std::vector<Seq_type_String> &get_ordered_types() const { return ordered_seq_types; }

    /** Return true if the given seq_type is registered. */
    bool contains(const Seq_type_String &seq_type) const { return type_to_index.count(seq_type) > 0; }

    /**
     * Return the zero-based position of seq_type in the 5'→3' order.
     * \throws std::out_of_range if seq_type is not registered.
     */
    size_t index_of(const Seq_type_String &seq_type) const
    {
        auto it = type_to_index.find(seq_type);
        if (it == type_to_index.end()) {
            throw std::out_of_range("Unknown seq_type in SeqTypeRegistry: " + seq_type);
        }
        return it->second;
    }

    /**
     * Return the seq_type immediately to the left (5' side) of the given seq_type,
     * or std::nullopt if seq_type is the leftmost segment or is not registered.
     */
    std::optional<Seq_type_String> get_left_neighbor(const Seq_type_String &seq_type) const
    {
        auto it = type_to_index.find(seq_type);
        if (it == type_to_index.end() || it->second == 0) {
            return std::nullopt;
        }
        return ordered_seq_types[it->second - 1];
    }

    /**
     * Return the seq_type immediately to the right (3' side) of the given seq_type,
     * or std::nullopt if seq_type is the rightmost segment or is not registered.
     */
    std::optional<Seq_type_String> get_right_neighbor(const Seq_type_String &seq_type) const
    {
        auto it = type_to_index.find(seq_type);
        if (it == type_to_index.end() || it->second + 1 >= ordered_seq_types.size()) {
            return std::nullopt;
        }
        return ordered_seq_types[it->second + 1];
    }

    bool empty() const { return ordered_seq_types.empty(); }
    size_t size() const { return ordered_seq_types.size(); }

private:
    std::vector<Seq_type_String> ordered_seq_types;
    std::unordered_map<Seq_type_String, size_t> type_to_index;
};
