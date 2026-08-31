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

#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * \brief Dense runtime handle for a registered sequence type.
 *
 * Ids are assigned consecutively from 0 as types are registered, so they index arrays
 * directly -- this is what lets the per-scenario maps be flat arrays rather than hash
 * lookups. The name remains the serialization identity (model files, event names,
 * diagnostics); the id is the runtime handle. See decision D1/D5 in
 * docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md.
 */
using SeqTypeId = std::uint16_t;

/// Returned by the neighbour lookups at the ends of the ordering, and held by an event
/// whose seq_type has not been resolved against a registry yet.
inline constexpr SeqTypeId kNoSeqType = std::numeric_limits<SeqTypeId>::max();

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
        throw_if_frozen("set_ordered_types");
        ordered_seq_types = types;
        type_to_index.clear();
        for (size_t i = 0; i < types.size(); ++i) {
            type_to_index[types[i]] = i;
        }

        //Assign ids in 5'->3' order, so that id order and ordering position coincide for
        //the common case of a registry built solely from an ordering.
        ordering_.clear();
        ordering_.reserve(types.size());
        for (const auto &name : types) {
            ordering_.push_back(register_type(name));
        }
        rebuild_neighbours();
    }

    /**
     * Pre-register the six legacy seq_types in Seq_type enum order, so that
     * SeqTypeId == Seq_type for all of them.
     *
     * This is what lets code still keyed by the Seq_type enum address a SeqTypeId-keyed
     * map correctly. Without it ids are assigned in ordering order, which happens to match
     * the enum for a VDJ model but not for a VJ one -- there J_gene_seq would land on id 2
     * while the enum says 4, silently aliasing two segments.
     *
     * Call before the ordering is applied. Idempotent, and harmless for models that use
     * only a subset: an unused standard type simply occupies an id nothing addresses.
     */
    void register_legacy_seq_types()
    {
        //Order matters: it must match the Seq_type enum declared in Utils.h.
        register_type("V_gene_seq");   // V_gene_seq  == 0
        register_type("VD_ins_seq");   // VD_ins_seq  == 1
        register_type("D_gene_seq");   // D_gene_seq  == 2
        register_type("DJ_ins_seq");   // DJ_ins_seq  == 3
        register_type("J_gene_seq");   // J_gene_seq  == 4
        register_type("VJ_ins_seq");   // VJ_ins_seq  == 5
    }

    /**
     * Register a seq_type name, returning its id. Idempotent: registering a name that is
     * already known returns the existing id rather than allocating a new one.
     * \throws std::logic_error if the registry is frozen.
     * \throws std::length_error if the id space (65535 types) is exhausted.
     */
    SeqTypeId register_type(const Seq_type_String &name)
    {
        auto it = name_to_id_.find(name);
        if (it != name_to_id_.end()) {
            return it->second;
        }
        throw_if_frozen("register_type");
        if (id_to_name_.size() >= static_cast<std::size_t>(kNoSeqType)) {
            throw std::length_error("SeqTypeRegistry: too many sequence types registered");
        }
        const auto id = static_cast<SeqTypeId>(id_to_name_.size());
        id_to_name_.push_back(name);
        name_to_id_.emplace(name, id);
        //Keep the neighbour tables id-indexable. A type registered outside set_ordered_types()
        //is not in the ordering, so it has no neighbours.
        left_.push_back(kNoSeqType);
        right_.push_back(kNoSeqType);
        return id;
    }

    /// Id of a registered name. \throws std::out_of_range if the name is unknown.
    SeqTypeId id(const Seq_type_String &name) const
    {
        auto it = name_to_id_.find(name);
        if (it == name_to_id_.end()) {
            throw std::out_of_range("Unknown seq_type in SeqTypeRegistry::id(): " + name);
        }
        return it->second;
    }

    /// Name behind an id. \throws std::out_of_range if the id was never allocated.
    const Seq_type_String &name(SeqTypeId type_id) const
    {
        if (type_id >= id_to_name_.size()) {
            throw std::out_of_range("Unknown SeqTypeId in SeqTypeRegistry::name(): "
                                    + std::to_string(type_id));
        }
        return id_to_name_[type_id];
    }

    /// Number of registered types; ids are exactly [0, total_count()).
    std::size_t total_count() const { return id_to_name_.size(); }

    /// The ordering as ids, 5'->3'. Types registered but absent from the ordering do not appear.
    const std::vector<SeqTypeId> &ordering() const { return ordering_; }

    /// Neighbour in the ordering, or kNoSeqType at the ends / for a type not in the ordering.
    SeqTypeId left_neighbor(SeqTypeId type_id) const
    {
        return type_id < left_.size() ? left_[type_id] : kNoSeqType;
    }
    SeqTypeId right_neighbor(SeqTypeId type_id) const
    {
        return type_id < right_.size() ? right_[type_id] : kNoSeqType;
    }

    /**
     * Seal the registry. Required before the per-scenario maps are sized from it, so that
     * the id space cannot grow underneath a map that has already allocated for it.
     * Idempotent.
     */
    void freeze() { frozen_ = true; }
    bool is_frozen() const { return frozen_; }

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
    void throw_if_frozen(const char *what) const
    {
        if (frozen_) {
            throw std::logic_error(std::string("SeqTypeRegistry::") + what
                                   + "() called on a frozen registry");
        }
    }

    /// Precompute id-indexed neighbour tables so traversal is an array lookup, not a hash.
    void rebuild_neighbours()
    {
        left_.assign(id_to_name_.size(), kNoSeqType);
        right_.assign(id_to_name_.size(), kNoSeqType);
        for (std::size_t pos = 0; pos != ordering_.size(); ++pos) {
            const SeqTypeId here = ordering_[pos];
            if (pos > 0) {
                left_[here] = ordering_[pos - 1];
            }
            if (pos + 1 < ordering_.size()) {
                right_[here] = ordering_[pos + 1];
            }
        }
    }

    std::vector<Seq_type_String> ordered_seq_types;              ///< ordering by name (5'->3')
    std::unordered_map<Seq_type_String, size_t> type_to_index;   ///< name -> position in the ordering

    std::vector<Seq_type_String> id_to_name_;                    ///< id -> name; size == total_count()
    std::unordered_map<Seq_type_String, SeqTypeId> name_to_id_;
    std::vector<SeqTypeId> ordering_;                            ///< ordering by id (5'->3')
    std::vector<SeqTypeId> left_, right_;                        ///< id-indexed neighbour tables
    bool frozen_ = false;
};
