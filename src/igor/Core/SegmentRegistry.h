/*
 * SegmentRegistry.h
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

//StdTypedefs.h rather than Utils.h: Utils.h pulls in DynamicSequenceMap.h, which needs
//this header, so including Utils.h here would be circular.
#include <igor/Core/Types.h>

#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * \brief A process-wide frozen registry holding just the six legacy segment names.
 *
 * Only tests use it, as a fixture for the VDJ layout. It was written for the junction length
 * bound precomputation, which named VD_ins_seq and DJ_ins_seq directly and has read the model's
 * own registry since B6/B7.
 *
 * Do not use it for anything that must respect a model's actual segments.
 */

namespace igor::core {

class SegmentRegistry;
const SegmentRegistry &legacy_seq_type_registry();

/**
 * \brief Dense runtime handle for a registered segment.
 *
 * Ids are assigned consecutively from 0 as types are registered, so they index arrays
 * directly -- this is what lets the per-scenario maps be flat arrays rather than hash
 * lookups. The name remains the serialization identity (model files, event names,
 * diagnostics); the id is the runtime handle. See decision D1/D5 in
 * docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md.
 */
using SegmentId = std::uint16_t;

/// Returned by the neighbour lookups at the ends of the ordering, and held by an event
/// whose segment has not been resolved against a registry yet.
inline constexpr SegmentId kNoSegment = std::numeric_limits<SegmentId>::max();

/**
 * \class SegmentRegistry SegmentRegistry.h
 * \brief The segments of a model: their names, their dense ids, and their ordering.
 *
 * Stores the left-to-right (5' to 3') ordering of sequence segment identifiers
 * (e.g. V_gene_seq, VD_ins_seq, D_gene_seq, ...) and provides neighbor lookup
 * for linked-list-style traversal used by the DinucMarkov model.
 *
 * For legacy (v1.x) model files the ordering is inferred automatically by
 * Model_Parms::finalize() and stored here so the rest of the code
 * has a single authoritative source.
 */
class SegmentRegistry
{
public:
    SegmentRegistry() = default;

    /**
     * Replace the registry contents with a new ordered list of seq_type names.
     * The order must reflect the 5'→3' direction of the constructed sequence.
     */
    void set_ordered_types(const std::vector<SegmentName> &types)
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
     * SegmentId == Seq_type for all of them.
     *
     * This is what lets code still keyed by the Seq_type enum address a SegmentId-keyed
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
        register_type("V_gene_seq"); // V_gene_seq  == 0
        register_type("VD_ins_seq"); // VD_ins_seq  == 1
        register_type("D_gene_seq"); // D_gene_seq  == 2
        register_type("DJ_ins_seq"); // DJ_ins_seq  == 3
        register_type("J_gene_seq"); // J_gene_seq  == 4
        register_type("VJ_ins_seq"); // VJ_ins_seq  == 5
    }

    /**
     * Register a seq_type name, returning its id. Idempotent: registering a name that is
     * already known returns the existing id rather than allocating a new one.
     * \throws std::logic_error if the registry is frozen.
     * \throws std::length_error if the id space (65535 types) is exhausted.
     */
    SegmentId register_type(const SegmentName &name)
    {
        auto it = name_to_id_.find(name);
        if (it != name_to_id_.end()) {
            return it->second;
        }
        throw_if_frozen("register_type");
        if (id_to_name_.size() >= static_cast<std::size_t>(kNoSegment)) {
            throw std::length_error("SegmentRegistry: too many sequence types registered");
        }
        const auto id = static_cast<SegmentId>(id_to_name_.size());
        id_to_name_.push_back(name);
        name_to_id_.emplace(name, id);
        //Keep the neighbour tables id-indexable. A type registered outside set_ordered_types()
        //is not in the ordering, so it has no neighbours.
        left_.push_back(kNoSegment);
        right_.push_back(kNoSegment);
        return id;
    }

    /// Id of a registered name. \throws std::out_of_range if the name is unknown.
    SegmentId id(const SegmentName &name) const
    {
        auto it = name_to_id_.find(name);
        if (it == name_to_id_.end()) {
            throw std::out_of_range("Unknown seq_type in SegmentRegistry::id(): " + name);
        }
        return it->second;
    }

    /// Name behind an id. \throws std::out_of_range if the id was never allocated.
    const SegmentName &name(SegmentId type_id) const
    {
        if (type_id >= id_to_name_.size()) {
            throw std::out_of_range("Unknown SegmentId in SegmentRegistry::name(): " + std::to_string(type_id));
        }
        return id_to_name_[type_id];
    }

    /// Number of registered types; ids are exactly [0, total_count()).
    std::size_t total_count() const { return id_to_name_.size(); }

    /// The ordering as ids, 5'->3'. Types registered but absent from the ordering do not appear.
    const std::vector<SegmentId> &ordering() const { return ordering_; }

    /// Neighbour in the ordering, or kNoSegment at the ends / for a type not in the ordering.
    SegmentId left_neighbor(SegmentId type_id) const { return type_id < left_.size() ? left_[type_id] : kNoSegment; }
    SegmentId right_neighbor(SegmentId type_id) const { return type_id < right_.size() ? right_[type_id] : kNoSegment; }

    /**
     * Seal the registry. Required before the per-scenario maps are sized from it, so that
     * the id space cannot grow underneath a map that has already allocated for it.
     * Idempotent.
     */
    void freeze() { frozen_ = true; }
    bool is_frozen() const { return frozen_; }

    /** Return the ordered list of seq_type names (5'→3'). */
    const std::vector<SegmentName> &get_ordered_types() const { return ordered_seq_types; }

    /** Return true if the given seq_type is registered. */
    bool contains(const SegmentName &seq_type) const { return type_to_index.count(seq_type) > 0; }

    /**
     * Return the zero-based position of seq_type in the 5'→3' order.
     * \throws std::out_of_range if seq_type is not registered.
     */
    size_t index_of(const SegmentName &seq_type) const
    {
        auto it = type_to_index.find(seq_type);
        if (it == type_to_index.end()) {
            throw std::out_of_range("Unknown seq_type in SegmentRegistry: " + seq_type);
        }
        return it->second;
    }

    /**
     * Return the seq_type immediately to the left (5' side) of the given seq_type,
     * or std::nullopt if seq_type is the leftmost segment or is not registered.
     */
    std::optional<SegmentName> get_left_neighbor(const SegmentName &seq_type) const
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
    std::optional<SegmentName> get_right_neighbor(const SegmentName &seq_type) const
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
            throw std::logic_error(std::string("SegmentRegistry::") + what + "() called on a frozen registry");
        }
    }

    /// Precompute id-indexed neighbour tables so traversal is an array lookup, not a hash.
    void rebuild_neighbours()
    {
        left_.assign(id_to_name_.size(), kNoSegment);
        right_.assign(id_to_name_.size(), kNoSegment);
        for (std::size_t pos = 0; pos != ordering_.size(); ++pos) {
            const SegmentId here = ordering_[pos];
            if (pos > 0) {
                left_[here] = ordering_[pos - 1];
            }
            if (pos + 1 < ordering_.size()) {
                right_[here] = ordering_[pos + 1];
            }
        }
    }

    std::vector<SegmentName> ordered_seq_types; ///< ordering by name (5'->3')
    std::unordered_map<SegmentName, size_t> type_to_index; ///< name -> position in the ordering

    std::vector<SegmentName> id_to_name_; ///< id -> name; size == total_count()
    std::unordered_map<SegmentName, SegmentId> name_to_id_;
    std::vector<SegmentId> ordering_; ///< ordering by id (5'->3')
    std::vector<SegmentId> left_, right_; ///< id-indexed neighbour tables
    bool frozen_ = false;
};

inline const SegmentRegistry &legacy_seq_type_registry()
{
    static const SegmentRegistry registry = [] {
        SegmentRegistry r;
        r.register_legacy_seq_types();
        r.freeze();
        return r;
    }();
    return registry;
}

} // namespace igor::core
