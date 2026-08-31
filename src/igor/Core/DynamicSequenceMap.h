/*
 * DynamicSequenceMap.h
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

#include <igor/Core/IntStr.h>
#include <igor/Core/LayeredArray.h>
#include <igor/Core/SeqTypeRegistry.h>

#include <cstddef>
#include <stdexcept>

/**
 * \brief Whether a stored value counts as an "absent" segment for ordered traversal.
 *
 * Only sequence segments can be actively absent. Offsets, probabilities and mismatch
 * vectors are meaningful as soon as they are written, so the default says "never empty"
 * and the traversal falls back to the written/unwritten distinction alone.
 */
template <typename V>
struct SeqSegmentEmptiness
{
    static bool is_empty(const V &) { return false; }
};

/// A constructed sequence segment is absent when it is null or zero-length. The
/// zero-length case is the one that matters: it is how an event says "this segment was
/// processed and is not there", as opposed to "not processed yet" (layer -1).
template <>
struct SeqSegmentEmptiness<Int_Str *>
{
    static bool is_empty(Int_Str *const &value) { return value == nullptr || value->empty(); }
};

/**
 * \class DynamicSequenceMap DynamicSequenceMap.h
 * \brief A LayeredArray keyed by SeqTypeId, with ordered traversal over the registry.
 *
 * This is the runtime-sized replacement for the Enum_fast_memory_map instantiations that
 * are keyed by Seq_type (constructed sequences, offsets, mismatch lists, downstream proba
 * bounds). It adds exactly one thing to LayeredArray: the ability to walk the registry's
 * 5'->3' ordering and find the nearest neighbouring segment that is actually present.
 *
 * See docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md, decisions D1/D5 and task B10.
 *
 * ### Three states, not two
 *
 * Traversal depends on distinguishing:
 *   - **not yet processed** -- layer -1, `exists()` is false;
 *   - **actively absent**   -- written, but the value is an empty segment;
 *   - **present**           -- written and non-empty.
 *
 * Both of the first two are skipped, but they are not the same thing, and conflating them
 * is what would let a tandem-D model seed a Markov chain from the wrong nucleotide. See
 * `occupied()`.
 *
 * ### Lifetime
 *
 * Holds a reference to the registry, which must outlive the map and must be frozen before
 * construction -- the map sizes itself from the id space, so that space must not grow
 * underneath it. In inference the registry lives in the thread-local Model_Parms copy,
 * whose lifetime spans the parallel region.
 */
template <typename V>
class DynamicSequenceMap : private LayeredArray<V>
{
    using Base = LayeredArray<V>;

public:
    /// \param registry        frozen registry; must outlive this map
    /// \param initial_layers  layers to allocate up front; grows on demand
    /// \param count           number of keys, defaulting to the whole id space. Pass a
    ///                        smaller value for maps that must not address the tail of the
    ///                        id space (planned use: excluding flank types, task B3).
    explicit DynamicSequenceMap(const SeqTypeRegistry &registry,
                                std::size_t initial_layers = 1,
                                std::size_t count = 0)
        : Base(count != 0 ? count : registry.total_count(), initial_layers), registry_(registry)
    {
        if (!registry.is_frozen()) {
            throw std::logic_error(
                "DynamicSequenceMap: the registry must be frozen before a map is sized from it");
        }
    }

    // The container surface, keyed by SeqTypeId.
    using Base::count;
    using Base::current_layer;
    using Base::current_layers;
    using Base::exists;
    using Base::get;
    using Base::init_first_layer;
    using Base::layer;
    using Base::layer_capacity;
    using Base::multiply_all;
    using Base::request_layer;
    using Base::reset;
    using Base::restore_layer;
    using Base::set;
    using Base::set_current_layer;

    const SeqTypeRegistry &registry() const { return registry_; }

    /**
     * True when this seq_type has been written **and** holds a present segment.
     *
     * This is the predicate the traversal skips on, and the reason the two "skip" cases
     * must stay distinguishable: an unwritten key means the event has not run yet, while a
     * written empty one means it ran and produced nothing.
     */
    bool occupied(SeqTypeId type_id) const
    {
        return Base::exists(type_id) && !SeqSegmentEmptiness<V>::is_empty(Base::get(type_id));
    }

    /**
     * Nearest occupied seq_type strictly to the 5' side, walking the registry ordering.
     * Returns kNoSeqType when the ordering runs out.
     *
     * This is how an event finds its left-hand neighbour without knowing the topology: an
     * insertion asks for the segment bounding it, a Dinucl_markov for the segment holding
     * its seed nucleotide. Absent segments are skipped, so a tandem-D model missing its
     * second D reaches the first one with no conditional logic.
     */
    SeqTypeId first_occupied_left(SeqTypeId from) const
    {
        for (SeqTypeId cur = registry_.left_neighbor(from); cur != kNoSeqType;
             cur = registry_.left_neighbor(cur)) {
            if (occupied(cur)) {
                return cur;
            }
        }
        return kNoSeqType;
    }

    /// Nearest occupied seq_type strictly to the 3' side. See first_occupied_left().
    SeqTypeId first_occupied_right(SeqTypeId from) const
    {
        for (SeqTypeId cur = registry_.right_neighbor(from); cur != kNoSeqType;
             cur = registry_.right_neighbor(cur)) {
            if (occupied(cur)) {
                return cur;
            }
        }
        return kNoSeqType;
    }

private:
    const SeqTypeRegistry &registry_;
};
