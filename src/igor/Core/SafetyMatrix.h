/*
 * SafetyMatrix.h
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

#include <igor/Core/LayeredArray.h>
#include <igor/Core/SeqTypeRegistry.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

/// A position in the 5'->3' ordering that the safety matrix cannot address: a seq_type the
/// model registered but left out of the ordering, or a cell nobody resolved.
inline constexpr std::uint8_t kNoOrderingPosition = 0xFF;

/**
 * \brief One pair of segments in the overlap-safety matrix, named by ordering position.
 *
 * `row` is the 5'-most of the two, `column` the 3'-most, both as positions in the model's
 * 5'->3' ordering -- not as `SeqTypeId`s, whose numbering is the legacy enum's and puts
 * `VJ_ins_seq` at 5 in a VJ model whose ordering position is 1. Resolved once, at
 * `initialize_event()`, and read in the hot loop; `SafetyMatrix::cell()` is what resolves it.
 */
struct SafetyCell {
    std::uint8_t row = kNoOrderingPosition;
    std::uint8_t column = kNoOrderingPosition;

    /// False for a default-constructed cell -- an event whose partner the model does not have.
    bool resolved() const noexcept { return row != kNoOrderingPosition; }

    bool operator==(const SafetyCell &) const = default;
};

/**
 * \class SafetyMatrix SafetyMatrix.h
 * \brief Which pairs of segments are known not to overlap, for the scenario being explored.
 *
 * Task S5 of docs/ITERATE_GENERIC_REWRITE_PLAN.md, section 2.3. It replaces
 * `LayeredArray<bool>` keyed by the three-valued `Event_safety` enum, which could name only
 * the three unordered pairs a single-D model has and therefore put a ceiling on tandem D.
 *
 * ### Storage: one word per row, not one slot per pair
 *
 * Storage is still the *n(n-1)/2* cells of a triangular matrix -- keying by the left member
 * alone is not enough, because the nearest chosen neighbour of a segment changes with
 * recursion depth and the pair identity would become implicit in a chosen set the reader
 * cannot see (section 2.3, "Why one slot per row is not enough"). What collapses is the
 * *work*: the cells of a row live in one `std::uint32_t`, so marking a whole row suffix is
 * one bitwise OR rather than a loop, and the per-event layer count stays at two.
 *
 * The cost is a limit of 32 positions in the ordering -- 5 for VDJ, 7 for tandem D -- which
 * the constructor rejects rather than silently truncating. Section 2.3 proposed checking it
 * at `SeqTypeRegistry::freeze()`; it lives here instead, because the limit belongs to this
 * container and not to the registry, and because the bound is on the *ordering*, which is a
 * subset of what the registry holds.
 *
 * ### The row word is the whole state of the row
 *
 * A write is read-modify-write: it starts from the row's word at its own layer if it has
 * already written there, and from the layer below otherwise, so a cell nobody touched at
 * this depth keeps the value the enclosing depth left. That is what makes `layer - 1` --
 * the read every consumer performs -- still mean "what the previous writer of this pair
 * left", even though the previous writer of the *row* may have been touching another cell.
 *
 * The corollary is a promise: an event that claims a layer on a row must write that row at
 * that layer before it hands off, exactly as under the per-slot container (section 7.9).
 * A claimed-but-unwritten layer is not readable, and `get()` says so rather than serving a
 * default.
 *
 * ### Row-suffix propagation
 *
 * Establishing that a pair `(A, B)` cannot overlap also marks every `(A, C)` with `C` right
 * of `B` -- section 2.3's corollary. This is *not* a claim that those pairs are separated:
 * it is the observation that if a completed scenario violated `C(A, C)`, it would also
 * violate `C(B, C)`, and `(B, C)` is checked for real by whichever of `B`, `C` is chosen
 * second. So propagation moves *where* such a scenario is discarded, never *whether* --
 * which is what makes replacing the container a bitwise-exact step.
 *
 * The argument rests on one property of the model, not of this container: **every pair
 * adjacent in chosen order is checked by someone**. Propagation only ever writes cells that
 * are *not* nearest-neighbour at the moment they are written, so it cannot displace such a
 * check; validating the property for an arbitrary ordering is B5's and phase C's.
 *
 * Marking a pair *unsafe* propagates nothing: it says the deciding event must look, and says
 * nothing about pairs further away.
 */
class SafetyMatrix
{
public:
    /// Positions addressable by one row word. Section 2.3: 5 for VDJ, 7 for tandem D.
    static constexpr std::size_t kMaxPositions = 32;

    /// \param registry        frozen registry whose ordering supplies the positions
    /// \param initial_layers  layers to allocate up front; grows on demand
    explicit SafetyMatrix(const SeqTypeRegistry &registry, std::size_t initial_layers = 1)
        : rows_(registry.ordering().size(), initial_layers),
          position_(registry.total_count(), kNoOrderingPosition)
    {
        const std::vector<SeqTypeId> &ordering = registry.ordering();
        if (ordering.size() > kMaxPositions) {
            throw std::length_error("SafetyMatrix: the ordering has " + std::to_string(ordering.size())
                                    + " segments, and a row bitmask addresses at most "
                                    + std::to_string(kMaxPositions));
        }
        for (std::size_t pos = 0; pos != ordering.size(); ++pos) {
            position_[ordering[pos]] = static_cast<std::uint8_t>(pos);
        }
    }

    /// Number of rows -- one per position in the ordering. Named `count()` so that the same
    /// generic layer introspection serves this map and every other layered one.
    std::size_t count() const noexcept { return rows_.count(); }

    /// Whether this segment has a position in the ordering, and so can be half of a pair.
    /// False for a type the model registered but left out -- `VJ_ins_seq` in a VDJ model,
    /// `D_gene_seq` in a VJ one -- which is what a caller enumerating legacy pairs must ask
    /// before naming one.
    bool addresses(SeqTypeId type_id) const noexcept
    {
        return static_cast<std::size_t>(type_id) < position_.size()
               && position_[type_id] != kNoOrderingPosition;
    }

    /**
     * The cell naming the pair `{a, b}`, whichever way round it is given.
     *
     * \throws std::out_of_range if either id is outside the registry this was built for, or
     *         is not in its ordering.
     * \throws std::invalid_argument if the two are the same segment: a segment does not
     *         overlap itself, and asking is a resolution bug rather than a query.
     */
    SafetyCell cell(SeqTypeId a, SeqTypeId b) const
    {
        const std::uint8_t pa = position_of(a);
        const std::uint8_t pb = position_of(b);
        if (pa == pb) {
            throw std::invalid_argument("SafetyMatrix::cell(): a segment cannot overlap itself "
                                        "(ordering position " + std::to_string(pa) + ")");
        }
        return pa < pb ? SafetyCell{pa, pb} : SafetyCell{pb, pa};
    }

    /// Claim the next layer on this cell's row. Two cells of the same row share one claim:
    /// they share the word, so an event checking both its flanks against segments on the
    /// same side requests once and uses the granted layer for both.
    void request_layer(SafetyCell cell) { rows_.request_layer(checked(cell).row); }

    /// The highest layer this cell's row has been claimed at, or -1.
    int claimed_layer(SafetyCell cell) const { return rows_.claimed_layer(checked(cell).row); }
    /// Row-keyed form, for the generic layer introspection the tests run over every map.
    int claimed_layer(std::size_t row) const { return rows_.claimed_layer(row); }

    /// The layer this cell's row currently stands at -- the last one written -- or -1.
    int current_layer(SafetyCell cell) const { return rows_.current_layer(checked(cell).row); }
    int current_layer(std::size_t row) const { return rows_.current_layer(row); }

    /// True once this cell's row has been written at least once.
    bool exists(SafetyCell cell) const { return rows_.exists(checked(cell).row); }

    /// Whether the pair is established as non-overlapping, at the row's current layer.
    bool get(SafetyCell cell) const { return bit(rows_.get(checked(cell).row), cell.column); }

    /// Whether the pair is established as non-overlapping, at an explicit layer.
    /// \throws std::out_of_range if that layer was claimed but never written.
    bool get(SafetyCell cell, std::size_t layer) const
    {
        return bit(rows_.get(checked(cell).row, layer), cell.column);
    }

    /**
     * Record a verdict for the pair at `layer`, carrying the rest of the row forward.
     *
     * `true` marks the whole row suffix from this column on -- see "Row-suffix propagation"
     * above. `false` marks this column alone.
     */
    void set(SafetyCell cell, bool is_safe, std::size_t layer)
    {
        checked(cell);
        std::uint32_t word = seed(cell.row, layer);
        if (is_safe) {
            word |= suffix_mask(cell.column);
        } else {
            word &= ~column_mask(cell.column);
        }
        rows_.set(cell.row, word, layer);
    }

    /// Every row all-unsafe at layer 0, marking them written. Mirrors the other maps'
    /// initialization, and gives the first writer on a row a defined word to build on.
    void init_first_layer() { rows_.init_first_layer(0u); }

private:
    static constexpr std::uint32_t column_mask(std::uint8_t column) noexcept
    {
        return std::uint32_t{1} << column;
    }

    /// Every column from `column` on: the pair itself, and the corollary's row suffix. The
    /// two halves are one mask because "safe at k" and "therefore safe past k" are written
    /// together; splitting them would let a caller record one without the other.
    static constexpr std::uint32_t suffix_mask(std::uint8_t column) noexcept
    {
        return ~(column_mask(column) - 1);
    }

    static constexpr bool bit(std::uint32_t word, std::uint8_t column) noexcept
    {
        return (word & column_mask(column)) != 0;
    }

    std::uint8_t position_of(SeqTypeId type_id) const
    {
        if (static_cast<std::size_t>(type_id) >= position_.size()
            || position_[type_id] == kNoOrderingPosition) {
            throw std::out_of_range("SafetyMatrix: seq_type id " + std::to_string(type_id)
                                    + " is not in the ordering this matrix was built for");
        }
        return position_[type_id];
    }

    const SafetyCell &checked(const SafetyCell &cell) const
    {
        if (not cell.resolved() || cell.row >= rows_.count() || cell.column >= rows_.count()
            || cell.row >= cell.column) {
            throw std::out_of_range("SafetyMatrix: cell (" + std::to_string(cell.row) + ", "
                                    + std::to_string(cell.column) + ") is not a pair of this "
                                    "matrix's " + std::to_string(rows_.count()) + " positions");
        }
        return cell;
    }

    /**
     * The row word a write at `layer` starts from.
     *
     * Above the row's data mark -- nothing written at this depth yet -- it is the enclosing
     * depth's word, so untouched cells carry over. At or below it, it is the row's own word
     * at that layer, so a second cell written at the same depth joins the first instead of
     * erasing it. Layer 0 on a never-written row starts from all-unsafe.
     */
    std::uint32_t seed(std::uint8_t row, std::size_t layer) const
    {
        const int current = rows_.current_layer(row);
        if (current >= static_cast<int>(layer)) {
            return rows_.get(row, layer);
        }
        if (layer == 0) {
            return 0u;
        }
        return rows_.get(row, layer - 1);
    }

    LayeredArray<std::uint32_t> rows_;
    /// SeqTypeId -> ordering position, kNoOrderingPosition for a type outside the ordering.
    std::vector<std::uint8_t> position_;
};
