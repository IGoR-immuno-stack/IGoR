/*
 * SpanAccumulator.h
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

#include <igor/Core/SeqTypeRegistry.h>

#include <cstddef>
#include <vector>

/**
 * \brief The segment lengths decided so far along one path of the junction-length fold.
 *
 * Replaces the `Seq_type_str_p_map` the fold used to be handed, which existed for one reason:
 * `Insertion` stashed a dummy `Int_Str` of the right length in it so that `Dinucl_markov`,
 * running later in the same traversal, could read `->size()` back out and raise its
 * per-nucleotide probability to that power. That was the
 * *"TODO constructed sequences should not be used but it is useful to compute the dinucl
 * contribution"* on Rec_Event.cpp -- a whole sequence map carried through the fold to move one
 * integer between two events.
 *
 * This carries the integer. Only events that **create** a segment publish here
 * (`SeqConstructionRole::Creates`), so each key has exactly one writer on any path and a
 * published length is always non-negative -- a `Deletion` contributes its negative delta to
 * the span total without ever touching this.
 *
 * Scoped to one fold, not to a scenario: the fold walks a tree of realizations, and an entry
 * is overwritten by the next realization of the same event rather than restored on backtrack,
 * because every event appears at most once on a path.
 */
class SpanAccumulator
{
public:
    SpanAccumulator() = default;

    /// Sized for `seq_type_count` ids, with nothing published anywhere.
    explicit SpanAccumulator(std::size_t seq_type_count) : lengths_(seq_type_count, kAbsent) {}

    /// Publish the length `id`'s creator just chose for it.
    void set(SeqTypeId id, int length)
    {
        if (id == kNoSeqType || static_cast<std::size_t>(id) >= lengths_.size()) {
            return; //an event whose seq_type never resolved publishes nothing
        }
        lengths_[id] = length;
    }

    /// Whether a creator has published a length for `id` on this path.
    bool has(SeqTypeId id) const
    {
        return id != kNoSeqType && static_cast<std::size_t>(id) < lengths_.size()
               && lengths_[id] != kAbsent;
    }

    /// The published length. Undefined unless has(id); callers guard, as the dinucl factor does.
    int length_of(SeqTypeId id) const { return lengths_[id]; }

    /// Forget every published length, for a fold that starts over.
    void reset() { lengths_.assign(lengths_.size(), kAbsent); }

private:
    /// Distinguishable from any real published length, which is a segment size and so >= 0.
    static constexpr int kAbsent = -1;

    std::vector<int> lengths_;
};
