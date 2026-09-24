/*
 * UnfilledSegmentLengths.h
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
 * \brief How many nucleotides each segment's offsets imply that nobody has chosen yet, along
 * one path of the junction-length fold.
 *
 * Replaces the `Seq_type_str_p_map` the fold used to be handed, which existed for one reason:
 * `Insertion` stashed a dummy `Int_Str` of the right length in it so that `Dinucl_markov`,
 * running later in the same traversal, could read `->size()` back out and raise its
 * per-nucleotide probability to that power. That was the
 * *"TODO constructed sequences should not be used but it is useful to compute the dinucl
 * contribution"* on Rec_Event.cpp -- a whole sequence map carried through the fold to move one
 * integer between two events.
 *
 * ### Why the content is "still to be chosen" and not "the length"
 *
 * `Dinucl_markov` raises its per-nucleotide probability to the published number. That is right
 * for a junction and wrong for a gene template, and while only creators of an all-placeholder
 * segment published, the two coincided and nothing said so (O12). The published quantity is
 * therefore **the number of positions still to be chosen**: `n` for an insertion whose offsets
 * are placed and whose sequence does not exist yet, **0** for a gene template, whose
 * nucleotides its creator fixed.
 *
 * That phrasing is what lets `Dinucl_markov` state its requirement as a property of the state
 * -- *offsets placed, sequence not yet created* -- rather than as a lookup from itself to its
 * `Insertion`, which is the coupling this refactor exists to remove.
 *
 * The publisher follows from it: **whoever creates a segment's offsets but not its sequence**,
 * which is `get_offset_role` and `get_seq_construction_role` together and needs no new
 * capability. Each key still has exactly one writer on any path, and a published value is
 * still non-negative -- a `Deletion` modifies an offset rather than creating one, and
 * contributes its negative delta to the span total without ever touching this.
 *
 * ### Still keyed, and the keying is load-bearing
 *
 * On a V->J span both insertion/dinucl pairs participate at once -- `affects_proba_of` answers
 * true for `VD_ins_seq` and `DJ_ins_seq` alike when the junction is `VJ_ins_seq` -- so a single
 * unkeyed scalar would alias the two.
 *
 * Scoped to one fold, not to a scenario: the fold walks a tree of realizations, and an entry
 * is overwritten by the next realization of the same event rather than restored on backtrack,
 * because every event appears at most once on a path.
 */
class UnfilledSegmentLengths
{
public:
    UnfilledSegmentLengths() = default;

    /// Sized for `seq_type_count` ids, with nothing published anywhere.
    explicit UnfilledSegmentLengths(std::size_t seq_type_count) : lengths_(seq_type_count, kAbsent) {}

    /// Publish how many of `id`'s positions are still to be chosen on this path.
    void set(SeqTypeId id, int length)
    {
        if (id == kNoSeqType || static_cast<std::size_t>(id) >= lengths_.size()) {
            return; //an event whose seq_type never resolved publishes nothing
        }
        lengths_[id] = length;
    }

    /// Whether anyone has published a count for `id` on this path.
    bool has(SeqTypeId id) const
    {
        return id != kNoSeqType && static_cast<std::size_t>(id) < lengths_.size()
               && lengths_[id] != kAbsent;
    }

    /// The published count. Undefined unless has(id); callers guard, as the dinucl factor does.
    int length_of(SeqTypeId id) const { return lengths_[id]; }

    /// Forget every published count, for a fold that starts over.
    void reset() { lengths_.assign(lengths_.size(), kAbsent); }

private:
    /// Distinguishable from any real published length, which is a segment size and so >= 0.
    static constexpr int kAbsent = -1;

    std::vector<int> lengths_;
};
