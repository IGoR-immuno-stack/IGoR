/*
 * FoldFrontier.h
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

#include <igor/Model/Legacy/UnfilledSegmentLengths.h>

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace igor::model::legacy {

/**
 * \brief One partial path of the junction-length fold, reduced to what a later participant can
 * read from it.
 *
 * Two paths that agree on all three fields are interchangeable for everything after them: every
 * later participant multiplies both by the same factors and adds the same lengths. That is what
 * lets the fold keep the better of the two and go on with one (see FoldFrontier).
 *
 * - `length`: the summed `length_delta()` of everything chosen so far -- the profile's key once
 *   every participant has chosen.
 * - `parent_offsets`: one entry per participant, where in its marginals a later participant reads,
 *   moved there by the realizations of the parents the fold has chosen (R6). An entry is cleared
 *   once its participant has read it, since nothing reads it again.
 * - `lengths`: what each segment's creator published on this path, for the participant whose
 *   factor depends on it -- a Dinucl_markov's `p^L`. Kept to the end of the fold, because the fold
 *   does not know which later participant reads which entry.
 */
struct FoldState {
    int length = 0;
    std::vector<int> parent_offsets;
    UnfilledSegmentLengths lengths;

    bool operator==(const FoldState &) const = default;
};

struct FoldStateHash {
    std::size_t operator()(const FoldState &state) const
    {
        std::size_t seed = std::hash<int>{}(state.length);
        const auto combine = [&seed](std::size_t value) {
            seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };
        for (const int offset : state.parent_offsets) {
            combine(std::hash<int>{}(offset));
        }
        combine(state.lengths.hash());
        return seed;
    }
};

/**
 * \brief Every partial path of a fold that has reached the same participant, merged by state.
 *
 * Per state only the best probability is kept, which loses nothing: a later participant
 * multiplies every path in a state by the same factor, and rounding a product to the nearest
 * double never reverses the order of two non-negative values, so the best path in a state stays
 * the best after any number of further factors. The profile the fold ends with is therefore the
 * one enumerating every path would produce, bit for bit, at the cost of the number of states
 * rather than the number of paths.
 *
 * A probability of 0 is kept like any other: SpanProfile::record() marks a length as reachable
 * whatever the probability recorded there.
 */
class FoldFrontier
{
public:
    void offer(FoldState state, double proba)
    {
        //try_emplace leaves `state` untouched when its key is already there.
        const auto [slot, inserted] = best_.try_emplace(std::move(state), proba);
        if (not inserted and proba > slot->second) {
            slot->second = proba;
        }
    }

    std::size_t size() const { return best_.size(); }
    auto begin() const { return best_.begin(); }
    auto end() const { return best_.end(); }

private:
    std::unordered_map<FoldState, double, FoldStateHash> best_;
};

} // namespace igor::model::legacy
