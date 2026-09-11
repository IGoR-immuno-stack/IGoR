/*
 * SpanProfile.h
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

#include <igor/Core/SegmentSpan.h>
#include <igor/Core/SeqTypeRegistry.h>

#include <map>
#include <optional>

/**
 * \brief What a SegmentSpan has instead of a length: every achievable distance, with the best
 * probability any completion of the scenario can reach at that distance.
 *
 * This is the `map<int,double>` the junction-length fold has always built, named. The naming is
 * what lets the bound stop being addressed by the `Seq_type` enum: a profile no longer *is*
 * `vd_length_best_proba_map`, it is the profile of some span, held by whichever event reads it
 * (JunctionBound below).
 *
 * ### It belongs to a span *and* to a reading position
 *
 * Not to a span alone. A boundary is a slot that a `Creates` event writes and a `Modifies` event
 * moves, so the same span has a different profile at each point of the priority ordering -- each
 * owner folds itself plus its suffix, and the contributors therefore differ. Section 2.5 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md has the derivation and the worked V->D example. Hence one
 * profile per consuming event rather than one per span.
 *
 * ### Reading it
 *
 * best_for() is the whole read interface, and it answers *value-or-absent* in one descent. The
 * `count()`-then-`at()` pattern it replaces cost two red-black descents of the same key at every
 * scenario node -- 18 such sites across Gene_choice and Deletion, plus one `at()` with no guard
 * at all (Insertion), which threw where the others discarded. See section 6.10, finding 6.
 */
class SpanProfile
{
public:
    using const_iterator = std::map<int, double>::const_iterator;

    /// Drop every entry, for an owner about to re-fold at the next EM iteration.
    void clear() { best_proba_by_distance_.clear(); }

    /**
     * \brief Offer \a proba as the bound at distance \a distance, keeping the better of the two.
     *
     * The fold's combine step -- max over scenarios reaching the same distance. One descent,
     * against the up-to-three the `count` / `at` / `operator[]` form took.
     */
    void record(int distance, double proba)
    {
        const auto [entry, inserted] = best_proba_by_distance_.try_emplace(distance, proba);
        if (not inserted and proba > entry->second) {
            entry->second = proba;
        }
    }

    /// The bound at \a distance, or nothing if no scenario reaches it -- in which case the
    /// caller's branch is dead and must be discarded rather than scored.
    std::optional<double> best_for(int distance) const
    {
        const auto entry = best_proba_by_distance_.find(distance);
        if (entry == best_proba_by_distance_.end()) {
            return std::nullopt;
        }
        return entry->second;
    }

    bool empty() const { return best_proba_by_distance_.empty(); }
    std::size_t size() const { return best_proba_by_distance_.size(); }

    /// Ordered by distance. For the consumers that enumerate rather than look up -- today only
    /// Gene_choice(D), building the retained decomposition of the span it splits.
    const_iterator begin() const { return best_proba_by_distance_.begin(); }
    const_iterator end() const { return best_proba_by_distance_.end(); }

private:
    std::map<int, double> best_proba_by_distance_;
};

/**
 * \brief One junction an event reads a bound from: which span, where its value goes, and the
 * profile itself.
 *
 * The replacement for the six `{vd,vj,dj}_length_best_proba_map` members and the three
 * `memory_layer_proba_map_junction*` scalars that went with them (decision O8). Everything about
 * *which* junction is resolved once, in `initialize_event()`; `iterate()` dereferences and never
 * asks a span question. That is the point: at 10^8-10^10 scenario nodes a span lookup would be a
 * pessimisation, so `(span, consumer)` is an **identity** here, not a runtime key.
 *
 * It also removes the enum ceiling. An event holds a fixed small number of these, whatever the
 * topology: nothing enumerates VD / DJ / VJ, so a tandem-D model needs no new member and no new
 * `Seq_type` value for `D1->D2`.
 *
 * \a proba_key and \a memory_layer address the slot of the downstream-bound array the value is
 * written to. Both were already resolved at init -- the layer explicitly, the key implicitly by
 * the `if (d_chosen) ... else if (j_chosen)` branch at each call site. This just keeps them
 * together with the profile they belong to.
 */
class JunctionBound
{
public:
    /// Whether initialize_Len_proba_bound() folds this junction's profile. A junction an event
    /// *splits* rather than measures carries no profile of its own: Gene_choice(D) writes the
    /// neutral 1.0 into the V->J slot and refines the two halves instead, which is the retained
    /// decomposition 5b generalises.
    enum class Fold { No, Yes };

    JunctionBound() = default;

    void resolve(SegmentSpan span, SeqTypeId proba_key, int memory_layer, Fold fold)
    {
        span_ = span;
        proba_key_ = proba_key;
        memory_layer_ = memory_layer;
        fold_ = fold;
    }

    /// False for a slot this event has no junction on -- a V gene choice has no left junction,
    /// and in a VJ model no event has a D-flanked one.
    bool resolved() const { return proba_key_ != kNoSeqType; }
    bool folded() const { return fold_ == Fold::Yes; }

    SegmentSpan span() const { return span_; }
    SeqTypeId proba_key() const { return proba_key_; }
    int memory_layer() const { return memory_layer_; }

    const SpanProfile &profile() const { return profile_; }
    SpanProfile &mutable_profile() { return profile_; }

private:
    SegmentSpan span_{};
    SeqTypeId proba_key_ = kNoSeqType;
    int memory_layer_ = -1;
    Fold fold_ = Fold::No;
    SpanProfile profile_{};
};
