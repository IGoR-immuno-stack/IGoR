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

#include <cstddef>
#include <iterator>
#include <optional>
#include <vector>

/**
 * \brief What a SegmentSpan has instead of a length: every achievable distance, with the best
 * probability any completion of the scenario can reach at that distance.
 *
 * A **dense array indexed by distance**, not an associative container. The distances a span can
 * take are the sumset of its contributors' realization ranges -- deletions, insertions and one
 * gene template, each a contiguous run of integers -- so the key space is a short contiguous
 * interval with no holes worth naming. Measured over the two ends of the model range:
 *
 *     Insertion(VD), demo TRB model   [0, 30]      31 entries     248 B
 *     Gene_choice(J), V->J, demo TRB  [-57, 85]   143 entries     1.1 kB
 *     the same, human BCR-heavy       [-104, 154] 259 entries     2.1 kB
 *
 * Every profile in every model IGoR ships therefore fits in L1, which is what makes indexing
 * the right answer and a tree the wrong one. The `std::map` this replaced spent seven or eight
 * *dependent* compares per lookup to address 2 kB: `best_for` was **4.6 % of total runtime and
 * 8 % of iterate()**, of which 87 % was `_M_lower_bound` alone. See section 6.13 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md for the profile.
 *
 * The access pattern makes it better still. `Deletion` walks its realizations in decreasing
 * deletion count, so the distance it asks for **decreases by one per iteration**: the loop is a
 * unit-stride backwards scan over this array, which a tree expressed as N independent descents.
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
 * best_for() is the whole read interface, and it answers *value-or-absent* in one bounds check
 * and one load. The `count()`-then-`at()` pattern it originally replaced cost two red-black
 * descents of the same key at every scenario node -- 18 such sites across Gene_choice and
 * Deletion, plus one `at()` with no guard at all (Insertion), which threw where the others
 * discarded. See section 6.10, finding 6.
 */
class SpanProfile
{
public:
    /**
     * One achievable distance and the bound at it. A named pair because the call site reads
     * better for it -- `vd.distance` over `vd_len_iter->first`.
     */
    struct Entry {
        int distance;
        double proba;
    };

    /**
     * Forward iterator over the *present* entries, in increasing distance. Absent slots are
     * skipped, so a walk visits exactly size() entries whatever the array's extent.
     */
    class const_iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Entry;
        using difference_type = std::ptrdiff_t;
        using reference = Entry;
        using pointer = void;

        const_iterator() = default;
        const_iterator(const SpanProfile *owner, std::size_t index) : owner_(owner), index_(index)
        {
            skip_absent();
        }

        Entry operator*() const
        {
            return Entry{owner_->min_distance_ + static_cast<int>(index_), owner_->slots_[index_]};
        }

        const_iterator &operator++()
        {
            ++index_;
            skip_absent();
            return *this;
        }

        const_iterator operator++(int)
        {
            const_iterator before = *this;
            ++*this;
            return before;
        }

        friend bool operator==(const const_iterator &l, const const_iterator &r)
        {
            return l.index_ == r.index_;
        }
        friend bool operator!=(const const_iterator &l, const const_iterator &r) { return not(l == r); }

    private:
        void skip_absent()
        {
            while (index_ < owner_->slots_.size() and is_absent(owner_->slots_[index_])) {
                ++index_;
            }
        }

        const SpanProfile *owner_ = nullptr;
        std::size_t index_ = 0;
    };

    /// Drop every entry, for an owner about to re-fold at the next EM iteration. Keeps the
    /// allocation: the next iteration's range is the same one, since it depends on the model's
    /// realization sets and not on its probabilities.
    void clear()
    {
        slots_.clear();
        present_ = 0;
        min_distance_ = 0;
    }

    /**
     * \brief Offer \a proba as the bound at \a distance, keeping the better of the two.
     *
     * The fold's combine step -- max over scenarios reaching the same distance -- and the only
     * operation that can grow the array. Growth happens a handful of times per fold while the
     * range fills in and never again, against the ~10^6 leaves a fold visits, so the memmove a
     * front-extension costs is not worth designing around.
     *
     * \a proba must be non-negative, which it is by construction: it is a product of marginal
     * entries. A negative value would be indistinguishable from an empty slot.
     */
    void record(int distance, double proba)
    {
        make_room_for(distance);
        double &slot = slots_[static_cast<std::size_t>(distance - min_distance_)];
        if (is_absent(slot)) {
            slot = proba;
            ++present_;
        } else if (proba > slot) {
            slot = proba;
        }
    }

    /**
     * \brief The bound at \a distance, or nothing if no scenario reaches it -- in which case the
     * caller's branch is dead and must be discarded rather than scored.
     *
     * The single unsigned comparison catches both ends: a distance below min_distance_ makes the
     * subtraction negative, which wraps to a value no smaller than the array's extent.
     */
    std::optional<double> best_for(int distance) const
    {
        const std::size_t index = static_cast<std::size_t>(distance - min_distance_);
        if (index >= slots_.size()) {
            return std::nullopt;
        }
        const double proba = slots_[index];
        if (is_absent(proba)) {
            return std::nullopt;
        }
        return proba;
    }

    /// Present entries, not the array's extent -- the two differ if the achievable distances
    /// have a hole, which the contributors' contiguous ranges make unlikely but not impossible.
    std::size_t size() const { return present_; }
    bool empty() const { return present_ == 0; }

    /// Ordered by distance, absent slots skipped. For the consumers that enumerate rather than
    /// look up -- today only Gene_choice(D), building the retained decomposition of the span it
    /// splits.
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end() const { return const_iterator(this, slots_.size()); }

private:
    /// Distinguishable from any bound, which is a probability and so >= 0. The same convention
    /// SpanAccumulator uses for an unpublished length, and for the same reason: zero is a
    /// legitimate value here. A scenario whose bound is 0.0 is pruned on probability, where an
    /// absent one is discarded outright -- and in Deletion::iterate those take different exits,
    /// `break` against `continue`, so conflating them would not be bitwise.
    static constexpr double kAbsent = -1.0;
    static bool is_absent(double proba) { return proba < 0.0; }

    /// Extend the array so that \a distance is addressable, leaving new slots absent.
    void make_room_for(int distance)
    {
        if (slots_.empty()) {
            min_distance_ = distance;
            slots_.assign(1, kAbsent);
            return;
        }
        if (distance < min_distance_) {
            slots_.insert(slots_.begin(), static_cast<std::size_t>(min_distance_ - distance), kAbsent);
            min_distance_ = distance;
            return;
        }
        const std::size_t index = static_cast<std::size_t>(distance - min_distance_);
        if (index >= slots_.size()) {
            slots_.resize(index + 1, kAbsent);
        }
    }

    /// The distance slots_[0] stands for. Meaningless while slots_ is empty.
    int min_distance_ = 0;
    std::size_t present_ = 0;
    std::vector<double> slots_;
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
