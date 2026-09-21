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

#include <algorithm>
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
    /// look up -- today only Rec_Event::build_retained_decomposition(), composing two folded
    /// halves into the SpanDecomposition of the span between them.
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
 * \brief The *retained* decomposition of a span, for the event that splits it by enumerating
 * its own placements inside it.
 *
 * `SpanProfile` is `⊗ᵐᵃˣ`: it keeps, per total distance, the best probability any composition
 * reaches, and forgets which composition that was. This is `⊗ᵉⁿᵘᵐ`, the variant that keeps its
 * arguments — section 2.5 of docs/ITERATE_GENERIC_REWRITE_PLAN.md. It cannot be rebuilt from two
 * folded profiles after the fact, which is why it is built alongside them rather than derived
 * on demand.
 *
 * ### Three components, whatever the topology
 *
 * An anchor `X` enumerating exhaustively between anchors `A` and `B` branches on exactly two
 * things — which realization, and where its 5' end sits — and its 3' remainder is then
 * arithmetic. Everything past the next anchor is already marginalised into `span(X,B)` by the
 * max-fold, because the events out there enumerate for themselves when they run. So a placement
 * is `(realization, len(span(A,X)), len(span(X,B)))` plus its bound, and stays that whether one
 * anchor or four sit between `X` and `B`. Section 2.5 records the earlier reading — that a
 * tandem D1 would need five components — and why it was a conflation of the fold that *builds*
 * a span with the decomposition *retained* in it.
 *
 * ### Sorted by decreasing probability, and that is load-bearing
 *
 * The consumer walks one total distance's placements and `break`s on the first prune. That is
 * exact only because the order is non-increasing in the bound, so 5a's sections pin it
 * (§6.15) and nothing here may reorder without saying so.
 *
 * ### Storage
 *
 * A dense array indexed by total distance, for the same reason `SpanProfile` is one: the
 * achievable totals are the sumset of contiguous realization ranges, so the key space is a
 * short contiguous interval. The `std::map<int, std::vector<std::tuple<std::string,int,int,
 * double>>>` this replaces paid a red-black descent per lookup and a string hash per candidate
 * *inside* the enumeration; the realization is an index here, resolved against the event's own
 * dense table.
 */
class SpanDecomposition
{
public:
    /// One way the enumerating segment can sit inside the span.
    struct Placement {
        int realization_index;   ///< index into the event's realizations, not a name
        int left_distance;       ///< len(span(A, X))
        int right_distance;      ///< len(span(X, B))
        double proba;            ///< the bound this composition reaches
    };

    /// Drop everything, for an owner about to rebuild at the next EM iteration.
    void clear()
    {
        buckets_.clear();
        min_total_ = 0;
        empty_ = true;
    }

    /// Add one composition at `total`. Order of calls is preserved within a bucket, which is
    /// what makes the sort below reproduce the one it replaces exactly: the comparator looks
    /// only at the probability, so ties keep their arrival order up to std::sort's permutation.
    void record(int total, Placement placement)
    {
        make_room_for(total);
        buckets_[static_cast<std::size_t>(total - min_total_)].push_back(placement);
        empty_ = false;
    }

    /// Order every bucket by decreasing probability. Called once, after the last record().
    void sort_by_decreasing_proba()
    {
        for (std::vector<Placement> &bucket : buckets_) {
            std::sort(bucket.begin(), bucket.end(),
                      [](const Placement &l, const Placement &r) { return l.proba > r.proba; });
        }
    }

    /// The placements reaching `total`, best first; empty when no composition reaches it.
    const std::vector<Placement> &at(int total) const
    {
        static const std::vector<Placement> kNone;
        const std::size_t index = static_cast<std::size_t>(total - min_total_);
        return index < buckets_.size() ? buckets_[index] : kNone;
    }

    /// True while nothing has been recorded -- the event has no retained decomposition, either
    /// because it does not enumerate or because one of its two halves is unreachable.
    bool empty() const { return empty_; }

private:
    void make_room_for(int total)
    {
        if (buckets_.empty()) {
            min_total_ = total;
            buckets_.resize(1);
            return;
        }
        if (total < min_total_) {
            buckets_.insert(buckets_.begin(), static_cast<std::size_t>(min_total_ - total),
                            std::vector<Placement>{});
            min_total_ = total;
            return;
        }
        const std::size_t index = static_cast<std::size_t>(total - min_total_);
        if (index >= buckets_.size()) {
            buckets_.resize(index + 1);
        }
    }

    int min_total_ = 0;
    bool empty_ = true;
    std::vector<std::vector<Placement>> buckets_;
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
    /// What initialize_Len_proba_bound() builds for this junction.
    ///
    ///  - `Yes`    -- fold the max-product profile, the pruning bound every consumer reads;
    ///  - `No`     -- nothing. The event measures this junction but does not bound it: it
    ///                writes the neutral 1.0 into the slot and lets its halves carry the bound;
    ///  - `Retain` -- the decomposition rather than the fold. For an event that *splits* this
    ///                junction by enumerating its own placements inside it, which needs to know
    ///                which composition reached each total and not only the best one. Built
    ///                from the two halves plus this event's realizations; see
    ///                SpanDecomposition and section 2.5.
    ///
    /// `Retain` also writes the neutral 1.0 into the slot at scenario time, exactly as `No`
    /// does -- the two differ in what initialization builds, not in what iterate() writes.
    enum class Fold { No, Yes, Retain };

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
    bool retained() const { return fold_ == Fold::Retain; }

    SegmentSpan span() const { return span_; }
    SeqTypeId proba_key() const { return proba_key_; }
    int memory_layer() const { return memory_layer_; }

    const SpanProfile &profile() const { return profile_; }
    SpanProfile &mutable_profile() { return profile_; }

    /// Only meaningful for a `Retain` junction; empty for every other.
    const SpanDecomposition &decomposition() const { return decomposition_; }
    SpanDecomposition &mutable_decomposition() { return decomposition_; }

private:
    SegmentSpan span_{};
    SeqTypeId proba_key_ = kNoSeqType;
    int memory_layer_ = -1;
    Fold fold_ = Fold::No;
    SpanProfile profile_{};
    SpanDecomposition decomposition_{};
};
