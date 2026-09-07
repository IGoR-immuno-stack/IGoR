/*
 * JunctionGeometry.h
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

#include <igor/Core/CoreEnums.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/StdTypedefs.h>
#include <igor/Core/Utils.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

/**
 * \brief The geometry shared by every overlap and junction-length decision in iterate().
 *
 * Task S2/S3 of docs/ITERATE_GENERIC_REWRITE_PLAN.md. What lives here is the part of the
 * scenario geometry that does not depend on which event is asking: how far an end can still
 * travel, and what that implies for the segment between two ends.
 */
namespace JunctionGeometry {

/**
 * \brief How far every segment end can still move once the undecided events have played out.
 *
 * Built once per event at `initialize_event()` and read in the hot loop. It replaces the
 * eight `*_min_del` / `*_max_del` scalars that `Deletion` and `Gene_choice` each carry, and
 * the two 60-line blocks that fill them by looking up the V-3', D-5', D-3' and J-5' deletion
 * events by hand (plan section 2.1). Because it asks the capability queries rather than
 * naming deletion events, a topology with more than one deletion per end -- tandem D -- needs
 * no new code here.
 *
 * **What "pending" means.** An event that is already in `processed_events` has had its
 * realization fixed further up the scenario tree, so it can no longer move anything and is
 * skipped. Everything still undecided contributes its full realization range. Note that the
 * event calling `rebuild()` is *not* itself processed yet -- the base
 * `Rec_Event::initialize_event()` inserts it only afterwards -- so its own travel is included,
 * exactly as in the code this replaces.
 *
 * **Two distinct situations collapse to `{0, 0}`**: an end whose modifier is already processed,
 * and an end that has no modifier in the model at all. The legacy scalars collapse them too,
 * and no consumer distinguishes them: both mean "this end will not move again". If a consumer
 * ever needs the difference it belongs in a separate query, not in a third state here.
 *
 * Deltas from several events *sum*. For every topology IGoR supports today at most one
 * deletion bears on a given end, so the sum has a single non-zero term and the result is
 * bitwise what the legacy lookup produced.
 */
/**
 * \brief The inclusive range of read positions one segment end can still occupy.
 *
 * `lo == hi` means the end is pinned: every pending modifier bearing on it has been decided.
 */
struct OffsetInterval {
    Seq_Offset lo = 0;
    Seq_Offset hi = 0;

    bool operator==(const OffsetInterval &) const = default;
};

/// The verdict of check_overlap(): what is still possible for a pair of neighbouring ends.
enum class Overlap {
    Infeasible,  ///< no combination of pending realizations avoids the overlap -- discard
    Safe,        ///< every combination avoids it -- no downstream event need re-check
    Undetermined ///< some do and some do not -- the deciding event must check again
};

class PendingModifierBounds
{
public:
    PendingModifierBounds() = default;

    /// Sized for `seq_type_count` ids with nothing pending anywhere.
    explicit PendingModifierBounds(std::size_t seq_type_count)
        : offset_(seq_type_count * 2), length_(seq_type_count)
    {
    }

    /**
     * Accumulate the bounds over every event of `events_map` not in `processed_events`.
     *
     * Sized from `registry.total_count()`, so ids and `SeqTypeId` agree with every other
     * id-keyed map in the traversal. Previous contents are discarded: calling this twice is
     * the same as building it fresh.
     */
    void rebuild(const SeqTypeRegistry &registry, const Events_map &events_map,
                 const std::unordered_set<Rec_Event_name> &processed_events)
    {
        const std::size_t count = registry.total_count();
        offset_.assign(count * 2, OffsetDelta{});
        length_.assign(count, LengthContribution{});

        for (const auto &[key, event] : events_map) {
            (void)key;
            if (!event || processed_events.count(event->get_name()) != 0) {
                continue;
            }
            for (std::size_t id = 0; id != count; ++id) {
                const auto type_id = static_cast<SeqTypeId>(id);
                for (const Seq_side side : {Five_prime, Three_prime}) {
                    const OffsetDelta delta = event->get_offset_delta_bounds(type_id, side);
                    check_ordered(delta.min, delta.max, *event, registry, type_id,
                                  "get_offset_delta_bounds");
                    OffsetDelta &accumulated = offset_[end_index(id, side)];
                    accumulated.min += delta.min;
                    accumulated.max += delta.max;
                }
                const LengthContribution contribution = event->get_length_contribution(type_id);
                check_ordered(contribution.min, contribution.max, *event, registry, type_id,
                              "get_length_contribution");
                LengthContribution &accumulated = length_[id];
                accumulated.min += contribution.min;
                accumulated.max += contribution.max;
            }
        }
    }

    /// Signed range by which `(type_id, side)` can still be shifted. `{0, 0}` if it cannot.
    OffsetDelta offset_delta(SeqTypeId type_id, Seq_side side) const
    {
        check_id(type_id, "offset_delta");
        return offset_[end_index(type_id, checked_side(side))];
    }

    /// Signed range of nucleotides still to be added to or removed from `type_id`.
    LengthContribution length(SeqTypeId type_id) const
    {
        check_id(type_id, "length");
        return length_[type_id];
    }

    /**
     * The interval `(type_id, side)` can still reach from where it currently sits.
     *
     * The sign flip is already absorbed by the delta -- a 3' end whose delta is `{-4, 0}`
     * reaches `[current - 4, current]`, a 5' end whose delta is `{0, 3}` reaches
     * `[current, current + 3]`, and a palindromic range spans both ways -- so this is a plain
     * translation. It does *not* re-sort the two bounds: `OffsetDelta` is contractually
     * ordered, and rebuild() rejects a provider that breaks that rather than papering over it.
     *
     * Degenerates to the point `[current, current]` exactly when nothing pending bears on the
     * end -- which is what makes one predicate serve both callers: in `Gene_choice::iterate`
     * the segment's own deletion is still pending and this is a proper interval, while in
     * `Deletion::iterate` the event being consumed *is* that modifier, so its own end has
     * collapsed by the time the check runs.
     */
    OffsetInterval reachable(SeqTypeId type_id, Seq_side side, Seq_Offset current) const
    {
        const OffsetDelta delta = offset_delta(type_id, side);
        return {current + delta.min, current + delta.max};
    }

    /// Number of ids this instance was sized for; valid ids are `[0, seq_type_count())`.
    std::size_t seq_type_count() const noexcept { return length_.size(); }

private:
    static std::size_t end_index(std::size_t type_id, Seq_side side)
    {
        return type_id * 2 + static_cast<std::size_t>(side);
    }

    static Seq_side checked_side(Seq_side side)
    {
        if (side != Five_prime && side != Three_prime) {
            throw std::invalid_argument("PendingModifierBounds: a segment end must be Five_prime "
                                        "or Three_prime, got " + to_string(side));
        }
        return side;
    }

    /// A provider returning `min > max` would make every interval derived from it nonsense --
    /// and, being merely *narrow* rather than obviously wrong, would survive a long way. This
    /// runs once per event per id at initialization, never in the hot loop.
    static void check_ordered(int min, int max, const Rec_Event &event, const SeqTypeRegistry &registry,
                              SeqTypeId type_id, const char *query)
    {
        if (min > max) {
            throw std::logic_error("PendingModifierBounds: " + event.get_name() + "::" + query
                                   + "(" + registry.name(type_id) + ") returned an unordered range ["
                                   + std::to_string(min) + ", " + std::to_string(max)
                                   + "]; min must not exceed max");
        }
    }

    void check_id(SeqTypeId type_id, const char *what) const
    {
        if (static_cast<std::size_t>(type_id) >= length_.size()) {
            throw std::out_of_range("PendingModifierBounds::" + std::string(what) + "(): seq_type id "
                                    + std::to_string(type_id) + " outside the "
                                    + std::to_string(length_.size()) + " ids this was built for");
        }
    }

    /// Indexed [id * 2 + side]; the two ends of a segment move independently.
    std::vector<OffsetDelta> offset_;
    /// Indexed [id].
    std::vector<LengthContribution> length_;
};

/**
 * Can what sits between two neighbouring ends still fit?
 *
 * `left_three_prime` is the reachable interval of the 3' end of the left segment,
 * `right_five_prime` that of the 5' end of the right one, and `gap` the minimum total number
 * of nucleotides that must sit strictly between them. The geometric constraint is
 * `left_3' + gap < right_5'`, so:
 *
 * - the best case is the left end as far left as it goes against the right end as far right,
 *   and if even that violates the constraint no realization can satisfy it -- `Infeasible`;
 * - the worst case is the reverse, and if even that satisfies it nothing downstream can break
 *   it -- `Safe`;
 * - otherwise the outcome depends on realizations not yet drawn -- `Undetermined`.
 *
 * This one predicate replaces the twelve hand-written comparisons in `Gene_choice::iterate`
 * and `Deletion::iterate` (plan section 2.2). The two callers differ only in whether the
 * moving end's own interval has collapsed to a point, which `reachable()` handles.
 *
 * `gap` is 0 at every site today -- every in-between segment can be empty, insertions
 * included -- and passing it explicitly is what keeps the predicate correct once a tandem-D
 * ordering puts a gene segment between two checked ends. It must be derived from a *minimum*
 * length; tightening it is a modelling change, not a refactor (plan section 7.2).
 */
inline Overlap check_overlap(OffsetInterval left_three_prime, OffsetInterval right_five_prime, int gap)
{
    if (left_three_prime.lo + gap >= right_five_prime.hi) {
        return Overlap::Infeasible;
    }
    if (left_three_prime.hi + gap < right_five_prime.lo) {
        return Overlap::Safe;
    }
    return Overlap::Undetermined;
}

} // namespace JunctionGeometry
