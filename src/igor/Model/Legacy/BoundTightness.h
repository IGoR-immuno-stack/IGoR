/*
 * BoundTightness.h
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

/**
 * \file
 * \brief How far the pruning bound sits above the probability a scenario actually realizes,
 *        why the walk explored subtrees that produced nothing, and which factor of the bound
 *        is responsible.
 *
 * Every scenario reaches its leaf carrying an *upper bound* on what it could still be worth, and
 * the walk discards it when that bound falls below the threshold. How much work that saves, and
 * how much a **tighter** bound would save on top, is the question section 6.10 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md leaves open and R6 has to answer: conditioning the bound on
 * the scenario's already-chosen parents costs storage, and nobody has measured what it buys.
 *
 * Four distributions are reported. See docs/PROBA_BOUND_MACHINERY.md section 9 for how to read
 * them and what the first runs said.
 *
 * **At leaves.** The ratio `bound / realized` for each completed scenario. This is what section
 * 6.10 proposed, and measured on the TRB corpus it reads **100 % within a factor of 1.8** -- by
 * the leaf every entry of the downstream bound map has been replaced by the factor the scenario
 * actually realized, so the two quantities are the same product and the ratio is one by
 * construction. Kept because it is a real invariant check, not because it discriminates.
 *
 * **At internal nodes.** The ratio of the bound a node carried to the best probability any leaf
 * *below* it turned out to reach. This is the quantity pruning acts on and the one a tighter
 * bound would improve: a node whose bound sits three decades above everything in its subtree is
 * three decades of scenarios explored for nothing. Reported per depth, so it is visible whether
 * the looseness is concentrated near the root -- where it is expensive -- or near the leaves.
 *
 * **Why a barren node was barren.** Most of the walk reaches no leaf at all, and "barren" alone
 * does not say whether a better bound could have avoided it. A node whose children were all
 * *probability-tested and rejected* is one a tighter bound at the node itself would have deleted;
 * a node whose children were never probability-tested at all died on geometry, safety or read
 * bounds, and no probability bound can reach it -- only a feasibility test can. The two are
 * counted apart, because they point at different work.
 *
 * **Which factor is loose.** The aggregate ratio says the product over-estimates without saying
 * which of its factors does. The decomposition is per *event*, and it telescopes: the bound at a
 * node divided by the bound at its parent is exactly how much the bound tightened when that
 * node's event ran, and the product of those steps down the path to the best leaf is the parent's
 * whole over-estimate. A step of 1 means the parent's bound already knew what this event would
 * contribute; a step of 10^-3 means three decades of the parent's slack were this event's, and
 * conditioning *its* bound is what would remove them. Taken along the winning path only, so what
 * is measured is slack and not the ordinary cost of choosing a realization.
 *
 * Not per *segment slot*: an event sets its slot in the downstream bound map to 1.0 once its
 * segment is resolved and multiplies what it realized into the scenario probability instead, so
 * by the leaf every slot is 1 and a slot-by-slot ratio against the leaf compares a bound against
 * nothing. That is also why the leaf ratio is 1 by construction. This was measured before it was
 * believed; see docs/PROBA_BOUND_MACHINERY.md section 9.
 *
 * ### It is compiled out unless asked for
 *
 * Every entry point is an empty inline function unless `IGOR_BOUND_INSTRUMENTATION` is defined, so
 * the default build pays nothing and the call sites need no `#ifdef`. Build a measuring binary
 * with
 *
 *     pixi run build_instrumented
 *
 * and run it as usual; the report goes to stderr at exit.
 *
 * ### What the instrumented build costs
 *
 * A `thread_local` pointer load and a handful of stores per node. That is not the hot path, but it
 * is not free either, and the instrumented build is for measuring the *bound*, never for timing
 * the walk.
 */

#include <cmath>
#include <cstdint>
#include <span>

#ifdef IGOR_BOUND_INSTRUMENTATION
#    include <algorithm>
#    include <cstdio>
#    include <cstring>
#    include <mutex>
#    include <string>
#    include <vector>
#endif

namespace BoundTightness {

/// Decades of over-estimation the histogram spans, and how finely each is cut. Quarter-decades
/// rather than whole ones because the answer turned out to live inside the first decade: a
/// histogram of whole decades reports "everything in bucket 0" and says nothing.
inline constexpr int kDecades = 20;
inline constexpr int kPerDecade = 4;
inline constexpr int kBuckets = kDecades * kPerDecade;

/// How far below 1 a ratio may fall before it counts as the bound being violated rather than
/// rounded. The bound and the realized probability are products of the same factors in different
/// associations, so an exact bound lands a few ULPs either side of equality; without this every
/// such leaf would be reported as unsound. One part in 10^9 is far above double's rounding and
/// far below any real violation.
inline constexpr long double kRoundingSlack = 1e-9L;

/// Bucket index for one leaf, as a pure function so it can be tested without the accumulator.
///
///  -1        the bound was below what the scenario realized, by more than rounding explains
///   0        the tightest quarter-decade: within a factor of 10^0.25
///   k        between 10^(k/kPerDecade) and 10^((k+1)/kPerDecade)
///   kBuckets everything above, including a realized probability of zero
///
/// A realized probability of zero has no finite ratio; it is counted in the last bucket rather
/// than dropped, because a scenario worth nothing that the bound let through is exactly the case
/// a tighter bound would have pruned.
inline int bucket_for(double bound, long double realized)
{
    if (realized <= 0.0L) {
        return kBuckets;
    }
    const long double ratio = static_cast<long double>(bound) / realized;
    if (ratio < 1.0L - kRoundingSlack) {
        return -1;
    }
    if (ratio <= 1.0L) {
        return 0;
    }
    const double steps = std::log10(static_cast<double>(ratio)) * kPerDecade;
    if (steps >= static_cast<double>(kBuckets)) {
        return kBuckets;
    }
    return static_cast<int>(steps);
}

/// Why a node's subtree produced no scenario. Pure, so the classification is testable without
/// running a walk -- the accumulator only counts what this returns.
enum class BarrenCause {
    /// A leaf was reached below this node. Not barren, whatever it was worth.
    NotBarren,
    /// Every realization the child event offered was rejected **before** any probability test:
    /// geometry, safety, read bounds, an empty realization set. A tighter probability bound
    /// cannot delete this node; only a feasibility test at the parent can.
    Starved,
    /// Every feasible child was probability-tested and fell below the cutoff. A tighter bound at
    /// *this* node would have deleted the node itself, and with it the whole enumeration below.
    Pruned,
    /// Some child was descended into and the waste is deeper down, where that child records its
    /// own cause. Counted so the three columns sum to the barren total, not as a finding.
    Hollow
};

/// \param reached_leaf whether any leaf at all was reached below the node
/// \param expanded     children the walk descended into
/// \param pruned       probability tests the node's child event failed
///
/// `Hollow` wins over `Pruned` when both apply: a node that expanded at least one child was not a
/// dead end at its own level, so its waste belongs to the child that was.
inline BarrenCause barren_cause(bool reached_leaf, unsigned expanded, unsigned pruned)
{
    if (reached_leaf) {
        return BarrenCause::NotBarren;
    }
    if (expanded != 0) {
        return BarrenCause::Hollow;
    }
    if (pruned != 0) {
        return BarrenCause::Pruned;
    }
    return BarrenCause::Starved;
}

#ifdef IGOR_BOUND_INSTRUMENTATION

namespace detail {

/// Deepest scenario path the per-depth node histogram distinguishes; anything deeper is folded
/// into the last row. The models IGoR ships have around ten events.
inline constexpr int kMaxDepth = 16;

/// One frame of the descent.
struct Frame {
    /// The bound that let the walk descend into this node.
    double bound = 0.0;
    /// The best probability any leaf below it turned out to reach.
    long double best_below = 0.0L;
    /// Children the walk descended into, and probability tests their event failed. Together with
    /// `reached_leaf` these are what `barren_cause()` reads.
    unsigned expanded = 0;
    unsigned pruned = 0;
    bool reached_leaf = false;
    /// The bound carried by whichever child produced `best_below` -- a child node's own bound, or
    /// a leaf's. Dividing this node's bound by it gives the slack that child's event resolved,
    /// and the product of those quotients down the winning path is this node's whole
    /// over-estimate. See the per-event table in report().
    double best_child_bound = 0.0;
};

/// One thread's tally. Never freed: the report reads it after the thread that owns it is gone.
struct Tally {
    std::uint64_t buckets[kBuckets + 1] = {};
    std::uint64_t unsound = 0;      ///< bucket -1: realized > bound
    double worst_unsound_ratio = 1.0; ///< the smallest bound/realized seen, when below 1
    std::uint64_t leaves = 0;

    /// The same histogram for internal nodes, split by depth -- and **signed**: a node whose
    /// bound sits below the best leaf under it goes into `node_under_buckets` at the decade of
    /// the shortfall, rather than out of the distribution. At the insertion depths every node is
    /// on that side (section 7.19), and a median taken over the handful that are not would
    /// describe nothing.
    std::uint64_t node_buckets[kMaxDepth][kBuckets + 1] = {};
    std::uint64_t node_under_buckets[kMaxDepth][kBuckets + 1] = {};
    std::uint64_t node_unsound[kMaxDepth] = {};
    std::uint64_t node_barren[kMaxDepth] = {};
    std::uint64_t nodes[kMaxDepth] = {};
    double node_worst_ratio[kMaxDepth] = {};

    /// The barren total, split by `BarrenCause`.
    std::uint64_t barren_starved[kMaxDepth] = {};
    std::uint64_t barren_pruned[kMaxDepth] = {};
    std::uint64_t barren_hollow[kMaxDepth] = {};

    /// The aggregate over-estimate, decomposed by the event that resolves it: the bound at a node
    /// of depth d-1 over the bound of the child of it that led to the best leaf, which is the
    /// slack the event at depth d was carrying. Signed like the node histogram, because an event
    /// whose bound *rises* on its own realization shows up here on the under side (section 7.19).
    std::uint64_t step_buckets[kMaxDepth][kBuckets + 1] = {};
    std::uint64_t step_under_buckets[kMaxDepth][kBuckets + 1] = {};
    std::uint64_t step_scored[kMaxDepth] = {};

    /// Prune decisions taken while no node was open, i.e. by the first event of the model. They
    /// belong to no frame; counted so they are not silently lost.
    std::uint64_t root_pruned = 0;

    /// The descent's open frames. Depth-first and per-thread, so a plain stack is enough.
    std::vector<Frame> stack;
};

inline std::mutex &registry_mutex()
{
    static std::mutex *m = new std::mutex();
    return *m;
}

inline std::vector<Tally *> &registry()
{
    static std::vector<Tally *> *r = new std::vector<Tally *>();
    return *r;
}

inline Tally &tally()
{
    thread_local Tally *mine = [] {
        Tally *fresh = new Tally();
        for (int d = 0; d != kMaxDepth; ++d) {
            fresh->node_worst_ratio[d] = 1.0;
        }
        std::lock_guard<std::mutex> lock(registry_mutex());
        registry().push_back(fresh);
        return fresh;
    }();
    return *mine;
}

/// Which event sits at each depth, for labelling the rows. The model ordering is the same for
/// every sequence and every thread, so the first name seen at a depth is the name.
///
/// **Copied**, never pointed at: the per-thread `Rec_Event` copies are destroyed before the
/// report runs, so a `const char *` kept from one of them dangles. (It did, and printed garbage.)
inline char (&depth_names())[kMaxDepth][40]
{
    static char names[kMaxDepth][40] = {};
    return names;
}

/// Name the row for `depth`, once. Races are benign: every writer writes the same bytes, and a
/// torn read can only mislabel a row of a report, never corrupt a count.
inline void remember_depth_name(int depth, const char *name)
{
    if (name == nullptr or depth < 0 or depth >= kMaxDepth or depth_names()[depth][0] != '\0') {
        return;
    }
    std::snprintf(depth_names()[depth], sizeof(depth_names()[depth]), "%s", name);
}

/// File one ratio into the signed pair of histograms: over 1 into `over` at its decade, under 1
/// into `under` at the decade of the shortfall -- which is the same bucketing applied to the
/// reciprocal, so one function decides both sides.
inline void tally_ratio(std::uint64_t *under, std::uint64_t *over, double numerator,
                        long double denominator)
{
    const int bucket = bucket_for(numerator, denominator);
    if (bucket < 0) {
        ++under[bucket_for(static_cast<double>(denominator),
                           static_cast<long double>(numerator))];
    } else {
        ++over[bucket];
    }
}

/// Decades at quantile `q` of a signed pair of histograms -- negative below 1 -- or NaN when the
/// pair is empty. Walking `under` from its far end inwards and then `over` outwards puts the
/// buckets in increasing order of ratio, which is what a quantile needs.
inline double signed_quantile(const std::uint64_t *under, const std::uint64_t *over, double q)
{
    std::uint64_t total = 0;
    for (int i = 0; i <= kBuckets; ++i) {
        total += under[i] + over[i];
    }
    if (total == 0) {
        return std::nan("");
    }
    const auto target = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(q * static_cast<double>(total)));
    std::uint64_t seen = 0;
    for (int i = kBuckets; i >= 0; --i) {
        seen += under[i];
        if (seen >= target) {
            return -static_cast<double>(i) / kPerDecade;
        }
    }
    for (int i = 0; i <= kBuckets; ++i) {
        seen += over[i];
        if (seen >= target) {
            return static_cast<double>(i) / kPerDecade;
        }
    }
    return static_cast<double>(kDecades);
}

/// A quantile, or a dash where there was nothing to take one of.
inline void print_decades(double decades)
{
    if (std::isnan(decades)) {
        std::fprintf(stderr, " %12s", "-");
    } else {
        std::fprintf(stderr, "      1e%-6.2f", decades);
    }
}

inline void report()
{
    Tally total;
    total.worst_unsound_ratio = 1.0;
    for (int d = 0; d != kMaxDepth; ++d) {
        total.node_worst_ratio[d] = 1.0;
    }
    {
        std::lock_guard<std::mutex> lock(registry_mutex());
        for (const Tally *t : registry()) {
            for (int i = 0; i <= kBuckets; ++i) {
                total.buckets[i] += t->buckets[i];
            }
            total.unsound += t->unsound;
            total.leaves += t->leaves;
            total.root_pruned += t->root_pruned;
            total.worst_unsound_ratio = std::min(total.worst_unsound_ratio, t->worst_unsound_ratio);
            for (int d = 0; d != kMaxDepth; ++d) {
                for (int i = 0; i <= kBuckets; ++i) {
                    total.node_buckets[d][i] += t->node_buckets[d][i];
                    total.node_under_buckets[d][i] += t->node_under_buckets[d][i];
                }
                total.node_unsound[d] += t->node_unsound[d];
                total.node_barren[d] += t->node_barren[d];
                total.barren_starved[d] += t->barren_starved[d];
                total.barren_pruned[d] += t->barren_pruned[d];
                total.barren_hollow[d] += t->barren_hollow[d];
                total.nodes[d] += t->nodes[d];
                for (int i = 0; i <= kBuckets; ++i) {
                    total.step_buckets[d][i] += t->step_buckets[d][i];
                    total.step_under_buckets[d][i] += t->step_under_buckets[d][i];
                }
                total.step_scored[d] += t->step_scored[d];
                if (t->node_worst_ratio[d] < total.node_worst_ratio[d]
                    or total.node_worst_ratio[d] == 0.0) {
                    total.node_worst_ratio[d] = t->node_worst_ratio[d];
                }
            }
        }
    }
    if (total.leaves == 0) {
        return;
    }
    std::fprintf(stderr, "\n=== bound / realized at scenario leaves ===\n");
    std::fprintf(stderr, "leaves            %llu\n",
                 static_cast<unsigned long long>(total.leaves));
    std::fprintf(stderr, "UNSOUND (< 1)     %llu", static_cast<unsigned long long>(total.unsound));
    if (total.unsound != 0) {
        std::fprintf(stderr, "   worst ratio %.6g", total.worst_unsound_ratio);
    }
    std::fprintf(stderr, "\n");

    std::uint64_t cumulative = 0;
    for (int i = 0; i <= kBuckets; ++i) {
        if (total.buckets[i] == 0) {
            continue;
        }
        cumulative += total.buckets[i];
        const double share = 100.0 * static_cast<double>(total.buckets[i])
                             / static_cast<double>(total.leaves);
        const double cum = 100.0 * static_cast<double>(cumulative)
                           / static_cast<double>(total.leaves);
        if (i == kBuckets) {
            std::fprintf(stderr, ">= 1e%-5.2f        %12llu  %6.2f%%  (cum %6.2f%%)\n",
                         static_cast<double>(kDecades),
                         static_cast<unsigned long long>(total.buckets[i]), share, cum);
        } else {
            std::fprintf(stderr, "[1e%-5.2f, 1e%-5.2f) %12llu  %6.2f%%  (cum %6.2f%%)\n",
                         static_cast<double>(i) / kPerDecade,
                         static_cast<double>(i + 1) / kPerDecade,
                         static_cast<unsigned long long>(total.buckets[i]), share, cum);
        }
    }

    //The distribution that matters: how far each node's bound sat above the best its subtree
    //turned out to hold. Summarised by decade rather than listed, because it is a table.
    std::fprintf(stderr, "\n=== bound / best realized below it, by depth ===\n");
    std::fprintf(stderr, "`barren` is a node no descendant of which ever reached a leaf -- the\n"
                         "bound admitted a subtree that produced nothing; the next table says\n"
                         "why. The quantiles cover the rest and are decades of over-estimation,\n"
                         "NEGATIVE where the bound sat below what the subtree reached -- which is\n"
                         "what `unsound` counts, and at the insertion depths is every node.\n\n");
    std::fprintf(stderr,
                 "depth        nodes       barren      unsound       median          p90          p99\n");
    for (int d = 0; d != kMaxDepth; ++d) {
        if (total.nodes[d] == 0) {
            continue;
        }
        std::fprintf(stderr, "%5d %12llu %12llu %12llu", d,
                     static_cast<unsigned long long>(total.nodes[d]),
                     static_cast<unsigned long long>(total.node_barren[d]),
                     static_cast<unsigned long long>(total.node_unsound[d]));
        for (const double q : {0.5, 0.9, 0.99}) {
            print_decades(signed_quantile(total.node_under_buckets[d], total.node_buckets[d], q));
        }
        std::fprintf(stderr, "\n");
    }
    for (int d = 0; d != kMaxDepth; ++d) {
        if (total.node_unsound[d] != 0) {
            std::fprintf(stderr,
                         "depth %d: %llu nodes whose bound sat BELOW the best leaf under them, "
                         "worst ratio %.6g\n",
                         d, static_cast<unsigned long long>(total.node_unsound[d]),
                         total.node_worst_ratio[d]);
        }
    }

    //Barren is most of the walk, and on its own it does not say what to fix. This does.
    std::fprintf(stderr, "\n=== why a barren node produced nothing, by depth ===\n");
    std::fprintf(stderr,
                 "starved  no child was ever probability-tested: every realization the child\n"
                 "         event offered was rejected on geometry, safety or read bounds. No\n"
                 "         probability bound can delete these; only a feasibility test can.\n"
                 "pruned   every feasible child was tested and fell below the cutoff. A tighter\n"
                 "         bound at this node would have deleted the node itself.\n"
                 "hollow   some child was expanded; the waste is deeper, and that child's own\n"
                 "         row records why.\n\n");
    std::fprintf(stderr,
                 "depth       barren      starved       pruned       hollow   starved%%  pruned%%\n");
    std::uint64_t all_barren = 0;
    std::uint64_t all_starved = 0;
    std::uint64_t all_pruned = 0;
    for (int d = 0; d != kMaxDepth; ++d) {
        if (total.nodes[d] == 0) {
            continue;
        }
        const double denom = total.node_barren[d] != 0
                                     ? static_cast<double>(total.node_barren[d])
                                     : 1.0;
        std::fprintf(stderr, "%5d %12llu %12llu %12llu %12llu   %7.2f%% %7.2f%%\n", d,
                     static_cast<unsigned long long>(total.node_barren[d]),
                     static_cast<unsigned long long>(total.barren_starved[d]),
                     static_cast<unsigned long long>(total.barren_pruned[d]),
                     static_cast<unsigned long long>(total.barren_hollow[d]),
                     100.0 * static_cast<double>(total.barren_starved[d]) / denom,
                     100.0 * static_cast<double>(total.barren_pruned[d]) / denom);
        all_barren += total.node_barren[d];
        all_starved += total.barren_starved[d];
        all_pruned += total.barren_pruned[d];
    }
    if (all_barren != 0) {
        std::fprintf(stderr,
                     "\nof %llu barren nodes, %.2f%% died with no probability test at all\n"
                     "and %.2f%% because every feasible child was below the cutoff.\n",
                     static_cast<unsigned long long>(all_barren),
                     100.0 * static_cast<double>(all_starved) / static_cast<double>(all_barren),
                     100.0 * static_cast<double>(all_pruned) / static_cast<double>(all_barren));
    }
    if (total.root_pruned != 0) {
        std::fprintf(stderr, "(%llu prune decisions were taken by the first event, outside any node)\n",
                     static_cast<unsigned long long>(total.root_pruned));
    }

    //Which event's bound carries the slack. This telescopes: multiply a depth's median down to
    //the leaf and you get the over-estimate the first table reports at that depth.
    std::fprintf(stderr, "\n=== the over-estimate, decomposed by the event that resolves it ===\n");
    std::fprintf(stderr,
                 "For each node on a path to the best leaf below it, the bound at its parent over\n"
                 "its own bound: how many decades the parent's bound was optimistic about THIS\n"
                 "event. 1e0.00 means the parent already knew what the event would contribute and\n"
                 "conditioning its bound would buy nothing. The steps multiply along the path, so\n"
                 "a row's decades are additive down to the leaf. NEGATIVE means the bound GREW on\n"
                 "the event's own realization, which is not something an upper bound may do.\n\n");
    std::fprintf(stderr,
                 "depth  event                              nodes       median          p90\n");
    for (int d = 0; d != kMaxDepth; ++d) {
        if (total.step_scored[d] == 0) {
            continue;
        }
        std::fprintf(stderr, "%5d  %-30.30s %11llu", d,
                     depth_names()[d][0] != '\0' ? depth_names()[d] : "(unnamed)",
                     static_cast<unsigned long long>(total.step_scored[d]));
        for (const double q : {0.5, 0.9}) {
            print_decades(signed_quantile(total.step_under_buckets[d], total.step_buckets[d], q));
        }
        std::fprintf(stderr, "\n");
    }

}

/// Prints at exit. The tallies are leaked on purpose, so they outlive every thread that wrote one.
struct Reporter {
    ~Reporter() { report(); }
};
inline Reporter reporter;

} // namespace detail

/// Record one probability test and its outcome, against the node whose child event took it.
inline void note_prune(bool pruned)
{
    if (not pruned) {
        return;
    }
    detail::Tally &mine = detail::tally();
    if (mine.stack.empty()) {
        ++mine.root_pruned;
        return;
    }
    ++mine.stack.back().pruned;
}

/// Record one leaf. `bound` is the upper bound that let this scenario through; `realized` is the
/// error-weighted probability it turned out to be worth.
inline void record(double bound, long double realized, const char *event_name)
{
    detail::Tally &mine = detail::tally();
    ++mine.leaves;
    detail::remember_depth_name(static_cast<int>(mine.stack.size()), event_name);
    if (not mine.stack.empty()) {
        detail::Frame &parent = mine.stack.back();
        //Reaching a leaf is what makes a node non-barren, even a leaf worth nothing: the walk
        //got all the way down, so the enumeration was not wasted in the sense this measures.
        parent.reached_leaf = true;
        if (realized > parent.best_below) {
            parent.best_below = realized;
            //A leaf's bound is its realized value, so this is the last step of the decomposition
            //and contributes nothing to the product -- but the event that took it is still named.
            parent.best_child_bound = bound;
        }
    }
    const int bucket = bucket_for(bound, realized);
    if (bucket < 0) {
        ++mine.unsound;
        const double ratio = static_cast<double>(static_cast<long double>(bound) / realized);
        if (ratio < mine.worst_unsound_ratio) {
            mine.worst_unsound_ratio = ratio;
        }
        return;
    }
    ++mine.buckets[bucket];
}

/// Open a frame for one internal node, carrying the bound that let the walk descend into it and
/// the slots that bound was the product of.
inline void enter(double bound, const char *event_name)
{
    detail::Tally &mine = detail::tally();
    detail::remember_depth_name(static_cast<int>(mine.stack.size()), event_name);
    mine.stack.emplace_back();
    mine.stack.back().bound = bound;
}

/// Close the frame `enter()` opened: record how far its bound sat above the best its subtree
/// reached, why it reached nothing if it reached nothing, which slot carried the slack, and carry
/// the best up to the parent.
inline void leave()
{
    detail::Tally &mine = detail::tally();
    if (mine.stack.empty()) {
        return;
    }
    const detail::Frame frame = mine.stack.back();
    mine.stack.pop_back();
    const int depth = std::min(static_cast<int>(mine.stack.size()), detail::kMaxDepth - 1);
    ++mine.nodes[depth];

    const BarrenCause cause = barren_cause(frame.reached_leaf, frame.expanded, frame.pruned);
    if (cause != BarrenCause::NotBarren) {
        //Nothing below this node ever reached a leaf. Counted apart from the ratio histogram,
        //because "the bound was enormously loose" and "there was nothing to be loose about" are
        //different observations and the second dominates.
        ++mine.node_barren[depth];
        switch (cause) {
        case BarrenCause::Starved:
            ++mine.barren_starved[depth];
            break;
        case BarrenCause::Pruned:
            ++mine.barren_pruned[depth];
            break;
        case BarrenCause::Hollow:
            ++mine.barren_hollow[depth];
            break;
        case BarrenCause::NotBarren:
            break;
        }
    } else {
        detail::tally_ratio(mine.node_under_buckets[depth], mine.node_buckets[depth], frame.bound,
                            frame.best_below);
        if (bucket_for(frame.bound, frame.best_below) < 0) {
            ++mine.node_unsound[depth];
            const double ratio = static_cast<double>(static_cast<long double>(frame.bound)
                                                     / frame.best_below);
            if (ratio < mine.node_worst_ratio[depth]) {
                mine.node_worst_ratio[depth] = ratio;
            }
        }

        //The same over-estimate, attributed to the one event that resolved part of it: how far
        //this node's bound sat above the bound of the child that led to the best leaf. Filed at
        //the child's depth, since it is that child's event whose bound was the optimistic one.
        if (frame.best_child_bound > 0.0 and depth + 1 < detail::kMaxDepth) {
            ++mine.step_scored[depth + 1];
            detail::tally_ratio(mine.step_under_buckets[depth + 1], mine.step_buckets[depth + 1],
                                frame.bound, static_cast<long double>(frame.best_child_bound));
        }
    }

    if (not mine.stack.empty()) {
        detail::Frame &parent = mine.stack.back();
        ++parent.expanded;
        if (frame.reached_leaf) {
            parent.reached_leaf = true;
        }
        if (frame.best_below > parent.best_below) {
            parent.best_below = frame.best_below;
            parent.best_child_bound = frame.bound;
        }
    }
}

#else

/// Compiled out. The call sites are unconditional so that enabling the measurement is a build
/// flag and never an edit.
inline void record(double, long double, const char *) {}
inline void enter(double, const char *) {}
inline void leave() {}
inline void note_prune(bool) {}

#endif

} // namespace BoundTightness
