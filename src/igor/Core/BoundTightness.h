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
 * \brief How far the pruning bound sits above the probability a scenario actually realizes.
 *
 * Every scenario reaches its leaf carrying an *upper bound* on what it could still be worth, and
 * the walk discards it when that bound falls below the threshold. How much work that saves, and
 * how much a **tighter** bound would save on top, is the question section 6.10 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md leaves open and R6 has to answer: conditioning the bound on
 * the scenario's already-chosen parents costs storage, and nobody has measured what it buys.
 *
 * Two distributions are reported, because the obvious one turns out not to answer the question.
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
 * Both also answer a question nobody had asked: whether the bound is a bound at all. A ratio
 * below one means the pruning is unsound -- scenarios are being discarded that should not be --
 * and the report counts those separately.
 *
 * ### It is compiled out unless asked for
 *
 * `record()` is an empty inline function unless `IGOR_BOUND_INSTRUMENTATION` is defined, so the
 * default build pays nothing and the call site needs no `#ifdef`. Build a measuring binary with
 *
 *     pixi run build_instrumented
 *
 * and run it as usual; the report goes to stderr at exit.
 *
 * ### What the instrumented build costs
 *
 * One `thread_local` pointer load, a `log10`, and an increment per leaf. Leaves are orders of
 * magnitude rarer than scenario nodes, so this is not the hot path -- but it is not free either,
 * and the instrumented build is for measuring the *bound*, never for timing the walk.
 */

#include <cmath>
#include <cstdint>

#ifdef IGOR_BOUND_INSTRUMENTATION
#    include <algorithm>
#    include <cstdio>
#    include <mutex>
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

#ifdef IGOR_BOUND_INSTRUMENTATION

namespace detail {

/// Deepest scenario path the per-depth node histogram distinguishes; anything deeper is folded
/// into the last row. The models IGoR ships have around ten events.
inline constexpr int kMaxDepth = 16;

/// One frame of the descent: the bound this node carried, and the best any leaf below it reached.
struct Frame {
    double bound = 0.0;
    long double best_below = 0.0L;
};

/// One thread's tally. Never freed: the report reads it after the thread that owns it is gone.
struct Tally {
    std::uint64_t buckets[kBuckets + 1] = {};
    std::uint64_t unsound = 0;      ///< bucket -1: realized > bound
    double worst_unsound_ratio = 1.0; ///< the smallest bound/realized seen, when below 1
    std::uint64_t leaves = 0;

    /// The same histogram for internal nodes, split by depth.
    std::uint64_t node_buckets[kMaxDepth][kBuckets + 1] = {};
    std::uint64_t node_unsound[kMaxDepth] = {};
    std::uint64_t node_barren[kMaxDepth] = {};
    std::uint64_t nodes[kMaxDepth] = {};
    double node_worst_ratio[kMaxDepth] = {};

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
            total.worst_unsound_ratio = std::min(total.worst_unsound_ratio, t->worst_unsound_ratio);
            for (int d = 0; d != kMaxDepth; ++d) {
                for (int i = 0; i <= kBuckets; ++i) {
                    total.node_buckets[d][i] += t->node_buckets[d][i];
                }
                total.node_unsound[d] += t->node_unsound[d];
                total.node_barren[d] += t->node_barren[d];
                total.nodes[d] += t->nodes[d];
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
                         "bound admitted a subtree that produced nothing. The quantiles cover\n"
                         "only the rest, and are decades of over-estimation.\n\n");
    std::fprintf(stderr,
                 "depth        nodes       barren     unsound   median    p90       p99\n");
    for (int d = 0; d != kMaxDepth; ++d) {
        if (total.nodes[d] == 0) {
            continue;
        }
        std::uint64_t scored = 0;
        for (int i = 0; i <= kBuckets; ++i) {
            scored += total.node_buckets[d][i];
        }
        const auto quantile = [&](double q) {
            if (scored == 0) {
                return -1.0;
            }
            const std::uint64_t target = static_cast<std::uint64_t>(q * static_cast<double>(scored));
            std::uint64_t seen = 0;
            for (int i = 0; i <= kBuckets; ++i) {
                seen += total.node_buckets[d][i];
                if (seen >= target) {
                    return static_cast<double>(i) / kPerDecade;
                }
            }
            return static_cast<double>(kDecades);
        };
        std::fprintf(stderr, "%5d %12llu %12llu %11llu   1e%-7.2f 1e%-7.2f 1e%-7.2f\n", d,
                     static_cast<unsigned long long>(total.nodes[d]),
                     static_cast<unsigned long long>(total.node_barren[d]),
                     static_cast<unsigned long long>(total.node_unsound[d]), quantile(0.5),
                     quantile(0.9), quantile(0.99));
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
}

/// Prints at exit. The tallies are leaked on purpose, so they outlive every thread that wrote one.
struct Reporter {
    ~Reporter() { report(); }
};
inline Reporter reporter;

} // namespace detail

/// Record one leaf. `bound` is the upper bound that let this scenario through; `realized` is the
/// error-weighted probability it turned out to be worth.
inline void record(double bound, long double realized)
{
    detail::Tally &mine = detail::tally();
    ++mine.leaves;
    if (not mine.stack.empty() and realized > mine.stack.back().best_below) {
        mine.stack.back().best_below = realized;
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

/// Open a frame for one internal node, carrying the bound that let the walk descend into it.
inline void enter(double bound)
{
    detail::Tally &mine = detail::tally();
    mine.stack.push_back(detail::Frame{bound, 0.0L});
}

/// Close the frame `enter()` opened, record how far its bound sat above the best its subtree
/// reached, and carry that best up to the parent.
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
    const int bucket = bucket_for(frame.bound, frame.best_below);
    if (frame.best_below <= 0.0L) {
        //Nothing below this node ever reached a leaf: every descendant was pruned. Counted
        //apart, because "the bound was enormously loose" and "there was nothing to be loose
        //about" are different observations and the second dominates.
        ++mine.node_barren[depth];
    } else if (bucket < 0) {
        ++mine.node_unsound[depth];
        const double ratio = static_cast<double>(static_cast<long double>(frame.bound)
                                                 / frame.best_below);
        if (ratio < mine.node_worst_ratio[depth]) {
            mine.node_worst_ratio[depth] = ratio;
        }
    } else {
        ++mine.node_buckets[depth][bucket];
    }
    if (not mine.stack.empty() and frame.best_below > mine.stack.back().best_below) {
        mine.stack.back().best_below = frame.best_below;
    }
}

#else

/// Compiled out. The call sites are unconditional so that enabling the measurement is a build
/// flag and never an edit.
inline void record(double, long double) {}
inline void enter(double) {}
inline void leave() {}

#endif

} // namespace BoundTightness
