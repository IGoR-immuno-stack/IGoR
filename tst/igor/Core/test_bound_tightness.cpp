/*
 * test_bound_tightness.cpp
 *
 *  Unit tests for the bound / realized-probability instrumentation (plan step 5a).
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

/**
 * The accumulator itself only runs in a build configured with
 * `-DENABLE_BOUND_INSTRUMENTATION=ON`, so what is tested here is the part that decides where a
 * leaf lands -- a pure function, compiled in either way. That keeps the classification honest
 * without making the default build carry the machinery.
 */

#include <igor/Core/BoundTightness.h>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("BoundTightness::bucket_for classifies a leaf by decades of over-estimation",
          "[bound_tightness]")
{
    SECTION("A bound equal to what was realized is the tightest possible")
    {
        CHECK(BoundTightness::bucket_for(0.5, 0.5L) == 0);
    }

    SECTION("The scale is quarter-decades, so the first factor of ten is resolved")
    {
        // Whole decades would report every leaf of the TRB corpus in one bucket -- the reason
        // the resolution is what it is. See section 6.10 of the plan.
        CHECK(BoundTightness::bucket_for(1.7, 1.0L) == 0);  // 10^0.23
        CHECK(BoundTightness::bucket_for(1.8, 1.0L) == 1);  // 10^0.26
        CHECK(BoundTightness::bucket_for(10.0, 1.0L) == BoundTightness::kPerDecade);
        CHECK(BoundTightness::bucket_for(100.0, 1.0L) == 2 * BoundTightness::kPerDecade);
    }

    SECTION("The bucket is the ratio's, not either value's")
    {
        // Same ratio, probabilities twenty orders of magnitude apart.
        CHECK(BoundTightness::bucket_for(1e-3, 1e-6L) == 3 * BoundTightness::kPerDecade);
        CHECK(BoundTightness::bucket_for(1e-23, 1e-26L) == 3 * BoundTightness::kPerDecade);
    }

    SECTION("A ratio a few ULPs below one is rounding, not a violated bound")
    {
        // The bound and the realized probability are products of the same factors in different
        // associations, so an exact bound lands either side of equality. Without the slack every
        // such leaf reads as unsound -- which is what the first run of this instrument reported,
        // for 18% of the corpus.
        CHECK(BoundTightness::bucket_for(1.0 - 1e-15, 1.0L) == 0);
        CHECK(BoundTightness::bucket_for(1.0 - 1e-6, 1.0L) == -1);
    }

    SECTION("A bound below what was realized is not a bound, and is reported apart")
    {
        // The case that would mean the pruning discards scenarios it should keep. It gets its
        // own bucket rather than being folded into the tightest one, because "tight" and
        // "wrong" are opposite findings.
        CHECK(BoundTightness::bucket_for(0.5, 0.6L) == -1);
        CHECK(BoundTightness::bucket_for(0.0, 1e-30L) == -1);
    }

    SECTION("A scenario worth nothing lands in the last bucket, not out of the count")
    {
        // No finite ratio exists, but a scenario the bound admitted and that turned out to be
        // worth zero is exactly what a tighter bound would have pruned -- so it is counted at
        // the far end rather than dropped.
        CHECK(BoundTightness::bucket_for(1e-9, 0.0L) == BoundTightness::kBuckets);
        CHECK(BoundTightness::bucket_for(0.0, 0.0L) == BoundTightness::kBuckets);
    }

    SECTION("Everything past the last decade saturates there")
    {
        CHECK(BoundTightness::bucket_for(1.0, 1e-30L) == BoundTightness::kBuckets);
        CHECK(BoundTightness::bucket_for(1.0, 1e-19L) == 19 * BoundTightness::kPerDecade);
    }
}

TEST_CASE("BoundTightness::barren_cause separates the waste a better bound could delete",
          "[bound_tightness]")
{
    using BoundTightness::BarrenCause;
    using BoundTightness::barren_cause;

    SECTION("A node with a leaf below it is not barren, whatever that leaf was worth")
    {
        // Including a leaf worth nothing: the walk reached the bottom, so the enumeration was
        // not wasted in the sense this measures. The ratio histogram takes it from there.
        CHECK(barren_cause(true, 0, 0) == BarrenCause::NotBarren);
        CHECK(barren_cause(true, 3, 7) == BarrenCause::NotBarren);
    }

    SECTION("No probability test at all means no probability bound could have helped")
    {
        // Every realization the child event offered was rejected on geometry, safety or read
        // bounds before any bound was compared. This is the count that decides whether
        // tightening the bound is worth doing at all -- if it dominates, it is not.
        CHECK(barren_cause(false, 0, 0) == BarrenCause::Starved);
    }

    SECTION("All feasible children tested and rejected is the case a tighter bound deletes")
    {
        // The node itself would never have been expanded had its own bound been tight enough to
        // fall below the cutoff -- so this is the volume a better bound is competing for.
        CHECK(barren_cause(false, 0, 1) == BarrenCause::Pruned);
        CHECK(barren_cause(false, 0, 400) == BarrenCause::Pruned);
    }

    SECTION("Expanding any child moves the finding down to that child")
    {
        // A node that descended into something was not a dead end at its own level, so it is
        // neither starved nor pruned however many of its other children were rejected. The
        // child that was the dead end carries the finding, at its own depth.
        CHECK(barren_cause(false, 1, 0) == BarrenCause::Hollow);
        CHECK(barren_cause(false, 1, 99) == BarrenCause::Hollow);
    }
}

TEST_CASE("BoundTightness entry points are callable whether or not the build measures",
          "[bound_tightness]")
{
    // The call sites in Rec_Event::iterate_wrap_up and ExplorationContext are unconditional;
    // this is what says so. In a default build they do nothing, in an instrumented one they
    // tally -- either way they link.
    BoundTightness::record(1.0, 0.5L, "V_gene");
    BoundTightness::record(0.0, 0.0L, nullptr);
    BoundTightness::note_prune(true);
    BoundTightness::note_prune(false);
    SUCCEED();
}
