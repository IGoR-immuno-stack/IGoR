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

TEST_CASE("BoundTightness::record is callable whether or not the build measures",
          "[bound_tightness]")
{
    // The call site in Rec_Event::iterate_wrap_up is unconditional; this is what says so. In a
    // default build it does nothing, in an instrumented one it tallies -- either way it links.
    BoundTightness::record(1.0, 0.5L);
    BoundTightness::record(0.0, 0.0L);
    SUCCEED();
}
