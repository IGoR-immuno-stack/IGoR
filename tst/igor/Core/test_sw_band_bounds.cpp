/*
 * test_sw_band_bounds.cpp
 *
 *  Unit tests for swalign::compute_band_bounds, the closed-form banded-alignment "cone"
 *  computation -- see its doc comment in AlignerInternal.h for the derivation.
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

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <igor/Core/AlignerInternal.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

using swalign::compute_band_bounds;
using swalign::SwBandBounds;

namespace {

// Independent brute-force reference for the PREFIX (seed-restricted) cone: for every column j and
// row i, enumerates every seed diagonal d0 in [d_lo,d_hi] directly (rather than using
// compute_band_bounds's peak-finding + binary search) and combines the forward (seed-to-(i,j))
// and backward ((i,j)-to-far-corner) terms exactly as compute_band_bounds's Combined(i,j) does,
// so it validates that shortcut -- and the peak-finding -- against a naive O(range) search over
// the exact same criterion. d_lo/d_hi are clamped to the matrix's own feasible diagonal span
// first, matching how compute_band_bounds's clamp(i-j,...) degenerates to an unconstrained
// "i-j itself" when given INT32_MIN/MAX bounds.
SwBandBounds brute_force_prefix(int n_rows, int n_cols, long long d_lo, long long d_hi, double m, int g,
                                double score_threshold)
{
    const long long lo = std::max(d_lo, static_cast<long long>(-(n_cols - 1)));
    const long long hi = std::min(d_hi, static_cast<long long>(n_rows - 1));

    SwBandBounds result;
    result.lo.assign(n_cols, 1);
    result.hi.assign(n_cols, 0);
    if (lo > hi) {
        return result;
    }

    for (int j = 1; j != n_cols; ++j) {
        for (int i = 1; i != n_rows; ++i) {
            double best = -std::numeric_limits<double>::infinity();
            for (long long d0 = lo; d0 <= hi; ++d0) {
                const long long a = i - std::max<long long>(0, d0);
                const long long b = j - std::max<long long>(0, -d0);
                if (a < 0 || b < 0) {
                    continue; // d0's own boundary seed point doesn't fit in the matrix at all
                }
                const double gap_moves = static_cast<double>(std::llabs(static_cast<long long>(i - j) - d0));
                const double forward =
                        m * static_cast<double>(std::min(a, b)) - static_cast<double>(g) * gap_moves;
                best = std::max(best, forward);
            }
            const double backward = m * std::min(n_rows - 1 - i, n_cols - 1 - j);
            if (std::isfinite(best) && best + backward >= score_threshold) {
                if (result.lo[j] > result.hi[j]) {
                    result.lo[j] = i;
                }
                result.hi[j] = i;
            }
        }
    }
    return result;
}

// Mirror-image brute force for the SUFFIX (target-diagonal-restricted) cone: enumerates every
// target diagonal d1 in [d_lo,d_hi] directly, validating compute_band_bounds's suffix peak-finding
// shortcut the same way brute_force_prefix validates the prefix one.
SwBandBounds brute_force_suffix(int n_rows, int n_cols, long long d_lo, long long d_hi, double m, int g,
                                double score_threshold)
{
    const long long lo = std::max(d_lo, static_cast<long long>(-(n_cols - 1)));
    const long long hi = std::min(d_hi, static_cast<long long>(n_rows - 1));

    SwBandBounds result;
    result.lo.assign(n_cols, 1);
    result.hi.assign(n_cols, 0);
    if (lo > hi) {
        return result;
    }

    for (int j = 1; j != n_cols; ++j) {
        for (int i = 1; i != n_rows; ++i) {
            double best = -std::numeric_limits<double>::infinity();
            const long long k = i - j;
            for (long long d1 = lo; d1 <= hi; ++d1) {
                const long long delta = d1 - k;
                const long long a = (n_rows - 1 - i) - std::max<long long>(0, delta);
                const long long b = (n_cols - 1 - j) - std::max<long long>(0, -delta);
                if (a < 0 || b < 0) {
                    continue; // d1's own far-corner target point doesn't fit in the matrix at all
                }
                const double forward = m * static_cast<double>(std::min(a, b))
                                      - static_cast<double>(g) * static_cast<double>(std::llabs(delta));
                best = std::max(best, forward);
            }
            const double backward = m * std::min(i, j);
            if (std::isfinite(best) && best + backward >= score_threshold) {
                if (result.lo[j] > result.hi[j]) {
                    result.lo[j] = i;
                }
                result.hi[j] = i;
            }
        }
    }
    return result;
}

std::string describe_bounds(const SwBandBounds &bounds)
{
    std::ostringstream out;
    for (size_t j = 0; j < bounds.lo.size(); ++j) {
        out << "  col " << j << ": [" << bounds.lo[j] << "," << bounds.hi[j] << "]"
            << (bounds.lo[j] > bounds.hi[j] ? " (empty)" : "") << "\n";
    }
    return out.str();
}

void require_bounds_equal(const SwBandBounds &actual, const SwBandBounds &expected)
{
    INFO("actual:\n" << describe_bounds(actual) << "expected:\n" << describe_bounds(expected));
    REQUIRE(actual.lo.size() == expected.lo.size());
    REQUIRE(actual.hi.size() == expected.hi.size());
    for (size_t j = 0; j < expected.lo.size(); ++j) {
        // Normalize "empty" representations: compute_band_bounds and the brute-force reference
        // may pick different but equally-valid lo>hi sentinel pairs for an empty column.
        const bool actual_empty = actual.lo[j] > actual.hi[j];
        const bool expected_empty = expected.lo[j] > expected.hi[j];
        INFO("column " << j);
        REQUIRE(actual_empty == expected_empty);
        if (!expected_empty) {
            REQUIRE(actual.lo[j] == expected.lo[j]);
            REQUIRE(actual.hi[j] == expected.hi[j]);
        }
    }
}

} // namespace

TEST_CASE("compute_band_bounds is unconstrained when offsets span INT32_MIN/INT32_MAX",
          "[aligner][sw][band_bounds]")
{
    const int n_rows = 8;
    const int n_cols = 11;
    const double m = 2.0;
    const int g = 1;

    SECTION("flip_seqs = false")
    {
        const SwBandBounds bounds = compute_band_bounds(n_rows, n_cols, std::numeric_limits<int>::min(),
                                                        std::numeric_limits<int>::max(), m, g,
                                                        -std::numeric_limits<double>::infinity(), false, 7, 10);
        for (int j = 1; j != n_cols; ++j) {
            INFO("column " << j);
            REQUIRE(bounds.lo[j] == 1);
            REQUIRE(bounds.hi[j] == n_rows - 1);
        }
    }

    SECTION("flip_seqs = true")
    {
        const SwBandBounds bounds = compute_band_bounds(n_rows, n_cols, std::numeric_limits<int>::min(),
                                                        std::numeric_limits<int>::max(), m, g,
                                                        -std::numeric_limits<double>::infinity(), true, 7, 10);
        for (int j = 1; j != n_cols; ++j) {
            INFO("column " << j);
            REQUIRE(bounds.lo[j] == 1);
            REQUIRE(bounds.hi[j] == n_rows - 1);
        }
    }
}

TEST_CASE("compute_band_bounds matches brute-force enumeration across randomized small configs",
          "[aligner][sw][band_bounds]")
{
    std::mt19937_64 rng(123456789ULL);
    std::uniform_int_distribution<int> dim_dist(2, 14);
    std::uniform_int_distribution<int> offset_dist(-12, 12);
    std::uniform_real_distribution<double> match_dist(0.5, 5.0);
    std::uniform_int_distribution<int> gap_dist(1, 6);
    // Mix -infinity (the SwDPConfig default -- exercises the "no score-based exclusion" path) with
    // real, sometimes-restrictive thresholds spanning comfortably-reachable to unreachable values.
    std::uniform_real_distribution<double> threshold_dist(-40.0, 40.0);
    std::uniform_int_distribution<int> use_real_threshold_dist(0, 1);

    for (int trial = 0; trial != 200; ++trial) {
        const int n_rows = dim_dist(rng);
        const int n_cols = dim_dist(rng);
        const int a = offset_dist(rng);
        const int b = offset_dist(rng);
        const int min_offset = std::min(a, b);
        const int max_offset = std::max(a, b);
        const double m = match_dist(rng);
        const int g = gap_dist(rng);
        const int data_seq_size = n_rows - 1;
        const int genomic_seq_size = n_cols - 1;
        const double score_threshold = use_real_threshold_dist(rng) ? threshold_dist(rng)
                                                                    : -std::numeric_limits<double>::infinity();

        INFO("trial=" << trial << " n_rows=" << n_rows << " n_cols=" << n_cols << " min_offset=" << min_offset
                       << " max_offset=" << max_offset << " m=" << m << " g=" << g
                       << " score_threshold=" << score_threshold);

        {
            const SwBandBounds actual = compute_band_bounds(n_rows, n_cols, min_offset, max_offset, m, g,
                                                             score_threshold, false, data_seq_size, genomic_seq_size);
            const SwBandBounds expected =
                    brute_force_prefix(n_rows, n_cols, min_offset, max_offset, m, g, score_threshold);
            require_bounds_equal(actual, expected);
        }
        {
            const SwBandBounds actual = compute_band_bounds(n_rows, n_cols, min_offset, max_offset, m, g,
                                                             score_threshold, true, data_seq_size, genomic_seq_size);
            // flip_seqs=true's suffix side uses target diagonal range [C-max_offset, C-min_offset]
            // with C = data_seq_size - genomic_seq_size (see compute_band_bounds's doc comment).
            const long long C = static_cast<long long>(data_seq_size) - static_cast<long long>(genomic_seq_size);
            const SwBandBounds expected_suffix =
                    brute_force_suffix(n_rows, n_cols, C - max_offset, C - min_offset, m, g, score_threshold);
            require_bounds_equal(actual, expected_suffix);
        }
    }
}

TEST_CASE("without a real score_threshold, offset restriction alone gives a rectangle, not a cone",
          "[aligner][sw][band_bounds]")
{
    // -infinity is SwDPConfig's default: Combined(i,j) is always finite, so ">= -infinity" holds
    // everywhere the offset-restricted seed diagonal is geometrically defined at all -- there is
    // no score-based narrowing. But the seed diagonal's own boundary point (i0=max(0,offset),
    // j0=0 here, since offset>=0) is a genuine geometric floor: no row before i0 is reachable at
    // all, in ANY column. So the expected shape is a RECTANGLE -- a single constant row cutoff at
    // i0, same for every column, unbounded above -- not the gradually-widening-from-row-1 shape a
    // "narrowing cone" would produce. This documents the "rectangle, not cone" finding directly:
    // offset alone restricts *which* seed diagonal is admissible and where its floor sits, not how
    // far a candidate may drift from it.
    const int n_rows = 40;
    const int n_cols = 40;
    const double m = 2.0;
    const int g = 1;
    const int offset = 3;

    const SwBandBounds bounds = compute_band_bounds(n_rows, n_cols, offset, offset, m, g,
                                                     -std::numeric_limits<double>::infinity(), false, n_rows - 1,
                                                     n_cols - 1);

    for (int j = 1; j != n_cols; ++j) {
        INFO("column " << j);
        REQUIRE(bounds.lo[j] == offset); // constant floor i0 = max(0,offset) = offset, every column
        REQUIRE(bounds.hi[j] == n_rows - 1); // unbounded above: no threshold-based narrowing
    }
}

TEST_CASE("with a real score_threshold, offset restriction DOES produce a genuine widening-then-"
          "narrowing cone",
          "[aligner][sw][band_bounds]")
{
    // The forward+backward Combined bound along the exact admissible diagonal itself is constant
    // (m * (n_rows-1-offset), see compute_band_bounds's doc comment): every column's "on-diagonal"
    // cell scores identically in the best case, so a threshold comfortably below that constant
    // still empties out both ends of the matrix (too little accumulated near the seed, too little
    // remaining room near the far corner) while leaving a genuinely widened region in the middle --
    // the cone shape the offset/score_threshold restriction is meant to produce.
    const int n_rows = 40;
    const int n_cols = 40;
    const double m = 2.0;
    const int g = 1;
    const int offset = 3;
    const double threshold = 0.9 * m * (n_rows - 1 - offset);

    const SwBandBounds bounds =
            compute_band_bounds(n_rows, n_cols, offset, offset, m, g, threshold, false, n_rows - 1, n_cols - 1);

    auto width_at = [&](int j) { return bounds.lo[j] > bounds.hi[j] ? 0 : bounds.hi[j] - bounds.lo[j] + 1; };

    const int width_first = width_at(1);
    const int width_last = width_at(n_cols - 1);
    int widest = 0;
    for (int j = 1; j != n_cols; ++j) {
        widest = std::max(widest, width_at(j));
    }

    INFO("width at column 1 = " << width_first << ", width at last column = " << width_last
                                 << ", widest = " << widest);
    REQUIRE(widest > width_first);
    REQUIRE(widest > width_last);
}

TEST_CASE("compute_band_bounds's score_threshold axis trims diagonals that can never qualify",
          "[aligner][sw][band_bounds]")
{
    const int n_rows = 12;
    const int n_cols = 12;
    const double m = 2.0;
    const int g = 1;

    // Full unconstrained offset range, so only score_threshold shapes the result.
    const int min_offset = std::numeric_limits<int>::min();
    const int max_offset = std::numeric_limits<int>::max();

    SECTION("a very high threshold that no diagonal can reach empties the whole band")
    {
        const double threshold = m * (n_rows + n_cols); // unreachable: exceeds the full-matrix max
        const SwBandBounds bounds =
                compute_band_bounds(n_rows, n_cols, min_offset, max_offset, m, g, threshold, false, n_rows - 1,
                                    n_cols - 1);
        for (int j = 1; j != n_cols; ++j) {
            INFO("column " << j);
            REQUIRE(bounds.lo[j] > bounds.hi[j]);
        }
    }

    SECTION("a moderate threshold narrows the band relative to the unconstrained case")
    {
        const double threshold = m * (n_rows - 2); // reachable only near the main diagonal
        const SwBandBounds narrow =
                compute_band_bounds(n_rows, n_cols, min_offset, max_offset, m, g, threshold, false, n_rows - 1,
                                    n_cols - 1);
        const SwBandBounds wide =
                compute_band_bounds(n_rows, n_cols, min_offset, max_offset, m, g,
                                    -std::numeric_limits<double>::infinity(), false, n_rows - 1, n_cols - 1);

        bool any_narrower = false;
        for (int j = 1; j != n_cols; ++j) {
            if (narrow.lo[j] > narrow.hi[j]) {
                continue;
            }
            REQUIRE(wide.lo[j] <= wide.hi[j]);
            REQUIRE(narrow.lo[j] >= wide.lo[j]);
            REQUIRE(narrow.hi[j] <= wide.hi[j]);
            if (narrow.lo[j] > wide.lo[j] || narrow.hi[j] < wide.hi[j]) {
                any_narrower = true;
            }
        }
        REQUIRE(any_narrower);
    }
}

TEST_CASE("compute_band_bounds handles degenerate tiny matrices", "[aligner][sw][band_bounds]")
{
    SECTION("n_rows == 2, n_cols == 2")
    {
        const SwBandBounds bounds = compute_band_bounds(2, 2, std::numeric_limits<int>::min(),
                                                        std::numeric_limits<int>::max(), 2.0, 1,
                                                        -std::numeric_limits<double>::infinity(), false, 1, 1);
        REQUIRE(bounds.lo[1] == 1);
        REQUIRE(bounds.hi[1] == 1);
    }

    SECTION("min_offset > max_offset (already-empty input range) never crashes and yields an empty band")
    {
        const SwBandBounds bounds =
                compute_band_bounds(6, 6, 5, -5, 2.0, 1, -std::numeric_limits<double>::infinity(), false, 5, 5);
        for (int j = 1; j != 6; ++j) {
            REQUIRE(bounds.lo[j] > bounds.hi[j]);
        }
    }
}
