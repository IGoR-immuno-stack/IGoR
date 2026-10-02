/*
 * test_error_rate_bound.cpp
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

#include <igor/Model/Legacy/Singleerrorrate.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

// ===========================================================================================
// The error-rate pruning bound's counts (plan section 7.18, repair R5a).
//
// The accessor took both counts as size_t. A credited length derived by subtraction could come
// out negative, wrap to a huge value, defeat the growth check by wrapping back to 9, and reach
// Matrix's `int` indices as -1 -- a read before the allocation, returning whatever sat there.
// The counts are `int` since R5a and a negative one is refused at the entry point, before any
// of that can happen. R0's Matrix bound is the second line; this is the first.
// ===========================================================================================

TEST_CASE("Error_rate::get_err_rate_upper_bound: the formula, including past the cached matrix",
          "[error_rate][bounds]")
{
    Single_error_rate error_rate(0.1);
    const double r = 0.1;

    CHECK(error_rate.get_err_rate_upper_bound(0, 0) == 1.0);
    CHECK_THAT(error_rate.get_err_rate_upper_bound(1, 7),
               Catch::Matchers::WithinRel(std::pow(r / 3.0, 1) * std::pow(1.0 - r, 7), 1e-12));
    // Far enough out that the matrix has to grow, which is the path the wrapped count used to
    // skip.
    CHECK_THAT(error_rate.get_err_rate_upper_bound(3, 250),
               Catch::Matchers::WithinRel(std::pow(r / 3.0, 3) * std::pow(1.0 - r, 250), 1e-12));
}

TEST_CASE("Error_rate::get_err_rate_upper_bound: a negative count is refused, not indexed with",
          "[error_rate][bounds]")
{
    // The measured 7.18 value: an error-free count of -1, which returned 5.31441e-07 on one
    // fixture and 0 on another depending on what preceded the allocation.
    Single_error_rate error_rate(0.1);
    CHECK_THROWS_AS(error_rate.get_err_rate_upper_bound(3, -1), std::invalid_argument);
    CHECK_THROWS_AS(error_rate.get_err_rate_upper_bound(-1, 3), std::invalid_argument);
    // The refusal leaves the accessor usable.
    CHECK(error_rate.get_err_rate_upper_bound(0, 0) == 1.0);
}
