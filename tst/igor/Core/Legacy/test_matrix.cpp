/*
 * test_matrix.cpp
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

#include <igor/Core/Legacy/Utils.h>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

// ===========================================================================================
// Matrix bounds (plan section 7.18, repair R0).
//
// operator()'s assert checked only the upper bound, so a negative index -- which the indices
// being `int` makes expressible, and which a credited length derived by subtraction can
// produce -- satisfied it and read before the allocation. The predicate is tested rather than
// the assert because the default build is RelWithDebInfo, where NDEBUG compiles the assert
// out; this is the pattern first_unfilled_segment() established for the same reason (7.14).
// ===========================================================================================

TEST_CASE("Matrix::in_range: rejects a negative index, not only an oversized one", "[matrix][bounds]")
{
    Matrix<double> m(3, 4);

    SECTION("Every cell of the matrix is in range")
    {
        for (int i = 0; i != 3; ++i) {
            for (int j = 0; j != 4; ++j) {
                CHECK(m.in_range(i, j));
            }
        }
    }

    SECTION("One past each end is out of range")
    {
        CHECK_FALSE(m.in_range(3, 0));
        CHECK_FALSE(m.in_range(0, 4));
    }

    SECTION("A negative index is out of range on either axis")
    {
        // The case the upper-bound-only check accepted. j == -1 is the measured one: the linear
        // index i + rows * j then lands `rows` elements before the allocation.
        CHECK_FALSE(m.in_range(0, -1));
        CHECK_FALSE(m.in_range(-1, 0));
        CHECK_FALSE(m.in_range(-1, -1));
    }

    SECTION("A huge unsigned count converted back to int is caught")
    {
        // How 7.18's value actually arrives: a negative credited length wraps in a size_t
        // parameter, and Matrix's `int` indices convert it back to -1.
        const std::size_t wrapped = static_cast<std::size_t>(-1);
        CHECK_FALSE(m.in_range(0, static_cast<int>(wrapped)));
    }
}

TEST_CASE("Matrix::get_field: throws on a negative index", "[matrix][bounds]")
{
    Matrix<double> m(3, 4);

    // get_field is the checked accessor, and it had the same omission as the assert. Its one
    // caller is the stream operator, which iterates in range, so this path is unreachable
    // today -- which is why widening the condition is bitwise.
    CHECK_THROWS_AS(m.get_field(0, -1), std::length_error);
    CHECK_THROWS_AS(m.get_field(-1, 0), std::length_error);
    CHECK_THROWS_AS(m.get_field(3, 0), std::length_error);
    CHECK_THROWS_AS(m.get_field(0, 4), std::length_error);
    CHECK_NOTHROW(m.get_field(2, 3));
}
