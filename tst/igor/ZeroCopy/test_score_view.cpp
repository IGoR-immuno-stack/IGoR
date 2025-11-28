/*
 * test_score_view.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ZeroCopy ScoreView
 */

#include <igor/ZeroCopy/ScoreView.h>

#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <vector>

using namespace igor;

TEST_CASE("ScoreView provides correct 2D view", "[zerocopy][score_view]") {
    // Create a flat buffer representing a 3x4 matrix (3 rows, 4 genes)
    // Row 0: 0.0, 1.0, 2.0, 3.0
    // Row 1: 4.0, 5.0, 6.0, 7.0
    // Row 2: 8.0, 9.0, 10.0, 11.0
    size_t num_rows = 3;
    size_t num_genes = 4;
    std::vector<double> data(num_rows * num_genes);
    std::iota(data.begin(), data.end(), 0.0);

    ScoreView view(data.data(), num_rows, num_genes);

    // Check individual elements
    REQUIRE(view.get_score(0, 0) == 0.0);
    REQUIRE(view.get_score(0, 3) == 3.0);
    REQUIRE(view.get_score(1, 0) == 4.0);
    REQUIRE(view.get_score(1, 2) == 6.0);
    REQUIRE(view.get_score(2, 3) == 11.0);

    // Check mdspan direct access
    auto matrix = view.get_matrix();
    REQUIRE(matrix.extent(0) == num_rows);
    REQUIRE(matrix.extent(1) == num_genes);
    REQUIRE(matrix(1, 1) == 5.0);
}

TEST_CASE("ScoreView handles single row", "[zerocopy][score_view]") {
    size_t num_rows = 1;
    size_t num_genes = 5;
    std::vector<double> data = {10.0, 20.0, 30.0, 40.0, 50.0};

    ScoreView view(data.data(), num_rows, num_genes);

    REQUIRE(view.get_score(0, 0) == 10.0);
    REQUIRE(view.get_score(0, 4) == 50.0);
}
