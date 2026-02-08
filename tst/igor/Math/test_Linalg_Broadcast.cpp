/*
 * test_Linalg_Broadcast.cpp
 *
 *  Created on: Feb 08, 2026
 *      Author: IGoR Agent
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>
// #include <cmath>
// #include <limits>

using namespace igor::math;
using Catch::Matchers::WithinAbs;

TEST_CASE("broadcast_to: 1D vector to 2D matrix", "[linalg][broadcast]") {
    // Create 1D vector [10, 20, 30]
    Tensor<double> vec({3});
    vec(0) = 10.0;
    vec(1) = 20.0;
    vec(2) = 30.0;

    // Broadcast to 4x3 matrix
    auto vec_view = vec.view<1>();
    auto broadcasted = linalg::broadcast_to(vec_view, std::array<size_t, 2>{4, 3});

    // Check shape
    REQUIRE(broadcasted.extent(0) == 4);
    REQUIRE(broadcasted.extent(1) == 3);

    // Check strides (key property!)
    REQUIRE(broadcasted.stride(0) == 0);  // No advance for row change
    REQUIRE(broadcasted.stride(1) == 1);  // Normal advance for column change

    // Verify data replication
    for (size_t row = 0; row < 4; ++row) {
        REQUIRE(broadcasted[row, 0] == 10.0);
        REQUIRE(broadcasted[row, 1] == 20.0);
        REQUIRE(broadcasted[row, 2] == 30.0);
    }

    // Verify zero-copy: modifying source affects view
    vec(1) = 99.0;
    for (size_t row = 0; row < 4; ++row) {
        REQUIRE(broadcasted[row, 1] == 99.0);
    }
}

TEST_CASE("broadcast_multiply: P(V,D,J) * P(J)", "[linalg][broadcast]") {
    // Setup: 3D probability tensor P(V,D,J)
    Tensor<double> p_vdj({4, 3, 5});  // Simplified from (64, 25, 12)
    for (size_t i = 0; i < p_vdj.size(); ++i) {
        p_vdj.data()[i] = 1.0;  // Uniform distribution
    }

    // Setup: 1D marginal P(J)
    Tensor<double> p_j({5});
    for (size_t j = 0; j < 5; ++j) {
        p_j(j) = 0.1 + j * 0.1;  // [0.1, 0.2, 0.3, 0.4, 0.5]
    }

    // Result tensor
    Tensor<double> result({4, 3, 5});

    // Perform broadcast multiply
    linalg::broadcast_multiply(p_vdj.view<3>(), p_j.view<1>(), result.view<3>());

    // Verify: Each (v,d) slice multiplied by p_j
    for (size_t v = 0; v < 4; ++v) {
        for (size_t d = 0; d < 3; ++d) {
            for (size_t j = 0; j < 5; ++j) {
                double expected = 1.0 * p_j(j);
                REQUIRE_THAT(result(v, d, j), WithinAbs(expected, 1e-10));
            }
        }
    }
}

TEST_CASE("broadcast workflow: normalize per-column", "[linalg][broadcast]") {
    // Matrix 3x4 where each column should be normalized independently
    Tensor<double> matrix({3, 4});
    matrix(0,0) = 1.0; matrix(0,1) = 2.0; matrix(0,2) = 1.0; matrix(0,3) = 4.0;
    matrix(1,0) = 2.0; matrix(1,1) = 2.0; matrix(1,2) = 2.0; matrix(1,3) = 2.0;
    matrix(2,0) = 3.0; matrix(2,1) = 2.0; matrix(2,2) = 3.0; matrix(2,3) = 2.0;

    // Step 1: Compute column sums
    Tensor<double> col_sums({4});
    for (size_t col = 0; col < 4; ++col) {
        col_sums(col) = 0.0;
        for (size_t row = 0; row < 3; ++row) {
            col_sums(col) += matrix(row, col);
        }
    }
    // col_sums = [6.0, 6.0, 6.0, 8.0]

    // Step 2: Broadcast and divide
    Tensor<double> normalized({3, 4});
    auto sums_broadcasted = linalg::broadcast_to(col_sums.view<1>(), std::array<size_t, 2>{3, 4});
    linalg::divide(matrix.view<2>(), sums_broadcasted, normalized.view<2>());

    // Verify: Each column sums to 1.0
    for (size_t col = 0; col < 4; ++col) {
        double sum = 0.0;
        for (size_t row = 0; row < 3; ++row) {
            sum += normalized(row, col);
        }
        REQUIRE_THAT(sum, WithinAbs(1.0, 1e-10));
    }

    // Verify first column
    REQUIRE_THAT(normalized(0, 0), WithinAbs(1.0/6.0, 1e-10));
    REQUIRE_THAT(normalized(1, 0), WithinAbs(2.0/6.0, 1e-10));
    REQUIRE_THAT(normalized(2, 0), WithinAbs(3.0/6.0, 1e-10));
}

TEST_CASE("broadcast_to: scalar to tensor", "[linalg][broadcast]") {
    // Create "scalar" as 1-element tensor
    Tensor<double> scalar({1});
    scalar(0) = 42.0;

    // Broadcast to 2x3 matrix
    auto scalar_view = scalar.view<1>();
    auto broadcasted = linalg::broadcast_to(scalar_view, std::array<size_t, 2>{2, 3});

    // All strides should be 0 (repeat same value)
    REQUIRE(broadcasted.stride(0) == 0);
    REQUIRE(broadcasted.stride(1) == 0);

    // Every element is 42.0
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            REQUIRE(broadcasted[i, j] == 42.0);
        }
    }
}

TEST_CASE("broadcast_to: incompatible shapes throw", "[linalg][broadcast]") {
    Tensor<double> vec({5});
    auto view = vec.view<1>();

    // Valid: [5] -> [3, 5] (prepend dimension)
    REQUIRE_NOTHROW(linalg::broadcast_to(view, std::array<size_t, 2>{3, 5}));

    // Invalid: [5] -> [3, 7] (last dimension mismatch)
    REQUIRE_THROWS_AS(
        linalg::broadcast_to(view, std::array<size_t, 2>{3, 7}),
        std::invalid_argument
    );

    // Invalid: [5] -> [5, 3] (broadcasting to smaller dimension)
    REQUIRE_THROWS_AS(
        linalg::broadcast_to(view, std::array<size_t, 2>{5, 3}),
        std::invalid_argument
    );
}
