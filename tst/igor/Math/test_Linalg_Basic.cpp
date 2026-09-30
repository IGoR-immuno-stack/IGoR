#include <catch2/catch_test_macros.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>
#include "TestHelpers.h"
#include <vector>

using namespace igor::math;

TEST_CASE("Linalg Basic Operations", "[Math][Linalg]") {

    SECTION("Element-wise Arithmetic") {
        Tensor<double> a({2, 2});
        Tensor<double> b({2, 2});
        Tensor<double> c({2, 2});

        test::fill_2d(a, [](size_t i, size_t j) { return 1.0 + i * 2 + j; });
        test::fill_2d(b, [](size_t i, size_t j) { return 4.0 - i * 2 - j; });

        // Add
        linalg::add(a.view<2>(), b.view<2>(), c.view<2>());
        REQUIRE(c(0, 0) == 5.0);
        REQUIRE(c(0, 1) == 5.0);
        REQUIRE(c(1, 0) == 5.0);
        REQUIRE(c(1, 1) == 5.0);

        // Subtract
        linalg::subtract(a.view<2>(), b.view<2>(), c.view<2>());
        REQUIRE(c(0, 0) == -3.0);

        // Multiply
        linalg::multiply(a.view<2>(), b.view<2>(), c.view<2>());
        REQUIRE(c(0, 0) == 4.0);
        REQUIRE(c(1, 1) == 4.0);

        // Divide
        linalg::divide(a.view<2>(), b.view<2>(), c.view<2>());
        REQUIRE(c(0, 0) == 0.25);
    }

    SECTION("Scale") {
        auto a = test::make_sequential<double>({3}, 1.0);
        Tensor<double> out({3});

        linalg::scale(a.view<1>(), 2.0, out.view<1>());
        REQUIRE(out(0) == 2.0);
        REQUIRE(out(1) == 4.0);
        REQUIRE(out(2) == 6.0);
    }

    SECTION("Reductions") {
        Tensor<double> t({2, 2});
        t(0, 0) = 1.0; t(0, 1) = 5.0;
        t(1, 0) = 2.0; t(1, 1) = -1.0;

        REQUIRE(linalg::sum(t.view<2>()) == 7.0);
        REQUIRE(linalg::max(t.view<2>()) == 5.0);
        REQUIRE(linalg::min(t.view<2>()) == -1.0);

        auto idx = linalg::argmax(t.view<2>());
        // Expecting {0, 1} flat index? No, array<size_t, 2>
        REQUIRE(idx[0] == 0);
        REQUIRE(idx[1] == 1);
    }

    SECTION("Tensor API (Runtime Dispatch)") {
        Tensor<double> a({2, 2});
        Tensor<double> b({2, 2});
        Tensor<double> c({2, 2});

        test::fill_2d(a, [](size_t i, size_t j) { return 1.0 + i * 2 + j; });
        test::fill_2d(b, [](size_t i, size_t j) { return 4.0 - i * 2 - j; });

        // Add
        linalg::add(a, b, c);
        REQUIRE(c(0, 0) == 5.0);
        REQUIRE(c(0, 1) == 5.0);
        REQUIRE(c(1, 1) == 5.0);

        // Subtract
        linalg::subtract(a, b, c);
        REQUIRE(c(0, 0) == -3.0);

        // Multiply
        linalg::multiply(a, b, c);
        REQUIRE(c(0, 1) == 6.0);

        // Scale
        Tensor<double> out({2, 2});
        linalg::scale(a, 2.0, out);
        REQUIRE(out(1, 1) == 8.0);

        // Sum
        REQUIRE(linalg::sum(a) == 10.0);

        // Max
        REQUIRE(linalg::max(a) == 4.0);

        // Argmax
        auto max_idx = linalg::argmax(a);
        REQUIRE(max_idx[0] == 1);
        REQUIRE(max_idx[1] == 1);

        // Rank Mismatch (always checked)
        Tensor<double> bad({4});
        REQUIRE_THROWS_AS(linalg::add(a, bad, c), std::invalid_argument);

#ifdef IGOR_MATH_DEBUG
        // Shape Mismatch (only checked in debug mode)
        Tensor<double> wrong_shape({3, 3});
        REQUIRE_THROWS_AS(linalg::add(a, wrong_shape, c), std::invalid_argument);
        REQUIRE_THROWS_AS(linalg::scale(a, 2.0, wrong_shape), std::invalid_argument);
#endif
    }
}
