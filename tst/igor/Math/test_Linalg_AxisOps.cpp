#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>


using namespace igor::math;
using Catch::Approx;

TEST_CASE("Linalg Axis Operations", "[Math][Linalg][Axis]") {

    SECTION("Sum Axis (Rank 2)") {
        Tensor<double> t({2, 3});
        // [ 1 2 3 ]
        // [ 4 5 6 ]
        t(0, 0) = 1.0; t(0, 1) = 2.0; t(0, 2) = 3.0;
        t(1, 0) = 4.0; t(1, 1) = 5.0; t(1, 2) = 6.0;

        // Sum along Axis 0 (Columns) -> [5, 7, 9]
        auto sum0 = linalg::sum_axis(t, 0);
        REQUIRE(sum0.ndim() == 1);
        REQUIRE(sum0.shape()[0] == 3);
        REQUIRE(sum0(0) == 5.0);
        REQUIRE(sum0(1) == 7.0);
        REQUIRE(sum0(2) == 9.0);

        // Sum along Axis 1 (Rows) -> [6, 15]
        auto sum1 = linalg::sum_axis(t, 1);
        REQUIRE(sum1.ndim() == 1);
        REQUIRE(sum1.shape()[0] == 2);
        REQUIRE(sum1(0) == 6.0);
        REQUIRE(sum1(1) == 15.0);
    }

    SECTION("Sum Axis (Rank 3)") {
        Tensor<double> t({2, 2, 2});
        // filled with 1.0
        for(size_t i=0; i<t.size(); ++i) t.data()[i] = 1.0 * i;

        // Sum Axis 0 -> [4, 6; 8, 10]
        auto sum0 = linalg::sum_axis(t, 0);
        REQUIRE(sum0.ndim() == 2);
        REQUIRE(sum0(0, 0) == 4.0);
        REQUIRE(sum0(0, 1) == 6.0);
        REQUIRE(sum0(1, 0) == 8.0);
        REQUIRE(sum0(1, 1) == 10.0);

        // Sum Axis 1 -> [4, 6; 8, 10]
        auto sum1 = linalg::sum_axis(t, 1);
        REQUIRE(sum1.ndim() == 2);
        REQUIRE(sum1(0, 0) == 2.0);
        REQUIRE(sum1(0, 1) == 4.0);
        REQUIRE(sum1(1, 0) == 10.0);
        REQUIRE(sum1(1, 1) == 12.0);

        // Sum Axis 2 -> [1, 5; 9, 13]
        auto sum2 = linalg::sum_axis(t, 2);
        REQUIRE(sum2.ndim() == 2);
        REQUIRE(sum2(0, 0) == 1.0);
        REQUIRE(sum2(0, 1) == 5.0);
        REQUIRE(sum2(1, 0) == 9.0);
        REQUIRE(sum2(1, 1) == 13.0);
    }

    SECTION("Normalize Axis (Rank 2)") {
        Tensor<double> t({2, 2});
        // [ 1, 3 ] -> Sum=4 -> [0.25, 0.75]
        // [ 2, 2 ] -> Sum=4 -> [0.5, 0.5]
        t(0, 0) = 1.0; t(0, 1) = 3.0;
        t(1, 0) = 2.0; t(1, 1) = 2.0;

        Tensor<double> out({2, 2});

        // Normalize rows (axis 1)
        linalg::normalize_axis(t, out, 1);

        REQUIRE(out(0, 0) == 0.25);
        REQUIRE(out(0, 1) == 0.75);
        REQUIRE(out(1, 0) == 0.5);
        REQUIRE(out(1, 1) == 0.5);

        // Sum of normalized rows should be 1.0
        auto check = linalg::sum_axis(out, 1);
        REQUIRE(check(0) == Approx(1.0));
        REQUIRE(check(1) == Approx(1.0));
    }

    SECTION("Normalize Axis (In-Place)") {
        Tensor<double> t({2, 2});
        t(0, 0) = 10.0; t(0, 1) = 10.0;
        t(1, 0) = 20.0; t(1, 1) = 20.0;

        // Normalize columns (axis 0)
        // Col 0: 10, 20 -> Sum 30 -> 1/3, 2/3
        linalg::normalize_axis(t, t, 0);

        REQUIRE(t(0, 0) == Approx(0.33333333));
        REQUIRE(t(1, 0) == Approx(0.66666666));
        REQUIRE(t(0, 1) == Approx(0.33333333));
        REQUIRE(t(1, 1) == Approx(0.66666666));
    }

    SECTION("Normalize Zero Sum") {
        Tensor<double> t({2});
        t(0) = 0.0; t(1) = 0.0;
        Tensor<double> out({2});

        // Should produce zeros and not NaN
        linalg::normalize_axis(t, out, 0);
        REQUIRE(out(0) == 0.0);
        REQUIRE(out(1) == 0.0);
    }

    SECTION("Rank 1 Special Case") {
        Tensor<double> t({4});
        t(0) = 1.0; t(1) = 2.0; t(2) = 3.0; t(3) = 4.0;

        // Sum Axis 0 (only axis)
        // Note: sum_axis on Rank 1 is not fully supported by generic logic unless we handled it.
        // But normalize_axis IS supported.
        Tensor<double> out({4});
        linalg::normalize_axis(t, out, 0);

        REQUIRE(out(0) == 0.1);
        REQUIRE(out(3) == 0.4);
    }

    SECTION("Error Handling") {
        Tensor<double> t({2, 2});

        REQUIRE_THROWS_AS(linalg::sum_axis(t, 2), std::invalid_argument);
        REQUIRE_THROWS_AS(linalg::normalize_axis(t, t, 2), std::invalid_argument);

        Tensor<double> wrong_shape({2, 3});
        REQUIRE_THROWS_AS(linalg::normalize_axis(t, wrong_shape, 0), std::invalid_argument);
    }

    SECTION("Randomized 3D Tensor Check") {
        constexpr size_t D0 = 4, D1 = 5, D2 = 3;
        Tensor<double> t({D0, D1, D2});

        // Fill with linear index: i*15 + j*3 + k
        for(size_t i=0; i<D0; ++i) {
            for(size_t j=0; j<D1; ++j) {
                for(size_t k=0; k<D2; ++k) {
                    t(i, j, k) = (double)(i*D1*D2 + j*D2 + k);
                }
            }
        }

        // Sum along Axis 2 (Last dim) -> Should be sum of 0,1,2 + offset
        auto sum2 = linalg::sum_axis(t, 2);

        for(size_t i=0; i<D0; ++i) {
            for(size_t j=0; j<D1; ++j) {
                double base = (double)(i*D1*D2 + j*D2);
                // Sum of base, base+1, base+2 = 3*base + 3
                REQUIRE(sum2(i, j) == Approx(3.0*base + 3.0));
            }
        }
        // This loop adds 4*5 = 20 assertions

        // Sum along Axis 1 -> Sum k over j=0..4
        auto sum1 = linalg::sum_axis(t, 1);
        for(size_t i=0; i<D0; ++i) {
            for(size_t k=0; k<D2; ++k) {
                // Sum of (i*15 + j*3 + k) for j=0..4
                // = 5*(i*15 + k) + 3*sum(0..4)
                // = 75*i + 5*k + 3*10 = 75*i + 5*k + 30
                double expected = 75.0*i + 5.0*k + 30.0;
                REQUIRE(sum1(i, k) == Approx(expected));
            }
        }
        // This loop adds 4*3 = 12 assertions
    }
}
