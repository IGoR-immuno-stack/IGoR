#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/Linalg.h>
#include <cmath>
#include <limits>

using namespace igor::math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Linalg Probabilistic Operations", "[Math][Linalg][Probabilistic]") {

    SECTION("Dot Product - Basic") {
        Tensor<double> v1({3}), v2({3});
        v1(0) = 1.0; v1(1) = 2.0; v1(2) = 3.0;
        v2(0) = 4.0; v2(1) = 5.0; v2(2) = 6.0;
        
        // mdspan interface
        auto result = linalg::dot(v1.view<1>(), v2.view<1>());
        REQUIRE(result == 32.0); // 1*4 + 2*5 + 3*6 = 32
        
        // Tensor interface
        auto result_tensor = linalg::dot(v1, v2);
        REQUIRE(result_tensor == 32.0);
    }
    
    SECTION("Dot Product - Dimension Mismatch") {
        Tensor<double> v1({3}), v2({4});
        REQUIRE_THROWS_AS(linalg::dot(v1.view<1>(), v2.view<1>()), std::invalid_argument);
        REQUIRE_THROWS_AS(linalg::dot(v1, v2), std::invalid_argument);
    }
    
    SECTION("Dot Product - Rank Mismatch") {
        Tensor<double> v1({3}), m2({2, 2});
        REQUIRE_THROWS_AS(linalg::dot(v1, m2), std::invalid_argument);
    }
    
    SECTION("Matrix Multiplication - Basic 2x2") {
        Tensor<double> A({2, 2}), B({2, 2}), C({2, 2});
        
        // A = [[1, 2],
        //      [3, 4]]
        A(0, 0) = 1.0; A(0, 1) = 2.0;
        A(1, 0) = 3.0; A(1, 1) = 4.0;
        
        // B = [[5, 6],
        //      [7, 8]]
        B(0, 0) = 5.0; B(0, 1) = 6.0;
        B(1, 0) = 7.0; B(1, 1) = 8.0;
        
        // mdspan interface
        linalg::matmul(A.view<2>(), B.view<2>(), C.view<2>());
        
        // C = [[1*5+2*7, 1*6+2*8],   = [[19, 22],
        //      [3*5+4*7, 3*6+4*8]]      [43, 50]]
        REQUIRE(C(0, 0) == 19.0);
        REQUIRE(C(0, 1) == 22.0);
        REQUIRE(C(1, 0) == 43.0);
        REQUIRE(C(1, 1) == 50.0);
        
        // Tensor interface
        Tensor<double> D({2, 2});
        linalg::matmul(A, B, D);
        REQUIRE(D(0, 0) == 19.0);
        REQUIRE(D(1, 1) == 50.0);
    }
    
    SECTION("Matrix Multiplication - Non-Square") {
        Tensor<double> A({2, 3}), B({3, 4}), C({2, 4});
        
        // Initialize A and B with simple pattern
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                A(i, j) = i * 3 + j + 1;
            }
        }
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                B(i, j) = i * 4 + j + 1;
            }
        }
        
        linalg::matmul(A, B, C);
        
        // Verify a few elements
        REQUIRE(C(0, 0) == 38.0);  // 1*1 + 2*5 + 3*9 = 38
        REQUIRE(C(1, 3) == 128.0); // 4*4 + 5*8 + 6*12 = 128
    }
    
    SECTION("Matrix Multiplication - Identity") {
        Tensor<double> A({3, 3}), I({3, 3}), C({3, 3});
        
        // A = some matrix
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                A(i, j) = i * 3 + j + 1;
            }
        }
        
        // I = identity
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                I(i, j) = (i == j) ? 1.0 : 0.0;
            }
        }
        
        linalg::matmul(A, I, C);
        
        // C should equal A
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                REQUIRE(C(i, j) == A(i, j));
            }
        }
    }
    
    SECTION("Matrix Multiplication - Dimension Mismatch") {
        Tensor<double> A({2, 3}), B({2, 4}), C({2, 4});
        // K dimension mismatch: A[2,3] x B[2,4] (should be [3,?])
        REQUIRE_THROWS_AS(linalg::matmul(A.view<2>(), B.view<2>(), C.view<2>()), std::invalid_argument);
        REQUIRE_THROWS_AS(linalg::matmul(A, B, C), std::invalid_argument);
    }
    
    SECTION("Matrix Multiplication - Output Size Mismatch") {
        Tensor<double> A({2, 3}), B({3, 4}), C({3, 3}); // C wrong size
        REQUIRE_THROWS_AS(linalg::matmul(A.view<2>(), B.view<2>(), C.view<2>()), std::invalid_argument);
    }
    
    SECTION("Matrix Multiplication - Rank Mismatch") {
        Tensor<double> A({2, 2}), B({2, 2}), C({2});
        REQUIRE_THROWS_AS(linalg::matmul(A, B, C), std::invalid_argument);
    }
    
    SECTION("Normalize - Basic") {
        Tensor<double> t({3}), out({3});
        t(0) = 1.0; t(1) = 2.0; t(2) = 3.0; // sum = 6
        
        // mdspan interface
        linalg::normalize(t.view<1>(), out.view<1>());
        REQUIRE_THAT(out(0), WithinAbs(1.0/6, 1e-10));
        REQUIRE_THAT(out(1), WithinAbs(2.0/6, 1e-10));
        REQUIRE_THAT(out(2), WithinAbs(3.0/6, 1e-10));
        REQUIRE_THAT(linalg::sum(out.view<1>()), WithinAbs(1.0, 1e-10));
        
        // Tensor interface
        Tensor<double> out2({3});
        linalg::normalize(t, out2);
        REQUIRE_THAT(linalg::sum(out2), WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Normalize - 2D Tensor") {
        Tensor<double> t({2, 3}), out({2, 3});
        for (size_t i = 0; i < 6; ++i) {
            t.data()[i] = i + 1; // 1,2,3,4,5,6 -> sum = 21
        }
        
        linalg::normalize(t, out);
        REQUIRE_THAT(linalg::sum(out), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(out(0, 0), WithinAbs(1.0/21, 1e-10));
        REQUIRE_THAT(out(1, 2), WithinAbs(6.0/21, 1e-10));
    }
    
    SECTION("Normalize - Zero Sum") {
        Tensor<double> zeros({3}), out({3});
        zeros(0) = 0.0; zeros(1) = 0.0; zeros(2) = 0.0;
        
        // Should throw exception
        REQUIRE_THROWS_AS(linalg::normalize(zeros.view<1>(), out.view<1>()), std::invalid_argument);
        REQUIRE_THROWS_AS(linalg::normalize(zeros, out), std::invalid_argument);
    }
    
    SECTION("Normalize - In-Place") {
        Tensor<double> t({3});
        t(0) = 2.0; t(1) = 3.0; t(2) = 5.0; // sum = 10
        
        linalg::normalize(t, t); // in-place
        REQUIRE_THAT(t(0), WithinAbs(0.2, 1e-10));
        REQUIRE_THAT(t(1), WithinAbs(0.3, 1e-10));
        REQUIRE_THAT(t(2), WithinAbs(0.5, 1e-10));
    }
}

TEST_CASE("Linalg Log-Space Operations", "[Math][Linalg][LogSpace]") {
    
    SECTION("Log Add Exp - Basic") {
        REQUIRE_THAT(linalg::log_add_exp(0.0, 0.0), WithinAbs(std::log(2.0), 1e-10));
        REQUIRE_THAT(linalg::log_add_exp(1.0, 2.0), WithinAbs(std::log(std::exp(1.0) + std::exp(2.0)), 1e-10));
    }
    
    SECTION("Log Add Exp - Infinity Handling") {
        double neg_inf = -std::numeric_limits<double>::infinity();
        REQUIRE(linalg::log_add_exp(neg_inf, 5.0) == 5.0);
        REQUIRE(linalg::log_add_exp(5.0, neg_inf) == 5.0);
        REQUIRE(linalg::log_add_exp(neg_inf, neg_inf) == neg_inf);
    }
    
    SECTION("Log Add Exp - Stability (Large Values)") {
        // These would overflow if computed naively as log(exp(1000) + exp(1001))
        double result = linalg::log_add_exp(1000.0, 1001.0);
        REQUIRE(std::isfinite(result));
        REQUIRE(result > 1001.0); // Should be slightly more than max
    }
    
    SECTION("Center - Basic") {
        Tensor<double> t({4}), out({4});
        t(0) = 1.0; t(1) = 5.0; t(2) = 3.0; t(3) = 2.0; // max = 5
        
        // mdspan interface
        linalg::center(t.view<1>(), out.view<1>());
        REQUIRE(out(0) == -4.0); // 1 - 5
        REQUIRE(out(1) == 0.0);  // 5 - 5
        REQUIRE(out(2) == -2.0); // 3 - 5
        REQUIRE(out(3) == -3.0); // 2 - 5
        
        // Tensor interface
        Tensor<double> out2({4});
        linalg::center(t, out2);
        REQUIRE(out2(1) == 0.0);
        REQUIRE(out2(0) == -4.0);
    }
    
    SECTION("Center - 2D Tensor") {
        Tensor<double> t({2, 2}), out({2, 2});
        t(0, 0) = 1.0; t(0, 1) = 5.0;
        t(1, 0) = 3.0; t(1, 1) = 2.0;
        
        linalg::center(t, out);
        REQUIRE(linalg::max(out) == 0.0); // max element becomes 0
        REQUIRE(out(0, 0) == -4.0);
    }
    
    SECTION("Log Sum Exp - Basic") {
        Tensor<double> t({3});
        t(0) = 0.0; t(1) = 0.0; t(2) = 0.0; // exp(0) + exp(0) + exp(0) = 3
        
        // mdspan interface
        auto result = linalg::log_sum_exp(t.view<1>());
        REQUIRE_THAT(result, WithinAbs(std::log(3.0), 1e-10));
        
        // Tensor interface
        auto result2 = linalg::log_sum_exp(t);
        REQUIRE_THAT(result2, WithinAbs(std::log(3.0), 1e-10));
    }
    
    SECTION("Log Sum Exp - Stability (Large Values)") {
        Tensor<double> large({3});
        large(0) = 1000.0; large(1) = 1001.0; large(2) = 1002.0;
        
        // Should not overflow (naive exp would overflow)
        auto lse = linalg::log_sum_exp(large);
        REQUIRE(std::isfinite(lse));
        REQUIRE(lse > 1002.0); // Should be slightly more than max
        REQUIRE(lse < 1003.0); // But not too much more
    }
    
    SECTION("Log Sum Exp - All Negative Infinity") {
        Tensor<double> t({3});
        double neg_inf = -std::numeric_limits<double>::infinity();
        t(0) = neg_inf; t(1) = neg_inf; t(2) = neg_inf;
        
        auto result = linalg::log_sum_exp(t);
        REQUIRE(result == neg_inf);
    }
    
    SECTION("Log Sum Exp - Mixed Values") {
        Tensor<double> t({4});
        t(0) = -1.0; t(1) = 0.0; t(2) = 1.0; t(3) = 2.0;
        
        // Manual calculation: log(exp(-1) + exp(0) + exp(1) + exp(2))
        double expected = std::log(std::exp(-1.0) + std::exp(0.0) + std::exp(1.0) + std::exp(2.0));
        auto result = linalg::log_sum_exp(t);
        REQUIRE_THAT(result, WithinAbs(expected, 1e-10));
    }
    
    SECTION("Log Normalize - Basic") {
        Tensor<double> t({3}), out({3});
        t(0) = 1.0; t(1) = 2.0; t(2) = 3.0;
        
        // mdspan interface
        linalg::log_normalize(t.view<1>(), out.view<1>());
        
        // exp(out) should sum to 1
        Tensor<double> exp_out({3});
        for (size_t i = 0; i < 3; ++i) {
            exp_out(i) = std::exp(out(i));
        }
        REQUIRE_THAT(linalg::sum(exp_out.view<1>()), WithinAbs(1.0, 1e-10));
        
        // Tensor interface
        Tensor<double> out2({3});
        linalg::log_normalize(t, out2);
        for (size_t i = 0; i < 3; ++i) {
            exp_out(i) = std::exp(out2(i));
        }
        REQUIRE_THAT(linalg::sum(exp_out), WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Log Normalize - Stability (Large Values)") {
        Tensor<double> large({3}), out({3});
        large(0) = 1000.0; large(1) = 1001.0; large(2) = 1002.0;
        
        linalg::log_normalize(large, out);
        
        // All values should be finite
        REQUIRE(std::isfinite(out(0)));
        REQUIRE(std::isfinite(out(1)));
        REQUIRE(std::isfinite(out(2)));
        
        // exp(out) should sum to 1
        double sum_exp = std::exp(out(0)) + std::exp(out(1)) + std::exp(out(2));
        REQUIRE_THAT(sum_exp, WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Log Normalize - 2D Tensor") {
        Tensor<double> t({2, 2}), out({2, 2});
        t(0, 0) = 0.0; t(0, 1) = 1.0;
        t(1, 0) = 2.0; t(1, 1) = 3.0;
        
        linalg::log_normalize(t, out);
        
        // exp(out) should sum to 1
        double sum_exp = 0.0;
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                sum_exp += std::exp(out(i, j));
            }
        }
        REQUIRE_THAT(sum_exp, WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Log Normalize - In-Place") {
        Tensor<double> t({3});
        t(0) = 1.0; t(1) = 2.0; t(2) = 3.0;
        
        linalg::log_normalize(t, t); // in-place
        
        // exp(t) should sum to 1
        double sum_exp = std::exp(t(0)) + std::exp(t(1)) + std::exp(t(2));
        REQUIRE_THAT(sum_exp, WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Log Normalize - Preserves Relative Ordering") {
        Tensor<double> t({3}), out({3});
        t(0) = 1.0; t(1) = 2.0; t(2) = 3.0;
        
        linalg::log_normalize(t, out);
        
        // Order should be preserved: out[0] < out[1] < out[2]
        REQUIRE(out(0) < out(1));
        REQUIRE(out(1) < out(2));
    }
}
