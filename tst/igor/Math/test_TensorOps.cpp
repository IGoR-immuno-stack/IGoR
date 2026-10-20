/*
 * test_TensorOps.cpp
 *
 *  Created on: Feb 11, 2026
 *
 *  Tests for Tensor copy semantics, compound assignment operators,
 *  exact equality, and NumPy-style allclose.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Math/Tensor.h>

using namespace igor::math;
using Catch::Matchers::WithinRel;

// ─── Copy Semantics ────────────────────────────────────────────────────

TEST_CASE("Tensor copy constructor deep copies", "[Math][Tensor][Copy]") {
    Tensor<double> original({3, 4});
    for (std::size_t i = 0; i < original.size(); ++i) {
        original.data()[i] = static_cast<double>(i);
    }

    Tensor<double> copy(original);

    SECTION("copy has same shape") {
        REQUIRE(copy.shape() == original.shape());
        REQUIRE(copy.ndim() == original.ndim());
        REQUIRE(copy.size() == original.size());
    }

    SECTION("copy has same values") {
        for (std::size_t i = 0; i < original.size(); ++i) {
            REQUIRE(copy.data()[i] == original.data()[i]);
        }
    }

    SECTION("copy is independent (different pointer)") {
        REQUIRE(copy.data() != original.data());
    }

    SECTION("copy is always owning") {
        REQUIRE(copy.is_owning());
    }

    SECTION("modifying copy does not affect original") {
        copy.data()[0] = 999.0;
        REQUIRE(original.data()[0] == 0.0);
    }
}

TEST_CASE("Tensor copy assignment deep copies", "[Math][Tensor][Copy]") {
    Tensor<double> original({2, 3});
    for (std::size_t i = 0; i < original.size(); ++i) {
        original.data()[i] = static_cast<double>(i * 10);
    }

    Tensor<double> copy({1});  // Different shape initially
    copy = original;

    SECTION("copy matches original") {
        REQUIRE(copy.shape() == original.shape());
        for (std::size_t i = 0; i < original.size(); ++i) {
            REQUIRE(copy.data()[i] == original.data()[i]);
        }
    }

    SECTION("copy is independent") {
        REQUIRE(copy.data() != original.data());
        REQUIRE(copy.is_owning());
    }

    SECTION("modifying copy does not affect original") {
        copy.data()[0] = 999.0;
        REQUIRE(original.data()[0] == 0.0);
    }
}

TEST_CASE("Tensor copy from borrowing creates owning copy", "[Math][Tensor][Copy]") {
    double raw_data[] = {1.0, 2.0, 3.0, 4.0};
    Tensor<double> borrowing(raw_data, {4});

    REQUIRE_FALSE(borrowing.is_owning());

    Tensor<double> copy(borrowing);

    SECTION("copy is owning") {
        REQUIRE(copy.is_owning());
    }

    SECTION("copy has same values") {
        for (std::size_t i = 0; i < 4; ++i) {
            REQUIRE(copy.data()[i] == raw_data[i]);
        }
    }

    SECTION("copy is independent from raw data") {
        copy.data()[0] = 999.0;
        REQUIRE(raw_data[0] == 1.0);
    }
}

TEST_CASE("Tensor self-assignment is safe", "[Math][Tensor][Copy]") {
    Tensor<double> t({3});
    t(0) = 1.0;
    t(1) = 2.0;
    t(2) = 3.0;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
#endif
    t = t;
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    REQUIRE(t(0) == 1.0);
    REQUIRE(t(1) == 2.0);
    REQUIRE(t(2) == 3.0);
}

// ─── Compound Assignment: Tensor ───────────────────────────────────────

TEST_CASE("Tensor operator+=", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    Tensor<double> b({3});
    b(0) = 10.0; b(1) = 20.0; b(2) = 30.0;

    a += b;

    REQUIRE(a(0) == 11.0);
    REQUIRE(a(1) == 22.0);
    REQUIRE(a(2) == 33.0);

    // b is unchanged
    REQUIRE(b(0) == 10.0);
}

TEST_CASE("Tensor operator-=", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 10.0; a(1) = 20.0; a(2) = 30.0;

    Tensor<double> b({3});
    b(0) = 1.0; b(1) = 2.0; b(2) = 3.0;

    a -= b;

    REQUIRE(a(0) == 9.0);
    REQUIRE(a(1) == 18.0);
    REQUIRE(a(2) == 27.0);
}

TEST_CASE("Tensor operator*= (element-wise)", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 2.0; a(1) = 3.0; a(2) = 4.0;

    Tensor<double> b({3});
    b(0) = 5.0; b(1) = 10.0; b(2) = 0.5;

    a *= b;

    REQUIRE(a(0) == 10.0);
    REQUIRE(a(1) == 30.0);
    REQUIRE(a(2) == 2.0);
}

TEST_CASE("Tensor operator/= (element-wise)", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 10.0; a(1) = 20.0; a(2) = 30.0;

    Tensor<double> b({3});
    b(0) = 2.0; b(1) = 5.0; b(2) = 10.0;

    a /= b;

    REQUIRE(a(0) == 5.0);
    REQUIRE(a(1) == 4.0);
    REQUIRE(a(2) == 3.0);
}

// ─── Compound Assignment: Scalar ───────────────────────────────────────

TEST_CASE("Tensor operator*= (scalar)", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    a *= 10.0;

    REQUIRE(a(0) == 10.0);
    REQUIRE(a(1) == 20.0);
    REQUIRE(a(2) == 30.0);
}

TEST_CASE("Tensor operator/= (scalar)", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 10.0; a(1) = 20.0; a(2) = 30.0;

    a /= 10.0;

    REQUIRE(a(0) == 1.0);
    REQUIRE(a(1) == 2.0);
    REQUIRE(a(2) == 3.0);
}

// ─── Multi-dimensional Operators ───────────────────────────────────────

TEST_CASE("Tensor operators work on 2D tensors", "[Math][Tensor][Ops]") {
    Tensor<double> a({2, 3});
    Tensor<double> b({2, 3});

    for (std::size_t i = 0; i < a.size(); ++i) {
        a.data()[i] = static_cast<double>(i);
        b.data()[i] = static_cast<double>(i * 2);
    }

    a += b;

    // a[i] = i + 2*i = 3*i
    for (std::size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a.data()[i] == static_cast<double>(i * 3));
    }
}

TEST_CASE("Tensor operators return *this for chaining", "[Math][Tensor][Ops]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    Tensor<double> b({3});
    b(0) = 1.0; b(1) = 1.0; b(2) = 1.0;

    // Chain: (a += b) *= 2.0
    (a += b) *= 2.0;

    REQUIRE(a(0) == 4.0);  // (1+1)*2
    REQUIRE(a(1) == 6.0);  // (2+1)*2
    REQUIRE(a(2) == 8.0);  // (3+1)*2
}

// ─── Exact Equality ────────────────────────────────────────────────────

TEST_CASE("Tensor operator== exact equality", "[Math][Tensor][Compare]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    Tensor<double> b({3});
    b(0) = 1.0; b(1) = 2.0; b(2) = 3.0;

    SECTION("equal tensors") {
        REQUIRE(a == b);
        REQUIRE_FALSE(a != b);
    }

    SECTION("different values") {
        b(1) = 2.0001;
        REQUIRE_FALSE(a == b);
        REQUIRE(a != b);
    }

    SECTION("different shapes") {
        Tensor<double> c({4});
        REQUIRE_FALSE(a == c);
    }

    SECTION("copy is equal") {
        Tensor<double> c(a);
        REQUIRE(a == c);
    }
}

TEST_CASE("Tensor operator== on 2D tensors", "[Math][Tensor][Compare]") {
    Tensor<double> a({2, 3});
    Tensor<double> b({2, 3});

    for (std::size_t i = 0; i < a.size(); ++i) {
        a.data()[i] = static_cast<double>(i);
        b.data()[i] = static_cast<double>(i);
    }

    REQUIRE(a == b);

    b(1, 2) = 999.0;
    REQUIRE(a != b);
}

// ─── allclose (NumPy Convention) ───────────────────────────────────────

TEST_CASE("allclose with identical tensors", "[Math][Tensor][Compare]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    Tensor<double> b(a);  // exact copy

    REQUIRE(allclose(a, b));
}

TEST_CASE("allclose with near-equal values", "[Math][Tensor][Compare]") {
    Tensor<double> a({3});
    a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;

    Tensor<double> b(a);
    b(0) += 1e-09;  // Within default tolerance

    SECTION("within default tolerance") {
        REQUIRE(allclose(a, b));
    }

    SECTION("outside tight tolerance") {
        REQUIRE_FALSE(allclose(a, b, 0.0, 1e-10));
    }
}

TEST_CASE("allclose with different shapes returns false", "[Math][Tensor][Compare]") {
    Tensor<double> a({3});
    Tensor<double> b({4});
    REQUIRE_FALSE(allclose(a, b));
}

TEST_CASE("allclose relative tolerance", "[Math][Tensor][Compare]") {
    Tensor<double> a({2});
    a(0) = 1.0;
    a(1) = 1000.0;

    Tensor<double> b({2});
    b(0) = 1.0;
    b(1) = 1000.01;

    // |1000.01 - 1000.0| = 0.01
    // rtol * |b[1]| = 1e-05 * 1000.01 = 0.01000...
    // atol = 1e-08
    // tol = 0.01000... + 1e-08 ≈ 0.01000
    // 0.01 <= 0.01000 → close enough
    REQUIRE(allclose(a, b));

    // With tighter tolerance
    REQUIRE_FALSE(allclose(a, b, 1e-06, 0.0));
}

TEST_CASE("allclose with long double", "[Math][Tensor][Compare]") {
    Tensor<long double> a({3});
    a(0) = 1.0L; a(1) = 2.0L; a(2) = 3.0L;

    Tensor<long double> b(a);
    b(0) += static_cast<long double>(1e-12);

    REQUIRE(allclose(a, b));
    REQUIRE_FALSE(allclose(a, b, 0.0, 1e-14));
}

TEST_CASE("allclose on 2D tensor", "[Math][Tensor][Compare]") {
    Tensor<double> a({3, 4});
    Tensor<double> b({3, 4});

    for (std::size_t i = 0; i < a.size(); ++i) {
        a.data()[i] = static_cast<double>(i) * 0.1;
        b.data()[i] = static_cast<double>(i) * 0.1 + 1e-10;
    }

    REQUIRE(allclose(a, b));
}
