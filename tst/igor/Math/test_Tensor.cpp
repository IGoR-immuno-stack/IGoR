/*
 * test_Tensor.cpp
 *
 *  Created on: Feb 06, 2026
 *      Author: IGoR Agent
 */

#include <catch2/catch_test_macros.hpp>
#include <igor/Math/Tensor.h>
<<<<<<< HEAD
#include "TestHelpers.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <span>
#include <numeric>     // For std::iota, std::accumulate
#include <type_traits> // For std::remove_cvref_t
=======
#include <vector>

>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)
using namespace igor::math;

TEST_CASE("Tensor Basic Operations", "[Math][Tensor]") {

    SECTION("2D Tensor Runtime Shape") {
        std::vector<size_t> shape{4, 3}; // 4 rows, 3 cols
        Tensor<double> tensor(shape);

<<<<<<< HEAD
        test::require_shape(tensor, {4, 3});
=======
        REQUIRE(tensor.size() == 12);
        REQUIRE(tensor.rank() == 2);
        REQUIRE(tensor.shape()[0] == 4);
        REQUIRE(tensor.shape()[1] == 3);
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)

        // Access via variadic operator()
        tensor(1, 1) = 42.0;

        REQUIRE(tensor.data()[1 * 3 + 1] == 42.0); // Assuming row-major

        // Also test span-based access
        std::vector<size_t> idx{1, 1};
        std::span<const size_t> idx_span{idx};
        REQUIRE(tensor(idx_span) == 42.0);
    }

    SECTION("3D Tensor Borrowing") {
        std::vector<int> raw(2 * 2 * 2, 0);
        std::vector<size_t> shape{2, 2, 2};

        Tensor<int> tensor(raw.data(), shape);

<<<<<<< HEAD
        test::require_shape(tensor, {2, 2, 2});
=======
        REQUIRE(tensor.size() == 8);
        REQUIRE(tensor.rank() == 3);
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)

        // Modify via view<3>
        auto v = tensor.view<3>();
        v[0, 1, 0] = 7;

        REQUIRE(raw[2] == 7); // 0*4 + 1*2 + 0 = 2
        REQUIRE(v[0, 1, 0] == 7);
    }

    SECTION("Runtime Rank Usage") {
        // Simulating runtime configuration (e.g. from file)
        size_t runtime_rows = 4;
        size_t runtime_cols = 3;
        std::vector<size_t> runtime_shape = {runtime_rows, runtime_cols};

        // Tensor construction (Rank is distinct from type)
        Tensor<double> t(runtime_shape);

<<<<<<< HEAD
        test::require_shape(t, {4, 3});
=======
        REQUIRE(t.rank() == 2);
        REQUIRE(t.size() == 12);
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)

        // Access via runtime indices (std::span)
        std::vector<size_t> idx = {1, 2};
        std::span<const size_t> idx_span{idx};
        t(idx_span) = 99.9;

        // Verify via raw pointer
        // Index: 1*3 + 2 = 5
        REQUIRE(t.data()[5] == 99.9);

        // Verify via runtime access
        REQUIRE(t(idx_span) == 99.9);
    }

    SECTION("Variadic Access") {
        Tensor<double> t({3, 4, 5});

        // Variadic operator()
        t(1, 2, 3) = 123.0;
        REQUIRE(t(1, 2, 3) == 123.0);

        // Also accessible via view
        auto v = t.view<3>();
        REQUIRE(v[1, 2, 3] == 123.0);
    }

    SECTION("Rank Mismatch Error") {
        Tensor<double> t({3, 4});

        // Wrong number of indices
        REQUIRE_THROWS_AS(t(0), std::runtime_error);
        REQUIRE_THROWS_AS(t(0, 1, 2), std::runtime_error);

        // Wrong rank in view request
        REQUIRE_THROWS_AS(t.view<3>(), std::runtime_error);
    }

    SECTION("Const Correctness") {
        std::vector<size_t> shape{3, 4};
        Tensor<double> t(shape);
        t(1, 2) = 42.0;

        const Tensor<double>& ct = t;

        // Const access should work
        REQUIRE(ct(1, 2) == 42.0);

        // Const view should work
        auto cv = ct.view<2>();
        REQUIRE(cv[1, 2] == 42.0);
    }

    SECTION("Rank > 5 Fallback") {
        // Create 6D tensor (exceeds variant range)
        Tensor<double> t({2, 2, 2, 2, 2, 2});

<<<<<<< HEAD
        test::require_shape(t, {2, 2, 2, 2, 2, 2});
=======
        REQUIRE(t.rank() == 6);
        REQUIRE(t.size() == 64);
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)

        // Variadic access falls back to manual strides
        t(0, 0, 0, 0, 0, 0) = 99.0;
        REQUIRE(t(0, 0, 0, 0, 0, 0) == 99.0);
        REQUIRE(t.data()[0] == 99.0);

        // Last element
        t(1, 1, 1, 1, 1, 1) = 77.0;
        REQUIRE(t.data()[63] == 77.0);

        // Runtime span access also works
        std::vector<size_t> idx{0, 1, 0, 1, 0, 1};
        std::span<const size_t> idx_span{idx};
        t(idx_span) = 55.0;
        REQUIRE(t(idx_span) == 55.0);
    }

    SECTION("Move Semantics") {
        Tensor<double> source({10, 10});
        source(5, 5) = 42.0;

        Tensor<double> dest = std::move(source);

<<<<<<< HEAD
        test::require_shape(dest, {10, 10});
=======
        REQUIRE(dest.size() == 100);
        REQUIRE(dest.rank() == 2);
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)
        REQUIRE(dest(5, 5) == 42.0);

        // Source should be in moved-from state (valid but unspecified)
        // We can at least check it doesn't crash
        REQUIRE(source.size() == 0);
    }
<<<<<<< HEAD

    SECTION("Initializer List Access") {
        Tensor<double> t({2, 2});
        t({0, 0}) = 10.0;
        t({0, 1}) = 20.0;

        REQUIRE(t({0, 0}) == 10.0);
        REQUIRE(t({0, 1}) == 20.0);

        // Verify via variadic
        REQUIRE(t(0, 0) == 10.0);
        REQUIRE(t(0, 1) == 20.0);

        // Mixed access patterns work correctly
        t({1, 0}) = 30.0;
        REQUIRE(t(1, 0) == 30.0);
    }

    SECTION("Iterators") {
        Tensor<int> t({2, 3});
        // Fill: 0, 1, 2, 3, 4, 5
        std::iota(t.begin(), t.end(), 0);

        // Verify order (Row-Major)
        REQUIRE(t(0, 0) == 0);
        REQUIRE(t(0, 1) == 1);
        REQUIRE(t(0, 2) == 2);
        REQUIRE(t(1, 0) == 3);
        REQUIRE(t(1, 1) == 4);
        REQUIRE(t(1, 2) == 5);

        // Verify const iterators
        const auto& ct = t;
        int sum = std::accumulate(ct.cbegin(), ct.cend(), 0);
        REQUIRE(sum == 15); // 0+1+2+3+4+5
    }

    SECTION("Visitor Pattern (apply)") {
        Tensor<double> t1({10});
        Tensor<double> t2({5, 5});

        t1(0) = 1.0;
        t2(0, 0) = 2.0;

        auto sum_first = [](auto&& view) {
             // Access first element regardless of rank
             // Need proper indexing based on rank?
             // Or just use flattened access if not available?
             // mdspan doesn't support flat access via [] without multi-index.
             // But we can check rank.
             if constexpr (std::remove_cvref_t<decltype(view)>::rank() == 1) {
                 return view[0];
             } else if constexpr (std::remove_cvref_t<decltype(view)>::rank() == 2) {
                 return view[0, 0];
             }
             return 0.0;
        };

        REQUIRE(t1.apply(sum_first) == 1.0);
        REQUIRE(t2.apply(sum_first) == 2.0);

        // Test generic mutation (Rank 1 specific)
        t1.apply([](auto&& v) {
            if constexpr (std::remove_cvref_t<decltype(v)>::rank() == 1) {
                v[0] = 99.0;
            }
        });
        REQUIRE(t1(0) == 99.0);
    }

    SECTION("Manual std::visit Usage") {
        Tensor<double> t({2, 2});
        t(0, 0) = 42.0;

        // Get the variant explicitly
        auto v = t.variant();

        bool visited = false;
        std::visit([&](auto&& view) {
            visited = true;
            using ViewType = std::remove_cvref_t<decltype(view)>;

            // We know t is Rank 2
            if constexpr (ViewType::rank() == 2) {
                REQUIRE(view[0, 0] == 42.0);
                // Verify we can mutate if non-const
                view[0, 1] = 84.0;
            } else {
                FAIL("Variant holds wrong rank view");
            }
        }, v);

        REQUIRE(visited);
        REQUIRE(t(0, 1) == 84.0);

        // Verify const variant
        const Tensor<double>& ct = t;
        std::visit([&](auto&& view) {
            using ViewType = std::remove_cvref_t<decltype(view)>;
            if constexpr (ViewType::rank() == 2) {
                REQUIRE(view[0, 1] == 84.0);
                // view[0, 0] = 0; // Should not compile (read-only view)
            }
        }, ct.variant());
    }

    SECTION("Shape and Size Info") {
        Tensor<float> t({2, 4, 3});
        test::require_shape(t, {2, 4, 3});

        REQUIRE(t.data() != nullptr);
    }
=======
>>>>>>> a9fb84b (feat(math): implement Tensor class with type-erased rank and mdspan views)
}
