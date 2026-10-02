#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/TensorCreation.h>
#include "TestHelpers.h"

using namespace igor::math;
namespace tn = igor::math::tensor;

TEST_CASE("tensor::zeros", "[Math][TensorCreation]") {
    SECTION("1D tensor") {
        auto t = tn::zeros({5});
        test::require_shape(t, {5});
        test::require_all_equal(t, 0.0);
    }

    SECTION("2D tensor") {
        auto t = tn::zeros({3, 4});
        test::require_shape(t, {3, 4});
        test::require_all_equal(t, 0.0);
    }

    SECTION("With explicit type") {
        auto t = tn::zeros<float>({2, 2});
        test::require_shape(t, {2, 2});
        test::require_all_equal(t, 0.0f);
    }
}

TEST_CASE("tensor::ones", "[Math][TensorCreation]") {
    SECTION("1D tensor") {
        auto t = tn::ones({4});
        test::require_shape(t, {4});
        test::require_all_equal(t, 1.0);
    }

    SECTION("3D tensor") {
        auto t = tn::ones({2, 3, 4});
        test::require_shape(t, {2, 3, 4});
        test::require_all_equal(t, 1.0);
    }
}

TEST_CASE("tensor::full", "[Math][TensorCreation]") {
    SECTION("Fill with custom value") {
        auto t = tn::full({3, 3}, 5.5);
        test::require_shape(t, {3, 3});
        test::require_all_equal(t, 5.5);
    }

    SECTION("Fill with integer") {
        auto t = tn::full<int>({2, 2}, 42);
        test::require_shape(t, {2, 2});
        test::require_all_equal(t, 42);
    }
}

TEST_CASE("tensor::arange", "[Math][TensorCreation]") {
    SECTION("Simple range [0, n)") {
        auto t = tn::arange(5.0);
        REQUIRE(t.ndim() == 1);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 0.0);
        REQUIRE(t.data()[1] == 1.0);
        REQUIRE(t.data()[2] == 2.0);
        REQUIRE(t.data()[3] == 3.0);
        REQUIRE(t.data()[4] == 4.0);
    }

    SECTION("Range with start and stop") {
        auto t = tn::arange(2.0, 7.0);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 2.0);
        REQUIRE(t.data()[4] == 6.0);
    }

    SECTION("Range with step") {
        auto t = tn::arange(0.0, 10.0, 2.0);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 0.0);
        REQUIRE(t.data()[1] == 2.0);
        REQUIRE(t.data()[2] == 4.0);
        REQUIRE(t.data()[3] == 6.0);
        REQUIRE(t.data()[4] == 8.0);
    }

    SECTION("Negative step") {
        auto t = tn::arange(10.0, 0.0, -2.0);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 10.0);
        REQUIRE(t.data()[4] == 2.0);
    }

    SECTION("Empty range") {
        auto t = tn::arange(5.0, 2.0, 1.0);  // start > stop with positive step
        REQUIRE(t.shape()[0] == 0);
    }

    SECTION("Integer range") {
        auto t = tn::arange<int>(0, 5);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 0);
        REQUIRE(t.data()[4] == 4);
    }

    SECTION("Zero step throws") {
        REQUIRE_THROWS_AS(tn::arange(0.0, 10.0, 0.0), std::invalid_argument);
    }
}

TEST_CASE("tensor::linspace", "[Math][TensorCreation]") {
    SECTION("Basic linspace with endpoint") {
        auto t = tn::linspace(0.0, 1.0, 5);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 0.0);
        REQUIRE_THAT(t.data()[1], Catch::Matchers::WithinRel(0.25, 1e-10));
        REQUIRE_THAT(t.data()[2], Catch::Matchers::WithinRel(0.5, 1e-10));
        REQUIRE_THAT(t.data()[3], Catch::Matchers::WithinRel(0.75, 1e-10));
        REQUIRE(t.data()[4] == 1.0);
    }

    SECTION("Linspace without endpoint") {
        auto t = tn::linspace(0.0, 1.0, 5, false);
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 0.0);
        REQUIRE_THAT(t.data()[1], Catch::Matchers::WithinRel(0.2, 1e-10));
        REQUIRE_THAT(t.data()[4], Catch::Matchers::WithinRel(0.8, 1e-10));
    }

    SECTION("Single point") {
        auto t = tn::linspace(5.0, 10.0, 1);
        REQUIRE(t.shape()[0] == 1);
        REQUIRE(t.data()[0] == 5.0);
    }

    SECTION("Empty linspace") {
        auto t = tn::linspace(0.0, 1.0, 0);
        REQUIRE(t.shape()[0] == 0);
    }
}

TEST_CASE("tensor::eye", "[Math][TensorCreation]") {
    SECTION("Square identity matrix") {
        auto I = tn::eye(3);
        test::require_shape(I, {3, 3});

        auto view = I.view<2>();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                REQUIRE(view[i, j] == (i == j ? 1.0 : 0.0));
            }
        }
    }

    SECTION("Rectangular identity") {
        auto E = tn::eye(3, 5);
        test::require_shape(E, {3, 5});

        auto view = E.view<2>();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                REQUIRE(view[i, j] == (i == j ? 1.0 : 0.0));
            }
        }
    }

    SECTION("Tall rectangular") {
        auto E = tn::eye(5, 3);
        test::require_shape(E, {5, 3});

        auto view = E.view<2>();
        // Only first 3 diagonal elements should be 1
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                REQUIRE(view[i, j] == (i == j ? 1.0 : 0.0));
            }
        }
    }
}

TEST_CASE("tensor::identity", "[Math][TensorCreation]") {
    auto I = tn::identity(4);
    test::require_shape(I, {4, 4});

    auto view = I.view<2>();
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            REQUIRE(view[i, j] == (i == j ? 1.0 : 0.0));
        }
    }
}

TEST_CASE("tensor::from_array", "[Math][TensorCreation]") {
    SECTION("1D tensor") {
        auto t = tn::from_array<double>({5}, {1.0, 2.0, 3.0, 4.0, 5.0});
        REQUIRE(t.shape()[0] == 5);
        REQUIRE(t.data()[0] == 1.0);
        REQUIRE(t.data()[4] == 5.0);
    }

    SECTION("2D tensor") {
        auto t = tn::from_array<double>(
            {2, 3},
            {1.0, 2.0, 3.0,
             4.0, 5.0, 6.0}
        );
        REQUIRE(t.shape()[0] == 2);
        REQUIRE(t.shape()[1] == 3);

        auto view = t.view<2>();
        REQUIRE(view[0, 0] == 1.0);
        REQUIRE(view[0, 2] == 3.0);
        REQUIRE(view[1, 0] == 4.0);
        REQUIRE(view[1, 2] == 6.0);
    }

    SECTION("3D tensor") {
        auto t = tn::from_array<int>(
            {2, 2, 2},
            {1, 2, 3, 4, 5, 6, 7, 8}
        );
        REQUIRE(t.ndim() == 3);
        REQUIRE(t.size() == 8);

        auto view = t.view<3>();
        REQUIRE(view[0, 0, 0] == 1);
        REQUIRE(view[1, 1, 1] == 8);
    }

    SECTION("Size mismatch throws") {
        REQUIRE_THROWS_AS(
            tn::from_array<double>({2, 3}, {1.0, 2.0, 3.0}),  // Only 3 elements, need 6
            std::invalid_argument
        );
    }
}

TEST_CASE("Combined usage example", "[Math][TensorCreation][Example]") {
    // Create probability distribution matrix
    auto probs = tn::zeros({3, 4});
    REQUIRE(probs.size() == 12);

    // Create indices
    auto v_idx = tn::arange<size_t>(50);
    REQUIRE(v_idx.shape()[0] == 50);

    // Create test matrix
    auto A = tn::from_array<double>(
        {2, 2},
        {1.0, 2.0,
         3.0, 4.0}
    );

    // Verify we can access it
    auto view = A.view<2>();
    REQUIRE(view[0, 1] == 2.0);
    REQUIRE(view[1, 0] == 3.0);
}

TEST_CASE("tensor::empty", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create empty tensor (uninitialized)
    auto tensor = tn::empty<double>({3, 4});

    REQUIRE(tensor.ndim() == 2);
    REQUIRE(tensor.shape()[0] == 3);
    REQUIRE(tensor.shape()[1] == 4);
    REQUIRE(tensor.size() == 12);

    // Test that we can write to it (won't segfault)
    auto view = tensor.view<2>();
    view[0, 0] = 42.0;
    REQUIRE(view[0, 0] == 42.0);

    // Test 1D
    auto vec = tn::empty<int>({100});
    REQUIRE(vec.ndim() == 1);
    REQUIRE(vec.size() == 100);
}

TEST_CASE("tensor::zeros_like", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create original tensor
    auto original = tn::ones<double>({3, 4});

    // Create zeros_like
    auto zeros = tn::zeros_like(original);

    test::require_shape(zeros, {3, 4});
    test::require_all_equal(zeros, 0.0);

    // Test 3D
    auto tensor3d = tn::full<float>({2, 3, 4}, 99.0f);
    auto zeros3d = tn::zeros_like(tensor3d);
    test::require_shape(zeros3d, {2, 3, 4});
    test::require_all_equal(zeros3d, 0.0f);
}

TEST_CASE("tensor::ones_like", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create original tensor
    auto original = tn::zeros<double>({2, 5});

    // Create ones_like
    auto ones = tn::ones_like(original);

    test::require_shape(ones, {2, 5});
    test::require_all_equal(ones, 1.0);
}

TEST_CASE("tensor::full_like", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create original tensor
    auto original = tn::zeros<int>({4, 3});

    // Create full_like with value 42
    auto filled = tn::full_like(original, 42);

    test::require_shape(filled, {4, 3});
    test::require_all_equal(filled, 42);

    // Test with negative value
    auto neg = tn::full_like(original, -7);
    REQUIRE(neg.view<2>()[0, 0] == -7);
}

TEST_CASE("tensor::diag - create from vector", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create diagonal matrix from values
    auto mat = tn::diag<double>({1.0, 2.0, 3.0});

    test::require_shape(mat, {3, 3});

    auto view = mat.view<2>();

    // Check diagonal
    REQUIRE(view[0, 0] == 1.0);
    REQUIRE(view[1, 1] == 2.0);
    REQUIRE(view[2, 2] == 3.0);

    // Check off-diagonal zeros
    REQUIRE(view[0, 1] == 0.0);
    REQUIRE(view[0, 2] == 0.0);
    REQUIRE(view[1, 0] == 0.0);
    REQUIRE(view[1, 2] == 0.0);
    REQUIRE(view[2, 0] == 0.0);
    REQUIRE(view[2, 1] == 0.0);

    // Test with single element
    auto single = tn::diag<double>({5.0});
    test::require_shape(single, {1, 1});
    REQUIRE(single.view<2>()[0, 0] == 5.0);
}

TEST_CASE("tensor::diag - extract from matrix", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create a matrix
    auto mat = tn::from_array<double>(
        {3, 3},
        {1.0, 2.0, 3.0,
         4.0, 5.0, 6.0,
         7.0, 8.0, 9.0}
    );

    // Extract diagonal
    auto diag_vec = tn::diag(mat);

    test::require_shape(diag_vec, {3});

    auto view = diag_vec.view<1>();
    REQUIRE(view[0] == 1.0);
    REQUIRE(view[1] == 5.0);
    REQUIRE(view[2] == 9.0);
    REQUIRE(view[1] == 5.0);
    REQUIRE(view[2] == 9.0);

    // Test with non-square matrix (take minimum dimension)
    auto rect = tn::from_array<double>(
        {2, 4},
        {1.0, 2.0, 3.0, 4.0,
         5.0, 6.0, 7.0, 8.0}
    );

    auto rect_diag = tn::diag(rect);
    REQUIRE(rect_diag.shape()[0] == 2);
    REQUIRE(rect_diag.view<1>()[0] == 1.0);
    REQUIRE(rect_diag.view<1>()[1] == 6.0);

    // Test error handling for non-2D tensor
    auto tensor3d = tn::zeros<double>({2, 3, 4});
    REQUIRE_THROWS_AS(tn::diag(tensor3d), std::invalid_argument);
}

TEST_CASE("tensor::diag - round trip", "[Math][Tensor][Creation]") {
    namespace tn = igor::math::tensor;

    // Create diagonal matrix, extract diagonal, recreate matrix
    auto original = tn::diag<double>({10.0, 20.0, 30.0});
    auto extracted = tn::diag(original);
    auto recreated = tn::diag<double>({10.0, 20.0, 30.0});

    REQUIRE(extracted.shape()[0] == 3);
    REQUIRE(extracted.view<1>()[0] == 10.0);
    REQUIRE(extracted.view<1>()[1] == 20.0);
    REQUIRE(extracted.view<1>()[2] == 30.0);

    // Verify recreated matches original
    auto orig_view = original.view<2>();
    auto rec_view = recreated.view<2>();
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            REQUIRE(orig_view[i, j] == rec_view[i, j]);
        }
    }
}

