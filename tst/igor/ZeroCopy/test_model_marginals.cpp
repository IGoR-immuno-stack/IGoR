/*
 * test_model_marginals.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  Unit tests for ZeroCopy Model_marginals
 */

#include <igor/ZeroCopy/Model_marginals.h>

#include <catch2/catch_test_macros.hpp>
#include <memory>

using namespace igor;

TEST_CASE("Model_marginals provides mdspan view", "[zerocopy][model_marginals]") {
    // Create a Model_marginals object
    Model_marginals marginals;

    // Manually setup internal state for testing
    // Note: marginal_array_smart_p and marginal_arr_size are public (or
    // accessible)
    size_t size = 10;
    marginals.marginal_array_smart_p = std::unique_ptr<long double[]>(new long double[size]);
    // We can't easily set marginal_arr_size because it's private in the original
    // class? Let's check the header again. In the modified header: size_t
    // get_length()const{return marginal_arr_size;}; size_t marginal_arr_size; is
    // private.

    // Set the size manually
    marginals.marginal_arr_size = size;

    // Initialize data
    for (size_t i = 0; i < size; ++i) {
        marginals.marginal_array_smart_p[i] = static_cast<long double>(i) * 1.5;
    }

    // Get the view
    auto view = marginals.get_view();

    // Check view properties
    REQUIRE(view.extent(0) == size);

    // Check data access via view
    REQUIRE(view(0) == 0.0);
    REQUIRE(view(1) == 1.5);
    REQUIRE(view(9) == 13.5);
}
