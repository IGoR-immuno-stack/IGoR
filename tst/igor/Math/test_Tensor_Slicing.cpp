#include <catch2/catch_test_macros.hpp>
#include <igor/Math/Tensor.h>
#include <igor/Math/MdspanCompat.h>
#include <iostream>

using namespace igor::math;

TEST_CASE("Tensor Slicing Example with submdspan", "[Math][Tensor][Example]") {
    // 1. Create a 3x4 Tensor
    Tensor<double> tensor({3, 4});
    
    // Fill with data:
    // 0  1  2  3
    // 4  5  6  7
    // 8  9 10 11
    int counter = 0;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            tensor.view<2>()[i, j] = (double)counter++;
        }
    }

#ifdef IGOR_NO_SUBMDSPAN
    WARN("submdspan not available - skipping example");
#else
    // 2. Get a view of the tensor (Rank 2)
    auto full_view = tensor.view<2>();

    // 3. Slice: Get the second row (index 1)
    // resulting view is rank 1, length 4
    // std::full_extent means "keep this dimension"
    auto row_slice = std::submdspan(full_view, 1, std::full_extent);

    // Verify properties of the slice
    REQUIRE(row_slice.rank() == 1);
    REQUIRE(row_slice.extent(0) == 4);
    
    // Verify data correctness (should be 4, 5, 6, 7)
    REQUIRE(row_slice[0] == 4.0);
    REQUIRE(row_slice[1] == 5.0);
    REQUIRE(row_slice[2] == 6.0);
    REQUIRE(row_slice[3] == 7.0);

    // 4. Slice: Get the third column (index 2)
    // resulting view is rank 1, length 3
    auto col_slice = std::submdspan(full_view, std::full_extent, 2);

    REQUIRE(col_slice.rank() == 1);
    REQUIRE(col_slice.extent(0) == 3);
    
    // Verify data (should be 2, 6, 10)
    REQUIRE(col_slice[0] == 2.0);
    REQUIRE(col_slice[1] == 6.0);
    REQUIRE(col_slice[2] == 10.0);

    // 5. Modification via slice affects original tensor
    row_slice[0] = 99.0;
    REQUIRE(tensor.view<2>()[1, 0] == 99.0);
#endif
}


TEST_CASE("Tensor::slice", "[Math][Tensor][Slicing]") {
#ifndef IGOR_NO_SUBMDSPAN
    // Setup 3x4 tensor
    Tensor<double> tensor({3, 4});
    // 0 1 2 3
    // 4 5 6 7
    // 8 9 10 11
    int c = 0;
    for(size_t i=0; i<3; ++i) 
        for(size_t j=0; j<4; ++j) 
            tensor.view<2>()[i, j] = (double)c++;

    // Test slice along dim 0 (row)
    auto row_view = tensor.slice<2>(0, 1); // Dim 0, Index 1 -> Row 1
    REQUIRE(row_view.rank() == 1);
    REQUIRE(row_view.extent(0) == 4);
    REQUIRE(row_view[0] == 4.0);

    // Test slice along dim 1 (col)
    auto col_view = tensor.slice<2>(1, 2); // Dim 1, Index 2 -> Col 2
    REQUIRE(col_view.rank() == 1);
    REQUIRE(col_view.extent(0) == 3);
    REQUIRE(col_view[0] == 2.0);
    REQUIRE(col_view[2] == 10.0);
#endif
}
