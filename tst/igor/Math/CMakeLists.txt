### CMakeLists.txt ---

set(MATH_TEST_SOURCES
    test_HybridBuffer.cpp
    test_Linalg_AxisOps.cpp
    test_Linalg_Basic.cpp
    test_Linalg_Broadcast.cpp
    test_Linalg_Probabilistic.cpp
    test_Tensor.cpp
    test_TensorCreation.cpp
    test_TensorOps.cpp
    test_Tensor_Slicing.cpp
)

add_executable(math_tests ${MATH_TEST_SOURCES})

set_target_properties(math_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(math_tests 
    PRIVATE 
        igor::Math 
        Catch2::Catch2WithMain
)

include(Catch)
# Labelled like tst/igor/Core so that CI's `ctest -L unit` runs this suite. Pure, no fixtures.
catch_discover_tests(math_tests
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

######################################################################
### CMakeLists.txt ends here
