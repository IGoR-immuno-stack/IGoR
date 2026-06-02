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
target_link_libraries(math_tests PRIVATE igor::Math Catch2::Catch2WithMain)
include(Catch)
catch_discover_tests(math_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES
        ENVIRONMENT "PATH=${PIXI_BIN_DIR};$ENV{PATH}"
)

######################################################################
### CMakeLists.txt ends here
