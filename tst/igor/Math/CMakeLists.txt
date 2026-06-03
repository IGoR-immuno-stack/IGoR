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

#igor_copy_runtime_dlls(math_tests)
include(Catch)
#set(_math_tests_dl_paths "$<TARGET_FILE_DIR:math_tests>")

catch_discover_tests(math_tests
    #DISCOVERY_MODE PRE_TEST
    #DL_PATHS ${_math_tests_dl_paths}
)

######################################################################
### CMakeLists.txt ends here
