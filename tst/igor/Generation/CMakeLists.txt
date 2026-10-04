### CMakeLists.txt ---

add_executable(generation_tests
    test_SamplingEngine.cpp
    test_SamplingHandlers.cpp
)

# Legacy tests: the tests of the code still under <Layer>/Legacy.
target_sources(generation_tests
    PRIVATE
        Legacy/test_fast_sampling.cpp
        Legacy/test_generation.cpp
)

set_target_properties(generation_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(generation_tests
    PRIVATE
        igor::Generation
        Catch2::Catch2WithMain
)

target_compile_definitions(generation_tests
    PRIVATE
        IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
)

include(Catch)
catch_discover_tests(generation_tests
    TEST_SPEC "~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(generation_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

catch_discover_tests(generation_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
