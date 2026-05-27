### CMakeLists.txt ---

set(CORE_TEST_SOURCES
    test_Contexts.cpp
    test_EventUtils.cpp
    test_fast_sampling.cpp
    test_generation.cpp
    test_inference.cpp
    test_sanity.cpp
    test_utils.cpp
    test_gene_choice_iterate.cpp
)

add_executable(igor_tests ${CORE_TEST_SOURCES})

target_link_libraries(igor_tests
    PRIVATE
        igor::Core
        Catch2::Catch2WithMain)

target_compile_definitions(igor_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

include(Catch)
catch_discover_tests(igor_tests
    TEST_SPEC "~[integration]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(igor_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

######################################################################
### CMakeLists.txt ends here
