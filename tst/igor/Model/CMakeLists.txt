### CMakeLists.txt ---

add_executable(model_tests
    test_EventFactory.cpp
    test_LegacyBridge.cpp
    test_RecombinationModel.cpp
    test_Topology.cpp
    test_bound_tightness.cpp
    Legacy/test_Contexts.cpp
    Legacy/test_deletion_iterate.cpp
    Legacy/test_dinucl_markov_iterate.cpp
    Legacy/test_error_rate_bound.cpp
    Legacy/test_event_capabilities.cpp
    Legacy/test_EventUtils.cpp
    Legacy/test_gene_choice_iterate.cpp
    Legacy/test_genetic_code.cpp
    Legacy/test_insertion_iterate.cpp
    test_junction_geometry.cpp
    Legacy/test_model_format_v2.cpp
    Legacy/test_model_json.cpp
    Legacy/test_proba_bound_benchmark.cpp
    test_safety_matrix.cpp
    test_span_profile.cpp
    Legacy/test_utils.cpp
)

set_target_properties(model_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(model_tests
    PRIVATE
        igor::Model
        Catch2::Catch2WithMain
)

target_compile_definitions(model_tests
    PRIVATE
        IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
)

include(Catch)
# CI's main job runs `ctest -L unit`, so an unlabelled discovery is a suite CI never runs.
catch_discover_tests(model_tests
    TEST_SPEC "~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(model_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

catch_discover_tests(model_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
