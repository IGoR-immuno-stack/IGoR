### CMakeLists.txt ---

add_executable(inference_tests
    test_InferenceEngine.cpp
    test_InferenceHandlers.cpp
    test_ModelInference.cpp
    test_NormalizationCompareLegacy.cpp
)

# Legacy tests: the tests of the code still under <Layer>/Legacy.
target_sources(inference_tests
    PRIVATE
        Legacy/test_inference.cpp
)

set_target_properties(inference_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(inference_tests
    PRIVATE
        igor::Inference
        igor::Generation
        Catch2::Catch2WithMain
)

target_compile_definitions(inference_tests
    PRIVATE
        IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
)

include(Catch)
# The tk [inference] cases drive a full EM through GenModel on generated data and go with
# [integration], except the ex-Core convergence case, which stays in its own group; CI runs `-L unit`.
catch_discover_tests(inference_tests
    TEST_SPEC "~[inference]~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(inference_tests
    TEST_SPEC "[integration],[inference]~[convergence]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

catch_discover_tests(inference_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
