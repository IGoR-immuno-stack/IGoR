### CMakeLists.txt ---

set(MODEL_TEST_SOURCES
    test_EventFactory.cpp
    test_InferenceEngine.cpp
    test_InferenceHandlers.cpp
    test_LegacyBridge.cpp
    test_ModelInference.cpp
    test_NormalizationCompareLegacy.cpp
    test_RecombinationModel.cpp
    test_SamplingEngine.cpp
    test_SamplingHandlers.cpp
    test_Topology.cpp
)

add_executable(model_tests ${MODEL_TEST_SOURCES})

set_target_properties(model_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(model_tests
    PRIVATE
        igor::Model
        igor::Core
        Catch2::Catch2WithMain
)

target_compile_definitions(model_tests
    PRIVATE
        IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models"
)

include(Catch)
# Same labelling as tst/igor/Core: CI's main job runs `ctest -L unit`, so an unlabelled
# discovery is a suite CI never runs. The two [inference] cases drive a full EM through
# GenModel on generated data (about 11 s together) and go with [integration].
catch_discover_tests(model_tests
    TEST_SPEC "~[inference]~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(model_tests
    TEST_SPEC "[inference],[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

######################################################################
### CMakeLists.txt ends here
