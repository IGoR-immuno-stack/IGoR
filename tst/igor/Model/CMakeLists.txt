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
catch_discover_tests(model_tests)

######################################################################
### CMakeLists.txt ends here
