### CMakeLists.txt ---

file(GLOB MODEL_TEST_SOURCES "*.cpp")

if(MODEL_TEST_SOURCES)
    add_executable(model_tests ${MODEL_TEST_SOURCES})
    set_target_properties(model_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(model_tests PRIVATE igor::Model igor::Core Catch2::Catch2WithMain)
    target_compile_definitions(model_tests PRIVATE IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models")
    igor_copy_runtime_dlls(model_tests)
    include(Catch)
    catch_discover_tests(model_tests
        DL_PATHS $<TARGET_FILE_DIR:model_tests> $ENV{CONDA_PREFIX}/Library/bin
        DISCOVERY_MODE PRE_TEST
    )
endif()

######################################################################
### CMakeLists.txt ends here
