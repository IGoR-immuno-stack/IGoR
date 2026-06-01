### CMakeLists.txt ---

file(GLOB MATH_TEST_SOURCES "*.cpp")

if(MATH_TEST_SOURCES)
    add_executable(math_tests ${MATH_TEST_SOURCES})
    set_target_properties(math_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(math_tests PRIVATE igor::Math Catch2::Catch2WithMain)
    igor_copy_runtime_dlls(math_tests)
    include(Catch)
    catch_discover_tests(math_tests
        DL_PATHS $<TARGET_FILE_DIR:math_tests> $ENV{CONDA_PREFIX}/Library/bin
        DISCOVERY_MODE PRE_TEST
    )
endif()

######################################################################
### CMakeLists.txt ends here
