### CMakeLists.txt for ModelIO Tests ---

add_executable(modelio_tests)

target_sources(modelio_tests
  PRIVATE
    test_model_structure.cpp
)

target_link_libraries(modelio_tests
  PRIVATE
    igor::ModelIO
    igor::Core
    Catch2::Catch2WithMain
)

# Register tests with CTest
catch_discover_tests(modelio_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  TEST_PREFIX "ModelIO::"
)

######################################################################
### CMakeLists.txt ends here
