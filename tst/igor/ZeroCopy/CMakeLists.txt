### CMakeLists.txt for ZeroCopy Tests ---

add_executable(zerocopy_tests
  test_sequence_view.cpp
  test_score_view.cpp
  test_model_marginals.cpp
)

target_link_libraries(zerocopy_tests
  PRIVATE
    igor::ZeroCopy
    Catch2::Catch2WithMain
    sparrow::sparrow
)

# Register tests
include(Catch)
catch_discover_tests(zerocopy_tests)
