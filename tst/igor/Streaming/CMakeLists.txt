### CMakeLists.txt ---

set(STREAMING_TEST_SOURCES
    StreamingTestUtils.cpp
    test_streaming_helpers.cpp
    test_parquet_reader.cpp
    test_parquet_writer.cpp
    test_airr_rearrangement_reader.cpp
    test_airr_rearrangement_writer.cpp
    test_airr_alignment.cpp
)

add_executable(streaming_tests ${STREAMING_TEST_SOURCES})

target_link_libraries(streaming_tests
    PRIVATE
        igor::Core
        igor::Streaming
        Catch2::Catch2WithMain
        sparrow::sparrow
        Arrow::arrow_shared
        Parquet::parquet_shared)

target_compile_definitions(streaming_tests PRIVATE
    IGOR_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data")

include(Catch)
# Labelled like tst/igor/Core so that CI's `ctest -L unit` runs this suite. Self-contained:
# temporary files written by the tests, plus the fixtures under IGOR_TEST_DATA_DIR.
catch_discover_tests(streaming_tests
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

######################################################################
### CMakeLists.txt ends here
