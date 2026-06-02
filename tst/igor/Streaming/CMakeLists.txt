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
catch_discover_tests(streaming_tests
    DISCOVERY_MODE PRE_TEST
)

######################################################################
### CMakeLists.txt ends here
