enable_testing()
# CMake configuration for tests

# Find Catch2
find_package(Catch2 3 REQUIRED)

# Add test sources (existing tests)
file(GLOB TEST_SOURCES "*.cpp")

if(TEST_SOURCES)
    add_executable(igor_tests ${TEST_SOURCES})
    target_link_libraries(igor_tests PRIVATE igor::Core Catch2::Catch2WithMain)
    include(Catch)
    catch_discover_tests(igor_tests)
endif()

# Add streaming tests
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/igor/Streaming")
    set(STREAMING_TEST_SOURCES
        igor/Streaming/StreamingTestUtils.cpp
        igor/Streaming/test_streaming_helpers.cpp
        igor/Streaming/test_parquet_reader.cpp
        igor/Streaming/test_parquet_writer.cpp
        igor/Streaming/test_airr_rearrangement_reader.cpp
        igor/Streaming/test_airr_rearrangement_writer.cpp
    )

    if(STREAMING_TEST_SOURCES)
        add_executable(streaming_tests ${STREAMING_TEST_SOURCES})

        target_link_libraries(streaming_tests
            PRIVATE
                igor::Core
                igor::Streaming
                Catch2::Catch2WithMain
                sparrow::sparrow
                Arrow::arrow_shared
                Parquet::parquet_shared)

        # Configure test data path at compile time
        target_compile_definitions(streaming_tests PRIVATE
            IGOR_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/igor/Streaming/test_data")

        # Register with CTest
        include(Catch)
        catch_discover_tests(streaming_tests)
    endif()
endif()
