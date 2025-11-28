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
    file(GLOB STREAMING_TEST_SOURCES "igor/Streaming/*.cpp")

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

        # Register with CTest
        include(Catch)
        catch_discover_tests(streaming_tests)
    endif()
endif()
