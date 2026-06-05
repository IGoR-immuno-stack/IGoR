### CMakeLists.txt ---

find_package(OpenSSL REQUIRED)

add_executable(igor_utils_tests
  test_sha256.cpp
)

target_include_directories(igor_utils_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/app/igor)  # find sha256_openssl.h

target_link_libraries(igor_utils_tests
    PRIVATE
        OpenSSL::Crypto
        Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(igor_utils_tests
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 30
        LABELS "unit"
)

######################################################################
### CMakeLists.txt ends here
