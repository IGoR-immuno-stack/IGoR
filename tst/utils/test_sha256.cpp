/*
 * test_sha256.cpp
 *
 *  Test suite for sha256_openssl.h — verifies igor_cli::sha256_file()
 *  against known SHA-256 reference values (NIST FIPS 180-4 vectors).
 */

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "sha256_openssl.h"

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

static fs::path write_tmp(const std::string &name, const std::string &content)
{
    fs::path p = fs::temp_directory_path() / name;
    std::ofstream f(p, std::ios::binary);
    f << content;
    return p;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("sha256_file: empty file", "[sha256][utils]")
{
    const fs::path p = write_tmp("igor_sha256_empty.bin", "");
    CHECK(igor_cli::sha256_file(p)
          == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    fs::remove(p);
}

TEST_CASE("sha256_file: known string 'abc'", "[sha256][utils]")
{
    // NIST FIPS 180-4 reference vector for "abc"
    const fs::path p = write_tmp("igor_sha256_abc.bin", "abc");
    //
    // echo -n 'abc' | sha256  -> ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
    //
    CHECK(igor_cli::sha256_file(p)
          == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    fs::remove(p);
}

TEST_CASE("sha256_file: known string 'hello world'", "[sha256][utils]")
{
    const fs::path p = write_tmp("igor_sha256_hello.bin", "hello world");
    //
    // echo -n 'hello world' | sha256 -> b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9
    //
    CHECK(igor_cli::sha256_file(p)
          == "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9");
    fs::remove(p);
}

TEST_CASE("sha256_file: binary content (256 bytes 0x00..0xff)", "[sha256][utils]")
{
    std::string bytes(256, '\0');
    for (int i = 0; i < 256; ++i) {
        bytes[i] = static_cast<char>(i);
    }
    const fs::path p = write_tmp("igor_sha256_binary.bin", bytes);
    CHECK(igor_cli::sha256_file(p)
          == "40aff2e9d2d8922e47afd4648e6967497158785fbd1da870e7110266bf944880");
    fs::remove(p);
}

TEST_CASE("sha256_file: missing file throws", "[sha256][utils]")
{
    const fs::path p = fs::temp_directory_path() / "igor_sha256_does_not_exist_xyz.bin";
    fs::remove(p);
    CHECK_THROWS_AS(igor_cli::sha256_file(p), std::runtime_error);
}

TEST_CASE("Sha256: streaming update matches file hash", "[sha256][utils]")
{
    const std::string content = "IGoR streaming SHA-256 test";
    const fs::path p = write_tmp("igor_sha256_stream.bin", content);
    const std::string expected = igor_cli::sha256_file(p);
    fs::remove(p);

    // hash via update() un octet à la fois
    igor_cli::Sha256 sha;
    for (unsigned char c : content) {
        sha.update(&c, 1);
    }
    CHECK(sha.final() == expected);
}
