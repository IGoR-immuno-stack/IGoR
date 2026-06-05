// sha256_openssl.h
// Drop-in replacement for sha256.h using OpenSSL EVP API (OpenSSL >= 1.1, 3.x compatible).
// Same interface: igor_cli::Sha256 (update/final) + igor_cli::sha256_file().

#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <openssl/evp.h>

namespace igor_cli {

class Sha256 {
public:
    Sha256() : ctx_(EVP_MD_CTX_new())
    {
        if (!ctx_) {
            throw std::runtime_error("Sha256: EVP_MD_CTX_new failed");
        }
        if (EVP_DigestInit_ex(ctx_.get(), EVP_sha256(), nullptr) != 1) {
            throw std::runtime_error("Sha256: EVP_DigestInit_ex failed");
        }
    }

    void update(const unsigned char *data, std::size_t len)
    {
        if (EVP_DigestUpdate(ctx_.get(), data, len) != 1) {
            throw std::runtime_error("Sha256: EVP_DigestUpdate failed");
        }
    }

    std::string final()
    {
        unsigned char hash[EVP_MAX_MD_SIZE];
        unsigned int hash_len = 0;
        if (EVP_DigestFinal_ex(ctx_.get(), hash, &hash_len) != 1) {
            throw std::runtime_error("Sha256: EVP_DigestFinal_ex failed");
        }

        std::ostringstream out;
        for (unsigned int i = 0; i < hash_len; ++i) {
            out << std::hex << std::setfill('0') << std::setw(2)
                << static_cast<int>(hash[i]);
        }
        return out.str();
    }

private:
    struct EvpCtxDeleter {
        void operator()(EVP_MD_CTX *p) const { EVP_MD_CTX_free(p); }
    };
    std::unique_ptr<EVP_MD_CTX, EvpCtxDeleter> ctx_;
};

inline std::string sha256_file(const std::filesystem::path &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open input for hashing: " + path.string());
    }
    Sha256 sha;
    std::array<unsigned char, 8192> buffer {};
    while (input) {
        input.read(reinterpret_cast<char *>(buffer.data()),
                   static_cast<std::streamsize>(buffer.size()));
        if (input.gcount() > 0) {
            sha.update(buffer.data(), static_cast<std::size_t>(input.gcount()));
        }
    }
    return sha.final();
}

} // namespace igor_cli
