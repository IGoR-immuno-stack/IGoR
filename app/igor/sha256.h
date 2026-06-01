#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace igor_cli {

class Sha256 {
public:
    void update(const unsigned char *data, std::size_t len)
    {
        for (std::size_t i = 0; i < len; ++i) {
            data_[datalen_++] = data[i];
            if (datalen_ == 64) {
                transform();
                bitlen_ += 512;
                datalen_ = 0;
            }
        }
    }

    std::string final()
    {
        std::size_t i = datalen_;
        data_[i++] = 0x80;
        if (i > 56) {
            while (i < 64) {
                data_[i++] = 0x00;
            }
            transform();
            i = 0;
        }
        while (i < 56) {
            data_[i++] = 0x00;
        }

        bitlen_ += datalen_ * 8;
        for (int j = 7; j >= 0; --j) {
            data_[63 - j] = static_cast<unsigned char>((bitlen_ >> (j * 8)) & 0xff);
        }
        transform();

        std::ostringstream out;
        for (auto word : state_) {
            out << std::hex << std::setfill('0') << std::setw(8) << word;
        }
        return out.str();
    }

private:
    static constexpr std::array<uint32_t, 64> k_ = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
            0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
            0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
            0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
            0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
            0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
            0xc67178f2};

    static uint32_t rotr(uint32_t x, uint32_t n) { return (x >> n) | (x << (32 - n)); }
    static uint32_t choose(uint32_t e, uint32_t f, uint32_t g) { return (e & f) ^ (~e & g); }
    static uint32_t majority(uint32_t a, uint32_t b, uint32_t c) { return (a & b) ^ (a & c) ^ (b & c); }
    static uint32_t sig0(uint32_t x) { return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3); }
    static uint32_t sig1(uint32_t x) { return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10); }
    static uint32_t ep0(uint32_t x) { return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22); }
    static uint32_t ep1(uint32_t x) { return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25); }

    void transform()
    {
        uint32_t m[64];
        for (int i = 0, j = 0; i < 16; ++i, j += 4) {
            m[i] = (static_cast<uint32_t>(data_[j]) << 24) | (static_cast<uint32_t>(data_[j + 1]) << 16)
                    | (static_cast<uint32_t>(data_[j + 2]) << 8) | static_cast<uint32_t>(data_[j + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            m[i] = sig1(m[i - 2]) + m[i - 7] + sig0(m[i - 15]) + m[i - 16];
        }

        uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int i = 0; i < 64; ++i) {
            uint32_t t1 = h + ep1(e) + choose(e, f, g) + k_[i] + m[i];
            uint32_t t2 = ep0(a) + majority(a, b, c);
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += h;
    }

    std::array<unsigned char, 64> data_ {};
    uint32_t datalen_ = 0;
    uint64_t bitlen_ = 0;
    std::array<uint32_t, 8> state_ {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
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
        input.read(reinterpret_cast<char *>(buffer.data()), static_cast<std::streamsize>(buffer.size()));
        if (input.gcount() > 0) {
            sha.update(buffer.data(), static_cast<std::size_t>(input.gcount()));
        }
    }
    return sha.final();
}

} // namespace igor_cli
