#pragma once

#include <catch2/generators/catch_generators.hpp>

enum class TestMode {
    Unit,
    Integration,
    Convergence,
    Regression,
    Benchmark
};

inline const char* to_string(TestMode mode) {
    switch (mode) {
    case TestMode::Unit: return "unit";
    case TestMode::Integration: return "integration";
    case TestMode::Convergence: return "convergence";
    case TestMode::Regression: return "regression";
    case TestMode::Benchmark: return "benchmark";
    }
    return "unknown";
}

inline void announce_mode(TestMode mode) {
    std::printf("Launching tests with '%s' option\n", to_string(mode));
}

#define GENERATE_ALL() \
    GENERATE(TestMode::Unit, TestMode::Integration, TestMode::Convergence, TestMode::Regression, TestMode::Benchmark)
