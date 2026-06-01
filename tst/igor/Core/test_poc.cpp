#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cstdio>

#include "igor_test_helpers.h"

TEST_CASE("TEST_CASE_1", "[unit][regression]") {
    auto mode = GENERATE_ALL();
    announce_mode(mode);

    std::printf("MODE_CASE_1: %s\n", to_string(mode));

    SECTION("1_unit") {
        if (mode != TestMode::Unit) {
            return;
        }
        std::printf("== Effectively running TEST_CASE_1 / 1_unit\n");
        SUCCEED();
    }

    SECTION("1_regression") {
        if (mode != TestMode::Regression) {
            return;
        }
        std::printf("== Effectively running TEST_CASE_1 / 1_regression\n");
        SUCCEED();
    }
}

TEST_CASE("TEST_CASE_2", "[unit][convergence]") {
    auto mode = GENERATE_ALL();
    announce_mode(mode);

    std::printf("MODE_CASE_2: %s\n", to_string(mode));

    SECTION("2_unit") {
        if (mode != TestMode::Unit) {
            return;
        }
        std::printf("== Effectively running TEST_CASE_2 / 2_unit\n");
        SUCCEED();
    }

    SECTION("2_convergence") {
        if (mode != TestMode::Convergence) {
            return;
        }
        std::printf("== Effectively running TEST_CASE_2 / 2_convergence\n");
        SUCCEED();
    }
}
