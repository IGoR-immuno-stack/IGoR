/*
 * test_InferenceHandlers.cpp
 *
 *  Created on: Feb 10, 2026
 *  Updated on: Feb 21, 2026
 *
 *  Unit tests for CategoricalInferenceHandler and MarkovInferenceHandler.
 *  Uses TEMPLATE_TEST_CASE to validate both double and long double.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <igor/Model/CategoricalInferenceHandler.h>
#include <igor/Model/MarkovInferenceHandler.h>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ─── Helpers ───────────────────────────────────────────────────────────

/// Create a 1D uniform tensor of the given size.
template <typename T>
static igor::math::Tensor<T> uniform_1d(std::size_t n) {
    igor::math::Tensor<T> t({n});
    T val = T(1) / static_cast<T>(n);
    std::fill(t.begin(), t.end(), val);
    return t;
}

/// Create a 2D tensor [rows x cols] with uniform rows (each row sums to 1).
template <typename T>
static igor::math::Tensor<T> uniform_2d(std::size_t rows, std::size_t cols) {
    igor::math::Tensor<T> t({rows, cols});
    T val = T(1) / static_cast<T>(cols);
    std::fill(t.begin(), t.end(), val);
    return t;
}

/// Create a 2D tensor [rows x cols] with uniform rows (for Markov transition matrix).
template <typename T>
static igor::math::Tensor<T> uniform_markov(std::size_t n_states) {
    return uniform_2d<T>(n_states, n_states);
}

// ─── CategoricalInferenceHandler Tests (1D — no parents) ──────────────

TEMPLATE_TEST_CASE("CategoricalInferenceHandler 1D construction",
                   "[model][handler]", double, long double) {
    auto weights = uniform_1d<TestType>(5);
    CategoricalInferenceHandler<TestType> handler("v_gene", 0, weights);

    SECTION("name is set") {
        REQUIRE(handler.name() == "v_gene");
    }
    SECTION("uid is set") {
        REQUIRE(handler.uid() == 0);
    }
    SECTION("realizationCount matches tensor dim 0") {
        REQUIRE(handler.realizationCount() == 5);
    }
    SECTION("weights() returns the borrowed tensor") {
        REQUIRE(&handler.weights() == &weights);
        REQUIRE(handler.weights().size() == 5);
    }
    SECTION("accumulator has same shape as weights") {
        REQUIRE(handler.accumulator().shape() == handler.weights().shape());
    }
    SECTION("accumulator starts at zero") {
        for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
            REQUIRE(handler.accumulator().data()[i] == TestType(0));
        }
    }
}

TEMPLATE_TEST_CASE("CategoricalInferenceHandler 1D maximizeLikelihood",
                   "[model][handler]", double, long double) {
    auto weights = uniform_1d<TestType>(4);
    CategoricalInferenceHandler<TestType> handler("j_gene", 1, weights);

    // Simulate accumulation: [1, 2, 3, 4] → normalise to [0.1, 0.2, 0.3, 0.4]
    handler.accumulator().data()[0] = TestType(1);
    handler.accumulator().data()[1] = TestType(2);
    handler.accumulator().data()[2] = TestType(3);
    handler.accumulator().data()[3] = TestType(4);

    handler.maximizeLikelihood();

    REQUIRE_THAT(double(handler.weights().data()[0]), WithinAbs(0.1, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[1]), WithinAbs(0.2, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[2]), WithinAbs(0.3, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[3]), WithinAbs(0.4, 1e-10));
}

TEMPLATE_TEST_CASE("CategoricalInferenceHandler M-step writes back into original tensor",
                   "[model][handler]", double, long double) {
    auto weights = uniform_1d<TestType>(3);
    CategoricalInferenceHandler<TestType> handler("d_gene", 2, weights);

    // Accumulate [3, 6, 9]
    handler.accumulator().data()[0] = TestType(3);
    handler.accumulator().data()[1] = TestType(6);
    handler.accumulator().data()[2] = TestType(9);

    handler.maximizeLikelihood();

    // The *original* tensor should have been updated
    REQUIRE_THAT(double(weights.data()[0]), WithinAbs(1.0 / 6.0, 1e-10));
    REQUIRE_THAT(double(weights.data()[1]), WithinAbs(2.0 / 6.0, 1e-10));
    REQUIRE_THAT(double(weights.data()[2]), WithinAbs(3.0 / 6.0, 1e-10));
}

TEMPLATE_TEST_CASE("CategoricalInferenceHandler resetAccumulator zeros out",
                   "[model][handler]", double, long double) {
    auto weights = uniform_1d<TestType>(3);
    CategoricalInferenceHandler<TestType> handler("v_gene", 0, weights);

    handler.accumulator().data()[0] = TestType(42);
    handler.accumulator().data()[1] = TestType(99);
    handler.resetAccumulator();

    for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
        REQUIRE(handler.accumulator().data()[i] == TestType(0));
    }
}

// ─── CategoricalInferenceHandler Tests (2D — with parents) ────────────

TEMPLATE_TEST_CASE("CategoricalInferenceHandler 2D construction",
                   "[model][handler]", double, long double) {
    // Shape: [3 realizations, 2 parent values]
    auto weights = uniform_2d<TestType>(3, 2);
    CategoricalInferenceHandler<TestType> handler("v_3_del", 3, weights);

    SECTION("shape is {3, 2}") {
        auto s = handler.shape();
        REQUIRE(s.size() == 2);
        REQUIRE(s[0] == 3);
        REQUIRE(s[1] == 2);
    }
    SECTION("realizationCount is first dim") {
        REQUIRE(handler.realizationCount() == 3);
    }
}

// ─── MarkovInferenceHandler Tests ─────────────────────────────────────

TEMPLATE_TEST_CASE("MarkovInferenceHandler construction",
                   "[model][handler]", double, long double) {
    auto weights = uniform_markov<TestType>(4);
    MarkovInferenceHandler<TestType> handler("vd_dinucl", 5, weights);

    SECTION("name is set") {
        REQUIRE(handler.name() == "vd_dinucl");
    }
    SECTION("uid is set") {
        REQUIRE(handler.uid() == 5);
    }
    SECTION("stateCount matches shape[0]") {
        REQUIRE(handler.stateCount() == 4);
    }
    SECTION("weights() returns the borrowed tensor") {
        REQUIRE(&handler.weights() == &weights);
    }
    SECTION("accumulator has same shape as weights") {
        REQUIRE(handler.accumulator().shape() == handler.weights().shape());
    }
    SECTION("accumulator starts at zero") {
        for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
            REQUIRE(handler.accumulator().data()[i] == TestType(0));
        }
    }
}

TEMPLATE_TEST_CASE("MarkovInferenceHandler maximizeLikelihood normalises rows",
                   "[model][handler]", double, long double) {
    // 2 states
    auto weights = uniform_markov<TestType>(2);
    MarkovInferenceHandler<TestType> handler("vd_dinucl", 5, weights);

    // Accumulate: row 0 = [1, 3], row 1 = [2, 2]
    handler.accumulator().data()[0] = TestType(1);
    handler.accumulator().data()[1] = TestType(3);
    handler.accumulator().data()[2] = TestType(2);
    handler.accumulator().data()[3] = TestType(2);

    handler.maximizeLikelihood();

    REQUIRE_THAT(double(handler.weights().data()[0]), WithinAbs(0.25, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[1]), WithinAbs(0.75, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[2]), WithinAbs(0.50, 1e-10));
    REQUIRE_THAT(double(handler.weights().data()[3]), WithinAbs(0.50, 1e-10));
}

TEMPLATE_TEST_CASE("MarkovInferenceHandler M-step writes back into original tensor",
                   "[model][handler]", double, long double) {
    auto weights = uniform_markov<TestType>(2);
    MarkovInferenceHandler<TestType> handler("vd_dinucl", 5, weights);

    handler.accumulator().data()[0] = TestType(1);
    handler.accumulator().data()[1] = TestType(3);
    handler.accumulator().data()[2] = TestType(2);
    handler.accumulator().data()[3] = TestType(2);

    handler.maximizeLikelihood();

    // The *original* tensor should have been updated
    REQUIRE_THAT(double(weights.data()[0]), WithinAbs(0.25, 1e-10));
    REQUIRE_THAT(double(weights.data()[1]), WithinAbs(0.75, 1e-10));
    REQUIRE_THAT(double(weights.data()[2]), WithinAbs(0.50, 1e-10));
    REQUIRE_THAT(double(weights.data()[3]), WithinAbs(0.50, 1e-10));
}

TEMPLATE_TEST_CASE("MarkovInferenceHandler resetAccumulator zeros out",
                   "[model][handler]", double, long double) {
    auto weights = uniform_markov<TestType>(4);
    MarkovInferenceHandler<TestType> handler("vd_dinucl", 5, weights);

    handler.accumulator().data()[0] = TestType(42);
    handler.accumulator().data()[3] = TestType(99);
    handler.resetAccumulator();

    for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
        REQUIRE(handler.accumulator().data()[i] == TestType(0));
    }
}

// ─── Polymorphism ─────────────────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceHandler polymorphism",
                   "[model][handler]", double, long double) {
    auto cat_weights = uniform_1d<TestType>(3);
    auto markov_weights = uniform_markov<TestType>(4);

    auto cat = std::make_unique<CategoricalInferenceHandler<TestType>>("v_gene", 0, cat_weights);
    auto markov = std::make_unique<MarkovInferenceHandler<TestType>>("vd_dinucl", 1, markov_weights);

    std::vector<std::unique_ptr<InferenceHandler<TestType>>> handlers;
    handlers.push_back(std::move(cat));
    handlers.push_back(std::move(markov));

    SECTION("all are InferenceHandler") {
        for (const auto& h : handlers) {
            REQUIRE(!h->name().empty());
            REQUIRE(h->weights().size() > 0);
            REQUIRE(h->accumulator().size() > 0);
        }
    }
    SECTION("resetAccumulator works polymorphically") {
        for (auto& h : handlers) {
            h->accumulator().data()[0] = TestType(123);
            h->resetAccumulator();
            REQUIRE(h->accumulator().data()[0] == TestType(0));
        }
    }
}
