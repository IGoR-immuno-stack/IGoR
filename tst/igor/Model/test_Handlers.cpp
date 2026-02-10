/*
 * test_Handlers.cpp
 *
 *  Created on: Feb 10, 2026
 *
 *  Unit tests for CategoricalHandler and MarkovHandler.
 *  Uses TEMPLATE_TEST_CASE to validate both double and long double.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <igor/Model/CategoricalHandler.h>
#include <igor/Model/MarkovHandler.h>

#include <sstream>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ─── CategoricalHandler Tests (1D — no parents) ───────────────────────

TEMPLATE_TEST_CASE("CategoricalHandler 1D construction", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("v_gene", 5);

    SECTION("name is set") {
        REQUIRE(handler.name() == "v_gene");
    }

    SECTION("n_realizations is set") {
        REQUIRE(handler.n_realizations() == 5);
    }

    SECTION("parameters initialized to uniform") {
        for (std::size_t i = 0; i < 5; ++i) {
            REQUIRE_THAT(static_cast<double>(handler.parameters()(i)),
                         WithinRel(0.2, 1e-10));
        }
    }

    SECTION("accumulator initialized to zero") {
        for (std::size_t i = 0; i < 5; ++i) {
            REQUIRE(handler.accumulator()(i) == TestType(0));
        }
    }

    SECTION("parameter tensor has correct shape") {
        REQUIRE(handler.parameters().ndim() == 1);
        REQUIRE(handler.parameters().shape()[0] == 5);
    }
}

TEMPLATE_TEST_CASE("CategoricalHandler 1D accumulate and normalize", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("d_gene", 3);

    handler.accumulator()(0) += TestType(10);
    handler.accumulator()(1) += TestType(20);
    handler.accumulator()(2) += TestType(30);

    SECTION("accumulator stores values") {
        REQUIRE_THAT(static_cast<double>(handler.accumulator()(0)),
                     WithinRel(10.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.accumulator()(1)),
                     WithinRel(20.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.accumulator()(2)),
                     WithinRel(30.0, 1e-10));
    }

    SECTION("normalize produces correct proportions") {
        handler.maximize_likelihood();

        // Sum to 1
        double total = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            total += static_cast<double>(handler.parameters()(i));
        }
        REQUIRE_THAT(total, WithinRel(1.0, 1e-10));

        // Proportional to accumulator
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0)),
                     WithinRel(1.0 / 6.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(1)),
                     WithinRel(2.0 / 6.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(2)),
                     WithinRel(3.0 / 6.0, 1e-10));
    }
}

TEMPLATE_TEST_CASE("CategoricalHandler zero accumulator stability", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("j_gene", 4);

    TestType initial = handler.parameters()(0);

    // M-step with zero accumulator should not change parameters
    handler.maximize_likelihood();

    for (std::size_t i = 0; i < 4; ++i) {
        REQUIRE(handler.parameters()(i) == initial);
    }
}

TEMPLATE_TEST_CASE("CategoricalHandler reset accumulator", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("v_gene", 3);

    handler.accumulator()(0) += TestType(100);
    handler.accumulator()(1) += TestType(200);
    handler.reset_accumulator();

    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE(handler.accumulator()(i) == TestType(0));
    }
}

TEMPLATE_TEST_CASE("CategoricalHandler I/O round-trip", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> writer("vd_ins", 4);
    writer.accumulator()(0) += TestType(1);
    writer.accumulator()(1) += TestType(2);
    writer.accumulator()(2) += TestType(3);
    writer.accumulator()(3) += TestType(4);
    writer.maximize_likelihood();

    std::stringstream ss;
    writer.write_parameters(ss);

    CategoricalHandler<TestType> reader("vd_ins", 4);
    reader.read_parameters(ss);

    for (std::size_t i = 0; i < 4; ++i) {
        REQUIRE_THAT(static_cast<double>(reader.parameters()(i)),
                     WithinRel(static_cast<double>(writer.parameters()(i)), 1e-10));
    }
}

// ─── CategoricalHandler Tests (2D — with parents) ─────────────────────

TEMPLATE_TEST_CASE("CategoricalHandler 2D construction", "[model][handler]",
                   double, long double) {
    // v_3_del conditioned on v_choice: shape = {n_del, n_v_genes}
    CategoricalHandler<TestType> handler("v_3_del", std::vector<std::size_t>{10, 5});

    SECTION("shape is correct") {
        REQUIRE(handler.shape().size() == 2);
        REQUIRE(handler.shape()[0] == 10);
        REQUIRE(handler.shape()[1] == 5);
    }

    SECTION("n_realizations is first dimension") {
        REQUIRE(handler.n_realizations() == 10);
    }

    SECTION("total size is product of dimensions") {
        REQUIRE(handler.parameters().size() == 50);
    }

    SECTION("parameters tensor is 2D") {
        REQUIRE(handler.parameters().ndim() == 2);
    }

    SECTION("initialized to uniform along axis 0") {
        // Each value should be 1/n_realizations = 1/10
        for (std::size_t i = 0; i < handler.parameters().size(); ++i) {
            REQUIRE_THAT(static_cast<double>(handler.parameters().data()[i]),
                         WithinRel(0.1, 1e-10));
        }
    }
}

TEMPLATE_TEST_CASE("CategoricalHandler 2D normalize per parent slice", "[model][handler]",
                   double, long double) {
    // Shape: {3 realizations, 2 parent values}
    CategoricalHandler<TestType> handler("del", std::vector<std::size_t>{3, 2});

    // Accumulate: for parent=0, counts are [10, 20, 30]; for parent=1, counts are [1, 1, 1]
    auto& acc = handler.accumulator();
    acc(0, 0) = TestType(10);
    acc(1, 0) = TestType(20);
    acc(2, 0) = TestType(30);
    acc(0, 1) = TestType(1);
    acc(1, 1) = TestType(1);
    acc(2, 1) = TestType(1);

    handler.maximize_likelihood();

    SECTION("parent=0 slice normalized") {
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0, 0)),
                     WithinRel(1.0 / 6.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(1, 0)),
                     WithinRel(2.0 / 6.0, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(2, 0)),
                     WithinRel(3.0 / 6.0, 1e-10));

        double slice_sum = 0.0;
        for (std::size_t r = 0; r < 3; ++r) {
            slice_sum += static_cast<double>(handler.parameters()(r, 0));
        }
        REQUIRE_THAT(slice_sum, WithinRel(1.0, 1e-10));
    }

    SECTION("parent=1 slice normalized (uniform)") {
        for (std::size_t r = 0; r < 3; ++r) {
            REQUIRE_THAT(static_cast<double>(handler.parameters()(r, 1)),
                         WithinRel(1.0 / 3.0, 1e-10));
        }
    }
}

// ─── MarkovHandler Tests (2D — no parents) ─────────────────────────────

TEMPLATE_TEST_CASE("MarkovHandler 2D construction", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> handler("vd_dinucl", 4);

    SECTION("name is set") {
        REQUIRE(handler.name() == "vd_dinucl");
    }

    SECTION("n_states is set") {
        REQUIRE(handler.n_states() == 4);
    }

    SECTION("parameters initialized to uniform rows") {
        for (std::size_t from = 0; from < 4; ++from) {
            double row_sum = 0.0;
            for (std::size_t to = 0; to < 4; ++to) {
                REQUIRE_THAT(static_cast<double>(handler.parameters()(from, to)),
                             WithinRel(0.25, 1e-10));
                row_sum += static_cast<double>(handler.parameters()(from, to));
            }
            REQUIRE_THAT(row_sum, WithinRel(1.0, 1e-10));
        }
    }

    SECTION("accumulator initialized to zero") {
        for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
            REQUIRE(handler.accumulator().data()[i] == TestType(0));
        }
    }

    SECTION("parameter tensor has correct shape") {
        REQUIRE(handler.parameters().ndim() == 2);
        REQUIRE(handler.parameters().shape()[0] == 4);
        REQUIRE(handler.parameters().shape()[1] == 4);
    }
}

TEMPLATE_TEST_CASE("MarkovHandler 2D row normalize", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> handler("vd_dinucl", 4);

    auto& acc = handler.accumulator();
    // Row 0: A→A=10, A→C=30
    acc(0, 0) = TestType(10);
    acc(0, 1) = TestType(30);
    // Row 2: equal counts
    acc(2, 0) = TestType(5);
    acc(2, 1) = TestType(5);
    acc(2, 2) = TestType(5);
    acc(2, 3) = TestType(5);

    handler.maximize_likelihood();

    SECTION("row 0 is normalized") {
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0, 0)),
                     WithinRel(0.25, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0, 1)),
                     WithinRel(0.75, 1e-10));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0, 2)),
                     WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(static_cast<double>(handler.parameters()(0, 3)),
                     WithinAbs(0.0, 1e-15));
    }

    SECTION("row 1 unchanged (zero accumulator)") {
        for (std::size_t to = 0; to < 4; ++to) {
            REQUIRE_THAT(static_cast<double>(handler.parameters()(1, to)),
                         WithinRel(0.25, 1e-10));
        }
    }

    SECTION("row 2 is uniform") {
        for (std::size_t to = 0; to < 4; ++to) {
            REQUIRE_THAT(static_cast<double>(handler.parameters()(2, to)),
                         WithinRel(0.25, 1e-10));
        }
    }
}

TEMPLATE_TEST_CASE("MarkovHandler reset accumulator", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> handler("vd_dinucl", 4);

    handler.accumulator()(0, 1) = TestType(100);
    handler.accumulator()(2, 3) = TestType(200);
    handler.reset_accumulator();

    for (std::size_t i = 0; i < handler.accumulator().size(); ++i) {
        REQUIRE(handler.accumulator().data()[i] == TestType(0));
    }
}

TEMPLATE_TEST_CASE("MarkovHandler I/O round-trip", "[model][handler]",
                   double, long double) {
    MarkovHandler<TestType> writer("vd_dinucl", 4);

    auto& acc = writer.accumulator();
    acc(0, 0) = TestType(1);
    acc(0, 1) = TestType(2);
    acc(0, 2) = TestType(3);
    acc(0, 3) = TestType(4);
    acc(1, 0) = TestType(5);
    acc(1, 1) = TestType(3);
    acc(1, 2) = TestType(1);
    acc(1, 3) = TestType(1);
    writer.maximize_likelihood();

    std::stringstream ss;
    writer.write_parameters(ss);

    MarkovHandler<TestType> reader("vd_dinucl", 4);
    reader.read_parameters(ss);

    for (std::size_t i = 0; i < writer.parameters().size(); ++i) {
        REQUIRE_THAT(static_cast<double>(reader.parameters().data()[i]),
                     WithinRel(static_cast<double>(writer.parameters().data()[i]),
                               1e-10));
    }
}

// ─── Polymorphic Access Tests ──────────────────────────────────────────

TEMPLATE_TEST_CASE("Handlers accessible through base pointer", "[model][handler]",
                   double, long double) {
    std::unique_ptr<MarginalHandler<TestType>> handler;

    SECTION("CategoricalHandler through base") {
        handler = std::make_unique<CategoricalHandler<TestType>>("v_gene", 3);

        REQUIRE(handler->name() == "v_gene");
        REQUIRE(handler->parameters().size() == 3);
        REQUIRE(handler->accumulator().size() == 3);

        handler->reset_accumulator();
        handler->maximize_likelihood();

        for (std::size_t i = 0; i < 3; ++i) {
            REQUIRE_THAT(static_cast<double>(handler->parameters()(i)),
                         WithinRel(1.0 / 3.0, 1e-10));
        }
    }

    SECTION("MarkovHandler through base") {
        handler = std::make_unique<MarkovHandler<TestType>>("vd_dinucl", 4);

        REQUIRE(handler->name() == "vd_dinucl");
        REQUIRE(handler->parameters().size() == 16);

        handler->reset_accumulator();
        handler->maximize_likelihood();

        for (std::size_t i = 0; i < 16; ++i) {
            REQUIRE_THAT(static_cast<double>(handler->parameters().data()[i]),
                         WithinRel(0.25, 1e-10));
        }
    }
}

// ─── EM Cycle Test ─────────────────────────────────────────────────────

TEMPLATE_TEST_CASE("CategoricalHandler full EM cycle", "[model][handler]",
                   double, long double) {
    CategoricalHandler<TestType> handler("v_gene", 3);

    // Iteration 1
    handler.reset_accumulator();
    handler.accumulator()(0) += TestType(10);
    handler.accumulator()(1) += TestType(20);
    handler.accumulator()(2) += TestType(30);
    handler.maximize_likelihood();

    REQUIRE_THAT(static_cast<double>(handler.parameters()(0)),
                 WithinRel(1.0 / 6.0, 1e-10));

    // Iteration 2 — strongly favor gene 2
    handler.reset_accumulator();
    handler.accumulator()(0) += TestType(1);
    handler.accumulator()(1) += TestType(1);
    handler.accumulator()(2) += TestType(98);
    handler.maximize_likelihood();

    REQUIRE_THAT(static_cast<double>(handler.parameters()(2)),
                 WithinRel(0.98, 1e-10));

    double total = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        total += static_cast<double>(handler.parameters()(i));
    }
    REQUIRE_THAT(total, WithinRel(1.0, 1e-10));
}
