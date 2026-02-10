/*
 * test_InferenceEngine.cpp
 *
 *  Created on: Feb 10, 2026
 *
 *  Unit and integration tests for InferenceEngine.
 *  Uses TEMPLATE_TEST_CASE to validate both double and long double.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <igor/Model/InferenceEngine.h>

#include <sstream>
#include <numeric>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ─── Helper: build a minimal VDJ-like descriptor set ───────────────────

static std::vector<EventDescriptor> make_vdj_descriptors() {
    // Mimics a simplified VDJ model:
    //   v_choice: 5 V genes (GeneChoice_t = 0), no parents
    //   d_gene:   3 D genes (GeneChoice_t = 0), parent=v_choice
    //   j_choice: 4 J genes (GeneChoice_t = 0), parent=d_gene
    //   v_3_del:  10 del values (Deletion_t = 2), parent=v_choice
    //   d_5_del:  8 del (Deletion_t = 2), parent=d_gene
    //   d_3_del:  8 del (Deletion_t = 2), parent=d_gene
    //   j_5_del:  12 del (Deletion_t = 2), parent=j_choice
    //   vd_ins:   20 insertion lengths (Insertion_t = 1), no parents
    //   dj_ins:   20 insertion lengths (Insertion_t = 1), no parents
    //   vd_dinucl: 4x4 Markov (Dinuclmarkov_t = 3), no parents
    //   dj_dinucl: 4x4 Markov (Dinuclmarkov_t = 3), no parents

    return {
        {"v_choice",  0, 0, 0, {5}},           // 1D
        {"d_gene",    0, 0, 0, {3, 5}},         // conditioned on v_choice
        {"j_choice",  0, 0, 0, {4, 3}},         // conditioned on d_gene
        {"v_3_del",   2, 0, 2, {10, 5}},        // conditioned on v_choice
        {"d_5_del",   2, 0, 0, {8, 3}},         // conditioned on d_gene
        {"d_3_del",   2, 0, 2, {8, 3}},         // conditioned on d_gene
        {"j_5_del",   2, 0, 0, {12, 4}},        // conditioned on j_choice
        {"vd_ins",    1, 0, 0, {20}},           // 1D
        {"dj_ins",    1, 0, 0, {20}},           // 1D
        {"vd_dinucl", 3, 0, 0, {4, 4}},         // 4x4 transition
        {"dj_dinucl", 3, 0, 0, {4, 4}},         // 4x4 transition
    };
}

// ─── Construction Tests ────────────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceEngine default construction", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;

    REQUIRE(engine.size() == 0);
    REQUIRE(engine.event_names().empty());
}

TEMPLATE_TEST_CASE("InferenceEngine construction from descriptors", "[model][engine]",
                   double, long double) {
    auto descs = make_vdj_descriptors();
    InferenceEngine<TestType> engine(descs);

    SECTION("correct number of handlers") {
        REQUIRE(engine.size() == 11);
    }

    SECTION("all events accessible by name") {
        for (const auto& desc : descs) {
            REQUIRE(engine.has_handler(desc.name));
        }
    }

    SECTION("event order preserved") {
        const auto& names = engine.event_names();
        REQUIRE(names[0] == "v_choice");
        REQUIRE(names[1] == "d_gene");
        REQUIRE(names[2] == "j_choice");
        REQUIRE(names[10] == "dj_dinucl");
    }

    SECTION("1D categorical handlers have correct size") {
        REQUIRE(engine.handler("v_choice").parameters().size() == 5);
        REQUIRE(engine.handler("vd_ins").parameters().size() == 20);
        REQUIRE(engine.handler("dj_ins").parameters().size() == 20);
    }

    SECTION("2D conditional handlers have correct size") {
        // d_gene: 3 × 5 = 15
        REQUIRE(engine.handler("d_gene").parameters().size() == 15);
        // v_3_del: 10 × 5 = 50
        REQUIRE(engine.handler("v_3_del").parameters().size() == 50);
        // j_5_del: 12 × 4 = 48
        REQUIRE(engine.handler("j_5_del").parameters().size() == 48);
    }

    SECTION("markov handlers have correct size") {
        // 4x4 = 16 elements
        REQUIRE(engine.handler("vd_dinucl").parameters().size() == 16);
        REQUIRE(engine.handler("dj_dinucl").parameters().size() == 16);
    }
}

// ─── Manual Registration Tests ─────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceEngine manual registration", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;

    engine.register_handler("v_choice",
        std::make_unique<CategoricalHandler<TestType>>("v_choice", 5));
    engine.register_handler("vd_dinucl",
        std::make_unique<MarkovHandler<TestType>>("vd_dinucl", 4));

    REQUIRE(engine.size() == 2);
    REQUIRE(engine.has_handler("v_choice"));
    REQUIRE(engine.has_handler("vd_dinucl"));
    REQUIRE_FALSE(engine.has_handler("nonexistent"));
}

TEMPLATE_TEST_CASE("InferenceEngine duplicate registration throws", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;

    engine.register_handler("v_choice",
        std::make_unique<CategoricalHandler<TestType>>("v_choice", 5));

    REQUIRE_THROWS_AS(
        engine.register_handler("v_choice",
            std::make_unique<CategoricalHandler<TestType>>("v_choice", 5)),
        std::runtime_error);
}

TEMPLATE_TEST_CASE("InferenceEngine access nonexistent throws", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;

    REQUIRE_THROWS_AS(engine.handler("nonexistent"), std::runtime_error);
}

// ─── EM Cycle Tests ────────────────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceEngine reset_accumulators", "[model][engine]",
                   double, long double) {
    auto descs = make_vdj_descriptors();
    InferenceEngine<TestType> engine(descs);

    // Accumulate some values directly into the accumulator tensor
    auto& v_acc = engine.handler("v_choice").accumulator();
    v_acc(0) += TestType(100);
    v_acc(1) += TestType(200);

    auto& m_acc = engine.handler("vd_dinucl").accumulator();
    m_acc(0, 1) += TestType(50);

    // Reset
    engine.reset_accumulators();

    // All accumulators should be zero
    engine.for_each_handler([](const std::string& /*name*/,
                               const MarginalHandler<TestType>& h) {
        for (std::size_t i = 0; i < h.accumulator().size(); ++i) {
            REQUIRE(h.accumulator().data()[i] == TestType(0));
        }
    });
}

TEMPLATE_TEST_CASE("InferenceEngine update_parameters", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;
    engine.register_handler("gene",
        std::make_unique<CategoricalHandler<TestType>>("gene", 3));
    engine.register_handler("markov",
        std::make_unique<MarkovHandler<TestType>>("markov", 4));

    // Accumulate into gene handler
    auto& gene_acc = engine.handler("gene").accumulator();
    gene_acc(0) = TestType(10);
    gene_acc(1) = TestType(20);
    gene_acc(2) = TestType(70);

    // Accumulate into markov handler
    auto& mk_acc = engine.handler("markov").accumulator();
    mk_acc(0, 0) = TestType(25);
    mk_acc(0, 1) = TestType(75);

    // M-step
    engine.update_parameters();

    SECTION("gene parameters normalized") {
        REQUIRE_THAT(static_cast<double>(engine.handler("gene").parameters()(0)),
                     WithinRel(0.1, 1e-10));
        REQUIRE_THAT(static_cast<double>(engine.handler("gene").parameters()(1)),
                     WithinRel(0.2, 1e-10));
        REQUIRE_THAT(static_cast<double>(engine.handler("gene").parameters()(2)),
                     WithinRel(0.7, 1e-10));
    }

    SECTION("markov parameters normalized per row") {
        REQUIRE_THAT(static_cast<double>(engine.handler("markov").parameters()(0, 0)),
                     WithinRel(0.25, 1e-10));
        REQUIRE_THAT(static_cast<double>(engine.handler("markov").parameters()(0, 1)),
                     WithinRel(0.75, 1e-10));
    }
}

TEMPLATE_TEST_CASE("InferenceEngine full EM cycle", "[model][engine]",
                   double, long double) {
    InferenceEngine<TestType> engine;
    engine.register_handler("gene",
        std::make_unique<CategoricalHandler<TestType>>("gene", 4));

    // Iteration 1
    engine.reset_accumulators();
    auto& acc1 = engine.handler("gene").accumulator();
    acc1(0) = TestType(25);
    acc1(1) = TestType(25);
    acc1(2) = TestType(25);
    acc1(3) = TestType(25);
    engine.update_parameters();

    // Should be uniform
    for (std::size_t i = 0; i < 4; ++i) {
        REQUIRE_THAT(static_cast<double>(engine.handler("gene").parameters()(i)),
                     WithinRel(0.25, 1e-10));
    }

    // Iteration 2 — skewed
    engine.reset_accumulators();
    auto& acc2 = engine.handler("gene").accumulator();
    acc2(0) = TestType(90);
    acc2(1) = TestType(5);
    acc2(2) = TestType(3);
    acc2(3) = TestType(2);
    engine.update_parameters();

    REQUIRE_THAT(static_cast<double>(engine.handler("gene").parameters()(0)),
                 WithinRel(0.9, 1e-10));

    // Sum should still be 1
    double total = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        total += static_cast<double>(engine.handler("gene").parameters()(i));
    }
    REQUIRE_THAT(total, WithinRel(1.0, 1e-10));
}

// ─── I/O Tests ─────────────────────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceEngine write/read round-trip", "[model][engine]",
                   double, long double) {
    // Build and configure engine
    InferenceEngine<TestType> writer;
    writer.register_handler("v_choice",
        std::make_unique<CategoricalHandler<TestType>>("v_choice", 3));
    writer.register_handler("vd_dinucl",
        std::make_unique<MarkovHandler<TestType>>("vd_dinucl", 4));

    // Set non-trivial parameters via accumulator
    auto& wv_acc = writer.handler("v_choice").accumulator();
    wv_acc(0) = TestType(10);
    wv_acc(1) = TestType(30);
    wv_acc(2) = TestType(60);
    writer.handler("v_choice").maximize_likelihood();

    auto& wm_acc = writer.handler("vd_dinucl").accumulator();
    wm_acc(0, 0) = TestType(1);
    wm_acc(0, 1) = TestType(2);
    wm_acc(0, 2) = TestType(3);
    wm_acc(0, 3) = TestType(4);
    writer.handler("vd_dinucl").maximize_likelihood();

    // Write
    std::stringstream ss;
    writer.write_marginals(ss);

    // Read into fresh engine with same structure
    InferenceEngine<TestType> reader;
    reader.register_handler("v_choice",
        std::make_unique<CategoricalHandler<TestType>>("v_choice", 3));
    reader.register_handler("vd_dinucl",
        std::make_unique<MarkovHandler<TestType>>("vd_dinucl", 4));

    reader.read_marginals(ss);

    // Compare
    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE_THAT(
            static_cast<double>(reader.handler("v_choice").parameters()(i)),
            WithinRel(
                static_cast<double>(writer.handler("v_choice").parameters()(i)),
                1e-10));
    }

    for (std::size_t from = 0; from < 4; ++from) {
        for (std::size_t to = 0; to < 4; ++to) {
            REQUIRE_THAT(
                static_cast<double>(reader.handler("vd_dinucl").parameters()(from, to)),
                WithinRel(
                    static_cast<double>(writer.handler("vd_dinucl").parameters()(from, to)),
                    1e-10));
        }
    }
}

// ─── for_each_handler Tests ────────────────────────────────────────────

TEMPLATE_TEST_CASE("InferenceEngine for_each_handler", "[model][engine]",
                   double, long double) {
    auto descs = make_vdj_descriptors();
    InferenceEngine<TestType> engine(descs);

    SECTION("const iteration counts all handlers") {
        std::size_t count = 0;
        engine.for_each_handler([&count](const std::string& /*name*/,
                                         const MarginalHandler<TestType>& /*h*/) {
            ++count;
        });
        REQUIRE(count == 11);
    }

    SECTION("mutable iteration can modify parameters") {
        // Zero out all parameters
        engine.for_each_handler([](const std::string& /*name*/,
                                   MarginalHandler<TestType>& h) {
            std::fill(h.parameters().begin(), h.parameters().end(), TestType(0));
        });

        // Verify
        engine.for_each_handler([](const std::string& /*name*/,
                                   const MarginalHandler<TestType>& h) {
            for (std::size_t i = 0; i < h.parameters().size(); ++i) {
                REQUIRE(h.parameters().data()[i] == TestType(0));
            }
        });
    }

    SECTION("iteration order matches event_names") {
        std::vector<std::string> visited;
        engine.for_each_handler([&visited](const std::string& name,
                                            const MarginalHandler<TestType>& /*h*/) {
            visited.push_back(name);
        });
        REQUIRE(visited == engine.event_names());
    }
}

// ─── Type Alias Tests ──────────────────────────────────────────────────

TEST_CASE("Engine and LegacyEngine aliases compile and work", "[model][engine]") {
    SECTION("Engine is InferenceEngine<double>") {
        Engine engine;
        engine.register_handler("test",
            std::make_unique<CategoricalHandler<double>>("test", 3));
        REQUIRE(engine.size() == 1);
    }

    SECTION("LegacyEngine is InferenceEngine<long double>") {
        LegacyEngine engine;
        engine.register_handler("test",
            std::make_unique<CategoricalHandler<long double>>("test", 3));
        REQUIRE(engine.size() == 1);
    }
}
