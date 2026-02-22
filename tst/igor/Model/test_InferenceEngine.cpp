/*
 * test_InferenceEngine.cpp
 *
 *  Created on: Feb 10, 2026
 *  Updated on: Feb 21, 2026
 *
 *  Integration tests for InferenceEngine.
 *  Uses the real Mouse TCR beta model loaded from files.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Model/InferenceEngine.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/InferenceHandler.h>

#include <numeric>

using namespace igor::model;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ─── Helper: load the mouse model ─────────────────────────────────────

static std::shared_ptr<RecombinationModel<double>> make_mouse_model() {
    const std::string parms_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_parms.txt";
    const std::string marginals_path =
        std::string(IGOR_MODELS_DIR) + "/mouse/tcr_beta/models/model_marginals.txt";

    try {
        auto m = recombination_model_from_files<double>(parms_path, marginals_path);
        return std::make_shared<RecombinationModel<double>>(std::move(m));
    } catch (...) {
        return nullptr;
    }
}

// ─── Construction ──────────────────────────────────────────────────────

TEST_CASE("InferenceEngine construction from RecombinationModel",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    SECTION("size matches topology") {
        REQUIRE(engine.size() == model->topology().size());
    }
    SECTION("each handler name matches topology event name") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            const auto& ev_name = model->topology().event(uid)->get_nickname();
            REQUIRE(engine.handler(uid).name() == ev_name);
        }
    }
    SECTION("handler weights point to model weights (same memory)") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            REQUIRE(&engine.handler(uid).weights() == &model->weight(uid));
        }
    }
}

// ─── Handler Lookup ────────────────────────────────────────────────────

TEST_CASE("InferenceEngine handler lookup by name and uid",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    SECTION("lookup by name matches lookup by uid") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            const auto& name = engine.handler(uid).name();
            REQUIRE(&engine.handler(name) == &engine.handler(uid));
        }
    }
    SECTION("hasHandler returns true for known events") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            REQUIRE(engine.hasHandler(engine.handler(uid).name()));
        }
    }
}

// ─── Reset Accumulators ────────────────────────────────────────────────

TEST_CASE("InferenceEngine resetAccumulators zeros all",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    // Dirty an accumulator
    engine.handler(0).accumulator().data()[0] = 42.0;
    engine.resetAccumulators();

    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(engine.size()); ++uid) {
        const auto& acc = engine.handler(uid).accumulator();
        for (std::size_t i = 0; i < acc.size(); ++i) {
            REQUIRE(acc.data()[i] == 0.0);
        }
    }
}

// ─── Update Parameters (M-step) ───────────────────────────────────────

TEST_CASE("InferenceEngine updateParameters writes back into model",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    // Set first handler's accumulator to known values
    auto& h0 = engine.handler(0);
    std::size_t n = h0.accumulator().size();
    for (std::size_t i = 0; i < n; ++i) {
        h0.accumulator().data()[i] = static_cast<double>(i + 1);
    }

    engine.updateParameters();

    // Weights should have been normalised
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += model->weight(0).data()[i];
    }
    REQUIRE_THAT(sum, WithinAbs(1.0, 1e-10));
}

// ─── Combine Accumulators ──────────────────────────────────────────────

TEST_CASE("InferenceEngine combineAccumulators merges correctly",
          "[model][engine]") {
    auto model1 = make_mouse_model();
    auto model2 = make_mouse_model();
    if (!model1 || !model2) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine1(model1);
    InferenceEngine<double> engine2(model2);

    // Set accumulators
    auto& acc1 = engine1.handler(0).accumulator();
    auto& acc2 = engine2.handler(0).accumulator();
    for (std::size_t i = 0; i < acc1.size(); ++i) {
        acc1.data()[i] = 1.0;
        acc2.data()[i] = 2.0;
    }

    engine1.combineAccumulators(engine2);

    for (std::size_t i = 0; i < acc1.size(); ++i) {
        REQUIRE_THAT(acc1.data()[i], WithinAbs(3.0, 1e-10));
    }
}

// ─── forEachHandler ──────────────────────────────────────────────────

TEST_CASE("InferenceEngine forEachHandler visits all",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    std::size_t count = 0;
    engine.forEachHandler([&](const std::string& name, const InferenceHandler<double>& h) {
        REQUIRE(!name.empty());
        REQUIRE(h.weights().size() > 0);
        ++count;
    });
    REQUIRE(count == engine.size());
}

// ─── orderedHandlers ───────────────────────────────────────────────────

TEST_CASE("InferenceEngine orderedHandlers visits all in topological order",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);
    auto ordered = engine.orderedHandlers();

    SECTION("orderedHandlers visits every handler") {
        REQUIRE(ordered.size() == engine.size());
    }
    SECTION("each element is a valid handler") {
        for (const auto& h : ordered) {
            REQUIRE(h != nullptr);
            REQUIRE(!h->name().empty());
        }
    }
    SECTION("parents are visited before children") {
        // Track which uids have been visited
        std::vector<bool> visited(engine.size(), false);
        for (const auto& h : ordered) {
            auto uid = h->uid();
            // All parents must have been visited already
            for (const auto& p : engine.parents(uid)) {
                REQUIRE(p != nullptr);
                REQUIRE(visited[p->uid()]);
            }
            visited[uid] = true;
        }
    }
}

// ─── range-based for ──────────────────────────────────────────────────

TEST_CASE("InferenceEngine supports range-based for",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    std::size_t count = 0;
    for (const auto& h : engine) {
        REQUIRE(h != nullptr);
        REQUIRE(!h->name().empty());
        ++count;
    }
    REQUIRE(count == engine.size());
}

// ─── parents / children ────────────────────────────────────────────────

TEST_CASE("InferenceEngine parents and children are consistent",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    SECTION("root events have no parents") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            auto p = engine.parents(uid);
            auto topo_parents = model->topology().parentsIds(uid);
            REQUIRE(p.size() == topo_parents.size());
        }
    }
    SECTION("parent/child relationship is symmetric") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            for (const auto& child_ptr : engine.children(uid)) {
                REQUIRE(child_ptr != nullptr);
                auto child_uid = child_ptr->uid();
                // uid must appear in child's parents
                bool found = false;
                for (const auto& p : engine.parents(child_uid)) {
                    if (p->uid() == uid) { found = true; break; }
                }
                REQUIRE(found);
            }
        }
    }
    SECTION("children navigators return valid handlers") {
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(engine.size()); ++uid) {
            for (const auto& c : engine.children(uid)) {
                REQUIRE(c != nullptr);
                REQUIRE(!c->name().empty());
            }
        }
    }
}

// ─── Type Aliases ──────────────────────────────────────────────────────

TEST_CASE("InferenceEngine type aliases", "[model][engine]") {
    static_assert(std::is_same_v<Engine, InferenceEngine<double>>);
    static_assert(std::is_same_v<LegacyEngine, InferenceEngine<long double>>);
}

// ─── run() ─────────────────────────────────────────────────────────────

TEST_CASE("InferenceEngine run() executes reset-eStep-update cycle",
          "[model][engine]") {
    auto model = make_mouse_model();
    if (!model) SKIP("Mouse TCR beta model files not found");

    InferenceEngine<double> engine(model);

    // Dirty an accumulator before run() — it should be reset inside run()
    engine.handler(0).accumulator().data()[0] = 999.0;

    bool e_step_called = false;

    engine.run([&](InferenceEngine<double>& eng) {
        e_step_called = true;

        // Verify accumulators were already reset
        REQUIRE(eng.handler(0).accumulator().data()[0] == 0.0);

        // Simulate accumulation: put known counts into handler 0
        auto& acc = eng.handler(0).accumulator();
        for (std::size_t i = 0; i < acc.size(); ++i) {
            acc.data()[i] = static_cast<double>(i + 1);
        }
    });

    REQUIRE(e_step_called);

    // After run(), the M-step should have normalised the weights
    double sum = 0.0;
    auto& w = model->weight(0);
    for (std::size_t i = 0; i < w.size(); ++i) {
        sum += w.data()[i];
    }
    REQUIRE_THAT(sum, WithinAbs(1.0, 1e-10));
}
