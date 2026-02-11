#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Errorrate.h>
#include <igor/Model/InferenceEngine.h>
#include <igor/Model/CategoricalHandler.h>

// Include the unit under test
#include "ParallelRecursion.h"

using namespace igor::model;
using namespace std;
using Catch::Matchers::WithinAbs;

// ─── Mocks ─────────────────────────────────────────────────────────────

class MockEvent : public Rec_Event {
public:
    MockEvent(Rec_Event_name name)
        : Rec_Event(Undefined_gene, Undefined_side) {
        this->name = name;
        this->nickname = name;
        this->type = Undefined_t;
        this->priority = 0;
    }

    // Copy required by clone
    MockEvent(const MockEvent& other) = default;

    std::shared_ptr<Rec_Event> copy() override {
        return std::make_shared<MockEvent>(*this);
    }

    int size() const override { return 1; }

    void iterate(double&, Downstream_scenario_proba_bound_map&, const std::string&, const Int_Str&,
                 Index_map&, const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>&,
                 std::shared_ptr<Next_event_ptr>&, Marginal_array_p& updated_marginals, const Marginal_array_p&,
                 const std::unordered_map<Gene_class, std::vector<Alignment_data>>&, Seq_type_str_p_map&,
                 Seq_offsets_map&, std::shared_ptr<Error_rate>&, std::map<size_t, std::shared_ptr<Counter>>&,
                 const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&,
                 Safety_bool_map&, Mismatch_vectors_map&, double&, double&) override {

        // Simulating accumulation
        // Hack for test: write to index 0 of updated_marginals.
        updated_marginals[0] += 2.0;
    }

    // Stubs

    // Correct signature: unordered_map<string, Event_realization>
    // But get_realizations_map is non-virtual accessor in base?
    // Rec_Event.h line 153: const ... get_realizations_map() const { return event_realizations; }
    // It returns the member variable. It is NOT virtual.
    // So I don't override it. I populate event_realizations in constructor if needed.
    // For this test, I don't need realizations if iterate doesn't use them.

    // Pure virtuals I must implement:

    std::queue<int> draw_random_realization(const Marginal_array_p&, std::unordered_map<Rec_Event_name, int>&,
                                          const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>&,
                                          std::unordered_map<Seq_type, std::string>&, std::mt19937_64&) const override {
        return {};
    }

    void write2txt(std::ofstream&) override {}

    void add_to_marginals(long double, Marginal_array_p &) const override {}

    bool has_effect_on(Seq_type) const override { return false; }

    void iterate_initialize_Len_proba(Seq_type considered_junction,
                                              std::map<int, double> &length_best_proba_map,
                                              std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                              double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                              Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                              int &seq_len) const override {}

    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                            const Marginal_array_p &model_parameters_point,
                                            Index_map &base_index_map) override {}

    // Note: Rec_Event destructor is virtual.
};

class MockErrorRate : public Error_rate {
public:
    MockErrorRate(double likelihood) {
        this->seq_likelihood = likelihood;
    }

    // Stubs
    double compare_sequences_error_prob(double, const std::string&, Seq_type_str_p_map&, const Seq_offsets_map&,
                                      const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>&,
                                      Mismatch_vectors_map&, double&, double&) override { return 0; }
    void update() override {}
    void add_to_norm_counter() override {}
    void clean_seq_counters() override {}
    void write2txt(std::ofstream&) override {}
    std::shared_ptr<Error_rate> copy() const override { return std::make_shared<MockErrorRate>(*this); }
    std::string type() const override { return "Mock"; }
    Error_rate* add_checked(Error_rate*) override { return this; }
    const double& get_err_rate_upper_bound(size_t, size_t) override { static double d=0; return d; }
    void build_upper_bound_matrix(size_t, size_t) override {}
    int get_number_non_zero_likelihood_seqs() const override { return 1; }
    std::queue<int> generate_errors(std::string&, std::mt19937_64&) const override { return {}; }
};


TEST_CASE("ParallelRecursionTest: CountScenarioToEngineAccumulatesCorrectly", "[integration][parallel]") {
    // 1. Setup Model Parameters with 1 Mock Event
    Model_Parms parms;
    auto mock_event = std::make_shared<MockEvent>("mock_event");
    parms.add_event(mock_event);

    auto mock_err = std::make_shared<MockErrorRate>(0.5); // Prob(data) = 0.5
    parms.set_error_ratep(mock_err);

    // 2. Setup Inference Engine
    InferenceEngine<long double> engine;
    // Register handler for "mock_event"
    // CategoricalHandler(std::string name, std::vector<std::size_t> shape)
    auto handler = std::make_unique<CategoricalHandler<long double>>("mock_event", std::vector<size_t>{1});
    engine.register_handler("mock_event", std::move(handler));

    // 3. Setup context maps
    // Model_marginals needed to generate index_map and provide params
    Model_marginals marginals(parms);
    // Initialize marginals (params) to something.
    // Note: marginals(parms) allocates array.
    if(marginals.get_length() > 0)
        marginals.marginal_array_smart_p[0] = 1.0;

    auto model_queue = parms.get_model_queue();
    auto index_map = marginals.get_index_map(parms, model_queue);
    auto offset_map = marginals.get_offsets_map(parms, model_queue);
    auto events_map = parms.get_events_map();

    // 4. Prepare inputs
    std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>> seq_tuple;
    std::get<0>(seq_tuple) = 0;
    std::get<1>(seq_tuple) = "ACGT";
    // Alignments empty

    // 5. Run the function (count_scenario_to_engine now creates its own next_event_ptr)
    double likelihood = count_scenario_to_engine(
        seq_tuple,
        engine,
        parms,
        marginals,
        events_map,
        index_map,
        offset_map,
        mock_err
    );

    // 6. Assertions
    REQUIRE(likelihood == 0.5);

    // Check accumulation
    // Mock adds 2.0. Likelihood is 0.5.
    // Result should be 2.0 * (1/0.5) = 4.0.
    auto& handler_ref = engine.handler("mock_event");
    auto& acc = handler_ref.accumulator();

    REQUIRE(acc.size() == 1);

    REQUIRE_THAT(static_cast<double>(acc.data()[0]), WithinAbs(4.0, 1e-9));
}
