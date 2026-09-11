/*
 * test_proba_bound_benchmark.cpp
 *
 *  Catch2 benchmarks for the junction-length probability bound precomputation --
 *  the reverse initialize_Len_proba_bound() sweep GenModel runs once per thread per
 *  EM iteration, before the per-sequence loop.
 *
 *  Why this is worth measuring rather than assuming: the cost is *model-only*, so it
 *  amortises to nothing over a large batch but is paid in full by a CLI user querying
 *  one sequence at a time. Measured on the VDJ TRB regression corpus it is ~0.7 s per
 *  thread per iteration, which dominates a single-sequence query by orders of
 *  magnitude. See docs/PROBA_BOUND_MACHINERY.md.
 *
 *  Hidden by default (tag starts with '.'); run explicitly with:
 *      igor_tests "[.benchmark][proba_bound]"
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>

#include <filesystem>
#include <forward_list>
#include <list>
#include <memory>
#include <queue>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef IGOR_SOURCE_DIR
#error "IGOR_SOURCE_DIR must be defined (set by CMake)"
#endif

namespace {

const std::string kModelsDir = std::string(IGOR_SOURCE_DIR) + "/models";
const std::string kRegressionDir = std::string(IGOR_SOURCE_DIR) + "/scripts/tests/data/input";

/**
 * Everything GenModel builds before the bound sweep, held so the benchmark times the
 * sweep alone rather than model loading.
 *
 * Marginals are the ones the parms imply rather than a loaded file: the sweep enumerates
 * every realization whatever the probabilities are, so its cost depends on the model's
 * *shape*, not its values.
 */
class BoundInitFixture
{
public:
    explicit BoundInitFixture(const std::string &parms_path)
        : parms_(), marginals_()
    {
        parms_.read_model_parms(parms_path);
        marginals_ = Model_marginals(parms_);
        marginals_.uniform_initialize(parms_);

        model_queue_ = parms_.get_model_queue();
        index_by_name_ = marginals_.get_index_map(parms_, model_queue_);

        events_ = parms_.get_event_list();
        index_map_ = Index_map(events_.size());
        for (auto &event : events_) {
            const int event_index = event->get_event_identifier();
            index_map_.request_layer(event_index);
            index_map_.set(event_index, index_by_name_.at(event->get_name()), 0);

            const std::size_t event_size = marginals_.get_event_size(event, parms_);
            event->set_event_marginal_size(event_size);
            event->set_crude_upper_bound_proba(index_by_name_.at(event->get_name()), event_size,
                                               marginals_.marginal_array_smart_p);
        }
        events_map_ = parms_.get_events_map();
        const auto offset_map = marginals_.get_offsets_map(parms_, model_queue_);

        //The forward initialize_event() pass, which the sweep depends on: it is where each
        //event resolves *which* junctions it reads a bound for. Skipping it used to be harmless
        //because the sweep picked its junctions from an enum switch; since S4c it would measure
        //an empty sweep. Part of the fixture rather than of the timed region, exactly as in
        //GenModel::infer_model, where it runs before the bound loop.
        Safety_bool_map safety_set(3);
        Seq_type_str_p_map constructed_sequences(parms_.get_seq_type_registry());
        Mismatch_vectors_map mismatches_lists(parms_.get_seq_type_registry());
        Seq_offsets_map seq_offsets(parms_.get_seq_type_registry());
        Downstream_scenario_proba_bound_map downstream_proba_map(parms_.get_seq_type_registry());
        downstream_proba_map.init_first_layer(1.0);
        std::shared_ptr<Error_rate> err_rate = parms_.get_err_rate_p();

        std::unordered_set<Rec_Event_name> processed_events;
        std::queue<std::shared_ptr<Rec_Event>> init_queue = model_queue_;
        while (!init_queue.empty()) {
            std::shared_ptr<Rec_Event> event = init_queue.front();
            init_queue.pop();
            event->initialize_event(processed_events, events_map_, offset_map, downstream_proba_map,
                                    constructed_sequences, safety_set, err_rate, mismatches_lists, seq_offsets,
                                    index_map_);
        }
    }

    /// The reverse sweep, exactly as GenModel::infer_model runs it.
    void run_bound_sweep()
    {
        std::stack<std::shared_ptr<Rec_Event>> init_stack;
        {
            std::queue<std::shared_ptr<Rec_Event>> walk = model_queue_;
            while (!walk.empty()) {
                init_stack.push(walk.front());
                walk.pop();
            }
        }

        double downstream_proba_bound = 1;
        std::forward_list<double *> updated_proba_list;
        while (!init_stack.empty()) {
            std::shared_ptr<Rec_Event> event = init_stack.top();
            init_stack.pop();

            //Each event is handed the queue *after* itself, which is what makes the
            //profiles a suffix family rather than copies of one another.
            std::queue<std::shared_ptr<Rec_Event>> remaining = model_queue_;
            while (remaining.front() != event) {
                remaining.pop();
            }
            remaining.pop();

            event->initialize_crude_scenario_proba_bound(downstream_proba_bound, updated_proba_list,
                                                         events_map_);
            event->initialize_Len_proba_bound(remaining, marginals_.marginal_array_smart_p, index_map_);
        }
    }

    std::size_t event_count() const { return events_.size(); }

private:
    Model_Parms parms_;
    Model_marginals marginals_;
    std::queue<std::shared_ptr<Rec_Event>> model_queue_;
    std::unordered_map<Rec_Event_name, int> index_by_name_;
    std::list<std::shared_ptr<Rec_Event>> events_;
    Index_map index_map_{0};
    Events_map events_map_;
};

} // namespace

TEST_CASE("Junction-length bound initialization", "[.benchmark][proba_bound]")
{
    struct Case {
        const char *label;
        std::string parms_path;
        const char *topology;
    };

    //Ordered cheapest first, so a partial run still says something.
    const std::vector<Case> cases{
            {"human TCR-alpha", kModelsDir + "/human/tcr_alpha/models/model_parms.txt", "VJ"},
            {"TRB regression corpus", kRegressionDir + "/TRB_model_parms.txt", "VDJ"},
            {"human TCR-beta", kModelsDir + "/human/tcr_beta/models/model_parms.txt", "VDJ"},
            {"human BCR-heavy", kModelsDir + "/human/bcr_heavy/models/model_parms.txt", "VDJ"},
    };

    for (const Case &model : cases) {
        if (!std::filesystem::exists(model.parms_path)) {
            WARN("skipping " << model.label << ": " << model.parms_path << " not found");
            continue;
        }

        std::ostringstream name;
        name << model.label << " (" << model.topology << ")";

        //Built outside the timed region: this measures the sweep, not model loading.
        BoundInitFixture fixture(model.parms_path);
        REQUIRE(fixture.event_count() > 0);

        BENCHMARK(name.str()) { return fixture.run_bound_sweep(); };
    }
}

//A tandem-D model cannot be benchmarked yet: iterate_initialize_Len_proba is reached
//through Seq_type-named maps, so a D1D2 junction throws before the sweep runs. Add the
//case with S4c, which is the step that removes that ceiling -- the benchmark is then the
//evidence that the generalisation did not cost throughput.
