#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <queue>
#include <stack>
#include <memory>
#include <forward_list>
#include <iostream>

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Utils.h>
#include <igor/Core/Counter.h>
#include <igor/Core/Errorrate.h>
#include <igor/Core/Genechoice.h> // For casting if needed
#include <igor/Core/Dinuclmarkov.h>

#include <igor/Model/InferenceEngine.h>
#include <igor/Core/LegacyBridge.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

/**
 * @brief Helper function to transfer accumulated marginals from a flat legacy array to the InferenceEngine.
 *
 * This mimics LegacyBridge::import_from_legacy but adds to the accumulator instead of replacing parameters.
 */
template <typename T>
void accumulate_from_legacy_flat(igor::model::InferenceEngine<T>& engine,
                                 const Model_marginals& marginals,
                                 const Model_Parms& parms) {

    // Get index map to locate events in the marginal array
    // Note: Model_marginals::get_index_map might be expensive if called per sequence.
    // Optimization: Should be passed in. But for validation speed isn't critical.
    // Actually, get_index_map takes parms and queue.
    auto model_queue = parms.get_model_queue();
    auto index_map = marginals.get_index_map(parms, model_queue);

    // Access the raw marginal array
    const long double* marginal_array = marginals.marginal_array_smart_p.get();

    // Iterate through each handler in the engine
    // We cannot use engine.for_each_handler easily because we need to modify non-const accumulator
    // But engine.get_handler returns a reference.

    // We iterate over the model events to match legacy order/indexing
    while(!model_queue.empty()){
        auto event = model_queue.front();
        model_queue.pop();

        std::string nickname = event->get_nickname();

        if (!engine.has_handler(nickname)) {
             continue;
        }

        auto& handler = engine.handler(nickname);

        std::string full_name = event->get_name();
        if (index_map.find(full_name) == index_map.end()) {
            std::cerr << "ERROR: Event not found in index_map: " << full_name << std::endl;
            std::cerr << "Available keys in index_map:" << std::endl;
            for (const auto& [key, val] : index_map) {
                std::cerr << "  " << key << " -> " << val << std::endl;
            }
            continue;
        }
        int base_index = index_map.at(full_name);

        // Get accumulator (T)
        auto& accumulator = handler.accumulator();
        size_t size = accumulator.size();

        // Add values
        T* acc_data = accumulator.data();
        for (size_t i = 0; i < size; ++i) {
            acc_data[i] += static_cast<T>(marginal_array[base_index + i]);
        }
    }
}


/**
 * @brief Test-only function to run the recursion for a single sequence and accumulate to the InferenceEngine.
 *
 * This mimics the body of the loop in GenModel::infer_model.
 * It uses the 'Feed Flat Array' strategy:
 * 1. Create a temporary Model_marginals for this sequence.
 * 2. Run the legacy Rec_Event::iterate() to fill it.
 * 3. Transfer the result to the InferenceEngine.
 *
 * Note: This function is NOT thread-safe if 'engine' is shared across threads.
 * For parallel execution, use thread-local engines and merge them using engine.combine_accumulators().
 */
inline double count_scenario_to_engine(
    const std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>& seq_tuple,
    igor::model::InferenceEngine<long double>& engine,
    const Model_Parms& parms,
    const Model_marginals& current_params_marginals, // Added parameter source
    // Context structures
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>>& events_map,
    const std::unordered_map<Rec_Event_name, int>& index_map,
    const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>& offset_map,
    std::shared_ptr<Error_rate> err_rate,
    double likelihood_threshold = 1e-25,
    double proba_threshold_factor = 0.001
) {
    using namespace std;


    // Unpack sequence
    string nt_seq = get<1>(seq_tuple);
    const auto& alignments = get<2>(seq_tuple);
    Int_Str int_sequence = nt2int(nt_seq);

    // Create separate Next_event_ptr array for engine path (isolated from legacy path)
    size_t num_events = parms.get_event_list().size();
    Next_event_ptr* raw_arr_engine = new Next_event_ptr[num_events];
    for(size_t k = 0; k < num_events; ++k) raw_arr_engine[k] = nullptr;
    shared_ptr<Next_event_ptr> next_event_ptr_arr_engine(
        raw_arr_engine,
        [](Next_event_ptr* p) { delete[] p; }
    );

    // Populate next_event_ptr_arr_engine
    auto model_queue = parms.get_model_queue();
    queue<shared_ptr<Rec_Event>> engine_queue = model_queue;
    while (!engine_queue.empty()) {
        shared_ptr<Rec_Event> event = engine_queue.front();
        engine_queue.pop();
        if (!engine_queue.empty()) {
            int id = event->get_event_identifier();
            next_event_ptr_arr_engine.get()[id] = engine_queue.front().get();
        } else {
            next_event_ptr_arr_engine.get()[event->get_event_identifier()] = nullptr;
        }
    }

    // 1. Initialize single seq marginals (accumulators)
    // We create a fresh accumulator array for this sequence
    Model_marginals single_seq_marginals(parms);
    // Zero out the array
    size_t arr_size = single_seq_marginals.get_length();
    for(size_t i=0; i<arr_size; ++i) single_seq_marginals.marginal_array_smart_p[i] = 0.0;

    double init_proba = 1.0;
    double max_proba_scenario = likelihood_threshold / proba_threshold_factor;

    // Initialize auxiliary maps
    Safety_bool_map safety_set(3);
    Seq_type_str_p_map constructed_sequences(6);
    Mismatch_vectors_map mismatches_lists(6);
    Seq_offsets_map seq_offsets(6, 3);

    // Downstream proba map
    Downstream_scenario_proba_bound_map downstream_proba_map(6);
    downstream_proba_map.init_first_layer(1.0);

    // Index Map wrapper
    auto events_list = parms.get_event_list();
    Index_map index_mapp_fast(events_list.size());
    for (const auto& event : events_list) {
         int event_identifiers = event->get_event_identifier();
         index_mapp_fast.request_memory_layer(event_identifiers);

         // Defensive: check if event name exists in index_map
         std::string event_name = event->get_name();
         if (index_map.find(event_name) == index_map.end()) {
             std::cerr << "WARNING: Event '" << event_name << "' not found in index_map, skipping" << std::endl;
             std::cerr << "Available keys in index_map:" << std::endl;
             for (const auto& [key, val] : index_map) {
                 std::cerr << "  '" << key << "' -> " << val << std::endl;
             }
             continue;
         }
         index_mapp_fast.set_value(event_identifiers, index_map.at(event_name), 0);
    }

    // Dummy pointers
    std::map<size_t, std::shared_ptr<Counter>> counters_dummy;

    // Get first event (already retrieved model_queue above)
    auto first_event = model_queue.front();

    // 2. Call iterate
    try {
        first_event->iterate(
            init_proba,
            downstream_proba_map,
            nt_seq,
            int_sequence,
            index_mapp_fast,
            // casting const ref to ref? iterate signature expects const ref for offset_map
            offset_map,
            next_event_ptr_arr_engine,  // Use engine-specific array
            single_seq_marginals.marginal_array_smart_p, // Accumulator
            current_params_marginals.marginal_array_smart_p, // Parameters (Read-only)
            alignments,
            constructed_sequences,
            seq_offsets,
            err_rate,
            counters_dummy,
            events_map,
            safety_set,
            mismatches_lists,
            max_proba_scenario,
            proba_threshold_factor
        );
    } catch (const std::exception& e) {
        // In validation, we might want to rethrow or just log
        // For now, rethrow to fail test
        throw;
    }
    // 3. Normalize error rate likelihoods?
    // GenModel does: err_rate->norm_weights_by_seq_likelihood(...)
    // This normalizes the added marginals by P(data).
    // The Engine expects Sum P(scenario|data).
    // iterate accumulates P(scenario, data).
    // So yes, we must normalize by P(data) = seq_likelihood.

    double seq_likelihood = err_rate->get_seq_likelihood();
    if (seq_likelihood > 0) {
        long double inv_lik = 1.0 / seq_likelihood;
        size_t size = single_seq_marginals.get_length(); // or arr_size
        for(size_t i=0; i<size; ++i) {
             single_seq_marginals.marginal_array_smart_p[i] *= inv_lik;
        }
    }

    // 4. Distribute to Engine
    accumulate_from_legacy_flat(engine, single_seq_marginals, parms);

    return seq_likelihood;
}

/**
 * @brief Run Parallel EM Validation (2 iterations)
 *
 * Replicates the E-M loop from GenModel::infer_model() but processes each sequence through
 * both the legacy Counter system and the new Parallel Recursion system.
 */
inline void infer_model_parallel_validation(
    const std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>& sequences,
    Model_Parms& model_parms,
    Model_marginals& legacy_marginals,
    igor::model::InferenceEngine<long double>& new_engine,
    int n_iterations,
    double mean_number_seq_err_thresh = INFINITY
) {
    using namespace std;
    using namespace igor::model;
    using Catch::Matchers::WithinAbs;

    // 1. Setup Logic
    std::queue<std::shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    std::unordered_map<Rec_Event_name, int> index_map = legacy_marginals.get_index_map(model_parms, model_queue);
    auto offsets_map = legacy_marginals.get_offsets_map(model_parms, model_queue);
    auto events_map = model_parms.get_events_map();

    // Copy legacy marginals to keep state tracking (shared read-only parameters for both paths)
    Model_marginals shared_parameters = legacy_marginals;

    // Build Index_map once (reused for all sequences)
    auto events_list = model_parms.get_event_list();
    Index_map shared_index_mapp(events_list.size());
    for (const auto& event : events_list) {
        int id = event->get_event_identifier();
        shared_index_mapp.request_memory_layer(id);
        shared_index_mapp.set_value(id, index_map.at(event->get_name()), 0);

        // Initialize event sizes and bounds (Crucial for iterate)
        size_t event_size = legacy_marginals.get_event_size(event, model_parms);
        event->set_event_marginal_size(event_size);
        // We skip crude upper bound for validation if not strictly needed, but GenModel does it.
        event->set_crude_upper_bound_proba(index_map.at(event->get_name()), event_size,
                                           legacy_marginals.marginal_array_smart_p);
    }

    // Initialize Rec_Event internal structures (GenModel::initialize_event)
    {
        std::unordered_set<Rec_Event_name> init_processed_events;
        Safety_bool_map safety_set(3);
        Seq_type_str_p_map constructed_sequences(6);
        Mismatch_vectors_map mismatches_lists(6);
        Seq_offsets_map seq_offsets(6, 3);
        Downstream_scenario_proba_bound_map downstream_proba_map(6);
        downstream_proba_map.init_first_layer(1.0);

        std::shared_ptr<Error_rate> tmp_err_rate = model_parms.get_err_rate_p();
        tmp_err_rate->initialize(events_map);
        tmp_err_rate->set_viterbi_run(false);

        // Initialize each event
        std::stack<std::shared_ptr<Rec_Event>> init_stack;
        std::queue<std::shared_ptr<Rec_Event>> init_q = model_queue;
        while(!init_q.empty()){
            auto event = init_q.front();
            init_q.pop();
            init_stack.push(event);
            event->initialize_event(init_processed_events, events_map, offsets_map,
                downstream_proba_map, constructed_sequences, safety_set,
                tmp_err_rate, mismatches_lists, seq_offsets, shared_index_mapp);
            event->set_viterbi_run(false);
        }

        // Initialize probability bounds (CRITICAL - GenModel does this)
        double downstream_proba_bound = 1.0;
        std::forward_list<double*> updated_proba_list;
        while(!init_stack.empty()) {
            auto event = init_stack.top();
            init_stack.pop();

            std::queue<std::shared_ptr<Rec_Event>> tmp_q = model_queue;
            while(tmp_q.front() != event) {
                tmp_q.pop();
            }
            tmp_q.pop();

            event->initialize_crude_scenario_proba_bound(downstream_proba_bound, updated_proba_list, events_map);
            event->initialize_Len_proba_bound(tmp_q, legacy_marginals.marginal_array_smart_p, shared_index_mapp);
        }

        // Update event internal probabilities (CRITICAL - GenModel does this)
        init_q = model_queue;
        while(!init_q.empty()) {
            init_q.front()->update_event_internal_probas(legacy_marginals.marginal_array_smart_p, index_map);
            init_q.pop();
        }
    }

    // 2. Loop Iterations
    for(int iter = 0; iter < n_iterations; ++iter) {

        // --- Legacy Setup ---
        Model_marginals new_legacy_marginals(model_parms);

        auto err_p = model_parms.get_err_rate_p();
        if(!err_p) std::cout << "ERROR: Error rate pointer is null!" << std::endl;
        std::shared_ptr<Error_rate> new_legacy_err_rate = err_p->copy();

        new_legacy_err_rate->initialize(events_map);

        // --- New Engine Setup ---
        new_engine.reset_accumulators();

        std::shared_ptr<Error_rate> new_engine_err_rate = model_parms.get_err_rate_p()->copy();
        new_engine_err_rate->initialize(events_map);

        // Debug info
        // cout << "Starting Iteration " << iter + 1 << endl;

        // 3. Loop Sequences
        // 3. Loop Sequences
        for (size_t i = 0; i < sequences.size(); ++i) {
            const auto& seq_tuple = sequences[i];

            // --- Common Prep ---
            // Prepare Next_event_ptr array (replicates GenModel init)
            // Note: Next_event_ptr is Rec_Event* (raw pointer to events owned by model_queue)
            // We allocate array of pointers but DON'T delete the pointed-to objects
            size_t num_events = model_parms.get_event_list().size();
            Next_event_ptr* raw_arr = new Next_event_ptr[num_events];
            for(size_t k = 0; k < num_events; ++k) raw_arr[k] = nullptr;  // Initialize all to nullptr
            shared_ptr<Next_event_ptr> next_event_ptr_arr(
                raw_arr,
                [](Next_event_ptr* p) { delete[] p; }  // Only delete the array, not the pointed-to objects
            );

            queue<shared_ptr<Rec_Event>> local_queue = model_queue;
            while (!local_queue.empty()) {
                shared_ptr<Rec_Event> event = local_queue.front();
                local_queue.pop();
                if (!local_queue.empty()) {
                    int id = event->get_event_identifier();
                    if(id < 0 || id >= (int)model_parms.get_event_list().size()) {
                         std::cout << "FATAL: ID out of bounds: " << id << std::endl;
                    }
                    next_event_ptr_arr.get()[id] = local_queue.front().get();
                } else {
                    next_event_ptr_arr.get()[event->get_event_identifier()] = nullptr;
                }
            }

            // --- LEGACY PATH ---
            // GenModel calls first_event->iterate(...) manually
            // Note: Keep these structures in scope until after engine path completes
            // to avoid memory corruption during early cleanup
            std::cout << "Seq " << i << ": Starting legacy iterate()..." << std::flush;
            queue<shared_ptr<Rec_Event>> legacy_q = model_queue;
            shared_ptr<Rec_Event> first_event = legacy_q.front();

            // Create fresh Index_map for this sequence (Index_map is stateful and gets modified during iterate)
            Index_map seq_index_map_legacy(events_list.size());
            for (const auto& event : events_list) {
                int id = event->get_event_identifier();
                seq_index_map_legacy.request_memory_layer(id);
                seq_index_map_legacy.set_value(id, index_map.at(event->get_name()), 0);
            }

            // Initialize single seq marginals using empty_copy() like GenModel does
            Model_marginals single_seq_marginals = shared_parameters.empty_copy();

            // Dummies
            map<size_t, shared_ptr<Counter>> dummy_counters;

            Seq_type_str_p_map constructed_sequences(6);

            Mismatch_vectors_map mismatches_lists(6);

            Seq_offsets_map seq_offsets(6, 3);

            Safety_bool_map safety_set(3);

            Downstream_scenario_proba_bound_map downstream_map(6);

            downstream_map.init_first_layer(1.0);

            double init_proba = 1.0;
            double max_proba_scenario = 1e-25 / 0.001;
            double proba_threshold_factor = 0.001;

            try {
                first_event->iterate(
                    init_proba, downstream_map, get<1>(seq_tuple), nt2int(get<1>(seq_tuple)),
                    seq_index_map_legacy, offsets_map, next_event_ptr_arr,
                    single_seq_marginals.marginal_array_smart_p,
                    shared_parameters.marginal_array_smart_p, // Use shared params
                    get<2>(seq_tuple), constructed_sequences, seq_offsets,
                    new_legacy_err_rate, dummy_counters, events_map, safety_set,
                    mismatches_lists, max_proba_scenario, proba_threshold_factor
                );
            } catch(const std::exception& e) {
                std::cerr << "\nException in iterate() for seq " << i << ": " << e.what() << std::endl;
                std::cerr << "Sequence: " << get<1>(seq_tuple) << std::endl;
                std::cerr << "Number of events in events_map: " << events_map.size() << std::endl;
                for (const auto& pair : events_map) {
                    std::cerr << "  Event: " << pair.second->get_nickname()
                              << " (type: " << pair.second->get_type() << ")" << std::endl;
                }
                std::cerr << "Alignments: " << std::endl;
                for (const auto& pair : get<2>(seq_tuple)) {
                    std::cerr << "  Gene " << pair.first << ": " << pair.second.size() << " alignments" << std::endl;
                }
                throw;
            }

            std::cout << " done" << std::endl;

            // Normalize and Add Legacy
            new_legacy_err_rate->norm_weights_by_seq_likelihood(
                single_seq_marginals.marginal_array_smart_p, single_seq_marginals.get_length()
            );

            double seq_lik = new_legacy_err_rate->get_seq_likelihood();
            std::cout << "Seq " << i << " legacy: seq_lik=" << seq_lik
                      << " mean_err=" << new_legacy_err_rate->get_seq_mean_error_number()
                      << " single_seq[0-4]=" << single_seq_marginals.marginal_array_smart_p[0]
                      << "," << single_seq_marginals.marginal_array_smart_p[1]
                      << "," << single_seq_marginals.marginal_array_smart_p[2]
                      << "," << single_seq_marginals.marginal_array_smart_p[3]
                      << "," << single_seq_marginals.marginal_array_smart_p[4] << std::endl;

            // Skip sequences with zero likelihood (no valid scenarios found)
            if(seq_lik <= 0) {
                std::cout << "  -> SKIPPED (zero likelihood)" << std::endl;
                continue;
            }

            if(new_legacy_err_rate->get_seq_mean_error_number() <= mean_number_seq_err_thresh) {
                new_legacy_err_rate->add_to_norm_counter();
                new_legacy_marginals += single_seq_marginals;
                std::cout << "  -> Added to legacy accumulator, new_legacy[0]=" << new_legacy_marginals.marginal_array_smart_p[0] << std::endl;
            } else {
                std::cout << "  -> SKIPPED (threshold exceeded)" << std::endl;
                continue;
            }

            // --- NEW PATH ---
            // count_scenario_to_engine handles normalization internally
            // Note: count_scenario creates its own separate next_event_ptr array
            std::cout << "Seq " << i << ": Starting engine count_scenario..." << std::flush;
            double lik_new = count_scenario_to_engine(
                seq_tuple,
                new_engine,
                model_parms,
                shared_parameters,
                events_map,
                index_map,
                offsets_map,
                new_engine_err_rate
            );
            std::cout << " done" << std::endl;
            std::cout << "Seq " << i << " engine: lik=" << lik_new << std::endl;

            // Apply same threshold filter as legacy path
            // Note: count_scenario_to_engine already accumulated to engine.
            // If we exceed threshold, we need to subtract it back.
            // Better approach: conditionally call count_scenario_to_engine
            // For validation, threshold should be INFINITY (disabled) to ensure identical behavior.
            // If threshold is finite, legacy may skip sequences that engine already processed.
            // To match legacy exactly when threshold is finite:
            if(new_engine_err_rate->get_seq_mean_error_number() > mean_number_seq_err_thresh) {
                // Legacy didn't accumulate this sequence - we shouldn't either
                // But we already did in count_scenario_to_engine!
                // Hacky fix: reset engine accumulators and rebuild without this sequence
                // Better: make count_scenario_to_engine return Model_marginals, add conditionally
                // For now: REQUIRE threshold is INFINITY for validation
                if(mean_number_seq_err_thresh != INFINITY) {
                    throw std::runtime_error("Finite error threshold not supported in validation - use INFINITY");
                }
            }
            new_engine_err_rate->add_to_norm_counter();

            // --- ASSERT PER SEQUENCE ---
            double lik_legacy = new_legacy_err_rate->get_seq_likelihood();
            REQUIRE_THAT(lik_new, WithinAbs(lik_legacy, 1e-12));
        }

        // 4. M-Step & Comparison

        // DEBUG: Check RAW accumulators BEFORE normalization
        std::cout << "=== BEFORE M-STEP (Iteration " << (iter + 1) << ") ===" << std::endl;
        std::cout << "Legacy raw accumulators (first 10): ";
        for(int k = 0; k < 10 && k < (int)new_legacy_marginals.get_length(); ++k) {
            std::cout << new_legacy_marginals.marginal_array_smart_p[k] << " ";
        }
        std::cout << std::endl;

        // Export engine accumulators to legacy format (WITHOUT normalizing)
        Model_marginals engine_raw_accumulators(model_parms);
        for(size_t k = 0; k < engine_raw_accumulators.get_length(); ++k) {
            engine_raw_accumulators.marginal_array_smart_p[k] = 0.0;
        }

        // Manually extract raw accumulators from engine
        auto events_list_debug = model_parms.get_event_list();
        for (const auto& event : events_list_debug) {
            std::string nickname = event->get_nickname();
            if (!new_engine.has_handler(nickname)) continue;

            auto& handler = new_engine.handler(nickname);
            auto& acc = handler.accumulator();

            std::string full_name = event->get_name();
            if (index_map.find(full_name) == index_map.end()) continue;
            int base_index = index_map.at(full_name);

            for (size_t i = 0; i < acc.size(); ++i) {
                engine_raw_accumulators.marginal_array_smart_p[base_index + i] = static_cast<long double>(acc.data()[i]);
            }
        }

        std::cout << "Engine raw accumulators (first 10): ";
        for(int k = 0; k < 10 && k < (int)engine_raw_accumulators.get_length(); ++k) {
            std::cout << engine_raw_accumulators.marginal_array_smart_p[k] << " ";
        }
        std::cout << std::endl;

        // Compare raw accumulators
        double max_raw_diff = 0.0;
        for(size_t k = 0; k < engine_raw_accumulators.get_length(); ++k) {
            double diff = std::abs(new_legacy_marginals.marginal_array_smart_p[k] - engine_raw_accumulators.marginal_array_smart_p[k]);
            max_raw_diff = std::max(max_raw_diff, diff);
        }
        std::cout << "Max difference in RAW accumulators: " << max_raw_diff << std::endl;

        // Normalize Legacy
        // Note: get_inverse_offset_map returns list, needed for normalize
        auto inv_offset_map = new_legacy_marginals.get_inverse_offset_map(model_parms, model_queue);
        new_legacy_marginals.normalize(inv_offset_map, index_map, model_queue);
        new_legacy_marginals.copy_fixed_events_marginals(shared_parameters, model_parms, index_map);

        // Maximize New Engine
        new_engine.update_parameters();

        // Update Error Rates (as GenModel does)
        new_legacy_err_rate->update();
        new_engine_err_rate->update();
        model_parms.set_error_ratep(new_legacy_err_rate);

        // Compare Parameters
        // Convert New Engine to Legacy format for easy comparison (and handle fixed events)
        Model_marginals new_engine_as_legacy(model_parms);
        igor::core::export_to_legacy(new_engine, new_engine_as_legacy, model_parms);

        // Handle Fixed Events for Engine Result (mimic Legacy behavior)
        new_engine_as_legacy.copy_fixed_events_marginals(shared_parameters, model_parms, index_map);

        // diff = new_legacy_marginals - new_engine_as_legacy
        Model_marginals diff = new_legacy_marginals - new_engine_as_legacy;
        size_t arr_size = diff.get_length();
        double max_diff = 0.0;
        for(size_t i=0; i<arr_size; ++i) {
            max_diff = std::max(max_diff, (double)std::abs(diff.marginal_array_smart_p[i]));
        }

        if(max_diff >= 1e-12) {
            std::cout << "MISMATCH at iteration " << (iter + 1) << ", max_diff = " << max_diff << std::endl;
        }
        REQUIRE(max_diff < 1e-12);

        // Update shared parameters for next iteration
        shared_parameters = new_legacy_marginals;
        legacy_marginals = new_legacy_marginals; // Update output arg

        // Update crude upper bounds for next iteration (GenModel does this)
        for (const auto& event : events_list) {
            event->set_crude_upper_bound_proba(index_map.at(event->get_name()),
                                              legacy_marginals.get_event_size(event, model_parms),
                                              shared_parameters.marginal_array_smart_p);
        }
    }
}
