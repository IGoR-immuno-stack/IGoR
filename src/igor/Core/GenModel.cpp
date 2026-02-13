/*
 * GenModel.cpp
 *
 *  Created on: 3 nov. 2014
 *      Author: Quentin Marcou
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 analyze and model immune receptors
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

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *      This class designs a generative model and supply all the methods to run
 a maximum likelihood estimate of the generative model
 */

#include <igor/Core/GenModel.h>
#include <igor/Core/SequenceTypes.h>
#include <igor/Core/DynamicSequenceMap.h>

using namespace std;

GenModel::GenModel(const Model_Parms &parms, const Model_marginals &marginals,
                   const map<size_t, shared_ptr<Counter>> &count_list)
    : model_parms(parms), model_marginals(marginals), counters_list(count_list)
{
}

GenModel::GenModel(const Model_Parms &parms, const Model_marginals &marginals)
    : GenModel(parms, marginals, map<size_t, shared_ptr<Counter>>())
{
}

GenModel::GenModel(const Model_Parms &parms)
    : GenModel(parms, Model_marginals(parms), map<size_t, shared_ptr<Counter>>())
{
}

GenModel::~GenModel()
{
    // TODO Auto-generated destructor stub
}

bool GenModel::infer_model(
        const vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> &sequences,
        const int iterations, const std::string path, bool fast_iter, double likelihood_threshold /*=1e-25*/,
        bool viterbi_like /*false*/)
{
    return this->infer_model(sequences, iterations, path, fast_iter, likelihood_threshold, viterbi_like, 0.001,
                             INFINITY);
}

bool GenModel::infer_model(
        const vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> &sequences,
        const int iterations, const std::string path, bool fast_iter, double likelihood_threshold,
        double proba_threshold_factor)
{
    return this->infer_model(sequences, iterations, path, fast_iter, likelihood_threshold, false,
                             proba_threshold_factor, INFINITY);
}

bool GenModel::infer_model(
        const vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> &sequences,
        const int iterations, const string path, bool fast_iter /*=true*/,
        double likelihood_threshold /*=1e-25 by default*/, bool viterbi_like /*=false*/,
        double proba_threshold_factor /*=0.001 by default*/,
        double mean_number_seq_err_thresh /*= INFINITY by default*/)
{

    // If viterbi like only the best scenario is of interest
    if (viterbi_like) {
        cerr << "******************************************************************" << endl;
        cerr << "*\t\t RUNNING \"VITERBI\" LIKE ALGORITHM \t\t *" << endl
             << "* \t(only the best scenario will be taken into account)\t *" << endl;
        cerr << "******************************************************************" << endl;
        proba_threshold_factor = 1.0;
    }

    if (likelihood_threshold > 1.0) {
        // throw invalid_argument("Likelihood threshold must be lesser or equal than one");
    }

    if (proba_threshold_factor > 1.0) {
        throw invalid_argument("Probability threshold ratio must be lesser or equal than one");
    }

    ofstream likelihood_file(path + "likelihoods.out");
    likelihood_file << "iteration;mean_log_Likelihood;n_seq" << endl;

    queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map = model_marginals.get_index_map(model_parms, model_queue);

    unordered_map<Rec_Event_name, list<pair<shared_ptr<const Rec_Event>, int>>> inv_offset_map =
            model_marginals.get_inverse_offset_map(model_parms, model_queue);
    int iteration_accomplished = 0;
    ofstream log_file(path + string("inference_logs.txt"));
    log_file << "iteration_n;seq_processed;seq_index;nt_sequence;n_V_aligns;n_J_"
                "aligns;seq_likelihood;seq_mean_n_errors;seq_n_scenarios;seq_"
                "best_scenario;time"
             << endl;
    ofstream general_logs(path + string("inference_info.out"));
    // Dump all inference parameters to file
    chrono::system_clock::time_point begin_time = chrono::system_clock::now();
    std::time_t tt;
    tt = chrono::system_clock::to_time_t(begin_time);

    general_logs << "Date: " << ctime(&tt) << endl;
    general_logs << "Max #iterations to be performed: " << iterations << endl;
    general_logs << "Path: " << path << endl;
    general_logs << "First iter fast(only best V and best J considered): " << fast_iter << endl;
    general_logs << "Min Likelihood threshold: " << likelihood_threshold << endl;
    general_logs << "Viterbi like (only keeps the best scenario): " << viterbi_like << endl;
    general_logs << "Proba threshold ratio: " << proba_threshold_factor
                 << "\t#(ratio between best scenario and current scenario needed "
                    "to explore/count the scenario)"
                 << endl;
    general_logs << "Mean #errors threshold: " << mean_number_seq_err_thresh
                 << "\t#Needs a very good reason to be set to another value than INFINITY" << endl;

    // Get the total number of sequences to process
    const double total_number_seqs = sequences.size(); // Use a double for float division afterwards

    /*
   * Get the list of fixed and inferred events and output them to the log file
   * Do it in a scope so the variables will be destroyed
   */
    {
        list<Rec_Event_name> fixed_events_list;
        list<Rec_Event_name> inferred_events_list;
        const list<shared_ptr<Rec_Event>> model_event_list = model_parms.get_event_list();
        for (list<shared_ptr<Rec_Event>>::const_iterator iter = model_event_list.begin();
             iter != model_event_list.end(); ++iter) {
            if (not(*iter)->is_fixed()) {
                inferred_events_list.emplace_back((*iter)->get_name());
            } else {
                fixed_events_list.emplace_back((*iter)->get_name());
            }
        }
        general_logs << endl;
        general_logs << "List of updated events: ";
        for (Rec_Event_name name : inferred_events_list) {
            general_logs << name << "\t";
        }
        general_logs << endl;
        general_logs << "List of fixed events: ";
        for (Rec_Event_name name : fixed_events_list) {
            general_logs << name << "\t";
        }
        general_logs << endl;
        general_logs << "Error model updated: " << model_parms.get_err_rate_p()->is_updated() << endl;
    }

    // Write initial condition to file
    this->model_marginals.write2txt(path + string("initial_marginals.txt"), this->model_parms);
    this->model_parms.write_model_parms(path + string("initial_model.txt"));

    /*
   * First initialization creates file streams
   * This is to make sure that the counter copies do not create new files each
   * time
   */
    for (map<size_t, shared_ptr<Counter>>::const_iterator iter = counters_list.begin(); iter != counters_list.end();
         ++iter) {
        (*iter).second->initialize_counter(model_parms, model_marginals);
    }

    // Loop over iterations
    while (iteration_accomplished != iterations) {

        Model_marginals new_marginals = Model_marginals(model_parms);
        shared_ptr<Error_rate> error_rate_copy = model_parms.get_err_rate_p()->copy();

        // Initialize error rate copy
        error_rate_copy->initialize(model_parms.get_events_map());

        // Initialize counters for the log file
        size_t sequences_processed = 0;

        const vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> *sequence_util_ptr;

        // Take only best alignments if fast_iter
        vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> fast_iter_sequences;
        if (fast_iter && iteration_accomplished == 0 && !sequences.empty()) {
            fast_iter_sequences = sequences;
            for (unordered_map<int, vector<Alignment_data>>::const_iterator gc_align_iter =
                         std::get<2>(sequences.at(0)).begin();
                 gc_align_iter != std::get<2>(sequences.at(0)).end(); ++gc_align_iter) {
                if ((*gc_align_iter).first == (int)D_gene)
                    continue;
                fast_iter_sequences = get_best_aligns(fast_iter_sequences, (*gc_align_iter).first);
            }
            sequence_util_ptr = &fast_iter_sequences;
        } else {
            sequence_util_ptr = &sequences;
        }

        cerr << "Performing Evaluate/Inference iteration " << iteration_accomplished + 1 << endl;

#pragma omp parallel shared(new_marginals, error_rate_copy, sequences_processed, sequence_util_ptr, sequences) \
        firstprivate(model_queue, proba_threshold_factor)
        {
            // Make single thread copies of objects for thread safety
            Model_Parms single_thread_model_parms(model_parms);
            queue<shared_ptr<Rec_Event>> single_thread_model_queue = single_thread_model_parms.get_model_queue();
            unordered_map<Rec_Event_name, int> single_thread_index_map =
                    model_marginals.get_index_map(single_thread_model_parms, single_thread_model_queue);
            Model_marginals single_thread_model_marginals(model_marginals);
            shared_ptr<Error_rate> single_thread_err_rate = single_thread_model_parms.get_err_rate_p();
            unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> events_map =
                    single_thread_model_parms.get_events_map();
            map<size_t, shared_ptr<Counter>> single_thread_counter_list;
            for (map<size_t, shared_ptr<Counter>>::const_iterator iter = this->counters_list.begin();
                 iter != this->counters_list.end(); ++iter) {
                // Copy only relevant counters for this iteration
                if ((not(*iter).second->is_last_iter_only()) or (iteration_accomplished == iterations - 1)) {
                    single_thread_counter_list.emplace((*iter).first, (*iter).second->copy());
                }
            }

            unordered_set<Rec_Event_name> init_processed_events;

            // Initialize Enum_fast_memory map and dual maps
            Safety_bool_map safety_set(SequenceTypeRegistry::get_instance().size());
            Seq_type_str_p_map constructed_sequences(SequenceTypeRegistry::get_instance().size()); // Dynamic size
            Mismatch_vectors_map mismatches_lists(SequenceTypeRegistry::get_instance().size());
            Seq_offsets_map seq_offsets(SequenceTypeRegistry::get_instance().size(), 3);

            // Initialize downstream probas to 1
            Downstream_scenario_proba_bound_map downstream_proba_map(SequenceTypeRegistry::get_instance().size());
            downstream_proba_map.init_first_layer(1.0);

            list<shared_ptr<Rec_Event>> events_list = single_thread_model_parms.get_event_list();
            Index_map index_mapp(events_list.size());

            // Initialize index_map
            for (list<shared_ptr<Rec_Event>>::iterator event_iter = events_list.begin();
                 event_iter != events_list.end(); ++event_iter) {
                int event_index = (*event_iter)->get_event_identifier();
                index_mapp.request_memory_layer(event_index);

                auto it = single_thread_index_map.find((*event_iter)->get_name());
                if (it == single_thread_index_map.end()) {
                    cerr << "!!! GenModel.cpp: Key not found in single_thread_index_map: " << (*event_iter)->get_name()
                         << " !!!" << endl;
                    cerr << "Available keys:" << endl;
                    for (auto const &[k, v] : single_thread_index_map)
                        cerr << "  " << k << endl;
                    throw out_of_range("Key not found in single_thread_index_map");
                }
                index_mapp.set_value(event_index, it->second, 0);

                // Get events probability upper bounds
                size_t event_size =
                        single_thread_model_marginals.get_event_size((*event_iter), single_thread_model_parms);

                (*event_iter)->set_event_marginal_size(event_size);
                (*event_iter)
                        ->set_crude_upper_bound_proba(it->second, event_size,
                                                      single_thread_model_marginals.marginal_array_smart_p);
            }

            queue<shared_ptr<Rec_Event>> init_single_thread_model_queue = single_thread_model_queue;
            unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> single_thread_offset_map =
                    model_marginals.get_offsets_map(single_thread_model_parms, single_thread_model_queue);

            stack<shared_ptr<Rec_Event>> init_single_thread_stack;

            // Initialize events
            while (!init_single_thread_model_queue.empty()) {
                shared_ptr<Rec_Event> first_init_event = init_single_thread_model_queue.front();
                init_single_thread_stack.push(first_init_event);
                init_single_thread_model_queue.pop();
                (*first_init_event)
                        .initialize_event(init_processed_events, events_map, single_thread_offset_map,
                                          downstream_proba_map, constructed_sequences, safety_set,
                                          single_thread_err_rate, mismatches_lists, seq_offsets, index_mapp);
                (*first_init_event).set_viterbi_run(viterbi_like);
            }

            shared_ptr<Next_event_ptr> next_event_ptr_arr(
                    new Next_event_ptr[single_thread_model_parms.get_event_list().size()]);
            init_single_thread_model_queue = single_thread_model_queue;
            while (!init_single_thread_model_queue.empty()) {
                shared_ptr<Rec_Event> first_init_event = init_single_thread_model_queue.front();
                init_single_thread_model_queue.pop();
                if (!init_single_thread_model_queue.empty()) {
                    next_event_ptr_arr.get()[first_init_event->get_event_identifier()] =
                            init_single_thread_model_queue.front().get();
                } else {
                    next_event_ptr_arr.get()[first_init_event->get_event_identifier()] = NULL;
                }
            }

            // Initialize error rate
            single_thread_err_rate->initialize(events_map);
            single_thread_err_rate->set_viterbi_run(viterbi_like);

            // Initialize Counters
            for (map<size_t, shared_ptr<Counter>>::iterator iter = single_thread_counter_list.begin();
                 iter != single_thread_counter_list.end(); ++iter) {
                (*iter).second->initialize_counter(single_thread_model_parms, single_thread_model_marginals);
                }
#pragma omp single nowait
            {
                cerr << "Initializing probability bounds..." << endl;
            }
            // Compute upper proba bounds for downstream scenarios for each event
            double downstream_proba_bound = 1;
            forward_list<double *> updated_proba_list;
            while (!init_single_thread_stack.empty()) {
                shared_ptr<Rec_Event> last_proba_init_event = init_single_thread_stack.top();
                queue<shared_ptr<Rec_Event>> tmp_init_proba_single_thread_model_queue = single_thread_model_queue;
                init_single_thread_stack.pop();
                while (tmp_init_proba_single_thread_model_queue.front() != last_proba_init_event) {
                    tmp_init_proba_single_thread_model_queue.pop();
                }
                tmp_init_proba_single_thread_model_queue.pop();
                last_proba_init_event->initialize_crude_scenario_proba_bound(downstream_proba_bound, updated_proba_list,
                                                                             events_map);

                last_proba_init_event->initialize_Len_proba_bound(tmp_init_proba_single_thread_model_queue,
                                                                  single_thread_model_marginals.marginal_array_smart_p,
                                                                  index_mapp);
            }
#pragma omp single nowait
            {
                cerr << "Initialization of probability bounds over." << endl;
            }

            // Now let all the events in the need of it get their own updated copy of
            // the marginals
            init_single_thread_model_queue = single_thread_model_queue;
            while (!init_single_thread_model_queue.empty()) {
                init_single_thread_model_queue.front()->update_event_internal_probas(
                        single_thread_model_marginals.marginal_array_smart_p, single_thread_index_map);
                init_single_thread_model_queue.pop();
            }

            chrono::system_clock::time_point single_seq_begin;
            chrono::duration<double> seq_time;

#pragma omp for schedule(dynamic) nowait
            for (auto i = 0; i < (*sequence_util_ptr).size(); ++i) {
                const auto &seq_it = (*sequence_util_ptr).at(i);

                single_seq_begin = chrono::system_clock::now();

                // Make a copy of the queue that can be modified in iterate
                queue<shared_ptr<Rec_Event>> model_queue_copy(single_thread_model_queue);

                // Get the first event from the queue
                shared_ptr<Rec_Event> first_event = model_queue_copy.front();
                model_queue_copy.pop();

                // Initialize single seq marginals
                Model_marginals single_seq_marginals = single_thread_model_marginals.empty_copy();
                double init_proba = 1;
                double max_proba_scenario = likelihood_threshold / proba_threshold_factor;

                Int_Str int_sequence = nt2int(get<1>(seq_it));

                try {

                    first_event->iterate(
                            init_proba, downstream_proba_map, get<1>(seq_it), int_sequence, index_mapp,
                            single_thread_offset_map, next_event_ptr_arr, single_seq_marginals.marginal_array_smart_p,
                            single_thread_model_marginals.marginal_array_smart_p, get<2>(seq_it), constructed_sequences,
                            seq_offsets, single_thread_err_rate, single_thread_counter_list, events_map, safety_set,
                            mismatches_lists, max_proba_scenario, proba_threshold_factor);

                } catch (exception &e) {
#pragma omp critical
                    {
                        cerr << "Exception caught during scenario exploration: " << e.what() << endl;
                        cerr << "Sequence index: " << get<0>(seq_it) << endl;
                        cerr << "Sequence: " << get<1>(seq_it) << endl;
                    }
                }

                single_seq_marginals.normalize(inv_offset_map, single_thread_index_map, single_thread_model_queue);
                single_seq_marginals.copy_fixed_events_marginals(single_thread_model_marginals,
                                                                 single_thread_model_parms, single_thread_index_map);

                single_thread_err_rate->norm_weights_by_seq_likelihood(single_seq_marginals.marginal_array_smart_p,
                                                                       single_seq_marginals.get_length());
                seq_time = chrono::system_clock::now() - single_seq_begin;
#pragma omp critical(dump_seq_info)
                {
                    ++sequences_processed;
                    log_file << iteration_accomplished << ";" << sequences_processed << ";" << get<0>(seq_it) << ";"
                             << get<1>(seq_it) << ";" << get<2>(seq_it).at(V_gene).size() << ";" << get<2>(seq_it).at(J_gene).size() << ";"
                             << single_thread_err_rate->get_seq_likelihood() << ";"
                             << single_thread_err_rate->get_seq_mean_error_number() << ";"
                             << single_thread_err_rate->debug_number_scenarios << ";" << max_proba_scenario << ";"
                             << seq_time.count() << endl;
                }
                for (map<size_t, shared_ptr<Counter>>::iterator iter = single_thread_counter_list.begin();
                     iter != single_thread_counter_list.end(); ++iter) {
                    iter->second->count_sequence(single_thread_err_rate->get_seq_likelihood(), single_seq_marginals,
                                                 single_thread_model_parms);
#pragma omp critical(dump_counters)
                    {
                        (*iter).second->dump_sequence_data(get<0>(seq_it), iteration_accomplished);
                    }
                }

                if (single_thread_err_rate->get_seq_mean_error_number() <= mean_number_seq_err_thresh) {
                    // Add weighed errors to the normalized error counter
                    single_thread_err_rate->add_to_norm_counter();

                    // Add the single_seq_marginals to the single thread marginals
                    single_thread_model_marginals += single_seq_marginals;
                } else {
                    // Erase seq specific counters so that it won't contribute to the
                    // error rate
                    single_thread_err_rate->clean_seq_counters();
                }

#pragma omp critical(update_progress_bar)
                {
                    if (sequences_processed % 100 == 0) {
                        // Output current progress to cerr
                        show_progress_bar(cerr, sequences_processed / total_number_seqs,
                                          "Iteration " + to_string(iteration_accomplished + 1), 50);
                    }
                }
            }

// Merge single thread error_rates and marginals
#pragma omp critical(merge_marginals_and_er)
            {
                new_marginals += single_thread_model_marginals;
                add_to_err_rate(error_rate_copy.get(), single_thread_err_rate.get());
                for (map<size_t, shared_ptr<Counter>>::iterator iter = single_thread_counter_list.begin();
                     iter != single_thread_counter_list.end(); ++iter) {
                    counters_list.at((*iter).first)->add_to_counter((*iter).second);
                }
            }
        }

        for (map<size_t, shared_ptr<Counter>>::const_iterator iter = counters_list.begin(); iter != counters_list.end();
             ++iter) {
            (*iter).second->dump_data_summary(iteration_accomplished);
        }

        likelihood_file << iteration_accomplished + 1 << ";"
                        << error_rate_copy->get_model_likelihood()
                        / error_rate_copy->get_number_non_zero_likelihood_seqs()
                        << ";" << error_rate_copy->get_number_non_zero_likelihood_seqs() << endl;
        error_rate_copy->update();
        this->model_parms.set_error_ratep(error_rate_copy);
        new_marginals.normalize(inv_offset_map, index_map, model_queue);
        new_marginals.copy_fixed_events_marginals(this->model_marginals, this->model_parms, index_map);
        this->model_marginals = new_marginals;
        ++iteration_accomplished;

        this->model_marginals.write2txt(
                path + string("iteration_") + to_string(iteration_accomplished) + string(".txt"), this->model_parms);
        this->model_parms.write_model_parms(path + string("iteration_") + to_string(iteration_accomplished)
                                            + string("_parms.txt"));

        // Close current iteration progress bar
        close_progress_bar(cerr, "Iteration " + to_string(iteration_accomplished), 50);
    }
    // Create a copy of the last iteration results with identifiable name
    this->model_marginals.write2txt(path + string("final_marginals.txt"), this->model_parms);
    this->model_parms.write_model_parms(path + string("final_parms.txt"));

    return 0;
}

vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>>
get_best_aligns(const vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> &sequences,
                int gc)
{
    vector<tuple<int, string, unordered_map<int, vector<Alignment_data>>>> best_aligns = sequences;
    for (auto &seq : best_aligns) {
        auto &align_map = get<2>(seq);
        if (align_map.count(gc) > 0 && !align_map.at(gc).empty()) {
            auto &aligns = align_map.at(gc);
            auto best_it =
                    max_element(aligns.begin(), aligns.end(),
                                [](const Alignment_data &a, const Alignment_data &b) { return a.score < b.score; });
            Alignment_data best_align = *best_it;
            aligns.clear();
            aligns.push_back(best_align);
        }
    }
    return best_aligns;
}

std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>>
GenModel::generate_sequences(int n_seq, bool output_realizations)
{
    std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> sequences;
    std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());

    queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map = model_marginals.get_index_map(model_parms, model_queue);
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map =
            model_marginals.get_offsets_map(model_parms, model_queue);

    for (int i = 0; i < n_seq; ++i) {
        sequences.emplace_front(
                generate_unique_sequence(model_queue, index_map, offset_map, generator, output_realizations));
    }

    return sequences;
}

void GenModel::generate_sequences(int n_seq, bool output_realizations, string seq_file_path, string real_file_path,
                                  list<pair<gen_seq_trans, shared_ptr<void>>> transformations, bool output_only_func,
                                  int seed)
{
    std::mt19937_64 generator;
    if (seed == -1) {
        generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    } else {
        generator.seed(seed);
    }

    queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map = model_marginals.get_index_map(model_parms, model_queue);
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map =
            model_marginals.get_offsets_map(model_parms, model_queue);

    ofstream seq_file(seq_file_path);
    seq_file << "seq_index;nt_sequence" << endl;
    ofstream real_file;
    if (output_realizations) {
        real_file.open(real_file_path);
        real_file << "seq_index";
        queue<shared_ptr<Rec_Event>> tmp_queue = model_queue;
        while (!tmp_queue.empty()) {
            real_file << ";" << tmp_queue.front()->get_name();
            tmp_queue.pop();
        }
        real_file << ";Errors" << endl;
    }

    for (int i = 0; i < n_seq; ++i) {
        pair<string, queue<queue<int>>> seq_and_real =
                generate_unique_sequence(model_queue, index_map, offset_map, generator, output_realizations);
        seq_file << i << ";" << seq_and_real.first << endl;
        if (output_realizations) {
            real_file << i;
            queue<queue<int>> tmp_real = seq_and_real.second;
            while (!tmp_real.empty()) {
                real_file << ";(";
                queue<int> event_real = tmp_real.front();
                tmp_real.pop();
                while (!event_real.empty()) {
                    real_file << event_real.front();
                    event_real.pop();
                    if (!event_real.empty())
                        real_file << ",";
                }
                real_file << ")";
            }
            real_file << endl;
        }

        for (auto const &trans : transformations) {
            trans.first(i, seq_and_real, trans.second);
        }
    }
}

pair<string, queue<queue<int>>> GenModel::generate_unique_sequence(
        queue<shared_ptr<Rec_Event>> model_queue, unordered_map<Rec_Event_name, int> index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        mt19937_64 &generator, bool output_realizations)
{
    queue<queue<int>> realizations;
    unordered_map<int, string> constructed_sequences;

    while (!model_queue.empty()) {
        shared_ptr<Rec_Event> event = model_queue.front();
        model_queue.pop();
        queue<int> event_real = event->draw_random_realization(model_marginals.marginal_array_smart_p, index_map,
                                                               offset_map, constructed_sequences, generator);
        if (output_realizations) {
            realizations.push(event_real);
        }
    }

    // Build the final sequence from constructed_sequences starting from V
    string final_seq = "";
    SequenceTypeRegistry::TypeId current = SequenceTypeRegistry::V_GENE_SEQ;
    unordered_set<SequenceTypeRegistry::TypeId> visited;

    while (visited.find(current) == visited.end()) {
        visited.insert(current);
        if (constructed_sequences.count(current)) {
            final_seq += constructed_sequences.at(current);
        }

        auto downstream = SequenceTypeRegistry::get_instance().get_downstream_neighbors(current);
        if (downstream.empty())
            break;

        // Follow the first path (standard for linear recombination)
        // If there's a junction, add it too
        if (downstream[0].junction_type != (SequenceTypeRegistry::TypeId)-1) {
            if (constructed_sequences.count(downstream[0].junction_type)) {
                final_seq += constructed_sequences.at(downstream[0].junction_type);
            }
        }
        current = downstream[0].neighbor_type;
    }

    return make_pair(final_seq, realizations);
}

bool GenModel::write2txt()
{
    this->model_marginals.write2txt("final_marginals.txt", this->model_parms);
    this->model_parms.write_model_parms("final_model.txt");
    return true;
}

bool GenModel::readtxt()
{
    // TODO
    return true;
}

void GenModel::write_seq2txt(std::string path, std::forward_list<std::string> sequences)
{
    // Dummy
}

void GenModel::write_seq_real2txt(std::string path, std::string batch_name,
                                  std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> sequences)
{
    // Dummy
}

bool GenModel::load_genmodel()
{
    return true;
}

void output_CDR3_gen_data(size_t i, std::pair<std::string, std::queue<std::queue<int>>> seq_and_real,
                          std::shared_ptr<void> func_data)
{
    // Dummy
}

//==============================================================================
// Fast Sequence Generation Implementation
//==============================================================================

void GenModel::generate_sequences_fast(size_t num_sequences, const string &seq_filename, const string &real_filename,
                                       size_t num_threads, int64_t seed, bool show_progress)
{

    // Get folder path for generation info
    string folder_path = seq_filename.substr(0, seq_filename.rfind("/") + 1);
    ofstream generation_infos_file(folder_path + "generation_info.out", fstream::out | fstream::app);

    // Log generation parameters
    chrono::system_clock::time_point begin_time = chrono::system_clock::now();
    std::time_t tt = chrono::system_clock::to_time_t(begin_time);

    generation_infos_file << endl << "================================================================" << endl;
    generation_infos_file << "FAST SEQUENCE GENERATION" << endl;
    generation_infos_file << "Generated sequences in file: " << seq_filename << endl;
    generation_infos_file << "Generated sequences realizations in file: " << real_filename << endl;
    generation_infos_file << "Date: " << ctime(&tt) << endl;
    generation_infos_file << "Number of sequences = " << num_sequences << endl;
    generation_infos_file << "Number of threads = "
                          << (num_threads == 0 ? igor::fast::get_optimal_thread_count() : num_threads) << endl;

    // Initialize fast generator if needed
    igor::fast::FastGenerator &generator = get_fast_generator();

    // Configure generation
    igor::fast::FastGeneratorConfig config;
    config.num_threads = num_threads;
    config.show_progress = show_progress;
    if (seed >= 0) {
        config.base_seed = static_cast<uint64_t>(seed);
    } else {
        config.base_seed = igor::fast::draw_random_seed();
    }
    generation_infos_file << "Seed = " << config.base_seed << endl;

    // Progress callback
    auto progress_callback = [show_progress](size_t completed, size_t total) {
        if (show_progress) {
            show_progress_bar(cerr, static_cast<double>(completed) / total, "Fast sequence generation", 50);
        }
    };

    // Generate sequences
    generator.generate_to_files(num_sequences, seq_filename, real_filename, config, progress_callback);

    if (show_progress) {
        close_progress_bar(cerr, "Fast sequence generation", 50);
    }

    // Log statistics
    auto stats = generator.get_stats();
    chrono::system_clock::time_point end_time = chrono::system_clock::now();
    double elapsed = chrono::duration<double>(end_time - begin_time).count();

    generation_infos_file << "Total time: " << elapsed << " seconds" << endl;
    generation_infos_file << "Sequences per second: " << stats.sequences_per_second << endl;
    generation_infos_file << "Bytes written: " << stats.bytes_written << endl;

    cerr << "Fast generation complete: " << num_sequences << " sequences in " << elapsed << "s ("
         << static_cast<size_t>(stats.sequences_per_second) << " seq/s)" << endl;
}

igor::fast::FastGenerator &GenModel::get_fast_generator()
{
    if (!fast_generator_) {
        fast_generator_ = make_unique<igor::fast::FastGenerator>();
        fast_generator_->initialize(model_parms, model_marginals);
    }
    return *fast_generator_;
}
