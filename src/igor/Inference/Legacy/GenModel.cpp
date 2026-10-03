/*
 * GenModel.cpp
 *
 *  Created on: 3 nov. 2014
 *      Author: Quentin Marcou
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

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *      This class designs a generative model and supply all the methods to run a maximum likelihood estimate of the generative model
 */

#include <igor/Inference/Legacy/GenModel.h>
#include <igor/Alignment/Legacy/Aligner.h>
#include <igor/Model/Legacy/QuerySequenceContext.h>
#include <igor/Model/Legacy/ModelContext.h>
#include <igor/Model/Legacy/ScenarioContext.h>
#include <igor/Model/Legacy/ExplorationContext.h>
#include <igor/Model/Legacy/AccumulationContext.h>


namespace igor::inference::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

using namespace std;

namespace {
/*
 * Memory budget for the per chunk accumulators of the deterministic reduction
 * performed in infer_model(). A larger budget yields finer chunks, hence a better
 * load balance on small datasets, at the cost of holding more partial marginals in
 * memory. The reduction result depends on the resulting chunk decomposition only,
 * never on the number of threads, so changing this constant changes the last bits
 * of the inferred parameters.
 */
constexpr size_t REDUCTION_MEMORY_BUDGET_BYTES = 64ull * 1024ull * 1024ull;
} // namespace

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
        const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> &sequences,
        const int iterations, const std::string path, bool fast_iter, double likelihood_threshold /*=1e-25*/,
        bool viterbi_like /*false*/)
{
    return this->infer_model(sequences, iterations, path, fast_iter, likelihood_threshold, viterbi_like, 0.001);
}

bool GenModel::infer_model(
        const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> &sequences,
        const int iterations, const std::string path, bool fast_iter /*=true*/, double likelihood_threshold /*=1e-25*/,
        double proba_threshold_factor /*=0.001*/)
{
    return this->infer_model(sequences, iterations, path, fast_iter, likelihood_threshold, false,
                             proba_threshold_factor, INFINITY);
}

bool GenModel::infer_model(
        const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> &sequences,
        const int iterations, const string path, bool fast_iter /*=true*/,
        double likelihood_threshold /*=1e-25 by default*/, bool viterbi_like /*=false*/,
        double proba_threshold_factor /*=0.001 by default*/,
        double mean_number_seq_err_thresh /*= INFINITY by default*/)
{

    //If viterbi like only the best scenario is of interest
    if (viterbi_like) {
        cerr << "******************************************************************" << endl;
        cerr << "*\t\t RUNNING \"VITERBI\" LIKE ALGORITHM \t\t *" << endl
             << "* \t(only the best scenario will be taken into account)\t *" << endl;
        cerr << "******************************************************************" << endl;
        proba_threshold_factor = 1.0;
    }

    if (likelihood_threshold > 1.0) {
        throw invalid_argument("Likelihood threshold must be lesser or equal than one");
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
    log_file << "iteration_n;seq_processed;seq_index;nt_sequence;n_V_aligns;n_J_aligns;seq_likelihood;seq_mean_n_"
                "errors;seq_n_scenarios;seq_best_scenario;time"
             << endl;
    ofstream general_logs(path + string("inference_info.out"));
    //Dump all inference parameters to file
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
                 << "\t#(ratio between best scenario and current scenario needed to explore/count the scenario)"
                 << endl;
    general_logs << "Mean #errors threshold: " << mean_number_seq_err_thresh
                 << "\t#Needs a very good reason to be set to another value than INFINITY" << endl;

    //Get the total number of sequences to process
    const double total_number_seqs = sequences.size(); //Use a double for float division afterwards

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

    //Write initial condition to file
    this->model_marginals.write2txt(path + string("initial_marginals.txt"), this->model_parms);
    this->model_parms.write_model_parms(path + string("initial_model.txt"));

    /*
	 * First initialization creates file streams
	 * This is to make sure that the counter copies do not create new files each time
	 */
    // Construct ModelContext for counter initialization
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map =
            model_marginals.get_offsets_map(model_parms, model_queue);
    Events_map events_map_ref = model_parms.get_events_map();
    ModelContext model_context(
        model_marginals.marginal_array_smart_p,
        offset_map,
        events_map_ref,
        model_queue
    );

    for (map<size_t, shared_ptr<Counter>>::const_iterator iter = counters_list.begin(); iter != counters_list.end();
         ++iter) {
        (*iter).second->initialize(model_context);
    }

    /*
 * Reduction using OpenMP 4.0 standards
 *
	#pragma omp declare reduction(+:Model_marginals:omp_out+=omp_in) initializer(omp_priv = omp_orig.empty_copy())
	//Note: since Error_rate is an abstract class, the following is very dirty, need to find a better solution
	#pragma omp declare reduction(+:shared_ptr<Error_rate>:add_to_err_rate(omp_out,omp_in)) initializer(omp_priv = omp_orig->copy())
*/

    //Loop over iterations
    while (iteration_accomplished != iterations) {

        //double proba_threshold_factor;
        Model_marginals new_marginals = Model_marginals(model_parms);
        shared_ptr<Error_rate> error_rate_copy = model_parms.get_err_rate_p()->copy();

        //Initialize error rate copy
        error_rate_copy->initialize(model_parms.get_events_map());

        //Initialize counters for the log file
        size_t sequences_processed = 0;

        new_marginals.debug_marg_name = "new_marginals";

        const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> *sequence_util_ptr;

        //Take only best alignments if fast_iter
        vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> fast_iter_sequences;
        if (fast_iter && iteration_accomplished == 0) {
            fast_iter_sequences = sequences;
            for (unordered_map<Gene_class, vector<Alignment_data>>::const_iterator gc_align_iter =
                         std::get<2>(sequences.at(0)).begin();
                 gc_align_iter != std::get<2>(sequences.at(0)).end(); ++gc_align_iter) {
                if ((*gc_align_iter).first == D_gene)
                    continue;
                fast_iter_sequences = get_best_aligns(fast_iter_sequences, (*gc_align_iter).first);
            }
            sequence_util_ptr = &fast_iter_sequences;
        } else {
            sequence_util_ptr = &sequences;
        }

        /*
		 * Deterministic reduction setup.
		 *
		 * Floating point addition is not associative, so the summed marginals depend on
		 * the order in which the per sequence contributions are added. Accumulating per
		 * thread makes that order depend on the number of threads and on the timing of
		 * the dynamic schedule, which changes the last bits of the inferred parameters
		 * from one run to the next; scenarios that are exactly degenerate under the
		 * model are then ranked differently by the counters.
		 *
		 * Instead the sequences are partitioned into fixed contiguous chunks, each
		 * chunk accumulating into its own slot. schedule(dynamic,chunk_size) hands out
		 * exactly those chunks, so a chunk is always summed by a single thread in
		 * increasing sequence order, whichever thread happens to pick it up. Merging the
		 * slots in chunk order after the parallel region therefore gives the same result
		 * for any number of threads.
		 *
		 * The number of slots is bounded by a memory budget: with few sequences every
		 * sequence gets its own slot, which is exactly how a plain dynamic schedule
		 * behaves; with many sequences the chunks grow instead.
		 */
        const size_t n_seqs_to_process = sequence_util_ptr->size();
        const size_t marginals_bytes = new_marginals.get_length() * sizeof(long double);
        const size_t max_reduction_slots =
                max<size_t>(1, REDUCTION_MEMORY_BUDGET_BYTES / max<size_t>(marginals_bytes, 1));
        const size_t chunk_size =
                max<size_t>(1, (n_seqs_to_process + max_reduction_slots - 1) / max_reduction_slots);
        const size_t n_reduction_slots = (n_seqs_to_process + chunk_size - 1) / chunk_size;
        vector<Model_marginals> chunk_marginals(n_reduction_slots, new_marginals.empty_copy());
        vector<shared_ptr<Error_rate>> chunk_err_rates(n_reduction_slots);
        vector<map<size_t, shared_ptr<Counter>>> thread_counter_lists(omp_get_max_threads());

        general_logs << "Iteration " << iteration_accomplished + 1 << ": reducing over " << n_reduction_slots
                     << " chunk(s) of " << chunk_size << " sequence(s)" << endl;

        cerr << "Performing Evaluate/Inference iteration " << iteration_accomplished + 1 << endl;

/* omp parallel declaration using OpenMP 4.0 standards
		 * #pragma omp parallel for schedule(dynamic) reduction(+:error_rate_copy,new_marginals) firstprivate(model_queue,index_map,offset_map,model_marginals_copy,events_map , processed_events , safety_set , write_index_list) //num_threads(6)
		 */

        //Where the thread that folds the junction-length bounds publishes its events for the
        //others to adopt from. Shared, and re-created every iteration because the bound follows
        //the marginals. Sized inside the region, by the thread that fills it.
        vector<shared_ptr<Rec_Event>> len_proba_bound_source_events;

//Declare variables to use OpenMP 3.1 standards
#pragma omp parallel shared(new_marginals, error_rate_copy, sequences_processed, sequence_util_ptr, sequences, \
        len_proba_bound_source_events) \
        firstprivate(model_queue, proba_threshold_factor) //num_threads(1)
        {
            //Make single thread copies of objects for thread safety
            Model_Parms single_thread_model_parms(model_parms);
            unordered_map<Rec_Event_name, int> single_thread_index_map =
                    model_marginals.get_index_map(model_parms, model_queue);
            Model_marginals single_thread_model_marginals(model_marginals);
            shared_ptr<Error_rate> single_thread_err_rate = single_thread_model_parms.get_err_rate_p();
            Events_map events_map =
                    single_thread_model_parms.get_events_map();
            map<size_t, shared_ptr<Counter>> single_thread_counter_list;
            for (map<size_t, shared_ptr<Counter>>::const_iterator iter = this->counters_list.begin();
                 iter != this->counters_list.end(); ++iter) {
                //Copy only relevant counters for this iteration
                if ((not(*iter).second->is_last_iter_only()) or (iteration_accomplished == iterations - 1)) {
                    single_thread_counter_list.emplace((*iter).first, (*iter).second->copy());
                }
            }

            single_thread_model_marginals.debug_marg_name = "single thread model marginals";

            unordered_set<Rec_Event_name> init_processed_events;

            //Initialize Enum_fast_memory map and dual maps
            SafetyMatrix safety_set(single_thread_model_parms.get_seq_type_registry());
            Seq_type_str_p_map constructed_sequences(single_thread_model_parms.get_seq_type_registry());
            Mismatch_vectors_map mismatches_lists(single_thread_model_parms.get_seq_type_registry());
            // Conservative pruning track. For exact NT queries it is identical to
            // mismatches_lists; it exists so that IUPAC/motif queries can prune on the
            // mismatch count that no choice of query branch can reduce.
            Pruning_mismatch_floor_map pruning_mismatch_floor(
                    single_thread_model_parms.get_seq_type_registry());
            Seq_offsets_map seq_offsets(single_thread_model_parms.get_seq_type_registry());

            //Initialize downstream probas to 1. Sized from the thread-local model's frozen
            //registry, which outlives the parallel region the map is used in.
            Downstream_scenario_proba_bound_map downstream_proba_map(
                    single_thread_model_parms.get_seq_type_registry());
            downstream_proba_map.init_first_layer(1.0);

            list<shared_ptr<Rec_Event>> events_list = single_thread_model_parms.get_event_list();
            Index_map index_mapp(events_list.size());

            //Initialize index_map
            for (list<shared_ptr<Rec_Event>>::iterator event_iter = events_list.begin();
                 event_iter != events_list.end(); ++event_iter) {
                int event_index = (*event_iter)->get_event_identifier();
                index_mapp.request_layer(event_index);
                index_mapp.set(event_index, single_thread_index_map.at((*event_iter)->get_name()), 0);
                //TODO update proba bound

                //Get events probability upper bounds
                size_t event_size =
                        single_thread_model_marginals.get_event_size((*event_iter), single_thread_model_parms);

                (*event_iter)->set_event_marginal_size(event_size);
                (*event_iter)
                        ->set_crude_upper_bound_proba(single_thread_index_map.at((*event_iter)->get_name()), event_size,
                                                      single_thread_model_marginals.marginal_array_smart_p);
            }

            queue<shared_ptr<Rec_Event>> single_thread_model_queue =
                    single_thread_model_parms.get_model_queue(); //single_thread_parms.get_model_queue();

            queue<shared_ptr<Rec_Event>> init_single_thread_model_queue = single_thread_model_queue;
            unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> single_thread_offset_map =
                    model_marginals.get_offsets_map(model_parms, single_thread_model_queue);

            stack<shared_ptr<Rec_Event>> init_single_thread_stack;

            //Initialize events
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

            /*
			 * Initialize the array of next event pointers
			 * This array replaces the formerly copied queue<shared_ptr<Rec_Event>> (was copied at each iterate_wrap_up call)
			 * Each event will access the pointer corresponding to its identifier address when calling iterate inside iterate_wrap_up
			 * The last event will point to null pointer enabling to call the error_rate
			 */
            shared_ptr<Next_event_ptr> next_event_ptr_arr(
                    new Next_event_ptr[single_thread_model_parms.get_event_list().size()],
                    std::default_delete<Rec_Event *[]>());
            init_single_thread_model_queue = single_thread_model_queue;
            while (!init_single_thread_model_queue.empty()) {
                shared_ptr<Rec_Event> first_init_event = init_single_thread_model_queue.front();
                init_single_thread_model_queue.pop();
                if (!init_single_thread_model_queue.empty()) {
                    next_event_ptr_arr.get()[first_init_event->get_event_identifier()] =
                            init_single_thread_model_queue.front().get();
                } else {
                    //This is the last event thus we emplace a shared null pointer
                    next_event_ptr_arr.get()[first_init_event->get_event_identifier()] = NULL; //Next_event_ptr(NULL);
                }
            }

            //Initialize error rate
            single_thread_err_rate->initialize(events_map);
            single_thread_err_rate->set_viterbi_run(viterbi_like);

            //Initialize Counters
            // Construct ModelContext for counter initialization
            ModelContext single_thread_model_context(
                single_thread_model_marginals.marginal_array_smart_p,
                single_thread_offset_map,
                events_map,
                single_thread_model_queue
            );

            for (map<size_t, shared_ptr<Counter>>::iterator iter = single_thread_counter_list.begin();
                 iter != single_thread_counter_list.end(); ++iter) {
                (*iter).second->initialize(single_thread_model_context);
            }
#pragma omp single nowait
            {
                cerr << "Initializing probability bounds..." << endl;
            }
            //The events in reverse queue order, which is the order the junction-length fold
            //below needs: each event's suffix is folded before anything reads it. (A per-thread
            //"crude" bound used to be initialized in this loop too; nothing read it at scenario
            //time, and R14 deleted it.)
            vector<shared_ptr<Rec_Event>> len_proba_bound_order;
            len_proba_bound_order.reserve(init_single_thread_stack.size());
            while (!init_single_thread_stack.empty()) {
                len_proba_bound_order.push_back(init_single_thread_stack.top());
                init_single_thread_stack.pop();
            }

            //The junction-length bound, by contrast, is a function of the marginals alone, so all
            //N threads used to compute the identical answer N times (section 2.5 finding 7 of
            //docs/ITERATE_GENERIC_REWRITE_PLAN.md). One thread folds it and the rest copy the
            //result, which is orders of magnitude below folding it again.
            bool folded_len_proba_bounds = false;
#pragma omp single
            {
                len_proba_bound_source_events.assign(single_thread_model_parms.get_event_list().size(), nullptr);
                for (const shared_ptr<Rec_Event> &len_proba_init_event : len_proba_bound_order) {
                    queue<shared_ptr<Rec_Event>> tmp_init_proba_single_thread_model_queue = single_thread_model_queue;
                    while (tmp_init_proba_single_thread_model_queue.front() != len_proba_init_event) {
                        tmp_init_proba_single_thread_model_queue.pop();
                    }
                    tmp_init_proba_single_thread_model_queue.pop();
                    len_proba_init_event->initialize_Len_proba_bound(
                            tmp_init_proba_single_thread_model_queue,
                            single_thread_model_marginals.marginal_array_smart_p, index_mapp);
                    len_proba_bound_source_events[len_proba_init_event->get_event_identifier()] = len_proba_init_event;
                }
                folded_len_proba_bounds = true;
            }
            //omp single's implicit barrier is what publishes the folds above to the readers below.
            if (not folded_len_proba_bounds) {
                for (const shared_ptr<Rec_Event> &len_proba_init_event : len_proba_bound_order) {
                    len_proba_init_event->adopt_Len_proba_bound(
                            *len_proba_bound_source_events[len_proba_init_event->get_event_identifier()]);
                }
            }
#pragma omp single nowait
            {
                cerr << "Initialization of probability bounds over." << endl;
            }

            //Now let all the events in the need of it get their own updated copy of the marginals
            init_single_thread_model_queue = single_thread_model_queue;
            while (!init_single_thread_model_queue.empty()) {
                init_single_thread_model_queue.front()->update_event_internal_probas(
                        single_thread_model_marginals.marginal_array_smart_p, index_map);
                init_single_thread_model_queue.pop();
            }

            chrono::system_clock::time_point single_seq_begin;
            chrono::duration<double> seq_time;

            //Loop over sequences in parallel, using the number of threads declared previously when declaring the parallel section
            //Chunks are handed out on demand to keep the load balanced, but each chunk
            //accumulates into its own slot so that the reduction does not depend on
            //which thread picked it up nor on how many threads there are
#pragma omp for schedule(dynamic, chunk_size) nowait
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n_seqs_to_process); ++i) {

                const size_t reduction_slot = static_cast<size_t>(i) / chunk_size;
                const auto &seq = (*sequence_util_ptr)[i]; // index access

                single_seq_begin = chrono::system_clock::now();

                //Make a copy of the queue that can be modified in iterate
                queue<shared_ptr<Rec_Event>> model_queue_copy(single_thread_model_queue);

                //Get the first event from the queue
                shared_ptr<Rec_Event> first_event = model_queue_copy.front();
                model_queue_copy.pop();

                //Initialize single seq marginals
                Model_marginals single_seq_marginals = single_thread_model_marginals.empty_copy();
                double init_proba = 1;
                //double init_tmp_err_w_proba = 1;
                double max_proba_scenario = likelihood_threshold / proba_threshold_factor;

                Int_Str int_sequence = nt2int(get<1>(seq));

                //cout<<int_sequence<<endl;

                single_seq_marginals.debug_marg_name = "single_seq_marginals";

                /*
				 * Call iterate on the first event
				 * The method will be called recursively for each event, this is equivalent to a nested loop and enumerates all possible scenarios
				 * The weight of each recombination scenario is added to the single_seq_marginals on the fly
				 *
				 * Context-based interface
				 */
                try {

                    // Construct context objects from existing variables
                    QuerySequenceContext query(
                        get<1>(seq),      // sequence
                        int_sequence,     // int_sequence
                        get<2>(seq)       // gene_alignments
                    );

                    ModelContext model(
                        single_thread_model_marginals.marginal_array_smart_p,  // model_parameters
                        single_thread_offset_map,                              // offset_map
                        events_map,                                        // events_map
                        single_thread_model_queue                              // model_queue
                    );

                    ScenarioContext scenario(
                        init_proba,              // scenario_proba
                        constructed_sequences,   // constructed_sequences
                        seq_offsets,            // seq_offsets
                        mismatches_lists        // mismatches_lists
                    );

                    ExplorationContext exploration(
                        downstream_proba_map,        // downstream_proba_map
                        max_proba_scenario,          // seq_max_prob_scenario
                        proba_threshold_factor,      // proba_threshold_factor
                        index_mapp,                  // index_map
                        next_event_ptr_arr,          // next_event_ptr_arr
                        safety_set,                  // safety_set
                        pruning_mismatch_floor       // pruning_mismatch_floor
                    );

                    AccumulationContext accumulation(
                        single_seq_marginals.marginal_array_smart_p,  // updated_marginals
                        single_thread_counter_list,                   // counters
                        single_thread_err_rate                        // error_rate
                    );

                    // Call context-based iterate
                    first_event->iterate(query, model, scenario, exploration, accumulation);

                }

                catch (exception &except) {
                    general_logs << "Exception caught calling iterate() on sequence:" << endl;
                    general_logs << get<1>(seq) << " with index " << get<0>(seq) << endl;
                    general_logs << "Exception caught after " << single_thread_err_rate->debug_number_scenarios
                                 << " scenarios explored" << endl;
                    general_logs << endl;
                    general_logs << "Throwing exception now..." << endl << endl;
                    general_logs << except.what() << endl;
                    throw;
                }

                //Normalize the weights on the single_seq_marginal so that each sequence has the same weight when merged to the single_thread_marginals
                single_thread_err_rate->norm_weights_by_seq_likelihood(single_seq_marginals.marginal_array_smart_p,
                                                                       single_seq_marginals.get_length());
                seq_time = chrono::system_clock::now() - single_seq_begin;
#pragma omp critical(dump_seq_info)
                {
                    ++sequences_processed;
                    //Output useful infos in the log file
                    //log_file<<iteration_accomplished<<";"<<sequences_processed<<";"<<(seq).first<<";"<<(seq).second.at(V_gene).size()<<";"<<(seq).second.at(D_gene).size()<<";"<<(seq).second.at(J_gene).size()<<";"<<single_thread_err_rate->get_seq_probability()<<";"<<single_thread_err_rate->get_seq_likelihood()<<";"<<single_thread_err_rate->debug_number_scenarios<<";"<<max_proba_scenario<<endl;
                    log_file << iteration_accomplished << ";" << sequences_processed << ";" << get<0>(seq) << ";"
                             << get<1>(seq) << ";" << get<2>(seq).at(V_gene).size() << ";"
                             << get<2>(seq).at(J_gene).size() << ";" << single_thread_err_rate->get_seq_likelihood()
                             << ";" << single_thread_err_rate->get_seq_mean_error_number() << ";"
                             << single_thread_err_rate->debug_number_scenarios << ";" << max_proba_scenario << ";"
                             << seq_time.count() << endl;
                }
                for (map<size_t, shared_ptr<Counter>>::iterator iter = single_thread_counter_list.begin();
                     iter != single_thread_counter_list.end(); ++iter) {
                    iter->second->count_sequence(single_thread_err_rate->get_seq_likelihood(), single_seq_marginals,
                                                 single_thread_model_parms);
#pragma omp critical(dump_counters)
                    {
                        (*iter).second->dump_sequence_data(get<0>(seq), iteration_accomplished);
                    }
                }

                if (single_thread_err_rate->get_seq_mean_error_number() <= mean_number_seq_err_thresh) {
                    //Add weighed errors to the normalized error counter
                    single_thread_err_rate->add_to_norm_counter();

                    //Add the single_seq_marginals to this chunk's marginals
                    chunk_marginals[reduction_slot] += single_seq_marginals;
                } else {
                    //Erase seq specific counters so that it won't contribute to the error rate
                    single_thread_err_rate->clean_seq_counters();
                }

                //Drain the error rate accumulators at the end of each chunk so that they
                //are summed in chunk order as well
                if (static_cast<size_t>(i) + 1 == min((reduction_slot + 1) * chunk_size, n_seqs_to_process)) {
                    chunk_err_rates[reduction_slot] = single_thread_err_rate->copy();
                    add_to_err_rate(chunk_err_rates[reduction_slot].get(), single_thread_err_rate.get());
                    single_thread_err_rate->clear_accumulators();
                }

#pragma omp critical(update_progress_bar)
                {
                    if (sequences_processed % 100 == 0) {
                        //Output current progress to cerr
                        show_progress_bar(cerr, sequences_processed / total_number_seqs,
                                          "Iteration " + to_string(iteration_accomplished + 1), 50);
                    }
                }
            }

            //Publish this thread's counters, they are merged in thread order below.
            //The team size can never exceed the omp_get_max_threads() value the vector
            //was sized with, so the index is always in range.
            thread_counter_lists[omp_get_thread_num()] = single_thread_counter_list;
        }

        //Merge the chunk accumulators in chunk order: this is what makes the inferred
        //parameters independent of the thread count and of the schedule timing
        for (size_t slot = 0; slot != n_reduction_slots; ++slot) {
            new_marginals += chunk_marginals[slot];
            if (chunk_err_rates[slot]) {
                add_to_err_rate(error_rate_copy.get(), chunk_err_rates[slot].get());
            }
        }

        /*
		 * FIXME counters are still accumulated per thread, so the summary they write
		 * (the coverage and error counters, the only ones summing over sequences) still
		 * depends on the sequence to thread assignment in the last bits. Merging in
		 * thread order only removes the dependency on the order in which the threads
		 * finish. Giving them chunk slots the way the marginals have would cost one
		 * counter copy per chunk, which the coverage counter is too large for.
		 *
		 * None of the counters feeds back into the model parameters, so this does not
		 * affect the inference itself. The proper fix is to make the accumulator
		 * containers themselves aware of the chunked reduction rather than to special
		 * case every accumulator here.
		 */
        for (const auto &thread_counters : thread_counter_lists) {
            for (const auto &counter_pair : thread_counters) {
                counters_list.at(counter_pair.first)->add_to_counter(counter_pair.second);
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

        //Close current iteration progress bar
        close_progress_bar(cerr, "Iteration " + to_string(iteration_accomplished), 50);
    }
    //Create a copy of the last iteration results with identifiable name
    this->model_marginals.write2txt(path + string("final_marginals.txt"), this->model_parms);
    this->model_parms.write_model_parms(path + string("final_parms.txt"));

    return 0;
}
/**
 * \deprecated This function used to store generated sequences in memory, and quickly overloaded it for large number of generated sequences.
 */

/*
 * Extract the best alignment for each sequence for a given gene class (used for the fast iter)
 */
vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>>
get_best_aligns(const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> &all_aligns,
                Gene_class gc)
{

    vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> all_aligns_copy(all_aligns);
    for (vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>>::iterator seq_iter =
                 all_aligns_copy.begin();
         seq_iter != all_aligns_copy.end(); ++seq_iter) {
        vector<Alignment_data> &align_vect = get<2>((*seq_iter)).at(gc); //TODO add exception

        //Get align best score
        double best_score = -1;
        for (vector<Alignment_data>::const_iterator align_iter = align_vect.begin(); align_iter != align_vect.end();
             ++align_iter) {
            if ((*align_iter).score > best_score) {
                best_score = (*align_iter).score;
            }
        }

        vector<Alignment_data> new_align_vect;
        for (vector<Alignment_data>::const_iterator align_iter = align_vect.begin(); align_iter != align_vect.end();
             ++align_iter) {
            if ((*align_iter).score == best_score) {
                new_align_vect.push_back((*align_iter));
            }
        }

        align_vect = new_align_vect;
    }

    return all_aligns_copy;
}

} // namespace igor::inference::legacy
