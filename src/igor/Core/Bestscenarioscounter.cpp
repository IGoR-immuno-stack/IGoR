/*
 * Bestscenarioscounter.cpp
 *
 *  Created on: Aug 19, 2016
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
 */

#include <igor/Core/Bestscenarioscounter.h>

#include <cmath>
#include <limits>

using namespace std;

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios) : Counter(), n_scenarios_counted(n_scenarios) { }

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios, bool is_last_iter_only)
    : Best_scenarios_counter(n_scenarios)
{
    this->last_iter_only = is_last_iter_only;
}

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios, string path)
    : Counter(path), n_scenarios_counted(n_scenarios)
{
}

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios, string path, bool is_last_iter_only)
    : Best_scenarios_counter(n_scenarios, path)
{
    this->last_iter_only = is_last_iter_only;
}

Best_scenarios_counter::Best_scenarios_counter() : Best_scenarios_counter(1, true)
{
    // TODO Auto-generated constructor stub
}

Best_scenarios_counter::~Best_scenarios_counter()
{
    // TODO Auto-generated destructor stub
}

double Best_scenarios_counter::tie_key(double proba) const
{
    if (this->tie_digits == 0) {
        return proba;
    }
    //Number of mantissa bits carrying the requested number of decimal digits
    const int kept_bits = static_cast<int>(ceil(this->tie_digits * 3.321928094887362));
    if (kept_bits >= numeric_limits<double>::digits) {
        //More precision was asked for than a double carries: nothing left to round
        return proba;
    }
    int exponent = 0;
    const double mantissa = frexp(proba, &exponent);
    return ldexp(nearbyint(ldexp(mantissa, kept_bits)), exponent - kept_bits);
}

bool Best_scenarios_counter::record_is_worse(const Scenario_record &lhs, const Scenario_record &rhs) const
{
    const double lhs_key = this->tie_key(get<0>(lhs));
    const double rhs_key = this->tie_key(get<0>(rhs));
    if (lhs_key != rhs_key) {
        return lhs_key < rhs_key;
    }
    //Tied to the requested precision: order on the realizations instead so that the
    //recorded scenarios do not depend on the order they were enumerated in, the lowest
    //realization indices ranking best
    if (get<1>(lhs) != get<1>(rhs)) {
        return get<1>(lhs) > get<1>(rhs);
    }
    return get<2>(lhs) > get<2>(rhs);
}

void Best_scenarios_counter::collect_realizations()
{
    for (forward_list<shared_ptr<const Rec_Event>>::const_iterator iter = this->event_fw_list.begin();
         iter != this->event_fw_list.end(); ++iter) {
        this->single_scenario_realizations.emplace_back((*iter)->get_current_realizations_index_vec());
    }
}

void Best_scenarios_counter::record_scenario(double scenario_seq_joint_proba)
{
    Scenario_record candidate(scenario_seq_joint_proba, this->single_scenario_realizations,
                              this->single_scenario_mismatches_list);

    auto jter = this->best_scenarios_vec.begin();
    while ((jter != this->best_scenarios_vec.end()) and this->record_is_worse(*jter, candidate)) {
        ++jter;
    }
    this->best_scenarios_vec.emplace(jter, move(candidate));

    if (this->best_scenarios_vec.size() <= this->n_scenarios_counted) {
        return;
    }

    //The list is ordered from the least to the most likely scenario, so the ones
    //falling out of it are the leading ones
    auto first_kept = this->best_scenarios_vec.end() - this->n_scenarios_counted;
    if (this->keep_ties) {
        //Walk back over the scenarios tied with the N-th best one: recording only part
        //of a degenerate group would make the reported member arbitrary
        const double cutoff_key = this->tie_key(get<0>(*first_kept));
        while ((first_kept != this->best_scenarios_vec.begin())
               and (this->tie_key(get<0>(*(first_kept - 1))) == cutoff_key)) {
            --first_kept;
        }
    }
    this->best_scenarios_vec.erase(this->best_scenarios_vec.begin(), first_kept);
}

// ===== CONTEXT-BASED INTERFACE =====

void Best_scenarios_counter::initialize(const ModelContext& model) {
    if (not fstreams_created) {
        output_scenario_file_ptr = shared_ptr<ofstream>(new ofstream);
        this->output_scenario_file_ptr->open(path_to_file + "best_scenarios_counts.csv");
        
        // Create the header with event names
        (*this->output_scenario_file_ptr.get()) << "seq_index;scenario_rank;scenario_proba_cond_seq";
        
        // Use model_queue which provides correct topological ordering
        // This matches the legacy behavior that used parms.get_model_queue()
        auto event_queue_copy = model.model_queue;  // Copy queue to iterate
        while (!event_queue_copy.empty()) {
            shared_ptr<Rec_Event> event_ptr = event_queue_copy.front();
            (*this->output_scenario_file_ptr.get()) << ";" << event_ptr->get_name();
            this->event_fw_list.emplace_front(event_ptr);
            event_queue_copy.pop();
        }
        event_fw_list.reverse();
        
        (*this->output_scenario_file_ptr.get()) << ";Mismatches" << endl;
        
        fstreams_created = true;
    } else {
        // File already created, just populate event_fw_list
        auto event_queue_copy = model.model_queue;  // Copy queue to iterate
        while (!event_queue_copy.empty()) {
            this->event_fw_list.emplace_front(event_queue_copy.front());
            event_queue_copy.pop();
        }
        event_fw_list.reverse();
    }
}

void Best_scenarios_counter::count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model)
{
    // Use scenario_error_w_proba (already normalized by error rate)
    const double scenario_seq_joint_proba = scenario.scenario_error_w_proba;

    //Discard the scenarios that cannot make the list before collecting their
    //realizations, which is the expensive part. Scenarios tied with the least likely
    //recorded one are kept: whether they make the list depends on their realizations.
    if ((this->best_scenarios_vec.size() >= this->n_scenarios_counted)
        and (this->tie_key(scenario_seq_joint_proba) < this->tie_key(get<0>(this->best_scenarios_vec[0])))) {
        return;
    }

    // Collect realization indices from events
    this->collect_realizations();

    // Get mismatches and add them to the mismatch list
    for (Seq_type seq_type : { V_gene_seq, D_gene_seq, J_gene_seq }) {
        if (scenario.mismatches[seq_type]) {
            const vector<size_t> &mismatch_list = *scenario.mismatches[seq_type];
            single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end(), mismatch_list.begin(),
                                                   mismatch_list.end());
        }
    }

    this->record_scenario(scenario_seq_joint_proba);

    // Clean up temporary containers
    this->single_scenario_realizations.clear();
    single_scenario_mismatches_list.clear();
}

// ===== LEGACY INTERFACE (DEPRECATED) =====

void Best_scenarios_counter::count_scenario(
        long double scenario_seq_joint_proba, double scenario_probability, const string &original_sequence,
        Seq_type_str_p_map &constructed_sequences, const Seq_offsets_map &seq_offsets,
        const Events_map &events_map,
        Mismatch_vectors_map &mismatches_lists)
{
    const double joint_proba = scenario_seq_joint_proba;

    //See the context based overload
    if ((this->best_scenarios_vec.size() >= this->n_scenarios_counted)
        and (this->tie_key(joint_proba) < this->tie_key(get<0>(this->best_scenarios_vec[0])))) {
        return;
    }

    this->collect_realizations();

    // Get mismatches and add them to the mismatch list
    for (Seq_type seq_type : { V_gene_seq, D_gene_seq, J_gene_seq }) {
        if (mismatches_lists.exist(seq_type)) {
            const vector<size_t> &mismatch_list = *mismatches_lists.at(seq_type);
            single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end(), mismatch_list.begin(),
                                                   mismatch_list.end());
        }
    }

    this->record_scenario(joint_proba);

    this->single_scenario_realizations.clear();
    single_scenario_mismatches_list.clear();
}

void Best_scenarios_counter::count_sequence(double seq_likelihood, const Model_marginals &single_seq_marginals,
                                            const Model_Parms &single_seq_model_parms)
{
    for (auto iter = this->best_scenarios_vec.begin();
         iter != this->best_scenarios_vec.end(); ++iter) {
        get<0>(*iter) /= seq_likelihood;
        //If an exception is thrown here there is a problem upstream
    }
}

void Best_scenarios_counter::initialize_counter(const Model_Parms &parms, const Model_marginals &marginals)
{
    if (not fstreams_created) {
        output_scenario_file_ptr = shared_ptr<ofstream>(new ofstream);
        this->output_scenario_file_ptr->open(path_to_file + "best_scenarios_counts.csv");
        //Create the header
        (*this->output_scenario_file_ptr.get()) << "seq_index;scenario_rank;scenario_proba_cond_seq";
        auto event_queue = parms.get_model_queue();

        while (not event_queue.empty()) {
            shared_ptr<Rec_Event> event_ptr = event_queue.front();
            (*this->output_scenario_file_ptr.get()) << ";" << event_ptr->get_name();
            this->event_fw_list.emplace_front(event_ptr);
            event_queue.pop();
        }
        event_fw_list.reverse();

        (*this->output_scenario_file_ptr.get()) << ";Mismatches" << endl;

        fstreams_created = true;
    } else {
        //Still need to fill in the fw list
        auto event_queue = parms.get_model_queue();
        while (not event_queue.empty()) {
            shared_ptr<Rec_Event> event_ptr = event_queue.front();
            this->event_fw_list.emplace_front(event_ptr);
            event_queue.pop();
        }
        event_fw_list.reverse();
    }
}

void Best_scenarios_counter::add_checked(shared_ptr<Counter> counter)
{
    return;
}

void Best_scenarios_counter::dump_sequence_data(int seq_index, int iteration_n)
{
    stringstream ss;
    size_t counter = 1;
    for (auto iter =
                 this->best_scenarios_vec.rbegin();
         iter != this->best_scenarios_vec.rend(); ++iter) {
        ss << seq_index << ";" << counter << ";" << get<0>(*iter);
        //Loop over events
        for (const vector<int> &real_vec : get<1>(*iter)) {
            ss << ";(";
            //Loop over event realizations
            for (vector<int>::const_iterator jter = real_vec.begin(); jter != real_vec.end(); ++jter) {
                if (jter != real_vec.begin()) {
                    ss << ",";
                }
                ss << (*jter);
            }
            ss << ")";
        }
        ss << ";(";
        //Loop over mismatches
        const vector<size_t> &mismatches_list = get<2>(*iter);
        for (auto kter = mismatches_list.begin(); kter != mismatches_list.end(); ++kter) {
            if (kter != mismatches_list.begin()) {
                ss << ",";
            }
            ss << (*kter);
        }
        ss << ")" << "\n";
        ++counter;
    }
    string str = ss.str();
#pragma omp critical(dump_best_scenarios)
    {
        (*this->output_scenario_file_ptr.get()) << str;
    }
    best_scenarios_vec.clear();
}

shared_ptr<Counter> Best_scenarios_counter::copy() const
{
    shared_ptr<Best_scenarios_counter> counter_copy_ptr(new Best_scenarios_counter(this->n_scenarios_counted));
    counter_copy_ptr->tie_digits = this->tie_digits;
    counter_copy_ptr->keep_ties = this->keep_ties;
    counter_copy_ptr->fstreams_created = this->fstreams_created;
    if (this->fstreams_created) {
        counter_copy_ptr->output_scenario_file_ptr = this->output_scenario_file_ptr;
    } else {
        throw runtime_error("Counters should not be copied before stream initalization");
    }
    return counter_copy_ptr;
}
