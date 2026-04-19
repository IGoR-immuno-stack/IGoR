/*
 * Counter.cpp
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

#include <igor/Core/Counter.h>

using namespace std;

/*
 * By default the counter will be output in /tmp
 * I hope you're using a UNIX based system
 */
Counter::Counter(const string &path /* = "/tmp/"*/, bool last_iter /* = false*/)
    : path_to_file(path), last_iter_only(last_iter), fstreams_created(false)
{
    //Make sure the given path is a folder
    if (path_to_file.back() != '/') {
        path_to_file.push_back('/');
    }
}

Counter::~Counter()
{
    // TODO Auto-generated destructor stub
}

// ===== NEW INTERFACE (Phase 4.1) - NEW→OLD Adapters =====
// Pattern: Following Phase 3 strategy, adapters delegate from NEW interface to OLD implementation
// This allows iterate_wrap_up to call new interface immediately while counters migrate gradually

void Counter::initialize(const ModelContext& model)
{
    // NEW→OLD adapter: Delegate to legacy initialize_counter()
    // This allows iterate_wrap_up to call new interface (Phase 4.2)
    // while counters remain unmigrated
    
    // TODO Phase 4.3: Complete this adapter
    // CHALLENGE: ModelContext only has model_parameters (Marginal_array_p)
    // but initialize_counter() requires Model_Parms and Model_marginals
    // 
    // Options:
    // 1. Extend ModelContext to include Model_Parms and Model_marginals references
    // 2. Make counters override new interface (no adapter)
    // 
    // For now: Temporary stub - counters must override new interface OR
    // we complete adapter in Phase 4.3 when ModelContext is extended
    
    // TEMPORARY: Do nothing - forces counter migration
    // Each counter subclass should override new interface during Phase 4.4 migration
}

void Counter::count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model)
{
    // NEW→OLD adapter: Delegate to legacy count_scenario() (Option A)
    // This allows unmigrated counters to continue working with empty implementations
    // while new counters can override this method directly.
    //
    // Pattern: Following Phase 3 strategy, adapter delegates from NEW → OLD
    // This enables gradual per-counter migration with independent validation
    //
    // TODO Phase 4.5: REMOVE this adapter after all counters migrated
    // This is temporary infrastructure for gradual migration only
    
    // Delegate to legacy interface (mutable references from Scenario's stored refs)
    this->count_scenario(
        scenario.scenario_error_w_proba,
        scenario.scenario_proba,
        query.sequence,
        scenario.constructed_sequences_ref,  // Mutable reference
        scenario.seq_offsets_ref,           // Const reference (OK)
        model.events_map,
        scenario.mismatches_lists_ref       // Mutable reference
    );
}

// ===== OLD INTERFACE (DEPRECATED) - Default implementations =====

void Counter::count_scenario(
        long double scenario_seq_joint_proba, double scenario_probability, const string &original_sequence,
        Seq_type_str_p_map &constructed_sequences, const Seq_offsets_map &seq_offsets,
        const unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        Mismatch_vectors_map &mismatches_lists)
{
    //Do nothing
    //This is a virtual method in case the counter does not have anything to count at the scenario level
}

void Counter::count_sequence(double seq_likelihood, const Model_marginals &single_seq_marginals,
                             const Model_Parms &single_seq_model_parms)
{
    //Do nothing
    //This is a virtual method in case the counter does not have anything to count at the sequence level
}

void Counter::add_to_counter(shared_ptr<Counter> other)
{
    if (this->type() != other->type()) {
        throw invalid_argument("Cannot add counters of type " + this->type() + " and " + other->type());
    } else {
        this->add_checked(other);
        return;
    }
}

/*
 * Dump sequence specific information to file
 * This method should also clean all the counters for the next sequence
 */
void Counter::dump_sequence_data(int seq_index, int iteration_n)
{
    //Do nothing
}

void Counter::dump_data_summary(int iteration_n)
{
    //Do nothing
}

void Counter::set_path_to_files(const string &new_path)
{
    if (fstreams_created) {
        throw runtime_error(
                "Cannot change file path when file streams have already been created (in Counter::set_path_to_files)");
    } else {
        this->path_to_file = new_path;
    }
}
