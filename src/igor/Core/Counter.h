/*
 * Counter.h
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

#pragma once

#include <igor/Core/Utils.h>
#include <string>
#include <memory>
#include <fstream>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/IntStr.h>

// Phase 4.1: Context-based Counter interface
#include <igor/Core/Scenario.h>
#include <igor/Core/QuerySequenceContext.h>
#include <igor/Core/ModelContext.h>

#include <igorCoreExport.h>

/**
 * \class Counter Counter.h
 * \brief Scenario statistics recording abstract class
 * \author Q.Marcou
 * \version 1.0
 *
 * The Counter abstract class provides an interface to collect individual scenarios statistics and aggregate them in various ways.
 * 
 * PHASE 4 REFACTORING: Counter interface now uses context objects for cleaner, more maintainable code.
 * - New interface: initialize(ModelContext) and count_scenario(Scenario, QuerySequenceContext, ModelContext)
 * - Legacy interface: Preserved for compatibility, marked deprecated
 */
class CORE_EXPORT Counter
{
public:
    Counter(const std::string &path = "/tmp/", bool last_iter = false);
    virtual ~Counter();

    virtual std::string type() const = 0;

    // ===== NEW INTERFACE (Phase 4.1 - context-based) =====
    
    /**
     * @brief Initialize counter with model context (NEW)
     * 
     * Uses ModelContext which provides access to:
     * - events_map: Event structure for querying event properties
     * - model_parameters: Probability parameters
     * - offset_map: Event dependencies
     * 
     * This is the preferred interface going forward.
     * 
     * @param model Model configuration context
     * 
     * NOTE: Default implementation does nothing (adapter incomplete - ModelContext needs extension).
     * Counter subclasses should override this method during Phase 4.4 migration.
     * TODO Phase 4.3: Complete NEW→OLD adapter when ModelContext extended with Model_Parms/Model_marginals.
     */
    virtual void initialize(const ModelContext& model);
    
    /**
     * @brief Count a completed scenario (NEW)
     * 
     * Uses context-based interface:
     * - Scenario: Flattened view of completed scenario (zero-copy)
     * - QuerySequenceContext: Input sequence data
     * - ModelContext: Model structure for event queries
     * 
     * This is the preferred interface going forward.
     * 
     * @param scenario Flattened scenario view (probabilities + sequences)
     * @param query Input sequence context (read sequence data)
     * @param model Model context (event lookup, structure queries)
     * 
     * NOTE: Default implementation throws (no adapter - legacy interface incomplete).
     * Counter subclasses MUST override this method during Phase 4.4 migration.
     */
    virtual void count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model
    );
    
    // ===== OLD INTERFACE (parameter-based) - DEPRECATED =====
    
    /**
     * @brief Initialize counter (DEPRECATED)
     * @deprecated Use initialize(const ModelContext&)
     */
    [[deprecated("Use initialize(const ModelContext&)")]]
    virtual void initialize_counter(const Model_Parms &, const Model_marginals &) = 0;

    /**
     * @brief Count scenario (DEPRECATED)
     * @deprecated Use count_scenario(const Scenario&, const QuerySequenceContext&, const ModelContext&)
     */
    [[deprecated("Use count_scenario(const Scenario&, const QuerySequenceContext&, const ModelContext&)")]]
    virtual void
    count_scenario(long double, double, const std::string &, Seq_type_str_p_map &, const Seq_offsets_map &,
                   const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>, std::shared_ptr<Rec_Event>> &,
                   Mismatch_vectors_map &);
                   
    // TODO Phase 4: Refactor count_sequence to use contexts
    // virtual void count_sequence(const Scenario&, const ModelContext&);
    virtual void count_sequence(double, const Model_marginals &, const Model_Parms &);

    virtual void add_to_counter(std::shared_ptr<Counter>);
    virtual void add_checked(std::shared_ptr<Counter>) = 0;

    virtual void dump_sequence_data(int, int);
    virtual void dump_data_summary(int);

    bool is_last_iter_only() const { return last_iter_only; }
    std::string get_path_to_files() const { return path_to_file; }
    void set_path_to_files(const std::string &new_path);

    virtual std::shared_ptr<Counter> copy() const = 0;

protected:
    std::string path_to_file;
    bool last_iter_only;
    bool fstreams_created;
    //TODO create a unique identifier of the counter? Make something up to prevent to have twice the same counter??
};
