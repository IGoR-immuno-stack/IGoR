/*
 * Errorrate.h
 *
 *  Created on: Jan 22, 2015
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
#include <igor/Core/IntStr.h>
#include <unordered_map>
#include <utility>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <random>
#include <queue>
#include <memory>

// Context objects for refactored iterate()
#include <igor/Core/QuerySequenceContext.h>
#include <igor/Core/ModelContext.h>
#include <igor/Core/ScenarioContext.h>
#include <igor/Core/ExplorationContext.h>
#include <igor/Core/AccumulationContext.h>

//Debug
#include <iostream>
#include <cmath>

#include <igorCoreExport.h>

//Forward declare Rec_event
class Rec_Event;

/**
 * \class Error_rate Error_rate.h
 * \brief Abstract class for generic error models behavior.
 * \author Q.Marcou
 * \version 1.0
 *
 * Base class for defining different error models such as additive or non-additive hypermutation models.
 * Errors are assessed when all RecEvent iterate have been processed (terminal leaf of the scenario tree)
 *
 * ARCHITECTURAL NOTE: Error_rate Position in Context Architecture
 *
 * Error_rate conceptually behaves like a terminal event in the recombination
 * graph (final probability computation after sequence construction), but has
 * key differences from Rec_Event:
 *
 * 1. MUTABLE ACCUMULATION STATE:
 *    Error_rate maintains per-sequence accumulators (seq_likelihood, seq_weighted_er)
 *    that are updated during iterate_wrap_up() calls. Rec_Events are stateless
 *    during iteration (parameters read from const ModelContext).
 *
 * 2. ACCUMULATOR + MODEL DUALITY:
 *    Error_rate is both:
 *    - A probabilistic model (has parameters: model_rate)
 *    - An accumulator (collects statistics across scenarios/sequences)
 *
 *    Rec_Events are only probabilistic models (parameters external in marginal_array).
 *
 * 3. CONTEXT PLACEMENT:
 *    - Rec_Events: Parameters in ModelContext (const - read-only model structure)
 *    - Error_rate: Lives in AccumulationContext (mutable - accumulation logic)
 *
 *    This placement respects const-correctness: ModelContext must stay const.
 *
 * 4. INTERFACE DESIGN:
 *    compute_scenario_error_probability() deliberately matches Rec_Event::iterate()
 *    signature pattern (contexts in, logic executed) to enable future conceptual
 *    unification while respecting const-correctness requirements.
 *
 * FUTURE UNIFICATION PATH:
 * If Rec_Event interface evolves to support terminal accumulators (mutable internal
 * state during iteration), Error_rate could inherit from Rec_Event. This would require:
 * - Rec_Event::iterate() accepting mutable self-reference or accumulation context
 * - Error_rate as last event in model_queue
 * - Rethinking Rec_Event contract to allow stateful iteration
 *
 * The current context-based interface makes this transition straightforward if/when needed.
 */
class CORE_EXPORT Error_rate
{
public:
    Error_rate();
    virtual ~Error_rate();

    // ===== CONTEXT-BASED INTERFACE =====

    /**
     * @brief Compute error-weighted scenario probability (context-based interface)
     *
     * This signature deliberately matches Rec_Event::iterate(contexts...) pattern
     * to enable future conceptual unification while respecting const-correctness.
     *
     * Error_rate uses ScenarioContext (full memory layer access) not flattened
     * Scenario view because:
     * - It's an accumulator (like a non-branching terminal Rec_Event), not a passive observer
     * - Future designs may use multiple Error_rate objects on different regions
     * - Needs flexibility to inspect memory layer maps
     *
     * Note: Error_rate stays in AccumulationContext (mutable) not ModelContext (const)
     * because it accumulates per-sequence state during iteration.
     *
     * @param query Input sequence context (const - read-only)
     * @param model Model structure context (const - read-only)
     * @param scenario Scenario state (mutable - scenario.scenario_error_w_proba updated)
     * @param exploration Exploration policy (mutable - threshold checking)
     * @return Error-weighted probability (scenario_proba * P(errors | model)), or 0 if below threshold
     *
     * MODIFIES:
     * - Internal per-sequence accumulators (seq_likelihood, seq_weighted_er, etc.)
     * - scenario.scenario_error_w_proba (stores computed value)
     */
    virtual double compute_scenario_error_probability(
        const QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration
    );

    // ===== LEGACY INTERFACE (DEPRECATED) =====

    /**
     * @brief Legacy error probability computation (DEPRECATED)
     * @deprecated Use compute_scenario_error_probability(contexts...) instead
     */
    [[deprecated("Use compute_scenario_error_probability(contexts...)")]]
    virtual double compare_sequences_error_prob(
            double, const std::string &, Seq_type_str_p_map &, const Seq_offsets_map &,
            const Events_map &,
            Mismatch_vectors_map &, double &, const double &) = 0;
    virtual void update() = 0;
    virtual void
    initialize(const Events_map &);
    bool is_updated() const { return updated; }
    void update_value(bool update_status) { updated = update_status; };
    virtual void add_to_norm_counter() = 0;
    virtual void clean_seq_counters() = 0;
    void norm_weights_by_seq_likelihood(Marginal_array_p &, const size_t, const double seq_weight = 1);
    virtual void write2txt(std::ofstream &) = 0;
    virtual std::shared_ptr<Error_rate> copy() const = 0;
    virtual std::string type() const = 0;
    virtual Error_rate *add_checked(Error_rate *) = 0;
    double get_model_likelihood() const { return model_log_likelihood; }
    double get_seq_likelihood() const { return seq_likelihood; }
    double get_seq_probability() const { return seq_probability; }
    double get_seq_mean_error_number() const;
    virtual const double &get_err_rate_upper_bound(size_t, size_t) = 0;
    virtual void build_upper_bound_matrix(size_t, size_t) = 0;
    virtual int get_number_non_zero_likelihood_seqs() const = 0;
    virtual std::queue<int> generate_errors(std::string &, std::mt19937_64 &) const = 0;
    void set_viterbi_run(bool viterbi_like) { viterbi_run = viterbi_like; }
    int debug_number_scenarios;

protected:
    bool updated;
    long double model_log_likelihood;
    int number_seq;
    long double seq_likelihood;
    double seq_mean_error_number;
    long double scenario_new_proba; //TODO rename this guy
    long double seq_probability; //Probability of generating one sequence without taking errors into account
    bool viterbi_run;
    Matrix<double> upper_bound_proba_mat; //Store the value of the error cost of i errors and j no errors
    size_t max_err;
    size_t max_noerr;
};

void add_to_err_rate(Error_rate *, Error_rate *);
