/*
 * Bestscenarioscounter.h
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

#include <igor/Core/Counter.h>
#include <queue>
#include <vector>
#include <list>

#include <igor/Core/Genechoice.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>

#include <igorCoreExport.h>

/**
 * \class Best_scenarios_counter Bestscenarioscounter.h
 * \brief Records the N best scenarios realizations and mismatches.
 * \author Q.Marcou
 * \version 1.0
 *
 * Implementation of the Counter abstract class.
 * Records the N most likely scenario realizations and mismatches and append it to a semicolon separated file.
 *
 * Recombination models are frequently degenerate: several genomic templates can be
 * indistinguishable over the region covered by a read, in which case the joint
 * probability of the corresponding scenarios is equal in exact arithmetic and differs
 * only by accumulated rounding. Comparing such probabilities bit for bit makes the
 * recorded scenarios depend on rounding noise, which is not reproducible across
 * platforms. Two settings control how that degeneracy is handled:
 *  - tie_digits: the number of significant digits the probabilities are compared on,
 *    scenarios agreeing to that precision being ordered by their realization indices
 *    instead. 0 compares the probabilities exactly.
 *  - keep_ties: whether a group of tied scenarios straddling the N-th rank is recorded
 *    in full, in which case more than N scenarios are written out for that sequence.
 *    Reporting an arbitrary member of the group instead would bias the output towards
 *    the lowest realization indices.
 */
class CORE_EXPORT Best_scenarios_counter : public Counter
{
public:
    Best_scenarios_counter();
    Best_scenarios_counter(size_t);
    Best_scenarios_counter(size_t, bool);
    Best_scenarios_counter(size_t, std::string);
    Best_scenarios_counter(size_t, std::string, bool);
    virtual ~Best_scenarios_counter();

    /**
     * \brief Number of significant digits scenario probabilities are compared on.
     * \param digits 0 to compare the probabilities exactly.
     */
    void set_tie_digits(size_t digits) { this->tie_digits = digits; };
    /**
     * \brief Whether tied scenarios straddling the N-th rank are all recorded.
     */
    void set_keep_ties(bool keep) { this->keep_ties = keep; };

    std::string type() const override { return "BestScenarioCounter"; }; //TODO return an enum

    // Context-based interface
    void initialize(const ModelContext& model) override;
    void count_scenario(
        const Scenario& scenario,
        const QuerySequenceContext& query,
        const ModelContext& model
    ) override;

    // LEGACY INTERFACE (DEPRECATED)
    void initialize_counter(const Model_Parms &, const Model_marginals &) override;

    void
    count_scenario(long double, double, const std::string &, Seq_type_str_p_map &, const Seq_offsets_map &,
                   const Events_map &,
                   Mismatch_vectors_map &) override;

    void count_sequence(double, const Model_marginals &, const Model_Parms &) override;

    void add_checked(std::shared_ptr<Counter>) override;

    void dump_sequence_data(int, int) override;

    std::shared_ptr<Counter> copy() const override;

    /**
     * A recorded scenario: its probability, the realization indices of every event in
     * model queue order, and the mismatch positions.
     */
    typedef std::tuple<double, std::vector<std::vector<int>>, std::vector<size_t>> Scenario_record;

    /// Default number of significant digits probabilities are compared on.
    static const size_t default_tie_digits = 12;

    size_t n_scenarios_counted;
    size_t tie_digits = default_tie_digits;
    bool keep_ties = true;
    std::shared_ptr<std::ofstream> output_scenario_file_ptr;

    std::vector<std::vector<int>> single_scenario_realizations;

    std::vector<size_t> single_scenario_mismatches_list;

    std::vector<Scenario_record> best_scenarios_vec;

    std::forward_list<std::shared_ptr<const Rec_Event>> event_fw_list;

private:
    /**
     * \brief Round a probability to tie_digits significant digits.
     *
     * Implemented by clearing the mantissa bits below the requested precision. Being a
     * pure function of the probability, the resulting comparison stays transitive,
     * which a tolerance based comparison would not be.
     */
    double tie_key(double) const;
    /// Strict weak ordering placing the least likely scenario first.
    bool record_is_worse(const Scenario_record &, const Scenario_record &) const;
    /// Insert a candidate at its place and drop the scenarios falling out of the list.
    void record_scenario(double);
    /// Collect the current realization indices and mismatches into the scratch members.
    void collect_realizations();
};
