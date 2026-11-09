/*
 * GenModel.h
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
 *
 *      This class designs a generative model and supply all the methods to run a maximum likelihood estimate of the generative model
 */

#pragma once

#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Counter.h>
#include <igor/Model/Legacy/Model_marginals.h>
#include <igor/Model/Legacy/Errorrate.h>
#include <igor/Core/Legacy/Utils.h>
#include <list>
#include <map>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <stdexcept>
#include <stack>
#include <memory>

#include <igor/Inference/Export.h>

//Make typedef for the function pointers

namespace igor::inference::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

/**
 * \class GenModel GenModel.h
 * \brief High level V(D)J generative model.
 * \author Q.Marcou
 * \version 1.0
 *
 * Highest level class to model the V(D)J recombination and subsequent processes.
 * It contains the model's graph structure (Model_Parms), the associated probability distribution (Model_Marginals).
 * The GenModel class provides high level functions to perform inference / sequence annotation as well as generating random sequences from the model.
 */
class INFERENCE_EXPORT GenModel
{
public:
    GenModel(const Model_Parms &);
    GenModel(const Model_Parms &, const Model_marginals &);
    GenModel(const Model_Parms &, const Model_marginals &, const std::map<size_t, std::shared_ptr<Counter>> &);
    //TODO: add all the necessary constructors: with just model_parms, with model_parms and marginals
    virtual ~GenModel();

    bool infer_model(
            const std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
                    &sequences,
            const int iterations, const std::string path, bool fast_iter, double likelihood_threshold = 1e-25,
            bool viterbi_like = false);
    bool infer_model(
            const std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
                    &sequences,
            const int iterations, const std::string path, bool fast_iter = true, double likelihood_threshold = 1e-25,
            double proba_threshold_factor = 0.001);
    bool infer_model(
            const std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
                    &sequences,
            const int iterations, const std::string path, bool fast_iter, double likelihood_threshold,
            bool viterbi_like, double proba_threshold_factor, double mean_number_seq_err_thresh = INFINITY);

    bool load_genmodel();
    bool write2txt();
    bool readtxt();
    const Model_marginals get_marginals() const { return this->model_marginals; };
    /// The model as this object holds it, after inference if infer_model() ran. What
    /// igor::generation::legacy::SequenceGenerator takes to generate from the inferred model.
    const Model_Parms &get_model_parms() const { return this->model_parms; }

    //write alignments, load alignments

private:
    Model_Parms model_parms;
    Model_marginals model_marginals;
    std::map<size_t, std::shared_ptr<Counter>>
            counters_list; //Size_t is a unique identifier for the Counter(useful for adding them up)
    Model_marginals compute_marginals(std::list<std::string> sequences);
    Model_marginals compute_seq_marginals(std::string sequence);
    Model_marginals compute_seq_marginals(std::string sequence, std::list<std::list<std::string>> allowed_scenarios);
};

INFERENCE_EXPORT std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
get_best_aligns(
        const std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>> &,
        Gene_class);

} // namespace igor::inference::legacy
