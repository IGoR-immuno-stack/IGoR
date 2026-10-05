/*
 * SequenceGenerator.h
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

#include <forward_list>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <igor/Core/Legacy/Utils.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Generation/Legacy/FastGenerator.h>
#include <igor/Generation/Export.h>

namespace igor::generation::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

typedef void (*gen_seq_trans)(size_t, std::pair<std::string, std::queue<std::queue<int>>>, std::shared_ptr<void>);

/**
 * Hardcode a data structure for the function extracting CDR3s in generated sequences
 */
struct gen_CDR3_data
{
    std::map<int, std::tuple<std::string, size_t, size_t, std::string>> v_anchors;
    size_t v_event_queue_position;
    std::map<int, std::tuple<std::string, size_t, size_t, std::string>> j_anchors;
    size_t j_event_queue_position;
    std::shared_ptr<std::ostream> output_stream;
    //Some config booleans
    bool output_nt_CDR3 = true;
    bool output_anchors_found = true;
    bool output_inframe = true;
    //FIXME
    //For now do not output aa CDR3 stats
    bool output_aa_CDR3 = false;
    bool output_productive = false;

    gen_CDR3_data(const std::unordered_map<std::string, size_t> &v_anchors_indices,
                  const std::unordered_map<std::string, Event_realization> &v_reals, size_t v_event_pos,
                  const std::unordered_map<std::string, size_t> &j_anchors_indices,
                  const std::unordered_map<std::string, Event_realization> &j_reals, size_t j_event_pos,
                  std::shared_ptr<std::ostream> output_stream_ptr =
                          std::shared_ptr<std::ostream>(&std::cout, null_delete<std::ostream>()))
        : v_event_queue_position(v_event_pos), j_event_queue_position(j_event_pos), output_stream(output_stream_ptr)
    {

        //First get all V anchors
        this->v_anchors.clear();
        for (const std::pair<std::string, Event_realization> v_real : v_reals) {
            size_t v_anchor_index;
            if (v_anchors_indices.count(v_real.second.name) > 0) {
                v_anchor_index = v_anchors_indices.at(v_real.second.name);
                v_anchors.emplace(v_real.second.index,
                                  std::make_tuple(v_real.second.name, v_anchor_index, v_real.second.value_str.size(),
                                                  v_real.second.value_str.substr(v_anchor_index, 3)));
            } else {
                v_anchor_index = 0;
                v_anchors.emplace(
                        v_real.second.index,
                        std::make_tuple(v_real.second.name, v_anchor_index, v_real.second.value_str.size(), ""));
            }
            /*try{
				v_anchor_index = v_anchors_indices.at(v_real.name);
			}
			catch (std::exception& e) {
				std::cerr<<"Could not find "<<v_real.name<<" in the V genes anchors map"<<std::endl;
				throw e;
			}*/
            //v_anchors.emplace(v_real.second.index,std::make_tuple(v_real.second.name,v_anchor_index,v_real.second.value_str.size(),v_real.second.value_str.substr(v_anchor_index,3)));
        }

        //Now get all J anchors
        this->j_anchors.clear();
        for (const std::pair<std::string, Event_realization> j_real : j_reals) {
            size_t j_anchor_index;
            /*try{
				j_anchor_index = j_anchors_indices.at(j_real.name);
			}
			catch (std::exception& e) {
				std::cerr<<"Could not find "<<j_real.name<<" in the J genes anchors map"<<std::endl;
				throw e;
			}*/
            if (j_anchors_indices.count(j_real.second.name) > 0) {
                j_anchor_index = j_anchors_indices.at(j_real.second.name);
                j_anchors.emplace(j_real.second.index,
                                  std::make_tuple(j_real.second.name, j_anchor_index, j_real.second.value_str.size(),
                                                  j_real.second.value_str.substr(j_anchor_index, 3)));
            } else {
                j_anchor_index = std::string::npos;
                j_anchors.emplace(
                        j_real.second.index,
                        std::make_tuple(j_real.second.name, j_anchor_index, j_real.second.value_str.size(), ""));
            }
            //j_anchors.emplace(j_real.second.index,std::make_tuple(j_real.second.name,j_anchor_index,j_real.second.value_str.size(),j_real.second.value_str.substr(j_anchor_index,3)));
        }

        //Write output file header
        *output_stream.get() << "seq_index";
        if (output_nt_CDR3) {
            *output_stream.get() << ",nt_CDR3";
        }
        if (output_anchors_found) {
            *output_stream.get() << ",anchors_found";
        }
        if (output_inframe) {
            *output_stream.get() << ",is_inframe";
        }
        if (output_aa_CDR3) {
            *output_stream.get() << ",aa_CDR3";
        }
        if (output_productive) {
            *output_stream.get() << ",is_productive";
        }
        *output_stream.get() << std::endl;
    }
};

/**
 * \class SequenceGenerator SequenceGenerator.h
 * \brief Generates random sequences from a model.
 *
 * The generation half of GenModel, split off in step 1c of doc/LAYER_REFACTORING_PROPOSAL.md so
 * that Inference no longer depends on Generation. The bodies are GenModel's, unchanged. It holds
 * its own copy of the model, as GenModel does; to generate from an inferred model, construct it
 * from GenModel::get_model_parms() and GenModel::get_marginals().
 */
class GENERATION_EXPORT SequenceGenerator
{
public:
    SequenceGenerator(const Model_Parms &, const Model_marginals &);
    virtual ~SequenceGenerator();

    /**
     * \brief Fast parallel sequence generation (100x+ speedup).
     *
     * Uses precomputed CDFs, binary search/alias sampling, multi-threading,
     * and batched I/O for high-performance sequence generation.
     *
     * \param num_sequences Number of sequences to generate
     * \param seq_filename Output file for sequences
     * \param real_filename Output file for realizations
     * \param num_threads Number of threads (0 = auto-detect)
     * \param seed Random seed (-1 = random)
     * \param show_progress Show progress bar
     */
    void generate_sequences_fast(size_t num_sequences, const std::string &seq_filename,
                                 const std::string &real_filename, size_t num_threads = 0, int64_t seed = -1,
                                 bool show_progress = true);

    /**
     * \brief Get the fast generator instance (for advanced usage).
     *
     * Initializes the fast generator if not already done.
     */
    igor::generation::legacy::fast::FastGenerator &get_fast_generator();

    std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>> generate_sequences(int, bool);
    void generate_sequences(int, bool, std::string, std::string,
                            std::list<std::pair<gen_seq_trans, std::shared_ptr<void>>> =
                                    std::list<std::pair<gen_seq_trans, std::shared_ptr<void>>>(),
                            bool output_only_func = false, int = -1);
    void write_seq2txt(std::string, std::forward_list<std::string>);
    void write_seq_real2txt(std::string, std::string,
                            std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>>);

private:
    Model_Parms model_parms;
    Model_marginals model_marginals;
    std::unique_ptr<igor::generation::legacy::fast::FastGenerator> fast_generator_;
    std::pair<std::string, std::queue<std::queue<int>>> generate_unique_sequence(
            std::queue<std::shared_ptr<Rec_Event>>, std::unordered_map<Rec_Event_name, int>,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::mt19937_64 &, bool = true);
};

GENERATION_EXPORT void output_CDR3_gen_data(size_t, std::pair<std::string, std::queue<std::queue<int>>> seq_and_real,
                                            std::shared_ptr<void> func_data);

} // namespace igor::generation::legacy
