/*
 * SequenceGenerator.cpp
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

#include <igor/Generation/Legacy/SequenceGenerator.h>
#include <igor/Alignment/Legacy/Aligner.h>

namespace igor::generation::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

using namespace std;

SequenceGenerator::SequenceGenerator(const Model_Parms &parms, const Model_marginals &marginals)
    : model_parms(parms), model_marginals(marginals)
{
}

SequenceGenerator::~SequenceGenerator() = default;

forward_list<pair<string, queue<queue<int>>>> SequenceGenerator::generate_sequences(int number_seq, bool generate_errors)
{

    queue<shared_ptr<Rec_Event>> model_queue = this->model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map = this->model_marginals.get_index_map(this->model_parms, model_queue);
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map =
            this->model_marginals.get_offsets_map(this->model_parms, model_queue);

    //Create seed for random generator
    //create a seed from timer
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point time = myclock::now();
    myclock::duration dur = (myclock::time_point::max)() - time;

    //Get a random seed
    uint64_t random_seed = draw_random_64bits_seed();
    //Instantiate random number generator
    mt19937_64 generator = mt19937_64(random_seed);
    forward_list<pair<string, queue<queue<int>>>> sequence_list = forward_list<pair<string, queue<queue<int>>>>();

    for (int seq = 0; seq != number_seq; ++seq) {
        pair<string, queue<queue<int>>> sequence =
                this->generate_unique_sequence(model_queue, index_map, offset_map, generator);
        if (generate_errors) {
            sequence.second.push(this->model_parms.get_err_rate_p()->generate_errors(sequence.first, generator));
        }
        sequence_list.push_front(sequence);
        if (seq % 1000 == 0) {
            //Output current progress to cerr
            show_progress_bar(cerr, seq / (double)number_seq, "Sequence generation", 50);
        }
    }
    close_progress_bar(cerr, "Sequence generation", 50);

    return sequence_list;
}
/*
 * Generate sequences in a memory efficient way
 */
void SequenceGenerator::generate_sequences(int number_seq, bool generate_errors, string filename_ind_seq,
                                  string filename_ind_real,
                                  list<pair<gen_seq_trans, shared_ptr<void>>>
                                          transform_func_and_data /*= list<pair<gen_seq_trans,shared_ptr<void>>>()*/,
                                  bool output_only_func /*= false*/, int seed /* =-1*/)
{
    ofstream outfile_ind_seq;
    ofstream outfile_ind_real;
    if (not output_only_func) {
        outfile_ind_seq.open(filename_ind_seq);
        outfile_ind_real.open(filename_ind_real);
    }
    string folder_path = filename_ind_seq.substr(0, filename_ind_seq.rfind("/") + 1); //Get the file path
    ofstream generation_infos_file(folder_path + "generation_info.out",
                                   fstream::out | fstream::app); //Opens the file in append mode

    //Create a header for the files
    queue<shared_ptr<Rec_Event>> model_queue = this->model_parms.get_model_queue();
    if (not output_only_func) {
        outfile_ind_seq << "seq_index;nt_sequence" << endl;
        outfile_ind_real << "seq_index";
        while (!model_queue.empty()) {
            outfile_ind_real << ";" << model_queue.front()->get_name();
            model_queue.pop();
        }
        outfile_ind_real << ";Errors" << endl;
    }
    model_queue = this->model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map = this->model_marginals.get_index_map(this->model_parms, model_queue);
    unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> offset_map =
            this->model_marginals.get_offsets_map(this->model_parms, model_queue);

    //Create seed for random generator
    //create a seed from timer if no seed was provided
    uint64_t random_seed;
    if (seed < 0) {
        //Get a random seed
        random_seed = draw_random_64bits_seed();
    } else {
        random_seed = seed;
    }
    clog << "Seed: " << random_seed << endl;
    //Instantiate random number generator
    mt19937_64 generator = mt19937_64(random_seed);

    chrono::system_clock::time_point begin_time = chrono::system_clock::now();
    std::time_t tt;
    tt = chrono::system_clock::to_time_t(begin_time);

    generation_infos_file << endl << "================================================================" << endl;
    generation_infos_file << "Generated sequences in file: " << filename_ind_seq << endl;
    generation_infos_file << "Generated sequences realizations in file: " << filename_ind_real << endl;
    generation_infos_file << "Date: " << ctime(&tt) << endl;
    generation_infos_file << "Number of sequences = " << number_seq << endl;
    generation_infos_file << "Generated with errors = " << generate_errors << endl;
    generation_infos_file << "Seed  = " << random_seed << endl;

    //Update events internal probas (e.g for dinucleotide ambiguous nucleotides)
    queue<shared_ptr<Rec_Event>> model_queue_copy = model_queue;
    while (not model_queue_copy.empty()) {
        model_queue_copy.front()->update_event_internal_probas(this->model_marginals.marginal_array_smart_p, index_map);
        model_queue_copy.pop();
    }

    for (size_t seq = 0; seq != number_seq; ++seq) {
        pair<string, queue<queue<int>>> sequence =
                this->generate_unique_sequence(model_queue, index_map, offset_map, generator, false);
        if (generate_errors) {
            sequence.second.push(this->model_parms.get_err_rate_p()->generate_errors(sequence.first, generator));
        }

        for (pair<gen_seq_trans, shared_ptr<void>> func_data_pair : transform_func_and_data) {
            func_data_pair.first(seq, sequence, func_data_pair.second);
        }

        if (not output_only_func) {
            outfile_ind_seq << seq << ";" << sequence.first << endl;
            outfile_ind_real << seq;
            queue<queue<int>> &realizations = sequence.second;
            while (!realizations.empty()) {
                outfile_ind_real << ";";
                queue<int> event_real = realizations.front();
                outfile_ind_real << "(";
                while (!event_real.empty()) {
                    outfile_ind_real << event_real.front();
                    event_real.pop();
                    if (!event_real.empty()) {
                        outfile_ind_real << ",";
                    }
                }
                outfile_ind_real << ")";
                realizations.pop();
            }
            outfile_ind_real << endl;
        }

        if (seq % 1000 == 0) {
            //Output current progress to cerr
            show_progress_bar(cerr, seq / (double)number_seq, "Sequence generation", 50);
        }
    }
    close_progress_bar(cerr, "Sequence generation", 50);
    return;
}

pair<string, queue<queue<int>>> SequenceGenerator::generate_unique_sequence(
        queue<shared_ptr<Rec_Event>> model_queue, unordered_map<Rec_Event_name, int> index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        mt19937_64 &generator, bool update_event_internal_proba /*= true*/)
{
    if (update_event_internal_proba) {
        queue<shared_ptr<Rec_Event>> model_queue_copy = model_queue;
        while (not model_queue_copy.empty()) {
            model_queue_copy.front()->update_event_internal_probas(this->model_marginals.marginal_array_smart_p,
                                                                   index_map);
            model_queue_copy.pop();
        }
    }

    unordered_map<Seq_type, string> *constructed_sequences_p = new unordered_map<Seq_type, string>;
    unordered_map<Seq_type, string> constructed_sequences = *constructed_sequences_p;
    queue<queue<int>> realizations;
    while (!model_queue.empty()) {
        realizations.push(model_queue.front()->draw_random_realization((this->model_marginals.marginal_array_smart_p),
                                                                       index_map, offset_map, constructed_sequences,
                                                                       generator));
        model_queue.pop();
    }
    //CAT strings
    string reconstructed_seq = constructed_sequences[V_gene_seq] + constructed_sequences[VJ_ins_seq]
            + constructed_sequences[VD_ins_seq] + constructed_sequences[D_gene_seq] + constructed_sequences[DJ_ins_seq]
            + constructed_sequences[J_gene_seq];
    delete constructed_sequences_p;
    return make_pair(reconstructed_seq, realizations);
}

void SequenceGenerator::write_seq2txt(string filename, forward_list<string> sequences)
{
    ofstream outfile(filename);
    for (forward_list<string>::const_iterator seq = sequences.begin(); seq != sequences.end(); ++seq) {
        outfile << (*seq) << endl;
    }
}

void SequenceGenerator::write_seq_real2txt(string filename_ind_seq, string filename_ind_real,
                                  forward_list<pair<string, queue<queue<int>>>> seq_and_realizations)
{
    ofstream outfile_ind_seq(filename_ind_seq);
    ofstream outfile_ind_real(filename_ind_real);

    //Create a header for the files
    outfile_ind_seq << "seq_index;nt_sequence" << endl;
    queue<shared_ptr<Rec_Event>> model_queue = this->model_parms.get_model_queue();
    outfile_ind_real << "Index";
    while (!model_queue.empty()) {
        outfile_ind_real << ";" << model_queue.front()->get_name();
        model_queue.pop();
    }
    outfile_ind_real << endl;

    size_t index = 0;
    for (forward_list<pair<string, queue<queue<int>>>>::const_iterator iter = seq_and_realizations.begin();
         iter != seq_and_realizations.end(); ++iter) {
        outfile_ind_seq << index << ";" << (*iter).first << endl;
        outfile_ind_real << index;
        queue<queue<int>> realizations = (*iter).second;
        while (!realizations.empty()) {
            outfile_ind_real << ";";
            queue<int> event_real = realizations.front();
            outfile_ind_real << "(";
            while (!event_real.empty()) {
                outfile_ind_real << event_real.front();
                event_real.pop();
                if (!event_real.empty()) {
                    outfile_ind_real << ",";
                }
            }
            outfile_ind_real << ")";
            realizations.pop();
        }
        outfile_ind_real << endl;
        index++;
    }
}

/**
 * FIXME for now the handling of non given anchors is very bad
 */
void output_CDR3_gen_data(size_t seq_index, std::pair<std::string, std::queue<std::queue<int>>> seq_and_real,
                          std::shared_ptr<void> func_data)
{
    gen_CDR3_data *func_data_cast = static_cast<gen_CDR3_data *>(func_data.get());

    tuple<string, size_t, size_t, string> *v_gene_anchors = nullptr;
    tuple<string, size_t, size_t, string> *j_gene_anchors = nullptr;

    size_t i = 0;
    while ((i != std::max(func_data_cast->v_event_queue_position, func_data_cast->j_event_queue_position) + 1)
           and (not seq_and_real.second.empty())) {
        if (i == func_data_cast->v_event_queue_position) {
            v_gene_anchors = &func_data_cast->v_anchors.at(seq_and_real.second.front().front());
            //There should be only one realization for v gene choice
        } else if (i == func_data_cast->j_event_queue_position) {
            j_gene_anchors = &func_data_cast->j_anchors.at(seq_and_real.second.front().front());
        }
        seq_and_real.second.pop();
        ++i;
    }

    //Compute the index of the last letter of the J anchor
    size_t tmp_index = seq_and_real.first.size() - get<2>(*j_gene_anchors) + get<1>(*j_gene_anchors) + 2;
    size_t tmp_v_index = get<1>(*v_gene_anchors);
    string nt_cdr3_seq = seq_and_real.first.substr(tmp_v_index, tmp_index - tmp_v_index + 1);

    /*	if( (nt_cdr3_seq.substr(0,3) == get<3>(*v_gene_anchors)) and (nt_cdr3_seq.substr(nt_cdr3_seq.size()-3,3) == get<3>(*j_gene_anchors))){
		func_data_cast->output_file<<nt_cdr3_seq<<",,"<<true<<",";
		if(nt_cdr3_seq.size()%3==0){
			func_data_cast->output_file<<true<<","<<endl;
		}
		else{
			func_data_cast->output_file<<false<<","<<endl;
		}
	}
	else{
		func_data_cast->output_file<<",,,"<<false<<","<<false<<","<<false<<endl;
	}
*/
    *func_data_cast->output_stream << seq_index;
    if (func_data_cast->output_nt_CDR3) {
        *func_data_cast->output_stream << "," << nt_cdr3_seq;
    }
    if (func_data_cast->output_anchors_found) {
        bool anchors_found = (nt_cdr3_seq.substr(0, 3) == get<3>(*v_gene_anchors))
                and (nt_cdr3_seq.substr(nt_cdr3_seq.size() - 3, 3) == get<3>(*j_gene_anchors));
        *func_data_cast->output_stream << "," << anchors_found;
    }
    if (func_data_cast->output_inframe) {
        bool is_inframe = nt_cdr3_seq.size() % 3 == 0;
        *func_data_cast->output_stream << "," << is_inframe;
    }
    if (func_data_cast->output_aa_CDR3) {
        //FIXME
        *func_data_cast->output_stream << "," << "";
    }
    if (func_data_cast->output_productive) {
        //FIXME
        *func_data_cast->output_stream << "," << "";
    }
    *func_data_cast->output_stream << endl;
}

//==============================================================================
// Fast Sequence Generation Implementation
//==============================================================================

void SequenceGenerator::generate_sequences_fast(size_t num_sequences, const string &seq_filename, const string &real_filename,
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
                          << (num_threads == 0 ? igor::generation::legacy::fast::get_optimal_thread_count() : num_threads) << endl;

    // Initialize fast generator if needed
    igor::generation::legacy::fast::FastGenerator &generator = get_fast_generator();

    // Configure generation
    igor::generation::legacy::fast::FastGeneratorConfig config;
    config.num_threads = num_threads;
    config.show_progress = show_progress;
    if (seed >= 0) {
        config.base_seed = static_cast<uint64_t>(seed);
    } else {
        config.base_seed = igor::generation::legacy::fast::draw_random_seed();
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

igor::generation::legacy::fast::FastGenerator &SequenceGenerator::get_fast_generator()
{
    if (!fast_generator_) {
        fast_generator_ = make_unique<igor::generation::legacy::fast::FastGenerator>();
        fast_generator_->initialize(model_parms, model_marginals);
    }
    return *fast_generator_;
}

} // namespace igor::generation::legacy
