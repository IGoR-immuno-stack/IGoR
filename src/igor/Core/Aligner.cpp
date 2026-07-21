/*
 * Aligner.cpp
 *
 *  Created on: Feb 16, 2015
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

#include <igor/Core/Aligner.h>
#include <igor/Core/AlignerInternal.h>

#include <cctype>
#include <unordered_set>
#include <cmath>

using namespace std;

SwAlignmentMode default_sw_alignment_mode_for_gene(Gene_class gene)
{
    switch (gene) {
    case V_gene:
        return { false, false, true, false, false };
    case D_gene:
        return { true, true, true, true, false };
    case J_gene:
        return { false, false, false, true, true };
    case Undefined_gene:
        return { true, true, true, true, false };
    default:
        throw runtime_error("Erroneous gene class for alignments");
    }
}

SwAlignmentMode effective_sw_mode_for_dp(const SwDPConfig &config)
{
    SwAlignmentMode effective_mode = config.alignment_mode;
    if (effective_mode.reverse_sequences) {
        // DP boundary initialization is expressed as leading edges; when input
        // sequences are reversed, leading/trailing semantics must be mirrored.
        std::swap(effective_mode.data_leading_free, effective_mode.data_trailing_free);
        std::swap(effective_mode.genomic_leading_free, effective_mode.genomic_trailing_free);
    }
    return effective_mode;
}

Aligner::Aligner()
{
    // TODO Auto-generated constructor stub
}
/**
 * Constructor for the Aligner class
 * @sub_mat : substitution matrix
 * @gap_pen : sets the gap penalty (the gap penalty is linear)
 * @gene : Gene class of the gene aligned. V gene allows for deletions on the 3' side of the genomic template, J gene on the 5' , D gene and undefined allow deletion on both sides
 */
Aligner::Aligner(Matrix<double> sub_mat, int gap_pen, Gene_class gene)
    : substitution_matrix(sub_mat), gap_penalty(gap_pen), gene(gene)
{
    switch (gene) {
    case V_gene:
        break;
    case D_gene:
        break;
    case J_gene:
        break;
    case Undefined_gene:
        break;
    default:
        throw runtime_error("Erroneous gene class for alignments");
        break;
    }
};

Aligner::~Aligner()
{
    // TODO Auto-generated destructor stub
}
/*
 * This method reads sequences in a fasta file and return a vector of indexed sequences
 */
vector<pair<const int, const string>> read_fasta(const string &filename)
{
    //TODO Check for \r,\n\s stuff
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    string seq_str;
    string temp_str;
    int seq_count = -1;
    vector<pair<const int, const std::string>> sequence_vect;

    while (getline(infile, temp_str)) {
        if (temp_str[temp_str.size() - 1] == '\r') {
            temp_str.erase(temp_str.size() - 1);
        }
        if (temp_str[0] == '>') {

            if (seq_count > (-1)) {
                //Read sequences in upper case
                transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
                sequence_vect.push_back(pair<const int, const string>(seq_count, seq_str));
            }

            seq_str = string();
            seq_count++;
        } else {
            seq_str += temp_str;
        }
    }
    if (seq_count > (-1)) {
        //Read sequences in upper case
        transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
        sequence_vect.push_back(pair<const int, const string>(seq_count, seq_str));
    }

    return sequence_vect;
}

/*
 * This method reads genomic templates from fasta files
 */
vector<pair<string, string>> read_genomic_fasta(const string &filename)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    string seq_str;
    string seq_name_str;
    string temp_str = "";
    int seq_count = -1;
    vector<pair<string, std::string>> sequence_vect;

    while (getline(infile, temp_str)) {
        if (!temp_str.empty() && temp_str.back() == '\r') {
            temp_str.pop_back();
        }
        if (temp_str[0] == '>') {

            if (seq_count > (-1)) {
                //Convert all strings to upper case
                transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
                //Get rid of gaps inserted by IMGT to align the different genes
                string::iterator letter = seq_str.begin();
                while (letter != seq_str.end()) {
                    if ((*letter) == '.') {
                        letter = seq_str.erase(letter);
                    } else if ((*letter) == '\r') {
                        //TODO Remove line return in sequences files in a proper way
                        seq_str.erase(letter);
                        break;
                    } else {
                        ++letter;
                    }
                }
                seq_name_str.erase(seq_name_str.begin()); //Get rid of the '>' at the beginning of the line
                sequence_vect.push_back(pair<string, string>(seq_name_str, seq_str));
            }
            seq_name_str = temp_str;
            seq_str = string();
            seq_count++;
        } else {
            seq_str += temp_str;
        }
    }
    if (seq_count > (-1)) {
        //Convert all strings to upper case
        transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
        //Get rid of gaps inserted by IMGT to align the different genes
        string::iterator letter = seq_str.begin();
        while (letter != seq_str.end()) {

            if ((*letter) == '.') {
                letter = seq_str.erase(letter);
            } else if ((*letter) == '\r') {
                //TODO Remove line return in sequences files in a proper way
                seq_str.erase(letter);
                break;
            } else {
                ++letter;
            }
        }
        seq_name_str.erase(seq_name_str.begin());
        sequence_vect.push_back(pair<string, string>(seq_name_str, seq_str));
    }

    return sequence_vect;
}

/*
 * This method reads a file with one sequence per line and returns a vector of indexed sequences
 */
vector<pair<const int, const string>> read_txt(const string &filename)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    int seq_count = 0;
    vector<pair<const int, const std::string>> sequence_vect;
    string seq_str;

    while (getline(infile, seq_str)) {
        if (seq_str[seq_str.size() - 1] == '\r') {
            seq_str.erase(seq_str.size() - 1);
        }
        if (!seq_str.empty()) {
            transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
            sequence_vect.push_back(pair<const int, const string>(seq_count, seq_str));
            seq_count++;
        }
    }

    return sequence_vect;
}
/*
 * This methods reads a file containing at each line the index of the sequence a semicolon and the actual sequence
 * @return : vector of indexed sequences vector<pair<const int , const string>>
 */
vector<pair<const int, const string>> read_indexed_csv(const string &filename)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    string line_str;
    vector<pair<const int, const std::string>> sequence_vect;
    getline(infile, line_str);
    while (getline(infile, line_str)) {
        size_t semi_col_index = line_str.find(";");
        int index = stoi(line_str.substr(0, semi_col_index));
        string seq_str = line_str.substr(semi_col_index + 1, string::npos);
        transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
        sequence_vect.push_back(pair<const int, const string>(index, seq_str));
    }
    return sequence_vect;
}

/**
 * \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                int min_offset, int max_offset, bool rev_offset_frame /*=false*/)
{
    //Create a map of offset bounds and call align_seq overloads
    //This is not very elegant, however this function will probably not be called anymore except for aligning a single sequence.
    return align_seq(nt_seq, score_threshold, best_align_only, build_genomic_bounds_map(min_offset, max_offset),
                     rev_offset_frame);
}
/**
 * \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                int min_offset, int max_offset, set<string> restricted_genomic_list,
                                                bool rev_offset_frame /*=false*/)
{
    //Create a map of offset bounds and call align_seq overloads
    //This is not very elegant, however this function will probably not be called anymore except for aligning a single sequence.
    return align_seq(nt_seq, score_threshold, best_align_only, build_genomic_bounds_map(min_offset, max_offset),
                     restricted_genomic_list, rev_offset_frame);
}
/**
 * \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                bool best_gene_only, int min_offset, int max_offset,
                                                bool rev_offset_frame /*=false*/)
{
    //Create a map of offset bounds and call align_seq overloads
    //This is not very elegant, however this function will probably not be called anymore except for aligning a single sequence.
    return align_seq(nt_seq, score_threshold, best_align_only, best_gene_only,
                     build_genomic_bounds_map(min_offset, max_offset), rev_offset_frame);
}
/**
 * \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                bool best_gene_only, int min_offset, int max_offset,
                                                set<string> restricted_genomic_list, bool rev_offset_frame /*=false*/)
{
    //Create a map of offset bounds and call align_seq overloads
    //This is not very elegant, however this function will probably not be called anymore except for aligning a single sequence.
    return align_seq(nt_seq, score_threshold, best_align_only, best_gene_only,
                     build_genomic_bounds_map(min_offset, max_offset), restricted_genomic_list, rev_offset_frame);
}
/**
 *  \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                unordered_map<string, pair<int, int>> genomic_offset_bounds,
                                                bool rev_offset_frame /*=false*/)
{
    // Call the align seq function enforcing alignments for all genes to be passed
    return align_seq(nt_seq, score_threshold, best_align_only, false, genomic_offset_bounds,
                     rev_offset_frame /*=false*/);
}
/**
 *  \overload
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                unordered_map<string, pair<int, int>> genomic_offset_bounds,
                                                set<string> restricted_genomic_list, bool rev_offset_frame /*=false*/)
{
    // Call the align seq function enforcing alignments for all genes to be passed
    return align_seq(nt_seq, score_threshold, best_align_only, false, genomic_offset_bounds, restricted_genomic_list,
                     rev_offset_frame /*=false*/);
}

forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                bool best_gene_only,
                                                unordered_map<string, pair<int, int>> genomic_offset_bounds,
                                                bool rev_offset_frame /*=false*/)
{
    set<string> all_genomic_names;
    for (pair<string, string> genomic_template : this->nt_genomic_sequences) {
        all_genomic_names.emplace(genomic_template.first);
    }
    return align_seq(nt_seq, score_threshold, best_align_only, best_gene_only, genomic_offset_bounds, all_genomic_names,
                     rev_offset_frame /*=false*/);
}

/**
 * \brief A function performing alignment of all genomic templates against a given sequence
 * \author Q.Marcou, M.Puelma Touzel
 * \version 1.2.1
 *
 * \param [in] nt_seq the nucleotide sequence to study
 * \param [in] score_threshold The SW alignment score threshold to record an alignment
 * \param [in] best_align_only Only retain the best alignment for each genomic template.
 * \param [in] best_gene_only Only retain the best gene/allele candidate (or best candidates if several have the same highest score).
 * \param [in] genomic_offset_bounds A hash map containing offsets lower and upper bounds for each genomic template. Keys of the map are the genomic templates names.
 * \param [in] restricted_genomic_list A set containing the names of the genes that should be aligned to the sequence.
 * \param [in] rev_offset_frame Are offsets bounds given reversed? (offset defined based on the last sequence nt instead of the first). Default is false.
 *
 * Call the SW alignment function for every genomic template aligning them against one target sequence.
 * There is possibility to pass different offset bounds for different genomic templates.
 * There is also a possibility to pass these offsets reversed (i.e defined from the last nucleotide of the target) in case it is more handy (e.g alignement of CDR3 sequences or J/C primer sequencing).
 *
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq, double score_threshold, bool best_align_only,
                                                bool best_gene_only,
                                                unordered_map<string, pair<int, int>> genomic_offset_bounds,
                                                set<string> restricted_genomic_list, bool rev_offset_frame /*=false*/)
{
    int min_offset;
    int max_offset;
    Int_Str int_seq = nt2int(nt_seq);
    size_t seqlen = int_seq.size();
    forward_list<Alignment_data> alignment_list; // = *(new forward_list<Alignment_data>());

    // Substitution matrix, gap penalty and alignment mode are invariant across all genomic
    // templates for this gene; only min_offset/max_offset vary per template below. Hoisting the
    // config out of the loop avoids re-deep-copying substitution_matrix on every template.
    SwDPConfig config{ score_threshold,
                       best_align_only,
                       /*min_offset=*/0,
                       /*max_offset=*/0,
                       this->substitution_matrix,
                       this->gap_penalty,
                       default_sw_alignment_mode_for_gene(gene),
                       this->enable_extension_ };

    for (forward_list<pair<string, Int_Str>>::const_iterator iter = int_genomic_sequences.begin();
         iter != int_genomic_sequences.end(); ++iter) {
        //If the gene must be aligned
        if (restricted_genomic_list.count((*iter).first) > 0) {
            // Extract min and max offset information from the offset bounds map
            try {
                min_offset = genomic_offset_bounds.at((*iter).first).first;
                max_offset = genomic_offset_bounds.at((*iter).first).second;
            } catch (exception &e) {
                cerr << "Exception caught trying to fetch template specific offset bounds in Aligner::align_seq"
                     << endl;
                throw runtime_error("Missing genomic offset bounds for genomic template \"" + iter->first + "\"");
            }

            // Reverse the offset if necessary (e.g for J CDR3 alignment or sequencing from J primer)
            min_offset += (rev_offset_frame) ? seqlen - 1
                                             : 0; //seqlen-1 correspond to the index of the last nt of the sequence
            max_offset += (rev_offset_frame) ? seqlen - 1 : 0;

            list<pair<int, Alignment_data>> alignments;
            try {
                config.min_offset = min_offset;
                config.max_offset = max_offset;
                alignments = sw_align(int_seq, (*iter).second, best_align_only, config);
            } catch (exception &e) {
                cerr << endl;
                cerr << "Exception caught calling sw_align() on genomic template:" << (*iter).first << endl;
                throw e;
            }
            //TODO quick and dirty fix for D genes alignments
            //alignment.second.gene_name = (*iter).first;
            //alignment_list.push_front(alignment.second);
            for (list<pair<int, Alignment_data>>::iterator jiter = alignments.begin(); jiter != alignments.end();
                 ++jiter) {
                (*jiter).second.gene_name = (*iter).first;
                alignment_list.push_front((*jiter).second);
            }
        }
    }

    if (best_gene_only) {
        //Only return alignments for the gene with best alignment
        alignment_list = extract_best_gene_alignments(alignment_list);
    }

    return alignment_list;
}

/*
 * Align sequences and hold them in memory
 */
unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only)
{
    unordered_map<int, forward_list<Alignment_data>> alignment_map =
            align_seqs(sequence_list, score_threshold, best_align_only, false, INT16_MIN, INT16_MAX);
    return alignment_map;
}

/*
 * Align sequences and hold them in memory
 */
unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only,
                    bool best_gene_only)
{
    unordered_map<int, forward_list<Alignment_data>> alignment_map =
            align_seqs(sequence_list, score_threshold, best_align_only, best_gene_only, INT16_MIN, INT16_MAX);
    return alignment_map;
}

/*
 * \brief A function performing alignment of all genomic templates against all provided sequences. Alignments are stored in memory.
 */
unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only,
                    int min_offset, int max_offset, bool rev_offset_frame /*=false*/)
{
    unordered_map<int, forward_list<Alignment_data>> alignment_map =
            align_seqs(sequence_list, score_threshold, best_align_only, false, min_offset, max_offset);
    return alignment_map;
}

unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only,
                    bool best_gene_only, int min_offset, int max_offset, bool rev_offset_frame /*=false*/)
{
    unordered_map<int, forward_list<Alignment_data>> alignment_map =
            align_seqs(sequence_list, score_threshold, best_align_only, best_gene_only,
                       build_genomic_bounds_map(min_offset, max_offset), rev_offset_frame);
    return alignment_map;
}

unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only,
                    unordered_map<string, pair<int, int>> genomic_offset_bounds, bool rev_offset_frame /*=false*/)
{
    unordered_map<int, forward_list<Alignment_data>> alignment_map =
            align_seqs(sequence_list, score_threshold, best_align_only, false, genomic_offset_bounds, rev_offset_frame);
    return alignment_map;
}

/*
 * \brief A function performing alignment of all genomic templates against all provided sequences. Alignments are stored in memory.
 */
unordered_map<int, forward_list<Alignment_data>>
Aligner::align_seqs(vector<pair<const int, const string>> sequence_list, double score_threshold, bool best_align_only,
                    bool best_gene_only, unordered_map<string, pair<int, int>> genomic_offset_bounds,
                    bool rev_offset_frame /*=false*/)
{
    unordered_map<int, forward_list<Alignment_data>>
            alignment_map; //= *(new unordered_map<int,forward_list<Alignment_data>>);

    int processed_seq_number = 0;
    double total_number_seqs = sequence_list.size(); //Use a double for float division afterwards

    /*
 * Declaring parellel loop using OpenMP 4.0 standards
	#pragma omp declare reduction (merge:unordered_map<int,forward_list<Alignment_data>>:omp_out.insert(omp_in.begin(),omp_in.end()))
	#pragma omp parallel for schedule(dynamic) reduction(merge:alignment_map) shared(processed_seq_number)
*/
    const auto n_seqs =
            sequence_list.size(); // Required for OpenMP with MSVC: can't use iterator for loop, need to use indexes

//Declare parallel loop using OpenMP 3.1 standards
#pragma omp parallel for schedule(dynamic) shared(processed_seq_number, alignment_map) //num_threads(1)
    for (auto i = 0; i < n_seqs; ++i) {

        const auto &seq_pair = sequence_list[i]; // au lieu de *seq_it

        forward_list<Alignment_data> seq_alignments =
                align_seq(seq_pair.second, score_threshold, best_align_only, best_gene_only, genomic_offset_bounds,
                          rev_offset_frame);
#pragma omp critical(emplace_seq_alignments)
        {
            alignment_map.emplace(seq_pair.first, seq_alignments);
            //cout<<"Seq "<<processed_seq_number<<" processed"<<endl;
            ++processed_seq_number;
        }

        if (processed_seq_number % 50 == 0) {
//Output current progress to cerr
#pragma omp critical(show_progress_align)
            {
                show_progress_bar(cerr, processed_seq_number / total_number_seqs, to_string(this->gene) + " alignments",
                                  50);
            }
        }
    }
    close_progress_bar(cerr, to_string(this->gene) + " alignments", 50);
    return alignment_map;
}
/**
 * \overload
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only, int min_offset, int max_offset, bool rev_offset_frame /*=false*/)
{
    return this->align_seqs(filename, sequence_list, score_threshold, best_align_only, false,
                            build_genomic_bounds_map(min_offset, max_offset), rev_offset_frame);
}

/**
 * \overload
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only, bool best_gene_only, int min_offset, int max_offset,
                         bool rev_offset_frame /*=false*/)
{
    return this->align_seqs(filename, sequence_list, score_threshold, best_align_only, best_gene_only,
                            build_genomic_bounds_map(min_offset, max_offset), rev_offset_frame);
}

/**
 * \overload
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only, unordered_map<string, pair<int, int>> genomic_offset_bounds,
                         bool rev_offset_frame /*=false*/)
{
    return this->align_seqs(filename, sequence_list, score_threshold, best_align_only, false, genomic_offset_bounds,
                            rev_offset_frame);
}

/**
 * \brief A function performing alignment of all genomic templates against all provided sequences. Output on file.
 * \author Q.Marcou, M.Puelma Touzel
 * \version 1.2.0
 *
 * \param [in] filename Path and filename for the ouput alignment file
 * \param [in] sequence_list A forward list containing pairs of nt sequence and the corresponding index
 * \param [in] nt_seq the nucleotide sequence to study
 * \param [in] score_threshold The SW alignment score threshold to record an alignment
 * \param [in] best_align_only Only retain the best alignment for each genomic template.
 * \param [in] best_gene_only Only retain the best gene/allele candidate (or best candidates if several have the same highest score).
 * \param [in] genomic_offset_bounds A hash map containing offsets lower and upper bounds for each genomic template. Keys of the map are the genomic templates names.
 * \param [in] rev_offset_frame Are offsets bounds given reversed? (offset defined based on the last sequence nt instead of the first). Default is false.
 *
 * Call the SW alignment function for every genomic template aligning them against all target sequences.
 * Alignments are all written on disk on the fly to avoid memory issues.
 * A summary file containing all alignments parameters and relevant information is created/appended in the directory.
 * There is possibility to pass different offset bounds for different genomic templates.
 * There is also a possibility to pass these offsets reversed (i.e defined from the last nucleotide of the target) in case it is more handy (e.g alignement of CDR3 sequences or J/C primer sequencing).
 *
 * \bug Summary file creation might not work on Windows systems
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only, bool best_gene_only,
                         unordered_map<string, pair<int, int>> genomic_offset_bounds, bool rev_offset_frame /*=false*/)
{

    unordered_map<int, forward_list<Alignment_data>>
            alignment_map; //= *(new unordered_map<int,forward_list<Alignment_data>>);

    string folder_path = filename.substr(0, filename.rfind("/") + 1); //Get the file path
    ofstream align_infos_file(folder_path + "aligns_info.out",
                              fstream::out | fstream::app); //Opens the file in append mode

    //Check if all templates have the same min and max offsets and compute min and max over all of them.
    tuple<bool, int, int> min_max_offsets = extract_min_max_genomic_templates_offsets(genomic_offset_bounds);

    // Start chronometer and get dates and time
    chrono::system_clock::time_point begin_time = chrono::system_clock::now();
    std::time_t tt;
    tt = chrono::system_clock::to_time_t(begin_time);

    align_infos_file << endl << "================================================================" << endl;
    align_infos_file << "Alignments in file: " << filename << endl;
    align_infos_file << "Date: " << ctime(&tt) << endl;
    align_infos_file << "Score threshold = " << score_threshold << endl;
    align_infos_file << "Best alignement per gene/allele only = " << best_align_only << endl;
    align_infos_file << "Best gene/allele candidate only = " << best_gene_only << endl;
    align_infos_file << "Min Offset = " << get<1>(min_max_offsets) << endl;
    align_infos_file << "Max Offset = " << get<2>(min_max_offsets) << endl;
    align_infos_file << "Using template specific offsets = " << get<0>(min_max_offsets) << endl;
    align_infos_file << "Using reversed offsets = " << rev_offset_frame << endl;
    align_infos_file << "Gap penalty = " << this->gap_penalty << endl;
    align_infos_file << "Substitution matrix:" << endl;
    align_infos_file << this->substitution_matrix << endl;
    align_infos_file << sequence_list.size() << " sequences processed in ";

    ofstream outfile(filename);
    outfile << "seq_index" << ";" << "gene_name" << ";" << "score" << ";" << "offset" << ";" << "insertions" << ";"
            << "deletions" << ";" << "mismatches" << ";" << "length" << ";5_p_align_offset;3_p_align_offset" << endl;

    int processed_seq_number = 0;
    double total_number_seqs = sequence_list.size(); //Use a double for float division afterwards

    /*
 * Declaring parallel loop using OpenMP 4.0 standards
	#pragma omp declare reduction (merge:unordered_map<int,forward_list<Alignment_data>>:omp_out.insert(omp_in.begin(),omp_in.end()))
	#pragma omp parallel for schedule(dynamic) reduction(merge:alignment_map) shared(processed_seq_number)
*/

    const auto n_seqs =
            sequence_list.size(); // Required for OpenMP with MSVC: can't use iterator for loop, need to use indexes

//Declare parallel loop using OpenMP 3.1 standards
#pragma omp parallel for schedule(dynamic) shared(processed_seq_number, alignment_map) //num_threads(1)
    for (auto i = 0; i < n_seqs; ++i) {

        const auto &seq_pair = sequence_list[i]; // au lieu de *seq_it
        try {
            forward_list<Alignment_data> seq_alignments =
                    align_seq(seq_pair.second, score_threshold, best_align_only, best_gene_only, genomic_offset_bounds,
                              rev_offset_frame);

#pragma omp critical(emplace_seq_alignments)
            {
                write_single_seq_alignment(outfile, seq_pair.first, seq_alignments);
                //cout<<"Seq "<<processed_seq_number<<" processed"<<endl;
                ++processed_seq_number;
            }

            //Output current progress to cerr
            if (processed_seq_number % 50 == 0) {
#pragma omp critical(show_progress_align)
                {
                    show_progress_bar(cerr, processed_seq_number / total_number_seqs,
                                      to_string(this->gene) + " alignments", 50);
                }
            }
        } catch (exception &except) {
            cerr << endl;
            cerr << "Exception caught calling align_seq() on sequence:" << endl;
            cerr << seq_pair.first << ";" << seq_pair.second << endl;
            cerr << endl;
            cerr << "Throwing exception now..." << endl << endl;
            cerr << except.what() << endl;
            throw except;
        }
    }
    close_progress_bar(cerr, to_string(this->gene) + " alignments", 50);

    chrono::duration<double> elapsed_time = chrono::system_clock::now() - begin_time;
    align_infos_file << elapsed_time.count() << " seconds" << endl;
}

/*
 * Align sequences and write them on disk on the fly (avoids memory issues)
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only)
{
    this->align_seqs(filename, sequence_list, score_threshold, best_align_only, INT16_MIN, INT16_MAX);
}
/*
 * Align sequences and write them on disk on the fly (avoids memory issues)
 */
void Aligner::align_seqs(string filename, vector<pair<const int, const string>> sequence_list, double score_threshold,
                         bool best_align_only, bool best_gene_only)
{
    this->align_seqs(filename, sequence_list, score_threshold, best_align_only, best_gene_only, INT16_MIN, INT16_MAX);
}

/**
 * \brief A small function to automatically build a hashmap containing genomic offset bounds from fixed bounds over genomic templates.
 */
std::unordered_map<std::string, std::pair<int, int>> Aligner::build_genomic_bounds_map(int min_offset,
                                                                                       int max_offset) const
{
    unordered_map<string, pair<int, int>> genomic_offset_bounds;
    for (forward_list<pair<string, Int_Str>>::const_iterator iter = this->int_genomic_sequences.begin();
         iter != this->int_genomic_sequences.end(); ++iter) {
        genomic_offset_bounds.emplace(iter->first, make_pair(min_offset, max_offset));
    }
    return genomic_offset_bounds;
}

/*
 * Writes the indexed sequences as semicolon separated files with 2 fields:
 * @seq_index
 * @sequence
 */
void write_indexed_seq_csv(const string &filename, const vector<pair<const int, const string>> &indexed_seq_list)
{
    ofstream outfile(filename);
    outfile << "seq_index" << ";" << "sequence" << endl;
    for (vector<pair<const int, const string>>::const_iterator iter = indexed_seq_list.begin();
         iter != indexed_seq_list.end(); ++iter) {
        outfile << (*iter).first << ";" << (*iter).second << endl;
    }
}
/*
 * Writes the alignment in a semicolon separated files with 5 fields:
 * @seq_index
 * @gene_name
 * @offset
 * @insertions: list of int coma separated surrounded by curly braces
 * @deletions: list of int coma separated surrounded by curly braces
 */
void Aligner::write_alignments_seq_csv(string filename,
                                       unordered_map<int, forward_list<Alignment_data>> indexed_alignments)
{
    ofstream outfile(filename);
    outfile << "seq_index" << ";" << "gene_name" << ";" << "score" << ";" << "offset" << ";" << "insertions" << ";"
            << "deletions" << ";" << "mismatches" << ";" << "length" << endl;

    for (unordered_map<int, forward_list<Alignment_data>>::const_iterator iter = indexed_alignments.begin();
         iter != indexed_alignments.end(); ++iter) {
        write_single_seq_alignment(outfile, (*iter).first, (*iter).second);
    }
}

/*
 * This method writes the alignments for one sequence in the given stream
 */
void write_single_seq_alignment(ofstream &outfile, int seq_index, forward_list<Alignment_data> seq_alignments)
{
    for (forward_list<Alignment_data>::const_iterator jiter = seq_alignments.begin(); jiter != seq_alignments.end();
         ++jiter) {
        outfile << seq_index << ";" << (*jiter).gene_name << ";" << (*jiter).score << ";" << (*jiter).offset << ";{";
        const vector<int>& insertions = (*jiter).get_all_insertions();
        for (size_t i = 0; i < insertions.size(); ++i) {
            if (i == 0) {
                outfile << insertions[i];
            } else {
                outfile << "," << insertions[i];
            }
        }
        outfile << "};{";
        const vector<int>& deletions = (*jiter).get_all_deletions();
        for (size_t i = 0; i < deletions.size(); ++i) {
            if (i == 0) {
                outfile << deletions[i];
            } else {
                outfile << "," << deletions[i];
            }
        }
        outfile << "};{"; //<<endl;
        const vector<int>& mismatches = (*jiter).get_all_mismatches();
        for (size_t i = 0; i < mismatches.size(); ++i) {
            if (i == 0) {
                outfile << mismatches[i];
            } else {
                outfile << "," << mismatches[i];
            }
        }
        outfile << "};" << (*jiter).align_length << ";" << (*jiter).five_p_offset << ";" << (*jiter).three_p_offset
                << endl;
    }
}

std::vector<std::pair<int, char>> parse_cigar(const std::string &cigar)
{
    std::vector<std::pair<int, char>> ops;
    if (cigar.empty()) {
        throw invalid_argument("empty CIGAR string");
    }
    size_t i = 0;
    while (i < cigar.size()) {
        if (!isdigit(static_cast<unsigned char>(cigar[i]))) {
            throw invalid_argument("CIGAR count expected");
        }
        long count = 0;
        while (i < cigar.size() && isdigit(static_cast<unsigned char>(cigar[i]))) {
            count = count * 10 + (cigar[i] - '0');
            ++i;
        }
        if (count <= 0 || i == cigar.size()) {
            throw invalid_argument("invalid CIGAR count");
        }
        char op = cigar[i++];
        const string valid_ops = "=XIDNSHPM";
        if (valid_ops.find(op) == string::npos) {
            throw invalid_argument("invalid CIGAR operation");
        }
        ops.push_back(make_pair(static_cast<int>(count), op));
    }
    return ops;
}

static void append_cigar_op(vector<pair<int, char>> &ops, char op)
{
    if (!ops.empty() && ops.back().second == op) {
        ++ops.back().first;
    } else {
        ops.push_back(make_pair(1, op));
    }
}

static string cigar_ops_to_string(const vector<pair<int, char>> &ops)
{
    string cigar_out;
    for (const auto &op : ops) {
        cigar_out += to_string(op.first);
        cigar_out += op.second;
    }
    return cigar_out;
}

static void append_cigar_run(vector<pair<int, char>> &ops, int count, char op)
{
    if (count <= 0)
        return;
    if (!ops.empty() && ops.back().second == op)
        ops.back().first += count;
    else
        ops.push_back(make_pair(count, op));
}





/**
 * Convert Alignment_data to core CIGAR string (1-parameter overload).
 * This function produces a CIGAR string for just the core alignment region
 * between five_p_offset and three_p_offset, using =, X, D, I operators.
 */
std::string alignment_data_to_core_cigar(const Alignment_data &aln)
{
    vector<int> insertions = aln.get_all_insertions();
    vector<int> deletions = aln.get_all_deletions();
    sort(insertions.begin(), insertions.end());
    sort(deletions.begin(), deletions.end());
    unordered_set<int> mismatches(aln.get_all_mismatches().begin(), aln.get_all_mismatches().end());

    size_t ins_i = 0;
    size_t del_i = 0;
    int t = static_cast<int>(aln.five_p_offset);
    int g = static_cast<int>(aln.five_p_offset) - aln.offset;
    const int t_end = static_cast<int>(aln.three_p_offset);
    vector<pair<int, char>> ops;

    while (t <= t_end) {
        while (del_i < deletions.size() && deletions[del_i] == g) {
            append_cigar_op(ops, 'D');
            ++g;
            ++del_i;
        }
        if (ins_i < insertions.size() && insertions[ins_i] == t) {
            append_cigar_op(ops, 'I');
            ++t;
            ++ins_i;
        } else {
            append_cigar_op(ops, mismatches.count(t) ? 'X' : '=');
            ++t;
            ++g;
        }
    }
    while (del_i < deletions.size() && deletions[del_i] == g) {
        append_cigar_op(ops, 'D');
        ++g;
        ++del_i;
    }

    return cigar_ops_to_string(ops);
}

/**
 * Convert Alignment_data to core CIGAR string (3-parameter version).
 * This was renamed from alignment_data_to_cigar_full_span as requested.
 * This function produces a CIGAR string that includes the full span with
 * AIRR-compliant N/S operators for end gaps.
 */
std::string alignment_data_to_core_cigar(const Alignment_data &aln, size_t sequence_length, size_t germline_length)
{
    int t = static_cast<int>(aln.five_p_offset);
    int g = static_cast<int>(aln.five_p_offset) - aln.offset;
    vector<pair<int, char>> ops;
    // Use AIRR-standard operators: N for reference-only gaps, S for query-only gaps
    append_cigar_run(ops, g, 'N');  // Leading reference gaps
    append_cigar_run(ops, t, 'S');  // Leading query gaps
    
    // Reuse the core version for the core region
    for (const auto &entry : parse_cigar(alignment_data_to_core_cigar(aln))) {
        append_cigar_run(ops, entry.first, entry.second);
        switch (entry.second) {
        case '=':
        case 'X':
        case 'M':
            t += entry.first;
            g += entry.first;
            break;
        case 'I':
        case 'S':
            t += entry.first;
            break;
        case 'D':
        case 'N':
            g += entry.first;
            break;
        default:
            break;
        }
    }
    // Use AIRR-standard operators: N for reference-only gaps, S for query-only gaps
    append_cigar_run(ops, static_cast<int>(germline_length) - g, 'N');  // Trailing reference gaps
    append_cigar_run(ops, static_cast<int>(sequence_length) - t, 'S');  // Trailing query gaps
    return cigar_ops_to_string(ops);
}

/**
 * Convert Alignment_data to extended CIGAR string including extended alignment regions.
 * This function creates a comprehensive CIGAR representation that includes:
 * - Leading reference gaps as N operators
 * - Leading query gaps as S operators  
 * - Core alignment with =, X, D, I operators (extended to include extended mismatches)
 * - Trailing reference gaps as N operators
 * - Trailing query gaps as S operators
 * 
 * Unlike alignment_data_to_core_cigar which only processes the core alignment region,
 * this function extends the alignment to include ALL mismatch positions from alignment_data.mismatches,
 * representing extended mismatches as X operators rather than as gaps.
 */
std::string alignment_data_to_extended_cigar(const Alignment_data &aln, size_t sequence_length, size_t germline_length)
{
    Alignment_data aln_copy = aln;
    // Update offset values to extend alignment bounds
    aln_copy.five_p_offset = max(0, aln.offset);
    size_t n_del = aln.get_all_deletions().size();
    size_t n_ins = aln.get_all_insertions().size();
    aln_copy.three_p_offset = min(aln_copy.offset + germline_length - 1 + n_ins - n_del, sequence_length - 1);
    aln_copy.align_length = aln_copy.three_p_offset - aln_copy.five_p_offset + n_del;
    return alignment_data_to_core_cigar(aln_copy, sequence_length, germline_length);
}

Alignment_data alignment_data_from_cigar(const std::string &gene_name, const std::string &cigar, int seq_start_1based,
                                         int seq_end_1based, int ref_start_1based, int /*ref_end_1based*/, double score)
{
    int offset = seq_start_1based - ref_start_1based;
    size_t five_p_offset = static_cast<size_t>(seq_start_1based - 1);
    size_t three_p_offset = static_cast<size_t>(seq_end_1based - 1);
    int t = seq_start_1based - 1;
    int g = ref_start_1based - 1;
    size_t align_length = 0;
    vector<int> insertions;
    vector<int> deletions;
    vector<int> mismatches;

    for (const auto &entry : parse_cigar(cigar)) {
        int count = entry.first;
        char op = entry.second;
        for (int i = 0; i < count; ++i) {
            switch (op) {
            case '=':
            case 'M':
                ++t;
                ++g;
                break;
            case 'X':
                mismatches.push_back(t);
                ++t;
                ++g;
                break;
            case 'I':
                insertions.push_back(t++);
                break;
            case 'D':
                deletions.push_back(g++);
                break;
            case 'N':
                deletions.push_back(g++);
                break;
            case 'S':
                ++t;
                break;
            case 'H':
            case 'P':
                break;
            default:
                throw invalid_argument("unsupported CIGAR operation");
            }
            ++align_length;
        }
    }
    return Alignment_data(gene_name, offset, five_p_offset, three_p_offset, align_length, insertions, deletions,
                          mismatches, score, 0, 0);
}

/**
 * Construct Alignment_data from core CIGAR and extended CIGAR strings.
 * 
 * Rules:
 * - offset: if extended CIGAR has leading N, offset = - (number of leading N)
 *           if extended CIGAR has leading S, offset = + (number of leading S)
 *           there should never be both leading N and S (invalid extended CIGAR)
 * - five_p_offset = number of leading "S" in the core CIGAR
 * - three_p_offset = five_p_offset + (number of =, X, M, I in core CIGAR) - 1
 * - align_length = count of =, X, M, I, D in core CIGAR
 * 
 * The extended_cigar is used to extract ALL mismatches (X operators).
 * The core_cigar is used to extract insertions (I) and deletions (D).
 * 
 * Validation: Extended CIGAR should not have both leading N and S,
 * and should not have both trailing N and S.
 */
Alignment_data alignment_data_from_cigar_and_extended(const std::string &gene_name, const std::string &core_cigar,
                                                      const std::string &extended_cigar, double score)
{
    // Parse extended CIGAR once
    const auto ext_entries = parse_cigar(extended_cigar);

    // Validate and compute offset from extended CIGAR
    int offset = 0;
    bool has_leading_n = false;
    bool has_leading_s = false;
    bool has_trailing_n = false;
    bool has_trailing_s = false;
    bool seen_non_gap = false;

    // Find last non-gap index and check leading gaps
    size_t last_non_gap_index = 0;

    // Merged pass: compute offset, track gaps, and populate position vectors
    std::vector<int> mismatches;
    std::vector<int> insertions;
    std::vector<int> deletions;
    // Keep track of both query and reference implied indices/positions
    int query_pos_ext = 0;
    int ref_pos = 0;

    for (size_t i = 0; i < ext_entries.size(); ++i) {
        int count = ext_entries[i].first;
        char op = ext_entries[i].second;

        // Check for leading gaps
        if (!seen_non_gap) {
            if (op == 'N') {
                has_leading_n = true;
                offset -= count;
            } else if (op == 'S') {
                has_leading_s = true;
                offset += count;
            } else {
                seen_non_gap = true;
            }
        }

        // Track last non-gap entry
        if (op != 'N' && op != 'S') {
            last_non_gap_index = i;
            seen_non_gap = true;
        }

        // Populate position vectors
        for (int j = 0; j < count; ++j) {
            switch (op) {
            case '=':
            case 'X':
            case 'M':
                if (op == 'X') {
                    mismatches.push_back(query_pos_ext);
                }
                query_pos_ext++;
                ref_pos++;
                break;
            case 'I':
                insertions.push_back(query_pos_ext);
                query_pos_ext++;
                break;
            case 'D':
                deletions.push_back(ref_pos);
                ref_pos++;
                break;
            case 'S':
                query_pos_ext++;
                break;
            case 'N':
                ref_pos++;
                break;
            default:
                break;
            }
        }
    }

    // Validate leading gaps
    if (has_leading_n && has_leading_s) {
        throw std::invalid_argument("Invalid extended CIGAR: has both leading N and S");
    }

    // Check trailing gaps (after last non-gap entry)
    for (size_t i = last_non_gap_index + 1; i < ext_entries.size(); ++i) {
        char op = ext_entries[i].second;
        if (op == 'N') {
            has_trailing_n = true;
        } else if (op == 'S') {
            has_trailing_s = true;
        }
    }

    // Validate trailing gaps
    if (has_trailing_n && has_trailing_s) {
        throw std::invalid_argument("Invalid extended CIGAR: has both trailing N and S");
    }

    // Parse core CIGAR once
    const auto core_entries = parse_cigar(core_cigar);

    // Parse core CIGAR to get alignment bounds, length, insertions, and deletions
    size_t five_p_offset = 0;
    size_t align_length = 0;
    size_t query_pos_core = 0; // counts =, X, M, I only (for three_p_offset)
    bool first_non_s_or_n = false;

    for (const auto &entry : core_entries) {
        int count = entry.first;
        char op = entry.second;

        // Count leading S for five_p_offset, account for possible leading N too.
        if(!first_non_s_or_n){
            switch (op)
            {
            case 'S':
                five_p_offset += count;
                break;
            case 'N':
                // Do nothing
                break;
            default:
                // First non S or N encountered
                first_non_s_or_n = true;
                break;
            }
        }

        for (int i = 0; i < count; ++i) {
            switch (op) {
            case '=':
            case 'X':
            case 'M':
                query_pos_core++;
                align_length++;
                break;
            case 'I':
                query_pos_core++;
                align_length++;
                break;
            case 'D':
                align_length++;
                break;
            case 'S':
                // S doesn't count towards alignment positions
                // (leading S are counted in five_p_offset, trailing S are ignored)
                break;
            case 'N':
                // N doesn't count in core alignment
                break;
            default:
                break;
            }
        }
    }

    // Compute three_p_offset
    // three_p_offset = five_p_offset + (number of =, X, M, I) - 1
    size_t three_p_offset = five_p_offset + query_pos_core - 1;

    // Sort mismatches for consistency
    std::sort(mismatches.begin(), mismatches.end());

    return Alignment_data(gene_name, offset, five_p_offset, three_p_offset, align_length,
                          insertions, deletions, mismatches, score, 0, 0);
}

int alignment_data_sequence_start(const Alignment_data &aln)
{
    return static_cast<int>(aln.five_p_offset) + 1;
}
int alignment_data_sequence_end(const Alignment_data &aln)
{
    return static_cast<int>(aln.three_p_offset) + 1;
}
int alignment_data_germline_start(const Alignment_data &aln)
{
    return static_cast<int>(aln.five_p_offset) - aln.offset + 1;
}
int alignment_data_germline_end(const Alignment_data &aln)
{
    return static_cast<int>(aln.three_p_offset) - aln.offset + 1;
}

/**
 * Extend alignment mismatches to regions outside the core alignment by comparing sequences directly.
 * This standalone function is useful for importing alignments from external software (igblast, mixcr, etc.)
 * where we only have the core CIGAR alignment data, not the DP matrices.
 * 
 * This function assumes no indels in the extended regions and simply compares sequences directly,
 * skipping positions that are insertions or deletions.
 * 
 * \param int_data_sequence     Query sequence (0-based).
 * \param int_genomic_sequence  Reference sequence (0-based).
 * \param offset               Alignment offset (index on the target sequence where the first nucleotide 
 *                             of the FULL genomic template aligns).
 * \param five_p_offset         5' position (0-based) of the first aligned nucleotide in the target sequence.
 * \param three_p_offset        3' position (0-based) of the last aligned nucleotide in the target sequence.
 * \param insertions            Indices on the TARGET of inserted nucleotides.
 * \param deletions            Indices on the GENOMIC TEMPLATE of deleted nucleotides.
 * \return Vector of mismatch positions (0-based query coordinates) found in the extended regions.
 */
vector<int> extend_alignment_mismatches(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                        int offset, size_t five_p_offset, size_t three_p_offset,
                                        const vector<int> &insertions, const vector<int> &deletions)
{
    vector<int> extended_mismatches;

    // Convert to unordered_sets for efficient lookup
    unordered_set<int> insertion_set(insertions.begin(), insertions.end());
    unordered_set<int> deletion_set(deletions.begin(), deletions.end());

    // 5' extension (left of alignment)
    // Start from beginning of sequences up to five_p_offset
    for (size_t i = 0; i < five_p_offset; ++i) {
        size_t ref_pos = i + offset; // Reference position accounting for offset

        // Skip if this position is an insertion or deletion
        if (insertion_set.count(static_cast<int>(i)) || deletion_set.count(static_cast<int>(ref_pos))) {
            continue;
        }

        // Check if we're within sequence bounds
        if (i < int_data_sequence.size() && ref_pos < int_genomic_sequence.size()) {
            if (!comp_nt_int(int_genomic_sequence[ref_pos], int_data_sequence[i])) {
                extended_mismatches.push_back(static_cast<int>(i));
            }
        }
    }

    // 3' extension (right of alignment)
    // Start from end of alignment to end of sequences
    for (size_t i = three_p_offset + 1; i < int_data_sequence.size(); ++i) {
        size_t ref_pos = i + offset; // Reference position accounting for offset

        // Skip if this position is an insertion or deletion
        if (insertion_set.count(static_cast<int>(i)) || deletion_set.count(static_cast<int>(ref_pos))) {
            continue;
        }

        // Check if we're within sequence bounds
        if (i < int_data_sequence.size() && ref_pos < int_genomic_sequence.size()) {
            if (!comp_nt_int(int_genomic_sequence[ref_pos], int_data_sequence[i])) {
                extended_mismatches.push_back(static_cast<int>(i));
            }
        }
    }

    return extended_mismatches;
}

/**
 * Compare two Alignment_data objects for equality.
 * Two alignments are considered equal if they have:
 * - Same gene name
 * - Same offset
 * - Same five_p_offset and three_p_offset
 * - Same insertions, deletions, mismatches (as sets, order-independent)
 * - Same align_length
 * - Same score (within tolerance)
 *
 * \param a First alignment to compare
 * \param b Second alignment to compare
 * \param score_tolerance Tolerance for score comparison (default: 1e-9)
 * \return true if alignments are considered equal, false otherwise
 */
bool alignment_data_equal(const Alignment_data &a, const Alignment_data &b, double score_tolerance /*= 1e-9*/)
{
    using namespace std;
    // Check basic fields
    if (a.gene_name != b.gene_name) return false;
    if (a.offset != b.offset) return false;
    if (a.five_p_offset != b.five_p_offset) return false;
    if (a.three_p_offset != b.three_p_offset) return false;
    if (a.align_length != b.align_length) return false;
    if (fabs(a.score - b.score) > score_tolerance) return false;
    
    // Check insertions (convert to sets for order-independent comparison)
    unordered_set<int> a_ins(a.get_all_insertions().begin(), a.get_all_insertions().end());
    unordered_set<int> b_ins(b.get_all_insertions().begin(), b.get_all_insertions().end());
    if (a_ins != b_ins) return false;
    
    // Check deletions
    unordered_set<int> a_del(a.get_all_deletions().begin(), a.get_all_deletions().end());
    unordered_set<int> b_del(b.get_all_deletions().begin(), b.get_all_deletions().end());
    if (a_del != b_del) return false;
    
    // Check mismatches (already sorted, but compare as sets to be safe)
    unordered_set<int> a_mis(a.get_all_mismatches().begin(), a.get_all_mismatches().end());
    unordered_set<int> b_mis(b.get_all_mismatches().begin(), b.get_all_mismatches().end());
    if (a_mis != b_mis) return false;
    
    return true;
}

/*
 * This method reads the indexed sequences from a given file(@filename)
 * The structure of the file is assumed to be the same as the one created by the Aligner::write_indexed_seq_csv method
 */
forward_list<pair<const int, const string>> read_indexed_seq_csv(string filename)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    string line_str;
    forward_list<pair<const int, const string>> indexed_seq_list;
    //get rid of the first line
    getline(infile, line_str);
    while (getline(infile, line_str)) {
        size_t scolon_index = line_str.find(';');
        indexed_seq_list.push_front(pair<const int, const string>(stoi(line_str.substr(0, scolon_index)),
                                                                  line_str.erase(0, (scolon_index + 1))));
    }
    return indexed_seq_list;
}

/*
 * This method reads the alignment data from a given file (@filename).
 * The structure of the file is assumed to be the same as the one created by the Aligner::write_alignments_seq_csv method
 */
unordered_map<int, vector<Alignment_data>> read_alignments_seq_csv(const string &filename, double score_threshold,
                                                                   bool allow_in_dels)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename);
    }
    string line_str;
    unordered_map<int, vector<Alignment_data>> indexed_alignments;
    //get rid of the first line
    getline(infile, line_str);
    while (getline(infile, line_str)) {
        auto align = parse_single_alignment_csv_line(line_str);
        if (align.second.score >= score_threshold) {
            if (allow_in_dels || (align.second.get_all_deletions().empty() && align.second.get_all_insertions().empty())) {
                indexed_alignments[align.first].push_back(align.second);
            }
        }
    }
    return indexed_alignments;
}

std::pair<int, Alignment_data> parse_single_alignment_csv_line(const string &line)
{
    size_t index_sep = line.find(';');
    size_t name_sep = line.find(';', index_sep + 1);
    size_t score_sep = line.find(';', name_sep + 1);
    size_t off_sep = line.find(';', score_sep + 1);
    size_t ins_sep = line.find(';', off_sep + 1);
    size_t del_sep = line.find(';', ins_sep + 1);
    size_t mism_sep = line.find(';', del_sep + 1);

    int index = stoi(line.substr(0, index_sep));
    string gene_name = line.substr((index_sep + 1), (name_sep - index_sep - 1));
    double score = stod(line.substr((name_sep + 1), (score_sep - name_sep - 1)));
    int offset = stoi(line.substr((score_sep + 1), (off_sep - score_sep - 1)));
    vector<int> insertions;
    vector<int> deletions;
    vector<int> mismatches;

    string ins_substr = line.substr((off_sep + 2), (ins_sep - off_sep - 3));
    size_t comma_index = ins_substr.find(',');
    if (comma_index != string::npos) {
        insertions.push_back(stoi(ins_substr.substr(0, comma_index)));
        while (comma_index != string::npos) {
            size_t next_comma_index = ins_substr.find(',', (comma_index + 1));
            insertions.push_back(stoi(ins_substr.substr((comma_index + 1), (next_comma_index - comma_index - 1))));
            comma_index = next_comma_index;
        }
    } else if (!ins_substr.empty()) {
        insertions.push_back(stoi(ins_substr));
    }

    string del_substr = line.substr((ins_sep + 2), (del_sep - ins_sep - 3));
    comma_index = del_substr.find(',');
    if (comma_index != string::npos) {
        deletions.push_back(stoi(del_substr.substr(0, comma_index)));
        while (comma_index != string::npos) {
            size_t next_comma_index = del_substr.find(',', (comma_index + 1));
            deletions.push_back(stoi(del_substr.substr((comma_index + 1), (next_comma_index - comma_index - 1))));
            comma_index = next_comma_index;
        }
    } else if (!del_substr.empty()) {
        deletions.push_back(stoi(del_substr));
    }

    string mismatch_substr;
    if (mism_sep == string::npos) {
        mismatch_substr = line.substr((del_sep + 2), (line.size() - del_sep - 3));
    } else {
        mismatch_substr = line.substr((del_sep + 2), (mism_sep - del_sep - 3));
    }

    comma_index = mismatch_substr.find(',');
    if (comma_index != string::npos) {
        mismatches.push_back(stoi(mismatch_substr.substr(0, comma_index)));
        while (comma_index != string::npos) {
            size_t next_comma_index = mismatch_substr.find(',', (comma_index + 1));
            mismatches.push_back(stoi(mismatch_substr.substr((comma_index + 1), (next_comma_index - comma_index - 1))));
            comma_index = next_comma_index;
        }
    } else if (!mismatch_substr.empty()) {
        mismatches.push_back(stoi(mismatch_substr));
    }

    size_t align_length = 0;
    size_t five_p_offset = 0;
    size_t three_p_offset = 0;

    if (mism_sep != string::npos) {
        size_t len_sep = line.find(';', mism_sep + 1);
        if (len_sep == string::npos) {
            string len_substr = line.substr(mism_sep + 1);
            if (!len_substr.empty()) {
                align_length = static_cast<size_t>(stoul(len_substr));
                five_p_offset = static_cast<size_t>(max(0, offset));
                size_t deletion_count = distance(deletions.begin(), deletions.end());
                three_p_offset = (align_length > deletion_count) ? five_p_offset + align_length - 1 - deletion_count
                                                                 : five_p_offset;
            }
        } else {
            string len_substr = line.substr(mism_sep + 1, len_sep - mism_sep - 1);
            if (!len_substr.empty()) {
                align_length = static_cast<size_t>(stoul(len_substr));
            }
            size_t five_sep = line.find(';', len_sep + 1);
            if (five_sep == string::npos) {
                five_p_offset = static_cast<size_t>(max(0, offset));
                size_t deletion_count = distance(deletions.begin(), deletions.end());
                three_p_offset = (align_length > deletion_count) ? five_p_offset + align_length - 1 - deletion_count
                                                                 : five_p_offset;
            } else {
                string five_substr = line.substr(len_sep + 1, five_sep - len_sep - 1);
                string three_substr = line.substr(five_sep + 1);
                five_p_offset = five_substr.empty() ? static_cast<size_t>(max(0, offset))
                                                    : static_cast<size_t>(stoul(five_substr));
                if (!three_substr.empty()) {
                    three_p_offset = static_cast<size_t>(stoul(three_substr));
                } else {
                    size_t deletion_count = distance(deletions.begin(), deletions.end());
                    three_p_offset = (align_length > deletion_count) ? five_p_offset + align_length - 1 - deletion_count
                                                                     : five_p_offset;
                }
            }
        }
    }

    return { index,
             Alignment_data(gene_name, offset, five_p_offset, three_p_offset, align_length, insertions, deletions,
                            mismatches, score, 0, 0) };
}

unordered_map<int, forward_list<Alignment_data>>
Aligner::read_alignments_seq_csv(string filename, double score_threshold, bool allow_in_dels)
{
    unordered_map<int, vector<Alignment_data>> parsed =
            ::read_alignments_seq_csv(filename, score_threshold, allow_in_dels);
    unordered_map<int, forward_list<Alignment_data>> converted;
    for (const auto &entry : parsed) {
        for (auto it = entry.second.rbegin(); it != entry.second.rend(); ++it) {
            converted[entry.first].push_front(*it);
        }
    }
    return converted;
}

unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>>
read_alignments_seq_csv(const string &filename, Gene_class aligned_gene, double score_threshold, bool allow_in_dels,
                        const vector<pair<const int, const string>> &indexed_sequences)
{
    unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> sorted_alignments;
    sorted_alignments = read_alignments_seq_csv(filename, aligned_gene, score_threshold, allow_in_dels,
                                                indexed_sequences, sorted_alignments);
    return sorted_alignments;
}

unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> read_alignments_seq_csv(
        const string &filename, Gene_class aligned_gene, double score_threshold, bool allow_in_dels,
        const vector<pair<const int, const string>> &indexed_sequences,
        unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> sorted_alignments)
{
    unordered_map<int, vector<Alignment_data>> alignments =
            read_alignments_seq_csv(filename, score_threshold, allow_in_dels);
    for (vector<pair<const int, const string>>::const_iterator seq_it = indexed_sequences.begin();
         seq_it != indexed_sequences.end(); ++seq_it) {
        sorted_alignments[(*seq_it).first].second[aligned_gene] = alignments[(*seq_it).first];
        sorted_alignments[(*seq_it).first].first = (*seq_it).second;
    }
    //sort(alignments.begin() , alignments.end() , align_compare);
    return sorted_alignments;
}

unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>>
read_alignments_seq_csv_score_range(const string &filename, Gene_class aligned_gene, double score_range,
                                    bool allow_in_dels, const vector<pair<const int, const string>> &indexed_sequences)
{
    unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> sorted_alignments;
    sorted_alignments = read_alignments_seq_csv_score_range(filename, aligned_gene, score_range, allow_in_dels,
                                                            indexed_sequences, sorted_alignments);
    return sorted_alignments;
}

unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> read_alignments_seq_csv_score_range(
        const string &filename, Gene_class aligned_gene, double score_range, bool allow_in_dels,
        const vector<pair<const int, const string>> &indexed_sequences,
        unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> sorted_alignments)
{
    unordered_map<int, vector<Alignment_data>> alignments = read_alignments_seq_csv(filename, 0, allow_in_dels);
    for (vector<pair<const int, const string>>::const_iterator seq_it = indexed_sequences.begin();
         seq_it != indexed_sequences.end(); ++seq_it) {
        vector<Alignment_data> &seq_alignments = alignments[(*seq_it).first];
        double max_score = -1;
        for (vector<Alignment_data>::const_iterator align_it = seq_alignments.begin(); align_it != seq_alignments.end();
             ++align_it) {
            if ((*align_it).score > max_score) {
                max_score = (*align_it).score;
            }
        }
        // use reverse iterator
        for (auto align_it = seq_alignments.rbegin(); align_it != seq_alignments.rend();) {
            if ((*align_it).score < (max_score - score_range)) {
                // Convert return value of erase back to reverse_iterator
                align_it = decltype(align_it)(seq_alignments.erase(std::next(align_it).base()));
            } else {
                ++align_it;
            }
        }
        sort(seq_alignments.begin(), seq_alignments.end(), align_compare);
        sorted_alignments[(*seq_it).first].second[aligned_gene] = alignments[(*seq_it).first];
        sorted_alignments[(*seq_it).first].first = (*seq_it).second;
    }
    return sorted_alignments;
}

vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>>
map2vect(unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>> alignments_map)
{
    vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>> alignmets_vect;
    for (unordered_map<int, pair<string, unordered_map<Gene_class, vector<Alignment_data>>>>::const_iterator seq_it =
                 alignments_map.begin();
         seq_it != alignments_map.end(); ++seq_it) {
        alignmets_vect.emplace_back((*seq_it).first, (*seq_it).second.first, (*seq_it).second.second);
    }
    return alignmets_vect;
}

void Aligner::set_genomic_sequences(vector<pair<string, string>> nt_genomic_seq)
{
    this->nt_genomic_sequences = forward_list<pair<string, string>>();
    this->int_genomic_sequences = forward_list<pair<string, Int_Str>>();
    for (vector<pair<string, string>>::const_iterator iter = nt_genomic_seq.begin(); iter != nt_genomic_seq.end();
         ++iter) {
        nt_genomic_sequences.emplace_front((*iter).first, (*iter).second);

        int_genomic_sequences.emplace_front((*iter).first, nt2int((*iter).second));
    }
}

/*
 * This method will incorporate gaps('-') at the places where insertion or deletion occured both in the data and genomic sequence
 * A deletion will correspond to a gap introduced in the data sequence
 * An insertion to a gap in the genomic sequence
 * <\code>prev_dels<code> will indicate the shift in the data_seq offset induced by the introduction of previous deletions
 *
 * The method returns the shift induced by introducing the deletions of this alignment
 */
int Aligner::incorporate_in_dels(string &data_seq, string &genomic_seq, const vector<int>,
                                 const vector<int>, int prev_dels)
{

    return prev_dels;
}

/*
 * Convert nucleotide alphabet sequence into int sequence
 * Conventions are the same as nt2int matlab function
 */
Int_Str nt2int(const string &nt_sequence)
{
    Int_Str int_seq;
    for (size_t i = 0; i != nt_sequence.size(); ++i) {
        if (nt_sequence[i] == 'A') {
            int_seq.append(int_A);
        } else if (nt_sequence[i] == 'C') {
            int_seq.append(int_C);
        } else if (nt_sequence[i] == 'G') {
            int_seq.append(int_G);
        } else if ((nt_sequence[i] == 'T') or (nt_sequence[i] == 'U')) {
            int_seq.append(int_T);
        } else if (nt_sequence[i] == 'R') {
            int_seq.append(int_R);
        } else if (nt_sequence[i] == 'Y') {
            int_seq.append(int_Y);
        } else if (nt_sequence[i] == 'K') {
            int_seq.append(int_K);
        } else if (nt_sequence[i] == 'M') {
            int_seq.append(int_M);
        } else if (nt_sequence[i] == 'S') {
            int_seq.append(int_S);
        } else if (nt_sequence[i] == 'W') {
            int_seq.append(int_W);
        } else if (nt_sequence[i] == 'B') {
            int_seq.append(int_B);
        } else if (nt_sequence[i] == 'D') {
            int_seq.append(int_D);
        } else if (nt_sequence[i] == 'H') {
            int_seq.append(int_H);
        } else if (nt_sequence[i] == 'V') {
            int_seq.append(int_V);
        } else if (nt_sequence[i] == 'N') {
            int_seq.append(int_N);
        } else {
            cerr << "print:" << nt_sequence << endl;
            cerr << i << endl;
            throw runtime_error("Unknown nucleotide: " + to_string(nt_sequence[i]) + "in string " + nt_sequence
                                + "in Aligner::nt2int");
        }
    }
    return int_seq;
}

/**
 * This function compares nucleotides and output a boolean if they do not necessarily imply an error (ambiguous nucleotides are thus treated in a loose sense).
 */
bool comp_nt_int(const int &nt_1, const int &nt_2)
{
    if (nt_1 != nt_2) {
        if ((nt_1 < 4) && (nt_2 < 4)) {
            return false;
        } else {
            switch (nt_1) {
            case int_A:
                switch (nt_2) {
                case int_R:
                case int_W:
                case int_M:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_C:
                switch (nt_2) {
                case int_Y:
                case int_S:
                case int_M:
                case int_B:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_G:
                switch (nt_2) {
                case int_R:
                case int_S:
                case int_K:
                case int_B:
                case int_D:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_T:
                switch (nt_2) {
                case int_Y:
                case int_W:
                case int_K:
                case int_B:
                case int_D:
                case int_H:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_R:
                switch (nt_2) {
                case int_A:
                case int_G:
                case int_S:
                case int_W:
                case int_K:
                case int_M:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_Y:
                switch (nt_2) {
                case int_C:
                case int_T:
                case int_S:
                case int_W:
                case int_K:
                case int_M:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_K:
                switch (nt_2) {
                case int_G:
                case int_T:
                case int_R:
                case int_Y:
                case int_S:
                case int_W:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_M:
                switch (nt_2) {
                case int_A:
                case int_C:
                case int_R:
                case int_Y:
                case int_S:
                case int_W:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_S:
                switch (nt_2) {
                case int_G:
                case int_C:
                case int_R:
                case int_Y:
                case int_K:
                case int_M:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_W:
                switch (nt_2) {
                case int_A:
                case int_T:
                case int_R:
                case int_Y:
                case int_K:
                case int_M:
                case int_B:
                case int_D:
                case int_H:
                case int_V:
                case int_N:
                    return true;
                    break;
                }
                break;
            case int_B:
                if (nt_2 != int_A) { //Of course this is in the hope that nt_2 is in the correct range of int
                    return true;
                }
                break;
            case int_D:
                if (nt_2 != int_C) { //Of course this is in the hope that nt_2 is in the correct range of int
                    return true;
                }
                break;
            case int_H:
                if (nt_2 != int_G) { //Of course this is in the hope that nt_2 is in the correct range of int
                    return true;
                }
                break;
            case int_V:
                if (nt_2 != int_T) { //Of course this is in the hope that nt_2 is in the correct range of int
                    return true;
                }
                break;
            case int_N:
                return true;
                break;
            default:
                throw runtime_error("Unknown nucleotide index: " + to_string(nt_1) + "in comp_nt_int()");
            }
            return false;
        }
    } else {
        return true;
    }
}

list<Int_nt> get_ambiguous_nt_list(const Int_nt &ambiguous_nt)
{
    list<Int_nt> nt_list;

    if ((ambiguous_nt == int_A) or (ambiguous_nt == int_C) or (ambiguous_nt == int_G) or (ambiguous_nt == int_T)) {
        nt_list.emplace_back(ambiguous_nt);
    } else {
        bool any_true = false;
        //Add an A for all cases implying an A	(do not add a break to allow for execution of other cases and possibly add more letters to the list)
        if ((ambiguous_nt == int_R) or (ambiguous_nt == int_W) or (ambiguous_nt == int_M) or (ambiguous_nt == int_D)
            or (ambiguous_nt == int_H) or (ambiguous_nt == int_V) or (ambiguous_nt == int_N)) {
            nt_list.emplace_back(int_A);
            any_true = true;
        }
        //Same for C
        if ((ambiguous_nt == int_Y) or (ambiguous_nt == int_S) or (ambiguous_nt == int_M) or (ambiguous_nt == int_B)
            or (ambiguous_nt == int_H) or (ambiguous_nt == int_V) or (ambiguous_nt == int_N)) {
            nt_list.emplace_back(int_C);
            any_true = true;
        }
        //Same for G
        if ((ambiguous_nt == int_R) or (ambiguous_nt == int_S) or (ambiguous_nt == int_K) or (ambiguous_nt == int_B)
            or (ambiguous_nt == int_D) or (ambiguous_nt == int_V) or (ambiguous_nt == int_N)) {
            nt_list.emplace_back(int_G);
            any_true = true;
        }
        //Same for T
        if ((ambiguous_nt == int_Y) or (ambiguous_nt == int_W) or (ambiguous_nt == int_K) or (ambiguous_nt == int_B)
            or (ambiguous_nt == int_D) or (ambiguous_nt == int_H) or (ambiguous_nt == int_N)) {
            nt_list.emplace_back(int_T);
            any_true = true;
        }

        if (not any_true) {
            throw runtime_error("Unknown nucleotide index: " + to_string(ambiguous_nt) + "in get_ambiguous_nt_list()");
        }
    }
    return nt_list;
}

/*
 * Randomly sample N indexed sequences from a given vector of indexed sequences
 */
vector<pair<const int, const string>> sample_indexed_seq(const vector<pair<const int, const string>> &indexed_seqs,
                                                         const size_t sample_size)
{

    //Return an error if trying to sample more than the number of available sequences
    if (sample_size > indexed_seqs.size()) {
        throw std::runtime_error("Trying to sample " + to_string(sample_size) + " sequences in a pool of "
                                 + to_string(indexed_seqs.size()) + " sequences in sample_indexed_seq()");
    }

    //Create seed for random generator
    //create a seed from timer
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point time = myclock::now();
    myclock::duration dur = (myclock::time_point::max)() - time;

    //Get a random seed
    uint64_t random_seed = draw_random_64bits_seed();
    //Instantiate random number generator
    mt19937_64 generator = mt19937_64(random_seed);

    //Need to make a copy because of the constness in indexed_seqs
    vector<pair<int, string>> indexed_seqs_copy(indexed_seqs.begin(), indexed_seqs.end());

    shuffle(indexed_seqs_copy.begin(), indexed_seqs_copy.end(), generator);
    return vector<pair<const int, const string>>(indexed_seqs_copy.begin(), indexed_seqs_copy.begin() + sample_size);
}

namespace swalign {

/*
 * SwCandidate, SwDPState, and SwPreparedInputs are defined in AlignerInternal.h
 * (moved there so whitebox tests can construct them directly).
 *
 * Internal workspace for one Smith-Waterman DP execution.
 *
 * Groups the four DP matrices and the candidate-tracking vector that are
 * allocated, mutated, and read together during a single sw_align call. Storing
 * them here avoids threading five separate references through every helper in
 * the SW pipeline.
 *
 * Coordinates: all matrices use the +1 padded convention (row 0 / col 0 are
 * initialization boundaries; sequence positions are 1-based inside the matrix).
 * candidates uses the same 1-based matrix row/column coordinates.
 *
 * candidates is a single vector<SwCandidate> rather than three parallel
 * vector<int> so the per-cell "update the best score for this candidate" access
 * (score/row/col together) touches one cache line instead of up to three.
 */

/**
   * \brief Coordinate conversion parameters for flipped sequences
   * 
   * Bundles the parameters needed to convert coordinates between original and flipped sequences.
   * This struct is used to isolate the coordinate conversion logic for easier testing and debugging.
   */
// Coordinate conversion functions for Smith-Waterman traceback
// These functions handle conversion from DP matrix coordinates to sequence coordinates,
// accounting for sequence reversal when flip_seqs = true.

/**
 * \brief Convert 1-based matrix row coordinate to 0-based query sequence coordinate.
 * 
 * Handles coordinate conversion from DP matrix row index to query sequence position.
 * When sequences are not flipped (normal case): row i (1-based) corresponds to query position i-1 (0-based).
 * When sequences are flipped: row i (1-based) corresponds to query position (data_seq_size - i).
 * 
 * \param i              1-based row coordinate from DP matrix
 * \param data_seq_size   Size of the original (unflipped) query sequence
 * \param flip_seqs      Whether sequences were flipped for alignment
 * 
 * \return 0-based query sequence coordinate
 */
size_t convert_matrix_row_to_query_pos(size_t i, size_t data_seq_size, bool flip_seqs)
{
    if (flip_seqs) {
        // When sequences are flipped, row i (1-based) maps to position (data_seq_size - i)
        return data_seq_size - static_cast<size_t>(i);
    } else {
        // Normal case: row i (1-based) maps to position i-1 (0-based)
        return static_cast<size_t>(i - 1);
    }
}

/**
 * \brief Convert 1-based matrix column coordinate to 0-based reference sequence coordinate.
 * 
 * Handles coordinate conversion from DP matrix column index to reference sequence position.
 * When sequences are not flipped (normal case): column j (1-based) corresponds to reference position j-1 (0-based).
 * When sequences are flipped: column j (1-based) corresponds to reference position (genomic_seq_size - j).
 * 
 * \param j                1-based column coordinate from DP matrix
 * \param genomic_seq_size Size of the original (unflipped) reference sequence
 * \param flip_seqs        Whether sequences were flipped for alignment
 * 
 * \return 0-based reference sequence coordinate
 */
size_t convert_matrix_col_to_ref_pos(size_t j, size_t genomic_seq_size, bool flip_seqs)
{
    if (flip_seqs) {
        // When sequences are flipped, column j (1-based) maps to position (genomic_seq_size - j)
        return genomic_seq_size - static_cast<size_t>(j);
    } else {
        // Normal case: column j (1-based) maps to position j-1 (0-based)
        return static_cast<size_t>(j - 1);
    }
}

struct SwReconstructionResult
{
    list<pair<int, Alignment_data>> alignments;
    double max_align_score;
};

/**
 * Prepare Smith-Waterman inputs before DP matrix allocation.
 *
 * Coordinates: sequences remain in 0-based nucleotide indexing, while the DP matrix
 * will use +1 row/column padding on top of these prepared strings.
 * Mutation: returns copied sequences so callers can safely reverse in place.
 *
 * \param config  Run policy; only config.flip_seqs is consulted here.
 */
SwPreparedInputs prepare_sw_inputs(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                   const SwDPConfig &config)
{
    SwPreparedInputs prepared{ int_data_sequence, int_genomic_sequence, 0 };
    if (config.alignment_mode.reverse_sequences) {
        reverse(prepared.data_sequence.begin(), prepared.data_sequence.end());
        reverse(prepared.genomic_sequence.begin(), prepared.genomic_sequence.end());
    }
    return prepared;
}

/**
 * Initialize score and traceback support matrices for a fresh DP run.
 *
 * Coordinates: the DP matrices use +1 padded dimensions where row 0 and column 0 are
 * initialization boundaries. Matrix values are written using this padded convention.
 * Mutation: fully initializes all matrix fields of dp; candidate vectors are left empty.
 *
 * \param dp       The DP workspace to initialize. n_rows and n_cols must already be set.
 * \param local_align  True for vanilla local (SW) alignment; false for semi-global.
 * \param gap_penalty  Linear gap penalty applied to column-0 initialization when semi-global.
 */
void initialize_sw_matrices(SwDPState &dp, const SwDPConfig &config)
{
    for (int i = 0; i != dp.n_rows; ++i) {
        if (config.alignment_mode.data_leading_free) {
            // free leading deletion in query
            // vanilla SW local alignment
            dp.score_matrix(i, 0) = 0;
        } else {
            // penalized leading deletion in query
            // akin to global alignment on the left/5'
            dp.score_matrix(i, 0) = -i * config.gap_penalty;
        }
        dp.col_memory_matrix(i, 0) = 0;
        dp.row_memory_matrix(i, 0) = 0;
        // Only the boundary needs -1: every interior cell (i,j >= 1) is unconditionally
        // overwritten by the fill routine before ever being read.
        dp.alignment_numb_tracker(i, 0) = -1;
    }

    for (int j = 0; j != dp.n_cols; ++j) {
        if (config.alignment_mode.genomic_leading_free) {
            // free leading insertion in query
            dp.score_matrix(0, j) = 0;
        } else {
            dp.score_matrix(0, j) = -j * config.gap_penalty;
        }
        dp.col_memory_matrix(0, j) = 0;
        dp.row_memory_matrix(0, j) = 0;
        dp.alignment_numb_tracker(0, j) = -1;
    }
}

/**
 * Extend alignment mismatches to the 5' (left) side assuming no indels in the extended region.
 * This function extends from the start of the core alignment (i, j) towards lower indices
 * by making diagonal moves only, checking for mismatches along the way.
 * Only applicable for local alignments where data_leading_free is true.
 *
 * \param prepared           Prepared (possibly flipped) copy sequences and offset_change.
 * \param i_start            Matrix row coordinate (1-based) of the start of the core alignment.
 * \param j_start            Matrix column coordinate (1-based) of the start of the core alignment.
 * \param data_seq_size      Size of the original (unflipped) query sequence.
 * \param genomic_seq_size   Size of the original (unflipped) reference sequence.
 * \param flip_seqs          Whether sequences were flipped for alignment.
 * \param matrix_n_rows      Total number of rows in the DP matrix.
 * \param matrix_n_cols      Total number of columns in the DP matrix.
 * \return Vector of mismatch positions (0-based query coordinates) found in the 5' extended region.
 */
vector<int> ungapped_extend_align_5p_from_dp(const SwPreparedInputs &prepared, int i_start, int j_start,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols)
{
    vector<int> extended_mismatches;

    // Extend 5' (towards lower indices) - always diagonal, no indels assumed
    // Start from the position just before the alignment start
    int i = i_start - 1;
    int j = j_start - 1;

    // Continue while we're within matrix boundaries (greater than 1 because matrices are +1 padded)
    while (i >= 1 && j >= 1) {
        if (!comp_nt_int(prepared.data_sequence.at(convert_matrix_row_to_query_pos(i, data_seq_size, false)),
                         prepared.genomic_sequence.at(convert_matrix_col_to_ref_pos(j, genomic_seq_size, false)))) {
            extended_mismatches.push_back(convert_matrix_row_to_query_pos(i, data_seq_size, flip_seqs));
        }
        --i;
        --j;
    }

    // Reverse to maintain order from 5' to alignment start (increasing query coordinates)
    reverse(extended_mismatches.begin(), extended_mismatches.end());
    return extended_mismatches;
}

/**
 * Extend alignment mismatches to the 3' (right) side assuming no indels in the extended region.
 * This function extends from the end of the core alignment (i_end, j_end) towards higher indices
 * by making diagonal moves only, checking for mismatches along the way.
 * Applicable for both local and semi-global alignments.
 *
 * \param prepared           Prepared (possibly flipped) copy sequences and offset_change.
 * \param i_end              Matrix row coordinate (1-based) of the end of the core alignment.
 * \param j_end              Matrix column coordinate (1-based) of the end of the core alignment.
 * \param data_seq_size      Size of the original (unflipped) query sequence.
 * \param genomic_seq_size   Size of the original (unflipped) reference sequence.
 * \param flip_seqs          Whether sequences were flipped for alignment.
 * \param matrix_n_rows      Total number of rows in the DP matrix.
 * \param matrix_n_cols      Total number of columns in the DP matrix.
 * \return Vector of mismatch positions (0-based query coordinates) found in the 3' extended region.
 */
vector<int> ungapped_extend_align_3p_from_dp(const SwPreparedInputs &prepared, int i_end, int j_end,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols)
{
    vector<int> extended_mismatches;

    // Extend 3' (towards higher indices) - always diagonal, no indels assumed
    // Start from the position just after the alignment end
    int i = i_end + 1;
    int j = j_end + 1;

    // Continue while we're within matrix boundaries
    while (i < matrix_n_rows && j < matrix_n_cols) {
        if (!comp_nt_int(prepared.data_sequence.at(convert_matrix_row_to_query_pos(i, data_seq_size, false)),
                         prepared.genomic_sequence.at(convert_matrix_col_to_ref_pos(j, genomic_seq_size, false)))) {
            extended_mismatches.push_back(convert_matrix_row_to_query_pos(i, data_seq_size, flip_seqs));
        }
        ++i;
        ++j;
    }

    return extended_mismatches;
}

/**
 * Merge and sort mismatch vectors from core alignment and extended regions.
 *
 * \param core_mismatches       Mismatches from the core alignment.
 * \param extended_mismatches   Mismatches from extended regions (5' and/or 3').
 * \return Sorted vector containing all mismatches (core + extended).
 */
vector<int> merge_and_sort_mismatches(const vector<int> &core_mismatches, const vector<int> &extended_mismatches)
{
    vector<int> all_mismatches = core_mismatches;
    all_mismatches.insert(all_mismatches.end(), extended_mismatches.begin(), extended_mismatches.end());
    sort(all_mismatches.begin(), all_mismatches.end());
    return all_mismatches;
}

/**
 * Merge and sort mismatch vectors from multiple sources (core + 5' extension + 3' extension).
 *
 * \param core_mismatches       Mismatches from the core alignment.
 * \param extended_5p_mismatches Mismatches from 5' extended region.
 * \param extended_3p_mismatches Mismatches from 3' extended region.
 * \return Sorted vector containing all mismatches (core + 5' extended + 3' extended).
 */
vector<int> merge_and_sort_mismatches(const vector<int> &core_mismatches, const vector<int> &extended_5p_mismatches,
                                      const vector<int> &extended_3p_mismatches)
{
    vector<int> all_mismatches = core_mismatches;
    all_mismatches.insert(all_mismatches.end(), extended_5p_mismatches.begin(), extended_5p_mismatches.end());
    all_mismatches.insert(all_mismatches.end(), extended_3p_mismatches.begin(), extended_3p_mismatches.end());
    sort(all_mismatches.begin(), all_mismatches.end());
    return all_mismatches;
}

/**
 * Trace back candidate alignments from max-score endpoints and build Alignment_data objects.
 *
 * Coordinates: dp matrix coordinates are 1-based padded indices; produced insertion/deletion/
 * mismatch coordinates keep the 0-based conventions currently used by Alignment_data.
 * Mutation: reads dp matrices and appends to the result list.
 *
 * \param int_data_sequence     Original (possibly un-flipped) query sequence, 0-based.
 * \param int_genomic_sequence  Original (possibly un-flipped) reference sequence, 0-based.
 * \param prepared              Prepared (possibly flipped) copy sequences and offset_change.
 * \param dp                    Completed DP workspace (read-only matrices and candidate vectors).
 */
SwReconstructionResult traceback_sw_alignments(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                               const SwPreparedInputs &prepared, const SwDPState &dp,
                                               const SwDPConfig &config)
{
    double score_threshold = config.score_threshold;
    const double best_score = std::max_element(dp.candidates.begin(), dp.candidates.end(),
                                               [](const SwCandidate &a, const SwCandidate &b) {
                                                   return a.score < b.score;
                                               })
                                       ->score;
    if (config.best_only && best_score >= config.score_threshold) {
        score_threshold = best_score;
    }
    const int min_offset = config.min_offset;
    const int max_offset = config.max_offset;
    const bool flip_seqs = config.alignment_mode.reverse_sequences;
    SwReconstructionResult output;
    output.max_align_score = 0;

    // Get sequence sizes for coordinate conversion
    const size_t data_seq_size = int_data_sequence.size();
    const size_t genomic_seq_size = int_genomic_sequence.size();

    for (size_t align = 0; align != dp.candidates.size(); ++align) {
        if (dp.candidates[align].score >= score_threshold) {

            vector<int> mismatches;
            vector<int> insertions;
            vector<int> deletions;
            size_t align_length = 0;

            bool end_of_alignment = false;

            size_t i = dp.candidates[align].row;
            size_t j = dp.candidates[align].col;
            // Save the original starting position for end offset calculation
            size_t i_end = i;
            size_t j_end = j;

            // TODO correct this to get the alignment until the end (not just until the best scoring nucl)
            while (!end_of_alignment) {
                if (dp.row_memory_matrix(i, j) == 0) {
                    // Deletion: use column coordinate (j) to get reference position
                    deletions.emplace_back(convert_matrix_col_to_ref_pos(j, genomic_seq_size, flip_seqs));
                } else if (dp.col_memory_matrix(i, j) == 0) {
                    // Insertion: use row coordinate (i) to get query position
                    insertions.emplace_back(convert_matrix_row_to_query_pos(i, data_seq_size, flip_seqs));
                } else {
                    if (!comp_nt_int(
                                prepared.data_sequence.at(convert_matrix_row_to_query_pos(i, data_seq_size, false)),
                                prepared.genomic_sequence.at(
                                        convert_matrix_col_to_ref_pos(j, genomic_seq_size, false)))) {
                        mismatches.emplace_back(convert_matrix_row_to_query_pos(i, data_seq_size, flip_seqs));
                    }
                }
                ++align_length;
                int i_temp = i;
                int j_temp = j;
                i -= dp.row_memory_matrix(i_temp, j);
                j -= dp.col_memory_matrix(i_temp, j);
                if ((dp.row_memory_matrix(i, j) == 0) && (dp.col_memory_matrix(i, j) == 0)) {
                    end_of_alignment = true;
                    // undo last move to remain away from initialization values
                    i += dp.row_memory_matrix(i_temp, j_temp);
                    j += dp.col_memory_matrix(i_temp, j_temp);
                    break;
                }
            }

            // Convert alignment boundaries from matrix coordinates to sequence coordinates
            // i, j are now at the start of the alignment (5' end after traceback)
            //int offset = convert_matrix_coords_to_offset(i, j, data_seq_size, genomic_seq_size, prepared.offset_change, flip_seqs);
            size_t begin_align_offset = convert_matrix_row_to_query_pos(i, data_seq_size, flip_seqs);
            size_t end_align_offset = convert_matrix_row_to_query_pos(i_end, data_seq_size, flip_seqs);
            int offset;
            if (flip_seqs) {
                // reverse offset order
                std::swap(begin_align_offset, end_align_offset);
                // assume that leading deletions (reverse trailing, hence j_end) would align 1 to 1 with the read.
                offset = static_cast<int>(begin_align_offset) - convert_matrix_col_to_ref_pos(j_end, genomic_seq_size, true);
            } else {
                /*
             * FIXME: this does not really make sense for local alignments. 
              It boils down to assuming that leading deletions would align 1 to 1 with the read.
              But it is what is expected by the legacy alignment data representation.
             * */
                offset = static_cast<int>(begin_align_offset) - convert_matrix_col_to_ref_pos(j, genomic_seq_size, flip_seqs);
                // reverse containers that have been filled via push back
                reverse(mismatches.begin(), mismatches.end());
                reverse(insertions.begin(), insertions.end());
                reverse(deletions.begin(), deletions.end());
            }

            if ((offset >= min_offset) && (offset <= max_offset)) {
                // TODO reduce computation time by truncating alignment from the beginning? = banded alignment
                // TODO change this and use incorporate_in_dels(), should probably change the list inside alignment data also to have the actual corresponding indices
                // TODO return the actual inserted/deleted sequences in the alignment data??

                // Extend alignment to capture mismatches in deleted V/J nucleotides
                // Only perform extension if enabled in config
                vector<int> extended_mismatches_5p;
                vector<int> extended_mismatches_3p;

                if (config.enable_extension) {
                    if(config.alignment_mode.is_local_alignment()){
                        // 5' extension: extend from the start of the alignment
                        extended_mismatches_5p = ungapped_extend_align_5p_from_dp(
                                prepared, i, j, data_seq_size, genomic_seq_size, flip_seqs, dp.n_rows, dp.n_cols);
                    }
                    
                    // 3' extension: extend from the end of the alignment  
                    extended_mismatches_3p = ungapped_extend_align_3p_from_dp(
                            prepared, i_end, j_end, data_seq_size, genomic_seq_size, flip_seqs, dp.n_rows, dp.n_cols);
                }

                // Merge core mismatches with extended mismatches and sort
                // TODO avoid sorting ops with proper design
                vector<int> all_mismatches =
                        merge_and_sort_mismatches(mismatches, extended_mismatches_5p, extended_mismatches_3p);

                if (dp.candidates[align].score > output.max_align_score) {
                    output.max_align_score = dp.candidates[align].score;
                }
                output.alignments.emplace_back(pair<int, Alignment_data>(
                        dp.candidates[align].score,
                        Alignment_data(offset, begin_align_offset, end_align_offset, align_length,
                                       insertions, deletions, all_mismatches,
                                       dp.candidates[align].score, int_data_sequence.size(),
                                       int_genomic_sequence.size())));
            }
        }
    }

    return output;
}

/**
 * Score a single DP cell and record the chosen predecessor direction.
 *
 * Computes the substitution, query-gap, and reference-gap candidate scores for
 * cell (i, j) using the Aligner's substitution matrix and gap penalty, writes
 * the winning score and predecessor flags into the raw DP matrix storage, and
 * updates candidates if a new high score is reached along the current
 * alignment path.
 *
 * Coordinates: i is the 1-based row index (query position i-1, 0-based);
 *              j is the 1-based column index (reference position j-1, 0-based).
 * score/row_mem/col_mem/numb_trk are the raw backing storage of dp.score_matrix,
 * dp.row_memory_matrix, dp.col_memory_matrix and dp.alignment_numb_tracker.
 * They share one linear index (i + n_rows*j, matching Matrix::operator()'s
 * formula) because all four matrices in SwDPState are constructed with
 * identical dimensions -- computing that index once per call here and reusing
 * it for every access, instead of each caller-side Matrix::operator() call
 * recomputing it independently, is why callers pass raw pointers/n_rows
 * rather than a SwDPState&. int_data_sequence/int_genomic_sequence are read
 * via operator[] instead of the bounds-checked at() since callers already
 * guarantee valid indices. Profiling identified the per-cell function-call
 * overhead, this repeated index arithmetic, and the bounds checking as the
 * dominant cost of sw_align() -- hence always_inline: this is called once per
 * DP cell from fill_sw_score_matrix's hot loop, and any un-inlined call would
 * reintroduce the per-cell prologue/epilogue and stack-protector overhead this
 * was written to eliminate.
 *
 * Move selection: the three candidate moves -- diagonal (match/mismatch),
 * "up" (query gap / deletion), "left" (reference gap / insertion) -- are
 * indexed 0/1/2 throughout, and their score/predecessor-cell/output-flag data
 * are laid out as parallel 3-element tables read by a single winner index
 * (max_idx) rather than as a nested if/else-if cascade choosing among three
 * differently-sized blocks of side effects. This is a rewrite of an
 * equivalent cascade, not a behavior change: it exists because the nested
 * form's side effects (writes, a conditional push_back) block the compiler
 * from generating branchless code for the part that is genuinely
 * data-dependent and hard to predict -- which of the three scores is
 * largest, which varies unpredictably cell to cell with sequence content.
 * Separating that decision (pure arithmetic, no side effects) from the
 * action (one shared block of writes parameterized by max_idx) lets the
 * decision compile to a couple of conditional moves instead of
 * mispredicted branches; profiling confirmed branch mispredictions as a
 * significant remaining cost once the surrounding per-cell overhead above
 * was removed.
 *
 * Why the nested cascade and the argmax-then-gate form are equivalent: in
 * local-alignment mode (reset_negative_scores) each original branch's guard
 * is "this move is (tied-for-)largest AND its own value is > 0". If the
 * true largest of the three values is <= 0, every other value is <= it and
 * therefore also <= 0 -- so no matter which value the argmax turns out to
 * be, checking that one value's sign reproduces the same pass/fail outcome
 * the original per-branch checks would have produced. Hence: compute the
 * argmax once (with the tie-break priority below), then gate the whole
 * decision on that single value's sign instead of repeating the check
 * inside every branch.
 *
 * Tie-break priority (matches the original's >= comparisons exactly -- see
 * the FIXME below, which this rewrite deliberately does not touch):
 * diagonal beats both gap moves on a tie; "up" beats "left" on a tie.
 * max_idx is computed with strict > comparisons in priority order, so a tie
 * leaves max_idx at the earlier (higher-priority) index.
 *
 * Candidate-tracking asymmetry: only the diagonal move may seed a brand-new
 * tracked candidate, and only when its predecessor cell is itself untracked
 * (numb_trk == -1); the two gap moves always just propagate whatever
 * tracker id their neighbor already has, even if that is -1. CAN_START_NEW
 * encodes this per move index instead of hardcoding it to index 0, so the
 * tracker-update code below is one block shared by all three move types.
 *
 * Mutation: writes score/row_mem/col_mem/numb_trk[idx], and may append to or
 *           update candidates.
 *
 * FIXME: the substitution score takes precedence over equal-scoring gap moves
 * (>= comparison), which collapses branching/convergent traceback paths into a
 * single ancestor. This will be fixed in Step 3 of the refactoring plan.
 */
IGOR_ALWAYS_INLINE void fill_sw_matrix_cell(int i, int j, int n_rows, double *score, int *row_mem, int *col_mem,
                                             int *numb_trk, vector<SwCandidate> &candidates,
                                             bool reset_negative_scores, const SwDPConfig &config,
                                             const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence)
{
    const int idx = i + n_rows * j;
    const int idx_up = idx - 1; // (i-1, j)
    const int idx_left = idx - n_rows; // (i, j-1)
    const int idx_diag = idx_left - 1; // (i-1, j-1)

    const int genomic_gap_score = static_cast<int>(score[idx_left] - config.gap_penalty);
    const int data_gap_score = static_cast<int>(score[idx_up] - config.gap_penalty);
    const int subs_score = static_cast<int>(
            score[idx_diag] + config.substitution_matrix(int_data_sequence[i - 1], int_genomic_sequence[j - 1]));

    // Move index convention used by every table below: 0 = diagonal, 1 = up, 2 = left.
    const int move_scores[3] = { subs_score, data_gap_score, genomic_gap_score };
    const int predecessor_idx[3] = { idx_diag, idx_up, idx_left };
    static constexpr int ROW_MOVE[3] = { 1, 1, 0 };
    static constexpr int COL_MOVE[3] = { 1, 0, 1 };
    static constexpr bool CAN_START_NEW[3] = { true, false, false };

    // argmax over move_scores with the tie-break priority documented above: strict `>`
    // means a tie leaves max_idx at the earlier, higher-priority index.
    // This manual version is much faster than using std::max_element
    int max_idx = 0;
    if (move_scores[1] > move_scores[max_idx]) {
        max_idx = 1;
    }
    if (move_scores[2] > move_scores[max_idx]) {
        max_idx = 2;
    }

    if ((!reset_negative_scores) || (move_scores[max_idx] > 0)) {
        // Keep the winning move (see the equivalence note above for why checking only
        // this one value's sign reproduces the original per-branch positivity checks).
        score[idx] = move_scores[max_idx];
        row_mem[idx] = ROW_MOVE[max_idx];
        col_mem[idx] = COL_MOVE[max_idx];

        const int predecessor_tracker = numb_trk[predecessor_idx[max_idx]];
        if ((predecessor_tracker == -1) && CAN_START_NEW[max_idx]) {
            // Only reachable for the diagonal move (see the candidate-tracking note above).
            numb_trk[idx] = static_cast<int>(candidates.size());
            candidates.push_back(SwCandidate{ move_scores[max_idx], i, j });
        } else {
            numb_trk[idx] = predecessor_tracker;
        }
    } else {
        score[idx] = 0;
        row_mem[idx] = 0;
        col_mem[idx] = 0;
        // TODO check this (local alignment reset)
        numb_trk[idx] = numb_trk[idx_diag];
    }

    // Keep max score in memory
    const int id = numb_trk[idx];
    if (id != -1) {
        SwCandidate &candidate = candidates[id];
        if (score[idx] > candidate.score) {
            candidate.score = static_cast<int>(score[idx]);
            candidate.row = i;
            candidate.col = j;
        }
    }
}

/**
 * Fill the Smith-Waterman score matrix and alignment trackers.
 *
 * Coordinates: this routine fills the +1 padded DP matrix starting at cell (1,1).
 * Mutation: updates dp.score_matrix, dp.row/col_memory_matrix, dp.alignment_numb_tracker,
 * and dp.candidates in place.
 *
 * is_local_alignment() and the four matrices' raw storage pointers are
 * computed once here (rather than once per cell) and passed to
 * fill_sw_matrix_cell, which is marked always_inline -- see its doc comment.
 *
 * Traversal order: cell (i,j) only depends on (i-1,j), (i,j-1) and (i-1,j-1),
 * so any traversal that visits both coordinates in increasing order is a
 * valid fill order. Matrix stores its backing array column-major
 * (array_p[i + rows*j], see Matrix::operator()), matching the idx = i +
 * n_rows*j indexing shared by all four matrices here.
 *
 * Column banding: a plain column-major nested loop (finish column j, then
 * column j+1, ...) is cache-friendly, but leaves only one thing in flight at
 * a time -- the "up" move's loop-carried dependency, cell (i,j) reading
 * score[(i-1,j)] written by the previous loop iteration. With nothing
 * independent to overlap it with, the core stalls on that chain's load-use
 * latency for the full column height, once per column (measured: this alone
 * dropped IPC from ~3.2 to ~2.0 versus the old expanding-square traversal,
 * more than offsetting its lower instruction count and far fewer L1 misses).
 * Interleaving BAND_WIDTH adjacent columns -- rows outermost, columns
 * innermost within the band -- turns each column's "up" chain into a
 * software-pipelined recurrence: by the time row i's cell for column j needs
 * score[(i-1,j)], BAND_WIDTH-1 other columns' cells have been computed in
 * between, giving that store time to retire. A band's working set (BAND_WIDTH
 * columns times the full column height, across all four matrices) stays
 * small enough to remain cache-resident, so this keeps the column-major
 * locality while restoring the instruction-level parallelism the old
 * expanding-square traversal used to provide incidentally, by interleaving a
 * row-chain and a column-chain.
 *
 * \param int_data_sequence     Prepared (possibly flipped) query sequence, 0-based.
 * \param int_genomic_sequence  Prepared (possibly flipped) reference sequence, 0-based.
 * \param dp  DP workspace whose matrices were already initialized by initialize_sw_matrices.
 */
void fill_sw_score_matrix(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence, SwDPState &dp,
                          const SwDPConfig &config)
{
    const bool reset_negative_scores = config.alignment_mode.is_local_alignment();
    const int n_rows = dp.n_rows;

    double *const score = dp.score_matrix.data();
    int *const row_mem = dp.row_memory_matrix.data();
    int *const col_mem = dp.col_memory_matrix.data();
    int *const numb_trk = dp.alignment_numb_tracker.data();

    // Width of the interleaved column band; see the traversal note above. #define'd rather than
    // constexpr so IGOR_UNROLL below can stringize the same value into a pragma -- keeping the
    // unroll count and the loop bound as a single source of truth.
#define IGOR_SW_BAND_WIDTH 8
    constexpr int BAND_WIDTH = IGOR_SW_BAND_WIDTH;

    // Always start at index 1 since first column and first row are initialization values
    int j = 1;
    for (; j + BAND_WIDTH <= dp.n_cols; j += BAND_WIDTH) {
        for (int i = 1; i != dp.n_rows; ++i) {
            IGOR_UNROLL(IGOR_SW_BAND_WIDTH)
            for (int b = 0; b != BAND_WIDTH; ++b) {
                fill_sw_matrix_cell(i, j + b, n_rows, score, row_mem, col_mem, numb_trk, dp.candidates,
                                    reset_negative_scores, config, int_data_sequence, int_genomic_sequence);
            }
        }
    }
    // Remaining columns too few to fill a whole band.
    for (; j != dp.n_cols; ++j) {
        for (int i = 1; i != dp.n_rows; ++i) {
            fill_sw_matrix_cell(i, j, n_rows, score, row_mem, col_mem, numb_trk, dp.candidates, reset_negative_scores,
                                config, int_data_sequence, int_genomic_sequence);
        }
    }
#undef IGOR_SW_BAND_WIDTH
}

} // namespace swalign

/**
 *\brief Performs Smith-Waterman alignment between two sequences (translated to int sequence as a prior)
 * Output:
 * int:Alignment score
 * Alignment_data: comprises offset, insertions and deletions locations.
 * Note: the gene_name field of the Alignment_data object is left blank and should be completed in a higher level method
 */
list<pair<int, Alignment_data>> sw_align(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                         bool best_only, SwDPConfig config)
{
    using namespace swalign;
    /*Convention:
        - data_sequence is the query, and is the vertical sequence in the matrix (i indexed)
        - genomic_sequence is the reference, and the horizontal sequence in the matrix (j indexed)
        - The alignment matrix and other utilities are of size sequence size + 1. The extra first row/column allows to initialize the algorithm (especially for the score matrix).
    */
    // config is owned by value here, so it can be mutated in place instead of being deep-copied
    // into a separate effective_config (substitution_matrix alone is a ~15x15 double matrix).
    // NOTE: the previous effective_config aggregate-init only listed 7 of SwDPConfig's 8 fields,
    // silently dropping back to enable_extension's default (true) regardless of the caller's
    // setting. That is preserved verbatim here (untested, latent behavior, tracked separately)
    // rather than fixed as a side effect of this perf-only change.
    config.enable_extension = true;
    config.alignment_mode = effective_sw_mode_for_dp(config);
    const SwPreparedInputs prepared_inputs = prepare_sw_inputs(int_data_sequence, int_genomic_sequence, config);
    const int n_rows = static_cast<int>(prepared_inputs.data_sequence.size()) + 1;
    const int n_cols = static_cast<int>(prepared_inputs.genomic_sequence.size()) + 1;

    SwDPState dp(n_rows, n_cols);
    initialize_sw_matrices(dp, config);
    fill_sw_score_matrix(prepared_inputs.data_sequence, prepared_inputs.genomic_sequence, dp, config);

    const SwReconstructionResult reconstruction =
            traceback_sw_alignments(int_data_sequence, int_genomic_sequence, prepared_inputs, dp, config);

    list<pair<int, Alignment_data>> seq_alignments_results = reconstruction.alignments;

    return seq_alignments_results;
}

//Compare alignments (sort by score)
bool align_compare(Alignment_data align1, Alignment_data align2)
{
    return align1.score > align2.score;
}

/**
 * \brief A dumb function to read CSV anchor gene indices
 */
unordered_map<string, size_t> read_gene_anchors_csv(const string &filename, string sep)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename + " in read_gene_anchors_csv()");
    }

    string temp_str;
    unordered_map<string, size_t> anchors_map;

    getline(infile, temp_str); //Ignore first header line

    vector<string> separated_strings;
    bool first_line = true;
    while (getline(infile, temp_str)) {
        separated_strings = extract_string_fields(temp_str, sep);
        if (first_line) {
            if (separated_strings.size() < 2) {
                throw runtime_error(
                        "Expected at least two fields in read_gene_anchors_csv(). Make sure file is separated by:\'"
                        + sep + "\'.");
            }

            try {
                auto unused = stoi(separated_strings.at(1));
            } catch (exception &e) {
                throw runtime_error("Expected an integer for the second field in read_gene_anchors_csv(), received:"
                                    + separated_strings.at(1));
            }

            first_line = false;
        }
        anchors_map.emplace(separated_strings.at(0), stoi(separated_strings.at(1)));
    }

    return anchors_map;
}

/**
 * \brief A dumb function to read CSV template specific offset bounds
 */
unordered_map<string, pair<int, int>> read_template_specific_offset_csv(const string &filename, string sep /*= ";"*/)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: " + filename + " in read_gene_anchors_csv()");
    }
    string temp_str;
    unordered_map<string, pair<int, int>> template_bounds_map;

    getline(infile, temp_str); //Ignore first header line

    vector<string> separated_strings;
    bool first_line = true;
    while (getline(infile, temp_str)) {
        separated_strings = extract_string_fields(temp_str, sep);
        if (first_line) {
            if (separated_strings.size() < 2) {
                throw runtime_error(
                        "Expected at least three fields in read_gene_anchors_csv(). Make sure file is separated by:\'"
                        + sep + "\'.");
            }

            try {
                auto unused = stoi(separated_strings.at(1));
            } catch (exception &e) {
                throw runtime_error(
                        "Expected an integer for the min offset (second field) in read_gene_anchors_csv(), received:"
                        + separated_strings.at(1));
            }

            try {
                auto unused = stoi(separated_strings.at(2));
            } catch (exception &e) {
                throw runtime_error(
                        "Expected an integer for the max offset (third field) in read_gene_anchors_csv(), received:"
                        + separated_strings.at(2));
            }

            first_line = false;
        }
        template_bounds_map.emplace(separated_strings.at(0),
                                    make_pair(stoi(separated_strings.at(1)), stoi(separated_strings.at(2))));
    }

    return template_bounds_map;
}

/**
 * Function reading a substitution matrix from a file
 * The matrix should be 4*4 or 15*15
 * Matrix can have a header or not
 */
Matrix<double> read_substitution_matrix(const string &filename, string sep /*=","*/)
{
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("File not found: \"" + filename + "\" in read_substitution_matrix()");
    }
    vector<double> tmp_vect;
    string temp_str;
    bool first_line = true;
    size_t line_size = 0;
    while (getline(infile, temp_str)) {
        vector<string> line_vect = extract_string_fields(temp_str, sep);
        if (first_line) {
            line_size = line_vect.size();
            if ((line_size != 4) and (line_size != 15)) {
                throw runtime_error("Substitution matrix should be 4*4 (A,C,G,T) or 15*15 "
                                    "(A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N) in read_substitution_matrix()");
            }
            first_line = false;
        }
        if (line_vect.size() != line_size) {
            throw runtime_error(
                    "Substitution matrix' rows length are inconsistent in input matrix for read_substitution_matrix()");
        }
        for (const string &field : line_vect) {
            tmp_vect.emplace_back(stod(field));
        }
    }
    if (line_size == 0 or (tmp_vect.size() / line_size != line_size)) {
        throw runtime_error("Substitution matrix' number of rows and columns are inconsistent in input matrix for "
                            "read_substitution_matrix()");
    }
    return Matrix<double>(line_size, line_size, tmp_vect);
}

/**
 * \brief Compute min and max offsets over all genomic templates and check if they are constant.
 * \author Q.Marcou
 * \version 1.2.0
 *
 * \param [in] genomic_offset_bounds A hash map containing offsets lower and upper bounds for each genomic template. Keys of the map are the genomic templates names.
 * \return A three component tuple: first: boolean true if offsets are template specific (not all the same for all templates, second and third the min and max offsets over all templates.
 */
tuple<bool, int, int>
extract_min_max_genomic_templates_offsets(const unordered_map<string, pair<int, int>> &genomic_offset_bounds)
{
    //Compute min and max offsets over all genomic templates and check if they are constant.
    int min_offset = INT32_MAX;
    int max_offset = INT32_MIN;
    bool template_specific_offsets =
            false; //Although silly we have to recover the fact that not all templates have the same offset.
    for (const pair<string, pair<int, int>> template_bounds : genomic_offset_bounds) {
        if (min_offset > template_bounds.second.first) {
            if (not template_specific_offsets and min_offset != INT32_MAX) {
                //If min_offset's values is no longer UINT64_MAX, it has been updated once and this is the second => templates have different bounds
                template_specific_offsets = true;
            }
            min_offset = template_bounds.second.first;
        }

        if (max_offset < template_bounds.second.second) {
            if (not template_specific_offsets and max_offset != INT32_MIN) {
                //If max_offset's values is no longer -UINT64_MAX, it has been updated once and this is the second => templates have different bounds
                template_specific_offsets = true;
            }
            max_offset = template_bounds.second.second;
        }
    }
    return tuple<bool, int, int>(template_specific_offsets, min_offset, max_offset);
}

/**
 * \brief Extract alignments of the gene/allele with best alignment score.
 * \author Q.Marcou
 * \version 1.2.1
 *
 * \param [in] all_aligns A forward list containing alignments for several genes/alleles.
 * \return A forward list containing only alignments of the best gene/allele candidate.
 */
forward_list<Alignment_data> extract_best_gene_alignments(const forward_list<Alignment_data> &all_aligns)
{
    //Find gene/allele whose alignments has best score
    set<string> best_genes_names;
    double best_align_score = 0;
    for (Alignment_data alignment : all_aligns) {
        if (alignment.score > best_align_score) {
            best_align_score = alignment.score;
            best_genes_names.clear();
            best_genes_names.emplace(alignment.gene_name);
        } else if (alignment.score == best_align_score) {
            //No need to make sure gene name is not already contained in the set before adding it, since emplace already check this
            best_genes_names.emplace(alignment.gene_name);
        }
    }

    //Extract needed alignments
    forward_list<Alignment_data> best_gene_aligns;
    for (Alignment_data alignment : all_aligns) {
        if (best_genes_names.count(alignment.gene_name) > 0) {
            best_gene_aligns.emplace_front(alignment);
        }
    }
    return best_gene_aligns;
}

// ============================================================================
// Alignment_data computed getter method implementations
// ============================================================================

std::vector<int> Alignment_data::get_core_mismatches() const {
    std::vector<int> core;
    for (int pos : mismatches) {
        if (pos >= five_p_offset && pos <= three_p_offset) {
            core.push_back(pos);
        }
    }
    return core;
}

std::vector<int> Alignment_data::get_5p_extended_mismatches() const {
    std::vector<int> result;
    for (int pos : extended_mismatches) {
        if (pos < static_cast<int>(five_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<int> Alignment_data::get_3p_extended_mismatches() const {
    std::vector<int> result;
    for (int pos : extended_mismatches) {
        if (pos > static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<int> Alignment_data::get_core_insertions() const
{
    std::vector<int> result;
    for (int pos : get_all_insertions()) {
        if (pos >= static_cast<int>(five_p_offset) && pos <= static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<int> Alignment_data::get_core_deletions() const
{
    std::vector<int> result;
    size_t n_ins_5p = get_5p_extended_insertions().size();
    size_t n_ins_all = get_5p_extended_insertions().size() + get_core_insertions().size();
    int ref_start = static_cast<int>(five_p_offset) - offset - n_ins_5p;
    int ref_end = static_cast<int>(three_p_offset) - offset - n_ins_all;
    for (int pos : get_all_deletions()) {
        if (pos < ref_start) {
            ref_start++;
            ref_end++;
        } else if (pos > ref_start && pos <= ref_end) {
            result.push_back(pos);
            ref_end++; // Account for reference coordinate shift
        } else {
            // Beyond core alignement
            break;
        }
    }
    return result;
}

// Extended insertion accessors (query coordinates)
std::vector<int> Alignment_data::get_5p_extended_insertions() const {
    std::vector<int> result;
    for (int pos : get_all_insertions()) {
        if (pos < static_cast<int>(five_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

std::vector<int> Alignment_data::get_3p_extended_insertions() const {
    std::vector<int> result;
    for (int pos : get_all_insertions()) {
        if (pos > static_cast<int>(three_p_offset)) {
            result.push_back(pos);
        }
    }
    return result;
}

// Extended deletion accessors (reference coordinates with dynamic threshold adjustment)
std::vector<int> Alignment_data::get_5p_extended_deletions() const
{
    std::vector<int> result;
    size_t n_ins_5p = get_5p_extended_insertions().size();
    int ref_start = static_cast<int>(five_p_offset) - offset - n_ins_5p;
    for (int pos : get_all_deletions()) {
        if (pos <= ref_start) {
            result.push_back(pos);
            ref_start++; // Account for reference coordinate shift
        }
    }
    return result;
}

std::vector<int> Alignment_data::get_3p_extended_deletions() const {
    std::vector<int> result;
    size_t n_ins_all = get_5p_extended_insertions().size() + get_core_insertions().size();
    int ref_end = static_cast<int>(three_p_offset) - offset - n_ins_all;
    
    // Second pass: collect 3' extended deletions (after adjusted ref_end)
    for (int pos : get_all_deletions()) {
        if (pos > ref_end) {
            result.push_back(pos);
        }
        else{
            ref_end++;
        }
    }
    return result;
}

std::string Alignment_data::core_cigar() const {
    if (query_length == 0 || germline_length == 0) {
        throw std::logic_error("query_length and germline_length required for CIGAR generation");
    }
    return alignment_data_to_core_cigar(*this, query_length, germline_length);
}

std::string Alignment_data::extended_cigar() const {
    if (query_length == 0 || germline_length == 0) {
        throw std::logic_error("query_length and germline_length required for CIGAR generation");
    }
    return alignment_data_to_extended_cigar(*this, query_length, germline_length);
}

bool Alignment_data::validate() const {
    // Check that mismatches are sorted
    for (size_t i = 1; i < mismatches.size(); ++i) {
        if (mismatches[i-1] > mismatches[i]) {
            return false; // Not sorted
        }
    }

    // Check that insertions/deletions are sorted (relied upon by
    // get_core_deletions/get_5p_extended_deletions/get_3p_extended_deletions)
    for (size_t i = 1; i < insertions.size(); ++i) {
        if (insertions[i-1] > insertions[i]) {
            return false; // Not sorted
        }
    }
    for (size_t i = 1; i < deletions.size(); ++i) {
        if (deletions[i-1] > deletions[i]) {
            return false; // Not sorted
        }
    }

    // Check that extended_mismatches are sorted
    for (size_t i = 1; i < extended_mismatches.size(); ++i) {
        if (extended_mismatches[i-1] > extended_mismatches[i]) {
            return false; // Not sorted
        }
    }
    
    // Check that extended_mismatches are actually outside core alignment
    for (int pos : extended_mismatches) {
        if (pos >= static_cast<int>(five_p_offset) && pos <= static_cast<int>(three_p_offset)) {
            return false; // Extended mismatch is within core alignment
        }
    }
    
    // Check that all mismatches in extended_mismatches are also in mismatches
    std::unordered_set<int> all_mismatches_set(mismatches.begin(), mismatches.end());
    for (int pos : extended_mismatches) {
        if (all_mismatches_set.count(pos) == 0) {
            return false; // Extended mismatch not in main mismatches
        }
    }
    
    // Check that core mismatches + extended mismatches = all mismatches
    std::vector<int> core = get_core_mismatches();
    std::vector<int> combined = core;
    combined.insert(combined.end(), extended_mismatches.begin(), extended_mismatches.end());
    std::sort(combined.begin(), combined.end());
    if (combined != mismatches) {
        return false; // Inconsistent mismatch categorization
    }
    
    return true;
}
