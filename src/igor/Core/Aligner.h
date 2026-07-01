/*
 * Aligner.h
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

#pragma once

#include <forward_list>
#include <list>
#include <unordered_map>
#include <set>
#include <utility>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <igor/Core/Utils.h>
#include <omp.h>
#include <stdexcept>
#include <random>
#include <chrono>

#include <igor/Core/IntStr.h>

#include <igorCoreExport.h>

/**
 * \class Alignment_data Aligner.h
 * \brief Stores information on the alignment of one genomic template against the target
 * \author Q.Marcou
 * \version 1.0
 *
 * Stores information on the alignment of one genomic template against the target.
 * It contains:
 * - the gene name
 * - the offset of the alignment (index on the target sequence on which the first letter of the FULL genomic template aligns (can be negative or lie outside the target)
 * - 5' and 3' offset positions of the best alignment first and last aligned nucleotide
 * - insertions : indices on the TARGET of inserted nucleotides
 * - deletions : indices on the GENOMIC TEMPLATE of deleted nucleotides
 * - alignment length
 * - list of mismatches (that lie even outside the best alignment to allow IGoR to know mismatch positions in advance while exploring different deletions numbers)
 *   NOTE: The mismatches vector may contain mismatches beyond the five_p_offset to three_p_offset range,
 *   representing mismatches in the extended alignment regions (e.g., in deleted V/J nucleotides).
 *   The vector is SORTED and contains ALL mismatches (both within the core alignment and in extended regions).
 * - the alignment score
 */
struct Alignment_data
{
    std::string gene_name;
    int offset;
    size_t five_p_offset;
    size_t three_p_offset;
    std::forward_list<int> insertions; //gap in the genomic sequence
    std::forward_list<int> deletions; //gap in the data sequence
    size_t align_length;
    mutable std::vector<int> mismatches;
    double score;

    Alignment_data(std::string gene, int off)
        : gene_name(gene),
          offset(off),
          insertions(*(new std::forward_list<int>)),
          deletions(*(new std::forward_list<int>)),
          score(0)
    {
    }
    Alignment_data(int off, size_t five_p_off, size_t three_p_off, size_t align_len, std::forward_list<int> ins,
                   std::forward_list<int> del, std::vector<int> mis, double alignment_score)
        : gene_name(std::string()),
          offset(off),
          five_p_offset(five_p_off),
          three_p_offset(three_p_off),
          insertions(ins),
          deletions(del),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score)
    {
    }
    Alignment_data(std::string gene, int off, size_t align_len, std::forward_list<int> ins, std::forward_list<int> del,
                   std::vector<int> mis, double alignment_score)
        : gene_name(gene),
          offset(off),
          insertions(ins),
          deletions(del),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score)
    {
    }
    Alignment_data(std::string gene, int off, size_t five_p_off, size_t three_p_off, size_t align_len,
                   std::forward_list<int> ins, std::forward_list<int> del, std::vector<int> mis, double alignment_score)
        : gene_name(gene),
          offset(off),
          five_p_offset(five_p_off),
          three_p_offset(three_p_off),
          insertions(ins),
          deletions(del),
          align_length(align_len),
          mismatches(mis),
          score(alignment_score)
    {
    }

    /*	bool operator<(const Alignment_data& align){
		//Hardcode to get the alignments in descending order using sort()
		return this->score > align.score;
	}*/
};

/**
 * \class Aligner Aligner.h
 * \brief A modified Smith-Waterman alignment class
 * \author Q.Marcou
 * \version 1.0
 *
 * The Aligner class allows to perform SW alignments according to the parameters (substitution matrix,gap penalty) supplied upon construction of the object.
 * The SW alignments matrix has been altered for V and J in order to allow for deletions on the deleted side only.
 * Alignments can be made in parallel using openMP
 *
 */

class CORE_EXPORT Aligner
{
public:
    Aligner();
    Aligner(Matrix<double>, int, Gene_class);
    virtual ~Aligner();

    // Single sequence alignments methods
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, int, int, bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, int, int, std::set<std::string>,
                                                bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, bool, int, int, bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, bool, int, int, std::set<std::string>,
                                                bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool,
                                                std::unordered_map<std::string, std::pair<int, int>>, bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool,
                                                std::unordered_map<std::string, std::pair<int, int>>,
                                                std::set<std::string>, bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, bool,
                                                std::unordered_map<std::string, std::pair<int, int>>, bool = false);
    std::forward_list<Alignment_data> align_seq(std::string, double, bool, bool,
                                                std::unordered_map<std::string, std::pair<int, int>>,
                                                std::set<std::string>, bool = false);

    // Multiple sequences alignments methods
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool);
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool, bool);
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool, int, int, bool = false);
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool, bool, int, int, bool = false);
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool,
               std::unordered_map<std::string, std::pair<int, int>>, bool = false);
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(std::vector<std::pair<const int, const std::string>>, double, bool, bool,
               std::unordered_map<std::string, std::pair<int, int>>, bool = false);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool, bool);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool, int, int,
                    bool = false);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool, bool, int, int,
                    bool = false);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool,
                    std::unordered_map<std::string, std::pair<int, int>>, bool = false);
    void align_seqs(std::string, std::vector<std::pair<const int, const std::string>>, double, bool, bool,
                    std::unordered_map<std::string, std::pair<int, int>>, bool = false);

    //I/O related methods
    void write_alignments_seq_csv(std::string, std::unordered_map<int, std::forward_list<Alignment_data>>);
    std::unordered_map<int, std::forward_list<Alignment_data>> read_alignments_seq_csv(std::string, double, bool);

    void set_genomic_sequences(std::vector<std::pair<std::string, std::string>>);
    int incorporate_in_dels(std::string &, std::string &, const std::forward_list<int>, const std::forward_list<int>,
                            int);

private:
    std::forward_list<std::pair<std::string, std::string>> nt_genomic_sequences;
    std::forward_list<std::pair<std::string, Int_Str>> int_genomic_sequences;
    Matrix<double> substitution_matrix;
    int gap_penalty;
    Gene_class gene;
    std::unordered_map<std::string, std::pair<int, int>> build_genomic_bounds_map(int, int) const;
};

CORE_EXPORT std::vector<std::pair<int, char>> parse_cigar(const std::string &cigar);
CORE_EXPORT std::string alignment_data_to_core_cigar(const Alignment_data &aln);
CORE_EXPORT std::string alignment_data_to_core_cigar(const Alignment_data &aln, size_t sequence_length,
                                                           size_t germline_length);
CORE_EXPORT std::string alignment_data_to_extended_cigar(const Alignment_data &aln, size_t sequence_length,
                                                           size_t germline_length);
CORE_EXPORT Alignment_data alignment_data_from_cigar(const std::string &gene_name, const std::string &cigar,
                                                     int seq_start_1based, int seq_end_1based, int ref_start_1based,
                                                     int ref_end_1based, double score);
CORE_EXPORT Alignment_data alignment_data_from_cigar_and_extended(const std::string &gene_name, const std::string &core_cigar,
                                                                 const std::string &extended_cigar, double score);
CORE_EXPORT int alignment_data_sequence_start(const Alignment_data &aln);
CORE_EXPORT int alignment_data_sequence_end(const Alignment_data &aln);
CORE_EXPORT int alignment_data_germline_start(const Alignment_data &aln);
CORE_EXPORT int alignment_data_germline_end(const Alignment_data &aln);
// Standalone function for external alignment import
std::vector<int> extend_alignment_mismatches(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                        const Alignment_data aln);
CORE_EXPORT bool alignment_data_equal(const Alignment_data &a, const Alignment_data &b, double score_tolerance = 1e-9);

std::pair<int, Alignment_data> parse_single_alignment_csv_line(const std::string &line);
CORE_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv(const std::string &, Gene_class, double, bool,
                        const std::vector<std::pair<const int, const std::string>> &);
CORE_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv(
        const std::string &, Gene_class, double, bool, const std::vector<std::pair<const int, const std::string>> &,
        std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
CORE_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv_score_range(const std::string &, Gene_class, double, bool,
                                    const std::vector<std::pair<const int, const std::string>> &);
CORE_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv_score_range(
        const std::string &, Gene_class, double, bool, const std::vector<std::pair<const int, const std::string>> &,
        std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
CORE_EXPORT std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
        map2vect(std::unordered_map<
                 int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
CORE_EXPORT std::forward_list<std::pair<const int, const std::string>> read_indexed_seq_csv(const std::string &);
CORE_EXPORT std::vector<std::pair<const int, const std::string>> read_indexed_csv(const std::string &);
CORE_EXPORT std::vector<std::pair<const int, const std::string>> read_fasta(const std::string &);
CORE_EXPORT std::vector<std::pair<std::string, std::string>> read_genomic_fasta(const std::string &);
CORE_EXPORT std::vector<std::pair<const int, const std::string>> read_txt(const std::string &);
CORE_EXPORT std::unordered_map<std::string, size_t> read_gene_anchors_csv(const std::string &,
                                                                          std::string separator = ";");
CORE_EXPORT std::unordered_map<std::string, std::pair<int, int>>
read_template_specific_offset_csv(const std::string &, std::string separator = ";");
CORE_EXPORT void write_indexed_seq_csv(const std::string &,
                                       const std::vector<std::pair<const int, const std::string>> &);
CORE_EXPORT Int_Str nt2int(const std::string &);
CORE_EXPORT bool comp_nt_int(const int &, const int &);
CORE_EXPORT std::list<Int_nt> get_ambiguous_nt_list(const Int_nt &);
CORE_EXPORT inline void write_single_seq_alignment(std::ofstream &, int, std::forward_list<Alignment_data>);
//Compare alignments (sort by score)
CORE_EXPORT bool align_compare(Alignment_data, Alignment_data);
CORE_EXPORT std::vector<std::pair<const int, const std::string>>
sample_indexed_seq(const std::vector<std::pair<const int, const std::string>> &, const size_t);
CORE_EXPORT Matrix<double> read_substitution_matrix(const std::string &, std::string sep = ",");
CORE_EXPORT std::tuple<bool, int, int> extract_min_max_genomic_templates_offsets(
        const std::unordered_map<std::string, std::pair<int, int>> &genomic_offset_bounds);
CORE_EXPORT std::forward_list<Alignment_data> extract_best_gene_alignments(const std::forward_list<Alignment_data> &);

struct SwAlignmentMode
{
    bool data_leading_free;
    bool data_trailing_free;
    bool genomic_leading_free;
    bool genomic_trailing_free;
    bool reverse_sequences;

    bool is_local_alignment() const
    {
        return data_leading_free && data_trailing_free && genomic_leading_free && genomic_trailing_free;
    }
};
/**
 * Run-policy for one Smith-Waterman alignment call.
 *
 * Bundles the scalar parameters and alignment mode that govern how a single
 * sw_align invocation prepares its inputs and filters its results, so they can
 * be passed as a unit instead of several separate arguments.
 *
 * Fields
 * ------
 * score_threshold  Minimum score an alignment must reach to be returned.
 * min_offset       Lower bound on the offset (genomic-vs-query position).
 * max_offset       Upper bound on the offset.
 * alignment_mode    Boundary and orientation policy for the DP run.
 */
struct SwDPConfig
{
    double score_threshold;
    bool best_only;
    int min_offset;
    int max_offset;
    Matrix<double> substitution_matrix;
    int gap_penalty;
    SwAlignmentMode alignment_mode;
};
std::list<std::pair<int, Alignment_data>> sw_align(const Int_Str &, const Int_Str &, bool, const SwDPConfig &);

namespace swalign {

// Forward declare internal structs
struct SwDPState;
struct SwPreparedInputs;

SwPreparedInputs prepare_sw_inputs(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                   const SwDPConfig &config);
// Coordinate conversion functions
size_t convert_matrix_row_to_query_pos(size_t i, size_t data_seq_size, bool flip_seqs);
size_t convert_matrix_col_to_ref_pos(size_t j, size_t genomic_seq_size, bool flip_seqs);
int convert_matrix_coords_to_offset(int i, int j, size_t data_seq_size, size_t genomic_seq_size, int offset_change,
                                    bool flip_seqs);
void fill_sw_score_matrix(const Int_Str &, const Int_Str &, SwDPState &, const SwDPConfig &);
void fill_sw_matrix_cell(const Int_Str &, const Int_Str &, int, int, SwDPState &, const SwDPConfig &);

// Alignment extension functions for capturing mismatches in extended regions
std::vector<int> ungapped_extend_align_5p_from_dp(const SwPreparedInputs &prepared, int i_start, int j_start,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols);

std::vector<int> ungapped_extend_align_3p_from_dp(const SwPreparedInputs &prepared, int i_end, int j_end,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols);

// Helper functions for merging and sorting mismatches
std::vector<int> merge_and_sort_mismatches(const std::vector<int> &core_mismatches, const std::vector<int> &extended_mismatches);
std::vector<int> merge_and_sort_mismatches(const std::vector<int> &core_mismatches, const std::vector<int> &extended_5p_mismatches,
                                      const std::vector<int> &extended_3p_mismatches);

// Extended mismatch identification functions (compute on-demand from existing fields)
std::vector<int> get_5p_extended_mismatches(const Alignment_data& aln);
std::vector<int> get_3p_extended_mismatches(const Alignment_data& aln);
std::vector<int> get_extended_mismatches(const Alignment_data& aln);
std::vector<int> get_core_mismatches(const Alignment_data& aln);
bool validate_mismatch_categorization(const Alignment_data& aln);

} // namespace swalign
