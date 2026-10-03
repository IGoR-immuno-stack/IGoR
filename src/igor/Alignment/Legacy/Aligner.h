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
#include <omp.h>
#include <stdexcept>
#include <random>
#include <chrono>
#include <limits>

#include <igor/Core/Legacy/IntStr.h>
#include <igor/Core/Legacy/Utils.h>
#include <igor/Core/Legacy/AlignmentData.h>
#include <igor/Alignment/Export.h>
#include <igor/Alignment/Legacy/TestingExport.h>

namespace igor::alignment::legacy {
using namespace igor::core::legacy;

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

class ALIGNMENT_EXPORT Aligner
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
    int incorporate_in_dels(std::string &, std::string &, const std::vector<size_t>, const std::vector<size_t>,
                            int);

    // Configuration for alignment extension
    void set_enable_extension(bool enable) { enable_extension_ = enable; }
    bool get_enable_extension() const { return enable_extension_; }

private:
    std::forward_list<std::pair<std::string, std::string>> nt_genomic_sequences;
    std::forward_list<std::pair<std::string, Int_Str>> int_genomic_sequences;
    Matrix<double> substitution_matrix;
    int gap_penalty;
    Gene_class gene;
    bool enable_extension_ = true; // Enable/disable alignment extension for capturing mismatches in extended regions
    std::unordered_map<std::string, std::pair<int, int>> build_genomic_bounds_map(int, int) const;
};

ALIGNMENT_EXPORT std::vector<std::pair<int, char>> parse_cigar(const std::string &cigar);
ALIGNMENT_EXPORT std::string alignment_data_to_core_cigar(const Alignment_data &aln);
ALIGNMENT_EXPORT std::string alignment_data_to_core_cigar(const Alignment_data &aln, size_t sequence_length,
                                                           size_t germline_length);
ALIGNMENT_EXPORT std::string alignment_data_to_extended_cigar(const Alignment_data &aln, size_t sequence_length,
                                                           size_t germline_length);
ALIGNMENT_EXPORT Alignment_data alignment_data_from_cigar(const std::string &gene_name, const std::string &cigar,
                                                     int seq_start_1based, int seq_end_1based, int ref_start_1based,
                                                     int ref_end_1based, double score);
ALIGNMENT_EXPORT Alignment_data alignment_data_from_cigar_and_extended(const std::string &gene_name, const std::string &core_cigar,
                                                                 const std::string &extended_cigar, double score);
ALIGNMENT_EXPORT int alignment_data_sequence_start(const Alignment_data &aln);
ALIGNMENT_EXPORT int alignment_data_sequence_end(const Alignment_data &aln);
ALIGNMENT_EXPORT int alignment_data_germline_start(const Alignment_data &aln);
ALIGNMENT_EXPORT int alignment_data_germline_end(const Alignment_data &aln);
// Standalone function for external alignment import
std::vector<size_t> extend_alignment_mismatches(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                        const Alignment_data aln);
ALIGNMENT_EXPORT bool alignment_data_equal(const Alignment_data &a, const Alignment_data &b, double score_tolerance = 1e-9);

ALIGNMENT_EXPORT std::pair<int, Alignment_data> parse_single_alignment_csv_line(const std::string &line);
ALIGNMENT_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv(const std::string &, Gene_class, double, bool,
                        const std::vector<std::pair<const int, const std::string>> &);
ALIGNMENT_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv(
        const std::string &, Gene_class, double, bool, const std::vector<std::pair<const int, const std::string>> &,
        std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
ALIGNMENT_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv_score_range(const std::string &, Gene_class, double, bool,
                                    const std::vector<std::pair<const int, const std::string>> &);
ALIGNMENT_EXPORT std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
read_alignments_seq_csv_score_range(
        const std::string &, Gene_class, double, bool, const std::vector<std::pair<const int, const std::string>> &,
        std::unordered_map<int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
ALIGNMENT_EXPORT std::vector<std::tuple<int, std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
        map2vect(std::unordered_map<
                 int, std::pair<std::string, std::unordered_map<Gene_class, std::vector<Alignment_data>>>>);
ALIGNMENT_EXPORT std::forward_list<std::pair<const int, const std::string>> read_indexed_seq_csv(const std::string &);
ALIGNMENT_EXPORT std::vector<std::pair<const int, const std::string>> read_indexed_csv(const std::string &);
ALIGNMENT_EXPORT std::vector<std::pair<const int, const std::string>> read_fasta(const std::string &);
ALIGNMENT_EXPORT std::vector<std::pair<std::string, std::string>> read_genomic_fasta(const std::string &);
ALIGNMENT_EXPORT std::vector<std::pair<const int, const std::string>> read_txt(const std::string &);
ALIGNMENT_EXPORT std::unordered_map<std::string, size_t> read_gene_anchors_csv(const std::string &,
                                                                          std::string separator = ";");
ALIGNMENT_EXPORT std::unordered_map<std::string, std::pair<int, int>>
read_template_specific_offset_csv(const std::string &, std::string separator = ";");
ALIGNMENT_EXPORT void write_indexed_seq_csv(const std::string &,
                                       const std::vector<std::pair<const int, const std::string>> &);
ALIGNMENT_EXPORT Int_Str nt2int(const std::string &);
ALIGNMENT_EXPORT bool comp_nt_int(const int &, const int &);
ALIGNMENT_EXPORT std::list<Int_nt> get_ambiguous_nt_list(const Int_nt &);
ALIGNMENT_EXPORT inline void write_single_seq_alignment(std::ofstream &, int, std::forward_list<Alignment_data>);
//Compare alignments (sort by score)
ALIGNMENT_EXPORT bool align_compare(Alignment_data, Alignment_data);
ALIGNMENT_EXPORT std::vector<std::pair<const int, const std::string>>
sample_indexed_seq(const std::vector<std::pair<const int, const std::string>> &, const size_t);
ALIGNMENT_EXPORT Matrix<double> read_substitution_matrix(const std::string &, std::string sep = ",");
ALIGNMENT_EXPORT std::tuple<bool, int, int> extract_min_max_genomic_templates_offsets(
        const std::unordered_map<std::string, std::pair<int, int>> &genomic_offset_bounds);
ALIGNMENT_EXPORT std::forward_list<Alignment_data> extract_best_gene_alignments(const std::forward_list<Alignment_data> &);

struct SwAlignmentMode
{
    bool data_leading_free = false;
    bool data_trailing_free = false;
    bool genomic_leading_free = false;
    bool genomic_trailing_free = false;
    bool reverse_sequences = false;

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
 * best_only        Retain only the alignment(s) reaching the best score for this call.
 * min_offset       Lower bound on the offset (genomic-vs-query position).
 * max_offset       Upper bound on the offset.
 * alignment_mode    Boundary and orientation policy for the DP run.
 */
struct SwDPConfig
{
    double score_threshold = -std::numeric_limits<double>::infinity();
    bool best_only = false;
    int min_offset = INT32_MIN;
    int max_offset = INT32_MAX;
    Matrix<double> substitution_matrix;
    int gap_penalty;
    SwAlignmentMode alignment_mode;
    bool enable_extension = true; // Enable/disable alignment extension for capturing mismatches in extended regions
};
ALIGNMENT_TESTING_EXPORT std::list<std::pair<int, Alignment_data>> sw_align(const Int_Str &, const Int_Str &, SwDPConfig);

namespace swalign {

// Forward declare internal structs
struct SwDPState;
struct SwPreparedInputs;

ALIGNMENT_TESTING_EXPORT SwPreparedInputs prepare_sw_inputs(const Int_Str &int_data_sequence, const Int_Str &int_genomic_sequence,
                                   const SwDPConfig &config);
// Coordinate conversion functions
size_t convert_matrix_row_to_query_pos(size_t i, size_t data_seq_size, bool flip_seqs);
size_t convert_matrix_col_to_ref_pos(size_t j, size_t genomic_seq_size, bool flip_seqs);
int convert_matrix_coords_to_offset(int i, int j, size_t data_seq_size, size_t genomic_seq_size, int offset_change,
                                    bool flip_seqs);
ALIGNMENT_TESTING_EXPORT void fill_sw_score_matrix(const Int_Str &, const Int_Str &, SwDPState &, const SwDPConfig &);

// Alignment extension functions for capturing mismatches in extended regions
std::vector<size_t> ungapped_extend_align_5p_from_dp(const SwPreparedInputs &prepared, int i_start, int j_start,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols);

std::vector<size_t> ungapped_extend_align_3p_from_dp(const SwPreparedInputs &prepared, int i_end, int j_end,
                                             size_t data_seq_size, size_t genomic_seq_size, bool flip_seqs,
                                             int matrix_n_rows, int matrix_n_cols);

// Helper functions for merging and sorting mismatches
std::vector<size_t> merge_and_sort_mismatches(const std::vector<size_t> &core_mismatches, const std::vector<size_t> &extended_mismatches);
std::vector<size_t> merge_and_sort_mismatches(const std::vector<size_t> &core_mismatches, const std::vector<size_t> &extended_5p_mismatches,
                                      const std::vector<size_t> &extended_3p_mismatches);

} // namespace swalign

} // namespace igor::alignment::legacy
