/*
 * AlignerInternal.h
 *
 *  Internal-only Smith-Waterman DP types and helpers, split out of Aligner.cpp
 *  so whitebox tests can construct a swalign::SwDPState and call
 *  swalign::fill_sw_score_matrix directly. Not part of the installed public
 *  API (not listed in Core's FILE_SET HEADERS): Aligner.h only forward-declares
 *  SwDPState/SwPreparedInputs for its own public declarations.
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
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <igor/Core/Aligner.h>
#include <vector>

namespace swalign {

/**
 * One tracked candidate local alignment: its best score so far and the DP matrix
 * coordinate (1-based, +1 padded convention) at which that best score was reached.
 */
struct SwCandidate
{
    int score;
    int row;
    int col;
};

/**
 * Internal workspace for one Smith-Waterman DP execution.
 *
 * Groups the four DP matrices and the candidate-tracking vector that are
 * allocated, mutated, and read together during a single sw_align call.
 *
 * Coordinates: all matrices use the +1 padded convention (row 0 / col 0 are
 * initialization boundaries; sequence positions are 1-based inside the matrix).
 */
struct SwDPState
{
    int n_rows;
    int n_cols;
    Matrix<double> score_matrix;
    Matrix<int> row_memory_matrix;
    Matrix<int> col_memory_matrix;
    Matrix<int> alignment_numb_tracker;
    std::vector<SwCandidate> candidates;

    SwDPState(int nr, int nc)
        : n_rows(nr),
          n_cols(nc),
          score_matrix(nr, nc),
          row_memory_matrix(nr, nc),
          col_memory_matrix(nr, nc),
          alignment_numb_tracker(nr, nc)
    {
    }
};

/**
 * Prepared (possibly sequence-reversed) integer-coded inputs for one sw_align call.
 * See swalign::prepare_sw_inputs.
 */
struct SwPreparedInputs
{
    Int_Str data_sequence;
    Int_Str genomic_sequence;
    int offset_change;
};

/**
 * Initialize score and traceback support matrices for a fresh DP run.
 *
 * \param dp      The DP workspace to initialize. n_rows and n_cols must already be set.
 * \param config  Run policy; alignment_mode's leading/trailing-free flags and gap_penalty are consulted.
 */
void initialize_sw_matrices(SwDPState &dp, const SwDPConfig &config);

} // namespace swalign
