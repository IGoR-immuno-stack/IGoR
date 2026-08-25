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
struct CORE_TESTING_EXPORT SwCandidate
{
    int score;
    // Running coordinate of alignment end.
    int row;
    int col;

    // Coordinates (1-based) of the cell where this candidate was first seeded ( Fixed for the candidate's whole lifetime, unlike row/col which track the
    // running-max position.
    const int start_row;
    const int start_col;

    SwCandidate &operator=(const SwCandidate &other)
    {
        score = other.score;
        row = other.row;
        col = other.col;
        return *this;
    }
};

/**
 * Column-indexed row range(s) restricting which DP cells are worth filling, derived in closed
 * form from an SwDPConfig's min_offset/max_offset/gap_penalty/substitution_matrix/score_threshold
 * -- see compute_band_bounds's doc comment for the derivation.
 *
 * lo/hi are sized n_cols; column 0 (the init boundary, always fully valid) is left at its default
 * {lo=1,hi=0} and never consulted. hi[j] < lo[j] means column j is empty: no cell in that column
 * can be part of any alignment admissible under min_offset/max_offset/score_threshold, so the
 * whole column can be skipped.
 */
struct CORE_TESTING_EXPORT SwBandBounds
{
    std::vector<int> lo;
    std::vector<int> hi;
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
struct CORE_TESTING_EXPORT SwDPState
{
    int n_rows;
    int n_cols;
    Matrix<double> score_matrix;
    Matrix<int> row_memory_matrix;
    Matrix<int> col_memory_matrix;
    Matrix<int> alignment_numb_tracker;
    std::vector<SwCandidate> candidates;
    // Column-indexed row bounds fill_sw_score_matrix restricts its fill to (see SwBandBounds).
    // Defaults to the full [1,n_rows-1] range for every column, i.e. unbanded, so constructing an
    // SwDPState without an explicit SwBandBounds (as existing whitebox tests do) behaves exactly
    // as before this field was added.
    std::vector<int> lo;
    std::vector<int> hi;

    SwDPState(int nr, int nc)
        : n_rows(nr),
          n_cols(nc),
          score_matrix(nr, nc),
          row_memory_matrix(nr, nc),
          col_memory_matrix(nr, nc),
          alignment_numb_tracker(nr, nc),
          lo(nc, 1),
          hi(nc, nr - 1)
    {
    }

    SwDPState(int nr, int nc, SwBandBounds bounds)
        : n_rows(nr),
          n_cols(nc),
          score_matrix(nr, nc),
          row_memory_matrix(nr, nc),
          col_memory_matrix(nr, nc),
          alignment_numb_tracker(nr, nc),
          lo(std::move(bounds.lo)),
          hi(std::move(bounds.hi))
    {
    }
};

/**
 * Prepared (possibly sequence-reversed) integer-coded inputs for one sw_align call.
 * See swalign::prepare_sw_inputs.
 */
struct CORE_TESTING_EXPORT SwPreparedInputs
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
CORE_TESTING_EXPORT void initialize_sw_matrices(SwDPState &dp, const SwDPConfig &config);

/**
 * Compute, for each column of the (n_rows x n_cols) prepared DP matrix, the contiguous range of
 * rows that can possibly belong to an alignment whose offset would pass config.min_offset/
 * max_offset (and, if set, config.score_threshold) -- i.e. the "cone" a banded fill can safely
 * restrict itself to, derived purely from the scoring scheme rather than a fixed diagonal +
 * tolerance.
 *
 * The math (see the per-call closed-form derivation in Aligner.cpp): a path using only diagonal
 * (score up to max_match_score per step) and gap (cost gap_penalty per step) moves can reach cell
 * (i,j) from a seed on diagonal d0 with score at most
 *   m * min(i - max(0,d0), j - max(0,-d0)) - g * |(i-j) - d0|
 * (front-load all diagonal moves, pay gap cost for the rest) -- a sound upper bound regardless of
 * where the true seed sits on that diagonal.
 *
 * Critically, this "forward" bound alone is NOT a valid pruning criterion on its own for
 * semi-global (non-reset) alignment modes (V/D/J-gene): unlike vanilla local alignment, a
 * candidate's tracked score is never reset when it would go negative, so a real candidate's
 * eventual peak can occur arbitrarily far downstream of a cell whose own best-case score is
 * unremarkable or negative -- excluding such a cell outright would silently drop legitimate
 * alignments. The correct criterion combines the forward bound with a "backward" bound -- the
 * best additional score achievable continuing from (i,j) to the matrix's far corner, extending
 * (i,j)'s own diagonal (the cheapest possible continuation, zero further gaps):
 *   Combined(i,j) = seed_reach_score(i,j,d0) + m * min(n_rows-1-i, n_cols-1-j)
 * compared against the *actual* config.score_threshold (not a hardcoded 0). If Combined(i,j) is
 * below score_threshold, no real path through (i,j) -- for any real sequence content -- could
 * ever contribute to a candidate whose eventual peak reaches threshold, so (i,j) is safe to skip.
 * With score_threshold left at -infinity (SwDPConfig's default), Combined(i,j) is always finite
 * and the comparison always holds, correctly degrading to "no score-based exclusion at all" --
 * only geometric reachability from an admissible offset diagonal restricts the band in that case
 * (a simple quadrant/rectangle cut, not a narrowing cone). With a real, finite threshold, the
 * combined bound genuinely narrows near the seed (too little accumulated yet) and again near the
 * far corner (too little room left to make up a bad start), producing the widening-then-plateauing
 * cone shape the offset/score_threshold restriction is meant to give.
 *
 * Both the forward and backward terms are concave/unimodal in i (for fixed j) and in d0 (for
 * fixed i,j, once optimally chosen), which is what makes an O(n_rows + n_cols) per-column sweep
 * possible instead of an O(n_rows * n_cols) scan.
 *
 * Critical asymmetry (verified against traceback_sw_alignments's own offset computation): when
 * flip_seqs is false (V/D-gene alignment), the traceback computes offset from the alignment's
 * *seed* (start) cell, so min_offset/max_offset restrict the seed's diagonal directly ("prefix
 * cone"): the forward term is offset-restricted, the backward term unrestricted. When flip_seqs is
 * true (J-gene alignment), offset is instead computed from the cell where the candidate's
 * *running-max* score was recorded, which can drift arbitrarily far from an *unrestricted* seed --
 * so min_offset/max_offset instead restrict the diagonal of that downstream cell ("suffix cone"):
 * the backward term is offset-restricted, the forward term unrestricted. These are not mirror
 * images of one another: which side is offset-restricted swaps with flip_seqs.
 *
 * \param n_rows            Prepared DP matrix row count (query length + 1).
 * \param n_cols             Prepared DP matrix column count (template length + 1).
 * \param min_offset        Same field as SwDPConfig::min_offset.
 * \param max_offset        Same field as SwDPConfig::max_offset.
 * \param max_match_score   Best (highest) entry of the substitution matrix in use for this call.
 * \param gap_penalty       Same field as SwDPConfig::gap_penalty.
 * \param score_threshold   Same field as SwDPConfig::score_threshold; -infinity disables this axis.
 * \param flip_seqs         config.alignment_mode.reverse_sequences, AFTER effective_sw_mode_for_dp.
 * \param data_seq_size     Original (unflipped) query sequence length.
 * \param genomic_seq_size  Original (unflipped) reference sequence length.
 */
CORE_TESTING_EXPORT SwBandBounds compute_band_bounds(int n_rows, int n_cols, int min_offset, int max_offset,
                                                     double max_match_score, int gap_penalty,
                                                     double score_threshold, bool flip_seqs, size_t data_seq_size,
                                                     size_t genomic_seq_size);

} // namespace swalign
