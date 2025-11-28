/*
 * ScoreView.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#include "ScoreView.h"

namespace igor {

ScoreView::ScoreView(double *data, size_t num_rows, size_t num_genes)
    : data_(data), num_rows_(num_rows), num_genes_(num_genes), matrix_(data, num_rows, num_genes) {}

ScoreView::ScoreMatrix ScoreView::get_matrix() const { return matrix_; }

double ScoreView::get_score(size_t row_idx, size_t gene_idx) const {
    // mdspan uses operator[] with multiple arguments in C++23
    // For the backport, it might use operator() or operator[] depending on
    // version Kokkos mdspan usually supports operator[]
    return matrix_(row_idx, gene_idx);
}

}  // namespace igor
