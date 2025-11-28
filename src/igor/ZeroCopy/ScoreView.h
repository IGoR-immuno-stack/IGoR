/*
 * ScoreView.h
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#pragma once

#include <igor/ZeroCopy/Export.h>

#include <experimental/mdspan>
#include <sparrow/record_batch.hpp>
#include <string>
#include <vector>

namespace igor {

/**
 * @brief Zero-copy view of a score grid using std::mdspan
 *
 * This class demonstrates how to use std::mdspan to provide a multi-dimensional
 * view over a flat Arrow buffer. It is designed to handle gene alignment scores
 * stored in a columnar format.
 *
 * It supports two layouts:
 * 1. Column-major (Separate columns for V, D, J scores)
 * 2. Contiguous (Single FixedSizeList column for all scores)
 */
class ZEROCOPY_EXPORT ScoreView {
public:
    // Define the mdspan type: 2D view over double, dynamic extent
    using ScoreMatrix = std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>>;

    /**
     * @brief Create a view from separate score columns (V, D, J)
     *
     * Note: Since Arrow columns are separate buffers, we can't create a single
     * contiguous mdspan over all of them. This constructor creates a "logical"
     * view where we might need to copy data to a contiguous buffer first,
     * OR we use a more complex strided layout if supported.
     *
     * For this POC, we will demonstrate the "Contiguous" use case (Use Case B)
     * where we assume the input is a single flat buffer (e.g. from a
     * FixedSizeList or a specific "scores" column).
     */
    ScoreView(double *data, size_t num_rows, size_t num_genes);

    /**
     * @brief Get the mdspan view
     */
    ScoreMatrix get_matrix() const;

    /**
     * @brief Get a score at a specific row and gene index
     */
    double get_score(size_t row_idx, size_t gene_idx) const;

private:
    double *data_;
    size_t num_rows_;
    size_t num_genes_;
    ScoreMatrix matrix_;
};

}  // namespace igor
