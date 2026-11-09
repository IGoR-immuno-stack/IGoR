/*
 * Matrix.h
 *
 *  Created on: Apr 9, 2015
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
 */

#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Out of Core/Legacy/Utils.h since step 1c of doc/LAYER_REFACTORING_PROPOSAL.md: a matrix is
// numerics, not vocabulary. Kept as is for its legacy users (the aligner, the error rates,
// Dinucl_markov); igor::math::Tensor replaces it as they are promoted.

namespace igor::math::legacy {

/*
 * Declare a simple matrix class with column major data ordering.
 *
 */
template <typename T>
struct Matrix
{
public:
    Matrix() : rows(0), cols(0), array_p(new T[0]) { }
    Matrix(int m, int n) : rows(m), cols(n), array_p(nullptr)
    {
        if (m * n > 0 and m > 0) {
            array_p = new T[m * n];
        }
    }
    Matrix(int m, int n, T arr[]) : rows(m), cols(n), array_p(new T[m * n])
    {
        for (size_t i = 0; i != m * n; i++) {
            array_p[i] = arr[i];
        }
    }
    Matrix(int m, int n, std::vector<T> vect) : rows(m), cols(n), array_p(new T[m * n])
    {
        for (size_t i = 0; i != m * n; i++) {
            array_p[i] = vect.at(i);
        }
    }
    Matrix(const Matrix<T> &other)
    {
        //Provides deep copy of a matrix
        this->rows = other.rows;
        this->cols = other.cols;
        this->array_p = new T[rows * cols];
        for (int i = 0; i != rows * cols; i++) {
            this->array_p[i] = other.array_p[i];
        }
    }
    ~Matrix() { delete[] array_p; }

    Matrix<T> &operator=(const Matrix &other)
    {
        delete[] array_p;
        this->rows = other.rows;
        this->cols = other.cols;
        this->array_p = new T[rows * cols];
        for (int i = 0; i != rows * cols; i++) {
            this->array_p[i] = other.array_p[i];
        }
        return *this;
    }

    Matrix(Matrix &&other) noexcept : rows(other.rows), cols(other.cols), array_p(other.array_p)
    {
        other.rows = 0;
        other.cols = 0;
        other.array_p = nullptr;
    }

    Matrix<T> &operator=(Matrix &&other) noexcept
    {
        if (this != &other) {
            delete[] array_p;
            this->rows = other.rows;
            this->cols = other.cols;
            this->array_p = other.array_p;
            other.rows = 0;
            other.cols = 0;
            other.array_p = nullptr;
        }
        return *this;
    }

    /**
     * \brief Whether (i, j) addresses a cell of this matrix.
     *
     * Exposed rather than left inside the assert so it can be tested without a debug build --
     * the same reason first_unfilled_segment() is a predicate (plan section 7.14).
     *
     * **The lower bound is not redundant.** The indices are `int`, and a caller that derives one
     * by subtraction can hand over a negative value. An upper-bound-only check accepts it, and
     * `array_p[i + rows * j]` then reads from before the allocation -- silently, since every
     * layer above fails to stop it too: an unsigned count wraps to a huge value, the growth
     * check adds to it and wraps back, and the conversion to `int` here brings it out as -1.
     * Plan section 7.18 measured that path returning two different answers for one arithmetic.
     * R5a corrects the derivation that produces the negative index; this is what stops the read.
     */
    bool in_range(const int &i, const int &j) const
    {
        return (i >= 0) && (j >= 0) && (i <= rows - 1) && (j <= cols - 1);
    }

    T &operator()(const int &i, const int &j)
    {
        assert(in_range(i, j));
        return array_p[i + rows * j];
    }

    const T &operator()(const int &i, const int &j) const
    {
        assert(in_range(i, j));
        return array_p[i + rows * j];
    }

    T get_field(const int &i, const int &j) const
    {
        if (not in_range(i, j)) {
            throw std::length_error("Cannot access indices [" + std::to_string(i) + "," + std::to_string(j)
                                    + "] with matrix dimensions [" + std::to_string(rows) + "," + std::to_string(cols)
                                    + "]");
            std::cout << "out_of range matrix coordinates: " << rows << "<" << i << " or " << cols << "<" << j
                      << std::endl;
        }
        return array_p[i + rows * j];
    }

    //Accessors
    const int &get_n_rows() const { return rows; }
    const int &get_n_cols() const { return cols; }

    // Raw storage access, for callers that need to compute a linear index once and share it
    // across several matrices with matching dimensions (see swalign::fill_sw_score_matrix).
    T *data() { return array_p; }
    const T *data() const { return array_p; }

    Matrix<T> transpose() const
    {
        Matrix<T> result(cols, rows);
        for (int i = 0; i != rows; ++i) {
            for (int j = 0; j != cols; ++j) {
                result(j, i) = array_p[i + rows * j];
            }
        }
        return result;
    }

    // Debug print
    void print(std::ostream &out = std::cout) const {
        out << rows << "x" << cols << " Matrix\n";
        if (rows == 0 || cols == 0) return;
        // Find max display width using string stream
        size_t max_width = 1;
        std::ostringstream oss;
        for (int k = 0; k < rows * cols; ++k) {
            oss.str("");
            oss.clear();
            oss << array_p[k];
            size_t len = oss.str().size();
            if (len > max_width) max_width = len;
        }
        // Print with fixed width
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (j > 0) out << " ";
                out << std::setw(static_cast<int>(max_width)) << std::right << array_p[i + j * rows];
            }
            out << "\n";
        }
    }

private:
    int rows;
    int cols;
    T *array_p;
};

template <typename T>
std::ostream &operator<<(std::ostream &stream, const Matrix<T> &mat)
{
    stream << mat.get_n_rows() << "x" << mat.get_n_cols() << " Matrix" << std::endl;
    for (int j = 0; j != mat.get_n_cols(); ++j) {
        for (int i = 0; i != mat.get_n_rows(); ++i) {
            if (i != 0) {
                stream << " ";
            }
            stream << mat.get_field(i, j);
        }
        stream << std::endl;
    }
    return stream;
}

} // namespace igor::math::legacy
