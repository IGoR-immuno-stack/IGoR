/*
 * Utils.h
 *
 *  Created on: Apr 9, 2015
 *      Author: Quentin Marcou
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 analyze and model immune receptors
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

#include <chrono>
#include <fstream>
#include <igor/Core/IntStr.h>
#include <iostream>
#include <list>
#include <memory>
#include <random>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

class Rec_Event;

enum Event_type { GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t, Undefined_t };
enum Event_safety { VD_safe = 0, DJ_safe = 1, VJ_safe = 2 };
enum Seq_side { Five_prime = 0, Three_prime = 1, Undefined_side = 2 };
enum Seq_type {
    V_gene_seq = 0,
    VD_ins_seq = 1,
    D_gene_seq = 2,
    DJ_ins_seq = 3,
    J_gene_seq = 4,
    VJ_ins_seq = 5
};
enum Gene_class {
    V_gene = 0,
    VD_genes = 1,
    D_gene = 2,
    DJ_genes = 3,
    J_gene = 4,
    VJ_genes = 5,
    VDJ_genes = 6,
    Undefined_gene = 7
};
enum Fileformat { CSV_f, FASTA_f, TXT_f, FASTQ_f };
enum Int_nt {
    int_A = 0,
    int_C = 1,
    int_G = 2,
    int_T = 3,
    int_R = 4,
    int_Y = 5,
    int_K = 6,
    int_M = 7,
    int_S = 8,
    int_W = 9,
    int_B = 10,
    int_D = 11,
    int_H = 12,
    int_V = 13,
    int_N = 14
};

Gene_class str2GeneClass(const std::string);
std::string to_string(const Gene_class);
Seq_side str2SeqSide(const std::string);
std::string to_string(const Seq_side);

std::ostream &operator<<(std::ostream &, Gene_class);
std::ostream &operator<<(std::ostream &, Seq_side);
std::string operator+(const std::string &, Gene_class);
std::string operator+(const std::string &, Seq_side);
std::string operator+(const std::string &, Event_type);

typedef std::unique_ptr<long double[]> Marginal_array_p;
typedef std::string Rec_Event_name;
typedef int Seq_Offset;
typedef Rec_Event *Next_event_ptr;

template <class T>
struct null_delete
{
    null_delete(void) = default;
    ~null_delete(void) = default;

    void operator()(T *) const { }
};

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

    T &operator()(const int &i, const int &j)
    {
        if ((i > rows - 1) || (j > cols - 1)) {
            throw std::length_error("Cannot access indices");
        }
        return array_p[i + rows * j];
    }

    const T &operator()(const int &i, const int &j) const
    {
        if ((i > rows - 1) || (j > cols - 1)) {
            throw std::length_error("Cannot access indices");
        }
        return array_p[i + rows * j];
    }

    T get_field(const int &i, const int &j) const
    {
        if ((i > rows - 1) || (j > cols - 1)) {
            throw std::length_error("Cannot access indices");
        }
        return array_p[i + rows * j];
    }

    const int &get_n_rows() const { return rows; }
    const int &get_n_cols() const { return cols; }

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

#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/FastMemoryMap.h>

typedef DynamicSequenceMap<Int_Str> Seq_type_str_p_map;
typedef Enum_fast_memory_map<int, bool> Safety_bool_map;
typedef Enum_fast_memory_map<int, std::vector<int> *> Mismatch_vectors_map;
typedef Enum_fast_memory_map<int, size_t> Index_map;
typedef Enum_fast_memory_map<int, double> Downstream_scenario_proba_bound_map;
typedef Enum_fast_memory_dual_key_map<int, Seq_side, Seq_Offset> Seq_offsets_map;

namespace std {
template <> struct hash<Seq_type> { std::size_t operator()(const Seq_type &seq_t) const { return hash<int>()(seq_t); } };
template <> struct hash<Gene_class> { std::size_t operator()(const Gene_class &gene) const { return hash<int>()(gene); } };
template <> struct hash<std::pair<Gene_class, Seq_side>> { std::size_t operator()(const pair<Gene_class, Seq_side> &gene_pair) const { return (hash<Gene_class>()(gene_pair.first) ^ (hash<int>()(gene_pair.second) << 1)) >> 1; } };
template <> struct hash<std::tuple<Event_type, Gene_class, Seq_side>> { std::size_t operator()(const std::tuple<Event_type, Gene_class, Seq_side> &event_triplet) const { Event_type ev_type; Gene_class g_class; Seq_side s_side; std::tie(ev_type, g_class, s_side) = event_triplet; return ((hash<int>()(ev_type) ^ (hash<int>()(g_class) << 1) >> 1) ^ (hash<int>()(s_side) << 1)); } };
template <> struct hash<std::tuple<Event_type, int, Seq_side>> { std::size_t operator()(const std::tuple<Event_type, int, Seq_side> &event_triplet) const { Event_type ev_type; int g_class; Seq_side s_side; std::tie(ev_type, g_class, s_side) = event_triplet; return ((hash<int>()(ev_type) ^ (hash<int>()(g_class) << 1) >> 1) ^ (hash<int>()(s_side) << 1)); } };
template <> struct hash<std::pair<Seq_type, Seq_side>> { std::size_t operator()(const std::pair<Seq_type, Seq_side> seq_pair) const { return (hash<int>()(seq_pair.first) ^ (hash<int>()(seq_pair.second) << 1)) >> 1; } };
template <> struct hash<Event_safety> { std::size_t operator()(const Event_safety ev_saf) const { return (hash<int>()(ev_saf)); } };
} // namespace std

struct D_position_comparator { bool operator()(std::tuple<std::string, int, int, double> position_1, std::tuple<std::string, int, int, double> position_2) { return std::get<3>(position_1) > std::get<3>(position_2); } };
struct inverse_offset_comparator { bool operator()(const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_1, const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_2) { return inv_offset_1.second < inv_offset_2.second; } };

std::vector<std::string> extract_string_fields(const std::string, const std::string);
void show_progress_bar(std::ostream &, double, std::string prefix_message = "", size_t progress_bar_size = 70);
void close_progress_bar(std::ostream &, std::string prefix_message = "", size_t progress_bar_size = 70);
uint64_t draw_random_64bits_seed();

typedef std::unordered_map<std::string, std::string> UMCodonTable;
std::string translate(std::string seq);
