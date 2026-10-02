/*
 * Utils.h
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
#include <cstdint>

#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <igor/Core/Legacy/IntStr.h>
#include <igor/Core/Legacy/CoreEnums.h>
#include <igor/Core/Legacy/StdTypedefs.h>
#include <igor/Core/Legacy/DynamicSequenceMap.h>
#include <igor/Core/Legacy/SeqOffsetsMap.h>
#include <igor/Core/Legacy/LayeredArray.h>
#include <memory>
#include <list>
#include <random>
#include <chrono>
#include <sys/types.h>
#include <igor/Core/Export.h>
#if defined(_WIN32)

#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif

#  include <process.h>
#  include <winsock2.h>
#  include <windows.h>
#  include <ws2tcpip.h>

inline int portable_getpid()
{
    return _getpid();
}

inline uint32_t portable_gethostid()
{
    // Version approximative: IP locale
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);

    char hostname[256];
    gethostname(hostname, sizeof(hostname));

    struct addrinfo hints{};
    hints.ai_family = AF_INET;

    struct addrinfo *info;
    if (getaddrinfo(hostname, nullptr, &hints, &info) != 0)
        return 0;

    uint32_t res = ((struct sockaddr_in *)info->ai_addr)->sin_addr.S_un.S_addr;

    freeaddrinfo(info);
    WSACleanup();
    return res;
}

#else

#  include <unistd.h>

inline int portable_getpid()
{
    return getpid();
}

inline uint32_t portable_gethostid()
{
    return gethostid();
}

#endif

#if defined(_MSC_VER)
#  include <intrin.h>
#endif

// Cross-platform population count (number of set bits).
// Used by the genetic-code utilities (GeneticCode.h) to size CodonMask sets.
#include <stdio.h>
#include <unordered_map>

namespace igor::core::legacy {

// The legacy to_string(Gene_class) overloads below would otherwise hide std::to_string
// from every legacy namespace that nominates this one.
using std::to_string;

inline int popcountll(uint64_t x)
{
#if defined(_MSC_VER)
    return __popcnt64(x);
#else
    return __builtin_popcountll(x);
#endif
}

// Force a function to be inlined at every call site, even at optimization levels where the
// compiler's own heuristics would otherwise leave it as a real call (e.g. -O2 on a function
// with several branches). Use sparingly, only where profiling has shown the call overhead
// itself (prologue/epilogue, stack-protector checks) to be a measurable cost.
#if defined(_MSC_VER)
#  define IGOR_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#  define IGOR_ALWAYS_INLINE __attribute__((__always_inline__)) inline
#else
#  define IGOR_ALWAYS_INLINE inline
#endif

// Force full unroll of a small, fixed-trip-count loop immediately below the pragma, so its
// per-iteration work stays visible to the scheduler as independent instructions instead of being
// serialized behind the loop's own increment/compare/branch (see fill_sw_score_matrix's
// column-banding note in Aligner.cpp for a case where that overhead was measurable). n must be an
// integer literal or an object-like macro -- it is textually stringized into a pragma at
// preprocessing time, so a constexpr variable will not work.
// Clang also parses "GCC unroll", including under clang-cl (which defines both __clang__ and
// _MSC_VER); check __clang__ before _MSC_VER so clang-cl takes this branch, not the MSVC one.
#define IGOR_UNROLL_STRINGIFY_(x) #x
#if defined(__clang__) || defined(__GNUC__)
#  define IGOR_UNROLL(n) _Pragma(IGOR_UNROLL_STRINGIFY_(GCC unroll n))
#elif defined(__INTEL_COMPILER)
// Classic (pre-oneAPI) Intel compiler; icx/dpcpp is LLVM-based and already covered by __clang__.
#  define IGOR_UNROLL(n) _Pragma(IGOR_UNROLL_STRINGIFY_(unroll(n)))
#else
// No portable equivalent on MSVC (or other unrecognized compilers): the loop is left to the
// compiler's own unrolling heuristics, or must be unrolled by hand for a guaranteed effect.
#  define IGOR_UNROLL(n)
#endif


/// Slim gene class: only the three fundamental gene types plus Undefined.
/// Use this everywhere in runtime code. Junction/compound values (VD, DJ, VJ, VDJ)
/// belong to Gene_class_legacy and are only used at legacy file I/O boundaries.
enum CORE_EXPORT Gene_class {
    V_gene = 0,
    D_gene = 1,
    J_gene = 2,
    Undefined_gene = 3
};

CORE_EXPORT Gene_class str2GeneClassNew(const std::string &);
CORE_EXPORT std::string to_string(const Gene_class);
CORE_EXPORT std::ostream &operator<<(std::ostream &, Gene_class);
CORE_EXPORT std::string operator+(const std::string &, Gene_class);

enum CORE_EXPORT Gene_class_legacy {
    V_gene_legacy = 0,
    VD_genes = 1,   ///< LEGACY: junction class, kept for legacy file read/write only
    D_gene_legacy = 2,
    DJ_genes = 3,   ///< LEGACY: junction class, kept for legacy file read/write only
    J_gene_legacy = 4,
    VJ_genes = 5,   ///< LEGACY: junction class, kept for legacy file read/write only
    VDJ_genes = 6,  ///< LEGACY: junction class, kept for legacy file read/write only
    Undefined_gene_legacy = 7
};

/// Convert legacy gene class (base gene only, not junction) to new Gene_class.
CORE_EXPORT Gene_class gene_class_legacy_to_new(Gene_class_legacy);
enum Fileformat { CSV_f, FASTA_f, TXT_f, FASTQ_f };
/**
 * \brief The IUPAC nucleotide codes, plus one value that is not a nucleotide.
 *
 * `int_A` through `int_N` are the fifteen codes a read can contain: four bases and the eleven
 * ambiguity codes. `int_undefined` is a **third state** on the same axis, and the distinction
 * it makes is easy to lose:
 *
 * - `int_N` means *this position is determined, and the read does not say which base it is*.
 *   It comes from the data, and every consumer handles it -- averaged over its underlying
 *   bases by Dinucl_markov, matched permissively by the aligner.
 * - `int_undefined` means *this position is not determined yet*. It comes from
 *   Dinucl_markov::iterate, which creates a junction of the right length from the offsets
 *   the Insertion placed, before anything knows its content, and then fills every position.
 *   **No read ever contains it**: nt2int() cannot produce it.
 *
 * It is a placeholder inside one scenario, not a value. A segment handed to an error rate, a
 * counter or an output file must contain none.
 *
 * It sits immediately *past* the real codes rather than at -1, so that using it as an index
 * is caught rather than silently wrapping: `dinuc_proba_matrix` is `kIntNtCount` square, so
 * `matrix(int_undefined, j)` trips its bounds assertion, while `matrix(-1, j)` would have read
 * one row before the array. The same choice keeps `< 4` tests meaning what they read as: an
 * undefined position is *not* one of the four bases, where -1 satisfied `x < 4` and sent an
 * undefined nucleotide down the unambiguous path with a negative marginal offset.
 */
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
    int_N = 14,
    int_undefined = 15 ///< allocated but not yet filled; never present in a read
};

/// Number of real nucleotide codes, i.e. every value a read can hold. Equal to
/// `int_undefined` by construction: anything that sizes a table by nucleotide gets a table
/// the placeholder cannot index into.
constexpr std::size_t kIntNtCount = static_cast<std::size_t>(int_undefined);

CORE_EXPORT Seq_type str2SeqType(const Seq_type_String &);
CORE_EXPORT Seq_type_String to_string(const Seq_type);
CORE_EXPORT Gene_class_legacy str2GeneClass(const std::string &);
CORE_EXPORT std::string to_string(const Gene_class_legacy);
CORE_EXPORT Seq_side str2SeqSide(const std::string &);
CORE_EXPORT std::string to_string(const Seq_side);

CORE_EXPORT std::ostream &operator<<(std::ostream &, Gene_class_legacy);
CORE_EXPORT std::ostream &operator<<(std::ostream &, Seq_side);
CORE_EXPORT std::string operator+(const std::string &, Gene_class_legacy);
CORE_EXPORT std::string operator+(const std::string &, Seq_side);
CORE_EXPORT std::string operator+(const std::string &, Event_type);

typedef Int_Str *Int_Str_ptr;

/**
 * \brief Declare a null_delete function
 * \author Q.Marcou
 * This function is not performing any task, it's purpose is to supply a "null_delete" function
 * to prevent shared pointer objects created when passing Rec_Event or Error_rate objects pointers to model_parms
 * to be destroyed when the model_parms object is destroyed itself and the rec_event and error_rate objects
 * might still be of used (and if not prevent from a segfault error by trying to delete them twice)
 */
template <class T>
struct null_delete
{
    null_delete(void) = default;
    ~null_delete(void) = default;

    void operator()(T *) const { }
};

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


//Keyed by SeqTypeId. This is the one map whose values can be *actively absent*: a
//zero-length segment means an event ran and produced nothing, which the ordered traversal
//distinguishes from a segment that has not been written yet (see SeqSegmentEmptiness).
typedef DynamicSequenceMap<Int_Str_ptr> Seq_type_str_p_map;

/**
 * The first segment still holding an `int_undefined` position, or `kNoSeqType` if none does.
 *
 * `int_undefined` is a placeholder inside one scenario, not a value: a `Dinucl_markov`
 * creates a junction of the right length holding it, and fills it before handing on. By the
 * time a scenario is complete every position must be determined, so anything that *consumes*
 * a finished scenario -- the error rate, the counters, the output writers -- may assume it.
 *
 * The assumption is checked at the one place it has to hold (`Rec_Event::iterate_wrap_up`'s
 * leaf branch) and only in a build with assertions enabled: the walk is linear in the
 * scenario's length and the default build is RelWithDebInfo, which defines `NDEBUG`. Exposed
 * rather than hidden in that assert so it can be unit-tested without a debug build.
 */
CORE_EXPORT SeqTypeId first_unfilled_segment(const Seq_type_str_p_map &constructed_sequences);

/**
 * The first segment end with no offset written, or `{kNoSeqType, Undefined_side}` if every
 * end is placed.
 *
 * The offsets counterpart of first_unfilled_segment(), and the other half of one invariant:
 *
 * > For every seq_type in the model, both a sequence **and** its offsets must have been
 * > created by the time a scenario reaches a leaf -- not necessarily by the same event.
 *
 * It could not be written before R3. While `Insertion` recorded no offsets at all, this would
 * have fired on every scenario that has a junction, which is every scenario; the check exists
 * because that stopped being true (plan section 2.5).
 *
 * The key is per (`SeqTypeId`, `Seq_side`) rather than per `SeqTypeId`, because
 * `get_offset_role` is side-taking and the two ends of a segment can in principle be created
 * by different events.
 *
 * **Swept over `registry.ordering()`, not over every registered id.** The registry pins all six
 * legacy names whatever the model is, so a VJ model carries `D_gene_seq` and both flanking
 * junctions as ids that no event in it will ever place. The ordering is the model's actual
 * segment layout, which is the set the invariant is about; unlike the content half this
 * predicate cannot skip what is absent, since absence is the thing it looks for.
 *
 * Checked at the same boundary and under the same conditions as the content half -- the leaf
 * branch of `Rec_Event::iterate_wrap_up`, assertions only. Exposed rather than hidden in that
 * assert so it can be unit-tested without a debug build.
 */
CORE_EXPORT std::pair<SeqTypeId, Seq_side> first_unplaced_segment_end(
        const Seq_offsets_map &seq_offsets);

//Keyed by SeqTypeId. An empty mismatch vector means "zero mismatches", which is present,
//not absent -- so the default SeqSegmentEmptiness (never empty) is the right one here.
typedef DynamicSequenceMap<std::vector<size_t> *> Mismatch_vectors_map;

/// NT-floor mismatch positions per Seq_type, used as a conservative pruning bound.
/// Same underlying type as Mismatch_vectors_map; the two carry different semantics:
///   Mismatch_vectors_map  - upper bound: position mismatches in at least one branch
///   Pruning_mismatch_floor_map - floor: position mismatches in every branch
/// For exact NT queries the two tracks are identical.
typedef DynamicSequenceMap<std::vector<size_t> *> Pruning_mismatch_floor_map;

//Index_map is keyed by event identifier, not by seq_type, so it uses the bare
//layered container rather than the registry-aware DynamicSequenceMap (plan D5).
typedef LayeredArray<size_t> Index_map;

//Keyed by SeqTypeId; the legacy Seq_type enum values are pinned to the same ids by
//Model_Parms::read_model_parms(), so enum-keyed call sites still address the right slot.
typedef DynamicSequenceMap<double> Downstream_scenario_proba_bound_map;


//Seq_offsets_map is a class in SeqOffsetsMap.h, holding one DynamicSequenceMap per end.


/*
 * Hash functions for the enums and tuples used as unordered_map keys
 */
} // namespace igor::core::legacy
namespace std {
template <>
struct hash<igor::core::legacy::Seq_type>
{
    std::size_t operator()(const igor::core::legacy::Seq_type &seq_t) const { return hash<int>()(seq_t); }
};

template <>
struct hash<igor::core::legacy::Gene_class>
{
    std::size_t operator()(const igor::core::legacy::Gene_class &gene) const { return hash<int>()(gene); }
};

template <>
struct hash<igor::core::legacy::Gene_class_legacy>
{
    std::size_t operator()(const igor::core::legacy::Gene_class_legacy &gene) const { return hash<int>()(gene); }
};

template <>
struct hash<std::tuple<igor::core::legacy::Event_type, igor::core::legacy::Seq_type_String, igor::core::legacy::Seq_side>>
{
    std::size_t operator()(const std::tuple<igor::core::legacy::Event_type, igor::core::legacy::Seq_type_String, igor::core::legacy::Seq_side> &event_triplet) const
    {
        igor::core::legacy::Event_type ev_type;
        igor::core::legacy::Seq_type_String seq_type_str;
        igor::core::legacy::Seq_side s_side;
        std::tie(ev_type, seq_type_str, s_side) = event_triplet;
        return ((hash<int>()(ev_type) ^ (hash<igor::core::legacy::Seq_type_String>()(seq_type_str) << 1) >> 1)
                ^ (hash<int>()(s_side) << 1));
    }
};

template <>
struct hash<std::pair<igor::core::legacy::Seq_type, igor::core::legacy::Seq_side>>
{
    std::size_t operator()(const std::pair<igor::core::legacy::Seq_type, igor::core::legacy::Seq_side> seq_pair) const
    {
        return (hash<int>()(seq_pair.first) ^ (hash<int>()(seq_pair.second) << 1)) >> 1;
    }
};
} // namespace std
namespace igor::core::legacy {


CORE_EXPORT std::vector<std::string> extract_string_fields(const std::string &, const std::string &);

CORE_EXPORT void show_progress_bar(std::ostream &, double, const std::string &prefix_message = "", size_t progress_bar_size = 70);
CORE_EXPORT void close_progress_bar(std::ostream &, const std::string &prefix_message = "", size_t progress_bar_size = 70);
CORE_EXPORT uint64_t draw_random_64bits_seed();


CORE_EXPORT std::string translate(const std::string &seq);

} // namespace igor::core::legacy
