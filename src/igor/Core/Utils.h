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
#include <igor/Core/IntStr.h>
#include <igor/Core/CoreEnums.h>
#include <igor/Core/StdTypedefs.h>
#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/SeqOffsetsMap.h>
#include <igor/Core/LayeredArray.h>
#include <igor/Core/LayeredMap.h>
#include <memory>
#include <list>
#include <random>
#include <chrono>
#include <sys/types.h>
#include <igorCoreExport.h>
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

// Marks a symbol that is internal to Core (not part of its installed public API) but still
// needs to cross the shared library boundary for whitebox tests to link against it directly
// (see AlignerInternal.h). Resolves to a real export/import only when CORE_TESTING_ENABLED is
// defined (see tst/igor/Core/CMakeLists.txt), so production builds keep these symbols hidden.
#ifdef CORE_TESTING_ENABLED
#  define CORE_TESTING_EXPORT CORE_EXPORT
#else
#  define CORE_TESTING_EXPORT CORE_NO_EXPORT
#endif

#include <stdio.h>
#include <unordered_map>

class Rec_Event;

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

//Typedef used for getting the next event ptr
//typedef std::shared_ptr<Rec_Event> Next_event_ptr; //Does not work for some reason
typedef Rec_Event *Next_event_ptr;

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

    T &operator()(const int &i, const int &j)
    {
        assert((i <= rows - 1) && (j <= cols - 1));
        return array_p[i + rows * j];
    }

    const T &operator()(const int &i, const int &j) const
    {
        assert((i <= rows - 1) && (j <= cols - 1));
        return array_p[i + rows * j];
    }

    T get_field(const int &i, const int &j) const
    {
        if ((i > rows - 1) || (j > cols - 1)) {
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

//Keyed by the Event_safety enum as an opaque dense integer, hence LayeredArray rather than
//DynamicSequenceMap: safety is not per-seq_type today. It is *pairwise* between gene
//segments -- VD_safe, DJ_safe and VJ_safe are the three unordered pairs of V, D and J, and
//VJ_safe is live in VDJ models where V and J are not adjacent. Re-expressing it per junction
//only becomes correct once B5/B6/B11 find neighbours dynamically; see the note under B5.
typedef LayeredArray<bool> Safety_bool_map;

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

/*
	template<> class Enum_fast_memory_map<Seq_type ,Str_ptr>{
	Enum_fast_memory_map():Enum_fast_memory_map<Seq_type,Str_ptr>(6){};
	size_t range = 6;
};
*/

/*template<>
class Seq_type_str_p_map : public Enum_fast_memory_map<Seq_type,Str_ptr>{

};*/


//Seq_offsets_map is a class in SeqOffsetsMap.h, holding one DynamicSequenceMap per end.

/**
 * String-keyed variant of Enum_fast_memory_dual_key_map.
 * First key is Seq_type_String (arbitrary seq_type name such as "D1_gene_seq"),
 * second key K2 is any type castable to int (e.g. Seq_side enum).
 * Supports the same layered-memory interface as Enum_fast_memory_dual_key_map.
 */
template <typename K2, typename V>
class Str_Dual_key_memory_map {
public:
    Str_Dual_key_memory_map() = default;

    V &at(const Seq_type_String &key1, const K2 &key2)
    {
        auto &entry = get_entry(key1, key2);
        if (entry.current_layer < 0)
            throw std::out_of_range("Trying to access uninitialized position in Str_Dual_key_memory_map::at()");
        return entry.layers[entry.current_layer];
    }

    const V &at(const Seq_type_String &key1, const K2 &key2) const
    {
        const auto &entry = get_entry_const(key1, key2);
        if (entry.current_layer < 0)
            throw std::out_of_range("Trying to access uninitialized position in Str_Dual_key_memory_map::at()");
        return entry.layers[entry.current_layer];
    }

    // Access (and rewind to) a specific memory layer.
    V &at(const Seq_type_String &key1, const K2 &key2, int memory_layer)
    {
        auto &entry = get_entry(key1, key2);
        if (memory_layer > entry.current_layer + 1)
            throw std::out_of_range("Trying to access uninitialized position in Str_Dual_key_memory_map::at()");
        if (memory_layer >= static_cast<int>(entry.layers.size()))
            entry.layers.resize(memory_layer + 1);
        entry.current_layer = memory_layer;
        return entry.layers[memory_layer];
    }

    const V &at(const Seq_type_String &key1, const K2 &key2, int memory_layer) const
    {
        const auto &entry = get_entry_const(key1, key2);
        if (memory_layer > entry.current_layer)
            throw std::out_of_range("Trying to access uninitialized position in Str_Dual_key_memory_map::at()");
        return entry.layers[memory_layer];
    }

    int get_current_memory_layer(const Seq_type_String &key1, const K2 &key2)
    {
        auto it1 = data_.find(key1);
        if (it1 == data_.end()) return -1;
        auto it2 = it1->second.find(static_cast<int>(key2));
        if (it2 == it1->second.end()) return -1;
        return it2->second.current_layer;
    }

    bool exist(const Seq_type_String &key1, const K2 &key2)
    {
        auto it1 = data_.find(key1);
        if (it1 == data_.end()) return false;
        auto it2 = it1->second.find(static_cast<int>(key2));
        if (it2 == it1->second.end()) return false;
        return it2->second.current_layer >= 0;
    }

    void request_memory_layer(const Seq_type_String &key1, const K2 &key2)
    {
        auto &entry = data_[key1][static_cast<int>(key2)];
        ++entry.current_layer;
        if (entry.current_layer >= static_cast<int>(entry.layers.size()))
            entry.layers.emplace_back();
    }

    void set_value(const Seq_type_String &key1, const K2 &key2, V value, int memory_layer)
    {
        auto &entry = data_[key1][static_cast<int>(key2)];
        if (memory_layer > entry.current_layer + 1)
            throw std::out_of_range(
                "Trying to access incorrect memory layer in Str_Dual_key_memory_map::set_value()");
        if (memory_layer >= static_cast<int>(entry.layers.size()))
            entry.layers.resize(memory_layer + 1);
        entry.layers[memory_layer] = value;
        entry.current_layer = memory_layer;
    }

private:
    struct LayeredValue {
        std::vector<V> layers;
        int current_layer = -1;
    };
    std::unordered_map<Seq_type_String, std::unordered_map<int, LayeredValue>> data_;

    LayeredValue &get_entry(const Seq_type_String &key1, const K2 &key2)
    {
        auto it1 = data_.find(key1);
        if (it1 == data_.end())
            throw std::out_of_range("Unknown key1 in Str_Dual_key_memory_map::at(): " + key1);
        auto it2 = it1->second.find(static_cast<int>(key2));
        if (it2 == it1->second.end())
            throw std::out_of_range("Unknown key2 in Str_Dual_key_memory_map::at()");
        return it2->second;
    }

    const LayeredValue &get_entry_const(const Seq_type_String &key1, const K2 &key2) const
    {
        auto it1 = data_.find(key1);
        if (it1 == data_.end())
            throw std::out_of_range("Unknown key1 in Str_Dual_key_memory_map::at(): " + key1);
        auto it2 = it1->second.find(static_cast<int>(key2));
        if (it2 == it1->second.end())
            throw std::out_of_range("Unknown key2 in Str_Dual_key_memory_map::at()");
        return it2->second;
    }
};

typedef Str_Dual_key_memory_map<Seq_side, Seq_Offset> Str_Seq_offsets_map;

/*
 * Defining a hash functions for Rec_Event, Gene_class_legacy and pair<Gene_class_legacy,Seq_side>
 */
namespace std {
/*
 	 template<>
	 struct hash<Rec_Event>{
		inline std::size_t operator()(const Rec_Event& event) const{ //TODO inline?
			return  (((hash<int>()(event.get_class())
					^(hash<int>()(event.get_side())<<1 )) >>1)
					^(hash<int>()(event.get_priority())<<1)>>1)
					^(hash<int>()(event.get_realizations_map().size())<<1);
			//Note : only consider the size of the realization map and not what it contains for speed purposes
			//this should be enough to ensure no collisions
		}
	 };
	 */

/*
	 template<>
	 struct hash<Rec_Event*>{
		 std::size_t operator()(const Rec_Event*& event_point) const{
			 return hash<Rec_Event>()(*event_point);
		 }
	 };
	 */

template <>
struct hash<Seq_type>
{
    std::size_t operator()(const Seq_type &seq_t) const { return hash<int>()(seq_t); }
};

template <>
struct hash<Gene_class>
{
    std::size_t operator()(const Gene_class &gene) const { return hash<int>()(gene); }
};

template <>
struct hash<Gene_class_legacy>
{
    std::size_t operator()(const Gene_class_legacy &gene) const { return hash<int>()(gene); }
};

template <>
struct hash<std::pair<Gene_class_legacy, Seq_side>>
{
    std::size_t operator()(const pair<Gene_class_legacy, Seq_side> &gene_pair) const
    {
        return (hash<Gene_class_legacy>()(gene_pair.first) ^ (hash<int>()(gene_pair.second) << 1)) >> 1;
    }
};

template <>
struct hash<std::tuple<Event_type, Gene_class_legacy, Seq_side>>
{
    std::size_t operator()(const std::tuple<Event_type, Gene_class_legacy, Seq_side> &event_triplet) const
    {
        Event_type ev_type;
        Gene_class_legacy g_class;
        Seq_side s_side;
        std::tie(ev_type, g_class, s_side) = event_triplet;
        return ((hash<int>()(ev_type) ^ (hash<int>()(g_class) << 1) >> 1) ^ (hash<int>()(s_side) << 1));
    }
};

template <>
struct hash<std::tuple<Event_type, Seq_type_String, Seq_side>>
{
    std::size_t operator()(const std::tuple<Event_type, Seq_type_String, Seq_side> &event_triplet) const
    {
        Event_type ev_type;
        Seq_type_String seq_type_str;
        Seq_side s_side;
        std::tie(ev_type, seq_type_str, s_side) = event_triplet;
        return ((hash<int>()(ev_type) ^ (hash<Seq_type_String>()(seq_type_str) << 1) >> 1)
                ^ (hash<int>()(s_side) << 1));
    }
};

template <>
struct hash<std::pair<Seq_type, Seq_side>>
{
    std::size_t operator()(const std::pair<Seq_type, Seq_side> seq_pair) const
    {
        return (hash<int>()(seq_pair.first) ^ (hash<int>()(seq_pair.second) << 1)) >> 1;
    }
};

template <>
struct hash<Event_safety>
{
    std::size_t operator()(const Event_safety ev_saf) const { return (hash<int>()(ev_saf)); }
};
} // namespace std

// v2.0 events map: keyed by (Event_type, seq_type string, Seq_side) so that
// multiple events of the same type and side but different seq_types (e.g. two
// D-gene deletions on D1 vs D2 in a tandem-D model) can coexist unambiguously.
typedef std::unordered_map<std::tuple<Event_type, Seq_type_String, Seq_side>,
                           std::shared_ptr<Rec_Event>>
        Events_map;

struct D_position_comparator
{
    bool operator()(std::tuple<std::string, int, int, double> position_1,
                    std::tuple<std::string, int, int, double> position_2)
    {
        return std::get<3>(position_1) > std::get<3>(position_2);
    }
};

struct inverse_offset_comparator
{
    bool operator()(const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_1,
                    const std::pair<std::shared_ptr<const Rec_Event>, int> &inv_offset_2)
    {
        return inv_offset_1.second < inv_offset_2.second;
    }
};

std::vector<std::string> extract_string_fields(const std::string &, const std::string &);

void show_progress_bar(std::ostream &, double, const std::string &prefix_message = "", size_t progress_bar_size = 70);
void close_progress_bar(std::ostream &, const std::string &prefix_message = "", size_t progress_bar_size = 70);
uint64_t draw_random_64bits_seed();


std::string translate(const std::string &seq);
