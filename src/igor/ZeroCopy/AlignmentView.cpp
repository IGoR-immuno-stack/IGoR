/*
 * AlignmentView.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#include "AlignmentView.h"

#include <algorithm>  // Required for std::find
#include <stdexcept>
#include <variant>  // Required for std::visit

namespace igor {

AlignmentView::AlignmentView(const sparrow::record_batch &batch, size_t row_index,
                             Gene_class gene_type)
    : batch_(batch), row_index_(row_index), gene_type_(gene_type) {}

std::string AlignmentView::get_column_name(const std::string &suffix) const {
    std::string prefix;
    switch (gene_type_) {
        case V_gene:
            prefix = "v_";
            break;
        case D_gene:
            prefix = "d_";
            break;
        case J_gene:
            prefix = "j_";
            break;
        default:
            return "";
    }
    return prefix + suffix;
}

std::string_view AlignmentView::get_gene_call() const {
    std::string col_name = get_column_name("call");
    if (col_name.empty()) return "";

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return "";
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> std::string_view {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<std::string_view>>) {
                    if (arg.has_value()) {
                        return arg.get();
                    }
                }
                return "";
            },
            value);
    } catch (...) {
        return "";
    }
}

double AlignmentView::get_score() const {
    std::string col_name = get_column_name("score");
    if (col_name.empty()) return 0.0;

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return 0.0;
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> double {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<const double &>> ||
                              std::is_same_v<T, sparrow::nullable<const float &>>) {
                    if (arg.has_value()) {
                        return static_cast<double>(arg.get());
                    }
                }
                return 0.0;
            },
            value);
    } catch (...) {
        return 0.0;
    }
}

int AlignmentView::get_sequence_start() const {
    std::string col_name = get_column_name("sequence_start");
    if (col_name.empty()) return 0;

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return 0;
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> int {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<const int64_t &>> ||
                              std::is_same_v<T, sparrow::nullable<const int32_t &>> ||
                              std::is_same_v<T, sparrow::nullable<int64_t>> ||
                              std::is_same_v<T, sparrow::nullable<int32_t>>) {
                    if (arg.has_value()) {
                        return static_cast<int>(arg.get());
                    }
                }
                return 0;
            },
            value);
    } catch (...) {
        return 0;
    }
}

int AlignmentView::get_sequence_end() const {
    std::string col_name = get_column_name("sequence_end");
    if (col_name.empty()) return 0;

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return 0;
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> int {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<const int64_t &>> ||
                              std::is_same_v<T, sparrow::nullable<const int32_t &>> ||
                              std::is_same_v<T, sparrow::nullable<int64_t>> ||
                              std::is_same_v<T, sparrow::nullable<int32_t>>) {
                    if (arg.has_value()) {
                        return static_cast<int>(arg.get());
                    }
                }
                return 0;
            },
            value);
    } catch (...) {
        return 0;
    }
}

int AlignmentView::get_germline_start() const {
    std::string col_name = get_column_name("germline_start");
    if (col_name.empty()) return 0;

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return 0;
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> int {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<const int64_t &>> ||
                              std::is_same_v<T, sparrow::nullable<const int32_t &>> ||
                              std::is_same_v<T, sparrow::nullable<int64_t>> ||
                              std::is_same_v<T, sparrow::nullable<int32_t>>) {
                    if (arg.has_value()) {
                        return static_cast<int>(arg.get());
                    }
                }
                return 0;
            },
            value);
    } catch (...) {
        return 0;
    }
}

int AlignmentView::get_germline_end() const {
    std::string col_name = get_column_name("germline_end");
    if (col_name.empty()) return 0;

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return 0;
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> int {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<const int64_t &>> ||
                              std::is_same_v<T, sparrow::nullable<const int32_t &>> ||
                              std::is_same_v<T, sparrow::nullable<int64_t>> ||
                              std::is_same_v<T, sparrow::nullable<int32_t>>) {
                    if (arg.has_value()) {
                        return static_cast<int>(arg.get());
                    }
                }
                return 0;
            },
            value);
    } catch (...) {
        return 0;
    }
}

std::string_view AlignmentView::get_cigar() const {
    std::string col_name = get_column_name("cigar");
    if (col_name.empty()) return "";

    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), col_name) == names.end()) {
            return "";
        }
        const auto &col = batch_.get_column(col_name);
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> std::string_view {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<std::string_view>>) {
                    if (arg.has_value()) {
                        return arg.get();
                    }
                }
                return "";
            },
            value);
    } catch (...) {
        return "";
    }
}

int AlignmentView::get_offset() const {
    // Legacy mapping logic: offset is 0-based index of start
    // AIRR sequence_start is 1-based
    int start = get_sequence_start();
    return start > 0 ? start - 1 : 0;
}

size_t AlignmentView::get_align_length() const {
    // Simplified length calculation from start/end
    // Real implementation should parse CIGAR string
    int start = get_sequence_start();
    int end = get_sequence_end();
    if (start > 0 && end >= start) {
        return static_cast<size_t>(end - start + 1);
    }
    return 0;
}

bool AlignmentView::is_valid() const { return !get_gene_call().empty(); }

}  // namespace igor
