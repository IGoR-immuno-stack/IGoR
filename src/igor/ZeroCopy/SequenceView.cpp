/*
 * SequenceView.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#include "SequenceView.h"

#include <algorithm>  // Required for std::find
#include <stdexcept>
#include <variant>  // Required for std::visit

namespace igor {

SequenceView::SequenceView(const sparrow::record_batch &batch, size_t row_index)
    : batch_(batch), row_index_(row_index) {}

std::string_view SequenceView::get_sequence_id() const {
    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), "sequence_id") == names.end()) {
            return "";
        }
        const auto &col = batch_.get_column("sequence_id");
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

std::string_view SequenceView::get_sequence() const {
    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), "sequence") == names.end()) {
            return "";
        }
        const auto &col = batch_.get_column("sequence");
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

bool SequenceView::productive() const {
    try {
        const auto &names = batch_.names();
        if (std::find(names.begin(), names.end(), "productive") == names.end()) {
            return false;
        }
        const auto &col = batch_.get_column("productive");
        auto value = col[row_index_];
        return std::visit(
            [](auto &&arg) -> bool {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, sparrow::nullable<bool>>) {
                    if (arg.has_value()) {
                        return arg.get();
                    }
                }
                return false;
            },
            value);
    } catch (...) {
        return false;
    }
}

AlignmentView SequenceView::get_alignment(Gene_class gene_type) const {
    return AlignmentView(batch_, row_index_, gene_type);
}

bool SequenceView::is_valid() const { return row_index_ < batch_.nb_rows(); }

}  // namespace igor
