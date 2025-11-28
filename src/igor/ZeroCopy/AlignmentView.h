/*
 * AlignmentView.h
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#pragma once

#include <igor/Core/Utils.h>
#include <igor/ZeroCopy/Export.h>

#include <optional>
#include <sparrow/record_batch.hpp>
#include <string_view>

namespace igor {

/**
 * @brief Zero-copy view of an alignment record from an Arrow batch
 *
 * This class provides a lightweight view over alignment data stored in Arrow
 * columns. It does not own any memory and performs no allocations. It
 * implements the AIRR standard naming conventions.
 */
class ZEROCOPY_EXPORT AlignmentView {
public:
    AlignmentView(const sparrow::record_batch &batch, size_t row_index, Gene_class gene_type);

    // AIRR Standard Fields
    std::string_view get_gene_call() const;  // v_call, d_call, j_call
    double get_score() const;                // v_score, d_score, j_score

    // Alignment coordinates (1-based, inclusive per AIRR)
    int get_sequence_start() const;  // v_sequence_start
    int get_sequence_end() const;    // v_sequence_end
    int get_germline_start() const;  // v_germline_start
    int get_germline_end() const;    // v_germline_end

    // CIGAR string for gaps/mismatches
    std::string_view get_cigar() const;  // v_cigar

    // Legacy IGoR compatibility (computed on the fly if possible, or mapped)
    int get_offset() const;           // Mapped from start/end
    size_t get_align_length() const;  // Mapped from cigar

    bool is_valid() const;

private:
    const sparrow::record_batch &batch_;
    size_t row_index_;
    Gene_class gene_type_;

    // Helper to resolve column names based on gene type (e.g., V_gene ->
    // "v_call")
    std::string get_column_name(const std::string &suffix) const;
};

}  // namespace igor
