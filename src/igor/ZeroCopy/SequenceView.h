/*
 * SequenceView.h
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#pragma once

#include <igor/Core/Utils.h>
#include <igor/ZeroCopy/AlignmentView.h>
#include <igor/ZeroCopy/Export.h>

#include <sparrow/record_batch.hpp>
#include <string_view>

namespace igor {

/**
 * @brief Zero-copy view of a sequence record from an Arrow batch
 *
 * This class provides a lightweight view over sequence data stored in Arrow
 * columns. It does not own any memory and performs no allocations. It
 * implements the AIRR standard naming conventions.
 */
class ZEROCOPY_EXPORT SequenceView {
public:
    SequenceView(const sparrow::record_batch &batch, size_t row_index);

    // AIRR Standard Fields
    std::string_view get_sequence_id() const;  // sequence_id
    std::string_view get_sequence() const;     // sequence
    bool productive() const;                   // productive

    // Access to alignments (AIRR: v_call, d_call, j_call)
    AlignmentView get_alignment(Gene_class gene_type) const;

    // Helper to check if row is valid
    bool is_valid() const;

private:
    const sparrow::record_batch &batch_;
    size_t row_index_;
};

}  // namespace igor
