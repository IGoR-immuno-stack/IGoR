/*
 * test_utils.cpp
 *
 *  Created on: Jan 21, 2026
 *      Author: IGoR Test Suite
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
 */

#include "test_utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <functional>

namespace IgorTestUtils {

Alignment_data create_mock_alignment_data(
    const std::string& gene_name,
    int offset,
    size_t five_p_offset,
    size_t three_p_offset,
    const std::vector<int>& mismatches,
    double score
) {
    // Calculate alignment length
    size_t align_length = three_p_offset - five_p_offset;
    
    // Create empty forward_lists for insertions and deletions
    std::forward_list<int> empty_insertions;
    std::forward_list<int> empty_deletions;
    
    // Use the appropriate constructor
    Alignment_data align_data(
        gene_name,
        offset,
        five_p_offset,
        three_p_offset,
        align_length,
        empty_insertions,
        empty_deletions,
        mismatches,
        score
    );
    
    return align_data;
}

} // namespace IgorTestUtils
