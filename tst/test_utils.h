/*
 * test_utils.h
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

#pragma once

#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <memory>
#include <string>
#include <vector>
#include <random>

namespace IgorTestUtils {

/**
 * @brief Create mock alignment data for testing
 * 
 * @param gene_name Name of the gene (e.g., "TRBV1", "TRBJ1-1")
 * @param offset Alignment position on target sequence
 * @param five_p_offset 5' alignment boundary
 * @param three_p_offset 3' alignment boundary
 * @param mismatches Vector of mismatch positions
 * @param score Alignment score
 * @return Alignment_data Mock alignment data structure
 */
Alignment_data create_mock_alignment_data(
    const std::string& gene_name,
    int offset,
    size_t five_p_offset,
    size_t three_p_offset,
    const std::vector<int>& mismatches = {},
    double score = 100.0
);

} // namespace IgorTestUtils
