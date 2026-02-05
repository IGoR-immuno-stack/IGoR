/*
 * ModelReader.h
 *
 *  Created on: Feb 5, 2026
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 *  analyze and model immune receptors generation, selection, mutation and all
 *  other processes.
 *   Copyright (C) 2026  IGoR Development Team
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
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file ModelReader.h
 * @brief Reader functions for IGoR model files
 *
 * This module provides functionality to read IGoR model files in both
 * the modern JSON format and legacy text format.
 *
 * **Supported Formats:**
 * - JSON format (model_struct.json, model_weights.json)
 * - Legacy format (model_parms.txt, model_marginals.txt)
 *
 * **Usage Example:**
 * @code
 *   using namespace igor::modelio;
 *
 *   // Read JSON format
 *   auto model = structure::read_json("model_struct.json");
 *   auto weights = weights::read_json("model_weights.json");
 *
 *   // Read legacy format
 *   auto legacy_model = structure::read_legacy("model_parms.txt");
 *   auto legacy_weights = weights::read_legacy("model_marginals.txt");
 * @endcode
 */

#pragma once

#include <igor/ModelIO/ModelIOCommon.h>

#include <string>

namespace igor::modelio::structure {

/**
 * @brief Read model structure from JSON file
 *
 * @param filepath Path to model_struct.json file
 * @return ModelData containing events, edges, and sequence types
 * @throws std::runtime_error if file cannot be read or parsed
 */
MODELIO_EXPORT
ModelData read_json(const std::string& filepath);

/**
 * @brief Read model structure from legacy model_parms.txt file
 *
 * @param filepath Path to model_parms.txt file
 * @return ModelData containing events and edges
 * @throws std::runtime_error if file cannot be read or parsed
 */
MODELIO_EXPORT
ModelData read_legacy(const std::string& filepath);

} // namespace igor::modelio::structure


namespace igor::modelio::weights {

/**
 * @brief Read model weights from JSON file
 *
 * @param filepath Path to model_weights.json file
 * @return WeightsData containing event weights and error rate
 * @throws std::runtime_error if file cannot be read or parsed
 */
MODELIO_EXPORT
WeightsData read_json(const std::string& filepath);

/**
 * @brief Read model weights from legacy model_marginals.txt file
 *
 * @param filepath Path to model_marginals.txt file
 * @return WeightsData containing event weights
 * @throws std::runtime_error if file cannot be read or parsed
 */
MODELIO_EXPORT
WeightsData read_legacy(const std::string& filepath);

} // namespace igor::modelio::weights
