/*
 * ModelWriter.h
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
 * @file ModelWriter.h
 * @brief Writer functions for IGoR model files
 *
 * This module provides functionality to write IGoR model files in both
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
 *   // Write JSON format
 *   structure::write_json("model_struct.json", model);
 *   weights::write_json("model_weights.json", weights);
 *
 *   // Write legacy format
 *   structure::write_legacy("model_parms.txt", model);
 *   weights::write_legacy("model_marginals.txt", weights);
 * @endcode
 */

#pragma once

#include <igor/ModelIO/ModelIOCommon.h>

#include <string>

namespace igor::modelio::structure {

/**
 * @brief Write model structure to JSON file
 *
 * @param filepath Path to output model_struct.json file
 * @param model ModelData to write
 * @param pretty Enable pretty-printing with indentation (default: true)
 * @throws std::runtime_error if file cannot be written
 */
MODELIO_EXPORT
void write_json(
    const std::string& filepath,
    const ModelData& model,
    bool pretty = true);

/**
 * @brief Write model structure to legacy model_parms.txt file
 *
 * @param filepath Path to output model_parms.txt file
 * @param model ModelData to write
 * @throws std::runtime_error if file cannot be written
 */
MODELIO_EXPORT
void write_legacy(
    const std::string& filepath,
    const ModelData& model);

} // namespace igor::modelio::structure


namespace igor::modelio::weights {

/**
 * @brief Write model weights to JSON file
 *
 * @param filepath Path to output model_weights.json file
 * @param weights WeightsData to write
 * @param pretty Enable pretty-printing with indentation (default: true)
 * @throws std::runtime_error if file cannot be written
 */
MODELIO_EXPORT
void write_json(
    const std::string& filepath,
    const WeightsData& weights,
    bool pretty = true);

/**
 * @brief Write model weights to legacy model_marginals.txt file
 *
 * @param filepath Path to output model_marginals.txt file
 * @param weights WeightsData to write
 * @throws std::runtime_error if file cannot be written
 */
MODELIO_EXPORT
void write_legacy(
    const std::string& filepath,
    const WeightsData& weights);

} // namespace igor::modelio::weights
