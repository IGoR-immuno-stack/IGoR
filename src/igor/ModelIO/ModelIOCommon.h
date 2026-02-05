/*
 * ModelIOCommon.h
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
 * @file ModelIOCommon.h
 * @brief Common types for Model I/O operations
 *
 * This header defines the data structures used for reading and writing
 * IGoR model files in both JSON and legacy formats.
 *
 * **File Formats:**
 * - model_struct.json: Immutable model structure (events, edges, sequence types)
 * - model_weights.json: Learned probability weights and error rate parameters
 * - model_parms.txt: Legacy structure format
 * - model_marginals.txt: Legacy weights format
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include <unordered_map>

#include <igor/ModelIO/Export.h>

namespace igor::modelio {

//==============================================================================
// Format version
//==============================================================================

/// Current format version for model files
inline constexpr const char* FORMAT_VERSION = "1.0.0";

//==============================================================================
// Realization data
//==============================================================================

/**
 * @brief Data for a single realization of an event
 *
 * For GeneChoice events, this represents a gene with name and sequence.
 * For Deletion/Insertion events, this represents a count value.
 */
struct MODELIO_EXPORT RealizationData {
    int index;                          ///< Index in the realization list
    std::string name;                   ///< Gene name or label
    std::string sequence;               ///< Gene sequence (empty for non-gene events)
    std::optional<int> value;           ///< Numeric value (for deletion/insertion counts)
};

//==============================================================================
// Event data
//==============================================================================

/**
 * @brief Complete description of a recombination event
 *
 * Events are the nodes in IGoR's Bayesian network. Each event has
 * a type (GeneChoice, Deletion, Insertion, DinucMarkov) and a set
 * of possible realizations.
 */
struct MODELIO_EXPORT EventData {
    std::string id;                     ///< Unique identifier (nickname)
    std::string legacy_name;            ///< Original name for legacy compatibility
    std::string type;                   ///< Event type: GeneChoice, Deletion, Insertion, DinucMarkov
    std::string gene_class;             ///< Gene class: V, D, J, or empty
    std::string seq_side;               ///< Sequence side: Five_prime, Three_prime, Undefined_side
    int priority;                       ///< Processing priority
    std::vector<RealizationData> realizations;  ///< Possible realizations
};

//==============================================================================
// Edge data
//==============================================================================

/**
 * @brief Edge in the Bayesian network
 *
 * Represents a parent-child dependency between two events.
 */
struct MODELIO_EXPORT EdgeData {
    std::string parent;                 ///< Parent event ID
    std::string child;                  ///< Child event ID
};

//==============================================================================
// Sequence type data
//==============================================================================

/**
 * @brief Definition of a sequence type
 */
struct MODELIO_EXPORT SequenceTypeDefinition {
    int id;                             ///< Numeric ID
    std::vector<std::string> aliases;   ///< Alternative names
    std::vector<std::string> parents;   ///< Parent sequence types
    std::vector<std::string> children;  ///< Child sequence types
};

/**
 * @brief Sequence type ordering and definitions
 */
struct MODELIO_EXPORT SequenceTypeData {
    std::vector<std::string> order;     ///< Ordered list of sequence type names
    std::unordered_map<std::string, SequenceTypeDefinition> definitions;
};

//==============================================================================
// Model metadata
//==============================================================================

/**
 * @brief Metadata about the model
 */
struct MODELIO_EXPORT ModelMetadata {
    std::string species;                ///< Species (e.g., "Homo sapiens")
    std::string chain;                  ///< Chain (e.g., "TRB", "IGH")
    std::string model_type;             ///< Model type (e.g., "VDJ", "VJ")
    std::string created;                ///< Creation timestamp (ISO 8601)
    std::string description;            ///< Human-readable description
};

//==============================================================================
// Model structure (model_struct.json)
//==============================================================================

/**
 * @brief Complete model structure data
 *
 * Contains all immutable information about a model:
 * events, edges, and sequence type definitions.
 */
struct MODELIO_EXPORT ModelData {
    std::string format_version;         ///< Format version string
    ModelMetadata metadata;             ///< Model metadata
    std::vector<EventData> events;      ///< List of events
    std::vector<EdgeData> edges;        ///< List of edges
    SequenceTypeData sequence_types;    ///< Sequence type information
};

//==============================================================================
// Event weights
//==============================================================================

/**
 * @brief Probability weights for a single event
 *
 * Contains the conditional probability tensor for an event,
 * with dimensions corresponding to the event's realizations
 * and any conditioning events.
 */
struct MODELIO_EXPORT EventWeights {
    std::string event_id;               ///< Event identifier
    std::vector<size_t> dims;           ///< Tensor dimensions
    std::vector<std::string> conditioning_events;  ///< Conditioning event IDs
    std::vector<double> values;         ///< Flattened probability values
    bool normalized;                    ///< Whether values are normalized
};

//==============================================================================
// Error rate data
//==============================================================================

/**
 * @brief Error rate model parameters
 *
 * Supports different error rate models:
 * - SingleErrorRate: Single error probability
 * - HypermutationGlobal: Global hypermutation model
 * - HypermutationFullNmer: Full N-mer hypermutation model
 */
struct MODELIO_EXPORT ErrorRateData {
    std::string type;                   ///< Error rate type
    double rate;                        ///< Error rate value
    std::string learn_on;               ///< Sequence type to learn on
    std::string apply_on;               ///< Sequence type to apply to
    std::optional<int> nmer_size;       ///< N-mer size (for hypermutation)
    std::vector<double> ei_contributions;    ///< Error-insertion contributions
    std::vector<double> mutation_probas;     ///< Mutation probabilities
};

//==============================================================================
// Model weights (model_weights.json)
//==============================================================================

/**
 * @brief Complete model weights data
 *
 * Contains all learned probability information:
 * event weights and error rate parameters.
 */
struct MODELIO_EXPORT WeightsData {
    std::string format_version;         ///< Format version string
    std::string struct_version;         ///< Compatible struct version
    std::vector<EventWeights> events;   ///< Event probability tensors
    ErrorRateData error_rate;           ///< Error rate parameters
};

} // namespace igor::modelio
