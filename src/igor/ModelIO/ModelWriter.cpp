/*
 * ModelWriter.cpp
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

#include <igor/ModelIO/ModelWriter.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>

using json = nlohmann::json;

namespace igor::modelio::detail {

//==============================================================================
// Helper functions
//==============================================================================

std::string generate_legacy_name(const EventData& event)
{
    std::stringstream ss;
    ss << event.type << "_";

    if (!event.gene_class.empty()) {
        ss << event.gene_class << "_gene_";
    }

    ss << event.seq_side << "_";
    ss << "prio" << event.priority << "_";
    ss << "size" << event.realizations.size();

    return ss.str();
}

} // namespace igor::modelio::detail

namespace igor::modelio::structure {

//==============================================================================
// JSON Writing
//==============================================================================

void write_json(
    const std::string& filepath,
    const ModelData& model,
    bool pretty)
{
    json j;

    // Format version
    j["format_version"] = model.format_version.empty() ? FORMAT_VERSION : model.format_version;

    // Metadata
    json meta;
    meta["species"] = model.metadata.species;
    meta["chain"] = model.metadata.chain;
    meta["model_type"] = model.metadata.model_type;
    meta["created"] = model.metadata.created;
    meta["description"] = model.metadata.description;
    j["metadata"] = meta;

    // Events
    json events = json::array();
    for (const auto& ev : model.events) {
        json event;
        event["id"] = ev.id;
        if (!ev.legacy_name.empty()) {
            event["legacy_name"] = ev.legacy_name;
        }
        event["type"] = ev.type;
        event["gene_class"] = ev.gene_class;
        event["seq_side"] = ev.seq_side;
        event["priority"] = ev.priority;

        json realizations = json::array();
        for (const auto& r : ev.realizations) {
            json real;
            real["index"] = r.index;
            real["name"] = r.name;
            if (!r.sequence.empty()) {
                real["sequence"] = r.sequence;
            }
            if (r.value.has_value()) {
                real["value"] = r.value.value();
            }
            realizations.push_back(real);
        }
        event["realizations"] = realizations;

        events.push_back(event);
    }
    j["events"] = events;

    // Edges
    json edges = json::array();
    for (const auto& e : model.edges) {
        json edge;
        edge["parent"] = e.parent;
        edge["child"] = e.child;
        edges.push_back(edge);
    }
    j["edges"] = edges;

    // Sequence types
    json st;
    st["order"] = model.sequence_types.order;

    json defs;
    for (const auto& [name, def] : model.sequence_types.definitions) {
        json d;
        d["id"] = def.id;
        if (!def.aliases.empty()) {
            d["aliases"] = def.aliases;
        }
        if (!def.parents.empty()) {
            d["parents"] = def.parents;
        }
        if (!def.children.empty()) {
            d["children"] = def.children;
        }
        defs[name] = d;
    }
    st["definitions"] = defs;
    j["sequence_types"] = st;

    // Write to file
    std::ofstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    if (pretty) {
        file << std::setw(2) << j << std::endl;
    } else {
        file << j << std::endl;
    }
}

//==============================================================================
// Legacy Writing
//==============================================================================

void write_legacy(
    const std::string& filepath,
    const ModelData& model)
{
    std::ofstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    // Write Event_list section
    file << "@Event_list\n";

    for (const auto& ev : model.events) {
        // Write event header
        std::string header = "#" + (ev.legacy_name.empty() ? detail::generate_legacy_name(ev) : ev.legacy_name);
        file << header << "\n";

        // Write realizations
        for (const auto& r : ev.realizations) {
            if (ev.type == "GeneChoice") {
                file << "%" << r.name << ";" << r.sequence << ";" << r.index << "\n";
            } else if (ev.type == "Deletion" || ev.type == "Insertion") {
                file << "%" << (r.value.has_value() ? std::to_string(r.value.value()) : r.name) << "\n";
            }
        }
    }

    // Write Edges section
    file << "@Edges\n";
    for (const auto& e : model.edges) {
        file << "%" << e.parent << ";" << e.child << "\n";
    }

    // Write ErrorRate section (placeholder - actual error rate in weights)
    file << "@ErrorRate\n";
    file << "#SingleErrorRate\n";
}

} // namespace igor::modelio::structure


namespace igor::modelio::weights {

//==============================================================================
// JSON Writing
//==============================================================================

void write_json(
    const std::string& filepath,
    const WeightsData& weights,
    bool pretty)
{
    json j;

    // Format version
    j["format_version"] = weights.format_version.empty() ? FORMAT_VERSION : weights.format_version;
    j["struct_version"] = weights.struct_version;

    // Event weights
    json events = json::array();
    for (const auto& ew : weights.events) {
        json event;
        event["event_id"] = ew.event_id;
        event["dims"] = ew.dims;
        event["conditioning_events"] = ew.conditioning_events;
        event["values"] = ew.values;
        event["normalized"] = ew.normalized;
        events.push_back(event);
    }
    j["events"] = events;

    // Error rate
    json er;
    er["type"] = weights.error_rate.type;
    er["rate"] = weights.error_rate.rate;
    er["learn_on"] = weights.error_rate.learn_on;
    er["apply_on"] = weights.error_rate.apply_on;

    if (weights.error_rate.nmer_size.has_value()) {
        er["nmer_size"] = weights.error_rate.nmer_size.value();
    }
    if (!weights.error_rate.ei_contributions.empty()) {
        er["ei_contributions"] = weights.error_rate.ei_contributions;
    }
    if (!weights.error_rate.mutation_probas.empty()) {
        er["mutation_probas"] = weights.error_rate.mutation_probas;
    }

    j["error_rate"] = er;

    // Write to file
    std::ofstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    if (pretty) {
        file << std::setw(2) << j << std::endl;
    } else {
        file << j << std::endl;
    }
}

//==============================================================================
// Legacy Writing
//==============================================================================

void write_legacy(
    const std::string& filepath,
    const WeightsData& weights)
{
    std::ofstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    for (const auto& ew : weights.events) {
        // Write event header
        file << "@" << ew.event_id << "\n";

        // Write dimensions
        file << "$Dim[";
        for (size_t i = 0; i < ew.dims.size(); ++i) {
            if (i > 0) file << ",";
            file << ew.dims[i];
        }
        file << "]\n";

        // Write values with conditioning headers
        if (ew.dims.size() == 1) {
            // Simple 1D case
            file << "#\n";
            file << "%";
            for (size_t i = 0; i < ew.values.size(); ++i) {
                if (i > 0) file << ",";
                file << ew.values[i];
            }
            file << "\n";
        } else if (ew.dims.size() == 2) {
            // 2D case with conditioning
            size_t rows = ew.dims[0];
            size_t cols = ew.dims[1];

            for (size_t i = 0; i < rows; ++i) {
                // Conditioning header
                if (!ew.conditioning_events.empty()) {
                    file << "#[" << ew.conditioning_events[0] << "," << i << "]\n";
                } else {
                    file << "#[row," << i << "]\n";
                }

                // Values for this row
                file << "%";
                for (size_t j = 0; j < cols; ++j) {
                    if (j > 0) file << ",";
                    file << ew.values[i * cols + j];
                }
                file << "\n";
            }
        } else {
            // Higher dimensional - flatten with row headers
            // This is a simplification; real implementation would need proper handling
            file << "#\n";
            file << "%";
            for (size_t i = 0; i < ew.values.size(); ++i) {
                if (i > 0) file << ",";
                file << ew.values[i];
            }
            file << "\n";
        }
    }
}

} // namespace igor::modelio::weights
