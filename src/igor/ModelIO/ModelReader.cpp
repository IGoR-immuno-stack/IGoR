/*
 * ModelReader.cpp
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

#include <igor/ModelIO/ModelReader.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

using json = nlohmann::json;

namespace igor::modelio::detail {

//==============================================================================
// Helper functions
//==============================================================================

std::string generate_nickname(const EventData& event)
{
    // Generate a short, unique nickname from event properties
    std::string nickname;

    if (!event.gene_class.empty()) {
        nickname = event.gene_class.substr(0, 1);  // V, D, or J
    }

    if (event.type == "GeneChoice") {
        nickname += "_choice";
    } else if (event.type == "Deletion") {
        if (event.seq_side == "Three_prime") {
            nickname += "_3_del";
        } else if (event.seq_side == "Five_prime") {
            nickname += "_5_del";
        } else {
            nickname += "_del";
        }
    } else if (event.type == "Insertion") {
        nickname = "VD_ins";  // or DJ_ins based on context
        if (event.gene_class == "D" && event.seq_side == "Five_prime") {
            nickname = "VD_ins";
        } else if (event.gene_class == "D" && event.seq_side == "Three_prime") {
            nickname = "DJ_ins";
        }
    } else if (event.type == "DinucMarkov") {
        nickname = "dinuc_markov";
    }

    return nickname;
}

} // namespace igor::modelio::detail

namespace igor::modelio::structure {

//==============================================================================
// JSON Reading
//==============================================================================

ModelData read_json(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }

    ModelData model;

    // Format version
    model.format_version = j.value("format_version", "1.0.0");

    // Metadata
    if (j.contains("metadata")) {
        const auto& meta = j["metadata"];
        model.metadata.species = meta.value("species", "");
        model.metadata.chain = meta.value("chain", "");
        model.metadata.model_type = meta.value("model_type", "");
        model.metadata.created = meta.value("created", "");
        model.metadata.description = meta.value("description", "");
    }

    // Events
    if (j.contains("events")) {
        for (const auto& ev : j["events"]) {
            EventData event;
            event.id = ev.value("id", "");
            event.legacy_name = ev.value("legacy_name", "");
            event.type = ev.value("type", "");
            event.gene_class = ev.value("gene_class", "");
            event.seq_side = ev.value("seq_side", "");
            event.priority = ev.value("priority", 0);

            if (ev.contains("realizations")) {
                for (const auto& real : ev["realizations"]) {
                    RealizationData r;
                    r.index = real.value("index", 0);
                    r.name = real.value("name", "");
                    r.sequence = real.value("sequence", "");
                    if (real.contains("value")) {
                        r.value = real["value"].get<int>();
                    }
                    event.realizations.push_back(r);
                }
            }

            model.events.push_back(event);
        }
    }

    // Edges
    if (j.contains("edges")) {
        for (const auto& edge : j["edges"]) {
            EdgeData e;
            e.parent = edge.value("parent", "");
            e.child = edge.value("child", "");
            model.edges.push_back(e);
        }
    }

    // Sequence types
    if (j.contains("sequence_types")) {
        const auto& st = j["sequence_types"];
        if (st.contains("order")) {
            model.sequence_types.order = st["order"].get<std::vector<std::string>>();
        }
        if (st.contains("definitions")) {
            for (const auto& [name, def] : st["definitions"].items()) {
                SequenceTypeDefinition d;
                d.id = def.value("id", 0);
                if (def.contains("aliases")) {
                    d.aliases = def["aliases"].get<std::vector<std::string>>();
                }
                if (def.contains("parents")) {
                    d.parents = def["parents"].get<std::vector<std::string>>();
                }
                if (def.contains("children")) {
                    d.children = def["children"].get<std::vector<std::string>>();
                }
                model.sequence_types.definitions[name] = d;
            }
        }
    }

    return model;
}

//==============================================================================
// Legacy Reading
//==============================================================================

ModelData read_legacy(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    ModelData model;
    model.format_version = FORMAT_VERSION;

    std::string line;
    std::string current_section;
    EventData* current_event = nullptr;

    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Section headers
        if (line[0] == '@') {
            current_section = line.substr(1);
            continue;
        }

        // Event headers
        if (line[0] == '#') {
            std::string event_header = line.substr(1);

            EventData event;
            event.legacy_name = event_header;

            // Parse event type from header
            if (event_header.find("GeneChoice") != std::string::npos) {
                event.type = "GeneChoice";
            } else if (event_header.find("Deletion") != std::string::npos) {
                event.type = "Deletion";
            } else if (event_header.find("Insertion") != std::string::npos) {
                event.type = "Insertion";
            } else if (event_header.find("DinucMarkov") != std::string::npos) {
                event.type = "DinucMarkov";
            }

            // Parse gene class
            if (event_header.find("_V_") != std::string::npos ||
                event_header.find("V_gene") != std::string::npos) {
                event.gene_class = "V";
            } else if (event_header.find("_D_") != std::string::npos ||
                       event_header.find("D_gene") != std::string::npos) {
                event.gene_class = "D";
            } else if (event_header.find("_J_") != std::string::npos ||
                       event_header.find("J_gene") != std::string::npos) {
                event.gene_class = "J";
            }

            // Parse seq_side
            if (event_header.find("Five_prime") != std::string::npos) {
                event.seq_side = "Five_prime";
            } else if (event_header.find("Three_prime") != std::string::npos) {
                event.seq_side = "Three_prime";
            } else {
                event.seq_side = "Undefined_side";
            }

            // Extract priority from header (e.g., "prio7")
            auto prio_pos = event_header.find("prio");
            if (prio_pos != std::string::npos) {
                auto num_start = prio_pos + 4;
                auto num_end = event_header.find('_', num_start);
                if (num_end == std::string::npos) {
                    num_end = event_header.length();
                }
                try {
                    event.priority = std::stoi(event_header.substr(num_start, num_end - num_start));
                } catch (...) {
                    event.priority = 0;
                }
            }

            // Generate nickname
            event.id = detail::generate_nickname(event);

            model.events.push_back(event);
            current_event = &model.events.back();
            continue;
        }

        // Realizations (lines starting with %)
        if (line[0] == '%' && current_event != nullptr) {
            std::string realization_line = line.substr(1);

            RealizationData r;
            r.index = static_cast<int>(current_event->realizations.size());

            // Parse realization based on event type
            if (current_event->type == "GeneChoice") {
                // Format: name;sequence;index
                auto first_semi = realization_line.find(';');
                auto last_semi = realization_line.rfind(';');

                if (first_semi != std::string::npos && last_semi != first_semi) {
                    r.name = realization_line.substr(0, first_semi);
                    r.sequence = realization_line.substr(first_semi + 1, last_semi - first_semi - 1);
                    // Index is at end, but we use our own indexing
                }
            } else if (current_event->type == "Deletion" || current_event->type == "Insertion") {
                // Format: value
                try {
                    r.value = std::stoi(realization_line);
                    r.name = realization_line;
                } catch (...) {
                    r.name = realization_line;
                }
            }

            current_event->realizations.push_back(r);
            continue;
        }

        // Edge definitions in @Edges section
        if (current_section == "Edges" && line[0] == '%') {
            std::string edge_line = line.substr(1);
            auto semicolon = edge_line.find(';');
            if (semicolon != std::string::npos) {
                EdgeData edge;
                edge.parent = edge_line.substr(0, semicolon);
                edge.child = edge_line.substr(semicolon + 1);
                model.edges.push_back(edge);
            }
        }
    }

    return model;
}

} // namespace igor::modelio::structure


namespace igor::modelio::weights {

//==============================================================================
// JSON Reading
//==============================================================================

WeightsData read_json(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }

    WeightsData weights;

    // Format version
    weights.format_version = j.value("format_version", "1.0.0");
    weights.struct_version = j.value("struct_version", "1.0.0");

    // Event weights
    if (j.contains("events")) {
        for (const auto& ev : j["events"]) {
            EventWeights ew;
            ew.event_id = ev.value("event_id", "");

            if (ev.contains("dims")) {
                ew.dims = ev["dims"].get<std::vector<size_t>>();
            }
            if (ev.contains("conditioning_events")) {
                ew.conditioning_events = ev["conditioning_events"].get<std::vector<std::string>>();
            }
            if (ev.contains("values")) {
                ew.values = ev["values"].get<std::vector<double>>();
            }
            ew.normalized = ev.value("normalized", true);

            weights.events.push_back(ew);
        }
    }

    // Error rate
    if (j.contains("error_rate")) {
        const auto& er = j["error_rate"];
        weights.error_rate.type = er.value("type", "SingleErrorRate");
        weights.error_rate.rate = er.value("rate", 0.0);
        weights.error_rate.learn_on = er.value("learn_on", "seq");
        weights.error_rate.apply_on = er.value("apply_on", "seq");

        if (er.contains("nmer_size")) {
            weights.error_rate.nmer_size = er["nmer_size"].get<int>();
        }
        if (er.contains("ei_contributions")) {
            weights.error_rate.ei_contributions = er["ei_contributions"].get<std::vector<double>>();
        }
        if (er.contains("mutation_probas")) {
            weights.error_rate.mutation_probas = er["mutation_probas"].get<std::vector<double>>();
        }
    }

    return weights;
}

//==============================================================================
// Legacy Reading
//==============================================================================

WeightsData read_legacy(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    WeightsData weights;
    weights.format_version = FORMAT_VERSION;

    std::string line;
    EventWeights* current_event = nullptr;

    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Event header (@nickname)
        if (line[0] == '@') {
            EventWeights ew;
            ew.event_id = line.substr(1);
            ew.normalized = true;
            weights.events.push_back(ew);
            current_event = &weights.events.back();
            continue;
        }

        // Dimensions ($Dim[...])
        if (line[0] == '$' && current_event != nullptr) {
            // Parse $Dim[n] or $Dim[n,m,...]
            auto bracket_start = line.find('[');
            auto bracket_end = line.find(']');
            if (bracket_start != std::string::npos && bracket_end != std::string::npos) {
                std::string dims_str = line.substr(bracket_start + 1, bracket_end - bracket_start - 1);
                std::stringstream ss(dims_str);
                std::string dim;
                while (std::getline(ss, dim, ',')) {
                    current_event->dims.push_back(std::stoull(dim));
                }
            }
            continue;
        }

        // Conditioning header (#[event,index])
        if (line[0] == '#' && current_event != nullptr) {
            // Parse conditioning events if present
            if (line.length() > 1 && line[1] == '[') {
                auto bracket_end = line.find(']');
                if (bracket_end != std::string::npos) {
                    std::string cond_str = line.substr(2, bracket_end - 2);
                    auto comma = cond_str.find(',');
                    if (comma != std::string::npos) {
                        std::string event_name = cond_str.substr(0, comma);
                        // Add to conditioning events if not already present
                        if (std::find(current_event->conditioning_events.begin(),
                                      current_event->conditioning_events.end(),
                                      event_name) == current_event->conditioning_events.end()) {
                            current_event->conditioning_events.push_back(event_name);
                        }
                    }
                }
            }
            continue;
        }

        // Probability values (%value,value,...)
        if (line[0] == '%' && current_event != nullptr) {
            std::string values_str = line.substr(1);
            std::stringstream ss(values_str);
            std::string val;
            while (std::getline(ss, val, ',')) {
                try {
                    current_event->values.push_back(std::stod(val));
                } catch (...) {
                    // Skip invalid values
                }
            }
            continue;
        }
    }

    return weights;
}

} // namespace igor::modelio::weights
