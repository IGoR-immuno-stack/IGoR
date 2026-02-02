/*
 * SequenceTypes.h
 *
 *  Created on: Nov 26, 2025
 *      Author: IGoR Team
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 * analyze and model immune receptors generation, selection, mutation and all
 * other processes. Copyright (C) 2017  Quentin Marcou
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
 *
 */

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

class SequenceTypeRegistry
{
public:
    using TypeId = uint16_t;

    // Predefined types for backward compatibility
    static constexpr TypeId V_GENE_SEQ = 0;
    static constexpr TypeId VD_INS_SEQ = 1;
    static constexpr TypeId D_GENE_SEQ = 2;
    static constexpr TypeId DJ_INS_SEQ = 3;
    static constexpr TypeId J_GENE_SEQ = 4;
    static constexpr TypeId VJ_INS_SEQ = 5;

    // Dynamic type allocation
    TypeId register_type(const std::string &name);
    TypeId get_type_id(const std::string &name) const;
    std::string get_type_name(TypeId id) const;

    // For tandem D support
    TypeId register_d_gene(int d_index); // D1, D2, etc.
    TypeId register_d_insertion(int from_d, int to_d); // D1D2_ins
    
    // For tandem D junctions
    TypeId register_junction_type(const std::string &name);
    TypeId register_junction_type(const std::string &name, TypeId upstream, TypeId downstream);

    size_t size() const { return id_to_name_.size(); }

    static SequenceTypeRegistry &get_instance();

    // Topology queries
    struct NeighborInfo
    {
        TypeId neighbor_type = -1;
        TypeId junction_type = -1;
        bool exists = false;
    };

    std::vector<NeighborInfo> get_upstream_neighbors(TypeId type_id) const;
    std::vector<NeighborInfo> get_downstream_neighbors(TypeId type_id) const;

    // Registration helpers that also build topology
    void register_connection(TypeId upstream, TypeId downstream, TypeId junction);

    // Helpers for insertion/junction handling
    bool is_junction_type(TypeId type_id) const;
    bool is_gene_type(TypeId type_id) const;
    TypeId get_junction_upstream(TypeId junction_type) const;
    TypeId get_junction_downstream(TypeId junction_type) const;

    struct JunctionNeighbor {
        TypeId neighbor_type;
        bool is_upstream;
    };
    std::vector<JunctionNeighbor> get_neighbors_for_junction(TypeId junction_type) const;

    // Safe type lookup (returns -1 if not found instead of throwing)
    int try_get_type_id(const std::string &name) const;

    // Direct access to internal structures for iteration
    struct TypeInfo {
        TypeId id;
        std::string name;
        bool is_gene;
        bool is_junction;
    };
    std::vector<TypeInfo> get_all_types() const;
    TypeInfo get_type_info(TypeId id) const;

private:
    SequenceTypeRegistry();

    std::unordered_map<std::string, TypeId> name_to_id_;
    std::vector<std::string> id_to_name_;
    TypeId next_id_ = 6; // Start after predefined types

    // Topology storage
    std::unordered_map<TypeId, std::vector<NeighborInfo>> upstream_neighbors_;
    std::unordered_map<TypeId, std::vector<NeighborInfo>> downstream_neighbors_;

    // Junction topology (junction -> upstream, downstream)
    std::unordered_map<TypeId, std::pair<TypeId, TypeId>> junction_endpoints_;
};

// Global registry accessor
inline SequenceTypeRegistry &get_sequence_type_registry()
{
    return SequenceTypeRegistry::get_instance();
}
