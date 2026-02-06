/*
 * SequenceTypes.cpp
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

#include <igor/Core/SequenceTypes.h>

SequenceTypeRegistry::SequenceTypeRegistry()
{
    // Initialize predefined types
    id_to_name_.resize(6);
    id_to_name_[V_GENE_SEQ] = "V_gene_seq";
    id_to_name_[VD_INS_SEQ] = "VD_ins_seq";
    id_to_name_[D_GENE_SEQ] = "D_gene_seq";
    id_to_name_[DJ_INS_SEQ] = "DJ_ins_seq";
    id_to_name_[J_GENE_SEQ] = "J_gene_seq";
    id_to_name_[VJ_INS_SEQ] = "VJ_ins_seq";

    name_to_id_["V_gene_seq"] = V_GENE_SEQ;
    name_to_id_["VD_ins_seq"] = VD_INS_SEQ;
    name_to_id_["D_gene_seq"] = D_GENE_SEQ;
    name_to_id_["DJ_ins_seq"] = DJ_INS_SEQ;
    name_to_id_["J_gene_seq"] = J_GENE_SEQ;
    name_to_id_["VJ_ins_seq"] = VJ_INS_SEQ;

    // Standard Topology
    register_connection(V_GENE_SEQ, D_GENE_SEQ, VD_INS_SEQ);
    register_connection(D_GENE_SEQ, J_GENE_SEQ, DJ_INS_SEQ);
    register_connection(V_GENE_SEQ, J_GENE_SEQ, VJ_INS_SEQ);
}

SequenceTypeRegistry &SequenceTypeRegistry::get_instance()
{
    static SequenceTypeRegistry instance;
    return instance;
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::register_type(const std::string &name)
{
    if (name_to_id_.count(name)) {
        return name_to_id_.at(name);
    }

    TypeId id = next_id_++;
    name_to_id_[name] = id;
    if (id_to_name_.size() <= id) {
        id_to_name_.resize(id + 1);
    }
    id_to_name_[id] = name;
    return id;
}

void SequenceTypeRegistry::register_alias(const std::string &alias, TypeId target_id)
{
    // Map the alias name to the same ID. Does NOT allocate a new ID.
    name_to_id_[alias] = target_id;
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::get_type_id(const std::string &name) const
{
    if (name_to_id_.count(name)) {
        return name_to_id_.at(name);
    }
    throw std::runtime_error("Unknown sequence type: " + name);
}

std::string SequenceTypeRegistry::get_type_name(TypeId id) const
{
    if (id < id_to_name_.size()) {
        return id_to_name_[id];
    }
    throw std::runtime_error("Unknown sequence type ID: " + std::to_string(id));
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::register_d_gene(int d_index)
{
    return register_type("D" + std::to_string(d_index) + "_gene_seq");
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::register_d_insertion(int from_d, int to_d)
{
    return register_type("D" + std::to_string(from_d) + "D" + std::to_string(to_d) + "_ins_seq");
}

void SequenceTypeRegistry::register_connection(TypeId upstream, TypeId downstream, TypeId junction)
{
    downstream_neighbors_[upstream].push_back({ downstream, junction, true });
    upstream_neighbors_[downstream].push_back({ upstream, junction, true });
    junction_endpoints_[junction] = { upstream, downstream };
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::register_junction_type(const std::string &name)
{
    if (name_to_id_.count(name)) {
        return name_to_id_[name];
    }
    TypeId id = id_to_name_.size();
    id_to_name_.push_back(name);
    name_to_id_[name] = id;
    // Set unknown endpoints for stand-alone registration
    junction_endpoints_[id] = { -1, -1 };
    return id;
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::register_junction_type(const std::string &name, TypeId upstream, TypeId downstream)
{
    TypeId id = register_junction_type(name);
    junction_endpoints_[id] = { upstream, downstream };
    return id;
}

std::vector<SequenceTypeRegistry::NeighborInfo>
SequenceTypeRegistry::get_upstream_neighbors(TypeId type_id) const
{
    if (upstream_neighbors_.count(type_id)) {
        return upstream_neighbors_.at(type_id);
    }
    return {};
}

std::vector<SequenceTypeRegistry::NeighborInfo>
SequenceTypeRegistry::get_downstream_neighbors(TypeId type_id) const
{
    if (downstream_neighbors_.count(type_id)) {
        return downstream_neighbors_.at(type_id);
    }
    return {};
}

bool SequenceTypeRegistry::is_junction_type(TypeId type_id) const
{
    return junction_endpoints_.count(type_id) > 0;
}

bool SequenceTypeRegistry::is_gene_type(TypeId type_id) const
{
    // A type is a gene if it's NOT a junction
    return !is_junction_type(type_id);
}

SequenceTypeRegistry::TypeInfo SequenceTypeRegistry::get_type_info(TypeId id) const
{
    if (id >= id_to_name_.size()) {
        throw std::out_of_range("SequenceTypeRegistry: invalid ID");
    }
    return { id, id_to_name_[id], is_gene_type(id), is_junction_type(id) };
}

std::vector<SequenceTypeRegistry::TypeInfo> SequenceTypeRegistry::get_all_types() const
{
    std::vector<TypeInfo> result;
    for (TypeId id = 0; id < id_to_name_.size(); ++id) {
        result.push_back(get_type_info(id));
    }
    return result;
}

SequenceTypeRegistry::TypeId SequenceTypeRegistry::get_junction_upstream(TypeId junction_type) const
{
    if (junction_endpoints_.count(junction_type)) {
        return junction_endpoints_.at(junction_type).first;
    }
    return static_cast<TypeId>(-1);
}

SequenceTypeRegistry::TypeId
SequenceTypeRegistry::get_junction_downstream(TypeId junction_type) const
{
    if (junction_endpoints_.count(junction_type)) {
        return junction_endpoints_.at(junction_type).second;
    }
    return static_cast<TypeId>(-1);
}

std::vector<SequenceTypeRegistry::JunctionNeighbor>
SequenceTypeRegistry::get_neighbors_for_junction(TypeId junction_type) const
{
    std::vector<JunctionNeighbor> result;
    if (junction_endpoints_.count(junction_type)) {
        auto endpoints = junction_endpoints_.at(junction_type);
        result.push_back({ endpoints.first, true });
        result.push_back({ endpoints.second, false });
    }
    return result;
}

int SequenceTypeRegistry::try_get_type_id(const std::string &name) const
{
    if (name_to_id_.count(name)) {
        return static_cast<int>(name_to_id_.at(name));
    }
    return -1;
}
