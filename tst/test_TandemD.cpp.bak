/*
 * test_TandemD.cpp
 *
 *  Unit tests for tandem D gene support.
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
 */

#include <catch2/catch_test_macros.hpp>
#include <igor/Core/Deletion.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/SequenceTypes.h>
#include <memory>
#include <string>

using namespace std;

TEST_CASE("Tandem D Type Registration", "[tandem_d]")
{
    auto &registry = SequenceTypeRegistry::get_instance();

    SECTION("Register D1 and D2 gene types")
    {
        auto d1_type = registry.register_type("D1_gene_seq");
        auto d2_type = registry.register_type("D2_gene_seq");

        REQUIRE(d1_type != d2_type);
        REQUIRE(registry.get_type_name(d1_type) == "D1_gene_seq");
        REQUIRE(registry.get_type_name(d2_type) == "D2_gene_seq");
    }

    SECTION("Register D1D2 insertion type")
    {
        auto d1_type = registry.register_type("D1_gene_seq");
        auto d2_type = registry.register_type("D2_gene_seq");
        auto d1d2_ins_type = registry.register_type("D1D2_ins_seq");

        // Register the connection
        registry.register_connection(d1_type, d2_type, d1d2_ins_type);

        // Verify topology
        REQUIRE(registry.is_junction_type(d1d2_ins_type));
        REQUIRE(registry.get_junction_upstream(d1d2_ins_type) == d1_type);
        REQUIRE(registry.get_junction_downstream(d1d2_ins_type) == d2_type);
    }

    SECTION("Register VD1 and D2J insertions")
    {
        auto d1_type = registry.register_type("D1_gene_seq");
        auto d2_type = registry.register_type("D2_gene_seq");
        auto vd1_ins_type = registry.register_type("VD1_ins_seq");
        auto d2j_ins_type = registry.register_type("D2J_ins_seq");

        registry.register_connection(SequenceTypeRegistry::V_GENE_SEQ, d1_type, vd1_ins_type);
        registry.register_connection(d2_type, SequenceTypeRegistry::J_GENE_SEQ, d2j_ins_type);

        REQUIRE(registry.get_junction_upstream(vd1_ins_type) == SequenceTypeRegistry::V_GENE_SEQ);
        REQUIRE(registry.get_junction_downstream(vd1_ins_type) == d1_type);
        REQUIRE(registry.get_junction_upstream(d2j_ins_type) == d2_type);
        REQUIRE(registry.get_junction_downstream(d2j_ins_type) == SequenceTypeRegistry::J_GENE_SEQ);
    }
}

TEST_CASE("Tandem D Topology", "[tandem_d]")
{
    auto &registry = SequenceTypeRegistry::get_instance();

    // Set up complete tandem D topology
    auto d1_type = registry.register_type("D1_gene_seq");
    auto d2_type = registry.register_type("D2_gene_seq");
    auto vd1_ins_type = registry.register_type("VD1_ins_seq");
    auto d1d2_ins_type = registry.register_type("D1D2_ins_seq");
    auto d2j_ins_type = registry.register_type("D2J_ins_seq");

    registry.register_connection(SequenceTypeRegistry::V_GENE_SEQ, d1_type, vd1_ins_type);
    registry.register_connection(d1_type, d2_type, d1d2_ins_type);
    registry.register_connection(d2_type, SequenceTypeRegistry::J_GENE_SEQ, d2j_ins_type);

    SECTION("D1 has correct neighbors")
    {
        auto d1_upstream = registry.get_upstream_neighbors(d1_type);
        auto d1_downstream = registry.get_downstream_neighbors(d1_type);

        // D1 should have V as upstream via VD1_ins
        bool has_v_upstream = false;
        for (const auto &neighbor : d1_upstream) {
            if (neighbor.neighbor_type == SequenceTypeRegistry::V_GENE_SEQ) {
                has_v_upstream = true;
                REQUIRE(neighbor.junction_type == vd1_ins_type);
            }
        }
        REQUIRE(has_v_upstream);

        // D1 should have D2 as downstream via D1D2_ins
        bool has_d2_downstream = false;
        for (const auto &neighbor : d1_downstream) {
            if (neighbor.neighbor_type == d2_type) {
                has_d2_downstream = true;
                REQUIRE(neighbor.junction_type == d1d2_ins_type);
            }
        }
        REQUIRE(has_d2_downstream);
    }

    SECTION("D2 has correct neighbors")
    {
        auto d2_upstream = registry.get_upstream_neighbors(d2_type);
        auto d2_downstream = registry.get_downstream_neighbors(d2_type);

        // D2 should have D1 as upstream via D1D2_ins
        bool has_d1_upstream = false;
        for (const auto &neighbor : d2_upstream) {
            if (neighbor.neighbor_type == d1_type) {
                has_d1_upstream = true;
                REQUIRE(neighbor.junction_type == d1d2_ins_type);
            }
        }
        REQUIRE(has_d1_upstream);

        // D2 should have J as downstream via D2J_ins
        bool has_j_downstream = false;
        for (const auto &neighbor : d2_downstream) {
            if (neighbor.neighbor_type == SequenceTypeRegistry::J_GENE_SEQ) {
                has_j_downstream = true;
                REQUIRE(neighbor.junction_type == d2j_ins_type);
            }
        }
        REQUIRE(has_j_downstream);
    }
}

TEST_CASE("Standard VDJ Backward Compatibility", "[tandem_d][regression]")
{
    auto &registry = SequenceTypeRegistry::get_instance();

    SECTION("Standard types are always available")
    {
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::V_GENE_SEQ) == "V_gene_seq");
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::VD_INS_SEQ) == "VD_ins_seq");
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::D_GENE_SEQ) == "D_gene_seq");
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::DJ_INS_SEQ) == "DJ_ins_seq");
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::J_GENE_SEQ) == "J_gene_seq");
        REQUIRE(registry.get_type_name(SequenceTypeRegistry::VJ_INS_SEQ) == "VJ_ins_seq");
    }

    SECTION("Standard topology is pre-configured")
    {
        // V -> D via VD_ins
        auto v_downstream = registry.get_downstream_neighbors(SequenceTypeRegistry::V_GENE_SEQ);
        bool has_d = false;
        for (const auto &neighbor : v_downstream) {
            if (neighbor.neighbor_type == SequenceTypeRegistry::D_GENE_SEQ) {
                has_d = true;
                REQUIRE(neighbor.junction_type == SequenceTypeRegistry::VD_INS_SEQ);
            }
        }
        REQUIRE(has_d);

        // D -> J via DJ_ins
        auto d_downstream = registry.get_downstream_neighbors(SequenceTypeRegistry::D_GENE_SEQ);
        bool has_j = false;
        for (const auto &neighbor : d_downstream) {
            if (neighbor.neighbor_type == SequenceTypeRegistry::J_GENE_SEQ) {
                has_j = true;
                REQUIRE(neighbor.junction_type == SequenceTypeRegistry::DJ_INS_SEQ);
            }
        }
        REQUIRE(has_j);

        // V -> J via VJ_ins (for VJ models)
        bool has_vj = false;
        for (const auto &neighbor : v_downstream) {
            if (neighbor.neighbor_type == SequenceTypeRegistry::J_GENE_SEQ) {
                has_vj = true;
                REQUIRE(neighbor.junction_type == SequenceTypeRegistry::VJ_INS_SEQ);
            }
        }
        REQUIRE(has_vj);
    }
}

