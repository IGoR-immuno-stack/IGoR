/*
 * ModelIOTestUtils.h
 *
 *  Created on: Feb 5, 2026
 *      Author: IGoR Development Team
 *
 *  Shared test utilities for ModelIO module tests
 */

#pragma once

#include <igor/ModelIO/ModelIOCommon.h>

#include <filesystem>
#include <string>
#include <vector>

namespace igor::modelio::test {

//==============================================================================
// RAII Test Directory Manager (replaces manual setup/cleanup)
//==============================================================================

/**
 * @brief RAII wrapper for test directory management
 *
 * Creates a unique temporary directory on construction,
 * removes it (and all contents) on destruction.
 * Ensures cleanup happens even if tests fail.
 */
class TestDirectory
{
public:
    explicit TestDirectory(const std::string& name = "igor_modelio_tests")
        : m_path(std::filesystem::temp_directory_path() / name)
    {
        std::filesystem::create_directories(m_path);
    }

    ~TestDirectory()
    {
        std::error_code ec;
        std::filesystem::remove_all(m_path, ec);
#ifndef NDEBUG
        if (ec) {
            std::cerr << "[TestDirectory] Warning: Failed to cleanup test directory '"
                      << m_path.string() << "': " << ec.message() << std::endl;
        }
#endif
    }

    // Non-copyable, non-movable
    TestDirectory(const TestDirectory&) = delete;
    TestDirectory& operator=(const TestDirectory&) = delete;
    TestDirectory(TestDirectory&&) = delete;
    TestDirectory& operator=(TestDirectory&&) = delete;

    /// Get full path to a file within the test directory
    [[nodiscard]] std::string path(const std::string& filename) const
    {
        return (m_path / filename).string();
    }

    /// Get the root test directory path
    [[nodiscard]] const std::filesystem::path& root() const { return m_path; }

private:
    std::filesystem::path m_path;
};

//==============================================================================
// Test Data Generators
//==============================================================================

/**
 * @brief Create a minimal valid ModelData for testing
 *
 * Creates a minimal model with:
 * - format_version: "1.0.0"
 * - metadata: chain="TRB", species="Homo sapiens"
 * - One v_choice event (no realizations)
 * - One edge (v_choice → v_del)
 *
 * @return A minimal but valid ModelData structure
 */
inline ModelData create_minimal_model()
{
    ModelData model;
    model.format_version = FORMAT_VERSION;
    model.metadata.chain = "TRB";
    model.metadata.species = "Homo sapiens";

    EventData event;
    event.id = "v_choice";
    event.type = "GeneChoice";
    event.gene_class = "V";
    event.seq_side = "Undefined_side";
    event.priority = 7;
    model.events.push_back(event);

    EdgeData edge;
    edge.parent = "v_choice";
    edge.child = "v_del";
    model.edges.push_back(edge);

    return model;
}

/**
 * @brief Create a ModelData with realizations
 *
 * Creates a model with:
 * - One v_choice event with `count` gene realizations
 * - Each realization has: TRBV{i+1}*01, sequence="ATCG...A"
 *
 * @param count Number of gene realizations to create
 * @return ModelData with populated realizations
 */
inline ModelData create_model_with_realizations(size_t count)
{
    ModelData model = create_minimal_model();
    model.events.clear();

    EventData event;
    event.id = "v_choice";
    event.type = "GeneChoice";
    event.gene_class = "V";
    event.seq_side = "Undefined_side";
    event.priority = 7;

    for (size_t i = 0; i < count; ++i) {
        RealizationData r;
        r.index = static_cast<int>(i);
        r.name = "TRBV" + std::to_string(i + 1) + "*01";
        r.sequence = "ATCG" + std::string(i + 1, 'A');
        event.realizations.push_back(r);
    }

    model.events.push_back(event);
    return model;
}

/**
 * @brief Create a complete VDJ recombination model
 *
 * Creates a realistic TRB model with:
 * - v_choice event (3 realizations)
 * - d_gene event (2 realizations)
 * - j_choice event (2 realizations)
 * - v_3_del, d_5_del, d_3_del, j_5_del deletion events
 * - vd_ins, dj_ins insertion events
 * - Full dependency graph (edges)
 * - Sequence types (V_gene, VD_ins, D_gene, DJ_ins, J_gene)
 *
 * @return A complete VDJ recombination ModelData
 */
inline ModelData create_vdj_model()
{
    ModelData model;
    model.format_version = FORMAT_VERSION;
    model.metadata.species = "Homo sapiens";
    model.metadata.chain = "TRB";
    model.metadata.model_type = "VDJ";
    model.metadata.description = "Test TRB recombination model";

    // V gene choice event with realizations
    EventData v_choice;
    v_choice.id = "v_choice";
    v_choice.type = "GeneChoice";
    v_choice.gene_class = "V";
    v_choice.seq_side = "Undefined_side";
    v_choice.priority = 7;

    for (int i = 0; i < 3; ++i) {
        RealizationData r;
        r.index = i;
        r.name = "TRBV" + std::to_string(i + 1) + "*01";
        r.sequence = "ATCGATCGATCG" + std::string(i + 1, 'A');
        v_choice.realizations.push_back(r);
    }
    model.events.push_back(v_choice);

    // D gene event with realizations
    EventData d_gene;
    d_gene.id = "d_gene";
    d_gene.type = "GeneChoice";
    d_gene.gene_class = "D";
    d_gene.seq_side = "Undefined_side";
    d_gene.priority = 5;

    for (int i = 0; i < 2; ++i) {
        RealizationData r;
        r.index = i;
        r.name = "TRBD" + std::to_string(i + 1) + "*01";
        r.sequence = "GGGAC";
        d_gene.realizations.push_back(r);
    }
    model.events.push_back(d_gene);

    // J gene choice event with realizations
    EventData j_choice;
    j_choice.id = "j_choice";
    j_choice.type = "GeneChoice";
    j_choice.gene_class = "J";
    j_choice.seq_side = "Undefined_side";
    j_choice.priority = 7;

    for (int i = 0; i < 2; ++i) {
        RealizationData r;
        r.index = i;
        r.name = "TRBJ" + std::to_string(i + 1) + "*01";
        r.sequence = "ACGTACGT";
        j_choice.realizations.push_back(r);
    }
    model.events.push_back(j_choice);

    // Deletion events
    EventData v_3_del;
    v_3_del.id = "v_3_del";
    v_3_del.type = "Deletion";
    v_3_del.gene_class = "V";
    v_3_del.seq_side = "Three_prime";
    v_3_del.priority = 6;
    model.events.push_back(v_3_del);

    EventData d_5_del;
    d_5_del.id = "d_5_del";
    d_5_del.type = "Deletion";
    d_5_del.gene_class = "D";
    d_5_del.seq_side = "Five_prime";
    d_5_del.priority = 4;
    model.events.push_back(d_5_del);

    EventData d_3_del;
    d_3_del.id = "d_3_del";
    d_3_del.type = "Deletion";
    d_3_del.gene_class = "D";
    d_3_del.seq_side = "Three_prime";
    d_3_del.priority = 4;
    model.events.push_back(d_3_del);

    EventData j_5_del;
    j_5_del.id = "j_5_del";
    j_5_del.type = "Deletion";
    j_5_del.gene_class = "J";
    j_5_del.seq_side = "Five_prime";
    j_5_del.priority = 6;
    model.events.push_back(j_5_del);

    // Insertion events
    EventData vd_ins;
    vd_ins.id = "vd_ins";
    vd_ins.type = "Insertion";
    vd_ins.seq_side = "VD_ins";
    vd_ins.priority = 3;
    model.events.push_back(vd_ins);

    EventData dj_ins;
    dj_ins.id = "dj_ins";
    dj_ins.type = "Insertion";
    dj_ins.seq_side = "DJ_ins";
    dj_ins.priority = 3;
    model.events.push_back(dj_ins);

    // Edges (dependency graph)
    model.edges = {
        {"v_choice", "v_3_del"},
        {"j_choice", "d_gene"},
        {"j_choice", "j_5_del"},
        {"d_gene", "d_5_del"},
        {"d_gene", "d_3_del"},
        {"v_3_del", "vd_ins"},
        {"d_5_del", "vd_ins"},
        {"d_3_del", "dj_ins"},
        {"j_5_del", "dj_ins"}
    };

    // Sequence types
    model.sequence_types.order = {"V_gene", "VD_ins", "D_gene", "DJ_ins", "J_gene"};

    SequenceTypeDefinition v_gene_def;
    v_gene_def.id = 0;
    v_gene_def.aliases = {"V"};
    model.sequence_types.definitions["V_gene"] = v_gene_def;

    SequenceTypeDefinition vd_ins_def;
    vd_ins_def.id = 1;
    vd_ins_def.parents = {"V_gene"};
    model.sequence_types.definitions["VD_ins"] = vd_ins_def;

    SequenceTypeDefinition d_gene_def;
    d_gene_def.id = 2;
    d_gene_def.parents = {"VD_ins"};
    d_gene_def.aliases = {"D"};
    model.sequence_types.definitions["D_gene"] = d_gene_def;

    SequenceTypeDefinition dj_ins_def;
    dj_ins_def.id = 3;
    dj_ins_def.parents = {"D_gene"};
    model.sequence_types.definitions["DJ_ins"] = dj_ins_def;

    SequenceTypeDefinition j_gene_def;
    j_gene_def.id = 4;
    j_gene_def.parents = {"DJ_ins"};
    j_gene_def.aliases = {"J"};
    model.sequence_types.definitions["J_gene"] = j_gene_def;

    return model;
}

/**
 * @brief Create minimal WeightsData for testing
 *
 * Creates weights for a single event (v_choice) with uniform probabilities.
 *
 * @return Minimal WeightsData structure
 */
inline WeightsData create_minimal_weights()
{
    WeightsData weights;
    weights.format_version = FORMAT_VERSION;

    EventWeights event_weights;
    event_weights.event_id = "v_choice";
    event_weights.dims = {3};
    event_weights.conditioning_events = {};
    event_weights.values = {0.333, 0.334, 0.333};
    event_weights.normalized = true;
    weights.events.push_back(event_weights);

    weights.error_rate.type = "SingleErrorRate";
    weights.error_rate.rate = 0.001;

    return weights;
}

/**
 * @brief Create complete WeightsData for VDJ model
 *
 * Creates weights for all events in a VDJ model:
 * - v_choice: [3] uniform
 * - d_gene: [2, 2] conditioned on j_choice
 * - j_choice: [2] uniform
 * - Deletion/insertion events with appropriate dimensions
 * - SingleErrorRate error model
 *
 * @return Complete WeightsData for VDJ model
 */
inline WeightsData create_vdj_weights()
{
    WeightsData weights;
    weights.format_version = FORMAT_VERSION;

    // V choice weights
    EventWeights v_choice;
    v_choice.event_id = "v_choice";
    v_choice.dims = {3};
    v_choice.conditioning_events = {};
    v_choice.values = {0.3, 0.4, 0.3};
    v_choice.normalized = true;
    weights.events.push_back(v_choice);

    // D gene weights (conditioned on J)
    EventWeights d_gene;
    d_gene.event_id = "d_gene";
    d_gene.dims = {2, 2};
    d_gene.conditioning_events = {"j_choice"};
    d_gene.values = {0.5, 0.5, 0.4, 0.6};
    d_gene.normalized = true;
    weights.events.push_back(d_gene);

    // J choice weights
    EventWeights j_choice;
    j_choice.event_id = "j_choice";
    j_choice.dims = {2};
    j_choice.conditioning_events = {};
    j_choice.values = {0.6, 0.4};
    j_choice.normalized = true;
    weights.events.push_back(j_choice);

    // Deletion weights (V 3' deletion: 0-10 nucleotides)
    EventWeights v_3_del;
    v_3_del.event_id = "v_3_del";
    v_3_del.dims = {11};
    v_3_del.conditioning_events = {};
    v_3_del.values = {0.05, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.03, 0.01, 0.005, 0.005};
    v_3_del.normalized = true;
    weights.events.push_back(v_3_del);

    // Error rate
    weights.error_rate.type = "SingleErrorRate";
    weights.error_rate.rate = 0.0015;
    weights.error_rate.learn_on = "seq";
    weights.error_rate.apply_on = "seq";

    return weights;
}

//==============================================================================
// Comparison Helpers
//==============================================================================

/**
 * @brief Compare two ModelData structures for equality
 *
 * Compares all fields including format_version, metadata, events, edges,
 * and sequence_types.
 *
 * @param a First model to compare
 * @param b Second model to compare
 * @return true if all fields are equal, false otherwise
 */
inline bool models_equal(const ModelData& a, const ModelData& b)
{
    if (a.format_version != b.format_version) return false;
    if (a.metadata.species != b.metadata.species) return false;
    if (a.metadata.chain != b.metadata.chain) return false;
    if (a.metadata.model_type != b.metadata.model_type) return false;

    if (a.events.size() != b.events.size()) return false;
    for (size_t i = 0; i < a.events.size(); ++i) {
        if (a.events[i].id != b.events[i].id) return false;
        if (a.events[i].type != b.events[i].type) return false;
        if (a.events[i].realizations.size() != b.events[i].realizations.size()) return false;
    }

    if (a.edges.size() != b.edges.size()) return false;
    for (size_t i = 0; i < a.edges.size(); ++i) {
        if (a.edges[i].parent != b.edges[i].parent) return false;
        if (a.edges[i].child != b.edges[i].child) return false;
    }

    return true;
}

/**
 * @brief Compare two WeightsData structures for equality
 *
 * Compares format_version, event weights, and error_rate.
 * Floating-point values are compared with tolerance.
 *
 * @param a First weights to compare
 * @param b Second weights to compare
 * @param tolerance Tolerance for floating-point comparison (default: 1e-9)
 * @return true if all fields are equal within tolerance, false otherwise
 */
inline bool weights_equal(
    const WeightsData& a,
    const WeightsData& b,
    double tolerance = 1e-9)
{
    if (a.format_version != b.format_version) return false;
    if (a.events.size() != b.events.size()) return false;

    for (size_t i = 0; i < a.events.size(); ++i) {
        if (a.events[i].event_id != b.events[i].event_id) return false;
        if (a.events[i].dims != b.events[i].dims) return false;
        if (a.events[i].values.size() != b.events[i].values.size()) return false;

        for (size_t j = 0; j < a.events[i].values.size(); ++j) {
            if (std::abs(a.events[i].values[j] - b.events[i].values[j]) > tolerance) {
                return false;
            }
        }
    }

    if (a.error_rate.type != b.error_rate.type) return false;
    if (std::abs(a.error_rate.rate - b.error_rate.rate) > tolerance) return false;

    return true;
}

} // namespace igor::modelio::test
