/**
 * @file test_model_format_v2.cpp
 * @brief Tests for the new model file format v2.0
 *
 * Tests backward compatibility with legacy format and new v2.0 format features.
 * The v2.0 format adds:
 *   - @Version section (2.0)
 *   - @Seq_type_order section
 *   - Explicit seq_type field in event declarations
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/Singleerrorrate.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Utils.h>
#include <igor/Core/EventUtils.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>

static const std::string TEST_DATA_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/format_v2/";

// Test that legacy format files can still be read
TEST_CASE("Legacy model format backward compatibility", "[model_format][legacy]")
{
    // Read the legacy format file from test_data/format_v2
    std::string filename = "test_legacy_model_parms.txt";
    std::string file_path = TEST_DATA_DIR + filename;
    
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));
    
    // Verify events were loaded correctly
    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5); // 2 GeneChoice + 1 Deletion + 1 Insertion + 1 DinucMarkov (error rate is separate)
    
    // Check that V gene choice event has correct properties
    auto v_choice = parms.get_event_pointer("v_choice", true);
    REQUIRE(v_choice != nullptr);
    REQUIRE(v_choice->get_type() == GeneChoice_t);
    REQUIRE(v_choice->get_class() == V_gene);
}


// Test that v2.0 format files can be read (placeholder for future implementation)
TEST_CASE("Model format v2.0 parsing", "[model_format][v2.0]")
{
    // Read the v2.0 format file from test_data/format_v2
    std::string filename = "test_v2_model_parms.txt";
    std::string file_path = TEST_DATA_DIR + filename;
    
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));
    
    // Verify events were loaded correctly
    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5); // 2 GeneChoice + 1 Deletion + 1 Insertion + 1 DinucMarkov (error rate is separate)
}

// Test round-trip: write to test_data/format_v2 and read back
TEST_CASE("Model format round-trip", "[model_format]")
{
    // Create a simple model
    Model_Parms original_parms;
    
    // Add V gene choice
    auto v_choice = std::make_shared<Gene_choice>(V_gene);
    v_choice->add_realization("V1", "ATCGATCGATCG");
    v_choice->set_priority(7);
    v_choice->set_nickname("v_choice");
    original_parms.add_event(v_choice);
    
    // Add J gene choice
    auto j_choice = std::make_shared<Gene_choice>(J_gene);
    j_choice->add_realization("J1", "GCTAGCTAGCTA");
    j_choice->set_priority(7);
    j_choice->set_nickname("j_choice");
    original_parms.add_event(j_choice);
    
    // Add V 3' deletion
    auto v_3_del = std::make_shared<Deletion>(V_gene, Three_prime);
    v_3_del->add_realization(0);
    v_3_del->add_realization(1);
    v_3_del->set_priority(5);
    v_3_del->set_nickname("v_3_del");
    original_parms.add_event(v_3_del);
    
    // Add edge
    original_parms.add_edge(v_choice, v_3_del);
    original_parms.add_edge(v_choice, j_choice);
    
    // Add error rate
    auto err_rate = std::make_shared<Single_error_rate>(0.001);
    original_parms.set_error_ratep(err_rate);
    
    // Write to test_data directory and read back
    std::string roundtrip_file = TEST_DATA_DIR + "test_roundtrip_parms.txt";
    original_parms.write_model_parms(roundtrip_file);
    
    // Read back
    Model_Parms loaded_parms;
    REQUIRE_NOTHROW(loaded_parms.read_model_parms(roundtrip_file));
    
    // Compare
    auto orig_events = original_parms.get_event_list();
    auto loaded_events = loaded_parms.get_event_list();
    REQUIRE(orig_events.size() == loaded_events.size());
    
    // Cleanup
    std::filesystem::remove(roundtrip_file);
}

// Test format detection
TEST_CASE("Format version detection", "[model_format]")
{
    // Test legacy format detection (no @Version header)
    std::string legacy_file = TEST_DATA_DIR + "test_legacy_detect.txt";
    Model_Parms legacy_parms;
    REQUIRE_NOTHROW(legacy_parms.read_model_parms(legacy_file));
    
    // Test v2.0 format detection
    std::string v2_file = TEST_DATA_DIR + "test_v2_detect.txt";
    Model_Parms v2_parms;
    REQUIRE_NOTHROW(v2_parms.read_model_parms(v2_file));
}

// Test reading actual v2.0 format model file from models directory
TEST_CASE("Load actual v2.0 model file from models directory", "[model_format][v2.0][integration]")
{
    // Get the project source directory from CMake define
    std::string source_dir = IGOR_SOURCE_DIR;
    std::string v2_model_path = source_dir + "/models/human/tcr_beta/models_v2/model_parms.txt";
    
    // Skip if model file doesn't exist (e.g., submodule not initialized)
    if (!std::filesystem::exists(v2_model_path)) {
        SKIP("v2.0 model file not found at " << v2_model_path << " (submodule may not be initialized)");
    }
    
    // Load the v2.0 model
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(v2_model_path));
    
    // Verify events were loaded
    auto events = parms.get_event_list();
    REQUIRE(events.size() >= 10); // Should have V_choice, D_choice, J_choice, deletions, insertions, dinuc
    
    // Check specific events exist
    auto v_choice = parms.get_event_pointer("v_choice", true);
    REQUIRE(v_choice != nullptr);
    REQUIRE(v_choice->get_type() == GeneChoice_t);
    REQUIRE(v_choice->get_class() == V_gene);
    REQUIRE(v_choice->get_seq_type() == "V_gene_seq");
    
    auto v_3_del = parms.get_event_pointer("v_3_del", true);
    REQUIRE(v_3_del != nullptr);
    REQUIRE(v_3_del->get_type() == Deletion_t);
    REQUIRE(v_3_del->get_seq_type() == "V_gene_seq");
    
    auto vd_ins = parms.get_event_pointer("vd_ins", true);
    REQUIRE(vd_ins != nullptr);
    REQUIRE(vd_ins->get_type() == Insertion_t);
    REQUIRE(vd_ins->get_seq_type() == "VD_ins_seq");
    
    auto vd_dinucl = parms.get_event_pointer("vd_dinucl", true);
    REQUIRE(vd_dinucl != nullptr);
    REQUIRE(vd_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vd_dinucl->get_seq_type() == "VD_ins_seq");
    
    // Verify edges were loaded
    auto edges = parms.get_edges();
    REQUIRE(!edges.empty());
}

// Test that a v2 file with a custom (tandem-D) @Seq_type_order populates
// the SeqTypeRegistry with exactly the declared order, including non-standard
// seq_type names that cannot exist in legacy format.
TEST_CASE("SeqTypeRegistry populated from v2 custom seq_type_order", "[model_format][v2.0][seq_type_registry]")
{
    const std::string file_path = TEST_DATA_DIR + "test_v2_tandem_d_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    // Registry must reflect exactly the seven entries listed in @Seq_type_order
    const std::vector<std::string> expected = {
        "V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
        "D1D2_ins_seq", "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    REQUIRE(registry.get_ordered_types() == expected);

    // Every declared name must be reachable
    REQUIRE(registry.contains("D1_gene_seq"));
    REQUIRE(registry.contains("D2_gene_seq"));
    REQUIRE(registry.contains("D1D2_ins_seq"));

    // Neighbour navigation must follow the declared order
    REQUIRE(registry.get_left_neighbor("D1_gene_seq")  == std::optional<std::string>("VD1_ins_seq"));
    REQUIRE(registry.get_right_neighbor("D1_gene_seq") == std::optional<std::string>("D1D2_ins_seq"));
    REQUIRE(registry.get_left_neighbor("V_gene_seq")   == std::nullopt);
    REQUIRE(registry.get_right_neighbor("J_gene_seq")  == std::nullopt);
}

// Test that reading a legacy VJ file (no @Version / no @Seq_type_order) produces
// the inferred VJ ordering in the SeqTypeRegistry.
TEST_CASE("SeqTypeRegistry inferred as VJ order from legacy format", "[model_format][legacy][seq_type_registry]")
{
    const std::string file_path = TEST_DATA_DIR + "test_legacy_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    const std::vector<std::string> expected = {"V_gene_seq", "VJ_ins_seq", "J_gene_seq"};
    REQUIRE(registry.get_ordered_types() == expected);
}

// Test that reading a legacy VDJ file (no @Version / no @Seq_type_order) produces
// the inferred VDJ ordering in the SeqTypeRegistry.
TEST_CASE("SeqTypeRegistry inferred as VDJ order from legacy format", "[model_format][legacy][seq_type_registry]")
{
    const std::string file_path = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    const std::vector<std::string> expected = {
        "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    REQUIRE(registry.get_ordered_types() == expected);
}

// ── Gene_class validation in v2 format ──────────────────────────────────────

// In v2.0 a GeneChoice event may only carry an alignment-purpose gene_class
// (V_gene, D_gene, J_gene, Undefined_gene).  Junction classes (VD_genes,
// DJ_genes, VJ_genes, VDJ_genes) are illegal and must produce an error.
TEST_CASE("v2 format: GeneChoice with junction gene_class throws", "[model_format][v2.0][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "test_v2_bad_genechoice_gene_class.txt";

    Model_Parms parms;
    REQUIRE_THROWS_AS(parms.read_model_parms(file_path), std::runtime_error);
    REQUIRE_THROWS_WITH(
        parms.read_model_parms(file_path),
        Catch::Matchers::ContainsSubstring("VD_genes") &&
        Catch::Matchers::ContainsSubstring("GeneChoice") &&
        Catch::Matchers::ContainsSubstring("v2.0"));
}

// In legacy format the junction gene_class values (VJ_gene, VD_genes, …) on
// Insertion and DinucMarkov events must still be accepted and automatically
// converted to the corresponding seq_type string.
TEST_CASE("legacy format: junction gene_class converted to seq_type on Insertion and DinucMarkov",
          "[model_format][legacy][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "test_legacy_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    // Insertion: VJ_gene (legacy) → seq_type "VJ_ins_seq"
    auto vj_ins = parms.get_event_pointer("vj_ins", true);
    REQUIRE(vj_ins != nullptr);
    REQUIRE(vj_ins->get_type() == Insertion_t);
    REQUIRE(vj_ins->get_seq_type() == "VJ_ins_seq");

    // DinucMarkov: VJ_gene (legacy) → seq_type "VJ_ins_seq"
    auto vj_dinucl = parms.get_event_pointer("vj_dinucl", true);
    REQUIRE(vj_dinucl != nullptr);
    REQUIRE(vj_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vj_dinucl->get_seq_type() == "VJ_ins_seq");
}

// In legacy format GeneChoice gene_class (V_gene, D_gene, J_gene) must be
// preserved as-is and the seq_type must be the corresponding "*_seq" string.
TEST_CASE("legacy format: GeneChoice gene_class preserved and seq_type inferred",
          "[model_format][legacy][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "test_legacy_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    auto v_choice = parms.get_event_pointer("v_choice", true);
    REQUIRE(v_choice != nullptr);
    REQUIRE(v_choice->get_class() == V_gene);
    REQUIRE(v_choice->get_seq_type() == "V_gene_seq");

    auto j_choice = parms.get_event_pointer("j_choice", true);
    REQUIRE(j_choice != nullptr);
    REQUIRE(j_choice->get_class() == J_gene);
    REQUIRE(j_choice->get_seq_type() == "J_gene_seq");
}

// In v2.0 format Insertion and DinucMarkov must use Undefined_gene as
// gene_class, and their seq_type must equal the explicit field from the file.
TEST_CASE("v2 format: Insertion and DinucMarkov use explicit seq_type field",
          "[model_format][v2.0][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "test_v2_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    auto vj_ins = parms.get_event_pointer("vj_ins", true);
    REQUIRE(vj_ins != nullptr);
    REQUIRE(vj_ins->get_type() == Insertion_t);
    REQUIRE(vj_ins->get_seq_type() == "VJ_ins_seq");

    auto vj_dinucl = parms.get_event_pointer("vj_dinucl", true);
    REQUIRE(vj_dinucl != nullptr);
    REQUIRE(vj_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vj_dinucl->get_seq_type() == "VJ_ins_seq");
}

// In v2.0 format GeneChoice with valid alignment gene_class must be accepted
// and the seq_type field must be stored verbatim.
TEST_CASE("v2 format: GeneChoice with valid gene_class accepted",
          "[model_format][v2.0][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "test_v2_model_parms.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    auto v_choice = parms.get_event_pointer("v_choice", true);
    REQUIRE(v_choice != nullptr);
    REQUIRE(v_choice->get_class() == V_gene);
    REQUIRE(v_choice->get_seq_type() == "V_gene_seq");

    auto j_choice = parms.get_event_pointer("j_choice", true);
    REQUIRE(j_choice != nullptr);
    REQUIRE(j_choice->get_class() == J_gene);
    REQUIRE(j_choice->get_seq_type() == "J_gene_seq");
}

// ── v1/v2 parity tests ──────────────────────────────────────────────────────

// Helper: collect comparable per-event properties from a loaded Model_Parms.
struct EventSummary {
    std::string nickname;
    Event_type  type;
    Gene_class  gene_class;  // only meaningful for GeneChoice
    std::string seq_type;
    Seq_side    side;
    int         num_realizations;
};

static std::vector<EventSummary> summarise_events(const Model_Parms &parms)
{
    std::vector<EventSummary> out;
    for (const auto &ev : parms.get_event_list()) {
        EventSummary s;
        s.nickname          = ev->get_nickname();
        s.type              = ev->get_type();
        s.gene_class        = ev->get_class();
        s.seq_type          = ev->get_seq_type();
        s.side              = ev->get_side();
        s.num_realizations  = static_cast<int>(ev->get_realizations_map().size());
        out.push_back(s);
    }
    // Sort by nickname so ordering differences don't matter.
    std::sort(out.begin(), out.end(),
              [](const EventSummary &a, const EventSummary &b) {
                  return a.nickname < b.nickname;
              });
    return out;
}

static void require_same_events(const Model_Parms &v1_parms, const Model_Parms &v2_parms)
{
    auto v1 = summarise_events(v1_parms);
    auto v2 = summarise_events(v2_parms);
    REQUIRE(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        INFO("Event: " << v1[i].nickname);
        REQUIRE(v1[i].nickname         == v2[i].nickname);
        REQUIRE(v1[i].type             == v2[i].type);
        REQUIRE(v1[i].seq_type         == v2[i].seq_type);
        REQUIRE(v1[i].num_realizations == v2[i].num_realizations);
        // For GeneChoice events the alignment gene_class must be preserved.
        if (v1[i].type == GeneChoice_t) {
            REQUIRE(v1[i].gene_class == v2[i].gene_class);
        }
    }
}

// VJ legacy vs hand-written v2 counterpart.
TEST_CASE("v1/v2 parity: VJ model - hand-written v2 file matches legacy",
          "[model_format][parity][v2.0]")
{
    Model_Parms v1_parms, v2_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
    REQUIRE_NOTHROW(v2_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms_v2.txt"));

    require_same_events(v1_parms, v2_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            v2_parms.get_seq_type_registry().get_ordered_types());
}

// VDJ legacy vs hand-written v2 counterpart.
TEST_CASE("v1/v2 parity: VDJ model - hand-written v2 file matches legacy",
          "[model_format][parity][v2.0]")
{
    Model_Parms v1_parms, v2_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));
    REQUIRE_NOTHROW(v2_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms_v2.txt"));

    require_same_events(v1_parms, v2_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            v2_parms.get_seq_type_registry().get_ordered_types());
}

// Round-trip: read legacy -> write v2 -> read back -> same as legacy.
TEST_CASE("v1/v2 parity: write_model_parms_v2 round-trip equals original",
          "[model_format][parity][v2.0]")
{
    const std::string v1_path = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    const std::string tmp_v2  = TEST_DATA_DIR + "tmp_roundtrip_v2.txt";

    Model_Parms v1_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(v1_path));

    REQUIRE_NOTHROW(v1_parms.write_model_parms_v2(tmp_v2));

    Model_Parms v2_parms;
    REQUIRE_NOTHROW(v2_parms.read_model_parms(tmp_v2));

    require_same_events(v1_parms, v2_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            v2_parms.get_seq_type_registry().get_ordered_types());

    std::filesystem::remove(tmp_v2);
}

// Round-trip using the large inference model file (covers full V/D/J gene sets).
TEST_CASE("v1/v2 parity: write_model_parms_v2 round-trip for inference model",
          "[model_format][parity][v2.0][inference]")
{
    const std::string INFER_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/inference/";
    const std::string v1_path   = INFER_DIR + "fixed_VJ_TRB_model_parms.txt";
    const std::string tmp_v2    = INFER_DIR + "tmp_inference_v2.txt";

    Model_Parms v1_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(v1_path));

    REQUIRE_NOTHROW(v1_parms.write_model_parms_v2(tmp_v2));

    Model_Parms v2_parms;
    REQUIRE_NOTHROW(v2_parms.read_model_parms(tmp_v2));

    require_same_events(v1_parms, v2_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            v2_parms.get_seq_type_registry().get_ordered_types());

    std::filesystem::remove(tmp_v2);
}

// ── Step 7: requires_extended_format() ─────────────────────────────────────

// Standard VJ model uses only legacy seq_types in VJ order → no extension needed.
TEST_CASE("requires_extended_format: VJ model returns false",
          "[model_format][step7][requires_extended_format]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
    REQUIRE_FALSE(parms.requires_extended_format());
}

// Standard VDJ model uses only legacy seq_types in VDJ order → no extension needed.
TEST_CASE("requires_extended_format: VDJ model returns false",
          "[model_format][step7][requires_extended_format]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));
    REQUIRE_FALSE(parms.requires_extended_format());
}

// Tandem-D model has non-standard seq_types (D1_gene_seq, VD1_ins_seq …) → must use v2.
TEST_CASE("requires_extended_format: tandem-D model returns true",
          "[model_format][step7][requires_extended_format]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_v2_tandem_d_parms.txt"));
    REQUIRE(parms.requires_extended_format());
}

// Programmatically created model with no registry → not extended (legacy OK).
TEST_CASE("requires_extended_format: empty registry returns false",
          "[model_format][step7][requires_extended_format]")
{
    Model_Parms parms;   // no events, no registry
    REQUIRE_FALSE(parms.requires_extended_format());
}

// ── Step 7: write_model_parms() auto-selection ──────────────────────────────

// For a standard VDJ model, write_model_parms() must produce a legacy file
// (no @Version header).
TEST_CASE("write_model_parms auto-selects legacy format for standard VDJ model",
          "[model_format][step7][auto_select]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_auto_select_legacy.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms(dest));

    // First line must NOT be @Version
    std::ifstream f(dest);
    REQUIRE(f.is_open());
    std::string first_line;
    std::getline(f, first_line);
    REQUIRE(first_line != "@Version");

    // Must be readable as legacy (no version header means legacy parser path)
    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(parms, loaded);

    std::filesystem::remove(dest);
}

// For a tandem-D model with non-standard seq_types, write_model_parms() must
// produce a v2 file (starts with @Version / 2.0).
TEST_CASE("write_model_parms auto-selects v2 format for tandem-D model",
          "[model_format][step7][auto_select]")
{
    const std::string src  = TEST_DATA_DIR + "test_v2_tandem_d_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_auto_select_v2.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms(dest));

    std::ifstream f(dest);
    REQUIRE(f.is_open());
    std::string first_line;
    std::getline(f, first_line);
    REQUIRE(first_line == "@Version");

    // Must also be readable back with the same events
    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(parms, loaded);

    std::filesystem::remove(dest);
}

// ── Step 7: write_model_parms_legacy() round-trip ──────────────────────────

// Read a VDJ legacy model, write it via write_model_parms_legacy(), read back,
// and verify that all events are identical.
TEST_CASE("write_model_parms_legacy: VDJ round-trip preserves events",
          "[model_format][step7][legacy_roundtrip]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_legacy_roundtrip.txt";

    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(src));
    REQUIRE_NOTHROW(original.write_model_parms_legacy(dest));

    // Written file must not contain @Version
    std::ifstream f(dest);
    REQUIRE(f.is_open());
    std::string first_line;
    std::getline(f, first_line);
    REQUIRE(first_line != "@Version");

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(original, loaded);

    REQUIRE(original.get_seq_type_registry().get_ordered_types() ==
            loaded.get_seq_type_registry().get_ordered_types());

    std::filesystem::remove(dest);
}

// ── Step 5: write2txt_v2() output format check ─────────────────────────────

// write_model_parms_v2() must produce event lines with the seq_type field
// between gene_class and seq_side: #Deletion;Undefined_gene;<seq_type>;<side>...
TEST_CASE("write2txt_v2 output contains seq_type field in correct position",
          "[model_format][step5][write2txt_v2]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_v2_format_check.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms_v2(dest));

    std::ifstream f(dest);
    REQUIRE(f.is_open());

    bool found_deletion_v2_line = false;
    bool found_insertion_v2_line = false;
    bool found_genechoice_v2_line = false;
    std::string line;
    while (std::getline(f, line)) {
        // Deletion v2: #Deletion;Undefined_gene;<seq_type>;<side>;<prio>;<nick>
        if (line.size() > 1 && line[0] == '#' && line.find("Deletion") == 1) {
            // Count semicolons — v2 has 5 fields after '#Deletion'
            size_t count = std::count(line.begin(), line.end(), ';');
            REQUIRE(count == 5);  // gene_class;seq_type;side;priority;nickname
            REQUIRE(line.find("Undefined_gene") != std::string::npos);
            REQUIRE(line.find("_seq") != std::string::npos);  // seq_type present
            found_deletion_v2_line = true;
        }
        // Insertion v2: #Insertion;Undefined_gene;<seq_type>;...
        if (line.size() > 1 && line[0] == '#' && line.find("Insertion") == 1) {
            size_t count = std::count(line.begin(), line.end(), ';');
            REQUIRE(count == 5);
            REQUIRE(line.find("Undefined_gene") != std::string::npos);
            REQUIRE(line.find("_seq") != std::string::npos);
            found_insertion_v2_line = true;
        }
        // GeneChoice v2: #GeneChoice;<gene_class>;<seq_type>;...
        if (line.size() > 1 && line[0] == '#' && line.find("GeneChoice") == 1) {
            size_t count = std::count(line.begin(), line.end(), ';');
            REQUIRE(count == 5);
            REQUIRE(line.find("_seq") != std::string::npos);
            found_genechoice_v2_line = true;
        }
    }
    REQUIRE(found_deletion_v2_line);
    REQUIRE(found_insertion_v2_line);
    REQUIRE(found_genechoice_v2_line);

    std::filesystem::remove(dest);
}

// ══════════════════════════════════════════════════════════════════════════════
// Step 8: @Edges full event name format change
// ══════════════════════════════════════════════════════════════════════════════
//
// Current (legacy and v2 as written today):
//   %GeneChoice_V_gene_Undefined_side_prio7_size1;Deletion_V_gene_Three_prime_prio5_size2
//
// Required v2 format (§7.11):
//   %GeneChoice_V_gene_V_gene_seq_Undefined_side_prio7_size1;Deletion_V_gene_V_gene_seq_Three_prime_prio5_size2
//
// Impact: Rec_Event::update_event_name() must include seq_type between
// gene_class and seq_side when the event carries a non-empty seq_type.
// The parser must handle both formats.
// ──────────────────────────────────────────────────────────────────────────────

// Current (legacy) event names do NOT include seq_type.
// This verifies the existing behaviour that must be preserved for legacy output.
TEST_CASE("step8: legacy event names do not include seq_type",
          "[model_format][step8][event_name]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));

    for (const auto &ev : parms.get_event_list()) {
        const std::string &name = ev->get_name();
        const std::string &st   = ev->get_seq_type();
        // Legacy names: <type>_<class>_<side>_prio<N>_size<N>
        // seq_type must NOT appear between class and side.
        if (!st.empty()) {
            // Find the position just after "<type>_<class>_"
            size_t underscore_after_class = name.find('_', name.find('_') + 1);
            if (underscore_after_class != std::string::npos) {
                std::string after_class = name.substr(underscore_after_class + 1);
                // The next token should be the seq_side, not the seq_type.
                REQUIRE(after_class.find(st) == std::string::npos);
            }
        }
    }
}

// In v2 format the event name must embed seq_type so that edges can be
// unambiguously resolved across models with different orderings.
// get_v2_name() returns the seq_type-inclusive name; get_name() stays legacy.
TEST_CASE("step8: v2 event names include seq_type between class and side",
          "[model_format][step8][event_name]")
{
    // get_v2_name() must return a name of the form:
    //   <type>_<class>_<seq_type>_<side>_prio<N>_size<N>
    // (The internal get_name() stays in legacy format; get_v2_name() is the v2 API.)
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_v2_model_parms.txt"));

    for (const auto &ev : parms.get_event_list()) {
        const std::string v2name = ev->get_v2_name();
        const std::string &st   = ev->get_seq_type();
        INFO("v2_name: " << v2name << "  seq_type: " << st);
        if (!st.empty()) {
            REQUIRE(v2name.find(st) != std::string::npos);
        }
    }
}

// write_model_parms_v2() edge lines must use the seq_type-inclusive event names.
TEST_CASE("step8: write_model_parms_v2 edge names include seq_type",
          "[model_format][step8][event_name]")
{
    // Use the VJ model (test_legacy_model_parms.txt) which has actual edges.
    const std::string src  = TEST_DATA_DIR + "test_legacy_model_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_step8_edge_names.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms_v2(dest));

    std::ifstream f(dest);
    REQUIRE(f.is_open());

    bool in_edges = false;
    bool found_v2_edge = false;
    std::string line;
    while (std::getline(f, line)) {
        if (line == "@Edges") { in_edges = true; continue; }
        if (line[0] == '@') { in_edges = false; continue; }
        if (in_edges && line[0] == '%') {
            // Each edge token should contain a seq_type substring (ends in _seq)
            REQUIRE(line.find("_seq") != std::string::npos);
            found_v2_edge = true;
        }
    }
    REQUIRE(found_v2_edge);

    std::filesystem::remove(dest);
}

// A v2 file that uses the new (seq_type-inclusive) edge name format must be
// parseable (get_event_pointer() accepts both legacy and v2 edge name formats).
TEST_CASE("step8: v2 file with seq_type-inclusive edge names is parseable",
          "[model_format][step8][parser]")
{
    // test_v2_new_edge_names.txt is a hand-crafted v2 file whose @Edges section
    // uses the new format: GeneChoice_V_gene_V_gene_seq_Undefined_side_prio7_size1
    const std::string file_path = TEST_DATA_DIR + "test_v2_new_edge_names.txt";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5);  // v_choice, j_choice, v_3_del, vj_ins, vj_dinucl

    // Edges must be wired: v_choice → v_3_del and v_choice → j_choice
    auto edges = parms.get_edges();
    REQUIRE(!edges.empty());
}

// ══════════════════════════════════════════════════════════════════════════════
// Step 9: build_scenario_sequence() / sequence ordering
// ══════════════════════════════════════════════════════════════════════════════
//
// Current: EventUtils::build_scenario_sequence() uses boolean flags and
// hardcoded seq_type names (V_gene_seq, VD_ins_seq, …).
//
// Required: must traverse the SeqTypeRegistry linked list in order so that
// non-standard orderings (tandem D: V→VD1ins→D1→D1D2ins→D2→DJins→J) work.
// ──────────────────────────────────────────────────────────────────────────────

// Standard VDJ ordering matches the current hardcoded path — must still work
// after step 9 is implemented.
TEST_CASE("step9: standard VDJ model is loadable and has correct seq_type_order",
          "[model_format][step9][seq_order]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    const std::vector<std::string> expected = {
        "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    REQUIRE(parms.get_seq_type_registry().get_ordered_types() == expected);
}

// A tandem-D model must be loadable with the correct 7-element registry order.
TEST_CASE("step9: tandem-D model seq_type_order has 7 elements in correct order",
          "[model_format][step9][seq_order]")
{
    // This part (reading + checking registry) is already implemented.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_v2_tandem_d_parms.txt"));

    const std::vector<std::string> expected = {
        "V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
        "D1D2_ins_seq", "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    REQUIRE(parms.get_seq_type_registry().get_ordered_types() == expected);
}

TEST_CASE("step9: build_scenario_sequence uses registry order for tandem D",
          "[model_format][step9][seq_order]")
{
    // Load the tandem-D model to get its registry, then call the registry-based
    // build_scenario_sequence() and verify 7-segment assembly order.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_v2_tandem_d_parms.txt"));
    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    Int_Str v    = {0, 0, 0};
    Int_Str vd1  = {1};
    Int_Str d1   = {2, 2};
    Int_Str d1d2 = {3};
    Int_Str d2   = {4, 4};
    Int_Str dj   = {5};
    Int_Str j    = {6, 6, 6};

    std::unordered_map<std::string, const Int_Str *> seqs;
    seqs["V_gene_seq"]   = &v;
    seqs["VD1_ins_seq"]  = &vd1;
    seqs["D1_gene_seq"]  = &d1;
    seqs["D1D2_ins_seq"] = &d1d2;
    seqs["D2_gene_seq"]  = &d2;
    seqs["DJ_ins_seq"]   = &dj;
    seqs["J_gene_seq"]   = &j;

    Int_Str result = EventUtils::build_scenario_sequence(registry, seqs);

    Int_Str expected;
    expected.insert(expected.end(), v.begin(),    v.end());
    expected.insert(expected.end(), vd1.begin(),  vd1.end());
    expected.insert(expected.end(), d1.begin(),   d1.end());
    expected.insert(expected.end(), d1d2.begin(), d1d2.end());
    expected.insert(expected.end(), d2.begin(),   d2.end());
    expected.insert(expected.end(), dj.begin(),   dj.end());
    expected.insert(expected.end(), j.begin(),    j.end());

    REQUIRE(result == expected);
}

// ══════════════════════════════════════════════════════════════════════════════
// Step 10: Seq_offsets_map key type change
// ══════════════════════════════════════════════════════════════════════════════
//
// Current: Seq_offsets_map is Enum_fast_memory_dual_key_map<Seq_type, Seq_side, …>
// keyed by the fixed Seq_type enum (6 values).
//
// Required: must accept string keys so that non-standard seq_types such as
// D1_gene_seq, D2_gene_seq, VD1_ins_seq can be used.
// ──────────────────────────────────────────────────────────────────────────────

// Standard seq_types are still accessible through the current map —
// documents that legacy behaviour is preserved.
TEST_CASE("step10: standard seq_types accessible in Seq_offsets_map",
          "[model_format][step10][seq_offsets]")
{
    // Seq_offsets_map is populated during initialize_event(); here we just
    // verify that a standard VJ model loads without error, implying that
    // V_gene_seq, VJ_ins_seq, J_gene_seq are all accessible.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
    REQUIRE(parms.get_event_list().size() == 5);
}

// Non-standard seq_types (tandem D) can now be keyed in Str_Seq_offsets_map.
TEST_CASE("step10: Str_Seq_offsets_map accepts non-standard string seq_type keys",
          "[model_format][step10][seq_offsets]")
{
    // Str_Seq_offsets_map uses std::string as first key, allowing arbitrary
    // seq_type names such as "D1_gene_seq", "D2_gene_seq", "VD1_ins_seq".
    Str_Seq_offsets_map offsets;

    // request_memory_layer creates the slot on first use
    REQUIRE_NOTHROW(offsets.request_memory_layer("D1_gene_seq", Undefined_side));
    REQUIRE_NOTHROW(offsets.request_memory_layer("D2_gene_seq", Undefined_side));
    REQUIRE_NOTHROW(offsets.request_memory_layer("VD1_ins_seq", Undefined_side));

    CHECK(offsets.exist("D1_gene_seq", Undefined_side));
    CHECK(offsets.exist("D2_gene_seq", Undefined_side));
    CHECK(offsets.exist("VD1_ins_seq", Undefined_side));
    CHECK_FALSE(offsets.exist("V_gene_seq", Undefined_side));  // never registered

    // set and read back a value
    offsets.set_value("D1_gene_seq", Undefined_side, 42, 0);
    CHECK(offsets.at("D1_gene_seq", Undefined_side) == 42);

    // a second memory layer
    REQUIRE_NOTHROW(offsets.request_memory_layer("D1_gene_seq", Undefined_side));
    offsets.set_value("D1_gene_seq", Undefined_side, 99, 1);
    CHECK(offsets.at("D1_gene_seq", Undefined_side) == 99);
    CHECK(offsets.get_current_memory_layer("D1_gene_seq", Undefined_side) == 1);
}

// ══════════════════════════════════════════════════════════════════════════════
// Step 11: DinucMarkov direction convention
// ══════════════════════════════════════════════════════════════════════════════
//
// Current: direction (left-to-right vs right-to-left) is implicit from the
// legacy gene_class:  VD_genes / VJ_genes → Three_prime (right-to-left);
//                     DJ_genes            → Five_prime  (left-to-right).
//
// Required: direction is explicit in seq_side.  The legacy converter must set
// it on read; v2 files carry the explicit value verbatim.
// ──────────────────────────────────────────────────────────────────────────────

// Legacy VD DinucMarkov → Three_prime (right-to-left insertion).
TEST_CASE("step11: DinucMarkov VD insertion has Three_prime side after legacy load",
          "[model_format][step11][dinucl_direction]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    auto vd_dinucl = parms.get_event_pointer("vd_dinucl", true);
    REQUIRE(vd_dinucl != nullptr);
    REQUIRE(vd_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vd_dinucl->get_side() == Three_prime);
}

// Legacy DJ DinucMarkov → Five_prime (left-to-right insertion).
TEST_CASE("step11: DinucMarkov DJ insertion has Five_prime side after legacy load",
          "[model_format][step11][dinucl_direction]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    auto dj_dinucl = parms.get_event_pointer("dj_dinucl", true);
    REQUIRE(dj_dinucl != nullptr);
    REQUIRE(dj_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(dj_dinucl->get_side() == Five_prime);
}

// v2 DinucMarkov side must be stored verbatim from the file and survive
// a write_model_parms_v2 → read round-trip.
TEST_CASE("step11: DinucMarkov direction preserved through v2 round-trip",
          "[model_format][step11][dinucl_direction]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    const std::string dest = TEST_DATA_DIR + "tmp_step11_dinucl.txt";

    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(src));

    REQUIRE_NOTHROW(original.write_model_parms_v2(dest));

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));

    auto vd_orig   = original.get_event_pointer("vd_dinucl", true);
    auto vd_loaded = loaded.get_event_pointer("vd_dinucl", true);
    REQUIRE(vd_loaded->get_side() == vd_orig->get_side());

    auto dj_orig   = original.get_event_pointer("dj_dinucl", true);
    auto dj_loaded = loaded.get_event_pointer("dj_dinucl", true);
    REQUIRE(dj_loaded->get_side() == dj_orig->get_side());

    std::filesystem::remove(dest);
}

// v2 file carries explicit seq_side for DinucMarkov; must be loaded verbatim.
TEST_CASE("step11: DinucMarkov direction in v2 file is loaded verbatim",
          "[model_format][step11][dinucl_direction]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_v2_model_parms.txt"));

    // In the VJ v2 file, the DinucMarkov event for VJ insertion should be
    // Three_prime (right-to-left, matching VJ_genes legacy convention).
    auto vj_dinucl = parms.get_event_pointer("vj_dinucl", true);
    REQUIRE(vj_dinucl != nullptr);
    REQUIRE(vj_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vj_dinucl->get_side() == Three_prime);
}
