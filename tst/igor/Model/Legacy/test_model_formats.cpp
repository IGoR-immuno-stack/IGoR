/**
 * @file test_model_formats.cpp
 * @brief The model file formats: v1, and JSON for a model v1 cannot hold.
 *
 * JSON replaced the v2 text format (plan T2), which no release read. Each v2 case has its JSON
 * counterpart here. What JSON adds over v1:
 *   - the segment order (`seq_type_order`);
 *   - each event's seq type, so any segment can be named;
 *   - edges as each event's parents, named by nickname rather than by a generated name.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "ScratchDir.h"

#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>
#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Model/Legacy/Genechoice.h>
#include <igor/Model/Legacy/Deletion.h>
#include <igor/Model/Legacy/Insertion.h>
#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Model/Legacy/Singleerrorrate.h>
#include <igor/Core/Legacy/SeqTypeRegistry.h>
#include <igor/Core/Legacy/Utils.h>
#include <igor/Model/Legacy/EventUtils.h>
#include <igor/Model/Legacy/ModelJson.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>

using namespace igor::core::legacy;
using namespace igor::alignment::legacy;
using namespace igor::model::legacy;

// Tests tagged [tandem_d] describe the intended tandem-D behaviour, the acceptance criteria
// for milestone 1 of the Rec_Event refactoring (see docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md).
// They were [!mayfail] until T1: no event could be built on a seq type outside the six legacy
// names, so a tandem-D model file could not be loaded.

namespace {

nlohmann::json read_json(const std::string &path)
{
    std::ifstream in(path);
    return nlohmann::json::parse(in);
}

std::string first_line_of(const std::string &path)
{
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    return line;
}

/// Each event's children, by nickname.
std::map<std::string, std::set<std::string>> children_by_nickname(const Model_Parms &parms)
{
    std::map<std::string, std::set<std::string>> graph;
    for (const auto &event : parms.get_event_list()) {
        auto &children = graph[event->get_nickname()];
        for (const auto &child : parms.get_children(event)) {
            children.insert(child->get_nickname());
        }
    }
    return graph;
}

}  // namespace
static const std::string TEST_DATA_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/model_formats/";

// Test that legacy format files can still be read
TEST_CASE("Legacy model format backward compatibility", "[model_format][legacy]")
{
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


// A JSON document can be read.
TEST_CASE("Model format JSON parsing", "[model_format][json]")
{
    std::string filename = "vj_model.json";
    std::string file_path = TEST_DATA_DIR + filename;
    
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));
    
    // Verify events were loaded correctly
    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5); // 2 GeneChoice + 1 Deletion + 1 Insertion + 1 DinucMarkov (error rate is separate)
}

// Test round-trip: write a model built in code and read it back
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
    auto v_3_del = std::make_shared<Deletion>("V_gene_seq", Three_prime);
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
    
    // Write to a scratch directory and read back
    ScratchDir scratch;
    std::string roundtrip_file = scratch.file("test_roundtrip_parms.txt");
    original_parms.write_model_parms(roundtrip_file);
    
    // Read back
    Model_Parms loaded_parms;
    REQUIRE_NOTHROW(loaded_parms.read_model_parms(roundtrip_file));
    
    // Compare
    auto orig_events = original_parms.get_event_list();
    auto loaded_events = loaded_parms.get_event_list();
    REQUIRE(orig_events.size() == loaded_events.size());
}

// Test format detection
TEST_CASE("Format version detection", "[model_format]")
{
    // Test legacy format detection (no @Version header)
    std::string legacy_file = TEST_DATA_DIR + "test_legacy_detect.txt";
    Model_Parms legacy_parms;
    REQUIRE_NOTHROW(legacy_parms.read_model_parms(legacy_file));
    
    // JSON is told apart by its content, not its file name.
    ScratchDir scratch;
    std::ifstream json_in(TEST_DATA_DIR + "vj_model.json");
    std::stringstream json_text;
    json_text << json_in.rdbuf();
    Model_Parms json_parms;
    REQUIRE_NOTHROW(json_parms.read_model_parms(scratch.file("named_like_v1_parms.txt", json_text.str())));
    REQUIRE(json_parms.get_event_list().size() == 5);
}

// What the schema does not describe is refused, rather than read as something else.
TEST_CASE("A JSON document the schema does not describe is refused", "[model_format][json]")
{
    using Catch::Matchers::ContainsSubstring;
    nlohmann::json doc = read_json(TEST_DATA_DIR + "vj_model.json");
    ScratchDir scratch;
    const auto read = [&](const nlohmann::json &edited) {
        Model_Parms parms;
        parms.read_model_parms(scratch.file("edited.json", edited.dump()));
    };

    SECTION("an unknown key") {
        doc["edges"] = nlohmann::json::array();
        REQUIRE_THROWS_WITH(read(doc), ContainsSubstring("unknown key \"edges\""));
    }
    SECTION("another schema version") {
        doc["schema_version"] = 2;
        REQUIRE_THROWS_WITH(read(doc), ContainsSubstring("schema version 2"));
    }
    SECTION("an error rate other than SingleErrorRate, whose JSON form is HJ") {
        doc["error_rate"] = {{"type", "HypermutationGlobalErrorRate"}, {"rate", 0.1}};
        REQUIRE_THROWS_WITH(read(doc), ContainsSubstring("HypermutationGlobalErrorRate") && ContainsSubstring("HJ"));
    }
    SECTION("no error rate") {
        doc.erase("error_rate");
        REQUIRE_THROWS_WITH(read(doc), ContainsSubstring("missing key \"error_rate\""));
    }
    SECTION("into a model that already has events") {
        Model_Parms parms;
        REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "vj_model.json"));
        REQUIRE_THROWS_AS(model_parms_from_json(doc, parms), std::logic_error);
    }
}

// The v2 text format was retired before any release read it.
TEST_CASE("A v2 text file is refused, and the message names JSON", "[model_format][json]")
{
    ScratchDir scratch;
    const std::string v2_file = scratch.file("json_parms.txt",
                                             "@Version\n2.0\n@Event_list\n"
                                             "#GeneChoice;V_gene;V_gene_seq;Undefined_side;7;v_choice\n%V1;ATCG;0\n"
                                             "@Edges\n@ErrorRate\n#SingleErrorRate\n0.001\n");
    Model_Parms parms;
    REQUIRE_THROWS_WITH(parms.read_model_parms(v2_file),
                        Catch::Matchers::ContainsSubstring("v2 text format")
                                && Catch::Matchers::ContainsSubstring("JSON"));
}

// A JSON document's seq_type_order populates the SeqTypeRegistry with exactly the declared
// order, including non-standard seq_type names that cannot exist in the v1 format.
TEST_CASE("SeqTypeRegistry populated from a JSON document's seq_type_order",
          "[model_format][json][seq_type_registry][tandem_d]")
{
    const std::string file_path = TEST_DATA_DIR + "tandem_d_model.json";

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(file_path));

    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    // Registry must reflect exactly the seven entries listed in seq_type_order
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

// ══════════════════════════════════════════════════════════════════════════════
// SeqTypeId: the runtime handle layer (plan decisions D1 / D5, B2 step 4)
// ══════════════════════════════════════════════════════════════════════════════

// The invariant the whole B8 migration rests on: for the six legacy seq_types the
// registry id equals the Seq_type enum value, so a map keyed by SeqTypeId and code still
// keyed by the enum address the same slot. Without the pre-registration this holds for a
// VDJ model by accident of ordering and fails for a VJ one -- J_gene_seq would take id 2
// while the enum says 4, silently aliasing two segments.
TEST_CASE("Loaded models pin legacy seq_type ids to the Seq_type enum",
          "[model_format][seq_type_id][b8]")
{
    auto check_enum_alignment = [](const Model_Parms &parms) {
        const SeqTypeRegistry &r = parms.get_seq_type_registry();
        REQUIRE(r.total_count() >= 6);
        CHECK(r.id("V_gene_seq")  == static_cast<SeqTypeId>(V_gene_seq));
        CHECK(r.id("VD_ins_seq")  == static_cast<SeqTypeId>(VD_ins_seq));
        CHECK(r.id("D_gene_seq")  == static_cast<SeqTypeId>(D_gene_seq));
        CHECK(r.id("DJ_ins_seq")  == static_cast<SeqTypeId>(DJ_ins_seq));
        CHECK(r.id("J_gene_seq")  == static_cast<SeqTypeId>(J_gene_seq));
        CHECK(r.id("VJ_ins_seq")  == static_cast<SeqTypeId>(VJ_ins_seq));
    };

    SECTION("VDJ model")
    {
        Model_Parms parms;
        REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));
        check_enum_alignment(parms);
    }

    SECTION("VJ model -- the case that breaks without pre-registration")
    {
        Model_Parms parms;
        REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
        check_enum_alignment(parms);
        // The ordering is still the VJ one; only the id space is pinned.
        const std::vector<std::string> vj = {"V_gene_seq", "VJ_ins_seq", "J_gene_seq"};
        CHECK(parms.get_seq_type_registry().get_ordered_types() == vj);
    }

    SECTION("JSON model")
    {
        Model_Parms parms;
        REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "vj_model.json"));
        check_enum_alignment(parms);
    }
}

TEST_CASE("SeqTypeRegistry: ids are dense and stable", "[seq_type_registry][seq_type_id]")
{
    SeqTypeRegistry registry;
    registry.set_ordered_types({"V_gene_seq", "VD_ins_seq", "D_gene_seq"});

    REQUIRE(registry.total_count() == 3);
    // Ids index arrays directly: exactly [0, total_count()), assigned in 5'->3' order.
    CHECK(registry.id("V_gene_seq") == 0);
    CHECK(registry.id("VD_ins_seq") == 1);
    CHECK(registry.id("D_gene_seq") == 2);

    CHECK(registry.name(0) == "V_gene_seq");
    CHECK(registry.name(2) == "D_gene_seq");

    const std::vector<SeqTypeId> expected_ordering = {0, 1, 2};
    CHECK(registry.ordering() == expected_ordering);

    CHECK_THROWS_AS(registry.id("nonexistent_seq"), std::out_of_range);
    CHECK_THROWS_AS(registry.name(7), std::out_of_range);
}

TEST_CASE("SeqTypeRegistry: register_type is idempotent", "[seq_type_registry][seq_type_id]")
{
    SeqTypeRegistry registry;
    const SeqTypeId first = registry.register_type("D1_gene_seq");
    const SeqTypeId again = registry.register_type("D1_gene_seq");
    CHECK(first == again);
    CHECK(registry.total_count() == 1);

    const SeqTypeId other = registry.register_type("D2_gene_seq");
    CHECK(other != first);
    CHECK(registry.total_count() == 2);
}

TEST_CASE("SeqTypeRegistry: id-keyed neighbour lookups", "[seq_type_registry][seq_type_id]")
{
    SeqTypeRegistry registry;
    registry.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
                                "D1D2_ins_seq", "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"});

    const SeqTypeId d1 = registry.id("D1_gene_seq");
    CHECK(registry.left_neighbor(d1)  == registry.id("VD1_ins_seq"));
    CHECK(registry.right_neighbor(d1) == registry.id("D1D2_ins_seq"));

    // The ends of the ordering report kNoSeqType rather than an optional.
    CHECK(registry.left_neighbor(registry.id("V_gene_seq")) == kNoSeqType);
    CHECK(registry.right_neighbor(registry.id("J_gene_seq")) == kNoSeqType);

    // A type registered but absent from the ordering has no neighbours.
    const SeqTypeId orphan = registry.register_type("left_flank_seq");
    CHECK(registry.left_neighbor(orphan) == kNoSeqType);
    CHECK(registry.right_neighbor(orphan) == kNoSeqType);
}

TEST_CASE("SeqTypeRegistry: freezing seals the id space", "[seq_type_registry][seq_type_id]")
{
    SeqTypeRegistry registry;
    registry.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
    REQUIRE_FALSE(registry.is_frozen());

    registry.freeze();
    CHECK(registry.is_frozen());

    // Re-registering a known name is still fine -- it allocates nothing.
    CHECK(registry.register_type("V_gene_seq") == registry.id("V_gene_seq"));
    // Anything that would grow or reorder the id space is refused.
    CHECK_THROWS_AS(registry.register_type("D1_gene_seq"), std::logic_error);
    CHECK_THROWS_AS(registry.set_ordered_types({"V_gene_seq"}), std::logic_error);

    registry.freeze();  // idempotent
    CHECK(registry.is_frozen());
}

TEST_CASE("Model_Parms::finalize resolves every event's SeqTypeId",
          "[model_format][seq_type_id]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    // read_model_parms() finalizes, so the registry is sealed on return.
    const SeqTypeRegistry &registry = parms.get_seq_type_registry();
    REQUIRE(registry.is_frozen());

    for (const auto &event : parms.get_event_list()) {
        const std::string name = event->get_seq_type();
        if (name.empty() || !registry.contains(name)) {
            continue;
        }
        INFO("event " << event->get_nickname() << " seq_type " << name);
        CHECK(event->get_seq_type_id() == registry.id(name));
        CHECK(event->get_seq_type_id() != kNoSeqType);
    }
}

TEST_CASE("Model_Parms: adding an event after finalize is refused",
          "[model_format][seq_type_id]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    auto extra = std::make_shared<Gene_choice>(V_gene);
    extra->set_nickname("late_addition");
    // It would carry an unresolved SeqTypeId, and might need a type the frozen registry
    // can no longer accept.
    CHECK_THROWS_AS(parms.add_event(extra), std::logic_error);
}

TEST_CASE("Rec_Event::copy carries the resolved SeqTypeId", "[model_format][seq_type_id]")
{
    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

    // Inference deep-copies Model_Parms per thread; the copies' events must keep their ids,
    // since the traversal reads the id, never the name.
    Model_Parms copy(original);
    // get_event_list() returns by value; bind before iterating.
    const auto orig_events = original.get_event_list();
    const auto copy_events = copy.get_event_list();
    REQUIRE(copy_events.size() == orig_events.size());

    auto orig_it = orig_events.begin();
    auto copy_it = copy_events.begin();
    for (; orig_it != orig_events.end(); ++orig_it, ++copy_it) {
        INFO("event " << (*orig_it)->get_nickname());
        CHECK((*copy_it)->get_seq_type() == (*orig_it)->get_seq_type());
        CHECK((*copy_it)->get_seq_type_id() == (*orig_it)->get_seq_type_id());
    }
}

// Model_Parms is deep-copied once per OpenMP thread in the inference loop, and the
// scenario maps are built from that copy. A copy that drops seq_type_registry leaves
// every worker with an empty ordering -- and silently downgrades write_model_parms()
// from JSON to the v1 format, since requires_extended_format() returns false on an empty
// registry.
TEST_CASE("Model_Parms copy constructor preserves the seq_type registry",
          "[model_format][seq_type_registry][copy]")
{
    SECTION("explicit JSON ordering")
    {
        Model_Parms original;
        REQUIRE_NOTHROW(original.read_model_parms(TEST_DATA_DIR + "tandem_d_model.json"));
        REQUIRE_FALSE(original.get_seq_type_registry().get_ordered_types().empty());

        Model_Parms copy(original);
        CHECK(copy.get_seq_type_registry().get_ordered_types()
              == original.get_seq_type_registry().get_ordered_types());
    }

    SECTION("ordering inferred from a legacy VDJ file")
    {
        Model_Parms original;
        REQUIRE_NOTHROW(original.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));

        Model_Parms copy(original);
        const std::vector<std::string> expected_vdj = {
            "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"
        };
        CHECK(copy.get_seq_type_registry().get_ordered_types() == expected_vdj);

        // Neighbour lookups exercise the registry's name->index map, not just the
        // ordered vector, so this also catches a partial copy.
        const SeqTypeRegistry &copied = copy.get_seq_type_registry();
        CHECK(copied.index_of("D_gene_seq") == 2);
        CHECK(copied.get_left_neighbor("D_gene_seq") == std::optional<std::string>("VD_ins_seq"));
        CHECK(copied.get_right_neighbor("D_gene_seq") == std::optional<std::string>("DJ_ins_seq"));
    }

    SECTION("copy keeps the same write_model_parms format selection")
    {
        Model_Parms original;
        REQUIRE_NOTHROW(original.read_model_parms(TEST_DATA_DIR + "tandem_d_model.json"));
        REQUIRE(original.requires_extended_format());

        Model_Parms copy(original);
        CHECK(copy.requires_extended_format());
    }
}

// ── Gene classes in the JSON format ────────────────────────────────────────────────

// A GeneChoice's gene_class names the gene library it aligns against: V_gene, D_gene,
// J_gene or Undefined_gene. The v1 junction classes (VD_genes, …) are refused.
TEST_CASE("JSON format: GeneChoice with junction gene_class throws", "[model_format][json][gene_class]")
{
    nlohmann::json doc = read_json(TEST_DATA_DIR + "vj_model.json");
    for (auto &node : doc.at("events")) {
        if (node.at("nickname") == "v_choice") {
            node["gene_class"] = "VD_genes";
        }
    }
    ScratchDir scratch;
    const std::string file_path = scratch.file("bad_gene_class.json", doc.dump());

    Model_Parms parms;
    REQUIRE_THROWS_WITH(parms.read_model_parms(file_path), Catch::Matchers::ContainsSubstring("VD_genes"));
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

// In the JSON format Insertion and DinucMarkov carry no gene_class, and their seq_type is
// the document's own field.
TEST_CASE("JSON format: Insertion and DinucMarkov use explicit seq_type field",
          "[model_format][json][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "vj_model.json";

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

// In the JSON format a GeneChoice with a valid alignment gene_class is accepted, and the
// seq_type field is stored verbatim.
TEST_CASE("JSON format: GeneChoice with valid gene_class accepted",
          "[model_format][json][gene_class]")
{
    const std::string file_path = TEST_DATA_DIR + "vj_model.json";

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

// ── v1/JSON parity tests ────────────────────────────────────────────────────

// Helper: collect comparable per-event properties from a loaded Model_Parms.
struct EventSummary {
    std::string nickname;
    Event_type  type;
    Gene_class  gene_class;  // only meaningful for GeneChoice
    Seq_type_String seq_type;
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

static void require_same_events(const Model_Parms &v1_parms, const Model_Parms &json_parms)
{
    auto v1 = summarise_events(v1_parms);
    auto json = summarise_events(json_parms);
    REQUIRE(v1.size() == json.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        INFO("Event: " << v1[i].nickname);
        REQUIRE(v1[i].nickname         == json[i].nickname);
        REQUIRE(v1[i].type             == json[i].type);
        REQUIRE(v1[i].seq_type         == json[i].seq_type);
        REQUIRE(v1[i].num_realizations == json[i].num_realizations);
        // For GeneChoice events the alignment gene_class must be preserved.
        if (v1[i].type == GeneChoice_t) {
            REQUIRE(v1[i].gene_class == json[i].gene_class);
        }
    }
}

// VJ legacy vs its JSON counterpart.
TEST_CASE("v1/JSON parity: VJ model - JSON file matches legacy",
          "[model_format][parity][json]")
{
    Model_Parms v1_parms, json_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
    REQUIRE_NOTHROW(json_parms.read_model_parms(TEST_DATA_DIR + "vj_model.json"));

    require_same_events(v1_parms, json_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            json_parms.get_seq_type_registry().get_ordered_types());
}

// VDJ legacy vs its JSON counterpart.
TEST_CASE("v1/JSON parity: VDJ model - JSON file matches legacy",
          "[model_format][parity][json]")
{
    Model_Parms v1_parms, json_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt"));
    REQUIRE_NOTHROW(json_parms.read_model_parms(TEST_DATA_DIR + "vdj_model.json"));

    require_same_events(v1_parms, json_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            json_parms.get_seq_type_registry().get_ordered_types());
}

// Round-trip: read legacy -> write JSON -> read back -> same as legacy.
TEST_CASE("v1/JSON parity: write_model_parms_json round-trip equals original",
          "[model_format][parity][json]")
{
    const std::string v1_path = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    ScratchDir scratch;
    const std::string tmp_json  = scratch.file("tmp_roundtrip.json");

    Model_Parms v1_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(v1_path));

    REQUIRE_NOTHROW(v1_parms.write_model_parms_json(tmp_json));

    Model_Parms json_parms;
    REQUIRE_NOTHROW(json_parms.read_model_parms(tmp_json));

    require_same_events(v1_parms, json_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            json_parms.get_seq_type_registry().get_ordered_types());
}

// Round-trip using the large inference model file (covers full V/D/J gene sets).
TEST_CASE("v1/JSON parity: write_model_parms_json round-trip for inference model",
          "[model_format][parity][json][inference]")
{
    const std::string INFER_DIR = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/inference/";
    const std::string v1_path   = INFER_DIR + "fixed_VJ_TRB_model_parms.txt";
    ScratchDir scratch;
    const std::string tmp_json    = scratch.file("tmp_inference.json");

    Model_Parms v1_parms;
    REQUIRE_NOTHROW(v1_parms.read_model_parms(v1_path));

    REQUIRE_NOTHROW(v1_parms.write_model_parms_json(tmp_json));

    Model_Parms json_parms;
    REQUIRE_NOTHROW(json_parms.read_model_parms(tmp_json));

    require_same_events(v1_parms, json_parms);

    REQUIRE(v1_parms.get_seq_type_registry().get_ordered_types() ==
            json_parms.get_seq_type_registry().get_ordered_types());
}

// A shipped model written as JSON and read back: its events, order and edges. Counterpart of
// the reduced v2 TRB model the models submodule carried.
TEST_CASE("A shipped model written as JSON reads back the same", "[model_format][json][integration]")
{
    const std::string source = std::string(IGOR_SOURCE_DIR) + "/models/human/tcr_beta/models/model_parms.txt";
    if (!std::filesystem::exists(source)) {
        SKIP("model file not found at " << source << " (submodule may not be initialized)");
    }
    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(source));

    ScratchDir scratch;
    const std::string written = scratch.file("human_trb.json");
    REQUIRE_NOTHROW(original.write_model_parms_json(written));

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(written));
    require_same_events(original, loaded);
    CHECK(loaded.get_seq_type_registry().get_ordered_types()
          == original.get_seq_type_registry().get_ordered_types());
    CHECK(children_by_nickname(loaded) == children_by_nickname(original));
    CHECK(loaded.get_event_pointer("vd_ins", true)->get_seq_type() == "VD_ins_seq");
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

// Tandem-D model has non-standard seq_types (D1_gene_seq, VD1_ins_seq …) → must use JSON.
TEST_CASE("requires_extended_format: tandem-D model returns true",
          "[model_format][step7][requires_extended_format][tandem_d]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "tandem_d_model.json"));
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

// For a standard VDJ model, write_model_parms() must produce a v1 file.
TEST_CASE("write_model_parms auto-selects legacy format for standard VDJ model",
          "[model_format][step7][auto_select]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_auto_select_legacy.txt");

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms(dest));

    REQUIRE(first_line_of(dest) == "@Event_list");
    CHECK(parms.model_file_extension() == ".txt");

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(parms, loaded);
}

// For a tandem-D model with non-standard seq_types, write_model_parms() must produce a JSON
// document, and its output files are named so.
TEST_CASE("write_model_parms auto-selects JSON for tandem-D model",
          "[model_format][step7][auto_select][tandem_d]")
{
    const std::string src  = TEST_DATA_DIR + "tandem_d_model.json";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_auto_select.json");

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms(dest));

    REQUIRE(first_line_of(dest) == "{");
    CHECK(parms.model_file_extension() == ".json");

    // Must also be readable back with the same events
    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(parms, loaded);
}

// ── Step 7: write_model_parms_legacy() round-trip ──────────────────────────

// Read a VDJ legacy model, write it via write_model_parms_legacy(), read back,
// and verify that all events are identical.
TEST_CASE("write_model_parms_legacy: VDJ round-trip preserves events",
          "[model_format][step7][legacy_roundtrip]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_legacy_roundtrip.txt");

    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(src));
    REQUIRE_NOTHROW(original.write_model_parms_legacy(dest));

    REQUIRE(first_line_of(dest) == "@Event_list");

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));
    require_same_events(original, loaded);

    REQUIRE(original.get_seq_type_registry().get_ordered_types() ==
            loaded.get_seq_type_registry().get_ordered_types());
}

// ── Step 5: the JSON event nodes ───────────────────────────────────────────

// Every event node of a written document carries its seq type; only a gene choice carries a
// gene class (counterpart of the v2 event line's seq_type field).
TEST_CASE("write_model_parms_json: every event carries its seq type", "[model_format][step5][json]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_json_nodes.json");

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms_json(dest));

    const nlohmann::json doc = read_json(dest);
    REQUIRE(doc.at("events").size() == parms.get_event_list().size());
    for (const auto &node : doc.at("events")) {
        const auto nickname = node.at("nickname").get<std::string>();
        INFO("event " << nickname);
        CHECK(node.at("seq_type") == parms.get_event_pointer(nickname, true)->get_seq_type());
        CHECK(node.contains("gene_class") == (node.at("type") == "GeneChoice"));
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// Step 8: how a file names an event in its edges
// ══════════════════════════════════════════════════════════════════════════════
//
// A v1 file's @Edges names events by their v1 name:
//   %GeneChoice_V_gene_Undefined_side_prio7_size1;Deletion_V_gene_Three_prime_prio5_size2
// which tells two events apart only where a v1 gene class does. A JSON document names each
// event's parents by nickname, the event's name (plan T2). The v2 format named them by a
// generated name carrying the seq type; it is retired.
// ──────────────────────────────────────────────────────────────────────────────

// v1 names do NOT include seq_type: they are what the v1 files carry.
TEST_CASE("step8: v1 event names do not include seq_type",
          "[model_format][step8][event_name]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));

    for (const auto &ev : parms.get_event_list()) {
        const std::string name = ev->get_v1_legacy_name();
        const std::string &st   = ev->get_seq_type();
        INFO("v1 name: " << name);
        REQUIRE(!st.empty());
        REQUIRE(name.find(st) == std::string::npos);
    }
}

// A JSON document names each event by its nickname, which is the event's name.
TEST_CASE("step8: a JSON document names each event by its name",
          "[model_format][step8][event_name][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "test_legacy_model_parms.txt"));
    const nlohmann::json doc = model_parms_to_json(parms);

    std::set<std::string> nicknames;
    for (const auto &node : doc.at("events")) {
        nicknames.insert(node.at("nickname").get<std::string>());
    }
    std::set<std::string> names;
    for (const auto &event : parms.get_event_list()) {
        names.insert(event->get_name());
    }
    CHECK(nicknames == names);
}

// write_model_parms_json() names each edge's parent by nickname.
TEST_CASE("step8: write_model_parms_json names parents by nickname",
          "[model_format][step8][event_name][json]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_model_parms.txt";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_step8_parents.json");

    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(src));
    REQUIRE_NOTHROW(parms.write_model_parms_json(dest));

    std::map<std::string, std::vector<std::string>> parents;
    for (const auto &node : read_json(dest).at("events")) {
        parents[node.at("nickname").get<std::string>()] = node.at("parents").get<std::vector<std::string>>();
    }
    CHECK(parents.at("j_choice") == std::vector<std::string>{"v_choice"});
    CHECK(parents.at("v_3_del") == std::vector<std::string>{"v_choice"});
    CHECK(parents.at("v_choice").empty());
}

// A JSON document's parents wire the model's edges.
TEST_CASE("step8: a JSON document's edges are wired",
          "[model_format][step8][parser][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "vj_model.json"));
    REQUIRE(parms.get_event_list().size() == 5);  // v_choice, j_choice, v_3_del, vj_ins, vj_dinucl

    const auto graph = children_by_nickname(parms);
    CHECK(graph.at("v_choice") == std::set<std::string>{"j_choice", "v_3_del"});
    CHECK(graph.at("j_choice").empty());
}

// A parent the document does not name is refused.
TEST_CASE("step8: a JSON parent no event is named is refused", "[model_format][step8][parser][json]")
{
    nlohmann::json doc = read_json(TEST_DATA_DIR + "vj_model.json");
    for (auto &node : doc.at("events")) {
        if (node.at("nickname") == "j_choice") {
            node["parents"] = {"v_gene"};
        }
    }
    ScratchDir scratch;
    Model_Parms parms;
    REQUIRE_THROWS_WITH(parms.read_model_parms(scratch.file("unknown_parent.json", doc.dump())),
                        Catch::Matchers::ContainsSubstring("unknown parent \"v_gene\""));
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
          "[model_format][step9][seq_order][tandem_d]")
{
    // This part (reading + checking registry) is already implemented.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "tandem_d_model.json"));

    const std::vector<std::string> expected = {
        "V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
        "D1D2_ins_seq", "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    REQUIRE(parms.get_seq_type_registry().get_ordered_types() == expected);
}

TEST_CASE("step9: build_scenario_sequence uses registry order for tandem D",
          "[model_format][step9][seq_order][tandem_d]")
{
    // Load the tandem-D model to get its registry, then call the registry-based
    // build_scenario_sequence() and verify 7-segment assembly order.
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "tandem_d_model.json"));
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

// Non-standard seq_types are addressable in the real Seq_offsets_map, which is keyed by
// SeqTypeId and sized from the registry. This supersedes the Str_Seq_offsets_map prototype
// (a nested hash map) that used to stand in for it here.
TEST_CASE("Seq_offsets_map accepts non-standard seq_types",
          "[model_format][seq_offsets][tandem_d]")
{
    SeqTypeRegistry registry;
    registry.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
                                "D1D2_ins_seq", "D2_gene_seq", "D2J_ins_seq", "J_gene_seq"});
    registry.freeze();

    Seq_offsets_map offsets(registry);
    const SeqTypeId d1 = registry.id("D1_gene_seq");
    const SeqTypeId d2 = registry.id("D2_gene_seq");

    CHECK_FALSE(offsets.exists(d1, Five_prime));

    offsets.request_layer(d1, Five_prime);
    offsets.set(d1, Five_prime, 42, 0);
    CHECK(offsets.exists(d1, Five_prime));
    CHECK(offsets.get(d1, Five_prime) == 42);
    // The two ends are independent maps, not two rows of one address space.
    CHECK_FALSE(offsets.exists(d1, Three_prime));
    CHECK_FALSE(offsets.exists(d2, Five_prime));

    offsets.request_layer(d1, Five_prime);
    offsets.set(d1, Five_prime, 99, 1);
    CHECK(offsets.get(d1, Five_prime) == 99);
    CHECK(offsets.claimed_layer(d1, Five_prime) == 1);
    offsets.restore_layer(d1, Five_prime);
    CHECK(offsets.get(d1, Five_prime) == 42);

    // An offset is one end of a segment or the other.
    CHECK_THROWS_AS(offsets.get(d1, Undefined_side), std::out_of_range);
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
// it on read; JSON documents carry the explicit value verbatim.
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

// A DinucMarkov's side is stored verbatim in a JSON document and survives a
// write_model_parms_json → read round-trip.
TEST_CASE("step11: DinucMarkov direction preserved through JSON round-trip",
          "[model_format][step11][dinucl_direction][json]")
{
    const std::string src  = TEST_DATA_DIR + "test_legacy_vdj_model_parms.txt";
    ScratchDir scratch;
    const std::string dest = scratch.file("tmp_step11_dinucl.json");

    Model_Parms original;
    REQUIRE_NOTHROW(original.read_model_parms(src));

    REQUIRE_NOTHROW(original.write_model_parms_json(dest));

    Model_Parms loaded;
    REQUIRE_NOTHROW(loaded.read_model_parms(dest));

    auto vd_orig   = original.get_event_pointer("vd_dinucl", true);
    auto vd_loaded = loaded.get_event_pointer("vd_dinucl", true);
    REQUIRE(vd_loaded->get_side() == vd_orig->get_side());

    auto dj_orig   = original.get_event_pointer("dj_dinucl", true);
    auto dj_loaded = loaded.get_event_pointer("dj_dinucl", true);
    REQUIRE(dj_loaded->get_side() == dj_orig->get_side());
}

// A JSON document carries a DinucMarkov's side explicitly; it is loaded verbatim, not derived.
TEST_CASE("step11: DinucMarkov direction in a JSON document is loaded verbatim",
          "[model_format][step11][dinucl_direction][json]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(TEST_DATA_DIR + "vj_model.json"));

    // As the VJ_genes v1 convention gives it: Three_prime, the chain seeded from V.
    auto vj_dinucl = parms.get_event_pointer("vj_dinucl", true);
    REQUIRE(vj_dinucl != nullptr);
    REQUIRE(vj_dinucl->get_type() == Dinuclmarkov_t);
    REQUIRE(vj_dinucl->get_side() == Three_prime);

    SECTION("whatever the convention would say") {
        nlohmann::json doc = read_json(TEST_DATA_DIR + "vj_model.json");
        for (auto &node : doc.at("events")) {
            if (node.at("nickname") == "vj_dinucl") {
                node["side"] = "Five_prime";
            }
        }
        ScratchDir scratch;
        Model_Parms edited;
        REQUIRE_NOTHROW(edited.read_model_parms(scratch.file("five_prime.json", doc.dump())));
        CHECK(edited.get_event_pointer("vj_dinucl", true)->get_side() == Five_prime);
    }
}
