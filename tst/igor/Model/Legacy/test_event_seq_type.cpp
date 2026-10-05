/**
 * @file test_event_seq_type.cpp
 * @brief An event carries the seq type its model names, whatever the name (plan T1).
 *
 * A deletion, an insertion and a Markov chain kept their seq type twice: as the name the model
 * gives, on Rec_Event, and as a member of their own typed by the six-value Seq_type enum. The
 * member was built through str2SeqType(), so a model naming any other seq type, a tandem D's
 * D1_gene_seq for instance, could not be loaded at all.
 *
 * The first case pins what the events of a legacy model are called and what they carry, through
 * every way an event is built: the text reader, its JSON node, copy(). The name is a key of the
 * marginals and a column of the outputs, so it must not move. The [tandem_d] cases are what
 * removing the member is for. They do not pin names: two events of a tandem model still get the
 * same generated name, which is T2.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Model/EventFactory.h>
#include <igor/Model/Legacy/ModelJson.h>
#include <igor/Model/Legacy/Model_Parms.h>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <map>
#include <string>

using namespace igor::core::legacy;
using namespace igor::model::legacy;

namespace {

const std::string kTestDataDir = std::string(IGOR_SOURCE_DIR) + "/tst/test_data/format_v2/";

/// What an event is called and what it carries: the parts a model file and the outputs see.
struct Identity {
    std::string name;
    Gene_class gene_class;
    Seq_type_String seq_type;
};

void require_identity(const Rec_Event &event, const Identity &expected, bool check_name)
{
    INFO("event " << event.get_nickname());
    if (check_name) {
        CHECK(event.get_name() == expected.name);
    }
    CHECK(event.get_class() == expected.gene_class);
    CHECK(event.get_seq_type() == expected.seq_type);
}

/// The model's events, each built three ways, against one table keyed by nickname.
void require_identities(const std::string &file, const std::map<std::string, Identity> &expected,
                        bool check_name)
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(kTestDataDir + file));
    REQUIRE(parms.get_event_list().size() == expected.size());

    SECTION("from the text reader") {
        for (const auto &event : parms.get_event_list()) {
            require_identity(*event, expected.at(event->get_nickname()), check_name);
        }
    }

    SECTION("from its JSON node") {
        const nlohmann::json doc = model_parms_to_json(parms);
        for (const auto &node : doc.at("events")) {
            const auto event = igor::model::event_factory::create(node);
            require_identity(*event, expected.at(event->get_nickname()), check_name);
        }
    }

    SECTION("through copy()") {
        const Model_Parms copied(parms);
        for (const auto &event : copied.get_event_list()) {
            require_identity(*event, expected.at(event->get_nickname()), check_name);
        }
    }
}

const std::map<std::string, Identity> kVdj = {
        {"v_choice", {"GeneChoice_V_gene_Undefined_side_prio7_size1", V_gene, "V_gene_seq"}},
        {"d_choice", {"GeneChoice_D_gene_Undefined_side_prio6_size1", D_gene, "D_gene_seq"}},
        {"j_choice", {"GeneChoice_J_gene_Undefined_side_prio7_size1", J_gene, "J_gene_seq"}},
        {"v_3_del", {"Deletion_V_gene_Three_prime_prio5_size1", V_gene, "V_gene_seq"}},
        {"d_5_del", {"Deletion_D_gene_Five_prime_prio5_size1", D_gene, "D_gene_seq"}},
        {"d_3_del", {"Deletion_D_gene_Three_prime_prio5_size1", D_gene, "D_gene_seq"}},
        {"j_5_del", {"Deletion_J_gene_Five_prime_prio5_size1", J_gene, "J_gene_seq"}},
        {"vd_ins", {"Insertion_VD_genes_Undefined_side_prio4_size1", Undefined_gene, "VD_ins_seq"}},
        {"dj_ins", {"Insertion_DJ_gene_Undefined_side_prio2_size1", Undefined_gene, "DJ_ins_seq"}},
        {"vd_dinucl", {"DinucMarkov_VD_genes_Undefined_side_prio3_size16", Undefined_gene, "VD_ins_seq"}},
        {"dj_dinucl", {"DinucMarkov_DJ_gene_Undefined_side_prio1_size16", Undefined_gene, "DJ_ins_seq"}},
};

const std::map<std::string, Identity> kVj = {
        {"v_choice", {"GeneChoice_V_gene_Undefined_side_prio7_size1", V_gene, "V_gene_seq"}},
        {"j_choice", {"GeneChoice_J_gene_Undefined_side_prio7_size1", J_gene, "J_gene_seq"}},
        {"v_3_del", {"Deletion_V_gene_Three_prime_prio5_size2", V_gene, "V_gene_seq"}},
        {"vj_ins", {"Insertion_VJ_gene_Undefined_side_prio4_size2", Undefined_gene, "VJ_ins_seq"}},
        {"vj_dinucl", {"DinucMarkov_VJ_gene_Undefined_side_prio1_size16", Undefined_gene, "VJ_ins_seq"}},
};

/// A deletion's gene class names its segment only for V, D and J. Anything else, a tandem D
/// included, is Undefined_gene, as for an insertion or a Markov chain: a gene class is the
/// gene choice's, which the deletion has no way to ask at construction.
const std::map<std::string, Identity> kTandemD = {
        {"v_choice", {"", V_gene, "V_gene_seq"}},
        {"d1_choice", {"", D_gene, "D1_gene_seq"}},
        {"d2_choice", {"", D_gene, "D2_gene_seq"}},
        {"j_choice", {"", J_gene, "J_gene_seq"}},
        {"v_3_del", {"", V_gene, "V_gene_seq"}},
        {"d1_5_del", {"", Undefined_gene, "D1_gene_seq"}},
        {"d1_3_del", {"", Undefined_gene, "D1_gene_seq"}},
        {"d2_5_del", {"", Undefined_gene, "D2_gene_seq"}},
        {"d2_3_del", {"", Undefined_gene, "D2_gene_seq"}},
        {"j_5_del", {"", J_gene, "J_gene_seq"}},
        {"vd1_ins", {"", Undefined_gene, "VD1_ins_seq"}},
        {"d1d2_ins", {"", Undefined_gene, "D1D2_ins_seq"}},
        {"dj_ins", {"", Undefined_gene, "DJ_ins_seq"}},
        {"vd1_dinucl", {"", Undefined_gene, "VD1_ins_seq"}},
        {"d1d2_dinucl", {"", Undefined_gene, "D1D2_ins_seq"}},
        {"dj_dinucl", {"", Undefined_gene, "DJ_ins_seq"}},
};

}  // namespace

TEST_CASE("The events of a legacy model keep their names, gene classes and seq types",
          "[model_format][seq_type]")
{
    SECTION("VDJ, legacy format") { require_identities("test_legacy_vdj_model_parms.txt", kVdj, true); }
    SECTION("VDJ, v2 format") { require_identities("test_legacy_vdj_model_parms_v2.txt", kVdj, true); }
    SECTION("VJ, legacy format") { require_identities("test_legacy_model_parms.txt", kVj, true); }
    SECTION("VJ, v2 format") { require_identities("test_v2_model_parms.txt", kVj, true); }
}

TEST_CASE("The events of a tandem-D model carry the seq types it names",
          "[model_format][seq_type][tandem_d][!shouldfail]")
{
    require_identities("test_v2_tandem_d_parms.txt", kTandemD, false);
}

TEST_CASE("A tandem-D event resolves its seq type and its neighbours from the model's order",
          "[model_format][seq_type][tandem_d][!shouldfail]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(kTestDataDir + "test_v2_tandem_d_parms.txt"));
    const SeqTypeRegistry &registry = parms.get_seq_type_registry();

    for (const auto &event : parms.get_event_list()) {
        INFO("event " << event->get_nickname());
        const SeqTypeId id = registry.id(event->get_seq_type());
        CHECK(event->get_seq_type_id() == id);
        CHECK(event->get_left_adjacent_id() == registry.left_neighbor(id));
        CHECK(event->get_right_adjacent_id() == registry.right_neighbor(id));
    }

    //The two D blocks, spelled out: what they sit between is what their deletions trim against
    //and what their insertions fill.
    const auto d1_3_del = parms.get_event_pointer("d1_3_del", true);
    CHECK(d1_3_del->get_left_adjacent_id() == registry.id("VD1_ins_seq"));
    CHECK(d1_3_del->get_right_adjacent_id() == registry.id("D1D2_ins_seq"));
    const auto d1d2_dinucl = parms.get_event_pointer("d1d2_dinucl", true);
    CHECK(d1d2_dinucl->get_left_adjacent_id() == registry.id("D1_gene_seq"));
    CHECK(d1d2_dinucl->get_right_adjacent_id() == registry.id("D2_gene_seq"));
}

TEST_CASE("A tandem-D model written back reads the same seq types",
          "[model_format][seq_type][tandem_d][!shouldfail]")
{
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(kTestDataDir + "test_v2_tandem_d_parms.txt"));

    const std::filesystem::path written =
            std::filesystem::temp_directory_path() / "igor_test_event_seq_type_tandem_d_parms.txt";
    parms.write_model_parms(written.string());

    Model_Parms read_back;
    REQUIRE_NOTHROW(read_back.read_model_parms(written.string()));
    std::filesystem::remove(written);
    REQUIRE(read_back.get_event_list().size() == kTandemD.size());
    for (const auto &event : read_back.get_event_list()) {
        require_identity(*event, kTandemD.at(event->get_nickname()), false);
    }
}
