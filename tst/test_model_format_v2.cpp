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
#include <igor/Core/Utils.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

// Test that legacy format files can still be read
TEST_CASE("Legacy model format backward compatibility", "[model_format][legacy]")
{
    // Create a minimal legacy format model file
    std::string temp_file = "/tmp/test_legacy_model_parms.txt";
    
    std::ofstream outfile(temp_file);
    outfile << "@Event_list" << std::endl;
    
    // GeneChoice event (legacy format: type;gene_class;seq_side;priority;nickname)
    outfile << "#GeneChoice;V_gene;Undefined_side;7;v_choice" << std::endl;
    outfile << "%V_test;ATCGATCG;0" << std::endl;
    
    outfile << "#GeneChoice;J_gene;Undefined_side;7;j_choice" << std::endl;
    outfile << "%J_test;GCTAGCTA;0" << std::endl;
    
    // Deletion event
    outfile << "#Deletion;V_gene;Three_prime;5;v_3_del" << std::endl;
    outfile << "%0;0" << std::endl;
    outfile << "%1;1" << std::endl;
    
    // Insertion event
    outfile << "#Insertion;VJ_gene;Undefined_side;4;vj_ins" << std::endl;
    outfile << "%0;0" << std::endl;
    outfile << "%1;1" << std::endl;
    
    // DinucMarkov event
    outfile << "#DinucMarkov;VJ_gene;Undefined_side;1;vj_dinucl" << std::endl;
    outfile << "%A;0" << std::endl;
    outfile << "%C;1" << std::endl;
    outfile << "%G;2" << std::endl;
    outfile << "%T;3" << std::endl;
    
    outfile << "@Edges" << std::endl;
    // Edge names use the internal event name format (legacy-compatible)
    outfile << "%GeneChoice_V_gene_Undefined_side_prio7_size1;Deletion_V_gene_Three_prime_prio5_size2" << std::endl;
    outfile << "%GeneChoice_V_gene_Undefined_side_prio7_size1;GeneChoice_J_gene_Undefined_side_prio7_size1" << std::endl;
    
    outfile << "@ErrorRate" << std::endl;
    outfile << "#SingleErrorRate" << std::endl;
    outfile << "0.001" << std::endl;
    outfile.close();
    
    // Read the legacy format file
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(temp_file));
    
    // Verify events were loaded correctly
    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5); // 2 GeneChoice + 1 Deletion + 1 Insertion + 1 DinucMarkov (error rate is separate)
    
    // Check that V gene choice event has correct properties
    auto v_choice = parms.get_event_pointer("v_choice", true);
    REQUIRE(v_choice != nullptr);
    REQUIRE(v_choice->get_type() == GeneChoice_t);
    REQUIRE(v_choice->get_class() == V_gene);
    
    // Cleanup
    std::filesystem::remove(temp_file);
}

// Test that v2.0 format files can be read (placeholder for future implementation)
TEST_CASE("Model format v2.0 parsing", "[model_format][v2.0][!mayfail]")
{
    // Create a v2.0 format model file
    std::string temp_file = "/tmp/test_v2_model_parms.txt";
    
    std::ofstream outfile(temp_file);
    outfile << "@Version" << std::endl;
    outfile << "2.0" << std::endl;
    
    outfile << "@Seq_type_order" << std::endl;
    outfile << "[V_gene_seq,VJ_ins_seq,J_gene_seq]" << std::endl;
    
    outfile << "@Event_list" << std::endl;
    
    // GeneChoice event (v2.0 format: type;gene_class;seq_type;seq_side;priority;nickname)
    outfile << "#GeneChoice;V_gene;V_gene_seq;Undefined_side;7;v_choice" << std::endl;
    outfile << "%V_test;ATCGATCG;0" << std::endl;
    
    outfile << "#GeneChoice;J_gene;J_gene_seq;Undefined_side;7;j_choice" << std::endl;
    outfile << "%J_test;GCTAGCTA;0" << std::endl;
    
    // Deletion event (v2.0: gene_class is Undefined_gene, seq_type identifies target)
    outfile << "#Deletion;Undefined_gene;V_gene_seq;Three_prime;5;v_3_del" << std::endl;
    outfile << "%0;0" << std::endl;
    outfile << "%1;1" << std::endl;
    
    // Insertion event
    outfile << "#Insertion;Undefined_gene;VJ_ins_seq;Undefined_side;4;vj_ins" << std::endl;
    outfile << "%0;0" << std::endl;
    outfile << "%1;1" << std::endl;
    
    // DinucMarkov event (v2.0: seq_side is explicit)
    outfile << "#DinucMarkov;Undefined_gene;VJ_ins_seq;Three_prime;1;vj_dinucl" << std::endl;
    outfile << "%A;0" << std::endl;
    outfile << "%C;1" << std::endl;
    outfile << "%G;2" << std::endl;
    outfile << "%T;3" << std::endl;
    
    outfile << "@Edges" << std::endl;
    // Edge names must match the internal event names (after v2.0 -> legacy conversion)
    // Deletion with seq_type=V_gene_seq gets gene_class=V_gene for internal use
    outfile << "%GeneChoice_V_gene_Undefined_side_prio7_size1;Deletion_V_gene_Three_prime_prio5_size2" << std::endl;
    outfile << "%GeneChoice_V_gene_Undefined_side_prio7_size1;GeneChoice_J_gene_Undefined_side_prio7_size1" << std::endl;
    
    outfile << "@ErrorRate" << std::endl;
    outfile << "#SingleErrorRate" << std::endl;
    outfile << "0.001" << std::endl;
    outfile.close();
    
    // Read the v2.0 format file
    Model_Parms parms;
    REQUIRE_NOTHROW(parms.read_model_parms(temp_file));
    
    // Verify events were loaded correctly
    auto events = parms.get_event_list();
    REQUIRE(events.size() == 5); // 2 GeneChoice + 1 Deletion + 1 Insertion + 1 DinucMarkov (error rate is separate)
    
    // Cleanup
    std::filesystem::remove(temp_file);
}

// Test round-trip: write and read back
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
    
    // Write to file
    std::string temp_file = "/tmp/test_roundtrip_parms.txt";
    original_parms.write_model_parms(temp_file);
    
    // Read back
    Model_Parms loaded_parms;
    REQUIRE_NOTHROW(loaded_parms.read_model_parms(temp_file));
    
    // Compare
    auto orig_events = original_parms.get_event_list();
    auto loaded_events = loaded_parms.get_event_list();
    REQUIRE(orig_events.size() == loaded_events.size());
    
    // Cleanup
    std::filesystem::remove(temp_file);
}

// Test format detection
TEST_CASE("Format version detection", "[model_format]")
{
    // Test legacy format detection (no @Version header)
    std::string legacy_file = "/tmp/test_legacy_detect.txt";
    std::ofstream legacy_out(legacy_file);
    legacy_out << "@Event_list" << std::endl;
    legacy_out << "#GeneChoice;V_gene;Undefined_side;7;v_choice" << std::endl;
    legacy_out << "%V1;ATCG;0" << std::endl;
    legacy_out << "@Edges" << std::endl;
    legacy_out << "@ErrorRate" << std::endl;
    legacy_out << "#SingleErrorRate" << std::endl;
    legacy_out << "0.001" << std::endl;
    legacy_out.close();
    
    Model_Parms legacy_parms;
    REQUIRE_NOTHROW(legacy_parms.read_model_parms(legacy_file));
    
    // Test v2.0 format detection
    std::string v2_file = "/tmp/test_v2_detect.txt";
    std::ofstream v2_out(v2_file);
    v2_out << "@Version" << std::endl;
    v2_out << "2.0" << std::endl;
    v2_out << "@Event_list" << std::endl;
    v2_out << "#GeneChoice;V_gene;V_gene_seq;Undefined_side;7;v_choice" << std::endl;
    v2_out << "%V1;ATCG;0" << std::endl;
    v2_out << "@Edges" << std::endl;
    v2_out << "@ErrorRate" << std::endl;
    v2_out << "#SingleErrorRate" << std::endl;
    v2_out << "0.001" << std::endl;
    v2_out.close();
    
    Model_Parms v2_parms;
    REQUIRE_NOTHROW(v2_parms.read_model_parms(v2_file));
    
    // Cleanup
    std::filesystem::remove(legacy_file);
    std::filesystem::remove(v2_file);
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
