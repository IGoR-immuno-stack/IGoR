#include <catch2/catch_test_macros.hpp>
#include <fstream> // Added for ofstream in MockEvent
#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/gene_to_seqtype_migr.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Rec_Event.h>
#include <memory>
#include <queue>  // Added for queue in MockEvent
#include <random> // Added for mt19937_64 in MockEvent
#include <tuple>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace EventUtils;

// Mock Rec_Event for testing
class MockEvent : public Rec_Event {
public:
  MockEvent(string name) : Rec_Event() { this->name = name; }
  void
  iterate(QuerySequenceContext& query,
          const ModelContext& model,
          ScenarioContext& scenario,
          ExplorationContext& exploration,
          AccumulationContext& accumulation) override {}

  queue<int> draw_random_realization(
      const Marginal_array_p &, unordered_map<Rec_Event_name, int> &,
      const unordered_map<Rec_Event_name,
                          vector<pair<shared_ptr<const Rec_Event>, int>>> &,
      unordered_map<Seq_type, string> &, mt19937_64 &) const override {
    return queue<int>();
  }
  void write2txt(ofstream &) override {}
  void write2txt_legacy(ofstream &) override {}
  void write2txt_v2(ofstream &) override {}
  void initialize_event(
      unordered_set<Rec_Event_name> &,
      const Events_map &,
      const unordered_map<Rec_Event_name,
                          vector<pair<shared_ptr<const Rec_Event>, int>>> &,
      Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &,
      Safety_bool_map &, shared_ptr<Error_rate>, Mismatch_vectors_map &,
      Seq_offsets_map &, Index_map &) override {}
  void add_to_marginals(long double, Marginal_array_p &) const override {}
  shared_ptr<Rec_Event> copy() override { return nullptr; }
  bool has_effect_on(Seq_type) const override { return false; }
  void iterate_initialize_Len_proba(Seq_type, map<int, double> &,
                                    queue<shared_ptr<Rec_Event>> &, double &,
                                    const Marginal_array_p &, Index_map &,
                                    Seq_type_str_p_map &,
                                    int &) const override {}
  void initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &,
                                  const Marginal_array_p &,
                                  Index_map &) override {}
};

TEST_CASE("EventUtils CheckGeneChoice", "[EventUtils]") {
  Events_map events_map;
  unordered_set<Rec_Event_name> processed_events;

  SECTION("Event does not exist") {
    auto status = check_gene_choice("V_gene_seq", events_map, processed_events);
    REQUIRE_FALSE(status.exists);
    REQUIRE_FALSE(status.chosen);
    REQUIRE(status.event_ptr == nullptr);
  }

  SECTION("Event exists but not processed") {
    auto v_event = make_shared<MockEvent>("V_choice");
    events_map[make_tuple(GeneChoice_t, string("V_gene_seq"), Undefined_side)] = v_event;

    auto status = check_gene_choice("V_gene_seq", events_map, processed_events);
    REQUIRE(status.exists);
    REQUIRE_FALSE(status.chosen);
    REQUIRE(status.event_ptr == v_event);
  }

  SECTION("Event exists and processed") {
    auto v_event = make_shared<MockEvent>("V_choice");
    events_map[make_tuple(GeneChoice_t, string("V_gene_seq"), Undefined_side)] = v_event;

    processed_events.insert(v_event->get_name());
    auto status = check_gene_choice("V_gene_seq", events_map, processed_events);
    REQUIRE(status.exists);
    REQUIRE(status.chosen);
    REQUIRE(status.event_ptr == v_event);
  }
}

TEST_CASE("EventUtils BuildScenarioSequence", "[EventUtils]") {
  Seq_type_str_p_map constructed_sequences(6);

  Int_Str v_seq = {0, 1, 2}; // A C G
  Int_Str d_seq = {3, 0};    // T A
  Int_Str j_seq = {1, 2};    // C G
  Int_Str vd_ins = {3};      // T
  Int_Str dj_ins = {0};      // A
  Int_Str vj_ins = {1};      // C

  constructed_sequences.set_value(V_gene_seq, &v_seq, 0);
  constructed_sequences.set_value(D_gene_seq, &d_seq, 0);
  constructed_sequences.set_value(J_gene_seq, &j_seq, 0);
  constructed_sequences.set_value(VD_ins_seq, &vd_ins, 0);
  constructed_sequences.set_value(DJ_ins_seq, &dj_ins, 0);
  constructed_sequences.set_value(VJ_ins_seq, &vj_ins, 0);

  SECTION("V-D-J with insertions") {
    Int_Str result = build_scenario_sequence(constructed_sequences, true, true,
                                             true, true, true, false);
    Int_Str expected = {0, 1, 2, 3, 3, 0, 0, 1, 2}; // V + VD + D + DJ + J
    REQUIRE(static_cast<std::vector<int>>(result) ==
            static_cast<std::vector<int>>(expected));
  }

  SECTION("V-J with insertion") {
    Int_Str result = build_scenario_sequence(constructed_sequences, true, false,
                                             true, false, false, true);
    Int_Str expected = {0, 1, 2, 1, 1, 2}; // V + VJ + J
    REQUIRE(static_cast<std::vector<int>>(result) ==
            static_cast<std::vector<int>>(expected));
  }
}

// Step 9: build_scenario_sequence must eventually use SeqTypeRegistry order
// rather than hardcoded boolean flags.  These tests document both the current
// (working) standard behaviour and the desired (not yet implemented) registry-
// based behaviour.
TEST_CASE("step9: build_scenario_sequence VJ produces V+VJ_ins+J",
          "[EventUtils][step9][seq_order]") {
    Seq_type_str_p_map constructed_sequences(6);
    Int_Str v_seq = {0, 1};   // A C
    Int_Str j_seq = {2, 3};   // G T
    Int_Str vj_ins = {0};     // A
    constructed_sequences.set_value(V_gene_seq,  &v_seq,  0);
    constructed_sequences.set_value(J_gene_seq,  &j_seq,  0);
    constructed_sequences.set_value(VJ_ins_seq,  &vj_ins, 0);

    Int_Str result = build_scenario_sequence(
            constructed_sequences, true, false, true, false, false, true);
    Int_Str expected = {0, 1, 0, 2, 3};  // V + VJ_ins + J
    REQUIRE(static_cast<std::vector<int>>(result) ==
            static_cast<std::vector<int>>(expected));
}

TEST_CASE("step9: build_scenario_sequence VDJ order is V+VD+D+DJ+J",
          "[EventUtils][step9][seq_order]") {
    Seq_type_str_p_map constructed_sequences(6);
    Int_Str v_seq  = {0};   Int_Str vd_ins = {1};
    Int_Str d_seq  = {2};   Int_Str dj_ins = {3};
    Int_Str j_seq  = {0};
    constructed_sequences.set_value(V_gene_seq,  &v_seq,  0);
    constructed_sequences.set_value(VD_ins_seq,  &vd_ins, 0);
    constructed_sequences.set_value(D_gene_seq,  &d_seq,  0);
    constructed_sequences.set_value(DJ_ins_seq,  &dj_ins, 0);
    constructed_sequences.set_value(J_gene_seq,  &j_seq,  0);

    Int_Str result = build_scenario_sequence(
            constructed_sequences, true, true, true, true, true, false);
    Int_Str expected = {0, 1, 2, 3, 0};  // V + VD + D + DJ + J
    REQUIRE(static_cast<std::vector<int>>(result) ==
            static_cast<std::vector<int>>(expected));
}

// Tandem-D ordering (V→VD1ins→D1→D1D2ins→D2→DJins→J) cannot be expressed
TEST_CASE("step9: build_scenario_sequence tandem-D uses registry-based API",
          "[EventUtils][step9][seq_order]")
{
    // The registry-based overload supports arbitrary orderings, including 7-segment
    // tandem D:   V → VD1_ins → D1 → D1D2_ins → D2 → DJ_ins → J
    SeqTypeRegistry registry;
    registry.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq",
                                 "D1D2_ins_seq", "D2_gene_seq", "DJ_ins_seq", "J_gene_seq"});

    Int_Str v    = {0, 0, 0};    // length 3
    Int_Str vd1  = {1};          // length 1
    Int_Str d1   = {2, 2};       // length 2
    Int_Str d1d2 = {3};          // length 1
    Int_Str d2   = {4, 4};       // length 2
    Int_Str dj   = {5};          // length 1
    Int_Str j    = {6, 6, 6};    // length 3

    std::unordered_map<std::string, const Int_Str *> seqs;
    seqs["V_gene_seq"]   = &v;
    seqs["VD1_ins_seq"]  = &vd1;
    seqs["D1_gene_seq"]  = &d1;
    seqs["D1D2_ins_seq"] = &d1d2;
    seqs["D2_gene_seq"]  = &d2;
    seqs["DJ_ins_seq"]   = &dj;
    seqs["J_gene_seq"]   = &j;

    Int_Str result = EventUtils::build_scenario_sequence(registry, seqs);

    // Expected: v + vd1 + d1 + d1d2 + d2 + dj + j (13 elements total)
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

TEST_CASE("EventUtils GetInsertionLenMax", "[EventUtils]") {
  Events_map events_map;

  // Create a mock event with specific len_max
  class MockInsertionEvent : public MockEvent {
  public:
    MockInsertionEvent(string name, int max_len) : MockEvent(name) {
      this->len_max = max_len;
    }
  };

  auto vd_ins = make_shared<MockInsertionEvent>("VD_ins", 10);
  auto dj_ins = make_shared<MockInsertionEvent>("DJ_ins", 15);
  auto vj_ins = make_shared<MockInsertionEvent>("VJ_ins", 20);

  events_map[make_tuple(Insertion_t, string("VD_ins_seq"), Undefined_side)] = vd_ins;
  events_map[make_tuple(Insertion_t, string("DJ_ins_seq"), Undefined_side)] = dj_ins;
  events_map[make_tuple(Insertion_t, string("VJ_ins_seq"), Undefined_side)] = vj_ins;

  SECTION("VD insertion") {
    REQUIRE(get_insertion_len_max("VD_ins_seq", events_map) == 10);
  }

  SECTION("DJ insertion") {
    REQUIRE(get_insertion_len_max("DJ_ins_seq", events_map) == 15);
  }

  SECTION("VJ insertion") {
    REQUIRE(get_insertion_len_max("VJ_ins_seq", events_map) == 20);
  }
}

TEST_CASE("EventUtils GeneClassToSeqType mapping", "[EventUtils]") {
  Seq_type seq_type;

  SECTION("Valid mappings") {
    REQUIRE(igor::migration::try_gene_class_to_gene_seq_type(V_gene, seq_type));
    REQUIRE(seq_type == V_gene_seq);

    REQUIRE(igor::migration::try_gene_class_to_gene_seq_type(D_gene, seq_type));
    REQUIRE(seq_type == D_gene_seq);

    REQUIRE(igor::migration::try_gene_class_to_gene_seq_type(J_gene, seq_type));
    REQUIRE(seq_type == J_gene_seq);
  }

  SECTION("Invalid mappings") {
    REQUIRE_FALSE(igor::migration::try_gene_class_to_gene_seq_type(VD_genes, seq_type));
    REQUIRE_FALSE(igor::migration::try_gene_class_to_gene_seq_type(Undefined_gene, seq_type));
  }
}

TEST_CASE("EventUtils InsertionGeneClassToSeqType mapping", "[EventUtils]") {
  Seq_type seq_type = VD_ins_seq;

  SECTION("Valid mappings") {
    REQUIRE(igor::migration::try_insertion_gene_class_to_seq_type(VD_genes, seq_type));
    REQUIRE(seq_type == VD_ins_seq);

    REQUIRE(igor::migration::try_insertion_gene_class_to_seq_type(DJ_genes, seq_type));
    REQUIRE(seq_type == DJ_ins_seq);

    REQUIRE(igor::migration::try_insertion_gene_class_to_seq_type(VJ_genes, seq_type));
    REQUIRE(seq_type == VJ_ins_seq);
  }

  SECTION("Invalid mappings") {
    REQUIRE_FALSE(igor::migration::try_insertion_gene_class_to_seq_type(V_gene, seq_type));
    REQUIRE_FALSE(igor::migration::try_insertion_gene_class_to_seq_type(Undefined_gene, seq_type));
  }
}

TEST_CASE("EventUtils HasInsertionSeqType", "[EventUtils]") {
  Events_map events_map;

  class MockInsertionEvent : public MockEvent {
  public:
    MockInsertionEvent(string name) : MockEvent(name) {}
  };

  auto vd_ins = make_shared<MockInsertionEvent>("VD_ins");
  auto vj_ins = make_shared<MockInsertionEvent>("VJ_ins");
  events_map[make_tuple(Insertion_t, string("VD_ins_seq"), Undefined_side)] = vd_ins;

  SECTION("Finds present insertion segments") {
    REQUIRE(has_insertion_seq_type(events_map, VD_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, DJ_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, VJ_ins_seq));
  }

  SECTION("Tracks a different insertion segment independently") {
    events_map[make_tuple(Insertion_t, string("VJ_ins_seq"), Undefined_side)] = vj_ins;
    REQUIRE(has_insertion_seq_type(events_map, VD_ins_seq));
    REQUIRE(has_insertion_seq_type(events_map, VJ_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, DJ_ins_seq));
  }
}

TEST_CASE("EventUtils TryGetEvent", "[EventUtils]") {
  Events_map events_map;

  auto v_event = make_shared<MockEvent>("V_choice");
  events_map[make_tuple(GeneChoice_t, string("V_gene_seq"), Undefined_side)] = v_event;

  SECTION("Returns a present event") {
    shared_ptr<Rec_Event> event_ptr;
    REQUIRE(try_get_event(events_map, GeneChoice_t, V_gene_seq, Undefined_side, event_ptr));
    REQUIRE(event_ptr == v_event);
  }

  SECTION("Returns false for a missing event") {
    shared_ptr<Rec_Event> event_ptr;
    REQUIRE_FALSE(try_get_event(events_map, GeneChoice_t, J_gene_seq, Undefined_side, event_ptr));
    REQUIRE(event_ptr == nullptr);
  }

  SECTION("Supports direct lookup") {
    auto vd_event = make_shared<MockEvent>("VD_dinuc");
    events_map[make_tuple(Dinuclmarkov_t, string("VD_ins_seq"), Undefined_side)] = vd_event;

    shared_ptr<Rec_Event> event_ptr;
    REQUIRE(try_get_event(events_map, Dinuclmarkov_t, VD_ins_seq, Undefined_side, event_ptr));
    REQUIRE(event_ptr == vd_event);
  }
}

TEST_CASE("EventUtils InsertionPriorityBridge - specific VD key found", "[EventUtils]") {
  // Validates that try_get_event resolves VD_ins_seq via the string-keyed Events_map.
  Events_map events_map;

  auto specific_vd = make_shared<MockEvent>("specific_VD_dinuc");
  auto vdj_shared  = make_shared<MockEvent>("VDJ_dinuc_shared");

  events_map[make_tuple(Dinuclmarkov_t, string("VD_ins_seq"), Undefined_side)] = specific_vd;
  events_map[make_tuple(Dinuclmarkov_t, string("DJ_ins_seq"), Undefined_side)] = vdj_shared;

  shared_ptr<Rec_Event> result;
  REQUIRE(EventUtils::try_get_event(events_map, Dinuclmarkov_t, VD_ins_seq, Undefined_side, result));
  REQUIRE(result == specific_vd);

  REQUIRE(EventUtils::try_get_event(events_map, Dinuclmarkov_t, DJ_ins_seq, Undefined_side, result));
  REQUIRE(result == vdj_shared);
}

TEST_CASE("EventUtils InsertionPriorityBridge - missing key returns false", "[EventUtils]") {
  // Validates that try_get_event returns false when the key is absent.
  Events_map events_map;

  auto vdj_shared = make_shared<MockEvent>("VDJ_dinuc_shared");
  events_map[make_tuple(Dinuclmarkov_t, string("VD_ins_seq"), Undefined_side)] = vdj_shared;

  shared_ptr<Rec_Event> result;
  REQUIRE(EventUtils::try_get_event(events_map, Dinuclmarkov_t, VD_ins_seq, Undefined_side, result));
  REQUIRE(result == vdj_shared);

  // DJ is absent
  REQUIRE_FALSE(EventUtils::try_get_event(events_map, Dinuclmarkov_t, DJ_ins_seq, Undefined_side, result));
}


TEST_CASE("EventUtils TryEventKeyToSeqKey", "[EventUtils]") {
  tuple<Event_type, Seq_type, Seq_side> seq_key;

  SECTION("GeneChoice and Deletion mapping") {
    REQUIRE(igor::migration::try_event_key_to_seq_key(GeneChoice_t, V_gene_legacy, Undefined_side, seq_key));
    REQUIRE(seq_key == make_tuple(GeneChoice_t, V_gene_seq, Undefined_side));

    REQUIRE(igor::migration::try_event_key_to_seq_key(Deletion_t, J_gene_legacy, Five_prime, seq_key));
    REQUIRE(seq_key == make_tuple(Deletion_t, J_gene_seq, Five_prime));
  }

  SECTION("Insertion and Dinuclmarkov mapping") {
    REQUIRE(igor::migration::try_event_key_to_seq_key(Insertion_t, VD_genes, Undefined_side, seq_key));
    REQUIRE(seq_key == make_tuple(Insertion_t, VD_ins_seq, Undefined_side));

    REQUIRE(igor::migration::try_event_key_to_seq_key(Dinuclmarkov_t, DJ_genes, Undefined_side, seq_key));
    REQUIRE(seq_key == make_tuple(Dinuclmarkov_t, DJ_ins_seq, Undefined_side));
  }

  SECTION("Invalid alias (VDJ_genes) mapping") {
    REQUIRE_FALSE(igor::migration::try_event_key_to_seq_key(Insertion_t, VDJ_genes, Undefined_side, seq_key));
    REQUIRE_FALSE(igor::migration::try_event_key_to_seq_key(Dinuclmarkov_t, VDJ_genes, Undefined_side, seq_key));
  }
}

class MockDinucEvent : public MockEvent {
public:
  double updated_flag = 0.0;
  MockDinucEvent(string name) : MockEvent(name) { 
    this->set_upper_bound_proba(1.0); // Required to prevent 0-multiplication
  }
  double* get_updated_ptr() override { return &updated_flag; }
};

TEST_CASE("Insertion Bridge Integration", "[Insertion]") {
  Insertion ins_vd(VD_ins_seq, std::make_pair(0, 10));

  SECTION("VD dinuc event found via string key") {
    Events_map events_map;
    auto vd_shared = make_shared<MockDinucEvent>("VD_dinuc");
    events_map[make_tuple(Dinuclmarkov_t, string("VD_ins_seq"), Undefined_side)] = vd_shared;

    double downstream_bound = 1.0;
    forward_list<double*> updated_list = { vd_shared->get_updated_ptr() };

    REQUIRE_NOTHROW(ins_vd.initialize_crude_scenario_proba_bound(downstream_bound, updated_list, events_map));

    // initialize_crude_scenario_proba_bound removes the dinuc updated ptr from the list
    REQUIRE(distance(updated_list.begin(), updated_list.end()) == 0);
  }

  SECTION("Throws if no compatible dinuc event found") {
    Events_map events_map;
    // Add a J dinuc — not compatible with VD insertion
    events_map[make_tuple(Dinuclmarkov_t, string("J_gene_seq"), Undefined_side)] = make_shared<MockDinucEvent>("J_dinuc");

    double downstream_bound = 1.0;
    forward_list<double*> updated_list;

    // VD_genes Insertion shouldn't find anything and will throw
    REQUIRE_THROWS_AS(ins_vd.initialize_crude_scenario_proba_bound(downstream_bound, updated_list, events_map), runtime_error);
  }
}

class DeletionTest {
public:
  static void test_initialization() {
    Deletion del_event;
    del_event.event_class = V_gene;

    // Create mock objects for initialize_event — use string-keyed Events_map
    std::unordered_set<Rec_Event_name> processed_events;
    Events_map events_map;

    // Create a mock gene choice in the map to simulate it being present
    std::shared_ptr<Rec_Event> mock_v = std::make_shared<MockEvent>("V_choice");
    events_map[std::make_tuple(GeneChoice_t, std::string("V_gene_seq"), Undefined_side)] = mock_v;

    std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> offset_map;
    Downstream_scenario_proba_bound_map downstream_proba_map(6);
    Seq_type_str_p_map constructed_sequences(6);
    Safety_bool_map safety_set(6);
    std::shared_ptr<Error_rate> error_rate_p;
    Mismatch_vectors_map mismatches_list(6);
    Seq_offsets_map seq_offsets(6, 2);
    Index_map index_map(6);

    // Force initial state to true to verify that initialize_event actually overwrites it
    del_event.v_chosen = true;
    del_event.d_chosen = true;
    del_event.j_chosen = true;

    del_event.initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                               constructed_sequences, safety_set, error_rate_p, mismatches_list,
                               seq_offsets, index_map);

    // GeneChoiceStatus sets chosen only if the event is in processed_events.
    // Since processed_events is empty, chosen will be false even though V exists.
    REQUIRE(del_event.v_chosen == false);
    REQUIRE(del_event.d_chosen == false);
    REQUIRE(del_event.j_chosen == false);
  }
};

TEST_CASE("Deletion Member Initialization Regression", "[Deletion]") {
  DeletionTest::test_initialization();
}
