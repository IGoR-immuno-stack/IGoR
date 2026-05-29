#include <catch2/catch_test_macros.hpp>
#include <fstream> // Added for ofstream in MockEvent
#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/EventUtils.h>
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
