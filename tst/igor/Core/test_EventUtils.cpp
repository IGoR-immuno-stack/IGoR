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
      const unordered_map<tuple<Event_type, Gene_class, Seq_side>,
                          shared_ptr<Rec_Event>> &,
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
  unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>>
      events_map;
  unordered_map<tuple<Event_type, Seq_type, Seq_side>, shared_ptr<Rec_Event>>
      seq_events_map;

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

  events_map[make_tuple(Insertion_t, VD_genes, Undefined_side)] = vd_ins;
  events_map[make_tuple(Insertion_t, DJ_genes, Undefined_side)] = dj_ins;
  events_map[make_tuple(Insertion_t, VJ_genes, Undefined_side)] = vj_ins;
  seq_events_map[make_tuple(Insertion_t, VD_ins_seq, Undefined_side)] = vd_ins;
  seq_events_map[make_tuple(Insertion_t, DJ_ins_seq, Undefined_side)] = dj_ins;
  seq_events_map[make_tuple(Insertion_t, VJ_ins_seq, Undefined_side)] = vj_ins;

  SECTION("VD insertion") {
    REQUIRE(get_insertion_len_max(VD_genes, events_map) == 10);
  }

  SECTION("DJ insertion") {
    REQUIRE(get_insertion_len_max(DJ_genes, events_map) == 15);
  }

  SECTION("VJ insertion") {
    REQUIRE(get_insertion_len_max(VJ_genes, events_map) == 20);
  }

  SECTION("VD insertion typed map") {
    REQUIRE(get_insertion_len_max(VD_ins_seq, seq_events_map) == 10);
  }

  SECTION("DJ insertion typed map") {
    REQUIRE(get_insertion_len_max(DJ_ins_seq, seq_events_map) == 15);
  }

  SECTION("VJ insertion typed map") {
    REQUIRE(get_insertion_len_max(VJ_ins_seq, seq_events_map) == 20);
  }
}

TEST_CASE("EventUtils GeneClassToSeqType mapping", "[EventUtils]") {
  Seq_type seq_type = V_gene_seq;

  SECTION("Gene class maps to gene sequence type") {
    REQUIRE(try_gene_class_to_gene_seq_type(V_gene, seq_type));
    REQUIRE(seq_type == V_gene_seq);

    REQUIRE(try_gene_class_to_gene_seq_type(D_gene, seq_type));
    REQUIRE(seq_type == D_gene_seq);

    REQUIRE(try_gene_class_to_gene_seq_type(J_gene, seq_type));
    REQUIRE(seq_type == J_gene_seq);
  }

  SECTION("Non-gene class fails gene sequence mapping") {
    REQUIRE_FALSE(try_gene_class_to_gene_seq_type(VD_genes, seq_type));
    REQUIRE_FALSE(try_gene_class_to_gene_seq_type(Undefined_gene, seq_type));
  }
}

TEST_CASE("EventUtils InsertionGeneClassToSeqType mapping", "[EventUtils]") {
  Seq_type seq_type = VD_ins_seq;

  SECTION("Insertion class maps to insertion sequence type") {
    REQUIRE(try_insertion_gene_class_to_seq_type(VD_genes, seq_type));
    REQUIRE(seq_type == VD_ins_seq);

    REQUIRE(try_insertion_gene_class_to_seq_type(DJ_genes, seq_type));
    REQUIRE(seq_type == DJ_ins_seq);

    REQUIRE(try_insertion_gene_class_to_seq_type(VJ_genes, seq_type));
    REQUIRE(seq_type == VJ_ins_seq);
  }

  SECTION("Non-insertion class fails insertion sequence mapping") {
    REQUIRE_FALSE(try_insertion_gene_class_to_seq_type(V_gene, seq_type));
    REQUIRE_FALSE(try_insertion_gene_class_to_seq_type(Undefined_gene, seq_type));
  }
}

TEST_CASE("EventUtils HasInsertionSeqType", "[EventUtils]") {
  unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>>
      events_map;

  class MockInsertionEvent : public MockEvent {
  public:
    MockInsertionEvent(string name) : MockEvent(name) {}
  };

  auto vd_ins = make_shared<MockInsertionEvent>("VD_ins");
  auto vj_ins = make_shared<MockInsertionEvent>("VJ_ins");
  events_map[make_tuple(Insertion_t, VD_genes, Undefined_side)] = vd_ins;

  SECTION("Finds present insertion segments") {
    REQUIRE(has_insertion_seq_type(events_map, VD_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, DJ_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, VJ_ins_seq));
  }

  SECTION("Tracks a different insertion segment independently") {
    events_map[make_tuple(Insertion_t, VJ_genes, Undefined_side)] = vj_ins;
    REQUIRE(has_insertion_seq_type(events_map, VD_ins_seq));
    REQUIRE(has_insertion_seq_type(events_map, VJ_ins_seq));
    REQUIRE_FALSE(has_insertion_seq_type(events_map, DJ_ins_seq));
  }
}

TEST_CASE("EventUtils TryGetEvent", "[EventUtils]") {
  unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>>
      events_map;

  auto v_event = make_shared<MockEvent>("V_choice");
  events_map[make_tuple(GeneChoice_t, V_gene, Undefined_side)] = v_event;

  SECTION("Returns a present event") {
    shared_ptr<Rec_Event> event_ptr;
    REQUIRE(try_get_event(events_map, GeneChoice_t, V_gene, Undefined_side, event_ptr));
    REQUIRE(event_ptr == v_event);
  }

  SECTION("Returns false for a missing event") {
    shared_ptr<Rec_Event> event_ptr;
    REQUIRE_FALSE(try_get_event(events_map, GeneChoice_t, J_gene, Undefined_side, event_ptr));
    REQUIRE(event_ptr == nullptr);
  }

  SECTION("Supports fallback lookup pattern") {
    auto vdj_event = make_shared<MockEvent>("VDJ_dinuc");
    events_map[make_tuple(Dinuclmarkov_t, VDJ_genes, Undefined_side)] = vdj_event;

    shared_ptr<Rec_Event> event_ptr;
    bool found_primary =
        try_get_event(events_map, Dinuclmarkov_t, VD_genes, Undefined_side, event_ptr);
    if (!found_primary) {
      found_primary =
          try_get_event(events_map, Dinuclmarkov_t, VDJ_genes, Undefined_side, event_ptr);
    }

    REQUIRE(found_primary);
    REQUIRE(event_ptr == vdj_event);
  }
}

TEST_CASE("EventUtils TryEventKeyToSeqKey", "[EventUtils]") {
  tuple<Event_type, Seq_type, Seq_side> seq_key;

  SECTION("Maps gene events to gene sequence keys") {
    REQUIRE(try_event_key_to_seq_key(GeneChoice_t, V_gene, Undefined_side, seq_key));
    REQUIRE(seq_key == make_tuple(GeneChoice_t, V_gene_seq, Undefined_side));

    REQUIRE(try_event_key_to_seq_key(Deletion_t, J_gene, Five_prime, seq_key));
    REQUIRE(seq_key == make_tuple(Deletion_t, J_gene_seq, Five_prime));
  }

  SECTION("Maps insertion-like events to insertion sequence keys") {
    REQUIRE(try_event_key_to_seq_key(Insertion_t, VD_genes, Undefined_side, seq_key));
    REQUIRE(seq_key == make_tuple(Insertion_t, VD_ins_seq, Undefined_side));

    REQUIRE(try_event_key_to_seq_key(Dinuclmarkov_t, DJ_genes, Three_prime, seq_key));
    REQUIRE(seq_key == make_tuple(Dinuclmarkov_t, DJ_ins_seq, Three_prime));
  }

  SECTION("Fails on incompatible class/type combinations") {
    REQUIRE_FALSE(try_event_key_to_seq_key(Insertion_t, V_gene, Undefined_side, seq_key));
    REQUIRE_FALSE(try_event_key_to_seq_key(GeneChoice_t, VD_genes, Undefined_side, seq_key));
    REQUIRE_FALSE(try_event_key_to_seq_key(Undefined_t, Undefined_gene, Undefined_side, seq_key));
  }
}
