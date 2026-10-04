/*
 * test_generation_construction.cpp
 *
 *  Characterization tests for the legacy generator's per-event step (generation plan, G0).
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
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

/**
 * Written against the *unmodified* `draw_random_realization()` overrides, before they are split
 * into a draw and a construction (docs/GENERATION_REWRITE_PLAN.md, G1) and before the V/D/J
 * switches go (G4). Follows ITERATE_TEST_GUIDE.md's rules: one TEST_CASE per pattern, one
 * SECTION per instance, numbers read off the running code, and a confirmed defect asserted at
 * its intended value under `[!shouldfail]`.
 *
 * Every draw puts all of an event's mass on one realization, so that what is chosen does not
 * depend on the order an `unordered_map` walks the realizations in. The Markov chain gets a
 * deterministic transition table (A->C->G->T->A) for the same reason. What the RNG is used for
 * is pinned separately, by comparing the generator against a copy advanced by the number of
 * uniforms the step should have consumed: that count is part of what makes a seed reproduce.
 *
 * Segments are addressed by seq_type *name* through GenerationState, never by the Seq_type
 * enum, so that the bodies below survive the container change (G3) and so that the tandem-D
 * cases can state what they expect before the code can deliver it.
 */

#include "test_utils.h"

#include <igor/Model/Legacy/EventUtils.h>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace IgorTestUtils;

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kMarginalSize = 1024;

using OffsetMap =
        std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>;

/**
 * What one generated sequence is built in, plus everything a draw reads.
 *
 * The only place that knows the generator's container. Everything else goes through the
 * seq_type name, so a change of container changes this class and nothing below it.
 */
class GenerationState {
public:
    explicit GenerationState(const SeqTypeRegistry &registry = vdj_seq_type_registry()) : registry_(registry) {}

    Marginal_array_p marginals{new long double[kMarginalSize]()};
    std::unordered_map<Rec_Event_name, int> index_map;
    OffsetMap offset_map;
    std::mt19937_64 rng{kSeed};

    /// Write a segment as an upstream event would have left it.
    void preset(const std::string &seq_type, const std::string &content)
    {
        sequences_[str2SeqType(seq_type)] = content;
    }

    bool has(const std::string &seq_type) const { return sequences_.count(str2SeqType(seq_type)) != 0; }

    std::string segment(const std::string &seq_type) const { return sequences_.at(str2SeqType(seq_type)); }

    /// Number of segments written so far, whatever their names.
    std::size_t written() const { return sequences_.size(); }

    /// The segments in the registry's 5'->3' order, skipping those nobody wrote.
    std::string assemble() const
    {
        std::string out;
        for (const auto &name : registry_.get_ordered_types()) {
            if (has(name)) {
                out += segment(name);
            }
        }
        return out;
    }

    /// All of `event`'s mass on realization `index`, in a block starting at `base`.
    void put_mass(const Rec_Event &event, int base, int index)
    {
        index_map[event.get_name()] = base;
        for (int i = 0; i != event.size(); ++i) {
            marginals[base + i] = 0.0L;
        }
        marginals[base + index] = 1.0L;
    }

    std::queue<int> draw(const Rec_Event &event)
    {
        return event.draw_random_realization(marginals, index_map, offset_map, sequences_, rng);
    }

private:
    const SeqTypeRegistry &registry_;
    std::unordered_map<Seq_type, std::string> sequences_;
};

std::vector<int> as_vector(std::queue<int> queue)
{
    std::vector<int> out;
    while (!queue.empty()) {
        out.push_back(queue.front());
        queue.pop();
    }
    return out;
}

/// A generator in the state GenerationState's starts in, advanced by `uniforms` draws.
std::mt19937_64 advanced_by(int uniforms)
{
    std::mt19937_64 rng{kSeed};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int i = 0; i != uniforms; ++i) {
        uniform(rng);
    }
    return rng;
}

int index_of_value(const Rec_Event &event, int value)
{
    for (const auto &[key, realization] : event.get_realizations_map()) {
        if (realization.value_int == value) {
            return realization.index;
        }
    }
    throw std::invalid_argument("no realization with value " + std::to_string(value));
}

int index_of_name(const Rec_Event &event, const std::string &name)
{
    return event.get_realizations_map().at(name).index;
}

/// Tell an event where it sits, as Model_Parms::finalize() does: its seq_type, the id the
/// registry gives that name, and the ids on either side of it in the ordering.
template <class Event>
std::shared_ptr<Event> place(std::shared_ptr<Event> event, const SeqTypeRegistry &registry,
                             const std::string &seq_type)
{
    event->set_seq_type(seq_type);
    const SeqTypeId id = registry.id(seq_type);
    event->set_seq_type_id(id);
    event->set_adjacent_segments(registry.left_neighbor(id), registry.right_neighbor(id));
    event->update_event_name();
    return event;
}

/// One-realization deletion of `value` on `seq_type`'s `side`, drawn on `content`.
std::string generated_trim(Seq_type target, Seq_side side, int value, const std::string &content)
{
    GenerationState state;
    auto deletion = make_deletion(target, side, value, value, /*id=*/0);
    const std::string name = EventUtils::seq_type_to_string(target);
    state.preset(name, content);
    state.put_mass(*deletion, 0, 0);
    state.draw(*deletion);
    return state.segment(name);
}

/// A Markov chain over A,C,G,T whose every transition is certain: A->C, C->G, G->T, T->A. Its
/// ambiguous rows follow from update_event_internal_probas(), as in a real run.
std::shared_ptr<Dinucl_markov> cyclic_chain(GenerationState &state, Seq_type target,
                                            const SeqTypeRegistry &registry, int base,
                                            const std::string &seq_type = "")
{
    auto chain = make_dinucl_markov(target, /*id=*/9);
    place(chain, registry, seq_type.empty() ? EventUtils::seq_type_to_string(target) : seq_type);
    state.index_map[chain->get_name()] = base;
    for (int previous = 0; previous != 4; ++previous) {
        for (int next = 0; next != 4; ++next) {
            state.marginals[base + previous * 4 + next] = (next == (previous + 1) % 4) ? 1.0L : 0.0L;
        }
    }
    chain->update_event_internal_probas(state.marginals, state.index_map);
    return chain;
}

/// V, two D slots and J: the milestone-1 layout.
const SeqTypeRegistry &tandem_registry()
{
    static const SeqTypeRegistry registry = [] {
        SeqTypeRegistry r;
        r.register_legacy_seq_types();
        r.set_ordered_types({"V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq", "D2_gene_seq",
                             "D2J_ins_seq", "J_gene_seq"});
        r.freeze();
        return r;
    }();
    return registry;
}

// For the comparison with inference: the read the segments are cut from, as in
// test_deletion_iterate.cpp, so that a palindrome is scored against what it mirrors.
const std::string kRead = "ACGTACGTACGTACGTACGTACGTACGT";

std::string read_run(Seq_Offset five, Seq_Offset three)
{
    if (three < five) { return {}; }
    return kRead.substr(static_cast<std::size_t>(five), static_cast<std::size_t>(three - five + 1));
}

/// What Deletion::iterate leaves in the segment for one deletion value, or nothing when it
/// drops the realization. The segment sits at [five, three] of kRead, with no neighbour chosen.
std::optional<std::string> inferred_trim(Seq_type target, Seq_side side, int value, Seq_Offset five,
                                         Seq_Offset three)
{
    IterateTestState state = create_iterate_state(kRead);
    auto deletion = make_deletion(target, side, value, value, /*id=*/0);
    state.preset_segment(target, five, three, read_run(five, three));
    for (std::size_t i = 0; i != 64; ++i) { state.set_marginal(i, 0.5L); }
    auto recorder = call_iterate_recording(deletion, state);
    if (recorder->call_count() == 0) {
        return std::nullopt;
    }
    return recorder->calls.at(0).sequences.at(target);
}

} // namespace

// =======================================================================================
// The draw: shared by Gene_choice, Deletion and Insertion
// =======================================================================================

TEST_CASE("Generation: a categorical draw follows the event's marginal row", "[generation]")
{
    GenerationState state;
    auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}, {"V3", "GGGG"}}, /*id=*/0);

    SECTION("whichever realization holds the mass is the one drawn")
    {
        for (const std::string name : {"V1", "V2", "V3"}) {
            GenerationState fresh;
            fresh.put_mass(*v, 0, index_of_name(*v, name));
            CHECK(as_vector(fresh.draw(*v)) == std::vector<int>{index_of_name(*v, name)});
        }
    }

    SECTION("the row starts at the event's entry in the index map")
    {
        //Mass at the start of the array would be read by an event that ignored its base.
        state.put_mass(*v, 100, index_of_name(*v, "V3"));
        state.marginals[0] = 1.0L;
        CHECK(as_vector(state.draw(*v)) == std::vector<int>{index_of_name(*v, "V3")});
    }

    SECTION("one uniform per draw, for each of the three scalar event types")
    {
        state.put_mass(*v, 0, 0);
        state.preset("V_gene_seq", "AAAA");
        state.draw(*v);
        CHECK(state.rng == advanced_by(1));

        auto deletion = make_deletion(V_gene_seq, Three_prime, 0, 2, /*id=*/1);
        state.put_mass(*deletion, 10, 1);
        state.draw(*deletion);
        CHECK(state.rng == advanced_by(2));

        auto insertion = make_insertion(VD_ins_seq, 0, 3, /*id=*/2);
        state.put_mass(*insertion, 20, 2);
        state.draw(*insertion);
        CHECK(state.rng == advanced_by(3));
    }
}

TEST_CASE("Generation: a drawn realization moves its children's rows", "[generation]")
{
    GenerationState state;
    auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}, {"V3", "GGGG"}}, /*id=*/0);
    auto j = make_gene_choice(J_gene, {{"J1", "TTTT"}, {"J2", "GGTT"}}, /*id=*/1);
    auto v_del = make_deletion(V_gene_seq, Three_prime, 0, 4, /*id=*/2);

    SECTION("each child moves by the realization's index times its stride")
    {
        state.put_mass(*v, 0, index_of_name(*v, "V3"));
        state.index_map[j->get_name()] = 100;
        state.index_map[v_del->get_name()] = 200;
        state.offset_map[v->get_name()] = {{j, 2}, {v_del, 5}};
        state.draw(*v);
        CHECK(state.index_map.at(j->get_name()) == 100 + index_of_name(*v, "V3") * 2);
        CHECK(state.index_map.at(v_del->get_name()) == 200 + index_of_name(*v, "V3") * 5);
    }

    SECTION("an event with no children moves nothing")
    {
        state.put_mass(*v, 0, index_of_name(*v, "V2"));
        state.index_map[j->get_name()] = 100;
        state.draw(*v);
        CHECK(state.index_map.at(j->get_name()) == 100);
        CHECK(state.index_map.at(v->get_name()) == 0);
    }

    SECTION("the child then draws from the row its parent selected")
    {
        //J's block holds one row of two per V realization; V3's row puts all its mass on J2.
        state.put_mass(*v, 0, index_of_name(*v, "V3"));
        state.index_map[j->get_name()] = 100;
        state.offset_map[v->get_name()] = {{j, 2}};
        for (int i = 0; i != 6; ++i) { state.marginals[100 + i] = 0.0L; }
        state.marginals[100 + index_of_name(*v, "V3") * 2 + index_of_name(*j, "J2")] = 1.0L;
        state.draw(*v);
        state.draw(*j);
        CHECK(state.segment("J_gene_seq") == "GGTT");
    }
}

// =======================================================================================
// Construction: what each event type writes
// =======================================================================================

TEST_CASE("Generation: a gene choice writes its template on its segment", "[generation]")
{
    struct Instance {
        Gene_class gene_class;
        const char *seq_type;
    };
    for (const Instance instance : {Instance{V_gene, "V_gene_seq"}, Instance{D_gene, "D_gene_seq"},
                                    Instance{J_gene, "J_gene_seq"}}) {
        DYNAMIC_SECTION(instance.seq_type)
        {
            GenerationState state;
            auto gene = make_gene_choice(instance.gene_class, {{"G1", "ACGTTGCA"}, {"G2", "TTTT"}}, /*id=*/0);
            state.put_mass(*gene, 0, index_of_name(*gene, "G1"));
            state.draw(*gene);
            CHECK(state.segment(instance.seq_type) == "ACGTTGCA");
            CHECK(state.written() == 1);
        }
    }
}

TEST_CASE("Generation: a deletion trims its segment", "[generation]")
{
    SECTION("V 3'") { CHECK(generated_trim(V_gene_seq, Three_prime, 3, "AACCGGTT") == "AACCG"); }
    SECTION("J 5'") { CHECK(generated_trim(J_gene_seq, Five_prime, 3, "AACCGGTT") == "CGGTT"); }
    SECTION("D 5'") { CHECK(generated_trim(D_gene_seq, Five_prime, 2, "AACCGGTT") == "CCGGTT"); }
    SECTION("D 3'") { CHECK(generated_trim(D_gene_seq, Three_prime, 2, "AACCGGTT") == "AACCGG"); }
    SECTION("zero leaves the segment as it was")
    {
        CHECK(generated_trim(V_gene_seq, Three_prime, 0, "AACCGGTT") == "AACCGGTT");
        CHECK(generated_trim(D_gene_seq, Five_prime, 0, "AACCGGTT") == "AACCGGTT");
    }
    SECTION("an internal segment may be deleted away, leaving it written and empty")
    {
        CHECK(generated_trim(D_gene_seq, Five_prime, 8, "AACCGGTT").empty());
        CHECK(generated_trim(D_gene_seq, Three_prime, 8, "AACCGGTT").empty());
    }
    SECTION("an anchored segment may be deleted away too")
    {
        //Observed, intent undecided (plan D3). Inference keeps at least one nucleotide of V and
        //J -- "Generation agrees with inference" below shows it dropping the realization -- so
        //a model whose marginals put mass here generates a sequence inference calls impossible.
        //An inferred model never does: the scenario is never counted.
        CHECK(generated_trim(V_gene_seq, Three_prime, 8, "AACCGGTT").empty());
        CHECK(generated_trim(J_gene_seq, Five_prime, 8, "AACCGGTT").empty());
    }
    SECTION("more than the segment holds")
    {
        //Observed, intent undecided (plan D2). Inference drops the realization; generation
        //clamps a 5' trim (erase(0, n) stops at the end) and throws on a 3' trim (erase(size - n)
        //wraps around). Which it does depends only on the side.
        CHECK(generated_trim(D_gene_seq, Five_prime, 10, "AACCGGTT").empty());
        CHECK(generated_trim(J_gene_seq, Five_prime, 10, "AACCGGTT").empty());
        CHECK_THROWS_AS(generated_trim(D_gene_seq, Three_prime, 10, "AACCGGTT"), std::out_of_range);
        CHECK_THROWS_AS(generated_trim(V_gene_seq, Three_prime, 10, "AACCGGTT"), std::out_of_range);
    }
    SECTION("trimming a segment nobody drew throws")
    {
        GenerationState state;
        auto deletion = make_deletion(D_gene_seq, Five_prime, 2, 2, /*id=*/0);
        state.put_mass(*deletion, 0, 0);
        CHECK_THROWS_AS(state.draw(*deletion), std::out_of_range);
    }
}

TEST_CASE("Generation: a negative deletion adds a palindrome", "[generation]")
{
    //The |k| nucleotides at the trimmed end, reversed and complemented, put back beyond it.
    SECTION("V 3'") { CHECK(generated_trim(V_gene_seq, Three_prime, -2, "AAACT") == "AAACTAG"); }
    SECTION("J 5'") { CHECK(generated_trim(J_gene_seq, Five_prime, -2, "GTCCC") == "ACGTCCC"); }
    SECTION("D 5'") { CHECK(generated_trim(D_gene_seq, Five_prime, -3, "GACTTT") == "GTCGACTTT"); }
    SECTION("D 3'") { CHECK(generated_trim(D_gene_seq, Three_prime, -3, "TTTGAC") == "TTTGACGTC"); }
    SECTION("longer than the template it mirrors")
    {
        //Observed, intent undecided (plan D2): inference drops it; generation mirrors the whole
        //template on a 5' end (substr(0, n) stops at the end) and throws on a 3' end.
        CHECK(generated_trim(D_gene_seq, Five_prime, -5, "GAC") == "GTCGAC");
        CHECK_THROWS_AS(generated_trim(D_gene_seq, Three_prime, -5, "GAC"), std::out_of_range);
    }
    SECTION("a template outside ACGT cannot be mirrored")
    {
        CHECK_THROWS_AS(generated_trim(V_gene_seq, Three_prime, -2, "AAANN"), std::runtime_error);
    }
}

TEST_CASE("Generation agrees with inference on a trimmed segment", "[generation]")
{
    //The arms a shipped model has, on segments placed as an aligner would place them. Every
    //value inference keeps must come out of generation identical; the values it drops are the
    //D2/D3 sections above.
    struct Arm {
        Seq_type target;
        Seq_side side;
        Seq_Offset five;
        Seq_Offset three;
        const char *label;
    };
    for (const Arm arm : {Arm{V_gene_seq, Three_prime, 0, 9, "V 3'"}, Arm{D_gene_seq, Five_prime, 12, 17, "D 5'"},
                          Arm{D_gene_seq, Three_prime, 12, 17, "D 3'"}, Arm{J_gene_seq, Five_prime, 20, 27, "J 5'"}}) {
        DYNAMIC_SECTION(arm.label)
        {
            int compared = 0;
            for (int value = -2; value <= 3; ++value) {
                const auto inferred = inferred_trim(arm.target, arm.side, value, arm.five, arm.three);
                if (!inferred) {
                    continue;
                }
                ++compared;
                INFO("deletion " << value);
                CHECK(generated_trim(arm.target, arm.side, value, read_run(arm.five, arm.three)) == *inferred);
            }
            CHECK(compared == 6);
        }
    }

    SECTION("inference trims a V 5' deletion at the 5' end and a J 3' one at the 3' end")
    {
        //The premise of the two [!shouldfail] cases below: the side the model names is the side
        //inference trims, whichever gene it is.
        CHECK(inferred_trim(V_gene_seq, Five_prime, 2, 0, 9) == read_run(2, 9));
        CHECK(inferred_trim(J_gene_seq, Three_prime, 2, 20, 27) == read_run(20, 25));
    }

    SECTION("inference keeps one nucleotide of an anchored segment, and drops a deletion past the end")
    {
        CHECK_FALSE(inferred_trim(V_gene_seq, Three_prime, 10, 0, 9).has_value());
        CHECK_FALSE(inferred_trim(D_gene_seq, Five_prime, 7, 12, 17).has_value());
    }
}

// Plan D1: the V and J arms of the generator never read the event's side. A confirmed defect --
// inference trims the side the model names (the section above), so the generator builds a
// sequence for a different scenario than the one it records. No shipped model has either
// deletion, which is why no golden output can show it.

TEST_CASE("Generation trims a V 5' deletion at the 5' end", "[generation][!shouldfail]")
{
    CHECK(generated_trim(V_gene_seq, Five_prime, 2, "AACCGGTT") == "CCGGTT");
}

TEST_CASE("Generation trims a J 3' deletion at the 3' end", "[generation][!shouldfail]")
{
    CHECK(generated_trim(J_gene_seq, Three_prime, 2, "AACCGGTT") == "AACCGG");
}

TEST_CASE("Generation: an insertion writes one placeholder per inserted nucleotide", "[generation]")
{
    for (const Seq_type target : {VD_ins_seq, DJ_ins_seq, VJ_ins_seq}) {
        const std::string name = EventUtils::seq_type_to_string(target);
        DYNAMIC_SECTION(name)
        {
            GenerationState state;
            auto insertion = make_insertion(target, 0, 5, /*id=*/0);
            state.put_mass(*insertion, 0, index_of_value(*insertion, 3));
            CHECK(as_vector(state.draw(*insertion)) == std::vector<int>{index_of_value(*insertion, 3)});
            CHECK(state.segment(name) == "III");
        }
    }

    SECTION("no insertion writes an empty segment, which is not the same as none")
    {
        GenerationState state;
        auto insertion = make_insertion(VD_ins_seq, 0, 5, /*id=*/0);
        state.put_mass(*insertion, 0, index_of_value(*insertion, 0));
        state.draw(*insertion);
        CHECK(state.has("VD_ins_seq"));
        CHECK(state.segment("VD_ins_seq").empty());
    }
}

TEST_CASE("Generation: a Markov chain fills its insertion from the anchor's end", "[generation]")
{
    SECTION("seeded from the 3' end of the segment on its left")
    {
        GenerationState state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "GGGA");
        state.preset("VD_ins_seq", "III");
        CHECK(as_vector(state.draw(*chain)) == std::vector<int>{1, 2, 3});
        CHECK(state.segment("VD_ins_seq") == "CGT");
        CHECK(state.rng == advanced_by(3));
    }

    SECTION("seeded from the 5' end of the segment on its right, and written reversed")
    {
        //The chain runs leftwards from J: C, then G, T, A. The draws are recorded in the order
        //they were made, and the segment is that chain read 5'->3'.
        GenerationState state;
        auto chain = cyclic_chain(state, DJ_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("J_gene_seq", "CAAA");
        state.preset("DJ_ins_seq", "III");
        CHECK(as_vector(state.draw(*chain)) == std::vector<int>{2, 3, 0});
        CHECK(state.segment("DJ_ins_seq") == "ATG");
    }

    SECTION("a VJ junction, seeded from V")
    {
        GenerationState state(vj_seq_type_registry());
        auto chain = cyclic_chain(state, VJ_ins_seq, vj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "TTTG");
        state.preset("VJ_ins_seq", "III");
        state.draw(*chain);
        CHECK(state.segment("VJ_ins_seq") == "TAC");
    }

    SECTION("an empty insertion draws nothing and consumes nothing")
    {
        GenerationState state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "GGGA");
        state.preset("VD_ins_seq", "");
        CHECK(as_vector(state.draw(*chain)).empty());
        CHECK(state.segment("VD_ins_seq").empty());
        CHECK(state.rng == advanced_by(0));
    }

    SECTION("only placeholders are drawn; any other character stays, and seeds the next draw")
    {
        //Never produced by an Insertion, which writes nothing but placeholders. Pinned because
        //the split (G1) has to reproduce it: a draw that returned a plain chain and a
        //construction that overwrote every position would differ exactly here.
        GenerationState state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "GGGA");
        state.preset("VD_ins_seq", "IAI");
        CHECK(as_vector(state.draw(*chain)) == std::vector<int>{1, 1});
        CHECK(state.segment("VD_ins_seq") == "CAC");
        CHECK(state.rng == advanced_by(2));
    }

    SECTION("an empty anchor throws")
    {
        //Plan D6: inference throws here too (§7.12, R2). What a chain should seed from when its
        //anchor is empty is B10's question, for both.
        GenerationState state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "");
        state.preset("VD_ins_seq", "II");
        CHECK_THROWS_AS(state.draw(*chain), std::out_of_range);
    }

    SECTION("a chain whose insertion was not drawn throws")
    {
        GenerationState state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "GGGA");
        CHECK_THROWS_AS(state.draw(*chain), std::out_of_range);
    }
}

// =======================================================================================
// Tandem D: what the generator should do on a V-D1-D2-J layout. None of it can today --
// every case below fails on the Seq_type enum -- and each tag comes off in the commit that
// makes its case pass (G4a, G4b, G4c).
// =======================================================================================

TEST_CASE("Generation on a tandem-D layout: each D slot's gene choice writes its own segment",
          "[generation][tandem_d][!shouldfail]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    GenerationState state(registry);
    auto d1 = place(make_gene_choice(D_gene, {{"D1", "GGGTTT"}}, /*id=*/0), registry, "D1_gene_seq");
    auto d2 = place(make_gene_choice(D_gene, {{"D1", "GGGTTT"}, {"D2", "CCCAAA"}}, /*id=*/1), registry,
                    "D2_gene_seq");
    state.put_mass(*d1, 0, 0);
    state.put_mass(*d2, 10, index_of_name(*d2, "D2"));
    state.draw(*d1);
    state.draw(*d2);
    CHECK(state.segment("D1_gene_seq") == "GGGTTT");
    CHECK(state.segment("D2_gene_seq") == "CCCAAA");
}

TEST_CASE("Generation on a tandem-D layout: an insertion writes its own junction",
          "[generation][tandem_d][!shouldfail]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    GenerationState state(registry);
    auto insertion = place(make_insertion(VD_ins_seq, 0, 3, /*id=*/0), registry, "D1D2_ins_seq");
    state.put_mass(*insertion, 0, index_of_value(*insertion, 2));
    state.draw(*insertion);
    CHECK(state.segment("D1D2_ins_seq") == "II");
}

TEST_CASE("Generation on a tandem-D layout: a deletion trims its own D slot",
          "[generation][tandem_d][!shouldfail]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    GenerationState state(registry);
    auto deletion = place(make_deletion(D_gene_seq, Five_prime, 2, 2, /*id=*/0), registry, "D2_gene_seq");
    state.preset("D1_gene_seq", "GGGTTT");
    state.preset("D2_gene_seq", "CCCAAA");
    state.put_mass(*deletion, 0, 0);
    state.draw(*deletion);
    CHECK(state.segment("D2_gene_seq") == "CAAA");
    CHECK(state.segment("D1_gene_seq") == "GGGTTT");
}

TEST_CASE("Generation on a tandem-D layout: the D1-D2 chain seeds from D1's 3' end",
          "[generation][tandem_d][!shouldfail]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    GenerationState state(registry);
    auto chain = cyclic_chain(state, VD_ins_seq, registry, 0, "D1D2_ins_seq");
    state.preset("D1_gene_seq", "GGGTTA");
    state.preset("D1D2_ins_seq", "II");
    state.draw(*chain);
    CHECK(state.segment("D1D2_ins_seq") == "CG");
}

TEST_CASE("Generation on a tandem-D layout: one sequence, end to end", "[generation][tandem_d][!shouldfail]")
{
    //Every event a milestone-1 model has, each with its mass on one realization, drawn in a
    //legacy queue order: genes, then deletions, then insertions, then the chains. Priorities
    //are distinct so that no two events share a generated name (T2).
    const SeqTypeRegistry &registry = tandem_registry();
    GenerationState state(registry);
    int base = 0;
    int priority = 30;
    std::vector<std::shared_ptr<Rec_Event>> queue;
    auto add = [&](std::shared_ptr<Rec_Event> event, int index) {
        event->set_priority(priority--);
        state.put_mass(*event, base, index);
        base += 64;
        queue.push_back(event);
    };
    auto v = place(make_gene_choice(V_gene, {{"V1", "AAAAAAAAGA"}}, 0), registry, "V_gene_seq");
    auto j = place(make_gene_choice(J_gene, {{"J1", "TTTTCC"}}, 1), registry, "J_gene_seq");
    auto d1 = place(make_gene_choice(D_gene, {{"D1", "GGGTTT"}}, 2), registry, "D1_gene_seq");
    auto d2 = place(make_gene_choice(D_gene, {{"D2", "CCCAAA"}}, 3), registry, "D2_gene_seq");
    add(v, 0);
    add(j, 0);
    add(d1, 0);
    add(d2, 0);
    auto v_3 = place(make_deletion(V_gene_seq, Three_prime, 2, 2, 4), registry, "V_gene_seq");
    auto d1_5 = place(make_deletion(D_gene_seq, Five_prime, 1, 1, 5), registry, "D1_gene_seq");
    auto d1_3 = place(make_deletion(D_gene_seq, Three_prime, 1, 1, 6), registry, "D1_gene_seq");
    auto d2_5 = place(make_deletion(D_gene_seq, Five_prime, 2, 2, 7), registry, "D2_gene_seq");
    auto d2_3 = place(make_deletion(D_gene_seq, Three_prime, 0, 0, 8), registry, "D2_gene_seq");
    auto j_5 = place(make_deletion(J_gene_seq, Five_prime, 0, 0, 9), registry, "J_gene_seq");
    for (const auto &deletion : {v_3, d1_5, d1_3, d2_5, d2_3, j_5}) {
        add(deletion, 0);
    }
    auto vd1 = place(make_insertion(VD_ins_seq, 2, 2, 10), registry, "VD1_ins_seq");
    auto d1d2 = place(make_insertion(VD_ins_seq, 1, 1, 11), registry, "D1D2_ins_seq");
    auto d2j = place(make_insertion(DJ_ins_seq, 2, 2, 12), registry, "D2J_ins_seq");
    for (const auto &insertion : {vd1, d1d2, d2j}) {
        add(insertion, 0);
    }
    for (const auto &[target, name] : {std::pair{VD_ins_seq, "VD1_ins_seq"}, std::pair{VD_ins_seq, "D1D2_ins_seq"},
                                       std::pair{DJ_ins_seq, "D2J_ins_seq"}}) {
        auto chain = cyclic_chain(state, target, registry, base, name);
        chain->set_priority(priority--);
        state.index_map[chain->get_name()] = base;
        chain->update_event_internal_probas(state.marginals, state.index_map);
        base += 64;
        queue.push_back(chain);
    }

    for (const auto &event : queue) {
        state.draw(*event);
    }
    //V AAAAAAAA | VD1 from A: CG | D1 GGTT | D1D2 from T: A | D2 CAAA | D2J from J's T: A, C,
    //reversed: CA | J TTTTCC
    CHECK(state.assemble() == "AAAAAAAACGGGTTACAAACATTTTCC");
}
