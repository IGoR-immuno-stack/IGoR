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
 * Written against the *unmodified* `draw_random_realization()` overrides (G0), before they were
 * split into a draw and a construction (docs/GENERATION_REWRITE_PLAN.md, G1) and before the
 * V/D/J switches go (G4). The cases on the split itself were added with it. Follows ITERATE_TEST_GUIDE.md's rules: one TEST_CASE per pattern, one
 * SECTION per instance, numbers read off the running code, and a confirmed defect asserted at
 * its intended value under `[!shouldfail]`.
 *
 * Every draw puts all of an event's mass on one realization, so that what is chosen does not
 * depend on the order an `unordered_map` walks the realizations in. The Markov chain gets a
 * deterministic transition table (A->C->G->T->A) for the same reason. What the RNG is used for
 * is pinned separately, by comparing the generator against a copy advanced by the number of
 * uniforms the step should have consumed: that count is part of what makes a seed reproduce.
 *
 * Segments are addressed by seq_type *name* through the Generator fixture, never by the Seq_type
 * enum, so that the bodies below survive the container change (G3) and so that the tandem-D
 * cases can state what they expect before the code can deliver it.
 */

#include "test_utils.h"

#include <igor/Model/Legacy/EventUtils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

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
 * A GenerationState, plus everything a draw reads.
 *
 * The only place that knows how the generator stores its segments. Everything else goes through
 * the seq_type name, so a change of container changes this class and nothing below it (G3 did).
 */
class Generator {
public:
    explicit Generator(const SeqTypeRegistry &registry = vdj_seq_type_registry())
        : registry_(registry), segments_(registry)
    {}

    Marginal_array_p marginals{new long double[kMarginalSize]()};
    std::unordered_map<Rec_Event_name, int> index_map;
    OffsetMap offset_map;
    std::mt19937_64 rng{kSeed};

    /// Write a segment as an upstream event would have left it.
    void preset(const std::string &seq_type, const std::string &content)
    {
        segments_.create(registry_.id(seq_type), content);
    }

    bool has(const std::string &seq_type) const { return segments_.built(registry_.id(seq_type)); }

    std::string segment(const std::string &seq_type) const { return segments_.read(registry_.id(seq_type)); }

    /// Number of segments built so far, whatever their names.
    std::size_t written() const
    {
        std::size_t count = 0;
        for (std::size_t id = 0; id != registry_.total_count(); ++id) {
            count += segments_.built(static_cast<SeqTypeId>(id)) ? 1 : 0;
        }
        return count;
    }

    /// The segments in the registry's 5'->3' order, skipping those nobody built.
    std::string assemble() const { return segments_.assemble(); }

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
        return event.draw_random_realization(marginals, index_map, offset_map, segments_, rng);
    }

    /// The draw alone: indices, and nothing written.
    std::vector<int> draw_indices(const Rec_Event &event)
    {
        return event.draw_realization(marginals, index_map, segments_, rng);
    }

    /// The construction alone, from indices however they were obtained.
    void construct(const Rec_Event &event, const std::vector<int> &indices)
    {
        event.construct_realization(indices, segments_);
    }

private:
    const SeqTypeRegistry &registry_;
    GenerationState segments_;
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

/// A random generator in the state Generator's starts in, advanced by `uniforms` draws.
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

/// One-realization deletion of `value` on `seq_type`'s `side`, drawn on `content`. The event is
/// placed in a VDJ ordering, as Model_Parms::finalize() places every event: where its segment
/// sits is what generation reads its anchoring from.
std::string generated_trim(Seq_type target, Seq_side side, int value, const std::string &content)
{
    Generator state;
    const std::string name = EventUtils::seq_type_to_string(target);
    auto deletion = place(make_deletion(target, side, value, value, /*id=*/0), vdj_seq_type_registry(), name);
    state.preset(name, content);
    state.put_mass(*deletion, 0, 0);
    state.draw(*deletion);
    return state.segment(name);
}

/// A Markov chain over A,C,G,T whose every transition is certain: A->C, C->G, G->T, T->A. Its
/// ambiguous rows follow from update_event_internal_probas(), as in a real run.
std::shared_ptr<Dinucl_markov> cyclic_chain(Generator &state, Seq_type target,
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
    Generator state;
    auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}, {"V3", "GGGG"}}, /*id=*/0);

    SECTION("whichever realization holds the mass is the one drawn")
    {
        for (const std::string name : {"V1", "V2", "V3"}) {
            Generator fresh;
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
    Generator state;
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

TEST_CASE("Generation: the draw writes nothing, and the construction draws nothing", "[generation]")
{
    //The split of G1. Each half is what a sampling engine and an event's `apply` will be: the
    //draw reads the marginals, the index map, the RNG and (for a chain) what has been built;
    //the construction reads the indices and nothing else.
    SECTION("gene choice")
    {
        Generator state;
        auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}}, /*id=*/0);
        state.put_mass(*v, 0, index_of_name(*v, "V2"));
        CHECK(state.draw_indices(*v) == std::vector<int>{index_of_name(*v, "V2")});
        CHECK(state.written() == 0);
        CHECK(state.rng == advanced_by(1));

        //Indices from anywhere -- here, not the ones the marginals favour.
        state.construct(*v, {index_of_name(*v, "V1")});
        CHECK(state.segment("V_gene_seq") == "AAAA");
        CHECK(state.rng == advanced_by(1));
        CHECK(state.index_map.at(v->get_name()) == 0);
    }

    SECTION("deletion")
    {
        Generator state;
        auto deletion = make_deletion(V_gene_seq, Three_prime, 0, 3, /*id=*/0);
        state.preset("V_gene_seq", "AACCGG");
        state.put_mass(*deletion, 0, index_of_value(*deletion, 2));
        CHECK(state.draw_indices(*deletion) == std::vector<int>{index_of_value(*deletion, 2)});
        CHECK(state.segment("V_gene_seq") == "AACCGG");
        state.construct(*deletion, {index_of_value(*deletion, 3)});
        CHECK(state.segment("V_gene_seq") == "AAC");
    }

    SECTION("insertion")
    {
        Generator state;
        auto insertion = make_insertion(DJ_ins_seq, 0, 5, /*id=*/0);
        state.put_mass(*insertion, 0, index_of_value(*insertion, 4));
        CHECK(state.draw_indices(*insertion) == std::vector<int>{index_of_value(*insertion, 4)});
        CHECK(state.written() == 0);
        state.construct(*insertion, {index_of_value(*insertion, 1)});
        CHECK(state.segment("DJ_ins_seq") == "I");
    }

    SECTION("Markov chain")
    {
        Generator state;
        auto chain = cyclic_chain(state, DJ_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("J_gene_seq", "CAAA");
        state.preset("DJ_ins_seq", "III");
        CHECK(state.draw_indices(*chain) == std::vector<int>{2, 3, 0});
        CHECK(state.segment("DJ_ins_seq") == "III");
        CHECK(state.rng == advanced_by(3));

        //A chain the table would never produce: construction writes what it is given.
        state.construct(*chain, {0, 0, 1});
        CHECK(state.segment("DJ_ins_seq") == "CAA");
        CHECK(state.rng == advanced_by(3));
    }

    SECTION("no realization drawn builds nothing")
    {
        Generator state;
        auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}}, /*id=*/0);
        state.construct(*v, {});
        CHECK(state.written() == 0);
    }
}

// Plan D4: the categorical walk stops at the first realization whose running sum is >= u.
// Confirmed defects, both pinned through pick_realization(), the walk as a function of u.

TEST_CASE("The categorical walk never draws a zero-mass realization", "[generation][!shouldfail]")
{
    //At u = 0 the first realization in walk order satisfies `0 >= 0`, whatever its mass. With
    //the mass on each realization in turn, at least two of the three have a zero-mass one
    //walked before them.
    auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}, {"V3", "GGGG"}}, /*id=*/0);
    for (const std::string name : {"V1", "V2", "V3"}) {
        Generator state;
        state.put_mass(*v, 0, index_of_name(*v, name));
        INFO("mass on " << name);
        CHECK(v->pick_realization(state.marginals, 0, 0.0) == index_of_name(*v, name));
    }
}

TEST_CASE("A row a rounding error short of the uniform still draws", "[generation][!shouldfail]")
{
    //A normalized row can sum to just under 1, and a uniform in [0, 1) can exceed the sum.
    //The walk then ends having chosen nothing: the event writes nothing, records "()", and the
    //next event to read its segment throws (G0 met this for real, through a misordered
    //marginals file).
    auto v = make_gene_choice(V_gene, {{"V1", "AAAA"}, {"V2", "CCCC"}}, /*id=*/0);
    Generator state;
    state.put_mass(*v, 0, index_of_name(*v, "V2"));
    state.marginals[index_of_name(*v, "V2")] = 1.0L - 1e-12L;
    CHECK(v->pick_realization(state.marginals, 0, 1.0 - 1e-13) == index_of_name(*v, "V2"));
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
            Generator state;
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
    SECTION("a segment on an end of the sequence keeps at least one nucleotide, on either side")
    {
        //Plan D3, through D2's refusal: inference keeps one nucleotide of a segment anchored on
        //the read, and drops the realization that would delete it away.
        CHECK(generated_trim(V_gene_seq, Three_prime, 7, "AACCGGTT") == "A");
        CHECK(generated_trim(J_gene_seq, Five_prime, 7, "AACCGGTT") == "T");
        for (const auto &[target, side] : {std::pair{V_gene_seq, Three_prime}, std::pair{V_gene_seq, Five_prime},
                                           std::pair{J_gene_seq, Five_prime}, std::pair{J_gene_seq, Three_prime}}) {
            CHECK_THROWS_WITH(generated_trim(target, side, 8, "AACCGGTT"),
                              Catch::Matchers::ContainsSubstring("keeps at least one nucleotide"));
        }
    }
    SECTION("more than the segment holds is refused, on either side")
    {
        //Plan D2: inference drops the realization. Generation used to clamp a 5' trim
        //(erase(0, n) stops at the end) and throw from erase() on a 3' one; both now refuse it.
        for (const auto &[target, side] : {std::pair{D_gene_seq, Five_prime}, std::pair{D_gene_seq, Three_prime},
                                           std::pair{V_gene_seq, Three_prime}, std::pair{J_gene_seq, Five_prime}}) {
            CHECK_THROWS_WITH(generated_trim(target, side, 9, "AACCGGTT"),
                              Catch::Matchers::ContainsSubstring("run past the end"));
        }
    }
    SECTION("trimming a segment nobody drew throws")
    {
        Generator state;
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
    SECTION("as long as the template it mirrors")
    {
        CHECK(generated_trim(D_gene_seq, Five_prime, -3, "GAC") == "GTCGAC");
        CHECK(generated_trim(D_gene_seq, Three_prime, -3, "GAC") == "GACGTC");
    }
    SECTION("longer than the template it mirrors is refused, on either side")
    {
        //Plan D2: inference drops it. Generation used to mirror the whole template on a 5' end
        //(substr(0, n) stops at the end) and throw from substr() on a 3' one.
        CHECK_THROWS_WITH(generated_trim(D_gene_seq, Five_prime, -4, "GAC"),
                          Catch::Matchers::ContainsSubstring("longer than the 3 nucleotides it mirrors"));
        CHECK_THROWS_WITH(generated_trim(D_gene_seq, Three_prime, -4, "GAC"),
                          Catch::Matchers::ContainsSubstring("longer than the 3 nucleotides it mirrors"));
    }
    SECTION("a template outside ACGT cannot be mirrored")
    {
        CHECK_THROWS_AS(generated_trim(V_gene_seq, Three_prime, -2, "AAANN"), std::runtime_error);
    }
}

TEST_CASE("Generation agrees with inference on a trimmed segment", "[generation]")
{
    //The arms a shipped model has, on segments placed as an aligner would place them, with no
    //neighbour chosen. Every value inference keeps comes out of generation identical, and every
    //value it drops -- here, only for want of template: plan D2 and D3 -- generation refuses.
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
            //Each arm runs from a palindrome one past its template to a deletion one past it.
            const int length = arm.three - arm.five + 1;
            const bool anchored = (arm.target != D_gene_seq);
            int kept = 0;
            for (int value = -(length + 1); value <= length + 1; ++value) {
                INFO("deletion " << value);
                const auto inferred = inferred_trim(arm.target, arm.side, value, arm.five, arm.three);
                const std::string template_run = read_run(arm.five, arm.three);
                if (inferred) {
                    ++kept;
                    CHECK(generated_trim(arm.target, arm.side, value, template_run) == *inferred);
                } else {
                    CHECK_THROWS_AS(generated_trim(arm.target, arm.side, value, template_run), std::out_of_range);
                }
            }
            //Every palindrome up to the template's length, and every deletion up to it, but the
            //last for a segment on an end of the sequence.
            CHECK(kept == 2 * length + (anchored ? 0 : 1));
        }
    }

    SECTION("inference trims a V 5' deletion at the 5' end and a J 3' one at the 3' end")
    {
        //The premise of the two D1 cases below: the side the model names is the side
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

// Plan D1: the V and J arms of the generator never read the event's side, so a V 5' or a J 3'
// deletion was trimmed at the end inference does not trim (the premise above). Fixed: a
// deletion trims the end its side names, whatever the segment. No model has either deletion,
// so no output moved.

TEST_CASE("Generation trims a V 5' deletion at the 5' end", "[generation]")
{
    CHECK(generated_trim(V_gene_seq, Five_prime, 2, "AACCGGTT") == "CCGGTT");
}

TEST_CASE("Generation trims a J 3' deletion at the 3' end", "[generation]")
{
    CHECK(generated_trim(J_gene_seq, Three_prime, 2, "AACCGGTT") == "AACCGG");
}

TEST_CASE("Generation: an insertion writes one placeholder per inserted nucleotide", "[generation]")
{
    for (const Seq_type target : {VD_ins_seq, DJ_ins_seq, VJ_ins_seq}) {
        const std::string name = EventUtils::seq_type_to_string(target);
        DYNAMIC_SECTION(name)
        {
            Generator state;
            auto insertion = make_insertion(target, 0, 5, /*id=*/0);
            state.put_mass(*insertion, 0, index_of_value(*insertion, 3));
            CHECK(as_vector(state.draw(*insertion)) == std::vector<int>{index_of_value(*insertion, 3)});
            CHECK(state.segment(name) == "III");
        }
    }

    SECTION("no insertion writes an empty segment, which is not the same as none")
    {
        Generator state;
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
        Generator state;
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
        Generator state;
        auto chain = cyclic_chain(state, DJ_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("J_gene_seq", "CAAA");
        state.preset("DJ_ins_seq", "III");
        CHECK(as_vector(state.draw(*chain)) == std::vector<int>{2, 3, 0});
        CHECK(state.segment("DJ_ins_seq") == "ATG");
    }

    SECTION("a VJ junction, seeded from V")
    {
        Generator state(vj_seq_type_registry());
        auto chain = cyclic_chain(state, VJ_ins_seq, vj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "TTTG");
        state.preset("VJ_ins_seq", "III");
        state.draw(*chain);
        CHECK(state.segment("VJ_ins_seq") == "TAC");
    }

    SECTION("an empty insertion draws nothing and consumes nothing")
    {
        Generator state;
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
        Generator state;
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
        Generator state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "");
        state.preset("VD_ins_seq", "II");
        CHECK_THROWS_AS(state.draw(*chain), std::out_of_range);
    }

    SECTION("a chain whose insertion was not drawn throws")
    {
        Generator state;
        auto chain = cyclic_chain(state, VD_ins_seq, vdj_seq_type_registry(), 0);
        state.preset("V_gene_seq", "GGGA");
        CHECK_THROWS_AS(state.draw(*chain), std::out_of_range);
    }
}

// =======================================================================================
// Tandem D: what the generator should do on a V-D1-D2-J layout. Each case was written in G0,
// failing on the Seq_type enum, and its tag came off in the commit that made it pass: gene
// choice and insertion in G4a; deletion in G4b; the chain and the end-to-end case in G4c.
// =======================================================================================

TEST_CASE("Generation on a tandem-D layout: each D slot's gene choice writes its own segment",
          "[generation][tandem_d]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    Generator state(registry);
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
          "[generation][tandem_d]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    Generator state(registry);
    auto insertion = place(make_insertion(VD_ins_seq, 0, 3, /*id=*/0), registry, "D1D2_ins_seq");
    state.put_mass(*insertion, 0, index_of_value(*insertion, 2));
    state.draw(*insertion);
    CHECK(state.segment("D1D2_ins_seq") == "II");
}

TEST_CASE("Generation on a tandem-D layout: a deletion trims its own D slot",
          "[generation][tandem_d]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    Generator state(registry);
    auto deletion = place(make_deletion(D_gene_seq, Five_prime, 2, 2, /*id=*/0), registry, "D2_gene_seq");
    state.preset("D1_gene_seq", "GGGTTT");
    state.preset("D2_gene_seq", "CCCAAA");
    state.put_mass(*deletion, 0, 0);
    state.draw(*deletion);
    CHECK(state.segment("D2_gene_seq") == "CAAA");
    CHECK(state.segment("D1_gene_seq") == "GGGTTT");
}

TEST_CASE("Generation on a tandem-D layout: the D1-D2 chain seeds from D1's 3' end",
          "[generation][tandem_d]")
{
    const SeqTypeRegistry &registry = tandem_registry();
    Generator state(registry);
    auto chain = cyclic_chain(state, VD_ins_seq, registry, 0, "D1D2_ins_seq");
    state.preset("D1_gene_seq", "GGGTTA");
    state.preset("D1D2_ins_seq", "II");
    state.draw(*chain);
    CHECK(state.segment("D1D2_ins_seq") == "CG");
}

TEST_CASE("Generation on a tandem-D layout: one sequence, end to end", "[generation][tandem_d]")
{
    //Every event a milestone-1 model has, each with its mass on one realization, drawn in a
    //legacy queue order: genes, then deletions, then insertions, then the chains. Priorities
    //are distinct so that no two events share a generated name (T2).
    const SeqTypeRegistry &registry = tandem_registry();
    Generator state(registry);
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
