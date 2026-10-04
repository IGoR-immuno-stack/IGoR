/*
 * test_generation.cpp
 *
 *  Characterization tests for the legacy generator end to end (generation plan, G0).
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
 * What SequenceGenerator::generate_unique_sequence() makes of a whole model: the queue order, the
 * conditioning between events, and the assembly of the segments into one sequence. The
 * per-event construction is pinned in tst/igor/Model/Legacy/test_generation_construction.cpp.
 *
 * The small models put all their mass on one realization per event, so their output does not
 * depend on the seed. The two seeded cases pin the shipped TRA (VJ) and TRB (VDJ) models: the
 * regression track `generate` covers a VDJ model only, and nothing covered VJ generation.
 * Their expected values were read off the unmodified generator; they hold as long as the
 * realizations are walked in the same order, which synthesis step 1 is the commit to change.
 */

#include <catch2/catch_test_macros.hpp>

#include <igor/Generation/Legacy/SequenceGenerator.h>
#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Model/Legacy/Model_marginals.h>
#include <igor/Model/Legacy/Rec_Event.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace igor::core::legacy;
using namespace igor::model::legacy;
using namespace igor::generation::legacy;

namespace {

/// A scratch directory, removed when the test is done with it.
class ScratchDir {
public:
    ScratchDir()
    {
        std::random_device device;
        path_ = std::filesystem::temp_directory_path() / ("igor_generation_" + std::to_string(device()));
        std::filesystem::create_directories(path_);
    }
    ~ScratchDir()
    {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }
    ScratchDir(const ScratchDir &) = delete;
    ScratchDir &operator=(const ScratchDir &) = delete;

    std::string file(const std::string &name, const std::string &content = "") const
    {
        const std::string full = (path_ / name).string();
        if (!content.empty()) {
            std::ofstream(full) << content;
        }
        return full;
    }

private:
    std::filesystem::path path_;
};

std::string read_file(const std::string &path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

std::uint64_t fnv1a(const std::string &text)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char c : text) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    return hash;
}

/// A dinucleotide table whose every transition is certain: A->C, C->G, G->T, T->A.
const char *const kCyclicDinucl = "%0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0\n";

const char *const kDinuclRealizations = "%A;0\n%C;1\n%G;2\n%T;3\n";

/// A VDJ model where J is conditioned on V, as in the shipped TRB model.
std::string vdj_parms(bool with_chains = true)
{
    std::string parms = "@Event_list\n"
                        "#GeneChoice;V_gene;Undefined_side;7;v_choice\n%V1;AAAAAAAAGA;0\n%V2;TTTTTTTTGC;1\n"
                        "#GeneChoice;D_gene;Undefined_side;6;d_gene\n%D1;GGGTTT;0\n"
                        "#GeneChoice;J_gene;Undefined_side;7;j_choice\n%J1;TTTTCC;0\n%J2;GGGGAA;1\n"
                        "#Deletion;V_gene;Three_prime;5;v_3_del\n%0;0\n%2;1\n"
                        "#Deletion;D_gene;Three_prime;5;d_3_del\n%0;0\n%1;1\n"
                        "#Deletion;D_gene;Five_prime;5;d_5_del\n%0;0\n%1;1\n"
                        "#Deletion;J_gene;Five_prime;5;j_5_del\n%0;0\n%1;1\n"
                        "#Insertion;VD_genes;Undefined_side;4;vd_ins\n%0;0\n%2;1\n"
                        "#Insertion;DJ_gene;Undefined_side;2;dj_ins\n%0;0\n%2;1\n";
    if (with_chains) {
        parms += std::string("#DinucMarkov;VD_genes;Undefined_side;3;vd_dinucl\n") + kDinuclRealizations
                 + "#DinucMarkov;DJ_gene;Undefined_side;1;dj_dinucl\n" + kDinuclRealizations;
    }
    parms += "@Edges\n"
             "%GeneChoice_V_gene_Undefined_side_prio7_size2;GeneChoice_J_gene_Undefined_side_prio7_size2\n"
             "@ErrorRate\n#SingleErrorRate\n0\n";
    return parms;
}

/// V2, then the J row V2 selects (J2), and the second realization of everything else.
///
/// In the model queue's order: Model_marginals::txt2marginals() reads the values one after the
/// other into the flat array and does not look at the `@` headers, so a block written out of
/// order lands on another event's row.
std::string vdj_marginals(bool with_chains = true)
{
    const std::string chain = std::string("$Dim[16]\n#\n") + kCyclicDinucl;
    return std::string("@v_choice\n$Dim[2]\n#\n%0,1\n"
                       "@j_choice\n$Dim[2,2]\n#[v_choice,0]\n%1,0\n#[v_choice,1]\n%0,1\n"
                       "@d_gene\n$Dim[1]\n#\n%1\n"
                       "@v_3_del\n$Dim[2]\n#\n%0,1\n"
                       "@d_3_del\n$Dim[2]\n#\n%0,1\n"
                       "@d_5_del\n$Dim[2]\n#\n%0,1\n"
                       "@j_5_del\n$Dim[2]\n#\n%0,1\n"
                       "@vd_ins\n$Dim[2]\n#\n%0,1\n")
           + (with_chains ? "@vd_dinucl\n" + chain : "") + "@dj_ins\n$Dim[2]\n#\n%0,1\n"
           + (with_chains ? "@dj_dinucl\n" + chain : "");
}

std::string vj_parms()
{
    return std::string("@Event_list\n"
                       "#GeneChoice;V_gene;Undefined_side;7;v_choice\n%V1;AAAAAAAAGA;0\n"
                       "#GeneChoice;J_gene;Undefined_side;6;j_choice\n%J1;TTTTCC;0\n"
                       "#Deletion;V_gene;Three_prime;5;v_3_del\n%0;0\n%2;1\n"
                       "#Deletion;J_gene;Five_prime;5;j_5_del\n%0;0\n%1;1\n"
                       "#Insertion;VJ_gene;Undefined_side;4;vj_ins\n%0;0\n%3;1\n"
                       "#DinucMarkov;VJ_gene;Undefined_side;3;vj_dinucl\n")
           + kDinuclRealizations + "@Edges\n@ErrorRate\n#SingleErrorRate\n0\n";
}

std::string vj_marginals()
{
    return std::string("@v_choice\n$Dim[1]\n#\n%1\n"
                       "@j_choice\n$Dim[1]\n#\n%1\n"
                       "@v_3_del\n$Dim[2]\n#\n%0,1\n"
                       "@j_5_del\n$Dim[2]\n#\n%0,1\n"
                       "@vj_ins\n$Dim[2]\n#\n%0,1\n"
                       "@vj_dinucl\n$Dim[16]\n#\n")
           + kCyclicDinucl;
}

/// One generated sequence and its realizations, keyed by the events' nicknames.
struct Generated {
    std::string sequence;
    std::map<std::string, std::vector<int>> realizations;
};

Generated generate_one(const ScratchDir &dir, const std::string &parms_text, const std::string &marginals_text)
{
    Model_Parms parms;
    parms.read_model_parms(dir.file("parms.txt", parms_text));
    Model_marginals marginals(parms);
    marginals.txt2marginals(dir.file("marginals.txt", marginals_text), parms);
    SequenceGenerator generator(parms, marginals);

    auto sequences = generator.generate_sequences(1, /*generate_errors=*/false);
    auto &[sequence, per_event] = sequences.front();

    Generated out{sequence, {}};
    auto queue = parms.get_model_queue();
    while (!queue.empty() && !per_event.empty()) {
        std::vector<int> indices;
        for (auto realization = per_event.front(); !realization.empty(); realization.pop()) {
            indices.push_back(realization.front());
        }
        out.realizations[queue.front()->get_nickname()] = indices;
        queue.pop();
        per_event.pop();
    }
    return out;
}

/// Generate `count` sequences from a shipped model with a fixed seed, as `igor generate` does,
/// and return the two files it writes.
std::pair<std::string, std::string> generate_shipped(const std::string &model, std::size_t count, int seed)
{
    const std::string dir = std::string(IGOR_MODELS_DIR) + "/human/" + model + "/models/";
    Model_Parms parms;
    parms.read_model_parms(dir + "model_parms.txt");
    Model_marginals marginals(parms);
    marginals.txt2marginals(dir + "model_marginals.txt", parms);
    SequenceGenerator generator(parms, marginals);

    ScratchDir scratch;
    const std::string seqs = scratch.file("seqs.csv");
    const std::string reals = scratch.file("reals.csv");
    generator.generate_sequences(static_cast<int>(count), /*generate_errors=*/false, seqs, reals, {},
                                 /*output_only_func=*/false, seed);
    return {read_file(seqs), read_file(reals)};
}

} // namespace

TEST_CASE("Generation: a VDJ model, one sequence", "[generation]")
{
    ScratchDir dir;
    const Generated generated = generate_one(dir, vdj_parms(), vdj_marginals());

    //V2 trimmed by 2 | VD: A, C from V's last T | D1 trimmed by 1 on each end | DJ: T, A from J's
    //first G, written reversed | J2, the row V2 selects, trimmed by 1.
    CHECK(generated.sequence == "TTTTTTTTACGGTTATGGGAA");
    CHECK(generated.realizations.at("v_choice") == std::vector<int>{1});
    CHECK(generated.realizations.at("j_choice") == std::vector<int>{1});
    CHECK(generated.realizations.at("d_gene") == std::vector<int>{0});
    CHECK(generated.realizations.at("vd_ins") == std::vector<int>{1});
    CHECK(generated.realizations.at("vd_dinucl") == std::vector<int>{0, 1});
    CHECK(generated.realizations.at("dj_dinucl") == std::vector<int>{3, 0});
}

TEST_CASE("Generation: a VJ model, one sequence", "[generation]")
{
    ScratchDir dir;
    const Generated generated = generate_one(dir, vj_parms(), vj_marginals());

    //V trimmed by 2 | VJ: C, G, T from V's last A | J trimmed by 1.
    CHECK(generated.sequence == "AAAAAAAACGTTTTCC");
    CHECK(generated.realizations.at("vj_dinucl") == std::vector<int>{1, 2, 3});
}

TEST_CASE("Generation: an insertion no Markov chain fills stays as placeholders", "[generation]")
{
    //Observed, intent undecided (plan D5): with no DinucMarkov, the insertion's 'I' markers
    //reach the output as if they were nucleotides.
    ScratchDir dir;
    const Generated generated = generate_one(dir, vdj_parms(false), vdj_marginals(false));
    CHECK(generated.sequence == "TTTTTTTTIIGGTTIIGGGAA");
}

TEST_CASE("Generation: the shipped TRA model with a fixed seed", "[generation]")
{
    const auto [seqs, reals] = generate_shipped("tcr_alpha", 20, 1234);
    CHECK(std::count(seqs.begin(), seqs.end(), '\n') == 21);
    CHECK(fnv1a(seqs) == 0x7253f96e34bc84a1ULL);
    CHECK(fnv1a(reals) == 0x47e291ca18b431c0ULL);
}

TEST_CASE("Generation: the shipped TRB model with a fixed seed", "[generation]")
{
    const auto [seqs, reals] = generate_shipped("tcr_beta", 20, 1234);
    CHECK(std::count(seqs.begin(), seqs.end(), '\n') == 21);
    CHECK(fnv1a(seqs) == 0x6061afd55c208e88ULL);
    CHECK(fnv1a(reals) == 0xee903fa28ca2b8d9ULL);
}
