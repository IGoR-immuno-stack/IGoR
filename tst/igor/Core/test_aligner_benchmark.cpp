/*
 * test_aligner_benchmark.cpp
 *
 *  Catch2 benchmarks for the Aligner, mimicking V/D/J-shaped alignments across
 *  sequencing-technology-realistic read lengths and offset/threshold constraints.
 *  Hidden by default (tag starts with '.'); run explicitly with:
 *      igor_tests "[.benchmark][aligner]"
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

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <igor/Core/Aligner.h>
#include "AlignerTestUtils.h"

#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace igor::test::align;

namespace {

std::string constant_match_sequence(size_t length)
{
    return std::string(length, 'A');
}

std::string constant_mismatch_sequence(size_t length)
{
    return std::string(length, 'G');
}

std::string random_sequence(size_t length, std::mt19937_64 &rng)
{
    static const char alphabet[] = { 'A', 'C', 'G', 'T' };
    std::uniform_int_distribution<int> dist(0, 3);
    std::string s(length, 'A');
    for (auto &c : s) c = alphabet[dist(rng)];
    return s;
}

double realistic_score_threshold(size_t read_length, double match_score)
{
    // TODO(user): tune the fraction; intent is "a real filter", not a precise quality model.
    return 0.5 * static_cast<double>(read_length) * match_score;
}

struct OffsetWindow
{
    int min_offset;
    int max_offset;
    bool rev_offset_frame;
};
const OffsetWindow UNCONSTRAINED{ INT16_MIN, INT16_MAX, false };

struct SequencingScenario
{
    const char *label;
    size_t read_length;
    OffsetWindow v_window; // plain offset, negative = primer inside V, positive = 5' RACE leader
    OffsetWindow j_window; // reversed offset (rev_offset_frame=true), negative = primer inside J,
                           // positive = 3' primer in constant region reading past J
};
const OffsetWindow D_WINDOW{ INT16_MIN, INT16_MAX, false }; // D is always fully spanned at the junction, all scenarios

const std::vector<SequencingScenario> SCENARIOS = {
    // label,      read_len, v_window (min,max,rev)  j_window (min,max,rev)
    { "multiplex", 150, { -300, -50, false }, { -30, -10, true } },  // Variable offset among genes
    { "MiSeq", 300, { -300, -50, false }, { 30, 30, true } },
    { "RACE", 500, { 20, 20, false }, { 30, 30, true } },  // Fixed offset for all genes
};

struct GeneVariant
{
    const char *label;
    Gene_class gene_class;
    size_t ref_length;
};
const std::vector<GeneVariant> GENE_VARIANTS = {
    { "V-like", V_gene, 300 }, // human IGHV ~ 290-300nt
    { "D-like", D_gene, 15 }, // human IGHD ~ 10-37nt
    { "J-like", J_gene, 50 }, // human IGHJ ~ 50-60nt
};

OffsetWindow window_for(const GeneVariant &gene, const SequencingScenario &scenario)
{
    if (gene.gene_class == V_gene) return scenario.v_window;
    if (gene.gene_class == J_gene) return scenario.j_window;
    return D_WINDOW;
}

} // namespace

TEST_CASE("Aligner alignment benchmarks", "[.benchmark][aligner]")
{
    const Matrix<double> matrix = build_test_score_matrix(7, -11);
    const double match_score = 7;
    const int gap_penalty = 13;
    std::mt19937_64 rng(12345); // fixed seed: reproducible benchmark runs

    for (const auto &scenario : SCENARIOS) {
        DYNAMIC_SECTION(scenario.label)
        {
            for (const auto &gene : GENE_VARIANTS) {
                DYNAMIC_SECTION(gene.label)
                {
                    const std::string reference = constant_match_sequence(gene.ref_length);
                    auto aligner = make_legacy_aligner(matrix, gap_penalty, gene.gene_class, { { "g1", reference } });
                    const OffsetWindow constrained = window_for(gene, scenario);
                    const double tight_threshold = realistic_score_threshold(scenario.read_length, match_score);

                    const std::string match_query = constant_match_sequence(scenario.read_length);
                    const std::string mismatch_query = constant_mismatch_sequence(scenario.read_length);
                    // Pre-generate so construction isn't counted inside meter.measure.
                    std::vector<std::string> random_queries;
                    for (int i = 0; i != 20; ++i) {
                        random_queries.push_back(random_sequence(scenario.read_length, rng));
                    }

                    struct ContentCase
                    {
                        const char *label;
                        const std::string *fixed_query; // nullptr => draw from random_queries instead
                    };
                    const std::vector<ContentCase> content_cases = {
                        { "constant match", &match_query },
                        { "constant mismatch", &mismatch_query },
                        { "random", nullptr },
                    };
                    const std::vector<std::pair<const char *, OffsetWindow>> offset_cases = {
                        { "unconstrained offset", UNCONSTRAINED },
                        { "constrained offset", constrained },
                    };
                    const std::vector<std::pair<const char *, double>> threshold_cases = {
                        { "permissive threshold", -1000.0 },
                        { "tight threshold", tight_threshold },
                    };

                    for (const auto &content : content_cases) {
                        for (const auto &offset_case : offset_cases) {
                            for (const auto &threshold_case : threshold_cases) {
                                std::ostringstream name;
                                name << content.label << ", " << offset_case.first << ", " << threshold_case.first;
                                const OffsetWindow &ow = offset_case.second;
                                const double threshold = threshold_case.second;

                                BENCHMARK_ADVANCED(name.str())(Catch::Benchmark::Chronometer meter)
                                {
                                    if (content.fixed_query != nullptr) {
                                        const std::string &query = *content.fixed_query;
                                        meter.measure([&] {
                                            return aligner.align_seq(query, threshold, true, true, ow.min_offset,
                                                                     ow.max_offset, ow.rev_offset_frame);
                                        });
                                    } else {
                                        meter.measure([&](int i) {
                                            return aligner.align_seq(random_queries[i % random_queries.size()],
                                                                     threshold, true, true, ow.min_offset,
                                                                     ow.max_offset, ow.rev_offset_frame);
                                        });
                                    }
                                };
                            }
                        }
                    }
                }
            }
        }
    }
}
