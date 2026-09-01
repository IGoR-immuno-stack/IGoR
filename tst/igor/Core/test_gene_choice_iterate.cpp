/*
 * test_gene_choice_iterate.cpp
 *
 *  Characterization tests for Gene_choice::iterate().
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

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * These are **characterization** tests: every expected value was read off the current
 * implementation, not derived from what it ought to compute. Where a pinned value is known
 * to be wrong, the section name says so and points at the plan.
 *
 * The file is organised by the generic pattern each branch belongs to rather than by
 * V/D/J, because V/D/J is the axis B11 removes -- see
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md sections 6.1 and 4.
 */

#include "test_utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace IgorTestUtils;

namespace {

/// Dump a recorder's observations. Used while deriving expected values; kept because the
/// next person to add a section will want it too.
void dump(const std::shared_ptr<RecordingEvent> &rec, const char *label)
{
    UNSCOPED_INFO(label << ": " << rec->call_count() << " call(s)");
    for (std::size_t i = 0; i != rec->calls.size(); ++i) {
        const ScenarioSnapshot &s = rec->calls[i];
        UNSCOPED_INFO("  [" << i << "] proba=" << s.scenario_proba);
        for (const auto &[seq_type, off] : s.offsets) {
            UNSCOPED_INFO("      offsets[" << seq_type << "] = (" << off.first << ", " << off.second
                                           << ")");
        }
        for (const auto &[seq_type, seq] : s.sequences) {
            UNSCOPED_INFO("      seq[" << seq_type << "] = \"" << seq << "\"");
        }
        for (const auto &[safety, value] : s.safety) {
            UNSCOPED_INFO("      safety[" << safety << "] = " << value);
        }
        for (const auto &[seq_type, bound] : s.downstream_bounds) {
            UNSCOPED_INFO("      bound[" << seq_type << "] = " << bound);
        }
    }
}

} // namespace

TEST_CASE("Gene_choice::iterate baseline writes (G4)", "[gene_choice][iterate][baseline]")
{
    SECTION("V: one alignment writes offsets, sequence and mismatches")
    {
        const std::string gene = "ACGTACGTACGT";
        const std::string read = gene + "TTTTTT";

        auto state = create_iterate_state(read);
        auto v_event = make_gene_choice(V_gene, {{"V1", gene}}, 0, /*fixed=*/false);
        state.set_alignments(V_gene, {create_perfect_alignment("V1", 0, gene.size())});
        state.set_marginal(0, 1.0L);

        auto rec = call_iterate_recording(v_event, state);
        dump(rec, "V baseline");

        REQUIRE(rec->call_count() == 1);
        const ScenarioSnapshot &s = rec->calls.at(0);
        CHECK(s.five_prime(V_gene_seq) == 0);
        CHECK(s.three_prime(V_gene_seq) == static_cast<Seq_Offset>(gene.size()) - 1);
        CHECK(s.sequences.at(V_gene_seq) == gene);
        CHECK(s.mismatches.at(V_gene_seq).empty());
        CHECK(s.scenario_proba == 1.0);
    }
}
