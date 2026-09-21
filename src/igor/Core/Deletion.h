/*
 * Deletion.h
 *
 *  Created on: Dec 9, 2014
 *      Author: Quentin Marcou
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

#pragma once

#include <igor/Core/JunctionGeometry.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <igor/Core/Errorrate.h>

#include <forward_list>
#include <unordered_map>
#include <string>
#include <list>
#include <queue>
#include <utility>
#include <random>
#include <math.h>

#include <igorCoreExport.h>

/**
 * \class Deletion Deletion.h
 * \brief Deletion recombination event
 * \author Q.Marcou
 * \version 1.0
 *
 * The Deletion RecEvent models deletions of one genomic fragment on a given side.
 * Deletions can be either positive or negative (= Palindromic insertions)
 *
 * By construction the corresponding GeneChoice must have been explored first.
 */
class CORE_EXPORT Deletion : public Rec_Event
{
    friend class Coverage_err_counter; //Grant friendship to access the current number of deletion
    friend class Hypermutation_global_errorrate; //Grant friendship to access the current number of deletion
    friend class Hypermutation_full_Nmer_errorrate; //Same
    friend class DeletionTest; // For unit testing private members

public:
        //Constructor
        Deletion();
        Deletion(Seq_type, Seq_side, std::pair<int, int>);
        Deletion(Seq_type, Seq_side);
        Deletion(Seq_type, Seq_side, std::unordered_map<std::string, Event_realization> &);
        ~Deletion() override;

    //Virtual methods
    std::shared_ptr<Rec_Event> copy() override;

    // Context-based iterate() interface
    inline void
    iterate(QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation) override;

    void add_realization(int);
    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const override;
    void write2txt(std::ofstream &) override;
    void write2txt_legacy(std::ofstream &) override;
    void write2txt_v2(std::ofstream &) override;
    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, SafetyMatrix &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;
    void add_to_marginals(long double, Marginal_array_p &) const override;

    //Capability queries (task A0)
    OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const override;
    LengthContribution get_length_contribution(SeqTypeId) const override;
    SeqConstructionRole get_seq_construction_role(SeqTypeId) const override;
    OffsetRole get_offset_role(SeqTypeId, Seq_side) const override;

    //Proba bound related computation methods
    bool affects_length_of(SegmentSpan) const override;
    int length_delta(const Event_realization &) const override;

private:
    /// Record \a span as the junction this deletion widens, in the slot matching the side it
    /// trims. Called from initialize_event() once the neighbour is known.
    void resolve_junction(SegmentSpan span, SeqTypeId proba_key, int memory_layer);

    /// Keep the mismatches that survive the trim. See the definition for the direction.
    void trim_mismatches(const std::vector<std::size_t> &previous_mismatches);

    Seq_type target_seq_type;

    inline void iterate_common(
            std::forward_list<Event_realization>::const_iterator &, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            const Marginal_array_p &);

    std::forward_list<Event_realization> int_value_and_index;

    /**
     * \brief One neighbouring gene segment this deletion's moving end is checked against.
     *
     * Resolved in initialize_event(), ordered 5'->3' by the partner's position -- the order the
     * four arms this replaces performed their comparisons in. A deletion is checked only
     * against the partners on the side it trims: the other side of its segment is not moving,
     * so nothing there can newly collide.
     */
    struct FlankCheck {
        SeqTypeId partner_id = kNoSeqType;
        Seq_side partner_side = Five_prime;   ///< the partner end facing this deletion
        SafetyCell safety_cell;               ///< the pair, by ordering position (S5)
        int safety_layer = -1;                ///< one per *row*, so same-row checks share it
        int partner_offset_layer = -1;
        bool partner_chosen = false;          ///< the partner is placed before this event

        /// How far the partner's facing end can still travel. `JunctionGeometry::
        /// legacy_offset_delta()`, which is §7.4's short answer rather than
        /// `PendingModifierBounds`'; R10 is the swap.
        OffsetDelta partner_delta{};

        //Per scenario, filled by the preamble of iterate():
        Seq_Offset offset = 0;                            ///< where the partner's facing end is
        JunctionGeometry::OffsetInterval reach{};         ///< ...and where it can still go
        bool active = false;                  ///< placed, and not already established safe
    };

    std::vector<FlankCheck> flank_checks_;


    /// Index into flank_checks_ of the neighbour bounding the junction this deletion widens --
    /// the nearest placed one on the trimmed side -- or -1 when there is none.
    int junction_partner_ = -1;

    /// `event_side == Three_prime`, cached: it selects the direction of every asymmetry in the
    /// body, from which end of the template is cut to which way the mismatch list is trimmed.
    bool trims_three_prime_ = false;

    /// Resolved from the ordering, not from the gene class. A segment anchored on an end of the
    /// read keeps at least one nucleotide of its template (§2.7), carries the dominated early
    /// prune stage (§6.14), and does *not* bounds-check its palindrome's read positions
    /// (§7.15); an internal one does the opposite in all three, and an internal 5' deletion
    /// additionally carries §7.3's surviving `//FIXME`.
    bool keep_one_nucleotide_ = false;
    bool early_prune_stage_ = false;
    bool guard_palindrome_positions_ = false;
    bool discard_offset_outside_read_ = false;

    /// Another unprocessed event still moves the opposite end of my segment, so nothing of the
    /// template is settled and its error bound stays neutral. Was `d_del_opposite_side_processed`,
    /// negated, and answered by an A0 query rather than by naming the D deletion on the far side.
    bool opposite_end_still_moves_ = false;

    /// This deletion's own end, for the scenario being explored.
    Seq_Offset my_offset = 0;
    Seq_Offset my_new_offset = 0;

    //Common variables
    mutable int base_index;
    double new_scenario_proba;
    double proba_contribution;
    int new_index;
    mutable Int_Str new_str;
    mutable Int_Str tmp_str;
    mutable std::string gen_new_str;
    mutable std::string gen_tmp_str;
    std::vector<std::size_t> mismatches_vector;
    std::vector<std::size_t>::const_iterator mis_iter;
    bool end_reached;
    int deletion_value;

    int memory_layer_cs;
    int memory_layer_mismatches;
    int memory_layer_offset_del;
    int memory_layer_proba_map_seq;
};

std::string &make_transversions(std::string &, bool);
Int_Str &make_transversions(Int_Str &);

bool del_numb_compare(const Event_realization &, const Event_realization &);
