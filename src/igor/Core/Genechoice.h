/*
 * Genechoice.h
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

#include <igor/Core/Rec_Event.h>
#include <igor/Core/JunctionGeometry.h>
#include <igor/Core/Utils.h>
#include <array>
#include <forward_list>
#include <unordered_map>
#include <string>
#include <list>
#include <queue>
#include <utility>
#include <igor/Core/Errorrate.h>
#include <random>
#include <igorCoreExport.h>

/**
 * \class Gene_choice Genechoice.h
 * \brief GeneChoice recombination event.
 * \author Q.Marcou
 * \version 1.0
 *
 * Models the gene choice recombination process.
 * The event realizations are explored based on the sequence alignments that were provdided to the inference.
 * Since D gene can be heavily deleted and might not be recognizable by sequence alignments, a special handling of the D gene choice exploring all D positions ranked by their likelihood has been implemented.
 */
class CORE_EXPORT Gene_choice : public Rec_Event
{
    friend class Coverage_err_counter; //Grant friendship to access current gene realization and offset
    friend class Hypermutation_global_errorrate; //Grant friendship to access current gene realization and offset
    friend class Hypermutation_full_Nmer_errorrate; //Same

public:

    //Constructors
    Gene_choice();
    Gene_choice(Gene_class);
    Gene_choice(Gene_class, std::unordered_map<std::string, Event_realization> &);
    Gene_choice(Gene_class, std::vector<std::pair<std::string, std::string>>);
    //Destructor
    ~Gene_choice() override;
    //Virtual methods overload
    std::shared_ptr<Rec_Event> copy() override;
    // Context-based iterate() interface
    inline void
    iterate(QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation) override;

    void add_realization(int);
    bool add_realization(std::string gene_name, std::string gene_sequence);
    void set_genomic_templates(const std::vector<std::pair<std::string, std::string>> &);
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

    //Proba bound related computation methods
    //Capability queries (task A0)
    OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const override;
    LengthContribution get_length_contribution(SeqTypeId) const override;
    SeqConstructionRole get_seq_construction_role(SeqTypeId) const override;
    OffsetRole get_offset_role(SeqTypeId, Seq_side) const override;

    bool affects_length_of(SegmentSpan) const override;
    int length_delta(const Event_realization &) const override;

private:
    inline double iterate_common(
            double, const int &, int, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            const Marginal_array_p &);

    /// Write the bound for every junction this placement fixes: the one on each flank, and the
    /// neutral element into the junction an internal segment splits. False when no scenario
    /// reaches one of those distances, which is the caller's cue to discard the placement.
    bool write_junction_bounds(Downstream_scenario_proba_bound_map &proba_map, Seq_Offset five_off,
                               Seq_Offset three_off) const;

    /**
     * \brief The error-free length credited to this placement, in nucleotides.
     *
     * **Section 7.1's defect, carried verbatim.** A segment with a movable end at each side
     * credits the core between them -- right in shape, one short of the inclusive count. A
     * segment anchored on a read end credits `size + travel` where at most `size - travel` can
     * survive, so its bound comes out too small and prunes harder than the model justifies.
     * Decision O4: step 3 reproduces both arithmetics rather than deriving one, and R5 replaces
     * them with `core_3 - core_5 + 1`.
     */
    int credited_core_length(Seq_Offset core_5, Seq_Offset core_3) const;

    /**
     * \brief One neighbouring segment end this gene's placement is checked against.
     *
     * Resolved in initialize_event(); iterate() reads it. The two entries a VDJ model produces
     * are the former memory_layer_safety_1 / _2 pair, in the same order -- the neighbour that
     * sits further 5' first -- which is what keeps the claimed layers where they were.
     */
    struct FlankCheck {
        SeqTypeId partner_id = kNoSeqType;
        Seq_side partner_side = Five_prime;   ///< the partner end facing this gene
        bool this_is_left = true;             ///< the partner sits 3' of this gene
        SafetyCell safety_cell;               ///< the pair, by ordering position (S5)
        int safety_layer = -1;                ///< one per *row*, so same-row checks share it
        int partner_offset_layer = -1;
        bool partner_exists = false;          ///< the model has a gene choice for that segment
        bool partner_chosen = false;          ///< ...and it is placed before this one
        bool active = false;                  ///< per scenario: the partner is placed, so check
    };

    /// Score one exhaustive placement against the read and write its error bound. Shared by
    /// the two sub-branches of the position scan; see the definition.
    void score_placement_against_read(const QuerySequenceContext &query,
                                      AccumulationContext &accumulation,
                                      ExplorationContext &exploration,
                                      JunctionGeometry::OffsetInterval five_reach,
                                      JunctionGeometry::OffsetInterval three_reach);

    /// Where this segment's own ends can still travel from a given placement, under §7.4's
    /// short bound. One line of R10 away from being `pending_.reachable()`.
    JunctionGeometry::OffsetInterval own_five_prime_reach(Seq_Offset five_off) const
    {
        return {static_cast<Seq_Offset>(five_off + own_five_prime_travel_.min),
                static_cast<Seq_Offset>(five_off + own_five_prime_travel_.max)};
    }
    JunctionGeometry::OffsetInterval own_three_prime_reach(Seq_Offset three_off) const
    {
        return {static_cast<Seq_Offset>(three_off + own_three_prime_travel_.min),
                static_cast<Seq_Offset>(three_off + own_three_prime_travel_.max)};
    }

    /// Section 7.1's two arithmetics. See credited_core_length().
    enum class EndogenousCore { Inflated, Truncated };

    std::vector<FlankCheck> flank_checks_;

    /// What each checked neighbour's facing end can still reach, and where it sits. Indexed by
    /// the neighbour's SeqTypeId; filled once per scenario by the preamble of iterate().
    std::vector<JunctionGeometry::OffsetInterval> neighbour_reach_;
    std::vector<Seq_Offset> neighbour_offset_;

    /// Every pending offset modifier in the model, rebuilt at initialize_event(). The four
    /// hand-written deletion lookups it replaces are the same query asked four times.
    JunctionGeometry::PendingModifierBounds pending_;

    /// Resolved from the ordering, not from event_class: a gene at an end of the constructed
    /// sequence is anchored by the read rather than by a neighbour, which is what makes its
    /// template able to overhang, its error-free core credited the other way, and its position
    /// not something to scan for when the aligner finds nothing.
    bool clip_template_before_read_ = false;
    bool clip_template_after_read_ = false;
    bool exhaustive_position_fallback_ = false;
    bool publishes_alignment_state_ = true;
    EndogenousCore endogenous_core_ = EndogenousCore::Inflated;

    /// This placement's own two ends, for the scenario being explored.
    Seq_Offset my_5_off = 0;
    Seq_Offset my_3_off = 0;

    //Inference variables
    //Bool checks

    //Offsets checks


    /// True until an alignment survives every check, which is what makes the exhaustive
    /// position scan run at all. Was `no_d_align`; the regression script and the plan keep
    /// that name for the *path*, which is the only D left in it.
    bool no_alignment_survived;

    /// The placement the exhaustive scan is currently scoring: where it sits, how long it is,
    /// and where it disagrees with the read.
    Seq_Offset placement_5_off;
    Seq_Offset placement_3_off;
    std::size_t template_size;
    std::vector<std::size_t> placement_mismatches;

    //Declare common variables
    mutable int base_index;
    double new_scenario_proba;
    double proba_contribution;
    Int_Str gene_seq;
    int new_index;
    const int *alignment_offset_p;
    /// Scratch for the position scan only; the alignment path counts its core inline.
    std::vector<std::size_t>::const_iterator placement_mism_iter;
    std::size_t endogeneous_mismatches;

    //Constants
    //Memory Layers
    int memory_layer_cs;
    int memory_layer_mismatches;
    int memory_layer_off_threep;
    int memory_layer_off_fivep;
    int memory_layer_proba_map_seq;

    /// How far this segment's own ends can still travel, for the exhaustive position scan.
    /// §7.4's short answer rather than `pending_`'s, carried until R10 -- see
    /// JunctionGeometry::legacy_offset_delta().
    OffsetDelta own_five_prime_travel_{};
    OffsetDelta own_three_prime_travel_{};

    /// This event's realizations by their index, so the retained decomposition can carry an
    /// index rather than a gene name -- §2.5's per-candidate hash lookup, removed. Points into
    /// event_realizations, which nothing mutates after initialize_event().
    std::vector<const Event_realization *> realizations_by_index_;

    /// The nearest *candidate* neighbour on each side -- placed or not -- as indices into
    /// flank_checks_, or -1. Distinct from the junction's endpoints, which are the nearest
    /// *placed* ones: the sliding window is bounded by a neighbour that has not been placed
    /// just as much as by one that has, through the read end the preamble records for it.
    int nearest_left_check_ = -1;
    int nearest_right_check_ = -1;
};
