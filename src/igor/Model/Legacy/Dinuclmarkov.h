/*
 * Dinuclmarkov.h
 *
 *  Created on: Mar 4, 2015
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
 *
 *      Please note that although this class inherits Rec_event properties, it is not currently designed to be parent
 *      of any other event in the model_parms graph
 */

#pragma once

#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Core/Legacy/Utils.h>
#include <igor/Math/Legacy/Matrix.h>
#include <forward_list>
#include <unordered_map>
#include <string>
#include <list>
#include <queue>
#include <utility>
#include <igor/Model/Legacy/Errorrate.h>
#include <random>

#include <igor/Model/Export.h>

/**
 * One junction a Dinucl_markov event fills, and where it takes its seed nucleotide from.
 *
 * `anchor_side` is the *anchor's* end facing the junction: `Three_prime` means the anchor sits
 * to the 5' side and the chain runs forward, `Five_prime` that it sits to the 3' side and the
 * chain runs backwards. It comes from the event's own `event_side`, which Model_Parms sets
 * from the model file -- it is a property of the Markov chain's direction, not of the
 * topology, so the registry supplies *which* segment is the neighbour and the event supplies
 * *which side* it seeds from.
 */

namespace igor::model::legacy {
using namespace igor::math::legacy;
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

struct DinuclTraversalSpec {
    SeqTypeId target_id = kNoSeqType;
    SeqTypeId anchor_id = kNoSeqType;
    Seq_side anchor_side = Undefined_side;
};

/**
 * \class Dinucl_markov Dinucl_markov.h
 * \brief Dinucleotide insertion Markov model.
 * \author Q.Marcou
 * \version 1.0
 *
 * Models a Markov chain dictating the identity of inserted nucleotides in the inserted region.
 * We assume a low error frequency and almost flat dinucleotide model regime such that we use an euristic to extract the most likely realization.
 * This choice has been made because the full handling through a forward algorithm would not be able to cope with e.g context dependent errors.
 *
 * By construction the Insertion event must have been explored first
 */
class MODEL_EXPORT Dinucl_markov : public Rec_Event
{
public:
        using Rec_Event::initialize_event;

    //Constructors
    Dinucl_markov(Seq_type); //TODO should be scalable on one side easily (mono di tri quadri nucl)
    /// From one event node of the JSON model schema. See ModelJson.h.
    explicit Dinucl_markov(const nlohmann::json &);
    //Destructor
    ~Dinucl_markov() override;

    //Accessors
    std::shared_ptr<Rec_Event> copy() override;
    /**
     * The junction this event fills and where it seeds from, derived rather than stored.
     *
     * Everything it needs is already on the event: its own seq_type id, the neighbours
     * Model_Parms::finalize() resolved, and `event_side` -- which says which way the Markov
     * chain runs and is therefore the event's own fact, not the topology's. Nothing to
     * resolve, nothing to keep in sync, nothing to make idempotent.
     *
     * `anchor_id` is `kNoSeqType` when the model does not define one: no ordering, no
     * declared direction, or nothing on the side the chain runs from. initialize_event()
     * refuses that; iterate() may then assume it.
     */
    DinuclTraversalSpec get_junction() const;
    int size() const override;
    /// {from, to}: a transition matrix, not a flat list. See Rec_Event::inherent_shape().
    std::vector<std::size_t> inherent_shape() const override;

    // Context-based iterate() interface
    inline void
    iterate(QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation) override;

    void add_realization(int);
    std::vector<int> draw_realization(const Marginal_array_p &, const std::unordered_map<Rec_Event_name, int> &,
                                      const GenerationState &, std::mt19937_64 &) const override;
    void construct_realization(const std::vector<int> &, GenerationState &) const override;
    /// A chain is sequence-valued: it has no single index to condition a child's row on, and it
    /// never moved one.
    void propagate_realization(
            const std::vector<int> &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &)
            const override
    {}
    void write2txt(std::ofstream &) override;
    void write2txt_legacy(std::ofstream &) override;
    void write2txt_v2(std::ofstream &) override;
    void ind_normalize(Marginal_array_p &, size_t) const override;
    
    void update_event_name() override;
    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, SafetyMatrix &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;
    void add_to_marginals(long double, Marginal_array_p &) const override;
    void update_event_internal_probas(const Marginal_array_p &, const std::unordered_map<Rec_Event_name, int> &) override;

    //Proba bound related computation methods
    //Capability queries (task A0)
    OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const override;
    LengthContribution get_length_contribution(SeqTypeId) const override;
    SeqConstructionRole get_seq_construction_role(SeqTypeId) const override;
    OffsetRole get_offset_role(SeqTypeId, Seq_side) const override;

    bool affects_proba_of(SegmentSpan, const SeqTypeRegistry &) const override;
    int length_delta(const Event_realization &) const override;
    double span_proba_factor(SegmentSpan, const UnfilledSegmentLengths &) const override;
    void prepare_span_proba_factor(const Marginal_array_p &, const Index_map &) override;

private:
    Matrix<double> dinuc_proba_matrix;

    int total_nucl_count;

    Int_Str previous_seq; //&
    /// The junction this event creates and fills, sized from the offsets its Insertion placed
    /// (O12 (a'), R1). Owned here rather than borrowed through the pointer the Insertion used
    /// to leave in the scenario's sequence map -- plan section 7.13.
    Int_Str junction_str;
    size_t previous_seq_size;
    int previous_nt_str;
    Int_Str data_seq_substr;

    mutable int base_index;
    int unmutable_base_index;
    double new_scenario_proba;
    double proba_contribution;

    Seq_type ins_seq_type;

    /// One entry per position this event actually filled in the last scenario: the marginal
    /// index credited there, or -1 where the pair involved an ambiguous nucleotide. Cleared
    /// and refilled per scenario; its capacity is reserved once at initialize_event() from the
    /// paired Insertion's longest realization, so the hot loop never allocates.
    std::vector<int> realization_indices;
    /// The longest junction this event can be asked to fill: its Insertion's longest realization.
    int longest_junction_ = 0;
    /// Per junction length L, the best probability the chain gives L nucleotides it has not seen
    /// (prepare_span_proba_factor()). Index 0 holds 1.
    std::vector<double> chain_bound_;
    /// Layer claimed in the downstream-proba map for the junction this event fills.
    int memory_layer_junction = -1;
    /// Layer claimed in the constructed-sequence map for that same junction, which this event
    /// now creates rather than filling in place.
    int memory_layer_seq = -1;

    //std::pair<Seq_type,Seq_side> v_5_pair = std::make_pair (V_gene_seq,Five_prime);
    //std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);

    //Iterate common
    int first_nt_index;
    int sec_nt_index;
    int offset;
    int realization_final_index;

    inline void iterate_common(std::vector<int> &, int &, Int_Str &, const Marginal_array_p &);
    std::vector<int> draw_chain(const std::string &, std::string &, std::uniform_real_distribution<double> &,
                                std::mt19937_64 &) const;
    inline double compute_nt_freq(int, const Marginal_array_p &) const;
};

} // namespace igor::model::legacy
