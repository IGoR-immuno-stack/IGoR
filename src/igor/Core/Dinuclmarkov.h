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

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
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
 * One junction a Dinucl_markov event fills, and where it takes its seed nucleotide from.
 *
 * `anchor_side` is the *anchor's* end facing the junction: `Three_prime` means the anchor sits
 * to the 5' side and the chain runs forward, `Five_prime` that it sits to the 3' side and the
 * chain runs backwards. It comes from the event's own `event_side`, which Model_Parms sets
 * from the model file -- it is a property of the Markov chain's direction, not of the
 * topology, so the registry supplies *which* segment is the neighbour and the event supplies
 * *which side* it seeds from.
 *
 * The per-junction scratch state lives here too, so that a model with more than one junction
 * per event needs no more named members.
 */
struct DinuclTraversalSpec {
    SeqTypeId target_id = kNoSeqType;
    SeqTypeId anchor_id = kNoSeqType;
    Seq_side anchor_side = Undefined_side;

    ///@{ \name Legacy enum handles, for the generation path only
    /// draw_random_realization() writes into an `unordered_map<Seq_type, string>`, which no
    /// non-legacy seq_type can key. Resolved when the names allow it and unused otherwise;
    /// B9 removes them with the rest of the generation path's enum dependency.
    Seq_type target_seq = VD_ins_seq;
    Seq_type anchor_seq = V_gene_seq;
    bool legacy_enums_valid = false;
    ///@}

    /// Marginal index per filled position, or -1 where the pair was ambiguous. Sized at
    /// initialize_event() from the paired Insertion's longest realization.
    std::vector<int> realization_indices;
    /// How much of `realization_indices` the last iterate() filled.
    std::size_t filled_size = 0;
    /// Layer this event claimed in the downstream-proba map for `target_id`.
    int memory_layer = -1;
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
class CORE_EXPORT Dinucl_markov : public Rec_Event
{
public:
        using Rec_Event::initialize_event;
        using Rec_Event::initialize_crude_scenario_proba_bound;

    //Constructors
    Dinucl_markov(Seq_type); //TODO should be scalable on one side easily (mono di tri quadri nucl)
    //Destructor
    ~Dinucl_markov() override;

    //Accessors
    std::shared_ptr<Rec_Event> copy() override;
    void resolve_topology(const SeqTypeRegistry &registry) override;

    /// The junctions this event fills, as resolved by resolve_topology(). Exposed for tests:
    /// a topology with no Seq_type enum entry cannot be checked through iterate() until the
    /// harness can build events by SeqTypeId.
    const std::vector<DinuclTraversalSpec> &get_traversal_specs() const { return traversal_specs; }
    int size() const override;

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
    void ind_normalize(Marginal_array_p &, size_t) const override;
    
    void update_event_name() override;
    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;
    void add_to_marginals(long double, Marginal_array_p &) const override;
    void update_event_internal_probas(const Marginal_array_p &, const std::unordered_map<Rec_Event_name, int> &) override;

    double *get_updated_ptr() override;
    void initialize_crude_scenario_proba_bound(
            double &, std::forward_list<double *> &,
            const Events_map &) override;

    //Proba bound related computation methods
    //Capability queries (task A0)
    OffsetDelta get_offset_delta_bounds(SeqTypeId, Seq_side) const override;
    LengthContribution get_length_contribution(SeqTypeId) const override;
    SeqConstructionRole get_seq_construction_role(SeqTypeId) const override;
    OffsetRole get_offset_role(SeqTypeId, Seq_side) const override;

    bool has_effect_on(Seq_type) const override;
    void iterate_initialize_Len_proba(Seq_type considered_junction, std::map<int, double> &length_best_proba_map,
                                      std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
                                      const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                      Seq_type_str_p_map &constructed_sequences, int &seq_len) const override;
    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                    const Marginal_array_p &model_parameters_point, Index_map &base_index_map) override;

private:
    double *updated_upper_bound_proba =
            nullptr; //This points to a double modified by the Insertion event given the number of insertion
    Matrix<double> dinuc_proba_matrix;

    int total_nucl_count;

    Int_Str previous_seq; //&
    size_t previous_seq_size;
    int previous_nt_str;
    Int_Str data_seq_substr;

    mutable int base_index;
    int unmutable_base_index;
    double new_scenario_proba;
    double proba_contribution;
    mutable bool correct_class;

    Seq_type ins_seq_type;
    std::vector<DinuclTraversalSpec> traversal_specs;

    //std::pair<Seq_type,Seq_side> v_5_pair = std::make_pair (V_gene_seq,Five_prime);
    //std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);

    //Iterate common
    int first_nt_index;
    int sec_nt_index;
    int offset;
    int realization_final_index;

    inline void iterate_common(int *, int &, Int_Str &, const Marginal_array_p &);
    inline std::queue<int> draw_random_common(const std::string &, std::string &, const Marginal_array_p &, int,
                                              std::uniform_real_distribution<double> &, std::mt19937_64 &) const;
    inline double compute_nt_freq(int, const Marginal_array_p &) const;
};
