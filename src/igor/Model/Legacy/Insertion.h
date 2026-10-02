/*
 * Insertion.h
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
 *
 */

#pragma once

#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Core/Legacy/Utils.h>
#include <forward_list>
#include <unordered_map>
#include <string>
#include <list>
#include <queue>
#include <utility>
#include <igor/Model/Legacy/Errorrate.h>
#include <random>
#include <map>

#include <igor/Model/Export.h>

/**
 * \class Insertion Insertion.h
 * \brief Insertion recombination events.
 * \author Q.Marcou
 * \version 1.0
 *
 *  The Insertion RecEvent models the distribution of junctional insertion length during the V(D)J recombination process.
 *
 */

namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

class MODEL_EXPORT Insertion : public Rec_Event
{
public:
    //Constructors
    Insertion();
    Insertion(Seq_type, std::pair<int, int>);
    Insertion(Seq_type, std::forward_list<int>);
    Insertion(Seq_type);
    Insertion(Seq_type, std::unordered_map<std::string, Event_realization> &);
    /// From one event node of the JSON model schema. See ModelJson.h.
    explicit Insertion(const nlohmann::json &);

    //Destructor
    ~Insertion() override;

    //virtual methods
    std::shared_ptr<Rec_Event> copy() override;
    // Context-based iterate() interface
    inline void
    iterate(QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation) override;

    bool add_realization(int);
    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const override;
    void write2txt(std::ofstream &) override;
    void write2txt_legacy(std::ofstream &) override;
    void write2txt_v2(std::ofstream &) override;

    void update_event_name() override;

    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, SafetyMatrix &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;
    void add_to_marginals(long double, Marginal_array_p &) const override;

    /// Throws unless `events_map` holds the Dinucl_markov that fills this event's junction.
    /// Called by initialize_event(); public so the check can be tested without the rest of it.
    void require_dinucl_markov(const Events_map &events_map) const;

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
            double, int, int, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            const Marginal_array_p &);

    std::map<int, Event_realization> ordered_realization_map;

    mutable int base_index;
    double new_scenario_proba;
    double proba_contribution;
    int insertions;
    int new_index;
    int previous_index;

    //Iterate common
    int realization_index;
    std::string insertions_str;

    /// Layers claimed for the two ends of the junction this event places (O12 (a'), R3).
    int memory_layer_offset_fivep = -1;
    int memory_layer_offset_threep = -1;

    Seq_type ins_seq_type;

    //Pre create pairs to call seq_offsets (otherwise cost of creating a pair at each call)
    //std::pair<Seq_type,Seq_side> d_5_pair = std::make_pair (D_gene_seq,Five_prime);
    //std::pair<Seq_type,Seq_side> v_3_pair = std::make_pair (V_gene_seq,Three_prime);
    //std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);
    //std::pair<Seq_type,Seq_side> d_3_pair = std::make_pair (D_gene_seq,Three_prime);
};

} // namespace igor::model::legacy
