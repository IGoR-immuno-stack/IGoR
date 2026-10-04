/*
 * Dinuclmarkov.cpp
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
 */

#include <igor/Model/Legacy/Dinuclmarkov.h>
#include <igor/Alignment/Legacy/Aligner.h>
#include <igor/Model/Legacy/EventUtils.h>
#include <igor/Model/Legacy/JsonDetail.h>

#include <algorithm>
#include <vector>

#include <cassert>


namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

using namespace std;

namespace {
/**
 * Resolve a Dinucl_markov's seq_type name to the legacy insertion enum.
 *
 * The Len_proba machinery is still Seq_type-keyed, so the string has to come back to the enum
 * somewhere; this is the one place, and the one place that rejects an unrecognised name. It
 * replaces the `correct_class` flag that used to do the same job inline in
 * Dinucl_markov::iterate_initialize_Len_proba().
 */
} // namespace

Dinucl_markov::Dinucl_markov(Seq_type seq_type) : Rec_Event(), total_nucl_count(0), ins_seq_type(seq_type)
{
    this->type = Event_type::Dinuclmarkov_t;
    //Same indexes as Aligner::nt2int

    event_realizations.emplace("A", Event_realization("A", INT16_MAX, "A", Int_Str(), 0));
    event_realizations.emplace("C", Event_realization("C", INT16_MAX, "C", Int_Str(), 1));
    event_realizations.emplace("G", Event_realization("G", INT16_MAX, "G", Int_Str(), 2));
    event_realizations.emplace("T", Event_realization("T", INT16_MAX, "T", Int_Str(), 3));

    dinuc_proba_matrix = Matrix<double>(kIntNtCount, kIntNtCount);
    this->update_event_name();
}

Dinucl_markov::Dinucl_markov(const nlohmann::json &node)
    : Dinucl_markov(str2SeqType(igor::model::legacy::json_detail::require(node, "seq_type").get<Seq_type_String>()))
{
    using namespace igor::model::legacy::json_detail;
    reject_unknown_keys(node, kEventKeys);
    expect_type(node, "DinucMarkov");

    this->set_seq_type(require(node, "seq_type").get<Seq_type_String>());
    this->set_priority(require(node, "priority").get<int>());
    this->set_nickname(require(node, "nickname").get<string>());

    // The call order used to matter here, and it no longer does: name_side() in Rec_Event.cpp
    // keeps a DinucMarkov's generated name on Undefined_side whatever its side, so the name
    // this constructor produces is the text reader's. The side itself is load-bearing beyond
    // the name — it is the traversal anchor, so Three_prime means the chain is seeded from the
    // segment on its left.
    this->set_event_side(str2SeqSide(require(node, "side").get<string>()));

    // The four nucleotides are self-initialized above, so the document's realizations are
    // verified rather than applied. Anything else means the document was not written by this
    // schema, and building a chain over a different alphabet would need more than a loop here.
    const auto realizations = realizations_in_index_order(node);
    if (realizations.size() != this->event_realizations.size())
        throw std::runtime_error("event json: a dinucleotide Markov chain has "
                                 + std::to_string(this->event_realizations.size())
                                 + " realizations, document has "
                                 + std::to_string(realizations.size()));
    for (const nlohmann::json *realization : realizations) {
        const auto name = require(*realization, "name").get<string>();
        const auto found = this->event_realizations.find(name);
        if (found == this->event_realizations.end()
            || found->second.index != require(*realization, "index").get<int>())
            throw std::runtime_error("event json: unexpected dinucleotide realization \"" + name
                                     + "\" at index "
                                     + std::to_string(require(*realization, "index").get<int>()));
    }
}

Dinucl_markov::~Dinucl_markov()
{
    // TODO delete realization indices
}

/**
 * Which junction this event fills, and which segment seeds it.
 *
 * Replaces a switch over the three legacy junctions *and* the resolution step that first
 * replaced it. Both facts are already on the event: Model_Parms::finalize() resolved the
 * neighbours, and `event_side` says which of the two the Markov chain runs from -- the
 * ordering is the model's fact, the direction is this event's.
 *
 * Deliberately the *ordering* neighbour and not DynamicSequenceMap's occupancy-skipping walk:
 * the two differ exactly when the anchor is empty, which today is plan section 7.12's crash.
 * Swapping them is that defect's fix and a behaviour change, so it lands separately (7.11).
 */
DinuclTraversalSpec Dinucl_markov::get_junction() const
{
    DinuclTraversalSpec spec;
    if (this->seq_type_id == kNoSeqType) {
        return spec;
    }
    spec.target_id = this->seq_type_id;
    spec.anchor_side = this->event_side;
    if (spec.anchor_side == Three_prime) {
        spec.anchor_id = this->left_adjacent_id;
    } else if (spec.anchor_side == Five_prime) {
        spec.anchor_id = this->right_adjacent_id;
    } else {
        //No direction declared: the model does not say which way this chain runs, and the two
        //neighbours are equally adjacent.
        return spec;
    }

    //The generation path is still keyed by the Seq_type enum. Resolve the handles when the
    //ids allow it and flag them otherwise, rather than leaving a plausible-looking default.
    if (spec.anchor_id != kNoSeqType
        && static_cast<std::size_t>(spec.target_id) < kLegacySeqTypeCount
        && static_cast<std::size_t>(spec.anchor_id) < kLegacySeqTypeCount) {
        spec.target_seq = static_cast<Seq_type>(spec.target_id);
        spec.anchor_seq = static_cast<Seq_type>(spec.anchor_id);
        spec.legacy_enums_valid = true;
    }
    return spec;
}

shared_ptr<Rec_Event> Dinucl_markov::copy()
{
    //TODO rewrite this by invoking a copy constructor?
    shared_ptr<Dinucl_markov> new_dinucl_markov_p = shared_ptr<Dinucl_markov>(new Dinucl_markov(this->ins_seq_type));
    new_dinucl_markov_p->priority = this->priority;
    new_dinucl_markov_p->nickname = this->nickname;
    new_dinucl_markov_p->fixed = this->fixed;
    new_dinucl_markov_p->update_event_name();
    new_dinucl_markov_p->set_event_identifier(this->event_index);
    new_dinucl_markov_p->set_seq_type(this->get_seq_type());
    new_dinucl_markov_p->set_seq_type_id(this->get_seq_type_id());
    new_dinucl_markov_p->set_event_side(this->get_side());
    //Adjacency is not carried here: Model_Parms' copy constructor re-runs finalize(), which
    //is the one place that resolves it.
    return new_dinucl_markov_p;
}

int Dinucl_markov::size() const
{
    return event_realizations.size() * event_realizations.size();
}

std::vector<std::size_t> Dinucl_markov::inherent_shape() const
{
    return { event_realizations.size(), event_realizations.size() };
}

/**
 * @brief Context-based iterate() implementation
 *
 * Unpacks 5 context objects into legacy parameters and delegates
 * to the existing iterate() implementation.
 *

 */
void Dinucl_markov::iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    base_index = exploration.index_map.get(this->event_index);
    proba_contribution = 1;

    //Clear all previous scenario realizations
    current_realizations_index_vec.clear();

    //No emptiness check on traversal_specs: initialize_event() refuses to leave them empty, so
    //the arm this replaces was unreachable through the production path -- the same dead
    //backstop B6 removed from Insertion::iterate.
    {
        const DinuclTraversalSpec spec = get_junction();
        previous_seq = (*scenario.constructed_sequences.get(spec.anchor_id));

        //The junction is *created* here, not filled here (O12 decision (a'), plan section
        //7.13). Its length comes from the offsets the Insertion placed, which is a property of
        //the scenario state -- "offsets placed, sequence not yet created" -- and not a lookup
        //from this event to that one. Which event established those offsets is exactly what
        //A0 exists to stop anyone needing to know.
        //
        //What this buys beyond the layer: no partially-constructed segment exists at any
        //hand-off, so `int_undefined` never leaves this event and the leaf invariant's
        //content half becomes global rather than conditional on somebody declaring `Fills`.
        //
        //The offsets are Insertion's own, so a negative width would mean it handed on a
        //realization its own realization map does not hold. It discards those instead
        //(iterate_common returns 0), hence the assert rather than a branch.
        const int junction_length = scenario.seq_offsets.get(spec.target_id, Three_prime)
                                    - scenario.seq_offsets.get(spec.target_id, Five_prime) + 1;
        assert(junction_length >= 0 && "junction offsets cross over");
        junction_str.assign(static_cast<std::size_t>(std::max(junction_length, 0)), int_undefined);
        Int_Str &target_seq = junction_str;
        scenario.constructed_sequences.set(spec.target_id, &junction_str, memory_layer_seq);

        //One entry per position filled, appended as the fill proceeds, where a fixed-width
        //array would keep the previous scenario's value at any position this one skips. There
        //are none to skip -- iterate_common() writes every position -- but the append is what
        //makes that a fact about the state rather than a coincidence the reader has to check
        //for.
        realization_indices.clear();
        //Capacity comes from the paired Insertion's longest realization, so the push_backs
        //below never reallocate. Worth stating: a junction that outgrew it would still be
        //correct, just quietly allocating in the hot loop.
        assert(realization_indices.capacity() >= target_seq.size()
               && "junction longer than its Insertion's longest realization");

        //The chain is seeded from the anchor's nucleotide facing the junction, and Int_Str is a
        //std::vector<int>, so front()/back() on a fully deleted anchor dereferences nullptr --
        //a segfault, not a wrong answer (plan section 7.12). The state is reachable on the
        //current corpus rather than only under tandem D: Deletion's V 3' branch has no size
        //guard at all, and its D 5' branch guards with a strict >, so deleting exactly the
        //whole segment is a legal realization and nothing downstream rejects the result.
        //
        //Throw, do not discard (decided Sep 8 2026). Such a scenario is geometrically
        //legitimate, so dropping it silently would remove probability mass the model should
        //account for and leave no trace that it happened; a scenario that cannot be scored is
        //a modelling error, and the user has to see it.
        //
        //Only where a seed is actually needed. An empty junction chooses no nucleotide, so it
        //never reads the anchor and stays perfectly scoreable -- which is what section 7.12
        //asks for: an anchor carrying no nucleotide must be rejected or skipped, *never read*.
        //
        //This is the throw half of section 7.12 alone. Its other half, the occupancy-skipping
        //walk that would seed from the next non-empty segment instead, is deferred to B10: per
        //section 7.11 it needs absence to have a defined meaning, and it needs the junction
        //offsets R3 writes, since skipping a gene segment lands the walk on a junction.
        //
        //anchor_side is the anchor's end facing the junction, so it also says which way the
        //chain runs: from a 3' anchor the read window follows it, from a 5' anchor it precedes
        //it and both the window and the filled segment are handled back to front.
        const bool reverse_traversal = (spec.anchor_side == Five_prime);

        if (target_seq.empty()) {
            //Nothing to choose, so nothing to seed from, so the anchor is not read at all --
            //which is what keeps an emptily-anchored scenario scoreable when its junction is
            //also empty. Bitwise for every case that worked before: iterate_common() was
            //already a no-op on an empty junction, and data_seq_substr is read nowhere else.
        } else if (previous_seq.empty()) {
            const SeqTypeRegistry &registry = scenario.constructed_sequences.registry();
            throw runtime_error("Dinucl_markov " + this->name + ": the anchor segment "
                                + registry.name(spec.anchor_id) + " carries no nucleotide to seed "
                                "the Markov chain from. This scenario deleted it entirely -- a "
                                "legal realization, but one that leaves " + registry.name(spec.target_id)
                                + " unscoreable.");
        } else if (reverse_traversal) {
            const size_t char_index =
                    scenario.seq_offsets.get(spec.anchor_id, spec.anchor_side) - target_seq.size();
            data_seq_substr = query.int_sequence.substr(char_index, target_seq.size());
            previous_nt_str = previous_seq.front();
            reverse(data_seq_substr.begin(), data_seq_substr.end());
            iterate_common(realization_indices, previous_nt_str, target_seq,
                           model.model_parameters);
            reverse(target_seq.begin(), target_seq.end());
        } else {
            const size_t start_index =
                    scenario.seq_offsets.get(spec.anchor_id, spec.anchor_side) + 1;
            data_seq_substr = query.int_sequence.substr(start_index, target_seq.size());
            previous_nt_str = previous_seq.back();
            iterate_common(realization_indices, previous_nt_str, target_seq,
                           model.model_parameters);
        }

        exploration.downstream_proba_map.set(spec.target_id, 1.0, memory_layer_junction);
    }

    scenario.scenario_proba *= proba_contribution;

    //Compute scenario downstream proba bound
    scenario_upper_bound_proba = exploration.compute_upper_bound(
        scenario.scenario_proba,
        current_downstream_proba_memory_layers
    );

    if (!exploration.should_prune(scenario_upper_bound_proba)) {
        Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
    }
}

namespace {

void require_legacy_enums(const DinuclTraversalSpec &spec, const string &name)
{
    if (!spec.legacy_enums_valid) {
        throw invalid_argument("Dinucl_markov " + name
                               + ": generation needs a junction and an anchor that the Seq_type "
                                 "enum names. Model_Parms::finalize() must have run, and the "
                                 "topology must be one of the legacy ones (see B9).");
    }
}

} // namespace

vector<int> Dinucl_markov::draw_realization(const Marginal_array_p &, const unordered_map<Rec_Event_name, int> &,
                                            const unordered_map<Seq_type, string> &constructed_sequences,
                                            mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);

    const DinuclTraversalSpec spec = get_junction();
    require_legacy_enums(spec, this->name);

    //Copies: the draw reads what has been built and writes nothing. The chain runs away from
    //its anchor, so a chain seeded from the anchor's 5' end reads that anchor backwards.
    string target_ins_seq = constructed_sequences.at(spec.target_seq);
    string anchor_seq = constructed_sequences.at(spec.anchor_seq);
    if (spec.anchor_side == Five_prime) {
        reverse(anchor_seq.begin(), anchor_seq.end());
    }
    return this->draw_chain(anchor_seq, target_ins_seq, distribution, generator);
}

/**
 * The chain, drawn into a copy of the insertion: one entry per placeholder, in chain order,
 * kNoRealization where the walk chose nothing. A position that is not a placeholder is not
 * drawn, and it seeds the next one like a drawn nucleotide does.
 */
vector<int> Dinucl_markov::draw_chain(const string &previous_seq, string &inserted_seq,
                                      uniform_real_distribution<double> &distribution, mt19937_64 &generator) const
{

    vector<int> chain;
    double prob_count;
    if (!inserted_seq.empty()) {
        double rand;
        if (inserted_seq[0] == 'I') {
            rand = distribution(generator);
            prob_count = 0;
            int drawn = kNoRealization;
            int prev_nt = nt2int(previous_seq.substr(previous_seq.size() - 1, 1)).at(0);
            for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
                 iter != event_realizations.end(); ++iter) {
                prob_count += this->dinuc_proba_matrix(prev_nt, (*iter).second.index);
                if (prob_count >= rand) {
                    inserted_seq[0] = (*iter).second.value_str[0];
                    drawn = (*iter).second.index;
                    break;
                }
            }
            chain.push_back(drawn);
        }
        for (size_t i = 1; i != inserted_seq.size(); ++i) {
            if (inserted_seq[i] == 'I') {
                int prev_nt = nt2int(inserted_seq.substr(i - 1, 1)).at(0);

                rand = distribution(generator);
                prob_count = 0;
                int drawn = kNoRealization;

                for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
                     iter != event_realizations.end(); ++iter) {
                    prob_count += this->dinuc_proba_matrix(prev_nt, (*iter).second.index);
                    if (prob_count >= rand) {
                        inserted_seq[i] = (*iter).second.value_str[0];
                        drawn = (*iter).second.index;
                        break;
                    }
                }
                chain.push_back(drawn);
            }
        }
    }
    return chain;
}

void Dinucl_markov::construct_realization(const vector<int> &chain,
                                          unordered_map<Seq_type, string> &constructed_sequences) const
{
    const DinuclTraversalSpec spec = get_junction();
    require_legacy_enums(spec, this->name);

    //One chain entry per placeholder, in chain order; the chain then reads 5'->3' once turned
    //around for a chain that ran from the right.
    string &target_ins_seq = constructed_sequences.at(spec.target_seq);
    size_t next = 0;
    for (char &position : target_ins_seq) {
        if (position != 'I') {
            continue;
        }
        if (next < chain.size() and chain[next] != kNoRealization) {
            position = this->realization_at(chain[next]).value_str[0];
        }
        ++next;
    }
    if (spec.anchor_side == Five_prime) {
        reverse(target_ins_seq.begin(), target_ins_seq.end());
    }
}

void Dinucl_markov::write2txt(ofstream &outfile)
{
    write2txt_legacy(outfile);
}

void Dinucl_markov::write2txt_legacy(ofstream &outfile)
{
    // Derive legacy gene_class from seq_type for backward compatibility
    // Use to_string() to match the exact strings expected by str2GeneClass()
    string legacy_gene_class;
    if (seq_type == "VD_ins_seq") legacy_gene_class = to_string(VD_genes);
    else if (seq_type == "DJ_ins_seq") legacy_gene_class = to_string(DJ_genes);
    else if (seq_type == "VJ_ins_seq") legacy_gene_class = to_string(VJ_genes);
    else legacy_gene_class = to_string(Undefined_gene);
    outfile << "#DinucMarkov;" << legacy_gene_class << ";" << Undefined_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_str << ";" << (*iter).second.index << endl;
    }
}

void Dinucl_markov::write2txt_v2(ofstream &outfile)
{
    outfile << "#DinucMarkov;Undefined_gene;" << seq_type << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_str << ";" << (*iter).second.index << endl;
    }
}


void Dinucl_markov::iterate_common(std::vector<int> &realization_indices, int &previous_assigned_nt,
                                   Int_Str &ins_seq, const Marginal_array_p &model_parameters_point)
{

    //Every position is written. The junction is created by iterate() for this scenario alone
    //(O12 (a')), so there is no earlier fill to preserve; the guard that skipped positions
    //already holding a nucleotide belonged to the buffer once shared with the Insertion, and
    //went with it (R11, plan section 7.13).
    if (!ins_seq.empty()) {

        //first_nt_index = event_realizations.at(previous_assigned_nt).index;
        //sec_nt_index = event_realizations.at(data_seq_substr.substr(0,1)).index;

        first_nt_index = previous_assigned_nt; //[0] -'0';
        sec_nt_index = data_seq_substr[0]; //-'0';

        current_realizations_index_vec.emplace_back(sec_nt_index);

        //For this Dinucl_Markov model the values on the marginal array represents the conditional probability of a couple of nucleotides (N2 | N1)
        if ((first_nt_index < 4) & (sec_nt_index < 4)) {
            offset = first_nt_index * event_realizations.size();
            realization_final_index = base_index + offset + sec_nt_index;
            proba_contribution *= model_parameters_point
                    [realization_final_index]; ///compute_nt_freq(base_index+offset , model_parameters_point);
            realization_indices.push_back(realization_final_index);
        } else {
            //If an ambiguous nucleotide is present we take the average probability over possible underlying nts
            proba_contribution *= dinuc_proba_matrix(first_nt_index, sec_nt_index);
            realization_indices.push_back(-1);
        }

        ins_seq.at(0) = data_seq_substr.at(0);
        total_nucl_count += 1;

        for (size_t i = 1; i != ins_seq.size(); ++i) {

            //first_nt_index = event_realizations.at(data_seq_substr.substr(i-1,1)).index;
            //sec_nt_index = event_realizations.at(data_seq_substr.substr(i,1)).index;

            first_nt_index = data_seq_substr[i - 1]; // -'0';
            sec_nt_index = data_seq_substr[i]; // -'0';

            current_realizations_index_vec.emplace_back(sec_nt_index);

            //For this Dinucl_Markov model the values on the marginal array represents the joint probability of a couple of nucleotides (N1 , N2)
            if ((first_nt_index < 4) & (sec_nt_index < 4)) {
                offset = first_nt_index * event_realizations.size();
                realization_final_index = base_index + offset + sec_nt_index;
                proba_contribution *= model_parameters_point
                        [base_index + offset
                         + sec_nt_index]; ///compute_nt_freq(base_index+offset , model_parameters_point);
                realization_indices.push_back(realization_final_index);
            } else {
                //If an ambiguous nucleotide is present we take the average probability over possible underlying nts
                proba_contribution *= dinuc_proba_matrix(first_nt_index, sec_nt_index);
                realization_indices.push_back(-1);
            }

            ins_seq.at(i) = data_seq_substr.at(i);
            total_nucl_count += 1;
        }
    }
}

/*
 * This way of proceeding is highly not optimal, should be modified
 */
double Dinucl_markov::compute_nt_freq(int index, const Marginal_array_p &model_marginals) const
{
    double nucl_freq = 0;
    for (size_t i = 0; i != event_realizations.size(); ++i) {
        nucl_freq += model_marginals[index + i];
    }
    return nucl_freq;
}

void Dinucl_markov::ind_normalize(Marginal_array_p &marginal_array_p, size_t base_index) const
{
    size_t numb_realizations = this->event_realizations.size();
    for (size_t i = 0; i != numb_realizations; ++i) {
        long double sum_marginals = 0;
        for (size_t j = 0; j != numb_realizations; ++j) {
            sum_marginals += marginal_array_p[base_index + i * numb_realizations + j];
        }
        if (sum_marginals != 0) {
            for (size_t j = 0; j != numb_realizations; ++j) {
                marginal_array_p[base_index + i * numb_realizations + j] /= sum_marginals;
            }
        }
    }
}

void Dinucl_markov::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const Events_map &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        SafetyMatrix &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{

    const DinuclTraversalSpec spec = get_junction();
    if (spec.anchor_id == kNoSeqType) {
        throw invalid_argument("Dinucl_markov " + this->name
                               + ": no junction to fill. Its seq_type must be in the registry "
                                 "(Model_Parms::finalize() resolves that), it must declare which "
                                 "side its chain runs from, and that side must have a neighbour.");
    }

    downstream_proba_map.request_layer(spec.target_id);
    memory_layer_junction = downstream_proba_map.claimed_layer(spec.target_id);

    //The segment itself, now that this event creates it (O12 (a')). Nothing else writes an
    //insertion's sequence, so the claimed layer is 0 and the write lands exactly where the
    //Insertion's set_current() used to put it -- the ownership is what changes, not the
    //storage.
    constructed_sequences.request_layer(spec.target_id);
    memory_layer_seq = constructed_sequences.claimed_layer(spec.target_id);

    //One index slot per position the junction can ever hold, taken from the Insertion that
    //allocates it, so the per-scenario push_backs never reallocate.
    const int longest = EventUtils::get_insertion_len_max(
            constructed_sequences.registry().name(spec.target_id), events_map);
    realization_indices.reserve(static_cast<std::size_t>(std::max(longest, 0)));
    longest_junction_ = std::max(longest, 0);

    index_map.set_current_layer(this->event_index, 0);
    unmutable_base_index = index_map.get(this->event_index);

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list, seq_offsets,
                                      index_map);
}

/**
 * \bug Will only count realizations of unambiguous nucleotides (realization indices>=0 since they are set to -1 in iterate_common)
 */
void Dinucl_markov::add_to_marginals(long double scenario_proba, Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        for (size_t i = 0; i != this->event_marginal_size; ++i) {
            updated_marginals[unmutable_base_index + i] = 0;
        }
    }

    for (const int index : realization_indices) {
        if (index >= 0) {
            updated_marginals[index] += scenario_proba;
        }
    }
}

/**
 * Update probability values contained in a matrix where coordinates >=4 indicates ambiguous nucleotides
 * We simply take the average probability over the different possible nucleotides.
 */
void Dinucl_markov::update_event_internal_probas(const Marginal_array_p &marginal_array,
                                                 const unordered_map<Rec_Event_name, int> &index_map)
{
    Int_nt const all_nt_vals[] = { int_A, int_C, int_G, int_T, int_R, int_Y, int_K, int_M,
                                   int_S, int_W, int_B, int_D, int_H, int_V, int_N };
    size_t event_index = index_map.at(this->get_name());
    for (size_t i = 0; i != kIntNtCount; ++i) {
        for (size_t j = 0; j != kIntNtCount; ++j) {
            list<Int_nt> previous_list = get_ambiguous_nt_list(all_nt_vals[i]);
            list<Int_nt> next_list = get_ambiguous_nt_list(all_nt_vals[j]);

            //Reset the dinuc proba matrix
            this->dinuc_proba_matrix(i, j) = 0;

            for (Int_nt prev_nt : previous_list) {
                for (Int_nt next_nt : next_list) {
                    this->dinuc_proba_matrix(i, j) +=
                            marginal_array[event_index + prev_nt * event_realizations.size() + next_nt];
                }
            }
            //By taking the average we assume all nucleotides underlying the ambiguous one are  equally probable
            this->dinuc_proba_matrix(i, j) /= (double)(previous_list.size() * next_list.size());
        }
    }
}

OffsetDelta Dinucl_markov::get_offset_delta_bounds(SeqTypeId, Seq_side) const
{
    return {};
}

LengthContribution Dinucl_markov::get_length_contribution(SeqTypeId) const
{
    //Nothing, even though this event now creates the segment: its realizations are single
    //nucleotides, and how many of them there are was decided by the offsets its Insertion
    //placed. Creating a segment and contributing length to a span are separate statements,
    //and this is the event that separates them. It is also the subclass whose len_min /
    //len_max were never set at all, and which therefore still carries their INT16 sentinels.
    return {};
}

SeqConstructionRole Dinucl_markov::get_seq_construction_role(SeqTypeId type_id) const
{
    //Creates, not Fills (O12 (a')). The segment it writes is one it allocated itself, at a
    //layer it claimed, sized from offsets that were already placed -- so a sibling scenario
    //can no longer find its nucleotides where it expects placeholders (plan section 7.13).
    return type_id == this->seq_type_id ? SeqConstructionRole::Creates : SeqConstructionRole::None;
}

OffsetRole Dinucl_markov::get_offset_role(SeqTypeId, Seq_side) const
{
    return OffsetRole::None;
}

bool Dinucl_markov::affects_proba_of(SegmentSpan span, const SeqTypeRegistry &registry) const
{
    //True exactly where the segment this model fills lies inside the span: the p^L factor scales
    //with that segment's length, so it belongs to every span containing it. The same ordering
    //question affects_length_of() asks, on the segment this event creates the sequence of.
    //
    //It used to validate the seq_type name as well, throwing on anything but the three legacy
    //insertions. That check belonged to the enum it resolved to; a model's seq_types are
    //validated where they are resolved, and initialize_event() refuses a Dinucl_markov with no
    //junction to fill.
    return lies_strictly_inside(registry, this->seq_type_id, span);
}

int Dinucl_markov::length_delta(const Event_realization &) const
{
    //Never reached: affects_length_of() is always false, so the fold takes the probability
    //branch and never enumerates this event. Stated rather than defaulted, so an event that
    //starts contributing length has to say so here too.
    return 0;
}

double Dinucl_markov::span_proba_factor(SegmentSpan, const UnfilledSegmentLengths &lengths) const
{
    //p^L over the segment this model creates, whose length was published by whoever placed its
    //offsets -- the Insertion, which under O12 (a') creates the offsets and not the sequence.
    //Absent means nobody placed them on this path, and the contribution is 1.
    const SeqTypeId filled = this->seq_type_id;
    if (not lengths.has(filled)) {
        return 1.0;
    }
    const int length = lengths.length_of(filled);
    if (length >= 0 and static_cast<std::size_t>(length) < chain_bound_.size()) {
        return chain_bound_[static_cast<std::size_t>(length)];
    }
    //Longer than any junction this event's Insertion can place, or no sweep has prepared the
    //chain: the cruder bound, every nucleotide at the best pair's probability, still holds.
    return pow(this->get_upper_bound_proba(), length);
}

/*
 * The best probability the chain can give L nucleotides it has not seen, for every L up to the
 * longest junction: a Viterbi pass over the fifteen codes a read position can hold -- the four
 * bases and the eleven ambiguity codes -- with each pair priced as iterate_common() prices it: the
 * marginal for two bases, the average over the bases they stand for otherwise
 * (update_event_internal_probas()). The seed is any code, since the anchor's nucleotide is not
 * known when the table is built. It replaces p^L, every nucleotide at the best pair's probability,
 * which it can only lower: every pair it multiplies is at most that best pair.
 *
 * Over the four bases alone it would not be a bound. An ambiguous position is priced by an
 * average, and a chain of averages can beat every chain of bases: with P(A|G) = 1 and P(T|C) = 1,
 * a junction reading M (A or C) then T after a G anchor is priced 0.5 x 0.625, where no chain of
 * two bases reaches more than 0.25 once the other rows are uniform.
 */
void Dinucl_markov::prepare_span_proba_factor(const Marginal_array_p &model_parameters, const Index_map &base_index_map)
{
    Int_nt const all_nt_vals[] = { int_A, int_C, int_G, int_T, int_R, int_Y, int_K, int_M,
                                   int_S, int_W, int_B, int_D, int_H, int_V, int_N };
    const int base = base_index_map.get(this->event_index, 0);
    const std::size_t nucleotides = event_realizations.size();
    //One block of nucleotides x nucleotides per configuration of the parents -- one, in every
    //model shipped -- and a pair has to be bounded whichever block the scenario reads.
    const std::size_t blocks = std::max<std::size_t>(1, this->event_marginal_size / this->size());

    double pair[kIntNtCount][kIntNtCount];
    for (std::size_t i = 0; i != kIntNtCount; ++i) {
        const list<Int_nt> previous_list = get_ambiguous_nt_list(all_nt_vals[i]);
        for (std::size_t j = 0; j != kIntNtCount; ++j) {
            const list<Int_nt> next_list = get_ambiguous_nt_list(all_nt_vals[j]);
            double best = 0.0;
            for (std::size_t block = 0; block != blocks; ++block) {
                double sum = 0.0;
                for (Int_nt prev_nt : previous_list) {
                    for (Int_nt next_nt : next_list) {
                        sum += model_parameters[base + block * this->size() + prev_nt * nucleotides + next_nt];
                    }
                }
                best = std::max(best, sum / (double)(previous_list.size() * next_list.size()));
            }
            pair[i][j] = best;
        }
    }

    chain_bound_.assign(1, 1.0);
    std::vector<double> ending(kIntNtCount, 1.0); //best chain of the current length ending in each code
    for (int length = 1; length <= longest_junction_; ++length) {
        std::vector<double> next(kIntNtCount, 0.0);
        for (std::size_t j = 0; j != kIntNtCount; ++j) {
            for (std::size_t i = 0; i != kIntNtCount; ++i) {
                next[j] = std::max(next[j], ending[i] * pair[i][j]);
            }
        }
        ending = std::move(next);
        chain_bound_.push_back(*std::max_element(ending.begin(), ending.end()));
    }
}

void Dinucl_markov::update_event_name()
{
    Seq_type_String seq_type_str;
    switch (ins_seq_type) {
    case VD_ins_seq: seq_type_str = "VD_genes"; break;
    case DJ_ins_seq: seq_type_str = "DJ_gene"; break;
    case VJ_ins_seq: seq_type_str = "VJ_gene"; break;
    default: seq_type_str = to_string(this->event_class); break;
    }
    // name_side() is what keeps this name on Undefined_side: the junction is already named by
    // seq_type_str just above, so the side would add a direction, not an identity. This override
    // exists for that class token, and it has to apply the base class's side rule.
    this->name = string() + this->type + "_" + seq_type_str + "_"
                 + to_string(name_side(this->type, this->event_side))
                 + "_prio" + to_string(priority) + "_size" + to_string(this->size());
}

} // namespace igor::model::legacy
