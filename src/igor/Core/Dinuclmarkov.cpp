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

#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/EventUtils.h>

#include <vector>

#include <cassert>

using namespace std;

Dinucl_markov::Dinucl_markov(Seq_type seq_type) : Rec_Event(), total_nucl_count(0), ins_seq_type(seq_type)
{
    this->type = Event_type::Dinuclmarkov_t;
    //Same indexes as Aligner::nt2int

    event_realizations.emplace("A", Event_realization("A", INT16_MAX, "A", Int_Str(), 0));
    event_realizations.emplace("C", Event_realization("C", INT16_MAX, "C", Int_Str(), 1));
    event_realizations.emplace("G", Event_realization("G", INT16_MAX, "G", Int_Str(), 2));
    event_realizations.emplace("T", Event_realization("T", INT16_MAX, "T", Int_Str(), 3));

    updated = true;
    updated_upper_bound_proba = new double;

    dinuc_proba_matrix = Matrix<double>(kIntNtCount, kIntNtCount);
    this->update_event_name();
}

Dinucl_markov::~Dinucl_markov()
{
    // TODO delete realization indices
    if (updated_upper_bound_proba)
        delete updated_upper_bound_proba;
}

/**
 * Turn "which junction do I fill, and who seeds it" into stored ids.
 *
 * Replaces a switch over the three legacy junctions. The ordering says which segment is
 * adjacent; `event_side` says which side this event's Markov chain runs from, which is a
 * property of the model and not of the topology -- Model_Parms derives it from the gene class
 * for legacy files and reads it directly for v2 ones.
 *
 * Deliberately the *ordering* neighbour and not DynamicSequenceMap's occupancy-skipping walk:
 * the two differ exactly when the anchor is empty, which today is plan section 7.12's crash.
 * Swapping them is that defect's fix and a behaviour change, so it lands separately (7.11).
 *
 * Idempotent: initialize_event() calls it again for models built in code rather than read from
 * a file, which never reach Model_Parms::finalize().
 */
void Dinucl_markov::resolve_topology(const SeqTypeRegistry &registry)
{
    traversal_specs.clear();
    if (this->seq_type_id == kNoSeqType
        || static_cast<std::size_t>(this->seq_type_id) >= registry.total_count()) {
        return;
    }

    DinuclTraversalSpec spec;
    spec.target_id = this->seq_type_id;
    spec.anchor_side = this->event_side;
    if (spec.anchor_side == Three_prime) {
        spec.anchor_id = registry.left_neighbor(spec.target_id);
    } else if (spec.anchor_side == Five_prime) {
        spec.anchor_id = registry.right_neighbor(spec.target_id);
    } else {
        //No direction declared: the model does not say which way this chain runs.
        return;
    }
    if (spec.anchor_id == kNoSeqType) {
        return;
    }

    //The generation path is still keyed by the Seq_type enum; resolve the handles when the
    //names allow it, and leave the spec usable for inference either way.
    const Seq_type_String &target_name = registry.name(spec.target_id);
    const Seq_type_String &anchor_name = registry.name(spec.anchor_id);
    try {
        spec.target_seq = str2SeqType(target_name);
        spec.anchor_seq = str2SeqType(anchor_name);
        spec.legacy_enums_valid = true;
    } catch (const std::runtime_error &) {
        spec.legacy_enums_valid = false;
    }

    traversal_specs.push_back(std::move(spec));
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
    //Carry the resolved topology: a per-thread copy is made after Model_Parms::finalize() has
    //run, and the generation path never initializes the copy.
    new_dinucl_markov_p->traversal_specs = this->traversal_specs;
    return new_dinucl_markov_p;
}

int Dinucl_markov::size() const
{
    return event_realizations.size() * event_realizations.size();
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
    for (auto &spec : this->traversal_specs) {
        previous_seq = (*scenario.constructed_sequences.get(spec.anchor_id));
        Int_Str &target_seq = *const_cast<Int_Str *>(scenario.constructed_sequences.get(spec.target_id));
        //One entry per position actually filled, appended as the fill proceeds: a position
        //that already held a nucleotide contributes no probability term and now records no
        //index either, where the fixed-width array kept the previous scenario's value there.
        spec.realization_indices.clear();
        //Capacity comes from the paired Insertion's longest realization, so the push_backs
        //below never reallocate. Worth stating: a junction that outgrew it would still be
        //correct, just quietly allocating in the hot loop.
        assert(spec.realization_indices.capacity() >= target_seq.size()
               && "junction longer than its Insertion's longest realization");

        //anchor_side is the anchor's end facing the junction, so it also says which way the
        //chain runs: from a 3' anchor the read window follows it, from a 5' anchor it precedes
        //it and both the window and the filled segment are handled back to front.
        const bool reverse_traversal = (spec.anchor_side == Five_prime);

        if (reverse_traversal) {
            const size_t char_index =
                    scenario.seq_offsets.get(spec.anchor_id, spec.anchor_side) - target_seq.size();
            data_seq_substr = query.int_sequence.substr(char_index, target_seq.size());
            previous_nt_str = previous_seq.front();
            reverse(data_seq_substr.begin(), data_seq_substr.end());
            iterate_common(spec.realization_indices, previous_nt_str, target_seq,
                           model.model_parameters);
            reverse(target_seq.begin(), target_seq.end());
        } else {
            const size_t start_index =
                    scenario.seq_offsets.get(spec.anchor_id, spec.anchor_side) + 1;
            data_seq_substr = query.int_sequence.substr(start_index, target_seq.size());
            previous_nt_str = previous_seq.back();
            iterate_common(spec.realization_indices, previous_nt_str, target_seq,
                           model.model_parameters);
        }

        exploration.downstream_proba_map.set(spec.target_id, 1.0, spec.memory_layer);
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

queue<int> Dinucl_markov::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        unordered_map<Seq_type, string> &constructed_sequences, mt19937_64 &generator) const
{

    uniform_real_distribution<double> distribution(0.0, 1.0);
    int index = index_map.at(this->get_name());
    queue<int> realization_queue;

    if (this->traversal_specs.empty()) {
        throw invalid_argument(std::string("Unknown seq_type for DinuclMarkov model: ") + to_string(this->ins_seq_type));
    }

    for (const auto &spec : this->traversal_specs) {
        string &target_ins_seq = constructed_sequences.at(spec.target_seq);
        string anchor_seq = constructed_sequences.at(spec.anchor_seq);
        bool reverse_traversal = (spec.anchor_side == Five_prime);
        if (reverse_traversal) {
            reverse(anchor_seq.begin(), anchor_seq.end());
        }

        queue<int> tmp = this->draw_random_common(anchor_seq, target_ins_seq, model_marginals_p, index, distribution,
                                                  generator);
        while (!tmp.empty()) {
            realization_queue.push(tmp.front());
            tmp.pop();
        }

        if (reverse_traversal) {
            reverse(target_ins_seq.begin(), target_ins_seq.end());
        }
    }

    return realization_queue;
}

queue<int> Dinucl_markov::draw_random_common(const string &previous_seq, string &inserted_seq,
                                             const Marginal_array_p &model_marginals_p, int index,
                                             uniform_real_distribution<double> &distribution,
                                             mt19937_64 &generator) const
{

    queue<int> realization_queue;
    double prob_count;
    if (!inserted_seq.empty()) {
        double rand;
        if (inserted_seq[0] == 'I') {
            rand = distribution(generator);
            prob_count = 0;
            /*THIS WAS REMOVED WHEN INTRODUCING AMBIGUOUS NUCLEOTIDES SUPPORT
			 * int offset;
			try{
				offset = event_realizations.at(previous_seq.substr(previous_seq.size()-1,1)).index*event_realizations.size();
			}
			catch(exception& except){
				cout<<"exception caught in DinucMarkov draw random common, key used: "<<previous_seq.substr(previous_seq.size()-1,1);
				throw except;
			}
			*/
            int prev_nt = nt2int(previous_seq.substr(previous_seq.size() - 1, 1)).at(0);
            for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
                 iter != event_realizations.end(); ++iter) {
                //prob_count += model_marginals_p[index + offset + (*iter).second.index];
                prob_count += this->dinuc_proba_matrix(prev_nt, (*iter).second.index);
                if (prob_count >= rand) {
                    inserted_seq[0] = (*iter).second.value_str[0];
                    realization_queue.push((*iter).second.index);
                    break;
                }
            }
        }
        for (size_t i = 1; i != inserted_seq.size(); ++i) {
            if (inserted_seq[i] == 'I') {
                /*THIS WAS REMOVED WHEN INTRODUCING AMBIGUOUS NUCLEOTIDES SUPPORT
				int offset;
				try{
					offset = event_realizations.at(inserted_seq.substr(i-1,1)).index*event_realizations.size();
				}
				catch(exception& except){
					cout<<"exception caught, key used: "<<inserted_seq.substr(i-1,1)<<",ins seq: "<<inserted_seq<<",i = "<<i<<", previous rand: "<<rand<<", previous prob_count: "<<prob_count<<endl;
					throw except;
				}
				*/
                int prev_nt = nt2int(inserted_seq.substr(i - 1, 1)).at(0);

                rand = distribution(generator);
                prob_count = 0;

                for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
                     iter != event_realizations.end(); ++iter) {
                    //prob_count += model_marginals_p[index + offset + (*iter).second.index];
                    prob_count += this->dinuc_proba_matrix(prev_nt, (*iter).second.index);
                    if (prob_count >= rand) {
                        inserted_seq[i] = (*iter).second.value_str[0];
                        realization_queue.push((*iter).second.index);
                        break;
                    }
                }
            }
        }
    }
    return realization_queue;
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

    if (!ins_seq.empty()) {
        if (ins_seq.at(0) == int_undefined) {

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
        }

        for (size_t i = 1; i != ins_seq.size(); ++i) {
            if (ins_seq.at(i) == int_undefined) {

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
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{

    //A model built in code never reaches Model_Parms::finalize(), so resolve here too. The
    //override is idempotent; for a model read from a file this recomputes the same specs.
    this->resolve_topology(constructed_sequences.registry());

    if (this->traversal_specs.empty()) {
        throw invalid_argument("Dinucl_markov " + this->name
                               + ": no junction to fill. Its seq_type must be in the registry, it "
                                 "must declare which side its chain runs from, and that side must "
                                 "have a neighbouring segment.");
    }

    for (auto &spec : this->traversal_specs) {
        downstream_proba_map.request_layer(spec.target_id);
        spec.memory_layer = downstream_proba_map.claimed_layer(spec.target_id);

        //One index slot per position the junction can ever hold, taken from the Insertion that
        //allocates it. Sized per spec, so a model with several junctions per event needs no
        //new members -- and freed with the event, unlike the three raw arrays this replaces.
        const int longest = EventUtils::get_insertion_len_max(
                constructed_sequences.registry().name(spec.target_id), events_map);
        spec.realization_indices.reserve(static_cast<std::size_t>(std::max(longest, 0)));
    }

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

    for (const auto &spec : this->traversal_specs) {
        for (const int index : spec.realization_indices) {
            if (index >= 0) {
                updated_marginals[index] += scenario_proba;
            }
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

double *Dinucl_markov::get_updated_ptr()
{
    return updated_upper_bound_proba;
}

void Dinucl_markov::initialize_crude_scenario_proba_bound(
        double &downstream_proba_bound, forward_list<double *> &updated_proba_list,
        const Events_map &events_map)
{
    this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
    this->updated_proba_bounds_list = updated_proba_list;
    updated_proba_list.push_front(this->updated_upper_bound_proba);
}

OffsetDelta Dinucl_markov::get_offset_delta_bounds(SeqTypeId, Seq_side) const
{
    return {};
}

LengthContribution Dinucl_markov::get_length_contribution(SeqTypeId) const
{
    //Fills placeholders that Insertion already allocated, so it adds no nucleotides of its
    //own. This is the subclass whose len_min / len_max were never set at all, and which
    //therefore still carries their INT16 sentinels.
    return {};
}

SeqConstructionRole Dinucl_markov::get_seq_construction_role(SeqTypeId type_id) const
{
    return type_id == this->seq_type_id ? SeqConstructionRole::Fills : SeqConstructionRole::None;
}

OffsetRole Dinucl_markov::get_offset_role(SeqTypeId, Seq_side) const
{
    return OffsetRole::None;
}

bool Dinucl_markov::has_effect_on(Seq_type seq_type_param) const
{
    const string &heo_dm_st = this->seq_type;
    if (heo_dm_st == "VD_ins_seq") {
        return (seq_type_param == VJ_ins_seq || seq_type_param == VD_ins_seq);
    } else if (heo_dm_st == "VJ_ins_seq") {
        return (seq_type_param == VJ_ins_seq);
    } else if (heo_dm_st == "DJ_ins_seq") {
        return (seq_type_param == VJ_ins_seq || seq_type_param == DJ_ins_seq);
    }
    return false;
}

void Dinucl_markov::iterate_initialize_Len_proba(Seq_type considered_junction,
                                                 std::map<int, double> &length_best_proba_map,
                                                 std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                                 double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                                 Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                                 int &seq_len /*=0*/) const
{
    base_index_map.set_current_layer(this->event_index, 0);
    base_index = base_index_map.get(this->event_index);

    correct_class = 0;
    const string &iilp_dm_st = this->seq_type;
    if (iilp_dm_st == "VD_ins_seq") {
        correct_class = 1;
        if (this->has_effect_on(considered_junction)) {
            if (constructed_sequences.exists(VD_ins_seq)) {
                scenario_proba *= pow(this->get_upper_bound_proba(), constructed_sequences.get(VD_ins_seq)->size());
            }
            //Otherwise the proba contribution is 1
        }
    }
    if (iilp_dm_st == "DJ_ins_seq") {
        correct_class = 1;
        if (this->has_effect_on(considered_junction)) {
            if (constructed_sequences.exists(DJ_ins_seq)) {
                scenario_proba *= pow(this->get_upper_bound_proba(), constructed_sequences.get(DJ_ins_seq)->size());
            }
            //Otherwise the proba contribution is 1
        }
    }
    if (iilp_dm_st == "VJ_ins_seq") {
        correct_class = 1;
        if (this->has_effect_on(considered_junction)) {
            if (constructed_sequences.exists(VJ_ins_seq)) {
                scenario_proba *= pow(this->get_upper_bound_proba(), constructed_sequences.get(VJ_ins_seq)->size());
            }
            //Otherwise the proba contribution is 1
        }
    }
    if (!correct_class) {
        throw invalid_argument(std::string("Unknown seq_type for DinuclMarkov model: ") + iilp_dm_st);
    }

    //TODO use a better proba bound for this dinucleotide markov model

    //Recursive call
    Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map, model_queue,
                                                    scenario_proba, model_parameters_point, base_index_map,
                                                    constructed_sequences, seq_len /*=0*/);
}

void Dinucl_markov::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                               const Marginal_array_p &model_parameters_point,
                                               Index_map &base_index_map)
{
    //Do nothing
    /*
	 * For now let's assume nothing can happen to the junction once the dinucleotide has been chosen
	 * =>no errors
	 * =>no in/dels
	 */
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
    this->name = string() + this->type + "_" + seq_type_str + "_" + to_string(this->event_side)
                 + "_prio" + to_string(priority) + "_size" + to_string(this->size());
}
