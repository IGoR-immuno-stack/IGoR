/*
 * Insertion.cpp
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

#include <igor/Core/Insertion.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/gene_to_seqtype_migr.h>

#include <algorithm>

#include <algorithm>

#include <limits>

using namespace std;

namespace {
/// Convert the seq_type string (e.g. "VD_ins_seq") to a Seq_type enum value.
/// Returns false when the string is not a recognised insertion seq_type.
///
/// Callers left are the ones feeding a genuinely Seq_type-keyed API: the generation path's
/// `unordered_map<Seq_type, string>` (B9) and the Len_proba machinery (G5/S4). What no longer
/// needs it is initialize_event(), which addresses the scenario maps -- those are keyed by
/// SeqTypeId, the only identity a non-legacy seq_type can have at all.
bool insertion_seq_type_str_to_enum(const Seq_type_String &seq_type_str, Seq_type &out)
{
    if (seq_type_str == "VD_ins_seq") { out = VD_ins_seq; return true; }
    if (seq_type_str == "DJ_ins_seq") { out = DJ_ins_seq; return true; }
    if (seq_type_str == "VJ_ins_seq") { out = VJ_ins_seq; return true; }
    return false;
}

/// The same conversion where there is no sensible way to continue without it.
///
/// Written as a return value rather than an out-parameter seeded with VD_ins_seq: that idiom
/// reads as "defaults to the VD junction", which is never what is meant -- the value is dead
/// on every path that does not throw.
Seq_type insertion_seq_type_or_throw(const Seq_type_String &seq_type_str, const char *where)
{
    Seq_type converted = VD_ins_seq;
    if (!insertion_seq_type_str_to_enum(seq_type_str, converted)) {
        throw std::runtime_error(std::string("Unknown insertion seq_type \"") + seq_type_str
                                 + "\" in " + where);
    }
    return converted;
}
} // namespace



Insertion::Insertion() : Insertion(VD_ins_seq)
{
    this->type = Event_type::Insertion_t;
    this->update_event_name();
}

Insertion::Insertion(Seq_type seq_type, pair<int, int> ins_range) : Insertion(seq_type)
{ //FIXME nonsense new

    int min_ins = std::min(ins_range.first, ins_range.second);
    int max_ins = std::max(ins_range.first, ins_range.second);

    this->type = Event_type::Insertion_t;
    this->len_max = max_ins;
    this->len_min = min_ins;

    for (int i = min_ins; i != max_ins + 1; ++i) {
        this->add_realization(i);
    }
    this->update_event_name();
}

Insertion::Insertion(Seq_type seq_type)
    : Rec_Event(Undefined_gene, Undefined_side),
      ins_seq_type(seq_type),
      proba_contribution(-1),
      previous_index(-1),
      insertions(-1),
      new_scenario_proba(-1),
      base_index(-1),
      new_index(-1),
      realization_index(-1),
      dinuc_updated_bound(NULL)
{
    this->type = Event_type::Insertion_t;
    //Keep the name and the enum in step from construction. An Insertion built directly from a
    //Seq_type used to carry an empty seq_type string until a caller happened to set one, which
    //left the two identities of the same fact disagreeing -- and the events_map is keyed by the
    //string.
    this->seq_type = EventUtils::seq_type_to_string(seq_type);
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > this->len_max) {
            this->len_max = (*iter).second.value_int;
        } else if ((*iter).second.value_int < this->len_min) {
            this->len_min = (*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Insertion::Insertion(Seq_type seq_type, unordered_map<string, Event_realization> &realizations) : Insertion(seq_type)
{
    this->event_realizations = realizations;

    this->type = Event_type::Insertion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > this->len_max) {
            this->len_max = (*iter).second.value_int;
        } else if ((*iter).second.value_int < this->len_min) {
            this->len_min = (*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Insertion::~Insertion()
{
    // TODO Auto-generated destructor stub
}

shared_ptr<Rec_Event> Insertion::copy()
{
    //TODO check this kind of copy for memory leak
    shared_ptr<Insertion> new_insertion_p =
            shared_ptr<Insertion>(new Insertion(this->ins_seq_type, this->event_realizations));
    new_insertion_p->priority = this->priority;
    new_insertion_p->nickname = this->nickname;
    new_insertion_p->fixed = this->fixed;
    new_insertion_p->update_event_name();
    new_insertion_p->set_event_identifier(this->event_index);
    new_insertion_p->set_seq_type(this->get_seq_type());
    new_insertion_p->set_seq_type_id(this->get_seq_type_id());
    return new_insertion_p;
}

bool Insertion::add_realization(int insertion_number)
{
    this->Rec_Event::add_realization(
            Event_realization(to_string(insertion_number), insertion_number, "", Int_Str(), this->size()));
    if (insertion_number > this->len_max) {
        this->len_max = insertion_number;
    } else if (insertion_number < this->len_min) {
        this->len_min = insertion_number;
    }
    this->update_event_name();
    return 0;
}

/**
 * @brief Context-based iterate() implementation
 *
 * Unpacks 5 context objects into legacy parameters and delegates
 * to the existing iterate() implementation.
 *
 * NOTE: Some const_casts are required because legacy iterate() signatures
 * accept non-const references even though they don't modify model data.
 * These casts are safe and will be eliminated when legacy interface is removed.
 */
void Insertion::iterate(
        QuerySequenceContext& query,
        const ModelContext& model,
        ScenarioContext& scenario,
        ExplorationContext& exploration,
        AccumulationContext& accumulation)
{
    base_index = exploration.index_map.get(this->event_index);
    proba_contribution = 1;

    //One junction, whoever its neighbours are. The three hardcoded seq_type comparisons this
    //replaces ran per scenario; the neighbour ids are resolved once at initialize_event().
    insertions = scenario.seq_offsets.get(right_neighbour_id, Five_prime)
                 - scenario.seq_offsets.get(left_neighbour_id, Three_prime) - 1;

    proba_contribution = iterate_common(proba_contribution, insertions, base_index, exploration.index_map,
                                        model.offset_map, model.model_parameters);

    if (proba_contribution != 0) {
        inserted_str.assign(insertions, int_undefined);
        //The VD and DJ arms used to re-derive this as
        //`event_realizations.at(to_string(insertions)).index` -- a string conversion and a hash
        //lookup in the hot loop, carrying its own FIXME. iterate_common() has already resolved
        //the same value; the two are pinned equal by test_insertion_iterate.cpp.
        new_index = base_index + realization_index;
        scenario.constructed_sequences.set_current(seq_type_id, &inserted_str);
        exploration.downstream_proba_map.set(seq_type_id, junction_length_best_proba_map.at(insertions),
                                             memory_layer_proba_map_junction);

        scenario.scenario_proba *= proba_contribution;
        //tmp_err_w_proba*=proba_contribution;
        (*dinuc_updated_bound) = upper_bound_per_ins.at(insertions);

        //Compute scenario downstream proba bound
        scenario_upper_bound_proba = exploration.compute_upper_bound(
            scenario.scenario_proba,
            current_downstream_proba_memory_layers
        );

        if (!exploration.should_prune(scenario_upper_bound_proba)) {
            Rec_Event::iterate_wrap_up(query, model, scenario, exploration, accumulation);
        }
    }
}

/*
 *This short method performs the iterate operations common to all Rec_event (modify index map and fetch realization probability)
 */
inline double Insertion::iterate_common(
        double scenario_proba, int insertions, int base_index, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        const Marginal_array_p &model_parameters_point)
{

    //insertions_str = to_string(insertions);
    //TODO just output proba contribution no need to take it as argument
    if (this->ordered_realization_map.count(insertions) > 0) {
        realization_index = this->ordered_realization_map.at(insertions).index;
        current_realizations_index_vec[0] = realization_index;
    } else {
        //discard out of range cases
        realization_index = INT32_MAX; //make sure the index called at the end is in the range
        //scenario_proba = 0;//discard the recombination scenario
        return 0;
    }

    /*	if (offset_map.count(this->name)!=0){
		for(vector<pair<const Rec_Event*,int>>::const_iterator jiter = offset_map.at(this->name).begin() ; jiter!= offset_map.at(this->name).end() ; jiter++){
			//modify index map using offset map
			base_index_map.at((*jiter).first->get_name()) += realization_index * ((*jiter).second);
		}
	}*/

    double prob = Rec_Event::iterate_common(
        realization_index, base_index, base_index_map, model_parameters_point);
    return scenario_proba * prob;
}

queue<int> Insertion::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        unordered_map<Seq_type, string> &constructed_sequences, mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            Seq_type seq_type = VD_ins_seq;
            if (insertion_seq_type_str_to_enum(this->seq_type, seq_type)) {
                constructed_sequences[seq_type] = string((*iter).second.value_int, 'I');
            }
            realization_queue.push((*iter).second.index);
            if (offset_map.count(this->get_name()) != 0) {
                for (vector<pair<shared_ptr<const Rec_Event>, int>>::const_iterator jiter =
                             offset_map.at(this->get_name()).begin();
                     jiter != offset_map.at(this->get_name()).end(); ++jiter) {
                    index_map.at((*jiter).first->get_name()) += (*iter).second.index * (*jiter).second;
                }
            }

            break;
        }
    }
    return realization_queue;
}

void Insertion::write2txt(ofstream &outfile)
{
    write2txt_legacy(outfile);
}

void Insertion::write2txt_legacy(ofstream &outfile)
{
    // Derive legacy gene_class from seq_type for backward compatibility
    // Use to_string() to match the exact strings expected by str2GeneClass()
    string legacy_gene_class;
    if (seq_type == "VD_ins_seq") legacy_gene_class = to_string(VD_genes);
    else if (seq_type == "DJ_ins_seq") legacy_gene_class = to_string(DJ_genes);
    else if (seq_type == "VJ_ins_seq") legacy_gene_class = to_string(VJ_genes);
    else legacy_gene_class = to_string(Undefined_gene);
    outfile << "#Insertion;" << legacy_gene_class << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Insertion::write2txt_v2(ofstream &outfile)
{
    outfile << "#Insertion;Undefined_gene;" << seq_type << ";" << event_side << ";" << priority << ";" << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Insertion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const Events_map &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>> &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map, Seq_type_str_p_map &constructed_sequences,
        Safety_bool_map &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{
    //An event carries its seq_type twice: as the serialization name and as the registry id
    //resolved from it. Everything downstream addresses the id, so a disagreement aliases this
    //event's whole scenario state onto another segment's keys and shows up as a wrong answer,
    //not as a failure -- the trap B2 hit with VJ. Checked once, here, where both are in hand.
    const SeqTypeRegistry &registry = constructed_sequences.registry();
    if (this->seq_type_id == kNoSeqType
        || static_cast<std::size_t>(this->seq_type_id) >= registry.total_count()
        || registry.name(this->seq_type_id) != this->seq_type) {
        throw runtime_error("Insertion " + this->name + ": seq_type \"" + this->seq_type
                            + "\" does not match the registry entry for its id");
    }

    downstream_proba_map.request_layer(this->seq_type_id);
    memory_layer_proba_map_junction = downstream_proba_map.claimed_layer(this->seq_type_id);

    //Resolve the junction's neighbours once, from the ordering rather than from the seq_type
    //name. This is what lets iterate() be topology-agnostic: a tandem-D D1D2_ins finds D1 and
    //D2 by the same two lookups that find V and D here.
    left_neighbour_id = registry.left_neighbor(this->seq_type_id);
    right_neighbour_id = registry.right_neighbor(this->seq_type_id);
    if (left_neighbour_id == kNoSeqType || right_neighbour_id == kNoSeqType) {
        throw runtime_error("Insertion " + this->name + " has no segment on one side: an "
                            "insertion is defined by the two segments it sits between");
    }

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map, downstream_proba_map,
                                      constructed_sequences, safety_set, error_rate_p, mismatches_list, seq_offsets,
                                      index_map);
}

void Insertion::add_to_marginals(long double scenario_proba, Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        updated_marginals[this->new_index] = scenario_proba;
    } else {
        updated_marginals[this->new_index] += scenario_proba;
    }
}

void Insertion::set_crude_upper_bound_proba(size_t base_index, size_t event_size, Marginal_array_p &marginal_array_p)
{

    size_t numb_realizations = this->size();
    upper_bound_per_ins.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        //Get the max proba for each number of insertions
        size_t real_index = (*iter).second.index;
        size_t j = 0;
        double max_proba = 0;

        while ((j * numb_realizations + real_index) < event_size) {
            if (marginal_array_p[base_index + ((j * numb_realizations + real_index))] > max_proba) {
                max_proba = marginal_array_p[base_index + ((j * numb_realizations + real_index))];
            }
            ++j;
        }
        upper_bound_per_ins[(*iter).second.value_int] = max_proba;
    }
}

void Insertion::initialize_crude_scenario_proba_bound(
        double &downstream_proba_bound, forward_list<double *> &updated_proba_list,
        const Events_map &events_map)
{
    this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
    this->updated_proba_bounds_list = updated_proba_list;
    this->event_upper_bound_proba = 0;
    shared_ptr<Rec_Event> dinuc_event_p;

    //TODO remove this and correct the way ordered realization map works
    ordered_realization_map.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter = (*this).event_realizations.begin();
         iter != (*this).event_realizations.end(); ++iter) {
        ordered_realization_map.emplace((*iter).second.value_int, (*iter).second);
    }

    //The switch this replaces converted ins_seq_type back into the very string the map is
    //keyed by. Dinucl_markov events are keyed with Undefined_side: the seq_type alone
    //identifies which junction they fill.
    const auto dinuc_entry = events_map.find(make_tuple(Dinuclmarkov_t, this->seq_type, Undefined_side));
    if (dinuc_entry != events_map.end()) {
        dinuc_event_p = dinuc_entry->second;
    }
    if (!dinuc_event_p) {
        throw runtime_error("Could not find associated Dinuclmarkov event for Insertion bounds");
    }
    double dinuc_upper_bound_proba = dinuc_event_p->get_upper_bound_proba();
    for (map<int, double>::iterator iter = upper_bound_per_ins.begin(); iter != upper_bound_per_ins.end(); ++iter) {
        //Compute joint upper bound of dinuc and insertion and store it as insertion upper bound
        (*iter).second *= pow(dinuc_upper_bound_proba, (*iter).first);
        if ((*iter).second > this->event_upper_bound_proba) {
            this->event_upper_bound_proba = (*iter).second;
        }
        //Only keep information about the dinucleotide probability to update the dinuc upperbound
        (*iter).second = pow(dinuc_upper_bound_proba, (*iter).first);
    }
    this->dinuc_updated_bound = dinuc_event_p->get_updated_ptr();
    //Remove the pointer from the list (otherwise dinuc upperbound is accounted for twice for events before insertion)
    updated_proba_list.remove(dinuc_updated_bound);

    //Apply the computed upper bound
    downstream_proba_bound *= event_upper_bound_proba;
}

OffsetDelta Insertion::get_offset_delta_bounds(SeqTypeId, Seq_side) const
{
    //An insertion writes no offsets at all: its span is *derived* from where its neighbours
    //already sit, which is exactly what makes the generic B6 rule possible.
    return {};
}

LengthContribution Insertion::get_length_contribution(SeqTypeId type_id) const
{
    if (type_id != this->seq_type_id || this->event_realizations.empty()) {
        return {};
    }
    int fewest = std::numeric_limits<int>::max();
    int most = std::numeric_limits<int>::min();
    for (const auto &[name, realization] : this->event_realizations) {
        (void)name;
        fewest = std::min(fewest, realization.value_int);
        most = std::max(most, realization.value_int);
    }
    return {fewest, most};
}

SeqConstructionRole Insertion::get_seq_construction_role(SeqTypeId type_id) const
{
    //Creates the segment, but as placeholders: Dinucl_markov supplies the nucleotides.
    return type_id == this->seq_type_id ? SeqConstructionRole::Creates : SeqConstructionRole::None;
}

OffsetRole Insertion::get_offset_role(SeqTypeId, Seq_side) const
{
    return OffsetRole::None;
}

bool Insertion::has_effect_on(Seq_type seq_type) const
{
    const string &heo_st = this->seq_type;
    if (heo_st == "VD_ins_seq") {
        return (seq_type == VJ_ins_seq || seq_type == VD_ins_seq);
    } else if (heo_st == "VJ_ins_seq") {
        return (seq_type == VJ_ins_seq);
    } else if (heo_st == "DJ_ins_seq") {
        return (seq_type == VJ_ins_seq || seq_type == DJ_ins_seq);
    } else {
        return false;
    }
}

void Insertion::iterate_initialize_Len_proba(Seq_type considered_junction, std::map<int, double> &length_best_proba_map,
                                             std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                             double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                             Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                             int &seq_len /*=0*/) const
{

    if (this->has_effect_on(considered_junction)) {

        base_index_map.set_current_layer(this->event_index, 0);
        base_index = base_index_map.get(this->event_index);

        //Insert sequence in the right constructed sequence. Still the enum: the whole
        //Len_proba machinery is Seq_type-keyed (Rec_Event::iterate_initialize_Len_proba), and
        //re-keying it is G5/S4's business, not B6's.
        const Seq_type seq_type = insertion_seq_type_or_throw(this->seq_type,
                                                              "iterate_initialize_Len_proba");

        for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {

            /*		//Update base index map
			for(forward_list<tuple<int,int,int>>::const_iterator jiter = memory_and_offsets.begin() ; jiter!=memory_and_offsets.end() ; ++jiter){
				//Get previous index for the considered event
				size_t previous_index = base_index_map.get(get<0>(*jiter),get<1>(*jiter)-1);
				//Update the index given the realization and the offset
				previous_index += iter->second.index *get<2>(*jiter);
				//Set the value
				base_index_map.set(get<0>(*jiter) , previous_index , get<1>(*jiter));
			}*/

            //Get the max proba for this realization (in case the event is child of another)
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()] > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index + i * this->size()];
                }
            }

            //Build an inserted sequence to let the Dinuc know about the number of insertions considered
            inserted_str.assign(iter->second.value_int, int_undefined);
            constructed_sequences.set_current(seq_type, &inserted_str);

            //Update the length and the probability within the recursive call
            Rec_Event::iterate_initialize_Len_proba_wrap_up(
                    considered_junction, length_best_proba_map, model_queue, scenario_proba * real_max_proba,
                    model_parameters_point, base_index_map, constructed_sequences, seq_len + (*iter).second.value_int);
        }
    } else {
        //Recursive call
        Rec_Event::iterate_initialize_Len_proba_wrap_up(considered_junction, length_best_proba_map, model_queue,
                                                        scenario_proba, model_parameters_point, base_index_map,
                                                        constructed_sequences, seq_len);
    }
}

void Insertion::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                           const Marginal_array_p &model_parameters_point, Index_map &base_index_map)
{
    //Still the enum, for the same reason as iterate_initialize_Len_proba above.
    const Seq_type seq_type = insertion_seq_type_or_throw(this->seq_type, "initialize_Len_proba_bound");

    //Scratch map for the junction length bound, which is still VDJ-hardcoded below;
    //see legacy_seq_type_registry(). It carries no ordering, so any event walked from here
    //that asks who its neighbours are gets kNoSeqType -- a landmine for B7.
    Seq_type_str_p_map constructed_sequences(legacy_seq_type_registry());

    junction_length_best_proba_map.clear();

    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        inserted_str.assign(iter->second.value_int, int_undefined);
        constructed_sequences.set_current(seq_type, &inserted_str);
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(seq_type, junction_length_best_proba_map, model_queue, init_proba,
                                                      model_parameters_point, base_index_map, constructed_sequences);
    }
}

void Insertion::update_event_name()
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
