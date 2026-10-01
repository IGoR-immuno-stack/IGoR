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
      realization_index(-1)
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
    insertions = scenario.seq_offsets.get(get_right_adjacent_id(), Five_prime)
                 - scenario.seq_offsets.get(get_left_adjacent_id(), Three_prime) - 1;

    proba_contribution = iterate_common(proba_contribution, insertions, base_index, exploration.index_map,
                                        model.offset_map, model.model_parameters);

    if (proba_contribution != 0) {
        //The VD and DJ arms used to re-derive this as
        //`event_realizations.at(to_string(insertions)).index` -- a string conversion and a hash
        //lookup in the hot loop, carrying its own FIXME. iterate_common() has already resolved
        //the same value; the two are pinned equal by test_insertion_iterate.cpp.
        new_index = base_index + realization_index;

        //Where the junction sits, and therefore how long it is. This event decides both and
        //creates neither the sequence nor a placeholder for it: under O12's decision (a') the
        //Dinucl_markov that follows creates the segment from exactly these two offsets. The
        //partially-constructed state -- a segment of int_undefined waiting for its filler --
        //stops existing anywhere in the program.
        //
        //Recording the offsets is the repair itself (R3, plan section 6.9). The old comment
        //beside get_offset_delta_bounds() ran two claims together: *the span is derived from
        //where the neighbours already sit* is true and is B6's whole point, *and therefore no
        //offsets are recorded* was the defect. A consumer is supposed to ask the segment where
        //it is rather than know which event produced it -- that is what A0 exists for.
        //
        //An empty junction is `3' == 5' - 1`, a pair of offsets rather than their absence,
        //which is the encoding B10 needs and the one the length below reads straight back.
        const Seq_Offset junction_five_prime =
                scenario.seq_offsets.get(get_left_adjacent_id(), Three_prime) + 1;
        scenario.seq_offsets.set(seq_type_id, Five_prime, junction_five_prime,
                                 memory_layer_offset_fivep);
        scenario.seq_offsets.set(seq_type_id, Three_prime, junction_five_prime + insertions - 1,
                                 memory_layer_offset_threep);

        const JunctionBound &junction = junction_bound(kEnclosingJunction);
        //One descent, where the other 18 consumption sites in Gene_choice and Deletion took two
        //(§6.10 finding 6). Unlike them this one has never had a guard, so an unreachable
        //distance throws rather than discarding the branch; turning that into a discard is a
        //behaviour change and belongs with the other Insertion defects in phase R.
        //
        //The table is conditioned on this event's realization (R13): read its profile for the
        //length chosen, at the gap between the neighbours -- which is where it stood before the
        //insertion, since an insertion fills the gap rather than moving its ends.
        const std::optional<double> junction_bound_proba =
                junction.profile_for(realization_index).best_for(insertions);
        if (not junction_bound_proba) {
            throw out_of_range("Insertion " + this->name + ": no junction bound for "
                               + to_string(insertions) + " insertions");
        }
        exploration.downstream_proba_map.set(junction.proba_key(), *junction_bound_proba,
                                             junction.memory_layer());

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
        SafetyMatrix &safety_set, shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
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

    //Both ends of the junction, because this event places both (O12 (a'), R3). Nothing stands
    //beneath them today -- no other event writes an insertion's offsets -- so the claimed
    //layer is 0 and the write lands where a set_current() would have. Requesting it anyway is
    //what makes the layer *owned*, which is the rule the harness holds every other event to.
    seq_offsets.request_layer(this->seq_type_id, Five_prime);
    this->memory_layer_offset_fivep = seq_offsets.claimed_layer(this->seq_type_id, Five_prime);
    seq_offsets.request_layer(this->seq_type_id, Three_prime);
    this->memory_layer_offset_threep = seq_offsets.claimed_layer(this->seq_type_id, Three_prime);

    //The neighbours come from Model_Parms::finalize(), which is the one place that reads the
    //ordering. A tandem-D D1D2_ins is told about D1 and D2 by the same pass that tells this
    //event about V and D.
    if (get_left_adjacent_id() == kNoSeqType || get_right_adjacent_id() == kNoSeqType) {
        throw runtime_error("Insertion " + this->name + " has no segment on one side: an "
                            "insertion is defined by the two segments it sits between, and "
                            "Model_Parms::finalize() must have resolved them");
    }

    //The junction this insertion fills -- the gap between the two segments it sits between, which
    //is the same pair iterate() measures the observed distance across. Resolved here and never
    //again: the profile is reached through the handle, not by asking which span this is.
    junction_bound(kEnclosingJunction)
            .resolve(SegmentSpan::gap(get_left_adjacent_id(), get_right_adjacent_id()), this->seq_type_id,
                     downstream_proba_map.claimed_layer(this->seq_type_id), JunctionBound::Fold::Yes);

    require_dinucl_markov(events_map);

    //iterate_common() finds a realization by its length. Built here rather than in the crude
    //bound initialization it used to share a function with, which R14 deleted: the map has
    //nothing to do with a bound. TODO the realization map itself could be keyed by length.
    ordered_realization_map.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter = this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        ordered_realization_map.emplace(iter->second.value_int, iter->second);
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

void Insertion::require_dinucl_markov(const Events_map &events_map) const
{
    //Dinucl_markov events are keyed with Undefined_side: the seq_type alone identifies which
    //junction they fill. Since O12 (a') it is the Dinucl_markov that creates the segment this
    //event only places, so an Insertion without one leaves a junction nobody fills.
    if (events_map.find(make_tuple(Dinuclmarkov_t, this->seq_type, Undefined_side)) == events_map.end()) {
        throw runtime_error("Insertion " + this->name + ": no Dinucl_markov event fills " + this->seq_type);
    }
}

OffsetDelta Insertion::get_offset_delta_bounds(SeqTypeId, Seq_side) const
{
    //Nothing, and for a reason that is *not* "this event has no offsets" -- it writes both of
    //them (R3). This query asks how far an end can still be *shifted* from where it already
    //stands, and an insertion never shifts one: its span is derived from where its neighbours
    //already sit, which is what makes the generic B6 rule possible, and the realization set's
    //effect on the 3' end is get_length_contribution()'s answer rather than this one's.
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

SeqConstructionRole Insertion::get_seq_construction_role(SeqTypeId) const
{
    //None, for any seq_type including its own. O12 (a'): this event decides the junction's
    //length and position, and the Dinucl_markov that follows creates the segment. Saying
    //`Creates` here while leaving every position undetermined was the placeholder regime.
    return SeqConstructionRole::None;
}

OffsetRole Insertion::get_offset_role(SeqTypeId type_id, Seq_side) const
{
    //Both ends, since iterate() derives and writes both. Plain `Creates`: the Anchors/Derives
    //split first proposed for this is the wrong answer, because whether a span is zero-width
    //is a *consumer* question and not a property of the segment (plan section 2.5).
    return type_id == this->seq_type_id ? OffsetRole::Creates : OffsetRole::None;
}

bool Insertion::affects_length_of(SegmentSpan span) const
{
    //An insertion segment lies strictly inside every span that brackets it, and adds its
    //realization's length to each. S4b generalises this off the enum, once the traversal
    //carries the registry ordering.
    const Seq_type junction = legacy_junction_of(span);
    const string &heo_st = this->seq_type;
    if (heo_st == "VD_ins_seq") {
        return (junction == VJ_ins_seq || junction == VD_ins_seq);
    } else if (heo_st == "VJ_ins_seq") {
        return (junction == VJ_ins_seq);
    } else if (heo_st == "DJ_ins_seq") {
        return (junction == VJ_ins_seq || junction == DJ_ins_seq);
    } else {
        return false;
    }
}

int Insertion::length_delta(const Event_realization &realization) const
{
    //The number of nucleotides inserted: the insertion creates the junction segment.
    return realization.value_int;
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
