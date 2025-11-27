/*
 * Insertion.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Quentin Marcou
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 analyze and model immune receptors
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
#include <igor/Core/SequenceTypes.h>

using namespace std;
Insertion::Insertion() : Insertion(Undefined_gene)
{
    this->type = Event_type::Insertion_t;
    this->update_event_name();
}

Insertion::Insertion(Gene_class genes, pair<int, int> ins_range) : Insertion(genes)
{ // FIXME nonsense new

    int min_ins = min(ins_range.first, ins_range.second);
    int max_ins = max(ins_range.first, ins_range.second);

    this->type = Event_type::Insertion_t;
    this->len_max = max_ins;
    this->len_min = min_ins;

    for (int i = min_ins; i != max_ins + 1; ++i) {
        this->add_realization(i);
    }
    this->update_event_name();
}

Insertion::Insertion(Gene_class gene)
    : Rec_Event(gene, Undefined_side),
      proba_contribution(-1),
      previous_index(-1),
      insertions(-1),
      new_scenario_proba(-1),
      base_index(-1),
      new_index(-1),
      realization_index(-1),
      dinuc_updated_bound(NULL),
      junction_type_id(-1),
      upstream_gene_type_id(-1),
      downstream_gene_type_id(-1)
{
    this->type = Event_type::Insertion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        if ((*iter).second.value_int > this->len_max) {
            this->len_max = (*iter).second.value_int;
        } else if ((*iter).second.value_int < this->len_min) {
            this->len_min = (*iter).second.value_int;
        }
    }
    this->update_event_name();
}

Insertion::Insertion(Gene_class gene, unordered_map<string, Event_realization> &realizations)
    : Insertion(gene)
{
    this->event_realizations = realizations;

    this->type = Event_type::Insertion_t;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
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
    // TODO check this kind of copy for memory leak
    shared_ptr<Insertion> new_insertion_p =
            shared_ptr<Insertion>(new Insertion(this->event_class, this->event_realizations));
    new_insertion_p->priority = this->priority;
    new_insertion_p->nickname = this->nickname;
    new_insertion_p->fixed = this->fixed;
    new_insertion_p->update_event_name();
    new_insertion_p->set_event_identifier(this->event_index);
    return new_insertion_p;
}

bool Insertion::add_realization(int insertion_number)
{
    this->Rec_Event::add_realization(Event_realization(
            to_string(insertion_number), insertion_number, "", Int_Str(), this->size()));
    if (insertion_number > this->len_max) {
        this->len_max = insertion_number;
    } else if (insertion_number < this->len_min) {
        this->len_min = insertion_number;
    }
    this->update_event_name();
    return 0;
}

void Insertion::iterate(
        double &scenario_proba, Downstream_scenario_proba_bound_map &downstream_proba_map,
        const string &sequence, const Int_Str &int_sequence, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        shared_ptr<Next_event_ptr> &next_event_ptr_arr, Marginal_array_p &updated_marginals_point,
        const Marginal_array_p &model_parameters_point,
        const unordered_map<Gene_class, vector<Alignment_data>> &allowed_realizations,
        Seq_type_str_p_map &constructed_sequences, Seq_offsets_map &seq_offsets,
        shared_ptr<Error_rate> &error_rate_p, map<size_t, shared_ptr<Counter>> &counters_list,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        Safety_bool_map &safety_set, Mismatch_vectors_map &mismatches_lists,
        double &seq_max_prob_scenario, double &proba_threshold_factor)
{
    base_index = base_index_map.at(this->event_index);
    new_scenario_proba = scenario_proba;
    proba_contribution = 1;

    // Use topology-based junction calculation
    // junction_type_id, upstream_gene_type_id, and downstream_gene_type_id
    // are set during initialize_event()

    // Calculate number of insertions: downstream_5' - upstream_3' - 1
    insertions = seq_offsets.at(downstream_gene_type_id, Five_prime)
            - seq_offsets.at(upstream_gene_type_id, Three_prime) - 1;

    proba_contribution =
            (*this).iterate_common(proba_contribution, insertions, base_index, base_index_map,
                                   offset_map, model_parameters_point);
    if (proba_contribution != 0) {
        inserted_str.assign(insertions, -1);
        new_index = base_index + this->event_realizations.at(to_string(insertions)).index;
        constructed_sequences[junction_type_id] = &inserted_str;
        downstream_proba_map.set_value(junction_type_id,
                                       junction_length_best_proba_map.at(insertions),
                                       memory_layer_proba_map_junction);
    }

    if (proba_contribution != 0) {
        new_scenario_proba *= proba_contribution;
        (*dinuc_updated_bound) = upper_bound_per_ins.at(insertions);

        // Compute scenario downstream proba bound
        scenario_upper_bound_proba = new_scenario_proba;
        // Multiply all downstream probas
        downstream_proba_map.multiply_all(scenario_upper_bound_proba,
                                          current_downstream_proba_memory_layers);

        if (scenario_upper_bound_proba >= (seq_max_prob_scenario * proba_threshold_factor)) {
            Rec_Event::iterate_wrap_up(
                    new_scenario_proba, downstream_proba_map, sequence, int_sequence,
                    base_index_map, offset_map, next_event_ptr_arr, updated_marginals_point,
                    model_parameters_point, allowed_realizations, constructed_sequences,
                    seq_offsets, error_rate_p, counters_list, events_map, safety_set,
                    mismatches_lists, seq_max_prob_scenario, proba_threshold_factor);
        }
    }
}

/*
 *This short method performs the iterate operations common to all Rec_event
 *(modify index map and fetch realization probability)
 */
inline double Insertion::iterate_common(
        double scenario_proba, int insertions, int base_index, Index_map &base_index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        const Marginal_array_p &model_parameters_point)
{

    // insertions_str = to_string(insertions);
    // TODO just output proba contribution no need to take it as argument
    if (this->ordered_realization_map.count(insertions) > 0) {
        realization_index = this->ordered_realization_map.at(insertions).index;
        current_realizations_index_vec[0] = realization_index;
    } else {
        // discard out of range cases
        realization_index = INT32_MAX; // make sure the index called at the end is in the range
        // scenario_proba = 0;//discard the recombination scenario
        return 0;
    }

    /*	if (offset_map.count(this->name)!=0){
                  for(vector<pair<const Rec_Event*,int>>::const_iterator jiter =
     offset_map.at(this->name).begin() ; jiter!= offset_map.at(this->name).end()
     ; jiter++){
                          //modify index map using offset map
                          base_index_map.at((*jiter).first->get_name()) +=
     realization_index * ((*jiter).second);
                  }
          }*/

    double prob = Rec_Event::iterate_common(realization_index, base_index, base_index_map,
                                            model_parameters_point);
    return scenario_proba * prob;
}

queue<int> Insertion::draw_random_realization(
        const Marginal_array_p &model_marginals_p, unordered_map<Rec_Event_name, int> &index_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        std::unordered_map<int, std::string> &constructed_sequences, mt19937_64 &generator) const
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand = distribution(generator);
    double prob_count = 0;
    queue<int> realization_queue;
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
        if (prob_count >= rand) {
            // Use topology-based junction type assignment
            constructed_sequences[junction_type_id] = string((*iter).second.value_int, 'I');

            realization_queue.push((*iter).second.index);
            if (offset_map.count(this->get_name()) != 0) {
                for (vector<pair<shared_ptr<const Rec_Event>, int>>::const_iterator jiter =
                             offset_map.at(this->get_name()).begin();
                     jiter != offset_map.at(this->get_name()).end(); ++jiter) {
                    index_map.at((*jiter).first->get_name()) +=
                            (*iter).second.index * (*jiter).second;
                }
            }

            break;
        }
    }
    return realization_queue;
}

void Insertion::write2txt(ofstream &outfile)
{
    outfile << "#Insertion;" << event_class << ";" << event_side << ";" << priority << ";"
            << nickname << endl;
    for (unordered_map<string, Event_realization>::const_iterator iter = event_realizations.begin();
         iter != event_realizations.end(); ++iter) {
        outfile << "%" << (*iter).second.value_int << ";" << (*iter).second.index << endl;
    }
}

void Insertion::initialize_event(
        unordered_set<Rec_Event_name> &processed_events,
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map,
        const unordered_map<Rec_Event_name, vector<pair<shared_ptr<const Rec_Event>, int>>>
                &offset_map,
        Downstream_scenario_proba_bound_map &downstream_proba_map,
        Seq_type_str_p_map &constructed_sequences, Safety_bool_map &safety_set,
        shared_ptr<Error_rate> error_rate_p, Mismatch_vectors_map &mismatches_list,
        Seq_offsets_map &seq_offsets, Index_map &index_map)
{

    // Determine junction type from nickname or event_class
    auto &registry = SequenceTypeRegistry::get_instance();

    // Try to get type from nickname first (for tandem D support: "VD1_ins", "D1D2_ins", etc.)
    int type_id = registry.try_get_type_id(this->nickname + "_seq");

    // If not found, use event_class for backward compatibility
    if (type_id < 0) {
        switch (this->event_class) {
        case VD_genes:
            type_id = SequenceTypeRegistry::VD_INS_SEQ;
            break;
        case DJ_genes:
            type_id = SequenceTypeRegistry::DJ_INS_SEQ;
            break;
        case VJ_genes:
            type_id = SequenceTypeRegistry::VJ_INS_SEQ;
            break;
        default:
            throw runtime_error("Unknown Gene_class for Insertion: "
                                + std::to_string(this->event_class));
        }
    }

    junction_type_id = type_id;

    // Get upstream and downstream gene types from registry topology
    upstream_gene_type_id = registry.get_junction_upstream(junction_type_id);
    downstream_gene_type_id = registry.get_junction_downstream(junction_type_id);

    // Request memory layer for this junction type
    downstream_proba_map.request_memory_layer(junction_type_id);
    memory_layer_proba_map_junction = downstream_proba_map.get_current_memory_layer(junction_type_id);

    this->Rec_Event::initialize_event(processed_events, events_map, offset_map,
                                      downstream_proba_map, constructed_sequences, safety_set,
                                      error_rate_p, mismatches_list, seq_offsets, index_map);
}

void Insertion::update_event_internal_probas(const Marginal_array_p &marginal_array,
                                             const unordered_map<Rec_Event_name, int> &index_map)
{
    // Update junction type from nickname if not yet done (for tandem D support)
    if (junction_type_id < 0 && !this->nickname.empty()) {
        auto &registry = SequenceTypeRegistry::get_instance();
        int type_id = registry.try_get_type_id(this->nickname + "_seq");
        if (type_id >= 0) {
            junction_type_id = type_id;
            upstream_gene_type_id = registry.get_junction_upstream(junction_type_id);
            downstream_gene_type_id = registry.get_junction_downstream(junction_type_id);
        }
    }

    // If still not set, fall back to event_class defaults
    if (junction_type_id < 0) {
        switch (this->event_class) {
        case VD_genes:
            junction_type_id = SequenceTypeRegistry::VD_INS_SEQ;
            upstream_gene_type_id = SequenceTypeRegistry::V_GENE_SEQ;
            downstream_gene_type_id = SequenceTypeRegistry::D_GENE_SEQ;
            break;
        case DJ_genes:
            junction_type_id = SequenceTypeRegistry::DJ_INS_SEQ;
            upstream_gene_type_id = SequenceTypeRegistry::D_GENE_SEQ;
            downstream_gene_type_id = SequenceTypeRegistry::J_GENE_SEQ;
            break;
        case VJ_genes:
            junction_type_id = SequenceTypeRegistry::VJ_INS_SEQ;
            upstream_gene_type_id = SequenceTypeRegistry::V_GENE_SEQ;
            downstream_gene_type_id = SequenceTypeRegistry::J_GENE_SEQ;
            break;
        default:
            // Default to VD for unknown
            junction_type_id = SequenceTypeRegistry::VD_INS_SEQ;
            upstream_gene_type_id = SequenceTypeRegistry::V_GENE_SEQ;
            downstream_gene_type_id = SequenceTypeRegistry::D_GENE_SEQ;
            break;
        }
    }
}

void Insertion::add_to_marginals(long double scenario_proba,
                                 Marginal_array_p &updated_marginals) const
{
    if (viterbi_run) {
        updated_marginals[this->new_index] = scenario_proba;
    } else {
        updated_marginals[this->new_index] += scenario_proba;
    }
}

void Insertion::set_crude_upper_bound_proba(size_t base_index, size_t event_size,
                                            Marginal_array_p &marginal_array_p)
{

    size_t numb_realizations = this->size();
    upper_bound_per_ins.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        // Get the max proba for each number of insertions
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
        const unordered_map<tuple<Event_type, int, Seq_side>, shared_ptr<Rec_Event>> &events_map)
{
    this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
    this->updated_proba_bounds_list = updated_proba_list;
    this->event_upper_bound_proba = 0;
    shared_ptr<Rec_Event> dinuc_event_p;

    // TODO remove this and correct the way ordered realization map works
    ordered_realization_map.clear();
    for (unordered_map<string, Event_realization>::const_iterator iter =
                 (*this).event_realizations.begin();
         iter != (*this).event_realizations.end(); ++iter) {
        ordered_realization_map.emplace((*iter).second.value_int, (*iter).second);
    }

    switch (this->event_class) {
        // TODO be careful in case there is both VDJ and VD/DJ (however this should
        // not happen)

    case VD_genes:
        if (events_map.count(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VD_genes, Undefined_side))) {
            dinuc_event_p = events_map.at(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VD_genes, Undefined_side));
        } else if (events_map.count(tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VDJ_genes,
                                                                     Undefined_side))) {
            dinuc_event_p = events_map.at(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VDJ_genes, Undefined_side));
        }
        break;

    case VJ_genes:
        if (events_map.count(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VJ_genes, Undefined_side))) {
            dinuc_event_p = events_map.at(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VJ_genes, Undefined_side));
        }
        break;

    case DJ_genes:
        if (events_map.count(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, DJ_genes, Undefined_side))) {
            dinuc_event_p = events_map.at(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, DJ_genes, Undefined_side));
        } else if (events_map.count(tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VDJ_genes,
                                                                     Undefined_side))) {
            dinuc_event_p = events_map.at(
                    tuple<Event_type, int, Seq_side>(Dinuclmarkov_t, VDJ_genes, Undefined_side));
        }
        break;
    default:
        throw runtime_error("Unknown Gene class for insertion in "
                            "initialize_scenario_proba_bound()");
        break;
    }
    double dinuc_upper_bound_proba = dinuc_event_p->get_upper_bound_proba();
    for (map<int, double>::iterator iter = upper_bound_per_ins.begin();
         iter != upper_bound_per_ins.end(); ++iter) {
        // Compute joint upper bound of dinuc and insertion and store it as
        // insertion upper bound
        (*iter).second *= pow(dinuc_upper_bound_proba, (*iter).first);
        if ((*iter).second > this->event_upper_bound_proba) {
            this->event_upper_bound_proba = (*iter).second;
        }
        // Only keep information about the dinucleotide probability to update the
        // dinuc upperbound
        (*iter).second = pow(dinuc_upper_bound_proba, (*iter).first);
    }
    this->dinuc_updated_bound = dinuc_event_p->get_updated_ptr();
    // Remove the pointer from the list (otherwise dinuc upperbound is accounted
    // for twice for events before insertion)
    updated_proba_list.remove(dinuc_updated_bound);

    // Apply the computed upper bound
    downstream_proba_bound *= event_upper_bound_proba;
}

bool Insertion::has_effect_on(int seq_type) const
{
    // An insertion has effect on a junction if this insertion IS that junction
    // or if they share a common endpoint (for length-proba bound calculations)

    // Direct match: this insertion is the junction being queried
    if (seq_type == junction_type_id) {
        return true;
    }

    // For backward compatibility with VJ_ins_seq which shares endpoints with VD and DJ
    // Check if this junction shares an endpoint with the queried junction
    const auto &registry = SequenceTypeRegistry::get_instance();

    if (registry.is_junction_type(seq_type)) {
        int query_upstream = registry.get_junction_upstream(seq_type);
        int query_downstream = registry.get_junction_downstream(seq_type);

        // If our upstream matches query's upstream or downstream,
        // or our downstream matches query's upstream or downstream
        if (upstream_gene_type_id == query_upstream
            || upstream_gene_type_id == query_downstream
            || downstream_gene_type_id == query_upstream
            || downstream_gene_type_id == query_downstream) {
            return true;
        }
    }

    return false;
}

void Insertion::iterate_initialize_Len_proba(
        int considered_junction, std::map<int, double> &length_best_proba_map,
        std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
        const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
        Seq_type_str_p_map &constructed_sequences, int &seq_len /*=0*/) const
{

    if (this->has_effect_on(considered_junction)) {

        base_index = base_index_map.at(this->event_index, 0);

        // Use topology-based junction type
        int seq_type = junction_type_id;

        for (unordered_map<string, Event_realization>::const_iterator iter =
                     this->event_realizations.begin();
             iter != this->event_realizations.end(); ++iter) {

            // Get the max proba for this realization (in case the event is child of
            // another)
            double real_max_proba = 0;
            for (size_t i = 0; i != this->event_marginal_size / this->size(); ++i) {
                if (model_parameters_point[base_index + (*iter).second.index + i * this->size()]
                    > real_max_proba) {
                    real_max_proba = model_parameters_point[base_index + (*iter).second.index
                                                            + i * this->size()];
                }
            }

            // Build an inserted sequence to let the Dinuc know about the number of
            // insertions considered
            inserted_str.assign(iter->second.value_int, -1);
            constructed_sequences[seq_type] = &inserted_str;

            // Update the length and the probability within the recursive call
            Rec_Event::iterate_initialize_Len_proba_wrap_up(
                    considered_junction, length_best_proba_map, model_queue,
                    scenario_proba * real_max_proba, model_parameters_point, base_index_map,
                    constructed_sequences, seq_len + (*iter).second.value_int);
        }
    } else {
        // Recursive call
        Rec_Event::iterate_initialize_Len_proba_wrap_up(
                considered_junction, length_best_proba_map, model_queue, scenario_proba,
                model_parameters_point, base_index_map, constructed_sequences, seq_len);
    }
}

void Insertion::initialize_Len_proba_bound(queue<shared_ptr<Rec_Event>> &model_queue,
                                           const Marginal_array_p &model_parameters_point,
                                           Index_map &base_index_map)
{
    // Use topology-based junction type
    int seq_type = junction_type_id;

    Seq_type_str_p_map constructed_sequences(SequenceTypeRegistry::get_instance().size());

    junction_length_best_proba_map.clear();

    for (unordered_map<string, Event_realization>::const_iterator iter =
                 this->event_realizations.begin();
         iter != this->event_realizations.end(); ++iter) {
        inserted_str.assign(iter->second.value_int, -1);
        constructed_sequences[seq_type] = &inserted_str;
        double init_proba = 1.0;
        this->Rec_Event::iterate_initialize_Len_proba(
                seq_type, junction_length_best_proba_map, model_queue, init_proba,
                model_parameters_point, base_index_map, constructed_sequences);
    }
}
