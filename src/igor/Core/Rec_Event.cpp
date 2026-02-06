/*
 * Rec_Event.cpp
 *
 *  Created on: 3 nov. 2014
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

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Counter.h>

using namespace std;

//std::ofstream log_file(std::string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/one_seq_comp/logs.txt"));

Rec_Event::Rec_Event(Gene_class gene , Seq_side side ): priority(0) , event_class(gene) , event_side(side) , name("Undefined_event_name") ,len_min(INT16_MAX) , len_max(INT16_MIN) , type(Undefined_t), event_index(INT16_MIN) , updated(false),fixed(false) , current_realizations_index_vec(vector<int>()) , scenario_downstream_upper_bound_proba(-1),event_upper_bound_proba(-1),scenario_upper_bound_proba(-1),current_realization_index(nullptr){} //FIXME why does this exist? anyway fix initilization


Rec_Event::Rec_Event(Gene_class gene , Seq_side side , map<string , Event_realization>& realizations): Rec_Event(gene,side)  {
	this->event_realizations = realizations;
}

Rec_Event::Rec_Event(): Rec_Event( Undefined_gene , Undefined_side ) {}



//TODO see this later
/*
Rec_Event::Rec_Event(list<Event_realization> realizations_list){
	for(list<Event_realization>::const_iterator iter = realizations_list.begin() ; iter!=realizations_list.end() ; iter++){
		this->add_realization((*iter));
	}
}


Rec_Event::Rec_Event(list<Event_realization> realization_list, int new_priority ) : Rec_Event(realization_list) {
	this->priority = new_priority;
}
*/


Rec_Event::~Rec_Event() {
	// TODO Auto-generated destructor stub
}

bool Rec_Event::operator ==(const Rec_Event& other)const {
	if(this->get_type() != other.get_type()) return 0;
	if( this->event_class != other.event_class) return 0;
	if( this->event_side != other.event_side) return 0;
	if( this->priority != other.priority) return 0;
	if( this->event_realizations.size() != other.event_realizations.size()) return 0;
	for(map< string,Event_realization >::const_iterator iter = this->event_realizations.begin() ; iter != this->event_realizations.end() ; ++iter){
		if(other.event_realizations.count((*iter).first) != 1 ) return 0;
	}
	return 1;
}

void Rec_Event::update_event_name(){
	this->name = string() + this->type +string("_")+ this->event_class + string("_") + this->event_side + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

void Rec_Event::add_realization(const Event_realization& realization){
	this->event_realizations.insert( make_pair ((realization).name,(realization)) );
	this->update_event_name();
}



bool Rec_Event::set_priority(int new_priority){
	this->priority = new_priority;
	this->update_event_name();
	return 1;
}


int Rec_Event::size()const{
	return event_realizations.size();
}

void Rec_Event::set_event_identifier(size_t identifier){
	this->event_index = identifier;
}

int Rec_Event::get_event_identifier() const {
	return event_index;
}


void Rec_Event::iterate_wrap_up(double& scenario_proba, Downstream_scenario_proba_bound_map& downstream_proba_bound_map, const std::string& sequence, const Int_Str& int_sequence, Index_map& index_map, const std::map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& map_of_events, std::shared_ptr<Next_event_ptr>& next_event_ptr_str, Marginal_array_p& updated_marginal_array_p, const Marginal_array_p& pdf_marginal_array_p, const std::map<Gene_class , std::vector<Alignment_data>>& alignments, Seq_type_str_p_map& constructed_sequences, Seq_offsets_map& seq_offsets, std::shared_ptr<Error_rate>& error_rate_p, std::map<size_t,std::shared_ptr<Counter>>& counters_list, const std::map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& events_map, Safety_bool_map& safety_bool_map, Mismatch_vectors_map& mismatch_vectors_map, double& seq_likelihood, double& model_log_likelihood){


	double scenario_error_w_proba = error_rate_p->compare_sequences_error_prob(scenario_proba, sequence, constructed_sequences,seq_offsets,events_map,mismatch_vectors_map, seq_likelihood, model_log_likelihood);
	scenario_error_w_proba *= scenario_proba;

	//Check for 0 proba
	if(scenario_error_w_proba==0){return;}

	//Update the marginals
	std::shared_ptr<const Rec_Event> parent_event_p = map_of_events.at(this->name).at(0).first;//TODO check this line
	//Loop over the events
	//while(parent_event_p!=nullptr){
	while(true){
		parent_event_p->add_to_marginals(scenario_error_w_proba,updated_marginal_array_p);
		if(parent_event_p->get_name() == "root" or parent_event_p->get_name() == "Root"){break;}
		//parent_event_p = map_of_events.at(parent_event_p->event_name).at(0).first;
		parent_event_p = events_map.at(std::tuple<Event_type,Gene_class,Seq_side>(parent_event_p->get_type(),parent_event_p->get_class(),parent_event_p->get_side()));
	}


	//Add to the counters
	for(std::map<size_t,std::shared_ptr<Counter>>::iterator iter=counters_list.begin(); iter!=counters_list.end(); ++iter){
		iter->second->count_scenario(scenario_proba,scenario_error_w_proba,sequence,constructed_sequences,seq_offsets,events_map,mismatch_vectors_map);
	}


	seq_likelihood+=scenario_error_w_proba;

	if (counters_list.count(BEST_SCENARIOS_COUNTER) > 0){
		double proba = error_rate_p->compare_sequences_error_prob(scenario_proba, sequence, constructed_sequences,seq_offsets,events_map,mismatch_vectors_map, seq_likelihood, model_log_likelihood);
		if (proba > 0) {
			model_log_likelihood+=scenario_error_w_proba*log10(scenario_error_w_proba);
		}
	}


}



void Rec_Event::initialize_event( unordered_set<Rec_Event_name>& processed_events , const map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , const map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>>& offset_map , Downstream_scenario_proba_bound_map& downstream_proba_map , Seq_type_str_p_map& constructed_sequences , Safety_bool_map& safety_set , shared_ptr<Error_rate> error_rate_p , Mismatch_vectors_map& mismatches_list , Seq_offsets_map& seq_offsets , Index_map& index_map){
	//No action performed on the event by default if the method is not overloaded
	//Need to call Rec_Event::initialize_event() to apply these common actions when the method is overloaded
	current_realizations_index_vec.push_back(-1);

	if(offset_map.count(this->get_name()) != 0){
		const vector<pair<shared_ptr<const Rec_Event>,int>>& offset_vector = offset_map.at(this->get_name());
		for(vector<pair<shared_ptr<const Rec_Event>,int>>::const_iterator iter = offset_vector.begin() ; iter != offset_vector.end() ; ++iter){
			//Request memory layer
			int event_identitfier = (*iter).first->get_event_identifier();
			index_map.request_memory_layer(event_identitfier);
			memory_and_offsets.emplace_front( event_identitfier , index_map.get_current_memory_layer(event_identitfier) , (*iter).second);
		}
	}

	downstream_proba_map.get_all_current_memory_layer(current_downstream_proba_memory_layers);

	processed_events.emplace(this->name);
	return;
}

void Rec_Event::ind_normalize(Marginal_array_p& marginal_array_p , size_t base_index) const{
	double sum_marginals = 0;
	for(int i =0 ; i != this->size() ; ++i){
		sum_marginals+= marginal_array_p[base_index + i];
	}
	if(sum_marginals!=0){
		for(int i =0 ; i != this->size() ; ++i){
			marginal_array_p[base_index + i] /= sum_marginals;
		}
	}
}


void Rec_Event::set_crude_upper_bound_proba( size_t base_index , size_t event_size , Marginal_array_p& marginal_array_p){

	double max_proba = 0;
	for(size_t i = 0 ; i!= event_size ; ++i){
		if(marginal_array_p[base_index + i] > max_proba){
			max_proba = marginal_array_p[base_index + i];
		}
	}
	this->event_upper_bound_proba = max_proba;
}

void Rec_Event::set_upper_bound_proba(double proba){
	this->event_upper_bound_proba = proba;
}

/**
 * Does nothing since in general events will not need to perform any operation on the marginal probabilities
 */
void Rec_Event::update_event_internal_probas(const Marginal_array_p& marginal_array , const map<Rec_Event_name,int>& index_map){
	//Do nothing
}

/*
 * This method initialize the scenario probability upper bound for each event
 * The point is to compute the upper bound probability (given the model) of the scenario for each event
 * This allows to discard scenarios with too low probability at early stages
 */
void Rec_Event::initialize_crude_scenario_proba_bound(double& downstream_proba_bound , forward_list<double*>& updated_proba_list , const map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
	this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
	this->updated_proba_bounds_list = updated_proba_list;
	if(!this->is_updated()){
		downstream_proba_bound*=this->event_upper_bound_proba;
	}
	else{
		throw logic_error("Updated events should overload Rec_event::initialize_scenario_proba_bound()");
	}
}


/*
 * Description??
 */
double* Rec_Event::get_updated_ptr(){
	throw logic_error("Updated events should overload Rec_event::get_updated_ptr()");
}


/*
 * Updates the value of scenario_upper_bound_proba according to the error weighted scenario and the upper bound of downstream scenarios
 */
void Rec_Event::compute_crude_upper_bound_scenario_proba( double& tmp_err_w_proba ) {
	scenario_upper_bound_proba = tmp_err_w_proba * scenario_downstream_upper_bound_proba;
	for (forward_list<double*>::const_iterator iter = updated_proba_bounds_list.begin() ; iter != updated_proba_bounds_list.end() ; ++iter){
		scenario_upper_bound_proba*=(*(*iter));
	}
}


void Rec_Event::iterate_initialize_Len_proba(Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences ) const{
	int seq_len = 0;
	this->iterate_initialize_Len_proba(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
}

/*
 * Called when iterating over all possible scenarios during initialization
 * Fills up the length-max_proba_bound for a given junction , and links the call to iterate_initialize_len_proba for two events
 *
 * TODO constructed sequences should not be used but it is useful to compute the dinucl contribution
 */
void Rec_Event::iterate_initialize_Len_proba_wrap_up(Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>> model_queue , double scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int seq_len ) const {

	if(not model_queue.empty()){
		std::shared_ptr<Rec_Event> next_event_p = model_queue.front();
		model_queue.pop();
		//TODO fix this and find a way not to loop over all events
		//if(next_event_p->has_effect_on(considered_junction)){
			// Explore realizations of this event
			next_event_p->iterate_initialize_Len_proba(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
		//}
		//else{
			// If this event has no effect on the junction skip it using a recursive call
			//next_event_p->iterate_initialize_Len_proba_wrap_up(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
		//}
	}
	else{
		// When all events with an effect on the junction have been processed update the length-proba map
		if(length_best_proba_map.count(seq_len)>0){
			if(scenario_proba>length_best_proba_map.at(seq_len)){
				//Keep the best proba for each length
				length_best_proba_map.at(seq_len) = scenario_proba;
			}
		}
		else{
			length_best_proba_map[seq_len] = scenario_proba;
		}
	}

}

