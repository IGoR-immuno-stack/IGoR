/*
 * Model_Parms.cpp
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
 *      This class defines a Model_Parms object. In our framework this will be used
 *      to easily define the inferred model, by using a directed graph structure.
 *      This helps defining conditional marginals in the model.
 */

#include <igor/Model/Legacy/Model_Parms.h>
#include <igor/Alignment/Legacy/Aligner.h>
#include <map>
#include <set>
#include <tuple>

#include <cassert>
#include <igor/Model/Legacy/gene_to_seqtype_migr.h>
#include <igor/Model/Legacy/EventUtils.h>
#include <igor/Model/Legacy/LegacyFormat.h>
#include <igor/Model/Legacy/ModelJson.h>

#include <nlohmann/json.hpp>

namespace igor::model::legacy {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

using namespace std;

/*
 * Default constructor, creates an empty Model_Parms
 */
Model_Parms::Model_Parms()
{
    /*	this->events = *(new list <shared_ptr<Rec_Event>>()); //FIXME nonsense new
	this->edges = *(new unordered_map <Rec_Event_name , Adjacency_list >());*/
}

/*
 * Construct a Model_Parms from a list of event and initialize connectivities/Adjacency_lists
 */
Model_Parms::Model_Parms(list<shared_ptr<Rec_Event>> event_list)
{
    for (const shared_ptr<Rec_Event> &event : event_list) {
        this->add_event(event);
    }
}

/*
 * Provide a deep copy constructor of the Model_Parms object
 * This method will make a copy of all the contained events and error rate as well
 */
Model_Parms::Model_Parms(const Model_Parms &other)
{
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = other.events.begin(); iter != other.events.end(); ++iter) {
        this->events.push_back((*iter)->copy());
    }

    //Each edge's events are found again by name, which is unique in a model. They used to be
    //found with Rec_Event::operator==, which ignored the seq type, so in a copy an edge to D2's
    //deletion could land on D1's.
    unordered_map<Rec_Event_name, shared_ptr<Rec_Event>> copy_named;
    for (const shared_ptr<Rec_Event> &event : this->events) {
        copy_named.emplace(event->get_name(), event);
    }
    for (const auto &[name, other_adjacency] : other.edges) {
        Adjacency_list adjacency_list;
        for (const shared_ptr<Rec_Event> &child : other_adjacency.children) {
            adjacency_list.children.push_back(copy_named.at(child->get_name()));
        }
        for (const shared_ptr<Rec_Event> &parent : other_adjacency.parents) {
            adjacency_list.parents.push_back(copy_named.at(parent->get_name()));
        }
        this->edges.emplace(name, adjacency_list);
    }
    this->error_rate = other.error_rate->copy();
    //The registry holds no owning pointers, so a plain copy is a deep copy.
    //It must be carried over: inference deep-copies Model_Parms once per OpenMP
    //thread, and every consumer of the seq_type ordering reads it from that copy.
    this->seq_type_registry = other.seq_type_registry;

    //The events above were rebuilt through copy(), which does not carry resolved ids or
    //adjacency. Re-running finalize() restores them; the registry is already frozen, so this
    //resolves and does not rebuild.
    this->finalize();
}

/*
 * Provide a deep copy of the model parameters
 */
/*Model_Parms::Model_Parms(const Model_Parms& other){
	this->error_rate = other.error_rate->copy();
	for(unordered_map<Rec_Event_name , Adjacency_list>::const_iterator iter=other.edges.begin() ; iter != other.edges.end();iter++){
		this->edges.emplace((*iter).first,*iter.)
	 }

	for(list<Rec_Event*>::const_iterator iter = other.events.begin() ; iter != other.events.end() ; iter++){
		this->events.emplace_back((*iter)->copy());
	}
}*/

/*
 * Avoid memory leaks by deleting all events
 */
Model_Parms::~Model_Parms()
{
    /*	for(list<Rec_Event*>::iterator iter = events.begin() ; iter != events.end() ; iter++){
		delete (*iter);
	}
	delete error_rate;*/
}

/*
 * Add an event to the event list, and initialize its connectivity
 */
bool Model_Parms::add_event(shared_ptr<Rec_Event> event_point)
{
    needs_finalize = true;
    //Adding an event after finalize() would leave it with an unresolved SeqTypeId, and
    //could require a seq_type the frozen registry no longer accepts.
    if (seq_type_registry.is_frozen()) {
        throw std::logic_error("Model_Parms::add_event(): the model has been finalized");
    }
    //The nickname is the event's name: the key of the edges and the marginals' maps, and what
    //the files and outputs name the event by (plan T2). An event without one, or with one the
    //model already has, would share every entry with another.
    if (event_point->get_nickname().empty()) {
        throw runtime_error("Model_Parms::add_event(): the " + (string() + event_point->get_type())
                            + " event on \"" + event_point->get_seq_type() + "\" has no nickname");
    }
    if (this->edges.count(event_point->get_name()) != 0) {
        throw runtime_error("Model_Parms::add_event(): the model already has an event named \""
                            + event_point->get_name() + "\"");
    }
    event_point->set_event_identifier(this->events.size());
    this->events.emplace_back(event_point);
    this->edges.emplace(event_point->get_name(), Adjacency_list());
    return 1;
}
bool Model_Parms::add_event(Rec_Event *event_point)
{
    return this->add_event(shared_ptr<Rec_Event>(event_point, null_delete<Rec_Event>()));
}

/*
 * Get the list of parents of the given Rec_event
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_parents(Rec_Event *event) const
{
    return this->get_parents(shared_ptr<Rec_Event>(event, null_delete<Rec_Event>()));
}
list<shared_ptr<Rec_Event>> Model_Parms::get_parents(shared_ptr<Rec_Event> event) const
{
    return this->get_parents(event->get_name());
}
list<shared_ptr<Rec_Event>> Model_Parms::get_parents(Rec_Event_name event_name) const
{
    if (edges.count(event_name) <= 0) {
        throw runtime_error("Model_Parms::get_parents(): event \"" + event_name + "\" does not exist in \"this\".");
    }
    return edges.at(event_name).parents;
}

/*
 * Get the list of children of the given Rec_event
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_children(Rec_Event *event) const
{
    return this->get_children(shared_ptr<Rec_Event>(event, null_delete<Rec_Event>()));
}
list<shared_ptr<Rec_Event>> Model_Parms::get_children(shared_ptr<Rec_Event> event) const
{
    return this->get_children(event->get_name());
}
list<shared_ptr<Rec_Event>> Model_Parms::get_children(Rec_Event_name event_name) const
{
    if (edges.count(event_name) <= 0) {
        throw runtime_error("Model_Parms::get_children(): event \"" + event_name + "\" does not exist in \"this\".");
    }
    return edges.at(event_name).children;
}

/**
 * \overload list<shared_ptr<Rec_Event>> Model_Parms::get_ancestors(Rec_Event_name event_name) const
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_ancestors(Rec_Event *event) const
{
    return this->get_ancestors(shared_ptr<Rec_Event>(event, null_delete<Rec_Event>()));
}
/**
 * \overload list<shared_ptr<Rec_Event>> Model_Parms::get_ancestors(Rec_Event_name event_name) const
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_ancestors(shared_ptr<Rec_Event> event) const
{
    return this->get_ancestors(event->get_name());
}
/**
 * \brief Get all ancestors of the supplied event
 * \author Q.Marcou
 * \version 1.0
 * \param[in] event_name The event whom we seek the ancestors
 * \return The list of the event ancestors
 *
 *  Get all ancestors of the supplied event. Ancestors are all events higher than the considered event
 *  in the genealogy and with link of any degree to it
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_ancestors(Rec_Event_name event_name) const
{
    set<shared_ptr<Rec_Event>> ancestor_set;
    list<shared_ptr<Rec_Event>> exploratory_list;
    list<shared_ptr<Rec_Event>> final_list;

    //Get the considered event parent
    if (edges.count(event_name) <= 0) {
        throw runtime_error("Model_Parms::get_parents(): event \"" + event_name + "\" does not exist in \"this\".");
    }
    exploratory_list = edges.at(event_name).parents;

    //Now explore all ancestors
    while (not exploratory_list.empty()) {
        //If the ancestor is not already in the list
        if (ancestor_set.count(exploratory_list.front()) <= 0) {
            final_list.emplace_back(exploratory_list.front());
            ancestor_set.emplace(exploratory_list.front());
            const list<shared_ptr<Rec_Event>> &parents = edges.at(exploratory_list.front()->get_name()).parents;
            ;
            exploratory_list.insert(exploratory_list.end(), parents.begin(), parents.end());
        }
        exploratory_list.pop_front();
    }

    return final_list;
}

/*
 * Add a directed edge between two events
 */
bool Model_Parms::add_edge(Rec_Event *parent_point, Rec_Event *child_point)
{
    shared_ptr<Rec_Event> parent_smart_p(parent_point, null_delete<Rec_Event>());
    shared_ptr<Rec_Event> child_smart_p(child_point, null_delete<Rec_Event>());
    return this->add_edge(parent_smart_p, child_smart_p);
}
/**
 * \bug this is not checking whether the exact same edge already exists
 */
bool Model_Parms::add_edge(shared_ptr<Rec_Event> parent_point, shared_ptr<Rec_Event> child_point)
{
    // check whether these events exist
    if (edges.count(parent_point->get_name()) <= 0) {
        throw runtime_error("Model_Parms::add_edge(): event \"" + parent_point->get_name()
                            + "\" does not exist in \"this\".");
    }
    if (edges.count(child_point->get_name()) <= 0) {
        throw runtime_error("Model_Parms::add_edge(): event \"" + child_point->get_name()
                            + "\" does not exist in \"this\".");
    }
    //Check whether it is creating a cycle
    //First check the silly self-loop cycle
    if (parent_point->get_name() == child_point->get_name()) {
        throw runtime_error("Trying to create an edge between " + parent_point->get_name() + "and itself");
    }

    //Else check if the new children is an ancestor of the new parent (this would close a cycle)
    list<shared_ptr<Rec_Event>> parent_ancestors = this->get_ancestors(parent_point);
    for (const shared_ptr<Rec_Event> event_ptr : parent_ancestors) {
        if (child_point->get_name() == event_ptr->get_name()) {
            throw runtime_error(child_point->get_name() + " is an ancestor of " + parent_point->get_name()
                                + " adding an edge would create a cycle in the graph, in Model_Parms::add_edge");
        }
    }

    this->edges.at(parent_point->get_name()).children.push_back(child_point);
    this->edges.at(child_point->get_name()).parents.push_back(parent_point);
    return 1;
}

bool Model_Parms::add_edge(Rec_Event_name parent_name, Rec_Event_name child_name)
{
    shared_ptr<Rec_Event> parent_smart_p = this->get_event_pointer(parent_name);
    shared_ptr<Rec_Event> child_smart_p = this->get_event_pointer(child_name);
    return this->add_edge(parent_smart_p, child_smart_p);
}

/*
 * Remove the corresponding edge if it exists
 */
bool Model_Parms::remove_edge(Rec_Event *parent_point, Rec_Event *child_point)
{
    shared_ptr<Rec_Event> parent_smart_p(parent_point, null_delete<Rec_Event>());
    shared_ptr<Rec_Event> child_smart_p(child_point, null_delete<Rec_Event>());
    return this->remove_edge(parent_smart_p, child_smart_p);
}
bool Model_Parms::remove_edge(shared_ptr<Rec_Event> parent_point, shared_ptr<Rec_Event> child_point)
{
    if (this->has_edge(parent_point, child_point)) {
        list<shared_ptr<Rec_Event>>::iterator iter;

        //Remove the child from the parent's children list
        list<shared_ptr<Rec_Event>> &children_list = this->edges.at(parent_point->get_name()).children;
        iter = std::find(children_list.begin(), children_list.end(), child_point);
        if (iter != children_list.end()) {
            children_list.erase(iter);
        }

        //Remove the parent from the child's parent list
        list<shared_ptr<Rec_Event>> &parents_list = this->edges.at(child_point->get_name()).parents;
        iter = std::find(parents_list.begin(), parents_list.end(), parent_point);
        if (iter != parents_list.end()) {
            parents_list.erase(iter);
        }
        return 1;
    } else {
        throw runtime_error("Model_Parms::remove_edge(): edge between \"" + parent_point->get_name() + "\" and \""
                            + child_point->get_name() + "\" does not exist.");
    }
}

bool Model_Parms::remove_edge(Rec_Event_name parent_name, Rec_Event_name child_name)
{
    shared_ptr<Rec_Event> parent_smart_p = this->get_event_pointer(parent_name);
    shared_ptr<Rec_Event> child_smart_p = this->get_event_pointer(child_name);
    return this->remove_edge(parent_smart_p, child_smart_p);
}

/*
 * Inverse the direction of an edge between two events
 * Detects automatically the edge direction invert it
 */
void Model_Parms::invert_edge(Rec_Event *ev1_point, Rec_Event *ev2_point)
{
    shared_ptr<Rec_Event> ev1_smart_p(ev1_point, null_delete<Rec_Event>());
    shared_ptr<Rec_Event> ev2_smart_p(ev2_point, null_delete<Rec_Event>());
    this->invert_edge(ev1_smart_p, ev2_smart_p);
}

void Model_Parms::invert_edge(shared_ptr<Rec_Event> ev1_point, shared_ptr<Rec_Event> ev2_point)
{

    //First check if ev2 is a child of ev1
    if (this->has_edge(ev1_point, ev2_point)) {
        //Invert the edge
        this->remove_edge(ev1_point, ev2_point);
        try {
            this->add_edge(ev2_point, ev1_point);
        } catch (exception &e) {
            //add_edge checks for cycle creation since 13/04/2017
            cerr << "Exception caught trying to invert an edge between " + ev1_point->get_name() + " and "
                            + ev2_point->get_name()
                            + " in Model_Parms::invert_edge , this is most likely creating a cycle, throwing exception "
                              "now..."
                 << endl;
            throw e;
        }

    }
    //Otherwise check if ev2 is a parent of ev1
    else if (this->has_edge(ev2_point, ev1_point)) {
        //Invert the edge
        this->remove_edge(ev2_point, ev1_point);
        try {
            this->add_edge(ev1_point, ev2_point);
        } catch (exception &e) {
            //add_edge checks for cycle creation since 13/04/2017
            cerr << "Exception caught trying to invert an edge between " + ev1_point->get_name() + " and "
                            + ev2_point->get_name()
                            + " in Model_Parms::invert_edge , this is most likely creating a cycle, throwing exception "
                              "now..."
                 << endl;
            throw e;
        }
    } else {
        //else: the edge do not exist and throw an exception
        throw runtime_error("In Model_Parms::invert_edge(): no edge exist between \"" + ev1_point->get_name()
                            + "\" and \"" + ev1_point->get_name() + "\" events in any direction");
    }
}

void Model_Parms::invert_edge(Rec_Event_name ev1_name, Rec_Event_name ev2_name)
{
    shared_ptr<Rec_Event> ev1_smart_p = this->get_event_pointer(ev1_name);
    shared_ptr<Rec_Event> ev2_smart_p = this->get_event_pointer(ev2_name);
    return this->invert_edge(ev1_smart_p, ev2_smart_p);
}

/*
 * Check if the model parms contain a given directed edge
 *
 * Note: assuming the edge has been correctly constructed it will test the existence in only one of the adjacency list
 */
bool Model_Parms::has_edge(Rec_Event *parent_point, Rec_Event *child_point) const
{
    shared_ptr<Rec_Event> parent_smart_p(parent_point, null_delete<Rec_Event>());
    shared_ptr<Rec_Event> child_smart_p(child_point, null_delete<Rec_Event>());
    return this->has_edge(parent_smart_p, child_smart_p);
}

bool Model_Parms::has_edge(shared_ptr<Rec_Event> parent_point, shared_ptr<Rec_Event> child_point) const
{
    return this->has_edge(parent_point->get_name(), child_point->get_name());
}

bool Model_Parms::has_edge(Rec_Event_name parent_name, Rec_Event_name child_name) const
{
    bool other_event_found = false;
    const Adjacency_list &adjacency_list = this->edges.at(parent_name);
    for (shared_ptr<Rec_Event> ev_ptr : adjacency_list.children) {
        if (ev_ptr->get_name() == child_name) {
            other_event_found = true;
        }
    }
    return other_event_found;
}

/*
 * This method returns all the roots of the tree (events with no parents).
 */
list<shared_ptr<Rec_Event>> Model_Parms::get_roots() const
{
    list<shared_ptr<Rec_Event>> root_list = list<shared_ptr<Rec_Event>>(); //FIXME nonsense new
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = this->events.begin(); iter != this->events.end(); ++iter) {
        if (this->edges.at((*iter)->get_name()).parents.empty()) {
            root_list.push_back(*iter);
        }
    }
    // root events are sorted for easier manipulation
    root_list.sort(Event_comparator());
    return root_list;
}

/*
 * This method output the queue (order in which events are processed) according to the structure of the model
 * This queue is used both for the order of iteration and the design of the marginal array.
 */
queue<shared_ptr<Rec_Event>> Model_Parms::get_model_queue() const
{

    list<shared_ptr<Rec_Event>> model_roots = this->get_roots();
    list<shared_ptr<Rec_Event>> events_copy_list = this->events;

    queue<shared_ptr<Rec_Event>> model_queue = queue<shared_ptr<Rec_Event>>();

    //if all events are root they are independent and thus the queue is sorted only by priority
    if (events_copy_list.size() == model_roots.size()) {
        //root list is already sorted by priority
        for (list<shared_ptr<Rec_Event>>::const_iterator iter = model_roots.begin(); iter != model_roots.end();
             ++iter) {
            model_queue.push(*iter);
        }
        return model_queue;
    }

    //The root with highest priority is the first event
    events_copy_list.remove(*(model_roots.begin()));
    model_queue.push(*(model_roots.begin()));

    //Keep track of the events already added to the queue
    unordered_map<Rec_Event_name, shared_ptr<Rec_Event>> *processed_events_point =
            new unordered_map<Rec_Event_name, shared_ptr<Rec_Event>>;
    (*processed_events_point).insert(make_pair((*model_roots.begin())->get_name(), *(model_roots.begin())));
    model_roots.pop_front();

    events_copy_list.sort(Event_comparator());
    list<shared_ptr<Rec_Event>>::iterator iter = events_copy_list.begin();

    while (!events_copy_list.empty()) {

        //Before an event is processed his parents must have been processed, this overrides the priority!

        list<shared_ptr<Rec_Event>> event_parents = this->edges.at((*iter)->get_name()).parents;
        bool parents_processed = 1;
        for (list<shared_ptr<Rec_Event>>::const_iterator jiter = event_parents.begin(); jiter != event_parents.end();
             ++jiter) {
            if ((*processed_events_point).count((*jiter)->get_name()) == 0) {
                parents_processed = 0;
            }
        }
        //If all parents have not been processed, tries the next highest priority event
        if (!parents_processed) {
            ++iter;
            continue;
        }

        model_queue.push(*iter);
        (*processed_events_point).insert(make_pair((*iter)->get_name(), (*iter)));
        events_copy_list.erase(iter);
        events_copy_list.sort(Event_comparator()); //TODO again check this sort function
        iter = events_copy_list.begin();
    }
    delete processed_events_point;
    return model_queue;
}

shared_ptr<Rec_Event> Model_Parms::get_event_pointer(const string &event_str, bool by_nickname) const
{
    //TODO find something better, might be a bottleneck
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = this->events.begin(); iter != this->events.end(); ++iter) {
        if (by_nickname) {
            if ((*iter)->get_nickname() == event_str) {
                return (*iter);
            }
        } else {
            if ((*iter)->get_name() == event_str) {
                return (*iter);
            }
        }
    }
    if (by_nickname) {
        throw runtime_error("Event pointer not found in Model_Parms::get_event_pointer for nickname:" + event_str);
    } else {
        throw runtime_error("Event pointer not found in Model_Parms::get_event_pointer for name:" + event_str);
    }
}

shared_ptr<Rec_Event> Model_Parms::get_event_pointer(const Rec_Event_name &event_name) const
{
    return this->get_event_pointer(event_name, false);
}

Events_map Model_Parms::get_events_map() const
{
    //Handing out a model whose events do not know their seq_type id or their neighbours is
    //the failure this flag exists to catch: it shows up far away, as an event addressing the
    //wrong segment. Debug-only -- the check is one bool, the fix is to call finalize().
    assert(!needs_finalize && "Model_Parms::finalize() must run after the last add_event()");
    Events_map events_map;
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = this->events.begin(); iter != this->events.end(); ++iter) {
        // DinucMarkov is keyed with Undefined_side: the seq_type already uniquely identifies
        // which junction it belongs to, and event_side carries direction info only (not identity).
        Seq_side map_side = ((*iter)->get_type() == Dinuclmarkov_t) ? Undefined_side : (*iter)->get_side();
        events_map.emplace(make_tuple((*iter)->get_type(), (*iter)->get_seq_type(), map_side), (*iter));
    }
    return events_map;
}

Events_map Model_Parms::get_events_map()
{
    //Handing out a model whose events do not know their seq_type id or their neighbours is
    //the failure this flag exists to catch: it shows up far away, as an event addressing the
    //wrong segment. Debug-only -- the check is one bool, the fix is to call finalize().
    assert(!needs_finalize && "Model_Parms::finalize() must run after the last add_event()");
    Events_map events_map;
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = this->events.begin(); iter != this->events.end(); ++iter) {
        Seq_side map_side = ((*iter)->get_type() == Dinuclmarkov_t) ? Undefined_side : (*iter)->get_side();
        events_map.emplace(make_tuple((*iter)->get_type(), (*iter)->get_seq_type(), map_side), (*iter));
    }
    return events_map;
}

void Model_Parms::finalize()
{
    //A name no longer stops two events of one type on one segment and side, and nothing in the
    //model forbids them. The legacy engine does: get_events_map() holds one event per such key,
    //and would drop the second without a word.
    map<tuple<Event_type, Seq_type_String, Seq_side>, shared_ptr<Rec_Event>> at_place;
    for (const shared_ptr<Rec_Event> &event : events) {
        const Seq_side map_side = (event->get_type() == Dinuclmarkov_t) ? Undefined_side : event->get_side();
        const auto [placed, inserted] =
                at_place.emplace(make_tuple(event->get_type(), event->get_seq_type(), map_side), event);
        if (not inserted) {
            throw runtime_error("Model_Parms::finalize(): \"" + placed->second->get_nickname() + "\" and \""
                                + event->get_nickname() + "\" are both " + (string() + event->get_type())
                                + " events on \"" + event->get_seq_type() + "\" at " + (string() + map_side)
                                + ": the legacy engine holds one event per type, segment and side");
        }
    }

    if (!seq_type_registry.is_frozen()) {
        //Pin the legacy ids first, so that an inferred ordering does not renumber them.
        seq_type_registry.register_legacy_seq_types();

        //Legacy files carry no @Seq_type_order, and programmatically built models set none;
        //infer it from whether the event graph has a D gene choice.
        if (seq_type_registry.empty()) {
            bool has_d = false;
            for (const auto &ev : events) {
                if (ev->get_type() == GeneChoice_t && ev->get_class() == D_gene) {
                    has_d = true;
                    break;
                }
            }
            if (has_d) {
                seq_type_registry.set_ordered_types(
                        {"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"});
            } else {
                seq_type_registry.set_ordered_types({"V_gene_seq", "VJ_ins_seq", "J_gene_seq"});
            }
        }
        seq_type_registry.freeze();
    }
    //Resolve each event's seq_type name to its runtime handle. Events whose seq_type is not
    //in the registry keep kNoSeqType: legacy models that never declared an ordering, and
    //error-rate-only parameter sets, are both legitimate in that state.
    for (const auto &event : events) {
        const Seq_type_String &name = event->get_seq_type();
        if (!name.empty() && seq_type_registry.contains(name)) {
            event->set_seq_type_id(seq_type_registry.id(name));
        }
    }

    //Second pass, once every id is known: tell each event what sits next to it. This is the
    //only place the ordering is turned into adjacency -- events are handed the answer rather
    //than given a registry to ask, so "who is my neighbour" has one owner and one moment.
    //Separate from the loop above because it reads ids the first pass may not have set yet.
    for (const auto &event : events) {
        const SeqTypeId id = event->get_seq_type_id();
        if (id == kNoSeqType) {
            event->set_adjacent_segments(kNoSeqType, kNoSeqType);
            continue;
        }
        event->set_adjacent_segments(seq_type_registry.left_neighbor(id),
                                     seq_type_registry.right_neighbor(id));
    }

    needs_finalize = false;
}

namespace {

/// Whether a v1 event line reads back as \a event. The reader derives a seq type from the line's
/// gene class: a gene choice's own class, the class legacy_gene_class_of() gives the others. An
/// event built in code with no seq type yet gets one from its class, so it reads back too.
bool v1_line_names_its_segment(const Rec_Event &event)
{
    if (event.get_seq_type().empty())
        return event.get_type() == GeneChoice_t;
    const Gene_class_legacy v1_class = legacy_gene_class_of(event.get_type(), event.get_seq_type());
    if (v1_class == Undefined_gene_legacy)
        return false;
    return event.get_type() != GeneChoice_t or gene_class_legacy_to_new(v1_class) == event.get_class();
}

}  // namespace

/**
 * Whether the model needs JSON, the v1 format being unable to write it and read it back as itself:
 * - an event whose v1 line would read back as another segment, or two events sharing a v1 name;
 * - a seq_type outside the six legacy ones;
 * - an ordering other than VDJ or VJ.
 */
bool Model_Parms::requires_extended_format() const
{
    static const std::unordered_set<Seq_type_String> legacy_seq_types = {
        "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq", "VJ_ins_seq"
    };
    static const std::vector<Seq_type_String> vdj_order = {
        "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"
    };
    static const std::vector<Seq_type_String> vj_order = {
        "V_gene_seq", "VJ_ins_seq", "J_gene_seq"
    };

    //A v1 file names each segment by an event line's gene class, and each event in @Edges by
    //its v1 name. So every event's line has to read back as that event, and no two events may
    //share a v1 name: an insertion on a gene segment has no v1 class, and two deletions on two
    //junctions would get one name.
    unordered_set<Rec_Event_name> v1_names;
    for (const auto &event : events) {
        if (not v1_line_names_its_segment(*event) or not v1_names.insert(event->get_v1_legacy_name()).second)
            return true;
    }

    const auto &ordered = seq_type_registry.get_ordered_types();
    if (ordered.empty())
        return false;
    for (const auto &st : ordered) {
        if (legacy_seq_types.find(st) == legacy_seq_types.end())
            return true;
    }
    return (ordered != vdj_order && ordered != vj_order);
}

void Model_Parms::write_model_parms(string filename)
{
    if (requires_extended_format())
        write_model_parms_json(filename);
    else
        write_model_parms_legacy(filename);
}

void Model_Parms::write_model_parms_legacy(string filename)
{
    ofstream outfile(filename);
    outfile << "@Event_list" << endl;
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = events.begin(); iter != events.end(); ++iter) {
        (*iter)->write2txt_legacy(outfile);
    }
    outfile << "@Edges" << endl;
    for (list<shared_ptr<Rec_Event>>::const_iterator iter = events.begin(); iter != events.end(); ++iter) {
        list<shared_ptr<Rec_Event>> children = edges[(*iter)->get_name()].children;
        for (list<shared_ptr<Rec_Event>>::const_iterator jiter = children.begin(); jiter != children.end(); ++jiter) {
            outfile << "%" << (*iter)->get_v1_legacy_name() << ";" << (*jiter)->get_v1_legacy_name() << endl;
        }
    }
    outfile << "@ErrorRate" << endl;
    error_rate->write2txt(outfile);
}

void Model_Parms::write_model_parms_json(string filename)
{
    ofstream outfile(filename);
    outfile << model_parms_to_json(*this).dump(2) << endl;
}

void Model_Parms::set_seq_type_order(const vector<Seq_type_String> &order)
{
    //As read_model_parms() does: the legacy ids pinned first (see there), then the order.
    seq_type_registry.register_legacy_seq_types();
    if (not order.empty()) {
        seq_type_registry.set_ordered_types(order);
    }
}

void Model_Parms::read_model_parms(string filename)
{
    ifstream infile(filename);
    if (!infile) {
        //Throw exception
        throw runtime_error("File not found : \"" + filename + "\"");
    }

    //A model the v1 format cannot hold is a JSON document (plan T2), told apart by its content.
    infile >> ws;
    if (infile.peek() == '{') {
        nlohmann::json doc;
        try {
            doc = nlohmann::json::parse(infile);
        } catch (const nlohmann::json::exception &e) {
            throw runtime_error("Model file \"" + filename + "\" is not valid JSON: " + e.what());
        }
        model_parms_from_json(doc, *this);
        return;
    }

    //Pin the six legacy seq_types to ids 0-5, matching the Seq_type enum, before any
    //ordering is applied. Code still keyed by the enum then addresses the same slots as
    //code keyed by SeqTypeId. See SeqTypeRegistry::register_legacy_seq_types().
    seq_type_registry.register_legacy_seq_types();

    string line_str;
    getline(infile, line_str);

    if (line_str == string("@Version")) {
        throw runtime_error("Model file \"" + filename + "\" is in the v2 text format, which was retired before "
                            "any release: a model the v1 format cannot hold is written as JSON");
    }

    // Now parse @Event_list section
    if (line_str == string("@Event_list")) {
        getline(infile, line_str);
        while (line_str[0] == '#') {
            size_t semicolon_index = line_str.find(";", 0);
            string event = line_str.substr(1, semicolon_index - 1);
            size_t next_semicolon_index = line_str.find(";", semicolon_index + 1);
            string event_class_str = line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1));
            Gene_class_legacy event_class;
            try {
                event_class = str2GeneClass(event_class_str);
            } catch (exception &e) {
                throw runtime_error("Unknown Gene_class_legacy\"" + event_class_str + "\" in model file: \"" + filename
                                    + "\"");
            }
            
            // event_type;gene_class;seq_side;priority;nickname
            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            string event_side_str = line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1));
            Seq_side event_side;
            try {
                event_side = str2SeqSide(event_side_str);
            } catch (exception &e) {
                throw runtime_error("Unknown Seq_side\"" + event_side_str + "\" in file: \"" + filename + "\"");
            }

            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            //The nickname is the event's name (plan T2). Early model files could leave it out; none
            //that IGoR ships does, and adding one is cheap.
            if (next_semicolon_index == string::npos
                or next_semicolon_index + 1 == line_str.size()) {
                throw runtime_error("Event line \"" + line_str + "\" in model file \"" + filename
                                    + "\" has no nickname: every event needs one, unique in its model");
            }
            const int priority = stoi(line_str.substr(semicolon_index + 1, next_semicolon_index - 1));
            const string nickname = line_str.substr(next_semicolon_index + 1, string::npos);

            cerr << event << " read" << endl;
            
            // The seq type comes from the gene class and the event type; a gene choice is also
            // built with its class. The other three events take only the seq type: a deletion
            // derives its gene class from it.
            const Gene_class effective_gene_class = gene_class_legacy_to_new(event_class);
            const Event_type evt_type = (event == "GeneChoice") ? GeneChoice_t
                                        : (event == "Deletion") ? Deletion_t
                                        : (event == "Insertion") ? Insertion_t
                                                                 : Dinuclmarkov_t;
            const Seq_type_String event_seq_type = legacy_gene_class_to_seq_type(event_class);
            // The two lines the event constructors used to refuse while turning the name back
            // into a Seq_type, refused here where the gene class still means something: a
            // class naming no segment, and a deletion on a junction. A gene choice was never
            // refused for either, and is not now.
            const bool junction_class =
                    (event_class == VD_genes or event_class == DJ_genes or event_class == VJ_genes);
            if ((evt_type != GeneChoice_t and event_seq_type == "Undefined_seq")
                or (evt_type == Deletion_t and junction_class)) {
                throw runtime_error("Gene class \"" + to_string(event_class) + "\" names no segment a "
                                    + event + " event can act on, in model file \"" + filename + "\"");
            }

            if (event == string("Insertion")) {
                unordered_map<string, Event_realization> event_realizations =
                        unordered_map<string, Event_realization>();
                getline(infile, line_str);
                while (line_str[0] == '%') {
                    semicolon_index = line_str.find(";", 0);
                    int value_int = stoi(line_str.substr(1, semicolon_index));
                    int index = stoi(line_str.substr(semicolon_index + 1, string::npos));
                    event_realizations.emplace(pair<string, Event_realization>(
                            to_string(value_int),
                            Event_realization(to_string(value_int), value_int, "", Int_Str(), index)));
                    getline(infile, line_str);
                }
                shared_ptr<Insertion> new_event_p =
                        shared_ptr<Insertion>(new Insertion(event_seq_type, event_realizations));
                new_event_p->set_priority(priority);
                new_event_p->set_nickname(nickname);
                this->add_event(new_event_p);
            } else if (event == string("Deletion")) {
                unordered_map<string, Event_realization> event_realizations =
                        unordered_map<string, Event_realization>();
                getline(infile, line_str);
                while (line_str[0] == '%') {
                    semicolon_index = line_str.find(";", 0);
                    int value_int = stoi(line_str.substr(1, semicolon_index));
                    int index = stoi(line_str.substr(semicolon_index + 1, string::npos));
                    event_realizations.emplace(pair<string, Event_realization>(
                            to_string(value_int),
                            Event_realization(to_string(value_int), value_int, "", Int_Str(), index)));
                    getline(infile, line_str);
                }
                shared_ptr<Deletion> new_event_p =
                        shared_ptr<Deletion>(new Deletion(event_seq_type, event_side, event_realizations));
                new_event_p->set_priority(priority);
                new_event_p->set_nickname(nickname);
                this->add_event(new_event_p);
            } else if (event == string("GeneChoice")) {
                unordered_map<string, Event_realization> event_realizations =
                        unordered_map<string, Event_realization>(); //FIXME nonsense new
                getline(infile, line_str);
                while (line_str[0] == '%') {
                    semicolon_index = line_str.find(";", 0);
                    string name = line_str.substr(1, semicolon_index - 1);
                    next_semicolon_index = line_str.find(";", semicolon_index + 1);
                    string value_str =
                            line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1));
                    string test = line_str.substr(next_semicolon_index + 1, string::npos);
                    int index = stoi(line_str.substr(next_semicolon_index + 1, string::npos));
                    event_realizations.emplace(pair<string, Event_realization>(
                            name, Event_realization(name, INT16_MAX, value_str, nt2int(value_str), index)));
                    getline(infile, line_str);
                }
                // For GeneChoice, use the effective_gene_class (Gene_class) derived from the legacy event_class
                shared_ptr<Gene_choice> new_event_p = shared_ptr<Gene_choice>(new Gene_choice(
                        effective_gene_class,
                        event_realizations));
                new_event_p->set_priority(priority);
                new_event_p->set_nickname(nickname);
                new_event_p->set_seq_type(event_seq_type);
                this->add_event(new_event_p);
            } else if (event == string("DinucMarkov")) {
                shared_ptr<Dinucl_markov> new_event_p = shared_ptr<Dinucl_markov>(new Dinucl_markov(event_seq_type));
                new_event_p->set_priority(priority);
                new_event_p->set_nickname(nickname);
                // The direction the chain runs in: a v1 line says Undefined_side, and the gene class
                // gives it.
                Seq_side dinucl_side = Undefined_side;
                if (event_class == VD_genes || event_class == VJ_genes) {
                    dinucl_side = Three_prime;
                } else if (event_class == DJ_genes) {
                    dinucl_side = Five_prime;
                }
                new_event_p->set_event_side(dinucl_side);
                this->add_event(new_event_p);
                getline(infile, line_str);
                while (line_str[0] == '%') {
                    getline(infile, line_str);
                }
            } else {
                throw runtime_error(event + " event is not implemented (thrown by Model_Parms::read_model_parms");
            }
        }
    } else {
        throw runtime_error("Unknown format for model_parms file");
    }

    if (line_str == string("@Edges")) {
        //An edge names its events by their v1 name. The token finds the event, and the edge is
        //keyed by the event's name, its nickname. A token naming two events is refused: the v1
        //name tells events apart only where a v1 gene class does.
        const auto event_named = [&](const string &token) {
            shared_ptr<Rec_Event> found;
            for (const shared_ptr<Rec_Event> &event : this->events) {
                if (event->get_v1_legacy_name() != token) {
                    continue;
                }
                if (found) {
                    throw runtime_error("Edge token \"" + token + "\" in model file \"" + filename
                                        + "\" names two events, \"" + found->get_nickname() + "\" and \""
                                        + event->get_nickname() + "\"");
                }
                found = event;
            }
            if (not found) {
                throw runtime_error("Edge token \"" + token + "\" in model file \"" + filename
                                    + "\" names no event");
            }
            return found;
        };
        getline(infile, line_str);
        while (line_str[0] == '%') {
            size_t semicolon_index = line_str.find(";", 0);
            const shared_ptr<Rec_Event> parent = event_named(line_str.substr(1, semicolon_index - 1));
            const shared_ptr<Rec_Event> child = event_named(line_str.substr(semicolon_index + 1, string::npos));
            edges.at(parent->get_name()).children.emplace_back(child);
            edges.at(child->get_name()).parents.emplace_back(parent);
            getline(infile, line_str);
        }
    } else {
        throw runtime_error("Unknown format for model file: " + filename);
    }
    if (line_str == string("@ErrorRate")) {
        getline(infile, line_str);
        size_t semicolon_index = line_str.find(";");
        string errrate = line_str.substr(0, semicolon_index);
        if (errrate == string("#SingleErrorRate")) {
            getline(infile, line_str);
            shared_ptr<Single_error_rate> err_rate_p =
                    shared_ptr<Single_error_rate>(new Single_error_rate(stod(line_str)));
            this->error_rate = err_rate_p;
        } else if (errrate == string("#Hypermutationglobalerrorrate")) {

            //Get the general attributes for the hypermutation model (size,learn_on,apply_on)
            size_t next_semicolon_index = line_str.find(";", semicolon_index + 1);
            size_t mutation_Nmer_size =
                    stoi(line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            Gene_class_legacy learn_on;
            try {
                learn_on = str2GeneClass(
                        line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            } catch (exception &e) {
                throw runtime_error("Unknown Gene_class_legacy\""
                                    + line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))
                                    + "\" for Hypermutationglobalerrorrate in model file: " + filename);
            }
            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            Gene_class_legacy apply_on;
            try {
                apply_on = str2GeneClass(
                        line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            } catch (exception &e) {
                throw runtime_error("Unknown Gene_class_legacy\""
                                    + line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))
                                    + "\" for Hypermutationglobalerrorrate in model file: " + filename);
            }

            //Get the global rate
            getline(infile, line_str);
            double R = stod(line_str);

            //Get the contributions
            getline(infile, line_str);
            semicolon_index = 0;
            next_semicolon_index = line_str.find(";");
            vector<double> ei_contributions;
            while (semicolon_index != string::npos) {
                if (semicolon_index == 0) {
                    ei_contributions.push_back(
                            stod(line_str.substr(semicolon_index, (next_semicolon_index - semicolon_index))));
                } else {
                    ei_contributions.push_back(
                            stod(line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))));
                }
                semicolon_index = next_semicolon_index;
                next_semicolon_index = line_str.find(";", semicolon_index + 1);
            }
            shared_ptr<Hypermutation_global_errorrate> err_rate_p = shared_ptr<Hypermutation_global_errorrate>(
                    new Hypermutation_global_errorrate(mutation_Nmer_size, learn_on, apply_on, R, ei_contributions));
            this->set_error_ratep(err_rate_p);
        } else if (errrate == string("#HypermutationfullNmererrorrate")) {
            //Get the general attributes for the hypermutation model (size,learn_on,apply_on)
            size_t next_semicolon_index = line_str.find(";", semicolon_index + 1);
            size_t mutation_Nmer_size =
                    stoi(line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            Gene_class_legacy learn_on;
            try {
                learn_on = str2GeneClass(
                        line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            } catch (exception &e) {
                throw runtime_error("Unknown Gene_class_legacy\""
                                    + line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))
                                    + "\" for Hypermutationglobalerrorrate in model file: " + filename);
            }
            semicolon_index = next_semicolon_index;
            next_semicolon_index = line_str.find(";", semicolon_index + 1);
            Gene_class_legacy apply_on;
            try {
                apply_on = str2GeneClass(
                        line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1)));
            } catch (exception &e) {
                throw runtime_error("Unknown Gene_class_legacy\""
                                    + line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))
                                    + "\" for Hypermutationglobalerrorrate in model file: " + filename);
            }

            //Get the mutation probas
            getline(infile, line_str);
            semicolon_index = 0;
            next_semicolon_index = line_str.find(";");
            vector<double> mutation_probas;
            while (semicolon_index != string::npos) {
                if (semicolon_index == 0) {
                    mutation_probas.push_back(
                            stod(line_str.substr(semicolon_index, (next_semicolon_index - semicolon_index))));
                } else {
                    mutation_probas.push_back(
                            stod(line_str.substr(semicolon_index + 1, (next_semicolon_index - semicolon_index - 1))));
                }
                semicolon_index = next_semicolon_index;
                next_semicolon_index = line_str.find(";", semicolon_index + 1);
            }
            shared_ptr<Hypermutation_full_Nmer_errorrate> err_rate_p = shared_ptr<Hypermutation_full_Nmer_errorrate>(
                    new Hypermutation_full_Nmer_errorrate(mutation_Nmer_size, learn_on, apply_on, mutation_probas));
            this->set_error_ratep(err_rate_p);
        } else {
            throw runtime_error("Unknown Error_rate type\"" + errrate + "\" in model file: " + filename);
        }

    } else {
        throw runtime_error("Unknown format for model file: " + filename);
    }

    //Seal the registry and resolve every event's SeqTypeId (plan decision D1b).
    this->finalize();
}

void Model_Parms::set_fixed_all_events(bool fix_bool_status)
{
    for (list<shared_ptr<Rec_Event>>::iterator iter = events.begin(); iter != events.end(); ++iter) {
        (*iter)->fix(fix_bool_status);
    }
}

} // namespace igor::model::legacy
