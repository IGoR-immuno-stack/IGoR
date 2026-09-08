/*
 * Rec_Event.h
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

#pragma once

#include <igor/Core/Errorrate.h>
#include <igor/Core/IntStr.h>
#include <igor/Core/Aligner.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/Utils.h>
#include <igorCoreExport.h>

// Context objects for refactored iterate()
#include <igor/Core/QuerySequenceContext.h>
#include <igor/Core/ModelContext.h>
#include <igor/Core/ScenarioContext.h>
#include <igor/Core/ExplorationContext.h>
#include <igor/Core/AccumulationContext.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <list>
#include <utility>
#include <forward_list>
#include <queue>
#include <random>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <memory>
#include <map>

class Counter;

//class Model_marginals; //forward declare model marginals to avoid circular inclusion

// value of event: struct: event identifier(name of Vgene), event value(sequence), event index(custom)
/**
 * \struct Event_realization Rec_Event.h
 * \brief Unit that stores an event realization name, value and index.
 * \author Q.Marcou
 * \version 1.0
 *
 *	Depending on the RecEvent type to which it belongs, the Event_realization must supply either a string (both std::string and IntStr) or an integer value.
 *	Integers values are e.g the number of deletions or insertions of Insertion or Deletion RecEvent
 *	String values are e.g realization of a GeneChoice Rec_Event, and stands for the gene sequence.
 *
 */
struct Event_realization
{
    std::string name;
    int value_int; //union? template? inheritance and reference? just use a virtual class containing two types of events:str and int
    std::string value_str;
    Int_Str value_str_int;
    int index; //Not defined by the user but at the creation of the event, not quite sure about the mutable

    Event_realization(std::string real_name, int val_int, std::string val_str, Int_Str val_str_int, int index_val)
        : name(real_name), value_int(val_int), value_str(val_str), value_str_int(val_str_int), index(index_val)
    {
    }
};

/**
 * \class Rec_Event Rec_Event.h
 * \brief Recombination event class (IGoR's graph nodes)
 * \author Q.Marcou
 * \version 1.0
 *
 * This class implements the recombination event object.
 * Rec_Events are the nodes in IGoR's Bayesian Network structure.
 * This is a purely abstract class and cannot be instanciated as is, only classes deriving from it and implementing the purely abstract methods can be.
 *
 * Rec_Events contain the different Event_realization associated to it in a hashmap.
 *
 * The RecEvents design is key to the way IGoR explore all possible scenarios (through the iterate method) and generate sequences (through the draw_random_realization)
 *
 */
/**
 * \brief Signed range by which an event shifts one end of a segment.
 *
 * A 3' end moves left as nucleotides are deleted, a 5' end moves right, so the sign is
 * carried here rather than being re-derived at every call site. `{0, 0}` means the event
 * cannot move that end at all.
 *
 * **Invariant: `min <= max`.** The pair is an ordered interval, not a (nearest, farthest)
 * pair -- an implementer flipping the sign of one end must reorder the two. Consumers rely
 * on this rather than re-sorting, so that a provider that gets it wrong surfaces as an error
 * instead of being silently normalised; JunctionGeometry::PendingModifierBounds::rebuild()
 * is where it is enforced.
 */
struct OffsetDelta {
    int min = 0;
    int max = 0;

    bool operator==(const OffsetDelta &) const = default;
};

/**
 * \brief Signed range of nucleotides an event contributes to a segment's length.
 *
 * Positive for an event that supplies sequence (a genomic template, an insertion), negative
 * for one that removes it (a deletion), zero for one that only fills placeholders.
 *
 * **Invariant: `min <= max`**, as for OffsetDelta.
 *
 * This is a *contribution*, not a length: several events compose onto one segment, and the
 * junction-length DP already sums them. It replaces `get_len_min()` / `get_len_max()`, whose
 * meaning differs per subclass -- signed delta on Deletion, a length on Insertion and
 * Gene_choice, never set on Dinucl_markov -- and whose value is accumulated by an
 * order-dependent `if / else if` over an unordered map. See
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md sections 2.1 and 7.4.
 */
struct LengthContribution {
    int min = 0;
    int max = 0;

    bool operator==(const LengthContribution &) const = default;
};

/// What an event does to a constructed sequence segment.
enum class SeqConstructionRole {
    None,     ///< does not touch this segment
    Creates,  ///< allocates it (possibly containing placeholders)
    Modifies, ///< truncates or extends an existing one
    Fills     ///< fills placeholder values in an existing one
};

/// What an event does to one end of a segment.
enum class OffsetRole {
    None,
    Creates, ///< sets the initial offset
    Modifies ///< shifts an existing one
};

class CORE_EXPORT Rec_Event
{
public:
    Rec_Event();
    Rec_Event(Gene_class, Seq_side);
    Rec_Event(Gene_class, Seq_side, std::unordered_map<std::string, Event_realization> &);
    virtual ~Rec_Event();
    virtual std::shared_ptr<Rec_Event> copy() = 0; //TODO make it const somehow
    virtual int size() const;
    //TODO get rid of deletion map and chosen gene map
    /**
     * @brief Context-based iterate() interface
     *
     * New signature using 5 context objects instead of 19 individual parameters.
     * Subclasses implement semantic iteration logic using these contexts.
     *
     * Contexts encapsulate:
     * - QuerySequenceContext: Input sequence and alignments
     * - ModelContext: Read-only model configuration
     * - ScenarioContext: Per-path mutable state
     * - ExplorationContext: Tree exploration policy
     * - AccumulationContext: Result accumulation
     *
     * This overload exists alongside the legacy signature during transition.
     */
    virtual void
    iterate(QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation) = 0;

    bool set_priority(int);

    //Accessors
    const Gene_class get_class() const { return event_class; };
    const Seq_side get_side() const { return event_side; };
    const std::unordered_map<std::string, Event_realization> get_realizations_map() const
    {
        return event_realizations;
    };
    const int get_priority() const { return priority; };
    const Rec_Event_name get_name() const { return name; };
    /// Returns the v2.0-format event name that includes seq_type between
    /// gene_class and seq_side.  Used in @Edges sections of v2.0 files.
    /// Falls back to get_name() when seq_type is empty.
    Rec_Event_name get_v2_name() const;
    /// Returns the legacy-format event name (without seq_type), even if the
    /// internal name has been updated to include seq_type.
    Rec_Event_name get_legacy_name() const;
    const std::string get_nickname() const { return nickname; };
    void set_nickname(std::string name) { nickname = name; }
    Event_type get_type() const { return this->type; }
    int get_len_max() const { return this->len_max; };
    int get_len_min() const { return this->len_min; };

    /**
     * \name Capability queries (Phase A, reduced -- task A0)
     *
     * Declarative properties, answered from the event's own instance data. Called during
     * model initialization, never inside the iterate hot loop, so the virtual dispatch is
     * free. Pure virtual rather than base members because several answers depend on
     * instance state a base constructor could not know: which seq_type an event targets,
     * which end it acts on, and its realization set.
     *
     * They exist so that the generic iterate() bodies (B5, B6, B11) can ask "how far can
     * this end still move" and "how long can what sits here still be" without knowing
     * which subclass answers. See docs/ITERATE_GENERIC_REWRITE_PLAN.md section 2.1.
     * @{
     */

    /// How far this event can still shift `(type_id, side)`. `{0, 0}` if it cannot.
    virtual OffsetDelta get_offset_delta_bounds(SeqTypeId type_id, Seq_side side) const = 0;

    /// How many nucleotides this event contributes to `type_id`. `{0, 0}` if none.
    virtual LengthContribution get_length_contribution(SeqTypeId type_id) const = 0;

    /// What this event does to the `type_id` segment.
    virtual SeqConstructionRole get_seq_construction_role(SeqTypeId type_id) const = 0;

    /// What this event does to one end of the `type_id` segment.
    virtual OffsetRole get_offset_role(SeqTypeId type_id, Seq_side side) const = 0;

    /**
     * The segments immediately 5' and 3' of this event's own, in the registry ordering.
     *
     * Set by Model_Parms::finalize() for every event, once per model, and `kNoSeqType` where
     * the ordering runs out or the event has no seq_type. **Topology is the model's fact, not
     * the event's**: events are told what is next to them rather than asking a registry, so
     * there is exactly one place that reads the ordering and one moment at which it is read.
     *
     * What an event *does* with its neighbours is its own business: an Insertion spans both,
     * a Dinucl_markov seeds from whichever side its Markov chain runs from, a Gene_choice
     * ignores them.
     */
    void set_adjacent_segments(SeqTypeId left, SeqTypeId right)
    {
        left_adjacent_id = left;
        right_adjacent_id = right;
    }
    SeqTypeId get_left_adjacent_id() const { return left_adjacent_id; }
    SeqTypeId get_right_adjacent_id() const { return right_adjacent_id; }

    /** @} */
    const Seq_type_String get_seq_type() const { return seq_type; };
    void set_seq_type(const Seq_type_String &st) { seq_type = st; }
    /// Runtime handle for seq_type, resolved against the model's frozen registry by
    /// Model_Parms::finalize(). kNoSeqType until then. Used in place of the name
    /// everywhere inside the scenario traversal, so the registry stays out of the hot path.
    SeqTypeId get_seq_type_id() const { return seq_type_id; }
    void set_seq_type_id(SeqTypeId id) { seq_type_id = id; }
    void set_event_side(Seq_side s) { event_side = s; }

    bool operator==(const Rec_Event &) const;
    virtual void update_event_name();
    virtual std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<Seq_type, std::string> &, std::mt19937_64 &) const = 0;
    virtual void write2txt(std::ofstream &) = 0;
    virtual void write2txt_legacy(std::ofstream &) = 0;
    virtual void write2txt_v2(std::ofstream &) = 0;
    virtual void ind_normalize(Marginal_array_p &, size_t) const;
    virtual void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const Events_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &);

private:
    void initialize_event_common(
            std::unordered_set<Rec_Event_name> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Index_map &);

public:
    virtual void initialize_crude_scenario_proba_bound(
            double &, std::forward_list<double *> &,
            const Events_map &);
    virtual void add_to_marginals(long double, Marginal_array_p &) const = 0;
    virtual void set_crude_upper_bound_proba(size_t, size_t, Marginal_array_p &);
    double iterate_common(int realization_index, int base_index,
                          Index_map &base_index_map,
                          const Marginal_array_p &model_parameters);

    /**
     * @brief Update parent realization tracking for marginal indexing
     *
     * Encapsulates the index_map update logic from iterate_common().
     *
     * Updates index_map based on this event's memory_and_offsets structure,
     * which tracks dependencies on parent event realizations.
     *
     * @param realization_index Current realization index
     * @param index_map Index map to update (from ExplorationContext)
     */
    void update_parent_tracking(int realization_index, Index_map& index_map) const {
        for (auto jiter = memory_and_offsets.begin();
             jiter != memory_and_offsets.end(); ++jiter) {
            size_t previous_index =
                index_map.get(std::get<0>(*jiter), std::get<1>(*jiter) - 1);
            previous_index += realization_index * std::get<2>(*jiter);
            index_map.set(std::get<0>(*jiter), previous_index,
                          std::get<1>(*jiter));
        }
    }

    void set_upper_bound_proba(double);
    double get_upper_bound_proba() const { return event_upper_bound_proba; };
    virtual void update_event_internal_probas(const Marginal_array_p &,
                                              const std::unordered_map<Rec_Event_name, int> &);
    //virtual double get_upper_bound_proba() const;
    void set_event_identifier(size_t);
    int get_event_identifier() const;
    void set_event_marginal_size(size_t ev_size) { this->event_marginal_size = ev_size; };
    bool is_updated() const { return updated; };
    void fix(bool fix_status) { fixed = fix_status; }
    bool is_fixed() const { return fixed; }
    void set_viterbi_run(bool viterbi_like) { viterbi_run = viterbi_like; }
    virtual double *get_updated_ptr();
    void compute_crude_upper_bound_scenario_proba(double &);
    const std::vector<int> &get_current_realizations_index_vec() const { return current_realizations_index_vec; };

    //Proba bound related computation methods
    virtual bool has_effect_on(Seq_type) const = 0;
    void iterate_initialize_Len_proba_wrap_up(Seq_type considered_junction,
                                              std::map<int, double> &length_best_proba_map,
                                              std::queue<std::shared_ptr<Rec_Event>> model_queue, double scenario_proba,
                                              const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                              Seq_type_str_p_map &constructed_sequences, int seq_len) const;
    virtual void iterate_initialize_Len_proba(Seq_type considered_junction,
                                              std::map<int, double> &length_best_proba_map,
                                              std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                              double &scenario_proba, const Marginal_array_p &model_parameters_point,
                                              Index_map &base_index_map, Seq_type_str_p_map &constructed_sequences,
                                              int &seq_len) const = 0;
    void iterate_initialize_Len_proba(Seq_type considered_junction, std::map<int, double> &length_best_proba_map,
                                      std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
                                      const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                      Seq_type_str_p_map &constructed_sequences) const;
    virtual void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                            const Marginal_array_p &model_parameters_point,
                                            Index_map &base_index_map) = 0;

protected:
    std::unordered_map<std::string, Event_realization> event_realizations;
    int priority;
    Gene_class event_class;
    Seq_side event_side;
    Rec_Event_name name; //Construct the name in a smart way so that it is unique
    std::string nickname;
    Seq_type_String seq_type; // Seq_type for v2.0 format (e.g., "V_gene_seq", "VD_ins_seq")
    SeqTypeId seq_type_id = kNoSeqType; // resolved from seq_type by Model_Parms::finalize()
    /// Neighbours in the registry ordering, also resolved by Model_Parms::finalize().
    SeqTypeId left_adjacent_id = kNoSeqType;
    SeqTypeId right_adjacent_id = kNoSeqType;
    int len_min;
    int len_max;
    Event_type type;
    int event_index;
    std::forward_list<std::tuple<int, int, int>> memory_and_offsets; //0: event identifier , 1: memory layer , 2: offset
    bool updated;
    bool viterbi_run;
    bool initialized;
    size_t event_marginal_size;
    bool fixed;
    double event_upper_bound_proba;
    double scenario_downstream_upper_bound_proba;
    double scenario_upper_bound_proba; // Used at runtime to store the upper bound probability of the whole scenario
    std::forward_list<double *> updated_proba_bounds_list;
    std::vector<int> current_realizations_index_vec;
    const int *current_realization_index;
    //Snapshot of the downstream proba map's per-key layers, taken at initialize_event().
    std::vector<int> current_downstream_proba_memory_layers;

    int compare_sequences(std::string, std::string); //TODO should probably not be a member functino
    void add_realization(const Event_realization &);

    /**
     * @brief Context-based iterate_wrap_up() interface
     *
     * New signature using 5 context objects instead of 18 individual parameters.
     * Called at leaf nodes to accumulate marginals and update error rate.
     */
    void iterate_wrap_up(
            QuerySequenceContext& query,
            const ModelContext& model,
            ScenarioContext& scenario,
            ExplorationContext& exploration,
            AccumulationContext& accumulation);
};

//bool compare_events(const Rec_Event*&, const Rec_Event*&);
struct Event_comparator
{
    bool operator()(std::shared_ptr<const Rec_Event> event_p1, std::shared_ptr<const Rec_Event> event_p2)
    {
        return event_p1->get_priority() > event_p2->get_priority();
    }
};
