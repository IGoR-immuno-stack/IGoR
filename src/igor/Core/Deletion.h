/*
 * Deletion.h
 */

#pragma once

#include <igor/Core/Errorrate.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>

#include <forward_list>
#include <list>
#include <math.h>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>

class Deletion : public Rec_Event
{
    friend class Coverage_err_counter;
    friend class Hypermutation_global_errorrate;
    friend class Hypermutation_full_Nmer_errorrate;
public:
    // Constructor
    Deletion();
    Deletion(Gene_class, Seq_side);
    Deletion(Gene_class, Seq_side, std::pair<int, int>); // bounds for deletion lengths
    Deletion(Gene_class, Seq_side, std::unordered_map<std::string, Event_realization> &);

    // Destructor
    virtual ~Deletion();

    // Virtual methods
    std::shared_ptr<Rec_Event> copy() override;

    void iterate(double &, Downstream_scenario_proba_bound_map &, const std::string &, const Int_Str &,
            Index_map &,
            const std::unordered_map<Rec_Event_name,
                                     std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>
                    &,
            std::shared_ptr<Next_event_ptr> &, Marginal_array_p &, const Marginal_array_p &,
            const std::unordered_map<Gene_class, std::vector<Alignment_data>> &,
            Seq_type_str_p_map &, Seq_offsets_map &, std::shared_ptr<Error_rate> &,
            std::map<size_t, std::shared_ptr<Counter>> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>,
                                     std::shared_ptr<Rec_Event>> &,
            Safety_bool_map &, Mismatch_vectors_map &, double &, double &) override;

    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name,
                                     std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>
                    &,
            std::unordered_map<int, std::string> &, std::mt19937_64 &) const override;

    void write2txt(std::ofstream &) override;

    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>
                    &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &,
            std::shared_ptr<Error_rate>, Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;

    void add_to_marginals(long double, Marginal_array_p &) const override;

    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                    const Marginal_array_p &model_parameters_point,
                                    Index_map &base_index_map) override;

    void set_nickname(std::string name) override;
    
    bool has_effect_on(int) const override;
    void iterate_initialize_Len_proba(
            int considered_junction, std::map<int, double> &length_best_proba_map,
            std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
            const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
            Seq_type_str_p_map &constructed_sequences, int &seq_len) const override;

    int get_sequence_type_id() const override { return sequence_type_id; }

    mutable int deletion_value;

private:
    int memory_layer_off_threep;
    int memory_layer_off_fivep;
    int memory_layer_safety_upstream;
    int memory_layer_safety_downstream;
    int memory_layer_off_check1;
    int memory_layer_off_check2;
    int memory_layer_cs;
    int memory_layer_proba_map_junction;
    int memory_layer_proba_map_junction_upstream;

    int sequence_type_id = -1;
    int upstream_seq_type;
    int downstream_seq_type;
    bool upstream_exists;
    bool downstream_exists;
    bool upstream_chosen;
    bool downstream_chosen;
    int upstream_ins_type;
    int downstream_ins_type;

    std::map<int, std::map<int, double>> junction_length_best_proba_maps;
};
