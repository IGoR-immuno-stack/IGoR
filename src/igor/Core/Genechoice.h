/*
 * Genechoice.h
 */

#pragma once

#include <igor/Core/Errorrate.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>

#include <igor/Core/Export.h>

class CORE_EXPORT Gene_choice : public Rec_Event
{
    friend class Coverage_err_counter;
    friend class Hypermutation_global_errorrate;
    friend class Hypermutation_full_Nmer_errorrate;

public:
    // Constructors
    Gene_choice();
    Gene_choice(int);
    Gene_choice(int, std::unordered_map<std::string, Event_realization> &);
    Gene_choice(int, std::vector<std::pair<std::string, std::string>>);
    // Destructor
    virtual ~Gene_choice();
    // Virtual methods overload
    std::shared_ptr<Rec_Event> copy() override;
    void
    iterate(double &, Downstream_scenario_proba_bound_map &, const std::string &, const Int_Str &, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::shared_ptr<Next_event_ptr> &, Marginal_array_p &, const Marginal_array_p &,
            const std::unordered_map<int, std::vector<Alignment_data>> &, Seq_type_str_p_map &,
            Seq_offsets_map &, std::shared_ptr<Error_rate> &, std::map<size_t, std::shared_ptr<Counter>> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
            Safety_bool_map &, Mismatch_vectors_map &, double &, double &) override;

    bool add_realization(std::string gene_name, std::string gene_sequence);
    void set_genomic_templates(const std::vector<std::pair<std::string, std::string>> &);
    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<int, std::string> &, std::mt19937_64 &) const override;
    void write2txt(std::ofstream &) override;
    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;
    void add_to_marginals(long double, Marginal_array_p &) const override;
    void update_event_internal_probas(const Marginal_array_p &,
                                      const std::unordered_map<Rec_Event_name, int> &) override;
    void set_nickname(std::string name) override;
    int get_sequence_type_id() const override { return sequence_type_id; }

    // Proba bound related computation methods
    bool has_effect_on(int) const override;
    void iterate_initialize_Len_proba(int considered_junction, std::map<int, double> &length_best_proba_map,
                                      std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
                                      const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                      Seq_type_str_p_map &constructed_sequences, int &seq_len) const override;
    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                    const Marginal_array_p &model_parameters_point, Index_map &base_index_map) override;

private:
    double iterate_common(
            double, const int &, int, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            const Marginal_array_p &);

    // Inference variables
    bool vd_check;
    bool vj_check;
    bool dj_check;

    Seq_Offset d_5_min_offset;
    Seq_Offset d_5_max_offset;
    Seq_Offset j_5_min_offset;
    Seq_Offset j_5_max_offset;
    Seq_Offset v_5_off;
    Seq_Offset v_3_off;
    Seq_Offset d_offset;
    Seq_Offset j_offset;
    Seq_Offset v_offset;
    Seq_Offset v_3_min_offset;
    Seq_Offset v_3_max_offset;
    Seq_Offset d_3_off;
    Seq_Offset d_5_off;
    Seq_Offset d_3_min_offset;
    Seq_Offset d_3_max_offset;
    Seq_Offset j_5_off;

    bool no_d_align;
    size_t d_size;
    Seq_Offset d_full_3_offset;

    mutable int base_index;
    double new_scenario_proba;
    double new_tmp_err_w_proba;
    double proba_contribution;
    int new_index;
    const int *alignment_offset_p;

    int memory_layer_cs;
    int memory_layer_mismatches;
    int memory_layer_safety_1;
    int memory_layer_safety_2;
    int memory_layer_off_threep;
    int memory_layer_off_fivep;
    int memory_layer_offset_check1;
    int memory_layer_offset_check2;
    int memory_layer_proba_map_seq;
    int memory_layer_proba_map_junction;
    int memory_layer_proba_map_junction_upstream;

    bool v_chosen;
    bool v_choice_exist;
    bool d_chosen;
    bool d_choice_exist;
    bool j_chosen;
    bool j_choice_exist;

    int sequence_type_id = -1;
    int upstream_seq_type = -1;
    int downstream_seq_type = -1;
    bool upstream_chosen = false;
    bool upstream_exists = false;
    bool downstream_chosen = false;
    bool downstream_exists = false;
    int upstream_ins_type = -1;
    int downstream_ins_type = -1;
    int memory_layer_safety_upstream = -1;
    int memory_layer_safety_downstream = -1;

    int d_5_max_del;
    int d_5_min_del;
    int d_5_real_max_del;
    int j_5_max_del;
    int j_5_min_del;
    int v_3_max_del;
    int v_3_min_del;
    int d_3_max_del;
    int d_3_min_del;

    std::vector<std::pair<int, int>> active_upstream_junctions;
    std::vector<std::pair<int, int>> active_downstream_junctions;
    std::map<int, std::map<int, double>> junction_length_best_proba_maps;
    std::map<int, std::vector<std::tuple<std::string, int, int, double>>> vj_length_d_position_proba;
};
