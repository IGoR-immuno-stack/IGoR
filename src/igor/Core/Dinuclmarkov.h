/*
 * Dinuclmarkov.h
 */

#pragma once

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>

#include <igor/Core/Export.h>

class CORE_EXPORT Dinucl_markov : public Rec_Event
{
public:
    Dinucl_markov();
    Dinucl_markov(int);
    Dinucl_markov(int, Seq_side);
    Dinucl_markov(int, Seq_side, std::unordered_map<std::string, Event_realization> &);

    virtual ~Dinucl_markov();

    std::shared_ptr<Rec_Event> copy() override;

    void
    iterate(double &, Downstream_scenario_proba_bound_map &, const std::string &, const Int_Str &, Index_map &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::shared_ptr<Next_event_ptr> &, Marginal_array_p &, const Marginal_array_p &,
            const std::unordered_map<int, std::vector<Alignment_data>> &, Seq_type_str_p_map &,
            Seq_offsets_map &, std::shared_ptr<Error_rate> &, std::map<size_t, std::shared_ptr<Counter>> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
            Safety_bool_map &, Mismatch_vectors_map &, double &, double &) override;

    std::queue<int> draw_random_realization(
            const Marginal_array_p &, std::unordered_map<Rec_Event_name, int> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            std::unordered_map<int, std::string> &, std::mt19937_64 &) const override;

    std::queue<int> draw_random_realization_with_engine(
            const igor::model::InferenceEngine<long double> &,
            std::unordered_map<int, std::string> &,
            std::unordered_map<std::string, std::size_t> &,
            const Model_Parms &,
            std::mt19937_64 &) const override;

    void write2txt(std::ofstream &) override;

    void initialize_event(
            std::unordered_set<Rec_Event_name> &,
            const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
            const std::unordered_map<Rec_Event_name, std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>> &,
            Downstream_scenario_proba_bound_map &, Seq_type_str_p_map &, Safety_bool_map &, std::shared_ptr<Error_rate>,
            Mismatch_vectors_map &, Seq_offsets_map &, Index_map &) override;

    void add_to_marginals(long double, Marginal_array_p &, const Int_Str &, Seq_offsets_map &) const;
    void add_to_marginals(long double, Marginal_array_p &) const override { }

    bool has_effect_on(int) const override { return false; }
    void iterate_initialize_Len_proba(int considered_junction, std::map<int, double> &length_best_proba_map,
                                      std::queue<std::shared_ptr<Rec_Event>> &model_queue, double &scenario_proba,
                                      const Marginal_array_p &model_parameters_point, Index_map &base_index_map,
                                      Seq_type_str_p_map &constructed_sequences, int &seq_len) const override;

    void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>> &model_queue,
                                    const Marginal_array_p &model_parameters_point, Index_map &base_map) override
    {
    }

    void ind_normalize(Marginal_array_p &, size_t) const override;
    void set_nickname(std::string name) override;
    int get_sequence_type_id() const override { return sequence_type_id; }
    int size() const override { return event_realizations.size() * event_realizations.size(); }

private:
    int memory_layer_off_threep;
    int memory_layer_off_fivep;
    int memory_layer_cs;
    int memory_layer_proba_map_junction;

    int sequence_type_id = -1;
    int upstream_seq_type;
    int downstream_seq_type;
    bool upstream_exists;
    bool downstream_exists;

    // Dinucleotide probability matrix (4x4 for ACGT combinations)
    double dinuc_proba_matrix[4][4] = {{0}};

    std::queue<int> draw_random_common(const std::string &previous_seq, std::string &inserted_seq,
                                        const Marginal_array_p &model_marginals_p, int index,
                                        std::uniform_real_distribution<double> &distribution,
                                        std::mt19937_64 &generator) const;
};
