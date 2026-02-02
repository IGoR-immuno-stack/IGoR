/*
 * HypermutationfullNmererrorrate.h
 */

#pragma once

#include <igor/Core/Errorrate.h>
#include <igor/Core/Utils.h>
#include <map>

class Hypermutation_full_Nmer_errorrate : public Error_rate
{
public:
    Hypermutation_full_Nmer_errorrate();
    Hypermutation_full_Nmer_errorrate(size_t, Gene_class, Gene_class, double, size_t = 0);
    Hypermutation_full_Nmer_errorrate(size_t, Gene_class, Gene_class, std::vector<double>, size_t = 0);
    Hypermutation_full_Nmer_errorrate(std::unordered_map<int, Matrix<double>> rates_map);
    virtual ~Hypermutation_full_Nmer_errorrate();

    std::shared_ptr<Error_rate> copy() const override;
    std::string type() const override { return "HypermutationFullNmerErrorrate"; }

    double compare_sequences_error_prob(double, const std::string &, Seq_type_str_p_map &,
                                 const Seq_offsets_map &,
                                 const std::unordered_map<std::tuple<Event_type, int, Seq_side>,
                                                          std::shared_ptr<Rec_Event>> &,
                                 Mismatch_vectors_map &, double &, double &) override;
    void update() override {}
    void initialize(const std::unordered_map<std::tuple<Event_type, int, Seq_side>,
                                             std::shared_ptr<Rec_Event>> &) override {}
    void add_to_norm_counter() override {}
    void clean_seq_counters() override {}
    void write2txt(std::ofstream &) override {}
    Error_rate *add_checked(Error_rate * other) override { return this; }
    const double &get_err_rate_upper_bound(size_t, size_t) override { static double d = 1.0; return d; }
    void build_upper_bound_matrix(size_t, size_t) override {}
    int get_number_non_zero_likelihood_seqs() const override { return 0; }
    std::queue<int> generate_errors(std::string &, std::mt19937_64 &) const override { return std::queue<int>(); }

private:
    std::unordered_map<int, Matrix<double>> error_rates_map;
};
