/*
 * HypermutationfullNmererrorrate.cpp
 */

#include <igor/Core/HypermutationfullNmererrorrate.h>
#include <igor/Core/SequenceTypes.h>

using namespace std;

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate() : Error_rate() { }

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t n, Gene_class g1, Gene_class g2, double d,
                                                                     size_t s)
    : Error_rate()
{
}

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t n, Gene_class g1, Gene_class g2,
                                                                     std::vector<double> v, size_t s)
    : Error_rate()
{
}

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(std::unordered_map<int, Matrix<double>> rates_map)
    : error_rates_map(rates_map)
{
}

Hypermutation_full_Nmer_errorrate::~Hypermutation_full_Nmer_errorrate() { }

shared_ptr<Error_rate> Hypermutation_full_Nmer_errorrate::copy() const
{
    return shared_ptr<Error_rate>(new Hypermutation_full_Nmer_errorrate(this->error_rates_map));
}

double Hypermutation_full_Nmer_errorrate::compare_sequences_error_prob(
        double scenario_proba, const std::string &sequence, Seq_type_str_p_map &constructed_sequences,
        const Seq_offsets_map &seq_offsets,
        const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Mismatch_vectors_map &mismatches_lists, double &seq_max_prob_scenario, double &proba_threshold_factor)
{
    double proba = 1.0;
    // Stub implementation
    return (long double)scenario_proba * proba;
}
