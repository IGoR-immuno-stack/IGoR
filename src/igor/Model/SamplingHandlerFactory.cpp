#include <igor/Model/SamplingHandlerFactory.h>
#include <igor/Model/CategoricalSamplingHandler.h>
#include <igor/Model/MarkovSamplingHandler.h>

namespace igor::model::sampling_handler_factory {

template void register_creator<double>(Event_type, Creator<double>);
template void register_creator<long double>(Event_type, Creator<long double>);

template HandlerPtr<double> create<double>(Event_type, EventPtr, const math::Tensor<double>&);
template HandlerPtr<long double> create<long double>(Event_type, EventPtr, const math::Tensor<long double>&);

bool is_registered(Event_type type)
{
    auto registered_double = detail::creators<double>.find(type) != detail::creators<double>.end();
    auto registered_long_double = detail::creators<long double>.find(type) != detail::creators<long double>.end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::model::CategoricalSamplingHandler<double>> categorial_registrar{
    GeneChoice_t, Deletion_t, Insertion_t
};
static Registrar<double, igor::model::MarkovSamplingHandler<double>> markov_registrar{
    Dinuclmarkov_t
};
}

}