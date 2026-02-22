#include <igor/Model/InferenceHandlerFactory.h>
#include <igor/Model/CategoricalInferenceHandler.h>
#include <igor/Model/MarkovInferenceHandler.h>

namespace igor::model::inference_handler_factory {

template void register_creator<double>(Event_type, Creator<double>);
template void register_creator<long double>(Event_type, Creator<long double>);

template HandlerPtr<double> create<double>(Event_type, EventPtr, math::Tensor<double>&);
template HandlerPtr<long double> create<long double>(Event_type, EventPtr, math::Tensor<long double>&);

bool is_registered(Event_type type)
{
    auto registered_double = detail::creators<double>.find(type) != detail::creators<double>.end();
    auto registered_long_double = detail::creators<long double>.find(type) != detail::creators<long double>.end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::model::CategoricalInferenceHandler<double>> categorical_registrar{
    GeneChoice_t, Deletion_t, Insertion_t
};
static Registrar<double, igor::model::MarkovInferenceHandler<double>> markov_registrar{
    Dinuclmarkov_t
};
}

}
