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
    auto registered_double = detail::get_creators<double>().find(type) != detail::get_creators<double>().end();
    auto registered_long_double = detail::get_creators<long double>().find(type) != detail::get_creators<long double>().end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::model::CategoricalInferenceHandler<double>> categorical_registrar{
    GeneChoice_t, Deletion_t, Insertion_t
};
static Registrar<double, igor::model::MarkovInferenceHandler<double>> markov_registrar{
    Dinuclmarkov_t
};
static Registrar<long double, igor::model::CategoricalInferenceHandler<long double>> categorical_registrar_ld{
    GeneChoice_t, Deletion_t, Insertion_t
};
static Registrar<long double, igor::model::MarkovInferenceHandler<long double>> markov_registrar_ld{
    Dinuclmarkov_t
};
}

}
