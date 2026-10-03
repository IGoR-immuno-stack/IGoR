#include <igor/Inference/InferenceHandlerFactory.h>
#include <igor/Inference/CategoricalInferenceHandler.h>
#include <igor/Inference/MarkovInferenceHandler.h>

namespace igor::inference::inference_handler_factory {

namespace detail {
template <typename T>
std::unordered_map<core::legacy::Event_type, Creator<T>>& get_creators() {
    static std::unordered_map<core::legacy::Event_type, Creator<T>> creators;
    return creators;
}
template INFERENCE_EXPORT std::unordered_map<core::legacy::Event_type, Creator<double>>& get_creators<double>();
template INFERENCE_EXPORT std::unordered_map<core::legacy::Event_type, Creator<long double>>& get_creators<long double>();
}

template void register_creator<double>(core::legacy::Event_type, Creator<double>);
template void register_creator<long double>(core::legacy::Event_type, Creator<long double>);

template HandlerPtr<double> create<double>(core::legacy::Event_type, EventPtr, math::Tensor<double>&);
template HandlerPtr<long double> create<long double>(core::legacy::Event_type, EventPtr, math::Tensor<long double>&);

bool is_registered(core::legacy::Event_type type)
{
    auto registered_double = detail::get_creators<double>().find(type) != detail::get_creators<double>().end();
    auto registered_long_double = detail::get_creators<long double>().find(type) != detail::get_creators<long double>().end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::inference::CategoricalInferenceHandler<double>> categorical_registrar{
    core::legacy::GeneChoice_t, core::legacy::Deletion_t, core::legacy::Insertion_t
};
static Registrar<double, igor::inference::MarkovInferenceHandler<double>> markov_registrar{
    core::legacy::Dinuclmarkov_t
};
static Registrar<long double, igor::inference::CategoricalInferenceHandler<long double>> categorical_registrar_ld{
    core::legacy::GeneChoice_t, core::legacy::Deletion_t, core::legacy::Insertion_t
};
static Registrar<long double, igor::inference::MarkovInferenceHandler<long double>> markov_registrar_ld{
    core::legacy::Dinuclmarkov_t
};
}

}
