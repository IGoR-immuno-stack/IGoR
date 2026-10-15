#include <igor/Inference/InferenceHandlerFactory.h>
#include <igor/Inference/CategoricalInferenceHandler.h>
#include <igor/Inference/MarkovInferenceHandler.h>

namespace igor::inference::inference_handler_factory {

namespace detail {
template <typename T>
std::unordered_map<core::EventType, Creator<T>>& get_creators() {
    static std::unordered_map<core::EventType, Creator<T>> creators;
    return creators;
}
template INFERENCE_EXPORT std::unordered_map<core::EventType, Creator<double>>& get_creators<double>();
template INFERENCE_EXPORT std::unordered_map<core::EventType, Creator<long double>>& get_creators<long double>();
}

template void register_creator<double>(core::EventType, Creator<double>);
template void register_creator<long double>(core::EventType, Creator<long double>);

template HandlerPtr<double> create<double>(core::EventType, EventPtr, math::Tensor<double>&);
template HandlerPtr<long double> create<long double>(core::EventType, EventPtr, math::Tensor<long double>&);

bool is_registered(core::EventType type)
{
    auto registered_double = detail::get_creators<double>().find(type) != detail::get_creators<double>().end();
    auto registered_long_double = detail::get_creators<long double>().find(type) != detail::get_creators<long double>().end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::inference::CategoricalInferenceHandler<double>> categorical_registrar{
    core::GeneChoice_t, core::Deletion_t, core::Insertion_t
};
static Registrar<double, igor::inference::MarkovInferenceHandler<double>> markov_registrar{
    core::Dinuclmarkov_t
};
static Registrar<long double, igor::inference::CategoricalInferenceHandler<long double>> categorical_registrar_ld{
    core::GeneChoice_t, core::Deletion_t, core::Insertion_t
};
static Registrar<long double, igor::inference::MarkovInferenceHandler<long double>> markov_registrar_ld{
    core::Dinuclmarkov_t
};
}

}
