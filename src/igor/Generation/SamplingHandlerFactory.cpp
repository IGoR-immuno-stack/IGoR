#include <igor/Generation/SamplingHandlerFactory.h>
#include <igor/Generation/CategoricalSamplingHandler.h>
#include <igor/Generation/MarkovSamplingHandler.h>

namespace igor::generation::sampling_handler_factory {

namespace detail {
template <typename T>
std::unordered_map<core::legacy::Event_type, Creator<T>>& get_creators() {
    static std::unordered_map<core::legacy::Event_type, Creator<T>> creators;
    return creators;
}
template GENERATION_EXPORT std::unordered_map<core::legacy::Event_type, Creator<double>>& get_creators<double>();
template GENERATION_EXPORT std::unordered_map<core::legacy::Event_type, Creator<long double>>& get_creators<long double>();
}

template void register_creator<double>(core::legacy::Event_type, Creator<double>);
template void register_creator<long double>(core::legacy::Event_type, Creator<long double>);

template HandlerPtr<double> create<double>(core::legacy::Event_type, EventPtr, const math::Tensor<double>&);
template HandlerPtr<long double> create<long double>(core::legacy::Event_type, EventPtr, const math::Tensor<long double>&);

bool is_registered(core::legacy::Event_type type)
{
    auto registered_double = detail::get_creators<double>().find(type) != detail::get_creators<double>().end();
    auto registered_long_double = detail::get_creators<long double>().find(type) != detail::get_creators<long double>().end();
    return registered_double || registered_long_double;
}

namespace {
static Registrar<double, igor::generation::CategoricalSamplingHandler<double>> categorial_registrar{
    core::legacy::GeneChoice_t, core::legacy::Deletion_t, core::legacy::Insertion_t
};
static Registrar<double, igor::generation::MarkovSamplingHandler<double>> markov_registrar{
    core::legacy::Dinuclmarkov_t
};
}

}