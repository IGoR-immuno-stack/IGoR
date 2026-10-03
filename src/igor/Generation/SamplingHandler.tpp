#pragma once

#include <igor/Core/Legacy/Typedef.h>
#include <vector>
#include <random>
#include <stdexcept>

#include <igor/Model/Forward.h>

namespace igor::generation {
using igor::model::RecombinationModel;
using igor::model::Navigator;
using igor::model::Topology;
using igor::model::SampledScenario;
using igor::model::SampledEvent;

template <typename T>
SamplingHandler<T>::SamplingHandler(std::string name, igor::core::legacy::index_type uid) : m_name(std::move(name)), m_uid(uid) 
{

}

template <typename T>
const std::string& SamplingHandler<T>::name(void) const 
{
     return m_name; 
}
    
template <typename T>
igor::core::legacy::index_type  SamplingHandler<T>::uid(void) const 
{
     return m_uid; 
}

template <typename T>
void  SamplingHandler<T>::setUid(igor::core::legacy::index_type id) 
{
     m_uid = id;  
}

// ─── Default sampleSequence: throws for non-Markov types ────────────────────

template <typename T>
std::vector<std::size_t> SamplingHandler<T>::sampleSequence(
    std::mt19937_64&,
    std::size_t,
    std::size_t,
    const std::vector<std::size_t>&) const
{
    throw std::logic_error("sampleSequence() is only valid for MarkovSamplingHandler ("
        + m_name + " is not a Markov handler)");
}

} // namespace igor::generation
