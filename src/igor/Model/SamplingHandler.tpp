#pragma once

#include "Typedef.h"
#include <vector>
#include <random>
#include <stdexcept>

namespace igor::model {

template <typename T>
SamplingHandler<T>::SamplingHandler(std::string name, igor::index_type uid) : m_name(std::move(name)), m_uid(uid) 
{

}

template <typename T>
const std::string& SamplingHandler<T>::name(void) const 
{
     return m_name; 
}
    
template <typename T>
igor::index_type  SamplingHandler<T>::uid(void) const 
{
     return m_uid; 
}

template <typename T>
void  SamplingHandler<T>::setUid(igor::index_type id) 
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

} // namespace igor::model
