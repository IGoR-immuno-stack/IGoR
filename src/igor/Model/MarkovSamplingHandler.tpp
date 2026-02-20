#pragma once

#include "igor/Model/MarkovSamplingHandler.h"
#include <algorithm>
#include <stdexcept>
#include <typeinfo>

namespace igor::model {

template <typename T>
MarkovSamplingHandler<T>::MarkovSamplingHandler(
    std::string name, igor::index_type uid, std::vector<std::size_t> shape)
    : SamplingHandler<T>(std::move(name), uid)
    , m_shape(std::move(shape))
    , m_state_count(m_shape[0])
    , m_parameters(m_shape)        // zero-initialised — caller fills before precompute_cdf()
{
    m_parent_slice_count = 1;
    for (std::size_t d = 2; d < m_shape.size(); ++d)
        m_parent_slice_count *= m_shape[d];
}

template <typename T>
const math::Tensor<T>& MarkovSamplingHandler<T>::parameters(void) const
{
    return m_parameters;
}

template <typename T>
math::Tensor<T>& MarkovSamplingHandler<T>::parameters(void)
{
    return m_parameters;
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::stateCount(void) const
{
    return m_state_count;
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::parentSliceCount(void) const
{
    return m_parent_slice_count;
}

template <typename T>
const std::vector<std::size_t>& MarkovSamplingHandler<T>::shape(void) const
{
    return m_shape;
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::parentSliceOffset(const std::vector<std::size_t>& parent_indices, std::size_t start_dim) const
{
    if (m_shape.size() <= 2) return 0;

    std::size_t stride = 1;
    std::size_t offset = 0;
    for (int d = static_cast<int>(m_shape.size()) - 1;
         d >= static_cast<int>(start_dim) + 1; --d)
    {
        std::size_t pi_idx = static_cast<std::size_t>(d) - start_dim;
        if (pi_idx < parent_indices.size())
            offset += parent_indices[pi_idx] * stride;
        stride *= m_shape[d];
    }
    return offset;
}

template <typename T>
void MarkovSamplingHandler<T>::precomputeCDF(void)
{
    const std::size_t n_rows = m_state_count * m_parent_slice_count;
    m_row_cdfs = math::Tensor<T>({n_rows, m_state_count});

    for (std::size_t from = 0; from < m_state_count; ++from) {
        for (std::size_t ps = 0; ps < m_parent_slice_count; ++ps) {
            const std::size_t row_idx = from * m_parent_slice_count + ps;
            T* cdf_row = m_row_cdfs.data() + row_idx * m_state_count;

            const T* row_ptr = m_parameters.data()
                               + from * m_state_count * m_parent_slice_count
                               + ps   * m_state_count;

            T cumsum = T(0);
            for (std::size_t to = 0; to < m_state_count; ++to) {
                cumsum += row_ptr[to];
                cdf_row[to] = cumsum;
            }
        }
    }

    // First-nucleotide marginal
    std::vector<T> marginal(m_state_count, T(0));
    for (std::size_t from = 0; from < m_state_count; ++from)
        for (std::size_t ps = 0; ps < m_parent_slice_count; ++ps) {
            const T* row_ptr = m_parameters.data()
                               + from * m_state_count * m_parent_slice_count
                               + ps   * m_state_count;
            for (std::size_t to = 0; to < m_state_count; ++to)
                marginal[from] += row_ptr[to];
        }

    T total = T(0);
    for (auto v : marginal) total += v;
    if (total > T(0))
        for (auto& v : marginal) v /= total;

    m_first_cdf.resize(m_state_count);
    T cumsum = T(0);
    for (std::size_t s = 0; s < m_state_count; ++s) {
        cumsum += marginal[s];
        m_first_cdf[s] = cumsum;
    }
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::sampleFromCDF(std::mt19937_64& gen, const T* cdf_row, std::size_t size) const
{
    std::uniform_real_distribution<T> dist(T(0), T(1));
    const T u = dist(gen);
    auto it = std::lower_bound(cdf_row, cdf_row + size, u);
    if (it == cdf_row + size) return size - 1;
    return static_cast<std::size_t>(std::distance(cdf_row, it));
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::sample(std::mt19937_64& gen, const std::vector<std::size_t>& parent_indices) const
{
    if (m_row_cdfs.size() == 0) {
        throw std::logic_error(
            "MarkovSamplingHandler \"" + this->m_name + "\": precomputeCDF() not called");
    }

    if (parent_indices.empty())
        return sampleFromCDF(gen, m_first_cdf.data(), m_first_cdf.size());

    const std::size_t from_state = parent_indices[0];
    if (from_state >= m_state_count)
        throw std::out_of_range(
            "MarkovSamplingHandler \"" + this->m_name + "\": from_state "
            + std::to_string(from_state) + " >= n_states " + std::to_string(m_state_count));

    const std::size_t ps      = parentSliceOffset(parent_indices, 1);
    const std::size_t row_idx = from_state * m_parent_slice_count + ps;
    const T* cdf_row = m_row_cdfs.data() + row_idx * m_state_count;
    return sampleFromCDF(gen, cdf_row, m_state_count);
}

template <typename T>
std::vector<std::size_t> MarkovSamplingHandler<T>::sampleSequence(std::mt19937_64& gen, std::size_t first_state, std::size_t n_steps, const std::vector<std::size_t>& parent_indices) const
{
    std::vector<std::size_t> chain;
    chain.reserve(n_steps + 1);
    chain.push_back(first_state);

    std::size_t from = first_state;
    for (std::size_t step = 0; step < n_steps; ++step) {
        std::vector<std::size_t> step_parents;
        step_parents.reserve(1 + parent_indices.size());
        step_parents.push_back(from);
        step_parents.insert(step_parents.end(), parent_indices.begin(), parent_indices.end());

        from = sample(gen, step_parents);
        chain.push_back(from);
    }
    return chain;
}

} // namespace igor::model
