#pragma once

#include <igor/Generation/MarkovSamplingHandler.h>
#include <algorithm>
#include <stdexcept>
#include <typeinfo>

#include <igor/Model/Forward.h>

namespace igor::generation {
using igor::model::RecombinationModel;
using igor::model::Navigator;
using igor::model::Topology;
using igor::model::SampledScenario;
using igor::model::SampledEvent;

template <typename T>
MarkovSamplingHandler<T>::MarkovSamplingHandler(
    std::string name, igor::core::index_type uid, const math::Tensor<T>& weights)
    : SamplingHandler<T>(std::move(name), uid)
    , m_weights(weights)
    , m_state_count(m_weights.shape()[m_weights.ndim() >= 2 ? m_weights.ndim() - 2 : 0])
{
    // Shape: [parent1, ..., from_state, to_state]
    // Parent dimensions are all except last 2
    m_parent_slice_count = 1;
    if (m_weights.ndim() > 2) {
        for (std::size_t d = 0; d + 2 < m_weights.shape().size(); ++d)
            m_parent_slice_count *= m_weights.shape()[d];
    }
}

template <typename T>
const math::Tensor<T>& MarkovSamplingHandler<T>::weights(void) const
{
    return m_weights;
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
    return m_weights.shape();
}

template <typename T>
std::size_t MarkovSamplingHandler<T>::parentSliceOffset(const std::vector<std::size_t>& parent_indices, std::size_t start_dim) const
{
    const auto& sh = m_weights.shape();
    if (sh.size() <= 2) return 0;

    // For shape [parent1, parent2, ..., from_state, to_state]
    // Parent dimensions are [0, ..., ndim-3]
    const std::size_t n_parent_dims = sh.size() - 2;

    // The caller's vector carries the from-state first and the parents after it, so the parent
    // coordinate of dimension d is parent_indices[start_dim + d]. Reading parent_indices[d]
    // instead, as this used to, took the from-state for the first parent's realization.
    if (parent_indices.size() != start_dim + n_parent_dims) {
        throw std::invalid_argument(
            "MarkovSamplingHandler \"" + this->m_name + "\": expected "
            + std::to_string(n_parent_dims) + " parent indices after the from-state, got "
            + std::to_string(parent_indices.size() - std::min(parent_indices.size(), start_dim)));
    }

    // Row-major offset over the parent dimensions only
    std::size_t stride = 1;
    std::size_t offset = 0;
    for (int d = static_cast<int>(n_parent_dims) - 1; d >= 0; --d) {
        offset += parent_indices[start_dim + static_cast<std::size_t>(d)] * stride;
        stride *= sh[d];
    }
    return offset;
}

template <typename T>
void MarkovSamplingHandler<T>::precomputeCDF(void)
{
    // Shape: [parent1, ..., from_state, to_state]
    // Row-major stride: to_state=1, from_state=n_to
    const std::size_t n_from = m_state_count;
    const std::size_t n_to = m_weights.shape().back();
    
    const std::size_t n_rows = n_from * m_parent_slice_count;
    m_row_cdfs = math::Tensor<T>({n_rows, n_to});

    for (std::size_t ps = 0; ps < m_parent_slice_count; ++ps) {
        for (std::size_t from = 0; from < n_from; ++from) {
            const std::size_t row_idx = ps * n_from + from;
            T* cdf_row = m_row_cdfs.data() + row_idx * n_to;
            
            // In row-major [parents..., from, to]:
            // Base offset = ps * (n_from * n_to) + from * n_to
            const T* row_ptr = m_weights.data() + ps * n_from * n_to + from * n_to;
            
            T cumsum = T(0);
            for (std::size_t to = 0; to < n_to; ++to) {
                cumsum += row_ptr[to];
                cdf_row[to] = cumsum;
            }
        }
    }

    // First-nucleotide marginal
    std::vector<T> marginal(n_from, T(0));
    for (std::size_t ps = 0; ps < m_parent_slice_count; ++ps) {
        for (std::size_t from = 0; from < n_from; ++from) {
            const T* row_ptr = m_weights.data() + ps * n_from * n_to + from * n_to;
            for (std::size_t to = 0; to < n_to; ++to)
                marginal[from] += row_ptr[to];
        }
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

    // Same row convention as precomputeCDF(): the parent slice is the outer index, the
    // from-state the inner one, as in the tensor layout [parents..., from, to]. The two used to
    // disagree (ca97599 moved the writer to parents-first and left this reader on the old
    // from-first order), which only showed with at least one parent.
    const std::size_t ps      = parentSliceOffset(parent_indices, 1);
    const std::size_t row_idx = ps * m_state_count + from_state;
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

} // namespace igor::generation
