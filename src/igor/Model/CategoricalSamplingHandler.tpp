#pragma once

#include <algorithm>
#include <stdexcept>
#include <typeinfo>

namespace igor::model {

// ─── Constructor (borrows tensor from RecombinationModel) ────────────────────

template <typename T>
CategoricalSamplingHandler<T>::CategoricalSamplingHandler(
    std::string name, igor::index_type uid, const math::Tensor<T>& weights)
    : SamplingHandler<T>(std::move(name), uid)
    , m_weights(weights)
    , m_realization_count(m_weights.shape()[0])
{
    m_parent_slice_count = 1;
    for (std::size_t d = 1; d < m_weights.shape().size(); ++d)
        m_parent_slice_count *= m_weights.shape()[d];
}

template <typename T>
const math::Tensor<T>& CategoricalSamplingHandler<T>::weights(void) const
{
    return m_weights;
}

template <typename T>
std::size_t CategoricalSamplingHandler<T>::realizationCount(void) const
{
    return m_realization_count;
}

template <typename T>
std::size_t CategoricalSamplingHandler<T>::parentSliceCount(void) const
{
    return m_parent_slice_count;
}

template <typename T>
const std::vector<std::size_t>& CategoricalSamplingHandler<T>::shape(void) const
{
    return m_weights.shape();
}

// ─── Parent slice offset ──────────────────────────────────────────────────────
//
// Maps parent_indices [p1, p2, ...] → flat index into parent dimensions.

template <typename T>
std::size_t CategoricalSamplingHandler<T>::parentSliceOffset(const std::vector<std::size_t>& parent_indices) const
{
    const auto& sh = m_weights.shape();
    if (sh.size() == 1) return 0;

    if (parent_indices.size() != sh.size() - 1) {
        throw std::invalid_argument(
            "CategoricalSamplingHandler \"" + this->m_name + "\": expected "
            + std::to_string(sh.size() - 1) + " parent indices, got "
            + std::to_string(parent_indices.size()));
    }

    std::size_t stride = 1;
    std::size_t offset = 0;
    for (int d = static_cast<int>(sh.size()) - 1; d >= 1; --d) {
        offset += parent_indices[d - 1] * stride;
        stride *= sh[d];
    }
    return offset;
}

// ─── CDF precomputation ───────────────────────────────────────────────────────

template <typename T>
void CategoricalSamplingHandler<T>::precomputeCDF(void)
{
    m_cdfs = math::Tensor<T>({m_parent_slice_count, m_realization_count});
    for (std::size_t s = 0; s < m_parent_slice_count; ++s) {
        const T* prob = m_weights.data() + s * m_realization_count;
        T* cdf_row = m_cdfs.data() + s * m_realization_count;
        T cumsum = T(0);
        for (std::size_t i = 0; i < m_realization_count; ++i) {
            cumsum     += prob[i];
            cdf_row[i]  = cumsum;
        }
    }
}

// ─── Sampling ─────────────────────────────────────────────────────────────────

template <typename T>
std::size_t CategoricalSamplingHandler<T>::sample(std::mt19937_64& generator, const std::vector<std::size_t>& parent_indices) const
{
    if (m_cdfs.size() == 0) {
        throw std::logic_error(
            "CategoricalSamplingHandler \"" + this->m_name + "\": precomputeCDF() not called");
    }

    const std::size_t s = parentSliceOffset(parent_indices);
    const T* cdf_row = m_cdfs.data() + s * m_realization_count;

    std::uniform_real_distribution<T> dist(T(0), T(1));
    const T u = dist(generator);

    auto it = std::lower_bound(cdf_row, cdf_row + m_realization_count, u);
    if (it == cdf_row + m_realization_count) return m_realization_count - 1;
    return static_cast<std::size_t>(std::distance(cdf_row, it));
}

} // namespace igor::model
