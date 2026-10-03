#pragma once

#include <igor/Generation/SamplingHandler.h>
#include <igor/Math/Tensor.h>

#include <vector>
#include <string>

#include <igor/Model/Forward.h>

namespace igor::generation {
using igor::model::RecombinationModel;
using igor::model::Navigator;
using igor::model::Topology;
using igor::model::SampledScenario;
using igor::model::SampledEvent;

// ─── CategoricalSamplingHandler<T> ───────────────────────────────────────────
//
// Handler for categorical distributions (Gene_choice, Deletion, Insertion).
// Tensor shape: [parent1_size, parent2_size, ..., n_realizations]  (parents first, own last)
//
// The handler borrows a const reference to the probability tensor stored in
// RecombinationModel. It only owns the precomputed CDF tables.
//
// ── Primary lifecycle ───────────────────────────────────────────────────────
//
//   const auto& w = model.weight(uid);          // tensor from RecombinationModel
//   auto sh = CategoricalSamplingHandler<T>(name, uid, w);
//   sh.precomputeCDF();
//   sh.sample(rng, parent_indices);

template <typename T = double>
class CategoricalSamplingHandler : public SamplingHandler<T>
{
public:
    CategoricalSamplingHandler(std::string name, igor::core::index_type uid, const math::Tensor<T>& weights);

    ~CategoricalSamplingHandler(void) = default;

    const math::Tensor<T>& weights(void) const override;

    void precomputeCDF(void) override;

    std::size_t sample(std::mt19937_64&                generator,
                       const std::vector<std::size_t>& parent_indices = {}) const override;

    std::size_t realizationCount(void) const;
    std::size_t parentSliceCount(void) const;
    const std::vector<std::size_t>& shape(void) const;

private:
    const math::Tensor<T>& m_weights; // borrowed from RecombinationModel

    std::size_t m_realization_count = 0;
    std::size_t m_parent_slice_count = 1; // product of parent dimensions

    // cdfs_[slice][j] = P(X <= j | parent_combination = slice)
    math::Tensor<T> m_cdfs;

    std::size_t parentSliceOffset(const std::vector<std::size_t>& parent_indices) const;
};

} // namespace igor::generation

#include <igor/Generation/CategoricalSamplingHandler.tpp>
