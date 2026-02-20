#pragma once

#include <igor/Model/SamplingHandler.h>
#include <igor/Math/Tensor.h>

#include <vector>
#include <string>

namespace igor::model {

// ─── CategoricalSamplingHandler<T> ───────────────────────────────────────────
//
// Handler for categorical distributions (Gene_choice, Deletion, Insertion).
// Tensor shape: [n_realizations, parent1_size, parent2_size, ...]
//
// ── Primary lifecycle (mirrors CategoricalHandler) ──────────────────────────
//
//   auto sh = CategoricalSamplingHandler<T>(name, uid, shape);
//   fill_from_model_marginals(sh.parameters(), marginals);   // same as MarginalHandler
//   sh.precomputeCDF();
//
// ── Validation shortcut ──────────────────────────────────────────────────────
//
//   auto sh = CategoricalSamplingHandler<T>(trained_cat_handler);
//   // parameters copied; caller must call precomputeCDF()
//   // make_sampling_handler(trained_cat_handler) does this automatically.

template <typename T = double>
class CategoricalSamplingHandler : public SamplingHandler<T>
{
public:
    CategoricalSamplingHandler(std::string name, igor::index_type uid, std::vector<std::size_t> shape);

    ~CategoricalSamplingHandler(void) = default;

    const math::Tensor<T>& parameters(void) const;
          math::Tensor<T>& parameters(void);

    void precomputeCDF(void) override;

    std::size_t sample(std::mt19937_64&                generator,
                       const std::vector<std::size_t>& parent_indices = {}) const override;

    std::size_t realizationCount(void) const;
    std::size_t parentSliceCount(void) const;
    const std::vector<std::size_t>& shape(void) const;

          T* rawData(void)       override { return m_parameters.data(); }
    const T* rawData(void) const override { return m_parameters.data(); }

    std::size_t rawDataSize(void) const override { return m_parameters.size(); }

private:
    std::vector<std::size_t> m_shape; // [n_realizations, p1, p2, ...]
    std::size_t m_realization_count = 0;
    std::size_t m_parent_slice_count = 1; // product of parent dimensions

    math::Tensor<T> m_parameters; // probability tensor

    // cdfs_[slice][j] = P(X <= j | parent_combination = slice)
    math::Tensor<T> m_cdfs;

    std::size_t parentSliceOffset(const std::vector<std::size_t>& parent_indices) const;
};

} // namespace igor::model

#include <igor/Model/CategoricalSamplingHandler.tpp>
