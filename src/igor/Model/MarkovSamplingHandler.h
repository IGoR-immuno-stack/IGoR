#pragma once

#include <igor/Model/SamplingHandler.h>
#include <igor/Math/Tensor.h>

#include <vector>
#include <string>

namespace igor::model {

// ─── MarkovSamplingHandler<T> ─────────────────────────────────────────────────
//
// Handler for Markov transition matrices (Dinucl_markov).
// Tensor shape: [n_states, n_states, parent1_size, ...]  (from, to, parents...)
//
// The handler borrows a const reference to the transition tensor stored in
// RecombinationModel. It only owns the precomputed CDF tables.
//
// ── Primary lifecycle ───────────────────────────────────────────────────────
//
//   const auto& w = model.weight(uid);
//   auto sh = MarkovSamplingHandler<T>(name, uid, w);
//   sh.precomputeCDF();
//   sh.sample(rng, parent_indices);

template <typename T = double>
class MarkovSamplingHandler : public SamplingHandler<T>
{
public:
    MarkovSamplingHandler(std::string name, igor::index_type uid, const math::Tensor<T>& weights);

    ~MarkovSamplingHandler(void) = default;

    const math::Tensor<T>& weights(void) const override;

    void precomputeCDF(void) override;

    std::size_t sample(std::mt19937_64& generator, const std::vector<std::size_t>& parent_indices = {}) const override;

    std::vector<std::size_t> sampleSequence(std::mt19937_64& generator,
        std::size_t first_state,
        std::size_t n_steps,
        const std::vector<std::size_t>& parent_indices = {}) const override;

    std::size_t stateCount(void) const;
    std::size_t parentSliceCount(void) const;
    const std::vector<std::size_t>& shape(void) const;

private:
    const math::Tensor<T>& m_weights; // borrowed from RecombinationModel

    std::size_t m_state_count = 0;
    std::size_t m_parent_slice_count = 1; // product of parent dimensions

    // row_cdfs_[from * n_parent_slices_ + parent_slice][to] = CDF over to-states
    math::Tensor<T> m_row_cdfs;

    // first_cdf_[s] = CDF for first nucleotide (marginal over from_states)
    std::vector<T> m_first_cdf;

    std::size_t parentSliceOffset(const std::vector<std::size_t>& parent_indices, std::size_t start_dim = 1) const;

    std::size_t sampleFromCDF(std::mt19937_64& gen, const T* cdf_row, std::size_t size) const;
};

} // namespace igor::model

#include <igor/Model/MarkovSamplingHandler.tpp>
