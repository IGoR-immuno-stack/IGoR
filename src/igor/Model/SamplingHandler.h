#pragma once

#include <igor/Core/Typedef.h>

#include <string>
#include <vector>
#include <memory>
#include <random>

namespace igor::model {

// ─── SamplingHandler<T> ───────────────────────────────────────────────────────
//
// Abstract base for all generation handlers.
// One instance per Topology node, owned by SamplingEngine.
//
// Satisfies HasUid: uid_ / uid() / setUid() mirror Rec_Event and MarginalHandler,
// enabling Adjacency<SamplingHandler<T>> for parent traversal in SamplingEngine.
//
// Lifecycle:
//   1. Construct from the corresponding trained MarginalHandler<T>.
//   2. Call precompute_cdf() once (after construction or after marginals are updated).
//   3. Call sample() / sample_sequence() as many times as needed.

template <typename T = double>
class SamplingHandler
{
public:
    using scalar_type = T;

    virtual ~SamplingHandler(void) = default;

    const std::string& name(void) const;

    igor::index_type uid(void) const;
    void setUid(igor::index_type id);

    // ── CDF Precomputation ────────────────────────────────────────────────────
    //
    // Build all CDF slices from the stored probability tensor.
    // Must be called once before any sample() calls.

    virtual void precomputeCDF(void) = 0;

    // ── Sampling ──────────────────────────────────────────────────────────────
    //
    // Draw one realization index given parent realization indices.
    //
    // For categorical events:
    //   parent_indices[k] = realization index of k-th parent event
    //
    // For Markov events — two modes:
    //   parent_indices empty   → sample first nucleotide (from stationary marginal)
    //   parent_indices[0]      → from_state; remaining are upstream parent realizations

    virtual std::size_t sample(
        std::mt19937_64&                 generator,
        const std::vector<std::size_t>&  parent_indices = {}) const = 0;

    // Drive a Markov chain of n steps, starting from first_state.
    // Returns [first_state, next1, next2, ..., nextN] — vector of length n+1.
    // Default: unsupported (throws). Override in MarkovSamplingHandler.

    virtual std::vector<std::size_t> sampleSequence(
        std::mt19937_64&                generator,
        std::size_t                     first_state,
        std::size_t                     n_steps,
        const std::vector<std::size_t>& parent_indices = {}) const;

    virtual       T* rawData(void) = 0;
    virtual const T* rawData(void) const = 0;

    virtual std::size_t rawDataSize(void) const = 0;

protected:
    explicit SamplingHandler(std::string name, igor::index_type uid);

    std::string      m_name;
    igor::index_type m_uid;
};

} // namespace igor::model

#include <igor/Model/SamplingHandler.tpp>
