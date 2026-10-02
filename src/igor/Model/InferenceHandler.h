// InferenceHandler.h ---

#pragma once

#include <igor/Math/Tensor.h>
#include <igor/Core/Typedef.h>

#include <string>

namespace igor::model {

// ─── InferenceHandler<T> ─────────────────────────────────────────────────────
//
// Abstract base for all EM inference handlers.
// One instance per Topology node, owned by InferenceEngine.
//
// Each handler holds a mutable reference to the probability tensor stored in
// RecombinationModel (the single source of truth). The handler only owns
// the accumulator tensor needed for the E-step.
//
// Lifecycle:
//   1. Construct with a mutable reference to the model's weight tensor.
//   2. Repeat per EM iteration:
//      a. reset_accumulator()              — zero the accumulator
//      b. ... E-step accumulates counts ...
//      c. maximize_likelihood()            — normalise accumulator → weights
//
// The M-step writes normalised values directly into the model's weight tensor,
// so all engines/handlers sharing the same RecombinationModel see the update.

template <typename T = double>
class InferenceHandler
{
public:
    using scalar_type = T;

    virtual ~InferenceHandler(void) = default;

    // ── Identity ─────────────────────────────────────────────────────────

    const std::string& name(void) const { return m_name; }

    igor::index_type uid(void) const { return m_uid; }
    void setUid(igor::index_type id) { m_uid = id; }

    // ── Tensor access ────────────────────────────────────────────────────

    virtual const math::Tensor<T>& weights(void) const = 0;
    virtual math::Tensor<T>& weights(void) = 0;

    virtual const math::Tensor<T>& accumulator(void) const = 0;
    virtual       math::Tensor<T>& accumulator(void) = 0;

    // ── EM operations ────────────────────────────────────────────────────

    virtual void resetAccumulator(void) = 0;
    virtual void maximizeLikelihood(void) = 0;

protected:
    explicit InferenceHandler(std::string name, igor::index_type uid)
        : m_name(std::move(name))
        , m_uid(uid) {}

    std::string      m_name;
    igor::index_type m_uid;
};

} // namespace igor::model

//
// InferenceHandler.h ends here
