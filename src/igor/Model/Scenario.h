#pragma once

#include <igor/Core/Typedef.h>

#include <vector>
#include <string>

namespace igor::model {

// ─── SampledEvent ─────────────────────────────────────────────────────────────
//
// Result of sampling one event node.
// indices holds:
//   - 1 element for categorical events (Gene_choice, Deletion, Insertion)
//   - N elements for Markov events (DinucMarkov): the full nucleotide chain
//     indices[0]   = first nucleotide (sampled from marginal)
//     indices[1..N] = subsequent nucleotides (sampled from transition rows)

struct SampledEvent {
    igor::index_type        event_id;  // Topology node index
    std::vector<std::size_t> indices;  // sampled realization indices
};

// ─── SampledScenario ──────────────────────────────────────────────────────────
//
// All events sampled in one pass, ordered by the Topology node index.
// events[i] corresponds to topology node i.

struct SampledScenario {
    explicit SampledScenario(std::size_t n_events = 0) : events(n_events) {}

    std::vector<SampledEvent> events;

    // Convenience: first index for node i (valid for categorical events)
    std::size_t index_of(igor::index_type i) const {
        return events[static_cast<std::size_t>(i)].indices[0];
    }
};

} // namespace igor::model
