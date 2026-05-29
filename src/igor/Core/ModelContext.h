#pragma once

#include <igor/Core/Utils.h>
#include <memory>
#include <unordered_map>
#include <vector>
#include <queue>

// Forward declarations
class Rec_Event;

/**
 * @brief Encapsulates read-only model configuration
 * 
 * ModelContext holds immutable model data:
 * - Model probability parameters
 * - Event topology (offset_map, events_map, model_queue, next_event pointers)
 * 
 * Crucially: ModelContext is NEVER modified, even across different
 * sequences in the same EM iteration. This enables:
 * - Thread-safe batch processing
 * - Clear const-correctness
 * - Model parameter caching
 * 
 * Note: offset_map is linked to marginals and will be removed in
 * future refactoring (separate branch).
 * 
 * Note: error_rate has been moved to AccumulationContext because it
 * has accumulation logic that modifies internal state during
 * iterate_wrap_up() at leaf nodes.
 */
struct ModelContext {
    // Model probability parameters (read-only)
    const Marginal_array_p& model_parameters;
    
    // Event dependency offsets (parent -> child index computation)
    const std::unordered_map<Rec_Event_name, 
        std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>& 
        offset_map;
    
    // All events in model (for querying neighbors), keyed by (Event_type, seq_type, Seq_side)
    const Events_map& events_map;

    // Topologically sorted event queue (for iteration order)
    const std::queue<std::shared_ptr<Rec_Event>>& model_queue;

    /**
     * @brief Constructor - binds const references
     */
    ModelContext(
        const Marginal_array_p& model_parameters_,
        const std::unordered_map<Rec_Event_name,
            std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>& offset_map_,
        const Events_map& events_map_,
        const std::queue<std::shared_ptr<Rec_Event>>& model_queue_
    ) : model_parameters(model_parameters_),
        offset_map(offset_map_),
        events_map(events_map_),
        model_queue(model_queue_)
    {}
    
    // Prevent copying and moving (const references)
    ModelContext(const ModelContext&) = delete;
    ModelContext& operator=(const ModelContext&) = delete;
    ModelContext(ModelContext&&) = delete;
    ModelContext& operator=(ModelContext&&) = delete;
};
