#pragma once

#include <igor/Core/Utils.h>
#include <memory>
#include <unordered_map>
#include <vector>

// Forward declarations
class Rec_Event;

/**
 * @brief Encapsulates read-only model configuration
 * 
 * ModelContext holds immutable model data:
 * - Model probability parameters
 * - Event topology (offset_map, events_map, next_event pointers)
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
    
    // Next event to call in iteration order
    const std::shared_ptr<Next_event_ptr>& next_event_ptr_arr;
    
    // All events in model (for querying neighbors)
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
        std::shared_ptr<Rec_Event>>& events_map;
    
    /**
     * @brief Constructor - binds const references
     */
    ModelContext(
        const Marginal_array_p& model_parameters_,
        const std::unordered_map<Rec_Event_name, 
            std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>>& offset_map_,
        const std::shared_ptr<Next_event_ptr>& next_event_ptr_arr_,
        const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
            std::shared_ptr<Rec_Event>>& events_map_
    ) : model_parameters(model_parameters_),
        offset_map(offset_map_),
        next_event_ptr_arr(next_event_ptr_arr_),
        events_map(events_map_)
    {}
    
    // Prevent copying and moving (const references)
    ModelContext(const ModelContext&) = delete;
    ModelContext& operator=(const ModelContext&) = delete;
    ModelContext(ModelContext&&) = delete;
    ModelContext& operator=(ModelContext&&) = delete;
};
