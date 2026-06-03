// InferenceHandlerFactory.h ---

#pragma once

#include <igor/Model/Export.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Model/InferenceHandler.h>
#include <igor/Math/Tensor.h>
#include <igor/Core/Utils.h>

#include <memory>
#include <functional>
#include <string>
#include <initializer_list>

namespace igor::model::inference_handler_factory {

using EventPtr = std::shared_ptr<Rec_Event>;

template <typename T>
using HandlerPtr = std::unique_ptr<InferenceHandler<T>>;

template <typename T>
using Creator = std::function<HandlerPtr<T>(EventPtr, math::Tensor<T>&)>;

namespace detail {
template <typename T>
MODEL_EXPORT std::unordered_map<Event_type, Creator<T>>& get_creators();
}

/**
 * @brief Register a creator function for an event type
 * @param type The Event_type enum value
 * @param func Function that creates a handler borrowing a mutable tensor reference
 */
template <typename T>
void register_creator(Event_type type, Creator<T> func);

/**
 * @brief Create a handler of the specified type
 * @param type The Event_type to create
 * @param event The Rec_Event for metadata (name, uid)
 * @param weights Mutable reference to the probability tensor from RecombinationModel
 * @return Unique pointer to the new handler
 * @throws std::runtime_error if type not registered
 */
template <typename T>
HandlerPtr<T> create(Event_type type, EventPtr event, math::Tensor<T>& weights);

/**
 * @brief Check if an event type has a registered creator
 * @param type The Event_type to check
 * @return true if registered, false otherwise
 */
MODEL_EXPORT bool is_registered(Event_type type);

} // namespace igor::model::inference_handler_factory

namespace igor::model {
    template <typename T> class RecombinationModel; // Forward declaration
}

namespace igor::model::inference_handler_factory {

/**
 * @brief Build a list of inference handlers borrowing tensors from a RecombinationModel
 * @param model The RecombinationModel providing topology and weight tensors
 * @return A vector of inference handlers, one per topology node
 *
 * Each handler holds a mutable reference to its corresponding tensor in the model.
 * The model must outlive the returned handlers.
 */
template <typename T>
std::vector<HandlerPtr<T>> build(igor::model::RecombinationModel<T>& model);

/**
 * @brief Template for automatic event registration using static initialization
 * @tparam T The scalar type (double, long double)
 * @tparam HandlerClass The concrete handler class
 *
 * Usage:
 *   static Registrar<double, CategoricalInferenceHandler<double>> cat{GeneChoice_t, Deletion_t};
 */
template<typename T, typename HandlerClass>
struct Registrar {
    Registrar(std::initializer_list<Event_type> types) {
        for (auto type : types) {
            register_creator<T>(type, [](EventPtr event, math::Tensor<T>& weights) {
                return std::make_unique<HandlerClass>(event->get_nickname(), event->uid(), weights);
            });
        }
    }
};

}

#include <igor/Model/InferenceHandlerFactory.tpp>

//
// InferenceHandlerFactory.h ends here
