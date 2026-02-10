#pragma once

namespace igor::model {

// ─── Construction from Descriptors ─────────────────────────────────────

template <typename T>
InferenceEngine<T>::InferenceEngine(const std::vector<EventDescriptor>& descriptors) {
    for (const auto& desc : descriptors) {
        handler_ptr handler;

        // Event_type values from Utils.h:
        // GeneChoice_t = 0, Insertion_t = 1, Deletion_t = 2, Dinuclmarkov_t = 3
        switch (desc.type) {
            case 0: // GeneChoice_t
            case 1: // Insertion_t
            case 2: // Deletion_t
                handler = std::make_unique<CategoricalHandler<T>>(
                    desc.name, desc.shape);
                break;
            case 3: // Dinuclmarkov_t
                handler = std::make_unique<MarkovHandler<T>>(
                    desc.name, desc.shape);
                break;
            default:
                handler = std::make_unique<CategoricalHandler<T>>(
                    desc.name, desc.shape);
                break;
        }

        register_handler(desc.name, std::move(handler));
    }
}

// ─── Handler Registration ──────────────────────────────────────────────

template <typename T>
void InferenceEngine<T>::register_handler(std::string name, handler_ptr handler) {
    if (handlers_.count(name)) {
        throw std::runtime_error(
            "InferenceEngine: duplicate handler name '" + name + "'");
    }
    handlers_[name] = std::move(handler);
    event_order_.push_back(std::move(name));
}

// ─── Handler Access ────────────────────────────────────────────────────

template <typename T>
MarginalHandler<T>& InferenceEngine<T>::handler(const std::string& name) {
    auto it = handlers_.find(name);
    if (it == handlers_.end()) {
        throw std::runtime_error(
            "InferenceEngine: no handler named '" + name + "'");
    }
    return *it->second;
}

template <typename T>
const MarginalHandler<T>& InferenceEngine<T>::handler(const std::string& name) const {
    auto it = handlers_.find(name);
    if (it == handlers_.end()) {
        throw std::runtime_error(
            "InferenceEngine: no handler named '" + name + "'");
    }
    return *it->second;
}

template <typename T>
bool InferenceEngine<T>::has_handler(const std::string& name) const {
    return handlers_.count(name) > 0;
}

// ─── EM Operations ─────────────────────────────────────────────────────

template <typename T>
void InferenceEngine<T>::reset_accumulators() {
    for (const auto& name : event_order_) {
        handlers_.at(name)->reset_accumulator();
    }
}

template <typename T>
void InferenceEngine<T>::update_parameters() {
    for (const auto& name : event_order_) {
        handlers_.at(name)->maximize_likelihood();
    }
}

// ─── I/O ───────────────────────────────────────────────────────────────

template <typename T>
void InferenceEngine<T>::write_marginals(std::ostream& out) const {
    for (const auto& name : event_order_) {
        handlers_.at(name)->write_parameters(out);
    }
}

template <typename T>
void InferenceEngine<T>::read_marginals(std::istream& in) {
    for (const auto& name : event_order_) {
        handlers_.at(name)->read_parameters(in);
    }
}

// ─── Iteration ─────────────────────────────────────────────────────────

template <typename T>
template <typename Func>
void InferenceEngine<T>::for_each_handler(Func&& func) const {
    for (const auto& name : event_order_) {
        func(name, *handlers_.at(name));
    }
}

template <typename T>
template <typename Func>
void InferenceEngine<T>::for_each_handler(Func&& func) {
    for (const auto& name : event_order_) {
        func(name, *handlers_.at(name));
    }
}

} // namespace igor::model
