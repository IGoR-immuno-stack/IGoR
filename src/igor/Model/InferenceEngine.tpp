#pragma once

namespace igor::model {

// ─── Construction from shared RecombinationModel ───────────────────────

template <typename T>
InferenceEngine<T>::InferenceEngine(std::shared_ptr<RecombinationModel<T>> model)
    : m_model(std::move(model))
    , m_handlers(inference_handler_factory::build(*m_model))
    , m_execution_order(m_model->topology().topologicalOrder())
{
}

// ─── Simple accessors ──────────────────────────────────────────────────

template <typename T>
std::size_t InferenceEngine<T>::size(void) const
{
    return m_handlers.size();
}

template <typename T>
const RecombinationModel<T>& InferenceEngine<T>::model(void) const
{
    return *m_model;
}

template <typename T>
RecombinationModel<T>& InferenceEngine<T>::model(void)
{
    return *m_model;
}

// ─── Handler Access ────────────────────────────────────────────────────

template <typename T>
InferenceHandler<T>& InferenceEngine<T>::handler(const std::string& name) {
    const auto& topology = m_model->topology();
    auto id = topology.eventId(name);
    if (id < 0 || id >= static_cast<igor::index_type>(m_handlers.size())) {
        throw std::runtime_error(
            "InferenceEngine: no handler named '" + name + "'");
    }
    return *m_handlers[id];
}

template <typename T>
const InferenceHandler<T>& InferenceEngine<T>::handler(const std::string& name) const {
    const auto& topology = m_model->topology();
    auto id = topology.eventId(name);
    if (id < 0 || id >= static_cast<igor::index_type>(m_handlers.size())) {
        throw std::runtime_error(
            "InferenceEngine: no handler named '" + name + "'");
    }
    return *m_handlers[id];
}

template <typename T>
InferenceHandler<T>& InferenceEngine<T>::handler(igor::index_type uid) {
    if (uid < 0 || uid >= static_cast<igor::index_type>(m_handlers.size())) {
        throw std::out_of_range("InferenceEngine: uid out of range");
    }
    return *m_handlers[uid];
}

template <typename T>
const InferenceHandler<T>& InferenceEngine<T>::handler(igor::index_type uid) const {
    if (uid < 0 || uid >= static_cast<igor::index_type>(m_handlers.size())) {
        throw std::out_of_range("InferenceEngine: uid out of range");
    }
    return *m_handlers[uid];
}

template <typename T>
bool InferenceEngine<T>::hasHandler(const std::string& name) const {
    const auto& topology = m_model->topology();
    auto id = topology.eventId(name);
    return id >= 0 && id < static_cast<igor::index_type>(m_handlers.size());
}

// ─── EM Operations ─────────────────────────────────────────────────────

template <typename T>
void InferenceEngine<T>::resetAccumulators() {
    for (auto& h : m_handlers) {
        h->resetAccumulator();
    }
}

template <typename T>
void InferenceEngine<T>::updateParameters() {
    for (auto& h : m_handlers) {
        h->maximizeLikelihood();
    }
}

template <typename T>
void InferenceEngine<T>::combineAccumulators(const InferenceEngine<T>& other) {
    for (std::size_t i = 0; i < m_handlers.size(); ++i) {
        const auto& name = m_handlers[i]->name();
        m_handlers[i]->accumulator() += other.handler(name).accumulator();
    }
}

// ─── Iteration ─────────────────────────────────────────────────────────

template <typename T>
auto InferenceEngine<T>::orderedHandlers(void) const -> OrderedList
{
    return OrderedList(m_handlers, m_execution_order);
}

template <typename T>
auto InferenceEngine<T>::parents(igor::index_type uid) const -> Adjacency_t
{
    return Adjacency_t(m_handlers, m_model->topology().parentsIds(uid));
}

template <typename T>
auto InferenceEngine<T>::children(igor::index_type uid) const -> Adjacency_t
{
    return Adjacency_t(m_handlers, m_model->topology().childrenIds(uid));
}

template <typename T>
template <typename Func>
void InferenceEngine<T>::forEachHandler(Func&& func) const {
    for (const auto& h : m_handlers) {
        func(h->name(), *h);
    }
}

template <typename T>
template <typename Func>
void InferenceEngine<T>::forEachHandler(Func&& func) {
    for (auto& h : m_handlers) {
        func(h->name(), *h);
    }
}

// ─── run() ─────────────────────────────────────────────────────────────

template <typename T>
template <typename Func>
void InferenceEngine<T>::run(Func&& eStep)
{
    resetAccumulators();
    std::forward<Func>(eStep)(*this);
    updateParameters();
}

} // namespace igor::model
