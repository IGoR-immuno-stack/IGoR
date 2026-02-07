#pragma once

namespace igor::math {

// Helper to compute size from vector dims
template <typename T>
typename Tensor<T>::size_type Tensor<T>::compute_size(const shape_type& e) {
    if (e.empty()) return 0;
    return std::accumulate(e.begin(), e.end(), size_type(1), std::multiplies<>());
}

// Helper to compute strides (Row-Major)
template <typename T>
typename Tensor<T>::shape_type Tensor<T>::compute_strides(const shape_type& e) {
    shape_type s(e.size());
    if (e.empty()) return s;
    
    size_type running_stride = 1;
    for (std::size_t i = e.size(); i > 0; --i) {
        s[i - 1] = running_stride;
        running_stride *= e[i - 1];
    }
    return s;
}

// Helper to Create View
template <typename T>
template<std::size_t Rank, std::size_t... Is>
auto Tensor<T>::make_view_impl(std::index_sequence<Is...>) {
    return std::mdspan<T, std::dextents<size_type, Rank>>(storage_.data(), dims_[Is]...);
}

template <typename T>
template<std::size_t Rank, std::size_t... Is>
auto Tensor<T>::make_view_const_impl(std::index_sequence<Is...>) const {
    return std::mdspan<const T, std::dextents<size_type, Rank>>(storage_.data(), dims_[Is]...);
}

// Owning Constructor
template <typename T>
Tensor<T>::Tensor(shape_type extents)
    : storage_(compute_size(extents)), 
      dims_(extents),
      strides_(compute_strides(extents)) {}

// Borrowing Constructor
template <typename T>
Tensor<T>::Tensor(T* ptr, shape_type extents)
    : storage_(ptr, compute_size(extents)), 
      dims_(extents),
      strides_(compute_strides(extents)) {}

// Variadic Accessor (tensor(i, j, k))
template <typename T>
template<typename... Indices>
T& Tensor<T>::operator()(Indices... indices) {
    if (sizeof...(Indices) != dims_.size()) {
        throw std::runtime_error(
            "Tensor rank mismatch: provided " + std::to_string(sizeof...(Indices)) +
            " indices but tensor has rank " + std::to_string(dims_.size())
        );
    }
#ifdef IGOR_DEBUG
    std::array<size_type, sizeof...(Indices)> idx{static_cast<size_type>(indices)...};
    for (size_t i = 0; i < sizeof...(Indices); ++i) {
        if (idx[i] >= dims_[i]) {
            throw std::out_of_range(
                "Index out of bounds: index[" + std::to_string(i) + "]=" +
                std::to_string(idx[i]) + " >= extent[" + std::to_string(i) + "]=" +
                std::to_string(dims_[i])
            );
        }
    }
#endif
    if constexpr (sizeof...(Indices) <= 5) {
        return view<sizeof...(Indices)>()[indices...];
    } else {
        // Fallback: Manual Strides (Rank > 5)
        std::array<size_type, sizeof...(Indices)> idx{static_cast<size_type>(indices)...};
        size_type offset = 0;
        for (size_type i = 0; i < sizeof...(Indices); ++i) {
            offset += idx[i] * strides_[i];
        }
        return storage_[offset];
    }
}

template <typename T>
template<typename... Indices>
const T& Tensor<T>::operator()(Indices... indices) const {
    if (sizeof...(Indices) != dims_.size()) {
        throw std::runtime_error(
            "Tensor rank mismatch: provided " + std::to_string(sizeof...(Indices)) +
            " indices but tensor has rank " + std::to_string(dims_.size())
        );
    }
#ifdef IGOR_DEBUG
    std::array<size_type, sizeof...(Indices)> idx{static_cast<size_type>(indices)...};
    for (size_t i = 0; i < sizeof...(Indices); ++i) {
        if (idx[i] >= dims_[i]) {
            throw std::out_of_range(
                "Index out of bounds: index[" + std::to_string(i) + "]=" +
                std::to_string(idx[i]) + " >= extent[" + std::to_string(i) + "]=" +
                std::to_string(dims_[i])
            );
        }
    }
#endif
    if constexpr (sizeof...(Indices) <= 5) {
        return view<sizeof...(Indices)>()[indices...];
    } else {
        // Fallback: Manual Strides (Rank > 5)
        std::array<size_type, sizeof...(Indices)> idx{static_cast<size_type>(indices)...};
        size_type offset = 0;
        for (size_type i = 0; i < sizeof...(Indices); ++i) {
            offset += idx[i] * strides_[i];
        }
        return storage_[offset];
    }
}

// Runtime Rank Accessor (Direct stride calculation)
template <typename T>
T& Tensor<T>::operator()(std::span<const size_type> indices) {
    if (indices.size() != dims_.size()) {
        throw std::runtime_error(
            "Tensor rank mismatch: provided " + std::to_string(indices.size()) +
            " indices but tensor has rank " + std::to_string(dims_.size())
        );
    }

    // Direct stride calculation (efficient for all ranks)
    size_type offset = 0;
    for (size_type i = 0; i < indices.size(); ++i) {
        offset += indices[i] * strides_[i];
    }
    return storage_[offset];
}

template <typename T>
const T& Tensor<T>::operator()(std::span<const size_type> indices) const {
    if (indices.size() != dims_.size()) {
        throw std::runtime_error(
            "Tensor rank mismatch: provided " + std::to_string(indices.size()) +
            " indices but tensor has rank " + std::to_string(dims_.size())
        );
    }

    // Direct stride calculation (efficient for all ranks)
    size_type offset = 0;
    for (size_type i = 0; i < indices.size(); ++i) {
        offset += indices[i] * strides_[i];
    }
    return storage_[offset];
}

template <typename T>
T& Tensor<T>::operator()(std::initializer_list<size_type> indices) {
    return operator()(std::span<const size_type>{indices.begin(), indices.size()});
}

template <typename T>
const T& Tensor<T>::operator()(std::initializer_list<size_type> indices) const {
    return operator()(std::span<const size_type>{indices.begin(), indices.size()});
}

// Variant Factories
template <typename T>
typename Tensor<T>::view_variant Tensor<T>::variant() {
    switch (dims_.size()) {
        case 1: return view<1>();
        case 2: return view<2>();
        case 3: return view<3>();
        case 4: return view<4>();
        case 5: return view<5>();
        default: throw std::runtime_error("Unsupported rank (max 5)");
    }
}

template <typename T>
typename Tensor<T>::const_view_variant Tensor<T>::variant() const {
    switch (dims_.size()) {
        case 1: return view<1>();
        case 2: return view<2>();
        case 3: return view<3>();
        case 4: return view<4>();
        case 5: return view<5>();
        default: throw std::runtime_error("Unsupported rank (max 5)");
    }
}

// Generic Dispatcher (Visitor Pattern)
template <typename T>
template<typename Func>
decltype(auto) Tensor<T>::apply(Func&& f) {
    return std::visit(std::forward<Func>(f), variant());
}

template <typename T>
template<typename Func>
decltype(auto) Tensor<T>::apply(Func&& f) const {
    return std::visit(std::forward<Func>(f), variant());
}

// View Accessor (Factory)
template <typename T>
template<std::size_t Rank>
auto Tensor<T>::view() {
    if (dims_.size() != Rank) {
        throw std::runtime_error(
            "Tensor rank mismatch: requested Rank=" + std::to_string(Rank) +
            " but tensor has rank=" + std::to_string(dims_.size())
        );
    }
    return make_view_impl<Rank>(std::make_index_sequence<Rank>{});
}

template <typename T>
template<std::size_t Rank>
auto Tensor<T>::view() const {
    if (dims_.size() != Rank) {
        throw std::runtime_error(
            "Tensor rank mismatch: requested Rank=" + std::to_string(Rank) +
            " but tensor has rank=" + std::to_string(dims_.size())
        );
    }
    return make_view_const_impl<Rank>(std::make_index_sequence<Rank>{});
}

} // namespace igor::math
