/*
 * Tensor.h
 *
 *  Created on: Feb 06, 2026
 *      Author: IGoR Agent
 */

#pragma once

#include "HybridBuffer.h"
#include <mdspan>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <utility> // for std::index_sequence
#include <span>
#include <variant>
#include <array>

namespace igor::math {

/**
 * \class Tensor
 * \brief N-dimensional array container with runtime rank and mdspan views.
 *
 * Supports dynamic rank determination at runtime (essential for model loading).
 */
template <typename T>
class Tensor {
public:
    using self_type = Tensor<T>;
    using storage_type = HybridBuffer<T>;
    using value_type = T;
    using size_type = std::size_t;
    using shape_type = std::vector<size_type>; // Dynamic shape storage

    // View Variant Types
    using view_variant = std::variant<
        std::mdspan<T, std::dextents<size_type, 1>>,
        std::mdspan<T, std::dextents<size_type, 2>>,
        std::mdspan<T, std::dextents<size_type, 3>>,
        std::mdspan<T, std::dextents<size_type, 4>>,
        std::mdspan<T, std::dextents<size_type, 5>>
    >;

    using const_view_variant = std::variant<
        std::mdspan<const T, std::dextents<size_type, 1>>,
        std::mdspan<const T, std::dextents<size_type, 2>>,
        std::mdspan<const T, std::dextents<size_type, 3>>,
        std::mdspan<const T, std::dextents<size_type, 4>>,
        std::mdspan<const T, std::dextents<size_type, 5>>
    >;

    // Iterator types
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

private:
    storage_type storage_;
    shape_type dims_;
    shape_type strides_; // Cached strides for fallback access (Rank > 5)

    // Helper to compute size from vector dims
    static size_type compute_size(const shape_type& e) {
        if (e.empty()) return 0;
        return std::accumulate(e.begin(), e.end(), size_type(1), std::multiplies<>());
    }

    // Helper to compute strides (Row-Major)
    static shape_type compute_strides(const shape_type& e) {
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
    template<std::size_t Rank, std::size_t... Is>
    auto make_view_impl(std::index_sequence<Is...>) {
        return std::mdspan<T, std::dextents<size_type, Rank>>(storage_.data(), dims_[Is]...);
    }

    template<std::size_t Rank, std::size_t... Is>
    auto make_view_const_impl(std::index_sequence<Is...>) const {
        return std::mdspan<const T, std::dextents<size_type, Rank>>(storage_.data(), dims_[Is]...);
    }

public:
    // Owning Constructor
    explicit Tensor(shape_type extents)
        : storage_(compute_size(extents)),
          dims_(extents),
          strides_(compute_strides(extents)) {}

    // Borrowing Constructor
    Tensor(T* ptr, shape_type extents)
        : storage_(ptr, compute_size(extents)),
          dims_(extents),
          strides_(compute_strides(extents)) {}

    // Default Constructor
    Tensor() = default;

    // Variadic Accessor (tensor(i, j, k))
    template<typename... Indices>
    T& operator()(Indices... indices) {
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
            // We need to unpack indices into an array/container to iterate with strides
            std::array<size_type, sizeof...(Indices)> idx{static_cast<size_type>(indices)...};
            size_type offset = 0;
            for (size_type i = 0; i < sizeof...(Indices); ++i) {
                offset += idx[i] * strides_[i];
            }
            return storage_[offset];
        }
    }

    template<typename... Indices>
    const T& operator()(Indices... indices) const {
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
    T& operator()(std::span<const size_type> indices) {
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

    const T& operator()(std::span<const size_type> indices) const {
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

    // Variant Factories
    view_variant variant() {
        switch (dims_.size()) {
            case 1: return view<1>();
            case 2: return view<2>();
            case 3: return view<3>();
            case 4: return view<4>();
            case 5: return view<5>();
            default: throw std::runtime_error("Unsupported rank (max 5)");
        }
    }

    const_view_variant variant() const {
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
    template<typename Func>
    decltype(auto) apply(Func&& f) {
        return std::visit(std::forward<Func>(f), variant());
    }

    template<typename Func>
    decltype(auto) apply(Func&& f) const {
        return std::visit(std::forward<Func>(f), variant());
    }

private:
    // Helper to unpack span to view
    template<typename View>
    static auto& apply_span_to_view(View&& v, std::span<const size_type> indices) {
        constexpr size_t Rank = std::remove_cvref_t<View>::rank();
        // We know indices.size() == Rank (checked in operator())
        return apply_indices_impl(std::forward<View>(v), indices, std::make_index_sequence<Rank>{});
    }

    template<typename View, size_t... Is>
    static auto& apply_indices_impl(View&& v, std::span<const size_type> indices, std::index_sequence<Is...>) {
        return v[indices[Is]...];
    }

public:
    // View Accessor (Factory)
    template<std::size_t Rank>
    auto view() {
        if (dims_.size() != Rank) {
            throw std::runtime_error(
                "Tensor rank mismatch: requested Rank=" + std::to_string(Rank) +
                " but tensor has rank=" + std::to_string(dims_.size())
            );
        }
        // Correctly construct extents from dynamic values
        // We use std::extents<size_type, dynamic...> constructor which takes dynamic values
        return make_view_impl<Rank>(std::make_index_sequence<Rank>{});
    }

    template<std::size_t Rank>
    auto view() const {
        if (dims_.size() != Rank) {
            throw std::runtime_error(
                "Tensor rank mismatch: requested Rank=" + std::to_string(Rank) +
                " but tensor has rank=" + std::to_string(dims_.size())
            );
        }
        return make_view_const_impl<Rank>(std::make_index_sequence<Rank>{});
    }

    // Direct Data Access
    T* data() { return storage_.data(); }
    const T* data() const { return storage_.data(); }
    size_type size() const { return storage_.size(); }
    const shape_type& shape() const { return dims_; }
    size_type rank() const { return dims_.size(); }

    // Iterators (Flat)
    iterator begin() { return storage_.begin(); }
    iterator end() { return storage_.end(); }
    const_iterator begin() const { return storage_.begin(); }
    const_iterator end() const { return storage_.end(); }
    const_iterator cbegin() const { return storage_.cbegin(); }
    const_iterator cend() const { return storage_.cend(); }
};

}
