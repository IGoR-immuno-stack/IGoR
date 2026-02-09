/*
 * Tensor.h
 *
 *  Created on: Feb 06, 2026
 *      Author: IGoR Agent
 */

#pragma once

#include <igor/Math/HybridBuffer.h>
#include <igor/Math/MdspanCompat.h>

#include <vector>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <span>
#include <variant>
#include <array>
#include <initializer_list>

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
    static size_type compute_size(const shape_type& e);

    // Helper to compute strides (Row-Major)
    static shape_type compute_strides(const shape_type& e);

    // Helper to Create View
    template<std::size_t Rank, std::size_t... Is>
    auto make_view_impl(std::index_sequence<Is...>);

    template<std::size_t Rank, std::size_t... Is>
    auto make_view_const_impl(std::index_sequence<Is...>) const;

public:
    // Owning Constructor
    explicit Tensor(shape_type extents);

    // Borrowing Constructor
    Tensor(T* ptr, shape_type extents);

    // Borrowing Constructor (Custom Strides)
    Tensor(T* ptr, shape_type extents, shape_type strides);

    // Default Constructor
    Tensor() = default;

    // Variadic Accessor (tensor(i, j, k))
    template<typename... Indices>
    T& operator()(Indices... indices);
    template<typename... Indices>
    const T& operator()(Indices... indices) const;

    // Runtime Rank Accessor (Direct stride calculation)
    T& operator()(std::span<const size_type> indices);
    const T& operator()(std::span<const size_type> indices) const;

    // Initializer list convenience overloads
    T& operator()(std::initializer_list<size_type> indices);
    const T& operator()(std::initializer_list<size_type> indices) const;

    // Variant Factories
    view_variant variant();
    const_view_variant variant() const;

    // Generic Dispatcher (Visitor Pattern)
    template<typename Func>
    decltype(auto) apply(Func&& f);
    template<typename Func>
    decltype(auto) apply(Func&& f) const;

    // View Accessor (Factory)
    template<std::size_t Rank>
    auto view();
    template<std::size_t Rank>
    auto view() const;

    // Direct Data Access
    T* data() { return storage_.data(); }
    const T* data() const { return storage_.data(); }
    size_type size() const { return storage_.size(); }
    const shape_type& shape() const { return dims_; }
    size_type ndim() const { return dims_.size(); }

    // Iterators (Flat)
    iterator begin() { return storage_.begin(); }
    iterator end() { return storage_.end(); }
    const_iterator begin() const { return storage_.begin(); }
    const_iterator end() const { return storage_.end(); }
    const_iterator cbegin() const { return storage_.cbegin(); }
    const_iterator cend() const { return storage_.cend(); }

#ifndef IGOR_NO_SUBMDSPAN
    /**
     * @brief Get a slice of the tensor along a dimension.
     * Returns a Strided mdspan view of Rank-1.
     */
    template<std::size_t Rank>
    auto slice(size_t dim, size_t index);

    template<std::size_t Rank>
    auto slice(size_t dim, size_t index) const;
#endif
};

} // namespace igor::math

#include <igor/Math/Tensor.tpp>
