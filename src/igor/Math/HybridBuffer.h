#pragma once

#include <cstddef>
#include <utility>

namespace igor::math {

/**
 * \class HybridBuffer
 * \brief A container that can either own its memory or view existing memory.
 *
 * solves the std::vector "ownership forced" problem.
 */
template <typename T>
class HybridBuffer {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using iterator = pointer;
    using const_iterator = const_pointer;

private:
    T* data_ = nullptr;
    size_type size_ = 0;
    bool is_owning_ = false;

public:
    // Default constructor
    HybridBuffer() = default;

    // Owning constructor
    explicit HybridBuffer(size_type n);

    // Borrowing constructor
    HybridBuffer(T* ptr, size_type n);

    // Destructor
    ~HybridBuffer();

    // Move Constructor
    HybridBuffer(HybridBuffer&& other) noexcept;

    // Move Assignment
    HybridBuffer& operator=(HybridBuffer&& other) noexcept;

    // Copy Constructor (always deep-copies into owning memory)
    HybridBuffer(const HybridBuffer& other);

    // Copy Assignment (always deep-copies into owning memory)
    HybridBuffer& operator=(const HybridBuffer& other);

    // Accessors
    pointer data() noexcept { return data_; }
    const_pointer data() const noexcept { return data_; }
    size_type size() const noexcept { return size_; }
    bool is_owning() const noexcept { return is_owning_; }

    // Iterators
    iterator begin() noexcept { return data_; }
    iterator end() noexcept { return data_ + size_; }
    const_iterator begin() const noexcept { return data_; }
    const_iterator end() const noexcept { return data_ + size_; }
    const_iterator cbegin() const noexcept { return data_; }
    const_iterator cend() const noexcept { return data_ + size_; }

    // Element Access
    reference operator[](size_type pos) { return data_[pos]; }
    const_reference operator[](size_type pos) const { return data_[pos]; }
};

} // namespace igor::math

#include "HybridBuffer.tpp"
