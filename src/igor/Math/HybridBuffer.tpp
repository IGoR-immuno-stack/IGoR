#pragma once

namespace igor::math {

template <typename T>
HybridBuffer<T>::HybridBuffer(size_type n)
    : data_(new T[n]()), size_(n), is_owning_(true) {}  // () after new ensures zero-init

template <typename T>
HybridBuffer<T>::HybridBuffer(T* ptr, size_type n)
    : data_(ptr), size_(n), is_owning_(false) {}

template <typename T>
HybridBuffer<T>::~HybridBuffer() {
    if (is_owning_ && data_) {
        delete[] data_;
    }
}

template <typename T>
HybridBuffer<T>::HybridBuffer(HybridBuffer&& other) noexcept
    : data_(other.data_), size_(other.size_), is_owning_(other.is_owning_) {
    other.data_ = nullptr;
    other.size_ = 0;
    other.is_owning_ = false;
}

template <typename T>
HybridBuffer<T>& HybridBuffer<T>::operator=(HybridBuffer&& other) noexcept {
    if (this != &other) {
        if (is_owning_ && data_) {
            delete[] data_;
        }
        data_ = other.data_;
        size_ = other.size_;
        is_owning_ = other.is_owning_;

        other.data_ = nullptr;
        other.size_ = 0;
        other.is_owning_ = false;
    }
    return *this;
}

// Copy Constructor — always deep-copies into new owning memory
template <typename T>
HybridBuffer<T>::HybridBuffer(const HybridBuffer& other)
    : data_(other.size_ > 0 ? new T[other.size_] : nullptr),
      size_(other.size_),
      is_owning_(true) {
    for (size_type i = 0; i < size_; ++i) {
        data_[i] = other.data_[i];
    }
}

// Copy Assignment — always deep-copies into new owning memory
template <typename T>
HybridBuffer<T>& HybridBuffer<T>::operator=(const HybridBuffer& other) {
    if (this != &other) {
        // Allocate new buffer first (strong exception safety)
        T* new_data = other.size_ > 0 ? new T[other.size_] : nullptr;
        for (size_type i = 0; i < other.size_; ++i) {
            new_data[i] = other.data_[i];
        }

        // Release old memory
        if (is_owning_ && data_) {
            delete[] data_;
        }

        data_ = new_data;
        size_ = other.size_;
        is_owning_ = true;
    }
    return *this;
}

} // namespace igor::math
