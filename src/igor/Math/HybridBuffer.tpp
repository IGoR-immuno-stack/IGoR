#pragma once

namespace igor::math {

template <typename T>
HybridBuffer<T>::HybridBuffer(size_type n)
    : data_(new T[n]), size_(n), is_owning_(true) {}

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

} // namespace igor::math
