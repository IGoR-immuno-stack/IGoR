/*
 * LayeredArray.h
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>


/// Error paths, kept out of line so that the checks in the hot accessors stay small
/// enough to inline. Building the message inline is enough to stop GCC inlining the
/// caller, which costs several times the price of the check itself.
namespace layered_array_detail {

#if defined(__GNUC__) || defined(__clang__)
#  define IGOR_LA_COLD [[gnu::cold, gnu::noinline]]
#else
#  define IGOR_LA_COLD
#endif

[[noreturn]] IGOR_LA_COLD inline void bad_key(std::size_t key, std::size_t count)
{
    throw std::out_of_range("LayeredArray: key " + std::to_string(key) + " out of range ("
                            + std::to_string(count) + " keys)");
}

[[noreturn]] IGOR_LA_COLD inline void unwritten_key(const char *what, std::size_t key)
{
    throw std::out_of_range(std::string(what) + ": key " + std::to_string(key)
                            + " has never been written");
}

[[noreturn]] IGOR_LA_COLD inline void bad_layer(const char *what, std::size_t key, std::size_t layer)
{
    throw std::out_of_range(std::string(what) + ": layer " + std::to_string(layer)
                            + " is not usable for key " + std::to_string(key));
}

} // namespace layered_array_detail

/**
 * \class LayeredArray LayeredArray.h
 * \brief Runtime-sized, integer-keyed array with a per-key stack of memory layers.
 *
 * Replacement for Enum_fast_memory_map (see docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md,
 * decision D5). It is a pure container: it knows nothing about sequence types or the
 * SeqTypeRegistry. DynamicSequenceMap adds the registry-aware ordered traversal on top,
 * while Index_map -- keyed by event index -- uses this class directly.
 *
 * ### Memory layers
 *
 * A layer is one level of the scenario traversal. Writing a key at layer n leaves its value
 * at layer n-1 intact, so backtracking is a decrement rather than a copy. Layers are tracked
 * **per key independently**: request_layer(k) affects only k.
 *
 * A key that has never been written sits at layer -1 and does not `exists()`. This is
 * deliberately distinct from a key written with an empty value: "not yet processed" versus
 * "actively absent" is the distinction the ordered traversal relies on.
 *
 * ### Storage layout and invariants
 *
 *   storage_[key + layer * count_],  size() == count_ * layer_capacity_
 *
 * Flat and contiguous: one integer multiply-add per access, no pointer chase. The layout is
 * identical to Enum_fast_memory_map's, so the migration is index-for-index.
 *
 *   - `layer_of_[key] == -1`  <=> key never written
 *   - `0 <= layer_of_[key] < layer_capacity_` otherwise
 *   - storage_ always holds exactly `count_ * layer_capacity_` elements
 *
 * ### Reads return by value
 *
 * `get()` returns V rather than V&, because growing the array reallocates and would dangle
 * any outstanding reference. Every value type stored here is small and trivially copyable
 * (pointers, offsets, doubles, bools), so this costs nothing. Writes go through set(), which
 * makes the target layer explicit at the call site.
 */
template <typename V>
class LayeredArray
{
public:
    /// \param count           number of keys; valid keys are [0, count)
    /// \param initial_layers  layers to allocate up front; grows on demand, never shrinks
    explicit LayeredArray(std::size_t count, std::size_t initial_layers = 1)
        : storage_(count * (initial_layers > 0 ? initial_layers : 1)),
          layer_of_(count, -1),
          count_(count),
          layer_capacity_(initial_layers > 0 ? initial_layers : 1)
    { }

    std::size_t count() const noexcept { return count_; }
    std::size_t layer_capacity() const noexcept { return layer_capacity_; }

    /// True once the key has been written at least once.
    bool exists(std::size_t key) const { check_key(key); return layer_of_[key] >= 0; }

    /// Current layer for this key, or -1 if never written.
    int current_layer(std::size_t key) const { check_key(key); return layer_of_[key]; }

    /// Value at the key's current layer.
    /// \throws std::out_of_range if the key is unknown or has never been written.
    V get(std::size_t key) const
    {
        check_key(key);
        const int layer = layer_of_[key];
        if (layer < 0) {
            layered_array_detail::unwritten_key("LayeredArray::get()", key);
        }
        return storage_[index(key, static_cast<std::size_t>(layer))];
    }

    /// Value at an explicit layer. Pure read: unlike Enum_fast_memory_map::at(key, layer),
    /// this does **not** rewind the key's current layer -- call set_current_layer() for that.
    V get(std::size_t key, std::size_t layer) const
    {
        check_key(key);
        if (layer >= layer_capacity_ || static_cast<int>(layer) > layer_of_[key]) {
            layered_array_detail::bad_layer("LayeredArray::get()", key, layer);
        }
        return storage_[index(key, layer)];
    }

    /// Write a value at an explicit layer and make it the key's current layer.
    /// Grows the storage when the layer lies beyond what is allocated.
    /// \throws std::out_of_range if the layer would leave a gap below it.
    void set(std::size_t key, const V &value, std::size_t layer)
    {
        check_key(key);
        // Layers must be filled bottom-up: a key at layer n may only be written at n+1 or
        // below. Writing higher means the caller is using another map's layer numbering.
        if (static_cast<int>(layer) > layer_of_[key] + 1) {
            layered_array_detail::bad_layer("LayeredArray::set()", key, layer);
        }
        ensure_layer(layer);
        storage_[index(key, layer)] = value;
        layer_of_[key] = static_cast<int>(layer);
    }

    /// Push a layer for this key. The new layer's content is **unspecified until written**:
    /// this mirrors Enum_fast_memory_map, whose callers always set() before reading, and
    /// avoids a read-plus-write of the value on every node of the traversal.
    void request_layer(std::size_t key)
    {
        check_key(key);
        const int next = layer_of_[key] + 1;
        ensure_layer(static_cast<std::size_t>(next));
        layer_of_[key] = next;
    }

    /// Pop a layer for this key. Popping layer 0 returns the key to the unwritten state.
    void restore_layer(std::size_t key)
    {
        check_key(key);
        if (layer_of_[key] < 0) {
            layered_array_detail::unwritten_key("LayeredArray::restore_layer()", key);
        }
        --layer_of_[key];
    }

    /// Rewind (or advance) a key to an already-written layer without reading it. This is the
    /// explicit half of Enum_fast_memory_map::at(key, layer), which conflated it with a read.
    void set_current_layer(std::size_t key, std::size_t layer)
    {
        check_key(key);
        if (static_cast<int>(layer) > layer_of_[key]) {
            layered_array_detail::bad_layer("LayeredArray::set_current_layer()", key, layer);
        }
        layer_of_[key] = static_cast<int>(layer);
    }

    /// One layer as a contiguous row over all keys.
    std::span<const V> layer(std::size_t l) const
    {
        if (l >= layer_capacity_) {
            layered_array_detail::bad_layer("LayeredArray::layer()", 0, l);
        }
        return std::span<const V>(storage_.data() + l * count_, count_);
    }

    /// Current layer index of every key, in key order. Replaces the raw int* out-parameter of
    /// Enum_fast_memory_map::get_all_current_memory_layer().
    std::span<const int> current_layers() const noexcept { return std::span<const int>(layer_of_); }

    /// Write `value` at layer 0 for every key, marking them all written.
    void init_first_layer(const V &value)
    {
        for (std::size_t k = 0; k != count_; ++k) {
            storage_[index(k, 0)] = value;
            layer_of_[k] = 0;
        }
    }

    /// Return every key that has been written to layer 0, keeping its value.
    void reset()
    {
        for (auto &l : layer_of_) {
            if (l > 0) {
                l = 0;
            }
        }
    }

private:
    std::size_t index(std::size_t key, std::size_t layer) const noexcept
    {
        return key + layer * count_;
    }

    void check_key(std::size_t key) const
    {
        if (key >= count_) {
            layered_array_detail::bad_key(key, count_);
        }
    }

    /// Make `layer` addressable. Capacity at least doubles, so repeated growth is amortised
    /// O(1); resize (not reserve) is what actually makes the elements exist.
    void ensure_layer(std::size_t layer)
    {
        if (layer < layer_capacity_) {
            return;
        }
        std::size_t new_capacity = layer_capacity_ * 2;
        while (new_capacity <= layer) {
            new_capacity *= 2;
        }
        storage_.resize(count_ * new_capacity);
        layer_capacity_ = new_capacity;
    }

    std::vector<V>   storage_;         ///< storage_[key + layer * count_]
    std::vector<int> layer_of_;        ///< current layer per key; -1 = never written
    std::size_t      count_;
    std::size_t      layer_capacity_;
};
