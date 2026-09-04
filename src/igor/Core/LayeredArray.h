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

#include <cassert>
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
 * ### Requested layers and written layers are tracked separately
 *
 * Two different questions get asked about a key's layers, and conflating them made the
 * container unable to answer either reliably:
 *
 *   - *which layer do I own?*  -- asked at initialization, by an event that has just called
 *     request_layer(), or that wants the layer its neighbour owns. Answered by
 *     `claimed_layer()` / `claimed_layers()`.
 *   - *where does this key's data currently stand?* -- asked at every read, and by
 *     `exists()`. Answered by `current_layer()`, which moves up and down as the traversal
 *     writes and backtracks.
 *
 * `request_layer()` raises only the first. That matters in three ways:
 *
 *   1. `get(key, layer)` on a layer that was requested but never written now always throws.
 *      When the two were one counter, requesting made the layer readable and it returned
 *      value-initialized storage; whether a missing write was caught depended on whether
 *      some *other* event's write happened to have pulled the counter back down first. See
 *      docs/ITERATE_GENERIC_REWRITE_PLAN.md section 7.9 for the bug that exposed this.
 *   2. `exists()` means written, not "requested or written". Its callers -- the Scenario
 *      view, Single_error_rate, the coverage counters, DynamicSequenceMap::occupied() --
 *      all guard dereferences with it, and a requested-but-unwritten key handed them a
 *      value-initialized null.
 *   3. The three-state distinction above becomes real. Under the old scheme a
 *      not-yet-processed key that some event had requested was indistinguishable from an
 *      actively-absent one, because both read back as a written default.
 *
 * Invariant: `current_layer_[key] <= claimed_layer_[key]`. Writing at a layer claims it,
 * so `set()` raises the requested mark when it has to.
 *
 * ### Storage layout and invariants
 *
 *   storage_[key + layer * count_],  size() == count_ * layer_capacity_
 *
 * Flat and contiguous: one integer multiply-add per access, no pointer chase. The layout is
 * identical to Enum_fast_memory_map's, so the migration is index-for-index.
 *
 *   - `current_layer_[key] == -1`  <=> key never written
 *   - `0 <= current_layer_[key] <= claimed_layer_[key] < layer_capacity_` otherwise
 *   - `current_layer()` is the data mark, `claimed_layer()` the ownership mark
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
          current_layer_(count, -1),
          claimed_layer_(count, -1),
          count_(count),
          layer_capacity_(initial_layers > 0 ? initial_layers : 1)
    { }

    std::size_t count() const noexcept { return count_; }
    std::size_t layer_capacity() const noexcept { return layer_capacity_; }

    /// True once the key has been **written** at least once. Requesting a layer does not
    /// make a key exist -- see "Requested layers and written layers" above.
    bool exists(std::size_t key) const { check_key(key); return current_layer_[key] >= 0; }

    /// The layer this key's data currently stands at -- the last one written -- or -1 if it
    /// has never been written. This is what get(key) reads and what a rewind moves.
    int current_layer(std::size_t key) const { check_key(key); return current_layer_[key]; }

    /// The highest layer this key has been *claimed* at, or -1. Ownership, not data: an
    /// event asks this at initialization to learn the layer request_layer() just granted it,
    /// or the layer its neighbour owns. Reading a claimed-but-unwritten layer throws.
    int claimed_layer(std::size_t key) const { check_key(key); return claimed_layer_[key]; }

    /// Value at the key's current layer.
    /// \throws std::out_of_range if the key is unknown or has never been written.
    V get(std::size_t key) const
    {
        check_key(key);
        const int layer = current_layer_[key];
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
        if (layer >= layer_capacity_ || static_cast<int>(layer) > current_layer_[key]) {
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
        // Layers are filled bottom-up: a key may be written at any layer it owns, or at one
        // above (which claims that layer). Writing higher means the caller is using another
        // map's layer numbering. Validated against the *requested* mark, since an event
        // writes at the layer it was granted however many writes happened below it.
        if (static_cast<int>(layer) > claimed_layer_[key] + 1) {
            layered_array_detail::bad_layer("LayeredArray::set()", key, layer);
        }
        ensure_layer(layer);
        storage_[index(key, layer)] = value;
        current_layer_[key] = static_cast<int>(layer);
        if (static_cast<int>(layer) > claimed_layer_[key]) {
            claimed_layer_[key] = static_cast<int>(layer);
        }
    }

    /**
     * Write at the key's current layer, or at layer 0 if it has never been written.
     *
     * This is the assignment form of Enum_fast_memory_map::operator[], which returned a
     * reference and lazily marked an unwritten key as written at layer 0. Spelled out as a
     * named write so the layer being targeted is visible at the call site.
     */
    void set_current(std::size_t key, const V &value)
    {
        check_key(key);
        const int layer = current_layer_[key] < 0 ? 0 : current_layer_[key];
        ensure_layer(static_cast<std::size_t>(layer));
        storage_[index(key, static_cast<std::size_t>(layer))] = value;
        current_layer_[key] = layer;
        if (layer > claimed_layer_[key]) {
            claimed_layer_[key] = layer;
        }
    }

    /// Claim the next layer for this key. The layer's content is **unspecified until
    /// written**, and it is not readable until then: requesting is a promise to write, not a
    /// write. An event that requests a layer and hands off without writing it leaves the
    /// next reader of `layer - 1` with nothing, which get() now reports rather than serving
    /// a default.
    void request_layer(std::size_t key)
    {
        check_key(key);
        const int next = claimed_layer_[key] + 1;
        ensure_layer(static_cast<std::size_t>(next));
        claimed_layer_[key] = next;
    }

    /// Release the top layer for this key, dropping any value written there.
    void restore_layer(std::size_t key)
    {
        check_key(key);
        if (claimed_layer_[key] < 0) {
            layered_array_detail::unwritten_key("LayeredArray::restore_layer()", key);
        }
        --claimed_layer_[key];
        if (current_layer_[key] > claimed_layer_[key]) {
            current_layer_[key] = claimed_layer_[key];
        }
    }

    /// Rewind (or advance) a key to an already-written layer without reading it. This is the
    /// explicit half of Enum_fast_memory_map::at(key, layer), which conflated it with a read.
    void set_current_layer(std::size_t key, std::size_t layer)
    {
        check_key(key);
        if (static_cast<int>(layer) > current_layer_[key]) {
            layered_array_detail::bad_layer("LayeredArray::set_current_layer()", key, layer);
        }
        current_layer_[key] = static_cast<int>(layer);
    }

    /// One layer as a contiguous row over all keys.
    std::span<const V> layer(std::size_t l) const
    {
        if (l >= layer_capacity_) {
            layered_array_detail::bad_layer("LayeredArray::layer()", 0, l);
        }
        return std::span<const V>(storage_.data() + l * count_, count_);
    }

    /// Claimed layer of every key, in key order -- the same ownership answer as
    /// claimed_layer(). Replaces the raw int* out-parameter of
    /// Enum_fast_memory_map::get_all_current_memory_layer(). Callers snapshot this at
    /// initialization and feed it to multiply_all(), so it must name the layers this event
    /// owns, not what happens to have been written when the snapshot is taken.
    std::span<const int> claimed_layers() const noexcept
    {
        return std::span<const int>(claimed_layer_);
    }

    /**
     * Multiply `acc` by one value per key, each read from the layer named in `layers`.
     *
     * `layers` is a snapshot previously taken from this same map (see claimed_layers()), so
     * the entries are known-good and no per-key checking is done -- this sits in the pruning
     * bound computation, which runs at every node of the traversal.
     */
    void multiply_all(V &acc, std::span<const int> layers) const
    {
        //The snapshot must cover every key; a short one would read past its end below.
        assert(layers.size() >= count_);
        for (std::size_t k = 0; k != count_; ++k) {
            acc *= storage_[index(k, static_cast<std::size_t>(layers[k]))];
        }
    }

    /// Write `value` at layer 0 for every key, marking them all written.
    void init_first_layer(const V &value)
    {
        for (std::size_t k = 0; k != count_; ++k) {
            storage_[index(k, 0)] = value;
            current_layer_[k] = 0;
            if (claimed_layer_[k] < 0) {
                claimed_layer_[k] = 0;
            }
        }
    }

    /// Return every key that has been written to layer 0, keeping its value.
    void reset()
    {
        for (auto &l : current_layer_) {
            if (l > 0) {
                l = 0;
            }
        }
        for (auto &l : claimed_layer_) {
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
    std::vector<int> current_layer_;   ///< last layer written per key; -1 = never written
    std::vector<int> claimed_layer_;   ///< highest layer claimed per key; -1 = never claimed
    std::size_t      count_;
    std::size_t      layer_capacity_;
};
