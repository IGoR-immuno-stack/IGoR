#pragma once

#include <igor/Core/Typedef.h>

#include <concepts>
#include <vector>
#include <memory>
#include <iterator>

namespace igor::model {

// ─── HasUid concept ──────────────────────────────────────────────────────────
//
// Any type used as a Topology / engine node must expose:
//   - uid()          → igor::index_type   (current assigned index)
//   - setUid(id)     → void               (called by Topology::addEvent / engine::add_handler)
//
// Rec_Event satisfies this already.
// InferenceHandler and SamplingHandler will satisfy it once uid_ is added.

template <typename T>
concept HasUid = requires(T& t, const T& ct, igor::index_type id) {
    { ct.uid()      } -> std::same_as<igor::index_type>;
    { t.setUid(id)  } -> std::same_as<void>;
};

// ─── Navigator<NodeType> ─────────────────────────────────────────────────────
//
// A lightweight, non-owning view over a subset of nodes in a parallel node
// vector. Constructed by Topology (for Rec_Event nodes) or by InferenceEngine /
// SamplingEngine (for their respective handler vectors).
//
// Usage:
//   for (auto& handler : engine.parents(i))   // Navigator<InferenceHandler<T>>
//       use(handler->uid(), ...);
//
//   for (auto& event : topology.parents(i))   // Navigator<Rec_Event>
//       use(event->get_nickname(), ...);
//
// NodeType must satisfy HasUid.

template <HasUid NodeType>
class Navigator {
public:
    // ── Iterator ─────────────────────────────────────────────────────────────

    class Iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = std::shared_ptr<NodeType>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const std::shared_ptr<NodeType>*;
        using reference         = std::shared_ptr<NodeType>;

        Iterator(const std::vector<std::shared_ptr<NodeType>>& nodes,
                 std::vector<igor::index_type>::const_iterator  it)
            : m_nodes(nodes), m_it(it) {}

        reference operator*()  const { return m_nodes[*m_it]; }
        reference operator[](difference_type n) const { return m_nodes[*(m_it + n)]; }
        pointer   operator->() const { return &m_nodes[*m_it]; }

        // Increment / Decrement
        Iterator& operator++()    { ++m_it; return *this; }
        Iterator  operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        Iterator& operator--()    { --m_it; return *this; }
        Iterator  operator--(int) { Iterator tmp = *this; --(*this); return tmp; }

        // Arithmetic
        Iterator& operator+=(difference_type n) { m_it += n; return *this; }
        Iterator& operator-=(difference_type n) { m_it -= n; return *this; }

        friend Iterator operator+(Iterator it, difference_type n) { return it += n; }
        friend Iterator operator+(difference_type n, Iterator it) { return it += n; }
        friend Iterator operator-(Iterator it, difference_type n) { return it -= n; }
        friend difference_type operator-(const Iterator& a, const Iterator& b) {
            return a.m_it - b.m_it;
        }

        // Comparison
        bool operator==(const Iterator& o) const { return m_it == o.m_it; }
        bool operator!=(const Iterator& o) const { return m_it != o.m_it; }
        bool operator< (const Iterator& o) const { return m_it <  o.m_it; }
        bool operator> (const Iterator& o) const { return m_it >  o.m_it; }
        bool operator<=(const Iterator& o) const { return m_it <= o.m_it; }
        bool operator>=(const Iterator& o) const { return m_it >= o.m_it; }

    private:
        const std::vector<std::shared_ptr<NodeType>>& m_nodes;
        std::vector<igor::index_type>::const_iterator  m_it;
    };

    // ── Construction ──────────────────────────────────────────────────────────

    /// View over a subset of `nodes` identified by `indices`.
    /// Both references must outlive this Navigator object.
    Navigator(const std::vector<std::shared_ptr<NodeType>>& nodes,
              const std::vector<igor::index_type>& indices)
              : m_nodes(nodes), m_indices(indices) {}

    // ── Range interface ───────────────────────────────────────────────────────

    Iterator begin(void) const { return Iterator(m_nodes, m_indices.begin()); }
    Iterator end(void)   const { return Iterator(m_nodes, m_indices.end());   }

    std::size_t size(void) const { return m_indices.size();  }
    bool        empty(void) const { return m_indices.empty(); }

    std::shared_ptr<NodeType> operator[](std::size_t n) const {
        return m_nodes[m_indices[n]];
    }

    /// Direct access to the raw index list (for engine parent resolution).
    const std::vector<igor::index_type>& indices(void) const { return m_indices; }

private:
    const std::vector<std::shared_ptr<NodeType>>& m_nodes;
    const std::vector<igor::index_type>&          m_indices;
};

} // namespace igor::model
