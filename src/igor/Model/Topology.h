#pragma once

#include <igor/Model/Export.h>
#include <igor/Model/Navigator.h>

#include <igor/Model/Legacy/Rec_Event.h>
#include <igor/Core/Legacy/Typedef.h>

#include <nlohmann/json_fwd.hpp>

#include <vector>
#include <unordered_map>
#include <memory>
#include <string>

namespace igor::model {

class MODEL_EXPORT Topology 
{
public:
    using Adjacency_t = Navigator<legacy::Rec_Event>;

    // Core Graph Construction
    core::legacy::index_type addEvent(std::shared_ptr<legacy::Rec_Event> event);
    void addEdge(core::legacy::index_type parent_id, core::legacy::index_type child_id);
    
    // Efficient Access
    std::shared_ptr<legacy::Rec_Event> event(core::legacy::index_type id) const;
    std::shared_ptr<legacy::Rec_Event> event(const std::string& name) const;
    core::legacy::index_type eventId(const std::string& name) const;
    std::string eventName(core::legacy::index_type id) const;
    bool hasEvent(const std::string& name) const;

    const std::vector<core::legacy::index_type>& childrenIds(core::legacy::index_type id) const;
    const std::vector<core::legacy::index_type>& parentsIds(core::legacy::index_type id) const;

    // Range-based iteration helpers (typed via Adjacency_t)
    Adjacency_t parents (core::legacy::index_type id) const { return Adjacency_t(m_events, m_parents [id]); }
    Adjacency_t children(core::legacy::index_type id) const { return Adjacency_t(m_events, m_children[id]); }
  
    // Graph Inspection
    bool hasEdge(core::legacy::index_type parent_id, core::legacy::index_type child_id) const;
    std::vector<core::legacy::index_type> roots() const;
    std::vector<core::legacy::index_type> ancestors(core::legacy::index_type id) const;

    // Topological ordering (Kahn's algorithm) — roots first, leaves last.
    // Required by InferenceEngine and SamplingEngine iteration.
    std::vector<core::legacy::index_type> topologicalOrder() const;

    // Graph Modification
    void removeEdge(core::legacy::index_type parent_id, core::legacy::index_type child_id);
    void invertEdge(core::legacy::index_type parent_id, core::legacy::index_type child_id);

    std::size_t size() const { return m_events.size(); }
    auto begin() const { return m_events.begin(); }
    auto end()   const { return m_events.end();   }

    /**
     * \brief Left-to-right order of the sequence segments this model builds.
     *
     * Model data, and not derivable from the graph: the graph says which event conditions
     * which, not which segment sits left of which. A tandem-D model has seven entries where a
     * VDJ model has five, and losing the order means an exported model gets the standard VDJ
     * one inferred back, silently wrong for tandem D.
     *
     * Empty when the Topology was built by hand rather than from a model.
     */
    const std::vector<core::legacy::Seq_type_String>& seqTypeOrder() const { return m_seq_type_order; }
    void setSeqTypeOrder(std::vector<core::legacy::Seq_type_String> order) { m_seq_type_order = std::move(order); }

private:
    std::vector<std::shared_ptr<legacy::Rec_Event>>     m_events;
    std::vector<std::vector<core::legacy::index_type>>         m_children;
    std::vector<std::vector<core::legacy::index_type>>         m_parents;
    std::unordered_map<std::string, core::legacy::index_type>  m_name_to_id;
    std::vector<core::legacy::Seq_type_String>                 m_seq_type_order;
};

/**
 * \brief Build a Topology from a model document, through the event factory.
 *
 * Two phases, because a graph has forward references: every node is created first, then the
 * edges are wired from each node's parent list. The first phase is the factory's, so this
 * function has no dispatch of its own and knows no concrete event class.
 *
 * \throws std::runtime_error if the document is malformed, a type is not registered, or an
 *         edge names an event the document does not define.
 */
MODEL_EXPORT std::shared_ptr<Topology> topology_from_json(const nlohmann::json& doc);

} // namespace igor::model