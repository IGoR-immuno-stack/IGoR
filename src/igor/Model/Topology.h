#pragma once

#include <igor/Model/Export.h>
#include <igor/Model/Navigator.h>

#include <igor/Core/Rec_Event.h>
#include <igor/Core/Typedef.h>

#include <vector>
#include <unordered_map>
#include <memory>

namespace igor::model {

class MODEL_EXPORT Topology 
{
public:
    using Adjacency_t = Navigator<Rec_Event>;

    // Core Graph Construction
    index_type addEvent(std::shared_ptr<Rec_Event> event);
    void addEdge(index_type parent_id, index_type child_id);
    
    // Efficient Access
    std::shared_ptr<Rec_Event> event(index_type id) const;
    std::shared_ptr<Rec_Event> event(const std::string& name) const;
    index_type eventId(const std::string& name) const;
    std::string eventName(index_type id) const;
    bool hasEvent(const std::string& name) const;

    const std::vector<index_type>& childrenIds(index_type id) const;
    const std::vector<index_type>& parentsIds(index_type id) const;

    // Range-based iteration helpers (typed via Adjacency_t)
    Adjacency_t parents (index_type id) const { return Adjacency_t(m_events, m_parents [id]); }
    Adjacency_t children(index_type id) const { return Adjacency_t(m_events, m_children[id]); }
  
    // Graph Inspection
    bool hasEdge(index_type parent_id, index_type child_id) const;
    std::vector<index_type> roots() const;
    std::vector<index_type> ancestors(index_type id) const;

    // Topological ordering (Kahn's algorithm) — roots first, leaves last.
    // Required by InferenceEngine and SamplingEngine iteration.
    std::vector<index_type> topologicalOrder() const;

    // Graph Modification
    void removeEdge(index_type parent_id, index_type child_id);
    void invertEdge(index_type parent_id, index_type child_id);

    std::size_t size() const { return m_events.size(); }
    auto begin() const { return m_events.begin(); }
    auto end()   const { return m_events.end();   }
    
private:
    std::vector<std::shared_ptr<Rec_Event>>     m_events;
    std::vector<std::vector<index_type>>         m_children;
    std::vector<std::vector<index_type>>         m_parents;
    std::unordered_map<std::string, index_type>  m_name_to_id;
};

MODEL_EXPORT std::shared_ptr<Topology> read_topology(const std::string& filename);

} // namespace igor::model