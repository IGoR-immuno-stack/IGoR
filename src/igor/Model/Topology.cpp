#include <igor/Model/Topology.h>

#include <igor/Model/EventFactory.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cassert>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace igor::model {

index_type Topology::addEvent(std::shared_ptr<Rec_Event> event)
{
    index_type id = m_events.size();
    event->setUid(id);
    m_events.push_back(event);
    m_children.emplace_back();
    m_parents.emplace_back();
    m_name_to_id[event->get_nickname()] = id;
    return id;
}

void Topology::addEdge(index_type parent_id, index_type child_id)
{
    assert(parent_id < m_events.size() && child_id < m_events.size() && "Event index out of bounds");
    assert(parent_id != child_id && "Self-loops are not allowed");
    assert(!hasEdge(parent_id, child_id) && "Edge already exists");

    m_children[parent_id].push_back(child_id);
    std::sort(m_children[parent_id].begin(), m_children[parent_id].end());
    
    m_parents[child_id].push_back(parent_id);
    std::sort(m_parents[child_id].begin(), m_parents[child_id].end());
}

std::shared_ptr<Rec_Event> Topology::event(index_type id) const
{
    assert(id >= 0 && id < m_events.size() && "Event index out of bounds");
    return m_events[id];
}

std::shared_ptr<Rec_Event> Topology::event(const std::string& name) const
{
    auto it = m_name_to_id.find(name);
    if (it != m_name_to_id.end()) {
        return m_events[it->second];
    }
    return nullptr;
}

index_type Topology::eventId(const std::string& name) const
{
    auto it = m_name_to_id.find(name);
    if (it != m_name_to_id.end()) {
        return it->second;
    }
    return -1;
}

std::string Topology::eventName(index_type id) const
{
    assert(id >= 0 && id < m_events.size() && "Event index out of bounds");
    return m_events[id]->get_nickname();
}

bool Topology::hasEvent(const std::string& name) const
{
    return m_name_to_id.count(name) > 0;
}

const std::vector<index_type>& Topology::childrenIds(index_type id) const
{
    assert(id >= 0 && id < m_events.size() && "Event index out of bounds");
    return m_children[id];
}

const std::vector<index_type>& Topology::parentsIds(index_type id) const
{
    assert(id >= 0 && id < m_events.size() && "Event index out of bounds");
    return m_parents[id];
}

bool Topology::hasEdge(index_type parent_id, index_type child_id) const
{
    assert(parent_id >= 0 && parent_id < m_events.size() && child_id >= 0 && child_id < m_events.size() && "Event index out of bounds");
    const auto& children = m_children[parent_id];
    return std::find(children.begin(), children.end(), child_id) != children.end();
}

std::vector<index_type> Topology::roots(void) const
{
    std::vector<index_type> roots;
    for (index_type i = 0; i < m_events.size(); ++i) {
        if (m_parents[i].empty()) {
            roots.push_back(i);
        }
    }
    return roots;
}

std::vector<index_type> Topology::ancestors(index_type id) const
{
    assert(id >= 0 && id < m_events.size() && "Event index out of bounds");
    std::vector<index_type> ancestors;
    std::vector<index_type> stack = m_parents[id];
    
    while (!stack.empty()) {
        index_type current = stack.back();
        stack.pop_back();
        
        // Avoid duplicates and self-loops
        if (std::find(ancestors.begin(), ancestors.end(), current) == ancestors.end()) {
            ancestors.push_back(current);
            // Add parents of current
            for (index_type parent : m_parents[current]) {
                stack.push_back(parent);
            }
        }
    }
    return ancestors;
}

void Topology::removeEdge(index_type parent_id, index_type child_id)
{
    assert(parent_id >= 0 && parent_id < m_events.size() && child_id >= 0 && child_id < m_events.size() && "Event index out of bounds");
    // Remove child from parent's children list
    auto& children = m_children[parent_id];
    children.erase(std::remove(children.begin(), children.end(), child_id), children.end());
    
    // Remove parent from child's parents list
    auto& parents = m_parents[child_id];
    parents.erase(std::remove(parents.begin(), parents.end(), parent_id), parents.end());
}

void Topology::invertEdge(index_type parent_id, index_type child_id)
{
    assert(parent_id >= 0 && parent_id < m_events.size() && child_id >= 0 && child_id < m_events.size() && "Event index out of bounds");
    if (hasEdge(parent_id, child_id)) {
        removeEdge(parent_id, child_id);
        addEdge(child_id, parent_id);
    }
}

std::vector<index_type> Topology::topologicalOrder() const
{
    const std::size_t n = m_events.size();

    // In-degree count for each node
    std::vector<std::size_t> in_degree(n, 0);
    for (index_type i = 0; i < static_cast<index_type>(n); ++i)
        in_degree[i] = m_parents[i].size();

    // Tie-break among ready nodes must match Core: Event_comparator orders by DESCENDING
    // priority ("event_p1->get_priority() > event_p2->get_priority()", Rec_Event.h), and
    // Model_Parms::get_model_queue() relies on that order both for iteration and for the
    // layout of the marginal array. The heap below is a min-heap, so the key holds the
    // NEGATED priority: the highest priority event pops first. On equal priorities the
    // smaller uid pops first, matching the stable sort Core applies to its event list.
    using Entry = std::pair<int, index_type>;   // { -priority, uid }
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> ready;

    for (index_type i = 0; i < static_cast<index_type>(n); ++i)
        if (in_degree[i] == 0)
            ready.push({ -m_events[i]->get_priority(), i });

    std::vector<index_type> order;
    order.reserve(n);

    while (!ready.empty()) {
        const index_type node = ready.top().second;
        ready.pop();
        order.push_back(node);
        for (index_type child : m_children[node]) {
            if (--in_degree[child] == 0)
                ready.push({ -m_events[child]->get_priority(), child });
        }
    }

    assert(order.size() == n && "Topology contains a cycle");
    return order;
}

std::shared_ptr<Topology> topology_from_json(const nlohmann::json& doc)
{
    const auto events = doc.find("events");
    if (events == doc.end() || !events->is_array())
        throw std::runtime_error("topology_from_json: the document has no \"events\" array");

    auto topo = std::make_shared<Topology>();

    const auto order = doc.find("seq_type_order");
    if (order != doc.end())
        topo->setSeqTypeOrder(order->get<std::vector<Seq_type_String>>());

    // Phase 1: every node, built complete by the factory. No dispatch here, and no concrete
    // event class named: adding an event type touches the factory's registrars, not this loop.
    for (const auto& node : *events)
        topo->addEvent(event_factory::create(node));

    // Phase 2: the edges, which is why this needs two phases at all. A parent may be named
    // before it is defined, so the wiring waits until every node exists.
    for (const auto& node : *events) {
        const auto child_nickname = node.at("nickname").get<std::string>();
        const index_type child_id = topo->eventId(child_nickname);
        if (child_id < 0)
            throw std::runtime_error("topology_from_json: event \"" + child_nickname
                                     + "\" was not registered by its own nickname");

        const auto parents = node.find("parents");
        if (parents == node.end())
            continue;
        for (const auto& parent : *parents) {
            const auto parent_nickname = parent.get<std::string>();
            const index_type parent_id = topo->eventId(parent_nickname);
            if (parent_id < 0)
                throw std::runtime_error("topology_from_json: event \"" + child_nickname
                                         + "\" names an unknown parent \"" + parent_nickname
                                         + "\"");
            if (!topo->hasEdge(parent_id, child_id))
                topo->addEdge(parent_id, child_id);
        }
    }

    return topo;
}

} // namespace igor::model