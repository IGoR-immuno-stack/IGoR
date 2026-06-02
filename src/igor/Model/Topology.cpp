#include <igor/Model/Topology.h>

#include <algorithm>
#include <cassert>
#include <queue>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <igor/Core/Utils.h>
#include <igor/Model/EventFactory.h>

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

    // Min-heap: (priority, index) — lowest priority value = highest priority event first,
    // matching Model_Parms::get_model_queue() which uses Event_comparator (sort by priority).
    using Entry = std::pair<int, index_type>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> ready;

    for (index_type i = 0; i < static_cast<index_type>(n); ++i)
        if (in_degree[i] == 0)
            ready.push({ m_events[i]->get_priority(), i });

    std::vector<index_type> order;
    order.reserve(n);

    while (!ready.empty()) {
        auto [prio, node] = ready.top(); ready.pop();
        order.push_back(node);
        for (index_type child : m_children[node]) {
            if (--in_degree[child] == 0)
                ready.push({ m_events[child]->get_priority(), child });
        }
    }

    assert(order.size() == n && "Topology contains a cycle");
    return order;
}

std::shared_ptr<Topology> read_topology(const std::string& filename)
{
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("File not found: \"" + filename + "\"");
    }

    auto topo = std::make_shared<Topology>();

    // Map from full event name -> topology id, needed for edge resolution
    std::unordered_map<std::string, index_type> name_to_id;

    std::string line;
    std::getline(infile, line);

    if (line != "@Event_list") {
        throw std::runtime_error("Unknown format for model_parms file: " + filename);
    }

    // ── Parse @Event_list section ──────────────────────────────────────
    std::getline(infile, line);
    while (line.size() > 0 && line[0] == '#') {
        // Parse header: #EventType;GeneClass;SeqSide;Priority[;Nickname]
        std::size_t sc1 = line.find(';', 0);
        std::string event_type_str = line.substr(1, sc1 - 1);

        std::size_t sc2 = line.find(';', sc1 + 1);
        std::string gene_class_str = line.substr(sc1 + 1, sc2 - sc1 - 1);
        Gene_class gene_class = str2GeneClass(gene_class_str);

        std::size_t sc3 = line.find(';', sc2 + 1);
        std::string seq_side_str = line.substr(sc2 + 1, sc3 - sc2 - 1);
        Seq_side seq_side = str2SeqSide(seq_side_str);

        std::size_t sc4 = line.find(';', sc3 + 1);
        int priority;
        std::string nickname;
        Event_type event_type;
        if (sc4 == std::string::npos) {
            priority = std::stoi(line.substr(sc3 + 1));
        } else {
            priority = std::stoi(line.substr(sc3 + 1, sc4 - sc3 - 1));
            nickname = line.substr(sc4 + 1);
        }

        if (event_type_str == "Insertion") {
            event_type = Insertion_t;
        } else if (event_type_str == "Deletion") {
            event_type = Deletion_t;
        } else if (event_type_str == "GeneChoice") {
            event_type = GeneChoice_t;
        } else if (event_type_str == "DinucMarkov") {
            event_type = Dinuclmarkov_t;
        } else {
            throw std::runtime_error(event_type_str + " event is not implemented (thrown by read_topology)");
        }

        // Phase 1: create via factory
        auto ev = event_factory::create(event_type);

        // Phase 2: configure
        ev->set_event_class(gene_class);
        ev->set_event_side(seq_side);
        ev->set_priority(priority);
        ev->set_nickname(nickname);

        // Parse realizations
        std::getline(infile, line);
        if (event_type == GeneChoice_t) {
            while (line.size() > 0 && line[0] == '%') {
                std::size_t si1 = line.find(';', 0);
                std::string name = line.substr(1, si1 - 1);
                std::size_t si2 = line.find(';', si1 + 1);
                std::string value_str = line.substr(si1 + 1, si2 - si1 - 1);
                int index = std::stoi(line.substr(si2 + 1));
                ev->add_realization(
                    Event_realization(name, INT16_MAX, value_str, nt2int(value_str), index));
                std::getline(infile, line);
            }
        } else if (event_type == Insertion_t || event_type == Deletion_t) {
            while (line.size() > 0 && line[0] == '%') {
                std::size_t si = line.find(';', 0);
                int value_int = std::stoi(line.substr(1, si));
                int index = std::stoi(line.substr(si + 1));
                ev->add_realization(
                    Event_realization(std::to_string(value_int), value_int, "", Int_Str(), index));
                std::getline(infile, line);
            }
        } else {
            // DinucMarkov: realizations are self-initialized, skip lines.
            // Also apply the legacy conversion for seq_side: in legacy files the
            // side is always Undefined_side, but the direction is implicit in the
            // gene_class (VD_genes/VJ_genes → Three_prime, DJ_genes → Five_prime).
            if (seq_side == Undefined_side) {
                if (gene_class == VD_genes || gene_class == VJ_genes || gene_class == VDJ_genes) {
                    ev->set_event_side(Three_prime);
                } else if (gene_class == DJ_genes) {
                    ev->set_event_side(Five_prime);
                }
            }
            while (line.size() > 0 && line[0] == '%') {
                std::getline(infile, line);
            }
        }

        index_type id = topo->addEvent(ev);
        name_to_id[ev->get_name()] = id;
    }

    // ── Parse @Edges section ───────────────────────────────────────────
    if (line != "@Edges") {
        throw std::runtime_error("Unknown format for model file: " + filename);
    }

    std::getline(infile, line);
    while (line.size() > 0 && line[0] == '%') {
        std::size_t si = line.find(';', 0);
        std::string parent_name = line.substr(1, si - 1);
        std::string child_name = line.substr(si + 1);

        auto pit = name_to_id.find(parent_name);
        auto cit = name_to_id.find(child_name);
        if (pit == name_to_id.end()) {
            throw std::runtime_error("read_topology: Unknown parent event \"" + parent_name + "\"");
        }
        if (cit == name_to_id.end()) {
            throw std::runtime_error("read_topology: Unknown child event \"" + child_name + "\"");
        }

        if (!topo->hasEdge(pit->second, cit->second)) {
            topo->addEdge(pit->second, cit->second);
        }

        std::getline(infile, line);
    }

    // @ErrorRate section is ignored

    return topo;
}

} // namespace igor::model