#include <igor/Model/Topology.h>

#include <algorithm>
#include <cassert>
#include <queue>
#include <stdexcept>
#include <fstream>
#include <igor/Core/Utils.h>
#include <igor/Core/SequenceTypes.h>
#include <igor/Model/EventFactory.h>

namespace {
void register_tandem_d_types_if_needed(const std::string &nickname)
{
    auto &registry = SequenceTypeRegistry::get_instance();
    if (nickname.find("D1") != std::string::npos || nickname.find("D2") != std::string::npos) {
        auto d1_gene = registry.register_type("D1_gene_seq");
        registry.register_type("D1_gene");
        auto d2_gene = registry.register_type("D2_gene_seq");
        registry.register_type("D2_gene");

        auto vd1_ins = registry.register_type("VD1_ins_seq");
        registry.register_type("VD1_ins");
        auto d1d2_ins = registry.register_type("D1D2_ins_seq");
        registry.register_type("D1D2_ins");
        auto d2j_ins = registry.register_type("D2J_ins_seq");
        registry.register_type("D2J_ins");

        registry.register_connection(SequenceTypeRegistry::V_GENE_SEQ, d1_gene, vd1_ins);
        registry.register_connection(d1_gene, d2_gene, d1d2_ins);
        registry.register_connection(d2_gene, SequenceTypeRegistry::J_GENE_SEQ, d2j_ins);
    }
}
}

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

std::shared_ptr<Topology> readTopology(const std::string& filename)
{
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("File not found : \"" + filename + "\"");
    }

    auto topo = std::make_shared<Topology>();

    std::string line_str;
    getline(infile, line_str);
    if (line_str == "@Event_list") {
        getline(infile, line_str);
        while (!line_str.empty() && line_str[0] == '#') {
            size_t semi = line_str.find(";");
            std::string event_type_str = line_str.substr(1, semi - 1);
            
            size_t next_semi = line_str.find(";", semi + 1);
            std::string event_class_str = line_str.substr(semi + 1, next_semi - semi - 1);
            int event_class;
            try {
                event_class = ::str2GeneClass(event_class_str);
            } catch (...) {
                throw std::runtime_error("Unknown Gene_class\"" + event_class_str + "\" in model file: \"" + filename + "\"");
            }

            semi = next_semi;
            next_semi = line_str.find(";", semi + 1);
            std::string event_side_str = line_str.substr(semi + 1, next_semi - semi - 1);
            ::Seq_side event_side;
            try {
                event_side = ::str2SeqSide(event_side_str);
            } catch (...) {
                throw std::runtime_error("Unknown Seq_side\"" + event_side_str + "\" in file: \"" + filename + "\"");
            }

            semi = next_semi;
            next_semi = line_str.find(";", semi + 1);
            int priority;
            std::string nickname;
            if (next_semi == std::string::npos) {
                priority = std::stoi(line_str.substr(semi + 1));
            } else {
                priority = std::stoi(line_str.substr(semi + 1, next_semi - semi - 1));
                nickname = line_str.substr(next_semi + 1);
            }

            if (!nickname.empty()) {
                ::register_tandem_d_types_if_needed(nickname);
            }

            ::Event_type type;
            if (event_type_str == "Insertion") type = ::Insertion_t;
            else if (event_type_str == "Deletion") type = ::Deletion_t;
            else if (event_type_str == "GeneChoice") type = ::GeneChoice_t;
            else if (event_type_str == "DinucMarkov") type = ::Dinuclmarkov_t;
            else throw std::runtime_error(event_type_str + " event is not implemented");

            auto ev = igor::model::event_factory::create(type);
            ev->set_event_class(event_class);
            ev->set_event_side(event_side);
            ev->set_priority(priority);
            ev->set_nickname(nickname);

            if (type == ::Insertion_t || type == ::Deletion_t) {
                getline(infile, line_str);
                while (!line_str.empty() && line_str[0] == '%') {
                    semi = line_str.find(";");
                    int val = std::stoi(line_str.substr(1, semi - 1));
                    int idx = std::stoi(line_str.substr(semi + 1));
                    ev->add_realization(::Event_realization(
                        std::to_string(val), val, "", ::Int_Str(), idx));
                    getline(infile, line_str);
                }
            } else if (type == ::GeneChoice_t) {
                getline(infile, line_str);
                while (!line_str.empty() && line_str[0] == '%') {
                    semi = line_str.find(";");
                    std::string name = line_str.substr(1, semi - 1);
                    next_semi = line_str.find(";", semi + 1);
                    std::string seq = line_str.substr(semi + 1, next_semi - semi - 1);
                    int idx = std::stoi(line_str.substr(next_semi + 1));
                    ev->add_realization(::Event_realization(
                        name, INT16_MAX, seq, ::nt2int(seq), idx));
                    getline(infile, line_str);
                }
            } else if (type == ::Dinuclmarkov_t) {
                getline(infile, line_str);
                while (!line_str.empty() && line_str[0] == '%') {
                    getline(infile, line_str); // Just skip the inner representations, handled in read_model_marginals
                }
            }
            topo->addEvent(ev);
        }
    } else {
        throw std::runtime_error("Unknown format for model file");
    }

    if (line_str == "@Edges") {
        getline(infile, line_str);
        while (!line_str.empty() && line_str[0] == '%') {
            size_t semi = line_str.find(";");
            std::string parent_name = line_str.substr(1, semi - 1);
            std::string child_name  = line_str.substr(semi + 1);

            size_t size_pos = parent_name.find("_size");
            if (size_pos != std::string::npos) parent_name = parent_name.substr(0, size_pos);
            size_pos = child_name.find("_size");
            if (size_pos != std::string::npos) child_name = child_name.substr(0, size_pos);

            index_type parent_id = -1;
            index_type child_id  = -1;
            for (const auto& ev : *topo) {
                if (ev->get_name() == parent_name || ev->get_nickname() == parent_name) {
                    parent_id = topo->eventId(ev->get_nickname());
                }
                if (ev->get_name() == child_name || ev->get_nickname() == child_name) {
                    child_id = topo->eventId(ev->get_nickname());
                }
            }

            if (parent_id >= 0 && child_id >= 0) {
                topo->addEdge(parent_id, child_id);
            }
            getline(infile, line_str);
        }
    }

    // Ignore @ErrorRate if it exists

    return topo;
}

} // namespace igor::model