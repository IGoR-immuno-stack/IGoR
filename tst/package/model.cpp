// Links against the Model layer: an empty topology has no event.
#include <igor/Model/Topology.h>

int main()
{
    igor::model::Topology topology;
    return static_cast<int>(topology.size());
}
