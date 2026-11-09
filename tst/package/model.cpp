// Links against the Model layer: an empty topology has no event.
#include <igor/Model.h>
#include <igor/Math.h>

int main()
{
    igor::model::Topology topology;
    return static_cast<int>(topology.size());
}
