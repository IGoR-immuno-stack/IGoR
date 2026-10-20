#pragma once

// Forward declarations of the Model types the engines (Inference, Generation) name in their
// headers. Including this header is enough to write `using igor::model::RecombinationModel;`
// before the full definition is visible.

namespace igor::model {

template <typename T> class RecombinationModel;
template <typename NodeType, typename PtrType> class Navigator;
class Topology;
struct SampledScenario;
struct SampledEvent;

} // namespace igor::model
