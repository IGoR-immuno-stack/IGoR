#pragma once

// Forward declarations of the Model types the engines (Inference, Generation) name in their
// headers, plus the legacy namespaces they qualify into. Including this header is enough to
// write `using igor::model::RecombinationModel;` or `model::legacy::Rec_Event` before the
// full definitions are visible.

namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy {}

namespace igor::model {

template <typename T> class RecombinationModel;
template <typename NodeType, typename PtrType> class Navigator;
class Topology;
struct SampledScenario;
struct SampledEvent;

} // namespace igor::model
