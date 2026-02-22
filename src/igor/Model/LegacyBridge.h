#pragma once

#include <igor/Core/Model_marginals.h>
#include <igor/Core/Model_Parms.h>

#include <memory>
#include <vector>

namespace igor::model {

class Topology;
template <typename T> class RecombinationModel;

/// Import the legacy Model_Parms structure into a modern Topology graph
std::shared_ptr<Topology> import_from_legacy(const Model_Parms& legacy_model);

/// Export a modern Topology graph back to a legacy Model_Parms
std::shared_ptr<Model_Parms> export_to_legacy(const Topology& topology);

/// Import marginal probabilities from legacy Model_marginals into a RecombinationModel.
/// The model must already be initialised with the correct Topology (which determines
/// tensor shapes); only the probability values are filled.
template <typename T = double>
void import_from_legacy(RecombinationModel<T>& model,
                        const Model_marginals& marginals);

} // namespace igor::model

#include <igor/Model/LegacyBridge.tpp>
