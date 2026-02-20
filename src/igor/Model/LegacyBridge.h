#pragma once

#include <igor/Core/Model_marginals.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Model/InferenceEngine.h>
#include <igor/Model/MarginalHandler.h>

#include <memory>
#include <vector>

namespace igor::model {

/// Build EventDescriptor vector from Model_Parms for InferenceEngine construction
std::vector<model::EventDescriptor> extract_event_descriptors(const Model_Parms& parms);

/// Import marginal probabilities from legacy Model_marginals into InferenceEngine
template <typename T = long double>
void import_from_legacy(model::InferenceEngine<T>& engine,
                       const Model_marginals& marginals,
                       const Model_Parms& parms);

/// Export marginal probabilities from InferenceEngine to legacy Model_marginals
template <typename T = long double>
void export_to_legacy(const model::InferenceEngine<T>& engine,
                     Model_marginals& marginals,
                     const Model_Parms& parms);

class Topology;

/// Import the legacy Model_Parms structure into a modern Topology graph
std::shared_ptr<Topology> import_from_legacy(const Model_Parms& legacy_model);

/// Export a modern Topology graph back to a legacy Model_Parms
std::shared_ptr<Model_Parms> export_to_legacy(const Topology& topology);

} // namespace igor::model

#include <igor/Model/LegacyBridge.tpp>
