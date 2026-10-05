#pragma once

// Model: the probabilistic model and its persistence. The engines that compute with it live in
// <igor/Inference.h> and <igor/Generation.h>.
#include <igor/Model/EventFactory.h>
#include <igor/Model/Topology.h>
#include <igor/Model/Navigator.h>
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/Scenario.h>

#include <igor/Model/LegacyBridge.h>

// Scenario geometry and pruning bounds, promoted out of Legacy/ in step 1c.
#include <igor/Model/SpanProfile.h>
#include <igor/Model/SafetyMatrix.h>
#include <igor/Model/UnfilledSegmentLengths.h>
#include <igor/Model/JunctionGeometry.h>
#include <igor/Model/BoundTightness.h>
