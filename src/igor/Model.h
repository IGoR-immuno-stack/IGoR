#pragma once

// Inference and marginal handlers
#include <igor/Model/MarginalHandler.h>
#include <igor/Model/CategoricalHandler.h>
#include <igor/Model/MarkovHandler.h>
#include <igor/Model/InferenceEngine.h>

// Sampling and generation
#include <igor/Model/SamplingHandler.h>
#include <igor/Model/CategoricalSamplingHandler.h>
#include <igor/Model/MarkovSamplingHandler.h>
#include <igor/Model/SamplingEngine.h>

// Architecture and graph topology
#include <igor/Model/EventFactory.h>
#include <igor/Model/SamplingHandlerFactory.h>
#include <igor/Model/Topology.h>
#include <igor/Model/Navigator.h>
#include <igor/Model/Scenario.h>

// Compatibility layer
#include <igor/Model/LegacyBridge.h>
