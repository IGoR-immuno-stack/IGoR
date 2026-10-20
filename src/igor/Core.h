#pragma once

// Core: the vocabulary every layer shares. The legacy helpers still under Core/Legacy
// (Utils.h) are not part of this umbrella; include them by their Legacy/ path.
#include <igor/Core/Platform.h>
#include <igor/Core/Types.h>
#include <igor/Core/IntStr.h>
#include <igor/Core/GeneticCode.h>
#include <igor/Core/AlignmentData.h>

#include <igor/Core/LayeredArray.h>
#include <igor/Core/SeqTypeRegistry.h>
#include <igor/Core/DynamicSequenceMap.h>
#include <igor/Core/SeqOffsetsMap.h>
#include <igor/Core/SegmentSpan.h>
