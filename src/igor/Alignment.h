#pragma once

// Alignment: algorithms on sequences. The whole layer is still under Alignment/Legacy, in
// igor::alignment::legacy; this umbrella is its public face until it is promoted.
// AlignerInternal.h is for the whitebox tests only and stays out.
#include <igor/Alignment/Legacy/Aligner.h>
#include <igor/Alignment/Legacy/CDR3SeqData.h>
#include <igor/Alignment/Legacy/ExtractFeatures.h>
#include <igor/Alignment/Legacy/JournaledQuery.h>
