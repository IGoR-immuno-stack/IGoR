/**
 * @file GeneticCode.h
 * @brief Standard genetic code utilities: codon indexing, CodonMask sets,
 *        translation, and OLGA-style AA motif parsing.
 *
 * Extracted from feature/AA_PGEN so that these utilities can live on develop
 * independently of the AA-level Pgen work, which is being re-derived on top of
 * the Rec_Event/iterate refactor.
 *
 * These declarations deliberately remain in `namespace EventUtils` (namespaces
 * are open) so that existing call sites compile unchanged, while living in
 * their own translation unit: EventUtils.h is under active rework for the
 * SeqTypeRegistry migration, and the genetic code has no reason to share its
 * merge surface.
 *
 * Nothing here depends on Rec_Event, Gene_class, or any event machinery.
 */

#pragma once

#include <igor/Core/GeneticCode.h>
#include <igor/Core/Legacy/IntStr.h>

// GeneticCode.h was promoted out of Legacy/ (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).
// This stub keeps the legacy include path and the legacy names for the code that has not
// been promoted yet; it goes when its last consumer switches.
namespace igor::core::legacy {

namespace genetic_code = igor::core::genetic_code;

} // namespace igor::core::legacy
