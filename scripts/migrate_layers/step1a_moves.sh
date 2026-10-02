#!/usr/bin/env bash
# Step 1a of doc/LAYER_REFACTORING_PROPOSAL.md: every ex-Core file goes to <Layer>/Legacy/,
# the tk engines leave Model for Inference and Generation. Pure `git mv`, no content change.
#
# Run from the repository root. Re-runnable: a source that no longer exists is skipped, so the
# script can be applied to another branch that already carries part of the move.
#
# Order matters for history: src/igor/Core is renamed to src/igor/Model/Legacy first, because
# the Rec_Event cluster is the bulk of it and is where the recent history lives; the other
# files are then moved out of Model/Legacy. git's rename detection follows both steps.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

mv_if_present() {
  local src="$1" dst="$2"
  if [ -e "$src" ]; then
    mkdir -p "$(dirname "$dst")"
    git mv "$src" "$dst"
  fi
}

# Stems of the ex-Core files, by target layer. A stem covers .h, .cpp and .tpp.
CORE_STEMS="CoreEnums StdTypedefs Typedef IntStr GeneticCode SeqTypeRegistry SegmentSpan
  LayeredArray DynamicSequenceMap SeqOffsetsMap Utils"
ALIGNMENT_STEMS="Aligner AlignerInternal ExtractFeatures CDR3SeqData JournaledQuery"
INFERENCE_STEMS="GenModel Pgencounter Bestscenarioscounter Coverageerrcounter Errorscounter"
GENERATION_STEMS="FastGenerator FastSampling"
# Everything else in the old Core stays in Model/Legacy: Rec_Event and the four events,
# EventUtils, gene_to_seqtype_migr, JsonDetail, ModelJson, the error rates, Model_Parms,
# Model_marginals, Counter, the five contexts, Scenario, JunctionGeometry, SpanProfile,
# SafetyMatrix, UnfilledSegmentLengths, BoundTightness.

# tk engine stems leaving src/igor/Model.
INFERENCE_ENGINE_STEMS="InferenceEngine InferenceHandler CategoricalInferenceHandler
  MarkovInferenceHandler InferenceHandlerFactory"
GENERATION_ENGINE_STEMS="SamplingEngine SamplingHandler CategoricalSamplingHandler
  MarkovSamplingHandler SamplingHandlerFactory"

move_stems() {
  local from="$1" to="$2" stems="$3" stem ext
  for stem in $stems; do
    for ext in h cpp tpp; do
      mv_if_present "$from/$stem.$ext" "$to/$stem.$ext"
    done
  done
}

# 1. The old Core becomes Model/Legacy.
if [ -d src/igor/Core ] && [ ! -d src/igor/Model/Legacy ]; then
  git mv src/igor/Core src/igor/Model/Legacy
fi

# 2. Core keeps its build files and its foundation headers.
mv_if_present src/igor/Model/Legacy/CMakeLists.txt src/igor/Core/CMakeLists.txt
mv_if_present src/igor/Model/Legacy/Config.h.in    src/igor/Core/Config.h.in
move_stems src/igor/Model/Legacy src/igor/Core/Legacy       "$CORE_STEMS"
move_stems src/igor/Model/Legacy src/igor/Alignment/Legacy  "$ALIGNMENT_STEMS"
move_stems src/igor/Model/Legacy src/igor/Inference/Legacy  "$INFERENCE_STEMS"
move_stems src/igor/Model/Legacy src/igor/Generation/Legacy "$GENERATION_STEMS"

# 3. The tk engines leave Model.
move_stems src/igor/Model src/igor/Inference  "$INFERENCE_ENGINE_STEMS"
move_stems src/igor/Model src/igor/Generation "$GENERATION_ENGINE_STEMS"

# 4. Tests follow their layer (see tst/igor/<Layer>/CMakeLists.txt for the executables).
#    Rule: the layer of the highest header a test includes.
move_tests() {
  local to="$1" names="$2" name
  for name in $names; do
    mv_if_present "tst/igor/Core/$name" "$to/$name"
  done
}
move_tests tst/igor/Core/Legacy "test_dynamic_sequence_map.cpp test_layered_array.cpp
  test_matrix.cpp test_sanity.cpp"
move_tests tst/igor/Alignment/Legacy "test_aligner.cpp test_aligner_benchmark.cpp
  test_alignment_cigar.cpp AlignerTestUtils.cpp AlignerTestUtils.h"
move_tests tst/igor/Model/Legacy "test_Contexts.cpp test_EventUtils.cpp test_deletion_iterate.cpp
  test_dinucl_markov_iterate.cpp test_gene_choice_iterate.cpp test_insertion_iterate.cpp
  test_event_capabilities.cpp test_junction_geometry.cpp test_safety_matrix.cpp
  test_span_profile.cpp test_bound_tightness.cpp test_error_rate_bound.cpp
  test_model_format_v2.cpp test_model_json.cpp test_proba_bound_benchmark.cpp
  test_genetic_code.cpp test_utils.cpp test_utils.h"
move_tests tst/igor/Inference/Legacy  "test_inference.cpp"
move_tests tst/igor/Generation/Legacy "test_generation.cpp test_fast_sampling.cpp"
mv_if_present tst/igor/Core/LayerClaim.h tst/igor/LayerClaim.h
mv_if_present tst/igor/Core/CMakeLists.txt tst/igor/Core/CMakeLists.txt.old

for name in test_InferenceEngine test_InferenceHandlers test_ModelInference \
            test_NormalizationCompareLegacy; do
  mv_if_present "tst/igor/Model/$name.cpp" "tst/igor/Inference/$name.cpp"
done
for name in test_SamplingEngine test_SamplingHandlers test_RecombinationModel; do
  mv_if_present "tst/igor/Model/$name.cpp" "tst/igor/Generation/$name.cpp"
done

echo "step1a_moves: done. Now run scripts/migrate_layers/step1a_includes.py"
