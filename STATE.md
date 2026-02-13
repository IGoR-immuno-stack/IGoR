================================================================================
FINAL STATUS - DEMO AT ISO FUNCTIONALITY WITH DEVELOP BRANCH
================================================================================

DATE: February 11, 2025

OVERALL STATUS:
===============

✓ Demo generation WORKS (VJ and VDJ models)
✓ Type consistency achieved across 30 files
✓ Event name generation fixed (backward compatible)
✓ All regression tests PASS
✓ Deterministic generation verified

KEY FIXES APPLIED:
==================

1. Event Priority Comparator (Rec_Event.h)
   - Changed from `<` back to `>` for priority ordering
   - Higher priority numbers are now processed first (matching develop)

2. DinucMarkov Initialization (Dinuclmarkov.cpp)
   - Added nucleotide realizations (A, C, G, T) in constructor
   - Added sequence_type_id initialization based on event_class
   - Fixed draw_random_realization() to use uniform random generation

3. Deletion Initialization (Deletion.cpp)
   - Added sequence_type_id initialization in constructor
   - Maps V_gene/D_gene/J_gene to correct SequenceTypeRegistry types
   - Fixed draw_random_realization() to handle negative deletion lengths
   - Added absolute value handling for deletion lengths (e.g., -3 -> 3)

4. Event Name Suffix (Rec_Event.cpp)
   - Fixed update_event_name() to append `_size` + to_string(this->size())
   - This ensures edges in model files match event names

5. Model Parsing (Model_Parms.cpp)
   - Event lookups now work with suffixed event names
   - Edges connect properly between events

FILES MODIFIED:
===============

Core Library:
- src/igor/Core/Rec_Event.h - Event_comparator fix
- src/igor/Core/Rec_Event.cpp - update_event_name() fix
- src/igor/Core/Dinuclmarkov.cpp - Initialization and generation fixes
- src/igor/Core/Dinuclmarkov.h - draw_random_common declaration
- src/igor/Core/Deletion.cpp - sequence_type_id and generation fixes
- src/igor/Core/Model_Parms.cpp - Event name lookups

Application:
- app/igor/main.cpp - Event name extraction for alignments

TIMING RESULTS:
===============

Simple VJ Model:
- 10,000 sequences: 0.122s
- ~82,000 sequences/second

TRB Complex VDJ Model:
- 10,000 sequences: 0.717s
- ~14,000 sequences/second

The TRB model is slower due to:
- More events (11 vs 5)
- DinucMarkov processing
- D gene handling (more complex recombination)

REGRESSION TESTS:
================

All tests PASS:
✓ test_regression_vdj.sh - VJ model generation and determinism
✓ test_determinism.sh - Identical outputs with same seed

VERIFICATION:
=============

Sample generated sequence (TRB model):
GATACTGGAGTCTCCCAGAACCCCAGACACAAGATCACAAAGAGGGGACAGAATGTAACTTTCAGGTGTGATCCAATTTCTGAACACAACCGCCTTTATTGGTACCGACAGACCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCCAGAATGAAGCTCAACTAGAAAAATCAAGGCTGCTCAGTGATCGGTTCTCTGCAGAGAGGCCTAAGGGATCTCTTTCCACCTTGGAGATCCAGCGCACAGAGCAGGGGGACTCGGCCATGTATCTCTGTGTGTATTTTGGCCCAGGCACCCGGCTGACAGTGCTCG

COMMITS:
========

1. Type consistency fixes (30 files) - Gene_class to int conversion
2. Event name generation fix - _size suffix handling
3. DinucMarkov initialization - realizations and sequence_type_id
4. Deletion fixes - sequence_type_id and negative length handling
5. Priority comparator fix - reverted to > for correct ordering

NEXT STEPS:
===========

1. Test inference with valid TRB data (not synthetic sequences)
2. Verify likelihood values match expected ranges
3. Consider adding debug mode for scenario rejection reasons
4. Performance optimization if needed

================================================================================
PREVIOUS INVESTIGATION NOTES (PRESERVED FOR REFERENCE)
================================================================================

PROBLEM:
The refactoring for the tandem feature changed events_map from Gene_class (static
enum) to int (dynamic runtime type). This was applied INCONSISTENTLY, causing:
- Type mismatches in event lookups
- Event name generation failures
- Broken inference (0 scenarios, NaN likelihoods)

SOLUTION APPLIED:
Completed the exhaustive conversion of Gene_class to int throughout the codebase
(see 30 files listed in git history).

KEY TECHNICAL NOTES:
- All events_map lookups now use tuple<Event_type, int, Seq_side>
- All alignment storage uses unordered_map<int, vector<Alignment_data>>
- Function parameters changed from Gene_class to int
- Member variables changed from Gene_class to int
- Added hash specialization for tuple<Event_type, int, Seq_side>
- Preserved to_string(Gene_class) for backward-compatible event names
- Fixed Rec_Event::update_event_name() to cast int back to Gene_class for naming

The type system is now consistent, and the demo achieves ISO functionality
with the develop branch for sequence generation.
