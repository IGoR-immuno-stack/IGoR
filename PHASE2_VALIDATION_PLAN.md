# Phase 2: Parallel Execution Validation Plan

**Objective**: Verify numerical equivalence between the legacy `Model_marginals` and the new `InferenceEngine` by running them simultaneously on the same data.

**Status**: Planned
**Target Verification Model**: Mouse TCRβ (Demo Dataset)

**Scope** (per reviewer feedback):
- **Single model**: Mouse TCRβ from demo (`demo/murugan_naive1_noncoding_demo_seqs.txt`)
- **2 iterations**: Sufficient to prove equivalence (no need for full convergence)
- **Real sequences**: Use actual demo workflow from `app/igor/main.cpp`
- **No multi-model testing**: TCRβ validation is adequate

---

## 1. Strategy: Parallel Execution (No Instrumentation)

**Key insight from reviewer**: All inference logic lives in `GenModel::infer_model()`. Instead of adding hooks to production code, we **clone the inference function** for test-only parallel execution.

We will create a **test-only function** `infer_model_parallel_validation()` that duplicates the core E-M loop from `GenModel::infer_model()` but processes each sequence through **both** systems:

1. **Legacy path**: Existing `Counter::count_scenario()` → accumulates to `Model_marginals`
2. **New path**: Duplicate recursion logic → accumulates to `InferenceEngine`

**Key insight from reviewer**: All inference logic lives in `GenModel::infer_model()`. We clone this function for validation rather than using hooks.

### Core Structure

**File**: `tst/igor/Integration/test_ParallelEM.cpp`

```cpp
// Clone of GenModel::infer_model() E-step loop
void infer_model_parallel_validation(
    const vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>>& sequences,
    Model_marginals& legacy_marginals,
    InferenceEngine& engine,
    const Model_Parms& model_parms,
    int n_iterations
) {
    // Prepare model queue and indices (same as GenModel::infer_model)
    queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    unordered_map<Rec_Event_name, int> index_map =
        legacy_marginals.get_index_map(model_parms, model_queue);

    for (int iter = 0; iter < n_iterations; ++iter) {
        // Reset accumulators
        Model_marginals new_legacy_marginals(model_parms);
        engine.reset_all_accumulators();

        long double legacy_ll = 0.0;
        long double engine_ll = 0.0;

        // E-step: Process each sequence through BOTH systems
        for (const auto& seq_tuple : sequences) {
            const string& sequence = std::get<1>(seq_tuple);
            const auto& alignments = std::get<2>(seq_tuple);

            // LEGACY PATH - real Counter::count_scenario
            Counter counter(&model_parms, legacy_marginals, ...);
            long double seq_ll_legacy = counter.count_scenario(sequence, alignments, ...);
            legacy_ll += seq_ll_legacy;
            // ^ Internally accumulates to new_legacy_marginals

            // NEW PATH - replicate recursion, accumulate to engine
            long double seq_ll_engine = count_scenario_to_engine(
                sequence, alignments, engine, model_parms, model_queue, index_map
            );
            engine_ll += seq_ll_engine;

            // Per-sequence validation
            REQUIRE(seq_ll_engine == Approx(seq_ll_legacy).epsilon(1e-12));
        }

        // M-step: Normalize
        new_legacy_marginals.normalize();
        engine.normalize_all();

        // Compare parameters
        REQUIRE(compare_marginals_to_engine(new_legacy_marginals, engine) < 1e-12);

        // Update for next iteration
        legacy_marginals = new_legacy_marginals;
    }
}
```

---

## 2. The Recursion Clone: `count_scenario_to_engine()`

**File**: `tst/igor/Integration/ParallelRecursion.h`

This helper function **duplicates the scenario exploration logic** from the legacy `Counter::count_scenario()` but routes accumulation to the `InferenceEngine` instead of `Model_marginals`.

**Key components** (simplified structure):

```cpp
long double count_scenario_to_engine(
    const string& sequence,
    const unordered_map<Gene_class, vector<Alignment_data>>& alignments,
    InferenceEngine& engine,
    const Model_Parms& parms,
    queue<shared_ptr<Rec_Event>> model_queue,
    const unordered_map<Rec_Event_name, int>& index_map
) {
    long double sequence_likelihood = 0.0;

    // Replicate the scenario exploration from Counter::count_scenario
    // This is a SIMPLIFIED outline - real implementation follows iterate() logic

    // 1. Initialize scenario exploration
    unordered_map<tuple<Event_type, Gene_class, Seq_side>, shared_ptr<Rec_Event>> events_map;
    // ... setup events from model_queue ...

    // 2. Iterate through model events (same order as legacy)
    for (auto& event : model_queue) {
        // Explore realizations for this event
        // ... iterate() logic ...
    }

    // 3. At scenario completion (same point as iterate_wrap_up)
    auto accumulate_scenario = [&](long double scenario_prob) {
        // Route each event realization to corresponding handler
        for (const auto& [key, event] : events_map) {
            string handler_name = event->get_nickname();
            auto handler = engine.find_handler(handler_name);
            if (!handler) continue; // Deletion events

            // Get realization (single index for categorical, multiple for Markov)
            int realization_idx = event->get_current_realization_index();

            // Special handling for Markov (multi-dimensional)
            if (event->get_type() == Dinuclmarkov_t) {
                // Access flat indices for all transitions
                auto* markov_event = dynamic_cast<Dinucl_markov*>(event.get());
                const auto& transition_indices = markov_event->get_realization_queue();
                for (int flat_idx : transition_indices) {
                    handler->accumulator().data()[flat_idx] += scenario_prob;
                }
            } else {
                handler->accumulate(realization_idx, scenario_prob);
            }
        }
        sequence_likelihood += scenario_prob;
    };

    // 4. Trigger accumulation for all explored scenarios
    // ... call accumulate_scenario(prob) for each valid scenario ...

    return sequence_likelihood;
}
```

**Note**: This is a **conceptual outline**. The real implementation will follow the structure of `Rec_Event::iterate()` and `iterate_wrap_up()` closely, ensuring identical scenario exploration order and probability calculations.

## 3. Test Fixture: Real Demo Workflow

**File**: `tst/igor/Integration/test_ParallelEM.cpp`

This integration test uses the **actual demo workflow** from `app/igor/main.cpp` to validate on real data.

**Test Structure**:

```cpp
TEST_CASE("Phase 2: Parallel EM Validation on Demo Dataset", "[phase2][integration]") {
    // 1. Setup - exactly as in demo code (lines 1512-1667 of main.cpp)

    // Load genomic templates
    vector<pair<string, string>> v_genomic =
        read_genomic_fasta(IGOR_DATA_DIR + "/demo/genomicVs_with_primers.fasta");
    vector<pair<string, string>> d_genomic =
        read_genomic_fasta(IGOR_DATA_DIR + "/demo/genomicDs.fasta");
    vector<pair<string, string>> j_genomic =
        read_genomic_fasta(IGOR_DATA_DIR + "/demo/genomicJs_all_curated.fasta");

    // Create TCRβ model structure
    Gene_choice v_choice(V_gene, v_genomic);
    v_choice.set_nickname("v_choice");
    // ... rest of model construction (deletions, insertions, markov) ...

    Model_Parms parms;
    parms.add_event(&v_choice);
    // ... add all events ...

    Model_marginals legacy_marginals(parms);
    legacy_marginals.uniform_initialize(parms);

    Single_error_rate error_rate(0.001);
    parms.set_error_ratep(&error_rate);

    // 2. Load sequences
    vector<pair<const int, const string>> indexed_seqlist =
        read_txt(IGOR_DATA_DIR + "/demo/murugan_naive1_noncoding_demo_seqs.txt");

    // 3. Align sequences (as in demo)
    Aligner v_aligner = Aligner(nuc44_sub_matrix, 50, V_gene);
    v_aligner.set_genomic_sequences(v_genomic);
    // ... align V, D, J ...

    auto sorted_alignments = read_alignments_seq_csv_score_range(...);
    vector<tuple<int, string, unordered_map<Gene_class, vector<Alignment_data>>>>
        sorted_alignments_vec = map2vect(sorted_alignments);

    // 4. Initialize InferenceEngine
    InferenceEngine engine(parms);
    LegacyBridge::import_from_legacy(engine, legacy_marginals, parms);

    // 5. Run parallel validation (2 iterations as reviewer suggested)
    SECTION("Two EM iterations produce identical results") {
        infer_model_parallel_validation(
            sorted_alignments_vec,
            legacy_marginals,
            engine,
            parms,
            2  // Only 2 iterations needed to prove equivalence
        );
        // If we reach here without assertion failures, equivalence is proven
    }
}
```

**Success Criteria** (per reviewer feedback):
- ✅ Both systems produce identical likelihoods per sequence (1e-12 tolerance)
- ✅ After 2 iterations, all parameters match (1e-12 tolerance)
- ✅ No need to test on multiple models - TCRβ demo is sufficient

---

## 4. Implementation Steps

**Per reviewer feedback**: Skip hooks/instrumentation entirely. Instead, clone the inference logic for test purposes.

1.  **Create test directory**: `tst/igor/Integration/`
2.  **Implement recursion clone**: `ParallelRecursion.h` with `count_scenario_to_engine()`
    - Study `Rec_Event::iterate()` and `iterate_wrap_up()` to understand scenario exploration
    - Replicate the same event queue processing order
    - Route accumulation to `InferenceEngine` handlers instead of `Model_marginals`
3.  **Create parallel validation function**: `infer_model_parallel_validation()`
    - Clone E-M loop structure from `GenModel::infer_model()`
    - Process each sequence through both legacy and new paths
    - Assert per-sequence and per-iteration equivalence
4.  **Write test case**: `test_ParallelEM.cpp`
    - Load demo model (TCRβ) exactly as `app/igor/main.cpp` does
    - Load demo sequences: `demo/murugan_naive1_noncoding_demo_seqs.txt`
    - Run 2 iterations (sufficient per reviewer)
5.  **Debug and iterate**: Fix any discrepancies in scenario exploration or index mapping

## 5. Status

I have completed **Steps 1-4** of the **Implementation Steps** (Section 4 of `PHASE2_VALIDATION_PLAN.md`).

**Completed:**
*   ✅ **Step 1: Create test directory**: Created `tst/igor/Integration/`.
*   ✅ **Step 2: Implement recursion clone**: Created `ParallelRecursion.h` with `count_scenario_to_engine()` and helper function `accumulate_from_legacy_flat()`. Validated with mock-based unit test.
*   ✅ **Step 3: Create parallel validation function**: Implemented `infer_model_parallel_validation()` that clones the E-M loop from `GenModel::infer_model()`, processes sequences through both legacy and new paths, and asserts per-sequence and per-iteration equivalence.
*   ✅ **Step 4: Write test case**: Created `test_ParallelEM.cpp` that loads the real Mouse TCRβ demo model with:
    - Full 15×15 nuc44 substitution matrix (handles ambiguous nucleotides)
    - Proper deletion/insertion ranges matching demo (not hardcoded 0)
    - Correct event priorities and topology
    - Demo alignment thresholds (V:55, D:35, J:10)
    - Real sequences from `demo/murugan_naive1_noncoding_demo_seqs.txt`
    - Validates 2 EM iterations with 5 sequences (sufficient to prove equivalence)

**Next (To Do):**
*   **Step 5: Debug and iterate** (**IN PROGRESS**):
    - ✅ Fixed memory management issue with `Next_event_ptr` array deleter
    - ✅ Added crude upper bound updates between iterations
    - ✅ Reduced to 2 sequences for faster iteration
    - 🔄 **Current Issue**: `unordered_map::at: key not found` on second sequence - investigating index_map/event name mismatch
    - 🔄 Need to ensure all event names in Model_Parms match keys in index_map
    - 🔄 Once passing, increase to 5 sequences and run full 2-iteration validation