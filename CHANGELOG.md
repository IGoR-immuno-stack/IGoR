Following are change highlights associated with official releases.  Important
bug fixes are all mentioned, but some internal enhancements are omitted here for
brevity.  Much more detail can be found in the git revision history:

    https://github.com/qmarcou/IGoR

* Unreleased


  Fixed:
  - Fix Pixi manifest parsing so build and regression tasks can run

  New features:
  - Build and Dev env:
    - Add a Pixi profiling task and BENCHMARK.md docs for perf/callgrind-based aligner profiling 
    - Add debugger quality-of-life config (.gdbinit STL pretty-printing, VS Code launch/debug settings)
  - CLI
    - Route the POSIX-style command interface through CLI11 for consistent help, version, and option validation
    - Add full configured pipeline replay and richer reproducibility manifests with runtime and compiler metadata
    - New output.scenario_tie_digits and output.scenario_keep_ties options, controlling
      the precision best scenario probabilities are compared on and whether scenarios
      tied with the last recorded one are all reported. Both are recorded in the run
      manifest. Sequences whose N-th best scenario is degenerate now get more than
      output.scenarios rows instead of an arbitrary member of the degenerate group
  - Sequence sampler:
    - Add fast parallel sequence generator (FastGenerator) with 100x+ speedup
      using precomputed CDFs, binary search/alias method sampling, and OpenMP
      multi-threading
    - New CLI flags --fast and --threads N for the -generate command
    - New sampling primitives: CategoricalSampler, ConditionalSampler,
      DinucleotideMarkovSampler with O(1) alias and O(log n) binary search
    - New GenModel::generate_sequences_fast() public method
  - Aligner: 
    - Add core and AIRR-compliant extended CIGAR string representations for
    alignments, with bidirectionnal conversion helpers
    - Improved Alignment_data struct and API (categorized mismatch/in/dels and offset getters)
    - Add a switch to enable/disable alignment extension (mismatch tracking
    beyond the core alignment bounds)
    - Support more leading/trailing free-end configurations in the
    Smith-Waterman aligner, and encapsulate its DP configuration/state into
    SwDPConfig/SwAlignmentMode/swalign:: helpers instead of long parameter
    lists

  Refactor:
    - Aligner: refactor core sw_align routine into dedicated namespace
    - Replace std::forward_list with std::vector for Alignment_data
    insertions/deletions (kept sorted), simplifying downstream consumers

  Performance:
    - Aligner: 
      - Only initialize the DP boundary trackers instead of the full matrix (O(N+M) instead of O(N*M))
      - Inline and restructure the per-cell Smith-Waterman fill loop (shared
        index arithmetic, packed candidate tracking, branchless argmax move
        selection); measured -23% CPU time / -20% cycles / -89% branch misses on a single-threaded profiling benchmark, and -24.7% wall-clock on a 1000-read pipeline benchmark
      - Improve alignment traceback and limit to in bounds and score acceptable alignments
      - Add 

  Bug fixes:
  - Fix uninitialized iterator bug in Model_marginals::swap_events_order
    (event_1_iterator / event_2_iterator were swapped)
  - Aligner:
    - fix first nucleotide mismatch detection
    - fix orphan duplicate alignment due to SW matrix filling pattern
    - fix best_only logic: reports the best alignment inside the offset bounds, instead of no alignment is the overall best alignment was outside the bounds.
    - Fix J gene alignment preset from local to semi-global
    - Fix offset computation bugs in reversed alignment  
  - Inference:
    - Make the EM reduction independent of the thread count and of the schedule timing,
      by accumulating into fixed per chunk slots merged in chunk order. Repeated runs of
      the same command no longer differ in the last bits of the inferred parameters,
      which used to make the recorded best scenarios flip between runs
    - Order scenarios that are degenerate under the model on their realization indices
      instead of on rounding noise
  - Counters: fix out of bounds reads in the best scenarios and errors counters
    insertion loops, and when dumping a scenario without mismatch

  Tests:
  - Add unit tests for all sampling primitives (16 test cases)
  - Add statistical convergence tests for VJ (TCRα) and VDJ (TCRβ) models
    using KL-divergence and entropy bounds
  - Add inference round-trip test (generate → infer → compare)
  - Add a comprehensive aligner/CIGAR test suite (test_aligner.cpp,
  test_alignment_cigar.cpp, AlignerTestUtils) covering V/D/J gene
  presets, getters, extended mismatches, and CIGAR round-tripping
  (803 assertions in 14 test cases) and micro benchmarks (test_aligner_benchmark.cpp).

  Miscellaneous:
  - Add submodule initialization step in CI workflow 

* 1.4.0 (April 11, 2020)

  A brief description of the release content's

  New features:
  - Add feature extracting CDR3s based on alignments (by @alfaceor)
  - New IGL and IGK models for humans (by @alfaceor)

  Bug fixes:
  - fix out of bound access for Dinucleotide markov model internal probabilities (by @Thopic)

  Miscellaneous :
  - New script igor-compute_pgen to calculate pgen for a single sequence (by @alfaceor)
  - Add error messages for unknown supplied chain or output subargument (by @qmarcou)
  - Add new option -getdatadir to output the path of default models directory (by @alfaceor)

* 1.3.0 (August 4, 2018)

  A brief description of the release content's

  New features:
  - Add ---best_gene_only boolean for alignments
  - Allow for template specific alignment offsets
  - Allow for reversed offsets use via command line
  - Pygor package installable via pip (by @penuts7644)
  - New documentation website available at: https://qmarcou.github.io/IGoR/

  Bug fixes:
  - Fix autogen creation of man and html README
  - Fix semi colon separated files for gene anchors (human BCR, mice TRB)
  - Fix pygor package broken imports and functions (by @penuts7644)
  - Fix IGoR directory creation and installation problems (by @smoe)
  - Fix bug on number of scenarios for the scenario counter CL
  - Fix alignment threshold for --ntCDR3 sequences

  Miscellaneous :
  - Split documentation in several files.
  - Create a proper manpage.
  - Improve error handling for alignment issues (report faulty genomic template)
  - Improve error handling for some input files reading
  - make the autogen.sh and build_release scripts executable (by @NickEngland)
  - Clean and make the Pygor package PEP8 compliant (by @penuts7644)

* 1.2.0 (April 18, 2018)

  A brief description of the release content's

  New features:
  - Enable CDR3 sequences alignment using defined gene anchors via --ntCDR3
  - better quality random seed generation
  - show a progress bar for alignment, inference/evaluation and generation.
  - add independent Nmer hypermutation models in the repository
  - add mouse beta chain model

  Bug fixes:
  - replace the use of default_random_engine by a 64bits mersenne twister

  Miscellaneous :
  - make a proper licensing under GNU-GPLv3
  - create a ChangeLog for better traceability
  - clean some Eclipse config file from the repo
  - remove events initialization message duplication
  - transition to asciidoc documentation
  - add juman BCR CDR3 anchors
