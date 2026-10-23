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
    - C++23 is now required. New dependencies: nlohmann_json, json-schema-validator (built
      from source on Windows) and mdspan
    - The installed CMake package works from a consumer project, with one target and one
      umbrella header per layer; CI checks it (`pixi run test_package`)
    - New ENABLE_BOUND_INSTRUMENTATION option (`pixi run build_instrumented`), building an
      IGoR that reports how far its pruning bound sits above the scenarios it admits
  - CLI
    - Route the POSIX-style command interface through CLI11 for consistent help, version, and option validation
    - Add full configured pipeline replay and richer reproducibility manifests with runtime and compiler metadata
    - New output.event_names option. The best scenarios and generated realizations files
      now head each event's column with its nickname; `legacy` writes the generated names
      IGoR used before, which repeat for the two D blocks of a tandem model
    - New output.scenario_tie_digits and output.scenario_keep_ties options, controlling
      the precision best scenario probabilities are compared on and whether scenarios
      tied with the last recorded one are all reported. Both are recorded in the run
      manifest. Sequences whose N-th best scenario is degenerate now get more than
      output.scenarios rows instead of an arbitrary member of the degenerate group
  - Model files:
    - A model the v1 text format cannot hold, a tandem D for instance, is read and written
      as a JSON document, recognised by its content. Inference writes such a model's
      parameters to .json files. The document holds the structure (segment layout, events
      keyed by nickname, conditioning graph) apart from the template library, and may carry
      comments. Its marginals are JSON too, inline or in a file of their own, with 17
      significant digits, and each value is matched to its realization by label
    - A JSON model names its JSON Schema in "$schema", so that an editor validates and
      completes it. The schemas are installed under share/igor/schemas/model/, and IGoR
      checks every document against them, naming each violation by its JSON pointer
    - Every event needs a nickname, unique in its model, which now names the event in the
      model and in its outputs. A v1 model file with an event line without one is refused
  - Recombination models:
    - Models with two D segments in tandem (V-D1-D2-J) can be loaded, generated from and
      inferred, as JSON documents. Every segment is present in every scenario: optional
      segments are not supported yet. Inference on such a model is costly: one iteration
      of a tandem model derived from human TRB takes about 250 times as long as TRB's, on
      matched reads
    - Before inference starts, IGoR refuses a model whose events cannot run in its order,
      or that the error model or a counter cannot read, naming the event or the reader.
      The hypermutation error models and the errors, coverage, best scenarios and Pgen
      counters read V-D-J and V-J models only
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
    - Split the sources into layers (Core, Alignment, Math, Model, Inference, Generation,
      Streaming), each with its own namespace and CMake target, the legacy code under
      <Layer>/Legacy (by T. Kloczko)
    - Scenario exploration no longer assumes a V-D-J or V-J layout. Segments are addressed
      by run-time ids from the model's segment registry, each event type has one iterate()
      body, overlap checks and pruning bounds follow from the events' declared roles and
      their order, and the scenario state lives in containers sized from the model.
      Outputs unchanged on the regression tracks. See docs/reference/ITERATE_INTERNALS.md
      and docs/reference/REC_EVENT_CONTRACT.md
    - The legacy sequence generator draws realization indices, then builds each segment by
      id and joins them in the model's segment order. Outputs unchanged

  Performance:
    - Aligner: 
      - Only initialize the DP boundary trackers instead of the full matrix (O(N+M) instead of O(N*M))
      - Inline and restructure the per-cell Smith-Waterman fill loop (shared
        index arithmetic, packed candidate tracking, branchless argmax move
        selection); measured -23% CPU time / -20% cycles / -89% branch misses on a single-threaded profiling benchmark, and -24.7% wall-clock on a 1000-read pipeline benchmark
      - Improve alignment traceback and limit to in bounds and score acceptable alignments
      - Add 
    - Inference:
      - The pruning bound's tables are built once per EM iteration and shared by all
        threads, and folded state by state rather than path by path: single-threaded, from
        about 10 s to 1 s on human IGH, and from 0.2 s to 0.055 s on TRB
      - The bound is read from dense arrays rather than tree maps at every scenario node
        (-9 % inference time on the 1000-read pipeline benchmark)
      - Costs: the containers sized from the model add about 9 % to inference time, and
        the corrected pruning bound (see Bug fixes) sums more scenarios, +13 to 20 % wall
        time on inferences from a uniform model

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
    - The pruning bound could fall below what a scenario realizes, pruning scenarios above
      the probability threshold: an insertion or a deletion counted its own probability
      twice, and the bound the J gene choice reads left out D's deletions. Every inference
      output moves: more scenarios are summed and likelihoods rise
    - Reads whose D is placed exhaustively, without a D alignment (no_d_align): D's 5' end
      was placed one nucleotide off, the D probability compounded over placements, and a
      deletion's reach was underestimated. Their outputs move
    - A Markov chain whose anchor gene is deleted entirely now throws instead of reading
      outside the sequence
    - The gene core credited to the error bound had its sign inverted for V and J and was
      one off for D, and could read outside the error-rate table. No shipped model's
      output moves
  - Counters: fix out of bounds reads in the best scenarios and errors counters
    insertion loops, and when dumping a scenario without mismatch
  - Sequence generation: a V 5' or J 3' deletion trims the end it names, as inference
    does. A deletion past its segment's end, a palindrome longer than its template, or a
    V or J deleted away is refused, as inference drops them. No shipped model reaches
    these cases

  Tests:
  - Add unit tests for all sampling primitives (16 test cases)
  - Add statistical convergence tests for VJ (TCRα) and VDJ (TCRβ) models
    using KL-divergence and entropy bounds
  - Add inference round-trip test (generate → infer → compare)
  - Add a comprehensive aligner/CIGAR test suite (test_aligner.cpp,
  test_alignment_cigar.cpp, AlignerTestUtils) covering V/D/J gene
  presets, getters, extended mismatches, and CIGAR round-tripping
  (803 assertions in 14 test cases) and micro benchmarks (test_aligner_benchmark.cpp).
  - Add a unit suite for each event's iterate() and for Single_error_rate. Its harness
    checks at every hand-off that each memory layer an event claimed was written, and that
    the event did what its capability queries declare
  - Add unit tests for the segment registry, the layered containers, the overlap safety
    matrix, the traversal plan and the bound profiles, a no_d_align regression track, and
    JSON model format tests
  - Tandem D: a V-D1-D2-J model derived from human TRB (scripts/tandem_d) and checked
    against it, 200 reads generated from it and frozen with their alignments, a smoke
    inference in the integration suite, a ten-iteration run hidden from CI, and a script
    comparing an iteration's cost with TRB's

  Miscellaneous:
  - Add submodule initialization step in CI workflow 
  - Add the design and migration plans under docs/, ordered in docs/IMPLEMENTATION_GRAPH.md

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
