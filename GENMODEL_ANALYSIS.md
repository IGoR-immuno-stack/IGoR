# GenModel Class - Detailed Analysis

## Overview

The **GenModel** class is the highest-level abstraction in the IGoR (Inference and Generation of Repertoires) framework. It represents a complete V(D)J recombination generative model capable of:
1. **Inferring model parameters** from observed immune receptor sequences
2. **Generating new synthetic sequences** from the learned probability distributions
3. **Annotating sequences** with recombination event realizations

---

## Class Hierarchy

### Direct Inheritance
- **GenModel** does NOT inherit from any class (standalone design)

### Component Dependencies (Composition)

GenModel is composed of several key components:

```
GenModel
├── Model_Parms (model_parms)
│   └── Bayesian Network structure of V(D)J recombination
├── Model_Marginals (model_marginals)
│   └── Probability distributions for event realizations
├── std::map<size_t, std::shared_ptr<Counter>> (counters_list)
│   └── Statistical counters for tracking sequence features
├── igor::fast::FastGenerator (fast_generator_)
│   └── Optimized sequence generation engine
└── Related Classes:
    ├── Rec_Event (abstract base for recombination events)
    ├── Rec_Event subclasses:
    │   ├── Gene_choice (V, D, J gene selection)
    │   ├── Deletion (P/N deletion events)
    │   ├── Insertion (P/N insertion events)
    │   ├── Dinucl_markov (dinucleotide markov model)
    │   └── Error_rate model
    └── Counter subclasses (various specialized counters)
```

### Related Key Classes

#### 1. **Model_Parms**
- **Purpose**: Encodes the Bayesian Network structure (directed acyclic graph)
- **Contains**:
  - Adjacency lists for event dependencies
  - List of all Rec_Event objects in the model
  - Error rate model (Error_rate object)
  - Methods for topological traversal

#### 2. **Model_Marginals**
- **Purpose**: Stores marginal probabilities for each event realization
- **Structure**: Giant array indexed by event indices and realization offsets
- **Contains**:
  - Probability distributions during inference (posterior frequencies)
  - Methods for normalization and arithmetic operations

#### 3. **Rec_Event** (Abstract Base Class)
- **Purpose**: Node in the Bayesian Network
- **Key Methods**:
  - `iterate()`: Explores all compatible scenarios for a given sequence
  - `draw_random_realization()`: Samples a random realization from marginals
  - `initialize_event()`: Sets up event-specific data structures

#### 4. **Event_realization**
- **Purpose**: Stores a single realization of a recombination event
- **Structure**:
  ```cpp
  struct Event_realization {
      std::string name;           // Gene/event name identifier
      int value_int;              // Integer value (e.g., deletion count)
      std::string value_str;      // String value (e.g., V gene sequence)
      Int_Str value_str_int;      // Hybrid representation
      int index;                  // Unique index in marginals array
  }
  ```

---

## Constructors

### Constructor Overloads

1. **GenModel(const Model_Parms &parms)**
   - Minimal initialization
   - Creates empty marginals and counter list

2. **GenModel(const Model_Parms &parms, const Model_marginals &marginals)**
   - Initializes with custom probability distributions

3. **GenModel(const Model_Parms &parms, const Model_marginals &marginals, const std::map<size_t, std::shared_ptr<Counter>> &counters)**
   - Full initialization with all components

---

## Core Algorithms

### 1. MODEL INFERENCE (infer_model)

**Purpose**: Maximum Likelihood Estimation of model parameters

**Function Signatures**:
```cpp
bool infer_model(
    const std::vector<std::tuple<int, std::string,
        std::unordered_map<int, std::vector<Alignment_data>>>> &sequences,
    const int iterations,
    const std::string path,
    bool fast_iter,
    double likelihood_threshold = 1e-25,
    bool viterbi_like = false,
    double proba_threshold_factor = 0.001,
    double mean_number_seq_err_thresh = INFINITY);
```

#### Input Data Structure
For each sequence:
- `int`: Sequence index
- `std::string`: Nucleotide sequence
- `unordered_map<int, vector<Alignment_data>>`:
  - Key: Gene class (V_gene, J_gene, D_gene)
  - Value: List of alignments with scores

#### Algorithm Overview

```
INFERENCE ALGORITHM:
└── Initialize:
    ├── Load model queue (topological order of events)
    ├── Create index and offset maps for marginal array navigation
    ├── Initialize probability bounds for pruning
    └── Create counter instances

└── Loop iterations = 1 to max_iterations:
    ├── Create new marginals (M^new = 0)
    ├── Create error rate copy

    ├── For each sequence in parallel (#pragma omp):
    │   ├── If fast_iter && iteration==0: Keep only best V and J alignments
    │   ├── Initialize single-sequence marginals to empty
    │   ├── Set max_proba_scenario = likelihood_threshold / proba_threshold_factor
    │   │
    │   ├── Call first_event->iterate():
    │   │   ├── For each realization of current event:
    │   │   │   ├── Compute compatibility with sequence
    │   │   │   ├── Calculate posterior (P(E|seq) ∝ P(seq|E) * P(E))
    │   │   │   ├── Prune if P(scenario) < max_proba_scenario
    │   │   │   ├── Update constructed sequences
    │   │   │   ├── Call next_event->iterate() recursively
    │   │   │   └── Update offset tracking (3'/5' positions)
    │   │   └── Process error realizations
    │   │
    │   ├── Normalize single-sequence marginals
    │   ├── Copy fixed events from previous iteration
    │   ├── Update error rate with sequence likelihood
    │   └── Add weighted marginals to M^new
    │
    ├── Merge thread-local results:
    │   ├── M_global += M^new
    │   ├── ER_global += ER^new
    │   └── Counter_global += Counter^new
    │
    ├── Normalize global marginals by conditioning probabilities
    ├── Update error rate parameters
    ├── Output iteration results to files
    └── Increment iteration counter

└── Return: true (success)
```

#### Key Algorithmic Features

1. **Scenario Pruning**:
   - Only explores scenarios where P(scenario) ≥ likelihood_threshold / proba_threshold_factor
   - Dramatically reduces computational complexity

2. **Viterbi Mode** (viterbi_like=true):
   - Only tracks the single best scenario per sequence
   - Sets proba_threshold_factor = 1.0 (only best path)
   - Fast but less accurate

3. **Fast First Iteration** (fast_iter=true):
   - First iteration only uses best V and J alignments
   - Reduces complexity in first pass
   - Full alignment exploration in subsequent iterations

4. **Thread Parallelization**:
   - Each thread processes independent sequences
   - Thread-local copies of marginals and error rates
   - Critical sections for merging results

5. **Posterior Computation**:
   ```
   P(scenario | sequence) ∝ P(sequence | scenario) × P(scenario)
   ```
   - Incorporates alignment scores
   - Incorporates error model likelihood
   - Incorporates model priors

#### Error Model Integration
- Error_rate object tracks sequence errors (mismatches)
- Weights contribution of each sequence by likelihood
- Updates error probabilities each iteration
- Filters sequences with too many errors (mean_number_seq_err_thresh)

#### Output Files
- `likelihoods.out`: Mean log-likelihood per iteration
- `inference_logs.txt`: Per-sequence details
- `inference_info.out`: Inference parameters and runtime info
- `initial_marginals.txt`, `iteration_N.txt`, `final_marginals.txt`: Model state
- `initial_parms.txt`, `iteration_N_parms.txt`, `final_parms.txt`: Model parameters

---

### 2. SEQUENCE GENERATION

#### 2.1 Standard Generation

**Function**:
```cpp
std::forward_list<std::pair<std::string, std::queue<std::queue<int>>>>
generate_sequences(int n_seq, bool output_realizations);
```

**Algorithm**:
```
├── Create random number generator with current timestamp seed
├── Get model queue (topological order)
├── Get index and offset maps from marginals
│
└── For i = 1 to n_seq:
    └── Call generate_unique_sequence():
        ├── For each event in model_queue:
        │   ├── Sample realization from marginal distribution
        │   ├── Update constructed_sequences map
        │   └── Store realization if requested
        │
        └── Build final sequence by traversing:
            ├── Start from V_GENE_SEQ
            ├── Follow downstream neighbors in SequenceTypeRegistry
            ├── Include junction sequences
            └── Concatenate all pieces -> final_seq
```

**Return Structure**:
```cpp
std::pair<
    std::string,                           // Final nucleotide sequence
    std::queue<std::queue<int>>            // Event realizations
>
```

#### 2.2 File-Based Generation

**Function**:
```cpp
void generate_sequences(
    int n_seq,
    bool output_realizations,
    std::string seq_file_path,
    std::string real_file_path,
    std::list<std::pair<gen_seq_trans, std::shared_ptr<void>>> transformations,
    bool output_only_func = false,
    int seed = -1);
```

**Features**:
- Writes sequences to CSV file with format: `seq_index;nt_sequence`
- Optional realizations file with all event choices
- Supports transformation functions (e.g., CDR3 extraction)
- Seeded randomization

#### 2.3 Fast Sequence Generation (100x+ speedup)

**Function**:
```cpp
void generate_sequences_fast(
    size_t num_sequences,
    const std::string &seq_filename,
    const std::string &real_filename,
    size_t num_threads = 0,
    int64_t seed = -1,
    bool show_progress = true);
```

**Optimization Strategies**:

1. **Precomputed CDFs**
   - Pre-calculate cumulative density functions for all events
   - Eliminates repeated probability calculations

2. **Binary Search / Alias Sampling**
   - Fast O(log n) realization sampling
   - Better cache locality than naive sampling

3. **Multi-threading**
   - Parallel sequence generation
   - Each thread independent (no synchronization needed)
   - Auto-detects optimal thread count

4. **Batched I/O**
   - Buffers multiple sequences before writing
   - Reduces system calls
   - Improves disk throughput

**Related Class**: `igor::fast::FastGenerator`
```cpp
class FastGenerator {
    void initialize(Model_Parms &, Model_marginals &);
    void generate_to_files(
        size_t num_sequences,
        const std::string &seq_filename,
        const std::string &real_filename,
        const FastGeneratorConfig &,
        ProgressCallback callback);
    FastGeneratorStats get_stats() const;
};
```

---

## Key Data Structures

### 1. Sequence Representation

**ntDNA (Integer Representation)**:
- 0 = A, 1 = C, 2 = G, 3 = T (or gap)
- Allows vectorized operations
- Compact bit representation

**Constructed Sequences Map**:
```cpp
std::unordered_map<int, std::string> constructed_sequences
// Key: SequenceTypeRegistry::TypeId (V, D, J, insertions, etc.)
// Value: Sequence fragment
```

### 2. Probability Bounds Map

**Purpose**: Dynamic pruning of low-probability scenarios

**Structure**:
```cpp
Downstream_scenario_proba_bound_map downstream_proba_map
// Multi-level map tracking upper probability bounds
// Allows early termination without exploring all scenarios
```

### 3. Offset Tracking

**Purpose**: Track insertion/deletion events and sequence positions

**Structure**:
```cpp
Seq_offsets_map seq_offsets
// Maps event to 5' and 3' offsets
// Prevents position overlap errors
```

### 4. Index Mapping

**Purpose**: Navigate marginals array efficiently

**Structures**:
- `Index_map`: Current position in marginals array
- `Inverse_offset_map`: Maps event to all its realizations
- `Offsets_map`: Realization to marginals array position

---

## Probability Computation

### Marginal Probability Formula

For a complete scenario (all events realized):
```
P(scenario | sequence) ∝ P(sequence | scenario) × P_marginal(scenario)

Where:
- P(sequence | scenario) = Product of alignment scores × error model likelihood
- P_marginal(scenario) = Product of marginal probabilities for each event
```

### Normalization

After processing all sequences in an iteration:
```
P(event_realization | all_sequences) =
    Sum of posterior contributions from all sequences /
    Normalization constant
```

The normalization is performed using the inverse_offset_map.

---

## Advanced Features

### 1. Error Rate Modeling

The GenModel integrates with error rate models:
- **SingleErrorrate**: Uniform error probability per position
- **HypermutationGlobalErrorrate**: Global error parameters
- **HypermutationfullNmerErrorrate**: N-mer specific errors
- **Dinuclmarkov**: Dinucleotide-based error model

**Integration Process**:
1. For each sequence scenario, compute mismatch positions
2. Update error model weights by sequence likelihood
3. Recompute error probabilities each iteration
4. Use updated probabilities in next iteration's inference

### 2. Counter System

Multiple specialized counters track:
- V/D/J gene usage frequencies
- Deletion/insertion length distributions
- Junction composition statistics
- Custom metrics

### 3. Conditional Probability Handling

The model respects Bayesian Network dependencies:
- Events marked as "fixed" retain previous iteration values
- Conditional events update only when parent events change
- Probability bounds updated respecting conditioning constraints

---

## Performance Characteristics

### Inference Scalability

**Time Complexity per Sequence**:
- Worst case: O(n_events × ∏(realization_counts))
- Pruned case: O(n_events × explored_scenarios)
- With good bounds: O(n_events × polylog(realization_space))

**Memory**:
- Marginal array: O(total_realizations × data_type_size)
- Per-thread overhead: O(n_sequences) marginals copies

### Generation Performance

**Standard Generation**:
- ~1,000-10,000 seq/sec per thread (model dependent)

**Fast Generation**:
- ~100,000-1,000,000 seq/sec per thread
- 100x+ speedup vs standard

### Optimization Parameters

Users can tune:
1. `likelihood_threshold`: Lower = more comprehensive, slower
2. `proba_threshold_factor`: Higher = more pruning, faster
3. `fast_iter`: Skip in iterations >0 for full accuracy
4. `mean_number_seq_err_thresh`: Filter low-quality sequences

---

## Workflow Example

### Typical Usage Pattern

```cpp
// 1. Create model parameters from configuration
Model_Parms params = load_model_from_file("model.txt");

// 2. Create GenModel with initial marginals
GenModel model(params);
model.uniform_initialize();  // Start with uniform distribution

// 3. Load sequences with alignments
auto sequences = load_and_align_sequences("data.fasta");

// 4. Inference
model.infer_model(
    sequences,
    iterations = 10,
    output_path = "./results/",
    fast_iter = true,
    likelihood_threshold = 1e-25,
    viterbi_like = false,
    proba_threshold_factor = 0.001
);

// 5. Generate new sequences
model.generate_sequences_fast(
    num_sequences = 100000,
    seq_filename = "generated_seqs.csv",
    real_filename = "generated_realizations.csv",
    num_threads = 0  // Auto-detect
);
```

---

## File Input/Output

### Parameters Format (parms.txt)
- Event definitions and parameters
- Error rate model specification
- Model topology (dependencies)

### Marginals Format (marginals.txt)
- Realization names and probabilities
- Normalized from previous iteration

### Sequence Format (generated_seqs.csv)
```
seq_index,nt_sequence
0,ATCGATCG...
1,GCTAGCTA...
```

### Realizations Format (generated_realizations.csv)
```
seq_index;V_choice;D_choice;J_choice;V_del;D_del_5;D_del_3;J_del;P_insert;N_insert;Errors
0;V1;D1;J1;(2);(1,0);(1);(3);(AT);(GC);...
```

---

## Key Methods Summary

| Method | Purpose | Complexity |
|--------|---------|-----------|
| `infer_model()` | Estimate parameters via MLE | O(iter × seqs × scenarios) |
| `generate_sequences()` | Generate and save sequences | O(n_seqs × n_events) |
| `generate_sequences_fast()` | Optimized generation | O(n_seqs × n_events / threads) |
| `get_fast_generator()` | Access optimization engine | O(1) with lazy init |
| `write2txt()` | Save model to disk | O(n_realizations) |
| `readtxt()` | Load model from disk | O(n_realizations) |

---

## Design Patterns

### 1. **Strategy Pattern** (Event Realization)
- Different Rec_Event subclasses implement event-specific logic
- GenModel treats all events polymorphically via shared_ptr<Rec_Event>

### 2. **Builder Pattern** (Model Construction)
- Model_Parms builds Bayesian Network structure
- Model_Marginals populates probability values
- GenModel orchestrates complete model

### 3. **Template Method Pattern** (Iteration)
- Rec_Event::iterate() provides the framework
- Subclasses implement realization-specific scoring

### 4. **Observer Pattern** (Counters)
- Counter objects register with GenModel
- Notified of sequence processing results
- Collect statistics asynchronously

### 5. **Lazy Initialization** (FastGenerator)
- FastGenerator created on first use via get_fast_generator()
- Avoids overhead if fast generation not needed

---

## Summary

The **GenModel** class is a sophisticated probabilistic framework for:
- **Inferring** complex multi-event generating processes from observed data
- **Generating** synthetic sequences respecting learned distributions
- **Integrating** error models and conditional dependencies

Its design leverages Bayesian Networks for interpretability while maintaining computational efficiency through:
- Adaptive scenario pruning
- Parallel processing
- Vectorized operations
- Multi-level caching

The class is the core engine of the IGoR software for understanding V(D)J immune receptor generation processes.
