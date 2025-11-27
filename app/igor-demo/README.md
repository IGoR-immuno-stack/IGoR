# igor-demo

Standalone demonstration application for IGoR (Inference and Generation of Repertoires).

## Overview

This executable extracts the demo code from the main IGoR application into a standalone program. It demonstrates the complete workflow of IGoR on a small TCRβ dataset.

## What it does

The demo performs the following steps:

1. **Read TCRβ genomic templates**
   - V genes from `genomicVs_with_primers.fasta`
   - D genes from `genomicDs.fasta`
   - J genes from `genomicJs_all_curated.fasta`

2. **Align sequences to templates**
   - Reads sequences from `murugan_naive1_noncoding_demo_seqs.txt`
   - Aligns against V, D, and J genes using the nuc44 substitution matrix
   - Writes alignment results to CSV files

3. **Construct a TCRβ model**
   - Gene choice events (V, D, J)
   - Deletion events (V 3', D 5'/3', J 5')
   - Insertion events (VD, DJ)
   - Dinucleotide Markov models for insertions
   - Simple error rate model

4. **Write and read model parameters**
   - Demonstrates serialization of model parameters and marginals

5. **Infer model from sequences**
   - Performs 4 iterations of Expectation-Maximization
   - Collects statistics: coverage, best scenarios, Pgen, errors

6. **Generate sequences**
   - Generates 100 sequences from the inferred model
   - Writes sequences and realizations to CSV files

## Usage

```bash
# Run with default output location (current directory)
# Creates: ./igor_demo/
igor-demo

# Run with custom output parent directory
# Creates: /path/to/output/igor_demo/
igor-demo /path/to/output

# Example: output to /tmp
# Creates: /tmp/igor_demo/
igor-demo /tmp
```

**Note:** The program automatically creates an `igor_demo/` subdirectory within the path you provide. If you want output in `/tmp/my_results/`, run `igor-demo /tmp` and you'll get `/tmp/igor_demo/` (not `/tmp/my_results/`).

## OpenMP Optimization

The demo uses OpenMP for parallel processing during model inference. You can control parallelization with environment variables:

```bash
# Set number of threads (default: number of CPU cores)
export OMP_NUM_THREADS=4
igor-demo

# Or set it inline
OMP_NUM_THREADS=8 igor-demo /tmp

# Disable nested parallelism (recommended for most cases)
export OMP_NESTED=false

# Set thread affinity for better cache performance
export OMP_PROC_BIND=true

# Control scheduling (dynamic is default in IGoR)
export OMP_SCHEDULE="dynamic"

# Combined example for optimal performance:
OMP_NUM_THREADS=8 OMP_PROC_BIND=true OMP_NESTED=false igor-demo
```

### Performance Tips

- **OMP_NUM_THREADS**: Set to the number of physical cores (not hyperthreads) for best performance
- **OMP_PROC_BIND=true**: Pins threads to cores, improving cache locality
- **OMP_NESTED=false**: Disables nested parallelism (not used in IGoR anyway)
- For small datasets (like this demo), using 4-8 threads is usually optimal
- For larger datasets, you can use all available cores

**Example on a 16-core machine:**
```bash
OMP_NUM_THREADS=16 OMP_PROC_BIND=true igor-demo /tmp
```

## Output

The demo creates the following directory structure:

```
<working_dir>/igor_demo/
├── murugan_naive1_noncoding_demo_seqs_alignments_V.csv
├── murugan_naive1_noncoding_demo_seqs_alignments_D.csv
├── murugan_naive1_noncoding_demo_seqs_alignments_J.csv
├── murugan_naive1_noncoding_demo_seqs_indexed_seq.csv
├── demo_write_model_parms.txt
├── demo_write_model_marginals.txt
├── run_demo/
│   ├── coverage/
│   ├── best_scenarios/
│   ├── pgen/
│   └── errors/
├── generated_seqs_indexed_demo.csv
└── generated_seqs_realizations_demo.csv
```

## Build

The executable is built automatically when you build IGoR:

```bash
mkdir build && cd build
cmake ..
make igor-demo
```

## Source Code

This is a direct extraction of the `if(run_demo)` section from the main IGoR application (`app/igor/main.cpp`). The code is unchanged to maintain exact behavior for reference and testing purposes.

## Note

This demo requires the IGoR data files to be available. By default, it looks for data files in `../../demo/` relative to the executable location. Make sure the demo data files are present:

- `../../demo/genomicVs_with_primers.fasta`
- `../../demo/genomicDs.fasta`
- `../../demo/genomicJs_all_curated.fasta`
- `../../demo/murugan_naive1_noncoding_demo_seqs.txt`

These files are typically installed with IGoR.
