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
# Run with default working directory (current directory)
igor-demo

# Run with custom working directory
igor-demo /path/to/working/directory/
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
