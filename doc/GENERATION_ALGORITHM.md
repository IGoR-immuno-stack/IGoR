# How IGoR Generates Immune Receptor Sequences

## Key Source Files

| File | Role |
|---|---|
| `src/igor/Core/GenModel.cpp` | `generate_unique_sequence()` — the main generation loop |
| `src/igor/Core/GenModel.h` | Generation method declarations |
| `src/igor/Core/Rec_Event.h` | Base class with `draw_random_realization_with_engine()` virtual method |
| `src/igor/Core/Genechoice.cpp` | Gene selection sampling (V, D, J) |
| `src/igor/Core/Deletion.cpp` | Deletion sampling (v_3_del, j_5_del, etc.) |
| `src/igor/Core/Insertion.cpp` | Insertion length sampling (vj_ins, vd_ins, etc.) |
| `src/igor/Core/Dinuclmarkov.cpp` | Dinucleotide Markov chain for junction filling |
| `src/igor/Model/CategoricalHandler.tpp` | Categorical sampling from conditional distributions |
| `src/igor/Model/MarkovHandler.tpp` | Markov transition sampling for nucleotide sequences |
| `src/igor/Model/InferenceEngine.h` | Holds all handlers, keyed by event nickname |
| `src/igor/Core/LegacyBridge.tpp` | Imports legacy `Model_Marginals` into the engine |
| `models/human/tcr_alpha/` | TCR alpha model (VJ, used in examples below) |
| `models/human/tcr_beta/` | TCR beta model (VDJ) |
| `tst/test_generation.cpp` | KL divergence validation tests for both generation paths |

## The Biology: V(D)J Recombination

Your immune system creates diverse antibodies/T-cell receptors by randomly assembling gene segments. For a **TCR beta** chain (VDJ model), the steps are:

1. Pick a V gene segment (e.g., TRBV7-2)
2. Pick a D gene segment (e.g., TRBD1)
3. Pick a J gene segment (e.g., TRBJ2-1)
4. Delete some nucleotides from the 3' end of V
5. Delete some nucleotides from the 5' end of D
6. Delete some nucleotides from the 3' end of D
7. Delete some nucleotides from the 5' end of J
8. Insert random nucleotides between V and D (VD junction)
9. Insert random nucleotides between D and J (DJ junction)
10. Fill junction nucleotides using dinucleotide Markov chain

The result is a single recombined DNA sequence. Each step has a **probability distribution** learned from real data.

## The Bayesian Network

These events aren't all independent. IGoR models them as a **Bayesian network** (directed acyclic graph):

```
TCR alpha (VJ):                    TCR beta (VDJ):

  v_choice                           v_choice    j_choice
     │                                  │           │
  j_choice (conditioned on V)        d_gene (conditioned on V, J)
   │    │                            │    │         │
v_3_del j_5_del                   v_3_del d_5_del d_3_del j_5_del
         │                                   │        │
      vj_ins                              vd_ins   dj_ins
         │                                   │        │
     vj_dinucl                           vd_dinucl dj_dinucl
```

So `j_choice` in TCR alpha has **shape `{68, 103}`** — 68 J genes × 103 V genes. The probability of picking J gene #5 depends on *which* V gene was chosen.

## The Algorithm

`GenModel::generate_unique_sequence()` proceeds as follows:

1. Get events in **topological order** (parents before children)
   → e.g., `[v_choice, j_choice, v_3_del, j_5_del, vj_ins, vj_dinucl]`

2. For each event:
   1. Look up which parents were already sampled (`sampled_indices` map)
   2. Use parent indices to slice the correct conditional distribution
   3. Sample a realization from that distribution
   4. Record the chosen index in `sampled_indices`
   5. Modify `constructed_sequences` accordingly

## Concrete Example (TCR alpha, engine path)

### Step 1: `v_choice` (no parents)

```
handler has shape {103}  → 1D tensor, 103 V genes
Sample from [0.015, 0.008, 0.032, ...]  → picks index 42
constructed_sequences[V_GENE_SEQ] = "ATGCAG...TCG"  (the V42 gene sequence)
sampled_indices["v_choice"] = 42
```

### Step 2: `j_choice` (parent: v_choice)

```
handler has shape {68, 103}  → 2D tensor
parent_indices = [42]  (looked up from sampled_indices["v_choice"])
Slice row for v_choice=42 → conditional distribution P(J | V=42)
Sample from that row → picks index 7
constructed_sequences[J_GENE_SEQ] = "GAAT...ACTT"
sampled_indices["j_choice"] = 7
```

### Step 3: `v_3_del` (parent: v_choice)

```
handler has shape {21, 103}  → 21 possible deletion lengths, conditioned on V
parent_indices = [42]
Sample from P(del | V=42) → picks index 3 (meaning delete 3 nucleotides)
constructed_sequences[V_GENE_SEQ] trimmed from 3' end by 3 nucleotides
sampled_indices["v_3_del"] = 3
```

### Step 4: `j_5_del` (parent: j_choice)

Same logic as `v_3_del`, but deletes from the 5' end of the J gene.

### Step 5: `vj_ins` (no parents, typically)

```
handler has shape {41}  → 0 to 40 inserted nucleotides
Sample → picks index 5 (insert 5 nucleotides)
constructed_sequences[VJ_INS_SEQ] = "NNNNN"  (placeholder, 5 chars)
sampled_indices["vj_ins"] = 5
```

### Step 6: `vj_dinucl` (Markov chain, fills the N's)

```
MarkovHandler has shape {4, 4}  → transition matrix P(nt_i | nt_{i-1})
For first nucleotide: sample from marginal distribution → "A"
For nucleotide 2: sample P(· | A) → "C"
For nucleotide 3: sample P(· | C) → "G"
For nucleotide 4: sample P(· | G) → "T"
For nucleotide 5: sample P(· | T) → "A"
constructed_sequences[VJ_INS_SEQ] = "ACGTA"
```

### Final Assembly

```
V gene (trimmed)  +  VJ insertion  +  J gene (trimmed)
"ATGCAG..."       +    "ACGTA"     +   "...ACTT"
```

That's one complete synthetic sequence. `generate_sequences_with_engine()` repeats this N times.

## Event Types and Their Handlers

| Event Type | Handler | Parameters | Example |
|---|---|---|---|
| `Gene_choice` | `CategoricalHandler` | P(gene) or P(gene \| parent) | v_choice, j_choice, d_gene |
| `Deletion` | `CategoricalHandler` | P(n_deleted \| gene) | v_3_del, j_5_del, d_5_del, d_3_del |
| `Insertion` | `CategoricalHandler` | P(n_inserted) | vj_ins, vd_ins, dj_ins |
| `Dinucl_markov` | `MarkovHandler` | P(nt_i \| nt_{i-1}) | vj_dinucl, vd_dinucl, dj_dinucl |

## Code Entry Points

- **Legacy path**: `GenModel::generate_sequences()` → `draw_random_realization()` per event
- **Engine path**: `GenModel::generate_sequences_with_engine()` → `draw_random_realization_with_engine()` per event
- **Model files**: `models/human/tcr_alpha/` and `models/human/tcr_beta/`
- **Test**: `tst/test_generation.cpp` — KL divergence validation for both paths
