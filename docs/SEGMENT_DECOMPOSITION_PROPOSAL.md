# Segment decomposition: an architecture proposal for IGoR

Status: superseded on 2026-10-02 by `ARCHITECTURE_SYNTHESIS.md` (Q. Marcou), which merges this
proposal with `SEGMENT_DECOMPOSITION_REVIEW.md` and the Model module design. Kept as the record of
the 2026-10-01 discussion. Where the two disagree, the synthesis holds: events keep their stated
attributes (targets, side, library), the per-segment rules are admissibility conditions of one
engine, the layout is a DAG, handlers are declared per node, and the per-segment dynamic
programming of §5 is one inference engine among several next to the legacy DFS.
Date: 2026-10-01.
Source: discussion between T. Kloczko and Q. Marcou on 2026-10-01. Drafted from that discussion.
Related: `REC_EVENT_CAPABILITY_REFACTORING_PLAN.md` (Phases A to E), `ITERATE_GENERIC_REWRITE_PLAN.md`,
`TK_REFACTORING_MERGE_ANALYSIS.md` (§0.1 the two graphs, §6.1 the sample/apply split, §7.2 Phase D).

## 1. Purpose

Restate IGoR's model and its two computations, inference and generation, on four concepts:
events, segments, a layout and a conditioning graph. The event stays the unit of
parameterization. Segments group events by the piece of sequence they act on. The layout
orders the segments. The conditioning graph carries the statistical dependencies between
events. Inference computes per segment and assembles along the layout. Generation is the same
structure run forward.

The probabilistic model does not change. What changes is how a scenario's probability is
computed: by a sum over the layout with per-segment tables, instead of the interleaved
exploration in `iterate()`. The data structures that exist for that exploration are the main
thing this proposal removes.

Alignment, the CLI and the streaming layer are out of scope.

## 2. Concepts

### 2.1 Event

An event has one role and one parameterization.

| Role | Meaning | Current events |
|---|---|---|
| Create | produces sequence content | `GeneChoice` (a germline template); `DinuclMarkov` (the nucleotides of a junction) |
| Modify | changes the extent of content that already exists | `Deletion` (trims or, with a negative value, extends one end of a gene segment); `Insertion` (sets the length of a junction) |
| Probe | measures the constructed sequence against the read, writes nothing | the error model, the counters. Not events today. |

The parameterization is the event's own axes in its tensor: `{n}` for a categorical event,
`{4, 4}` for a transition kernel. This is `inherent_shape()` in the current code; `size()` is
its product. Realizations are indexed `0..n-1` and the index is the tensor coordinate.

An event no longer carries a gene class, a side or a generated name. Its identity is its
position in a segment (§2.2). The current `GeneChoice_V_gene_Undefined_side_prio7_size89`
style of name has no equivalent.

### 2.2 Segment

A segment is an ordered composite of events that act on one piece of the receptor. It is what
the current code calls a seq_type.

| Segment | Events | Kind |
|---|---|---|
| `V` | `GeneChoice V`, `Deletion V 3'` | gene |
| `D` | `GeneChoice D`, `Deletion D 5'`, `Deletion D 3'` | gene |
| `J` | `GeneChoice J`, `Deletion J 5'` | gene |
| `VD` | `Insertion VD`, `DinuclMarkov VD` | junction |
| `DJ` | `Insertion DJ`, `DinuclMarkov DJ` | junction |

A gene segment has one Create event and at most one Modify event per end. A junction segment
has a length (Modify) and a content (Create); its content depends on its two neighbours and on
nothing else.

A segment exposes to its neighbours, for each of its ends: an offset (in the read during
inference, in the constructed sequence during generation) and the nucleotides at that end.
Every cross-segment need goes through this interface (§5.5).

### 2.3 Layout

The layout is the ordered list of segments, left to right.

```
VDJ       [V, VD, D, DJ, J]
VJ        [V, VJ, J]
tandem D  [V, VD1, D1, D1D2, D2, DJ, J]
```

Adjacency is position plus or minus one. A junction segment's two inputs are its two
neighbours in the layout. This replaces `@Seq_type_order`, `set_adjacent_segments()` on each
event, and the `(Event_type, seq_type, side)` key of `Events_map`.

### 2.4 Conditioning graph

A sparse directed acyclic graph over events. An edge `parent -> child` means the child's
distribution is conditioned on the parent's realization.

The shipped models, edges compacted (`Del` for `Deletion`, prio and size suffixes removed):

| Model | Edges |
|---|---|
| human TRB | `V->Del_V3'`, `V->D`, `V->J`, `D->Del_D3'`, `D->Del_D5'`, `J->Del_J5'`, `J->D`, `Del_D5'->Del_D3'` |
| mouse TRB | `V->Del_V3'`, `D->Del_D3'`, `D->Del_D5'`, `J->Del_J5'`, `J->D`, `Del_D5'->Del_D3'` |
| human IGK, human TRA | `V->J`, `V->Del_V3'`, `J->Del_J5'` |

Three properties follow from this table and fix the status of the graph:

- It is not derivable from the layout. `V->J` exists in every model and V and J are not
  adjacent; `V->D` exists in human TRB and not in mouse TRB.
- It varies per model and a user can edit it. The tandem D fixture ships with an empty edge
  section because that structure has not been decided.
- It determines the shape of every tensor and the layout of the EM accumulators.

The graph is therefore model data. What the code needs from it: the parents of each event, a
topological order, and an acyclicity check at load time. The edge-inversion and ancestor
machinery of the current `Topology` class is not needed at runtime.

### 2.5 Parameters

One tensor per event. Its shape is the concatenation of the parents' own axes, in topological
order, then the event's own axes: `[parent_1..., parent_k..., own...]`. In row-major storage
one parent combination is one contiguous row.

Two handler kinds, selected by the rank of the own axes and not by the event type: categorical
(rank 1) normalizes the last axis; Markov (rank 2) normalizes the last axis per from-state.
This is the existing `Model` and `Math` layer and it is unchanged by this proposal.

## 3. Model file

One JSON document. Human TRB, abbreviated:

```json
{
  "schema_version": 3,
  "layout": ["V", "VD", "D", "DJ", "J"],
  "segments": {
    "V":  { "create": "v_choice", "modify": { "3'": "v_3_del" } },
    "D":  { "create": "d_gene",   "modify": { "5'": "d_5_del", "3'": "d_3_del" } },
    "J":  { "create": "j_choice", "modify": { "5'": "j_5_del" } },
    "VD": { "length": "vd_ins", "content": "vd_dinucl" },
    "DJ": { "length": "dj_ins", "content": "dj_dinucl" }
  },
  "conditioning": {
    "j_choice": ["v_choice"],
    "d_gene":   ["v_choice", "j_choice"],
    "v_3_del":  ["v_choice"],
    "d_5_del":  ["d_gene"],
    "d_3_del":  ["d_gene", "d_5_del"],
    "j_5_del":  ["j_choice"]
  },
  "events": {
    "v_choice": { "type": "GeneChoice", "realizations": [ { "name": "TRBV5-8*01", "seq": "GAGG...", "index": 57 } ] },
    "v_3_del":  { "type": "Deletion",   "realizations": [ { "value": -4, "index": 0 } ] },
    "vd_ins":   { "type": "Insertion",  "realizations": [ { "value": 0, "index": 0 } ] },
    "vd_dinucl":{ "type": "DinucMarkov" }
  },
  "error_model": { "type": "SingleErrorRate", "rate": 0.000846 }
}
```

- `layout` replaces `@Seq_type_order`.
- `conditioning` replaces `@Edges`. It is keyed by nickname; generated names do not exist.
- An event's identity is its place in `segments`. The nickname is a handle, not an identity.
- v1 and v2 text files are read by an importer and are not written any more. The current
  `model_parms_to_json()` already produces the `events` and `conditioning` parts; `layout` and
  `segments` are derived from the seq_type order and from each event's seq_type and side.
- Marginals go in the same document or in a sibling file, as one array per event in the
  tensor layout of §2.5. This replaces the positional `model_marginals.txt`.

## 4. Generation

Generation runs the model forward in two loops.

1. Outer loop, over the conditioning graph in topological order: draw each event's realization
   from the row of its tensor selected by its parents' realizations. This covers gene choices,
   deletions and insertion lengths.
2. Inner loop, along the layout: build each gene segment from its template and its deletions;
   for each junction, draw the chain of the drawn length, seeded by the last nucleotide of the
   left neighbour; concatenate.

Draws are made from CDFs indexed by realization index. The current `draw_random_realization()`
walks an `unordered_map`, so the realization a seed selects depends on the standard library's
hash order; this is why the generation goldens only reproduce on Linux. The indexed sampler
removes that. It changes the output for a given seed once, and the goldens are regenerated
once.

What exists: `SamplingEngine` and its handlers perform step 1. What is missing: the assembler
of step 2. `FastGenerator` and `draw_random_realization()` are retired together.

## 5. Inference

### 5.1 Factorization

For one read `R` and the human TRB conditioning graph, the likelihood is

```
P(R) = sum_V sum_J sum_D   P(V) P(J|V) P(D|V,J)
       * sum_delV sum_delD5,delD3 sum_delJ
             L_V(V, delV)
           * G_VD(end_V, start_D, last_nt_V)
           * L_D(D, delD5, delD3) P(delD5|D) P(delD3|D, delD5)
           * G_DJ(end_D, start_J, last_nt_D)
           * L_J(J, delJ)
           * P(delV|V) P(delJ|J)
```

- `L_seg` is the local likelihood of a gene segment against the read: the match and mismatch
  counts of its aligned part after deletion, turned into a probability by the error model.
- `G_gap` is the probability of a junction: the insertion length, and the Markov chain over
  the nucleotides the read contains between the two offsets, seeded by the left neighbour's
  last nucleotide.

The outer sum runs over the clique of gene choices and is structured by the conditioning
graph. The inner sum, at fixed gene choices, is a chain along the layout: `V - VD - D - DJ - J`.
A chain is summed by a forward-backward pass over its boundary states. There is no
enumeration of scenarios and no pruning.

### 5.2 Segment tables

For each read and each gene segment, a table indexed by (allele alignment, deletion
realizations) holding: the offset of each end in the read, the nucleotides at each end, and
the match and mismatch counts of the aligned part.

The table depends on the read and the alignments. It does not depend on the model parameters.
It is computed once per read and reused by every EM iteration. With `SingleErrorRate`, the
local likelihood is `(1-e)^matches * e^mismatches`, recomputed from the stored counts when `e`
changes.

For each junction and each pair of neighbour offsets within the read, the inserted sequence is
the read's content between the two offsets. The chain probability is a product over that
sequence, recomputed per iteration from the Markov matrix. The nucleotide path is fixed per
read.

### 5.3 Assembly

The outer loop runs over the tuples of gene choices that have alignments in the read. The
inner pass is a forward-backward over boundary states along the layout. A boundary state is
(offset, boundary nucleotide); the nucleotide is needed to seed the next junction's chain.

The same pass serves three computations by changing the semiring:

| Semiring | Result | Replaces |
|---|---|---|
| (+, x) | likelihood and posterior marginals | the E-step of `iterate()` |
| (max, x) with backtracking | best scenarios | `best_scenarios` output |
| (+, x) with extended state | counters (coverage, Pgen per scenario) | the `Counter` hierarchy |

Cost, human TRB, one read, one EM iteration. Gene choices with alignments: about 10 V, 4 J, 3 D,
120 tuples. Inner chain: 21 V deletions, 21 x 21 D deletions, 23 J deletions, about 20,000
operations per tuple. About 2.5 million operations per read. The segment tables are not
recomputed.

The result is exact. The current implementation prunes with `infer.likelihood_threshold` and
`infer.probability_ratio_threshold`; the two heuristics disappear.

### 5.4 EM

E-step: the posterior marginals from the forward-backward pass are accumulated into one
accumulator per event, at the index of the realization, in the tensor layout. M-step: each
handler normalizes its last axis. This is the existing `maximizeLikelihood()`.

The deterministic reduction of #68 keeps its structure: one accumulator per read, merged into
chunk accumulators in a fixed order, chunks merged in order after the parallel region.

### 5.5 Segment boundaries

Everything a segment needs from a neighbour crosses the interface of §2.2: an offset and the
boundary nucleotides. Today one nucleotide is needed, to seed a Markov chain. Context-dependent
error models (`Hypermutation_full_Nmer_errorrate`, used by the `bcr_heavy` supplementary
models) read a window of `k` nucleotides that can straddle a boundary; for them the table
exposes `k` boundary nucleotides and the assembler evaluates the window at the join. The
information is available in the tables; the width of the interface is a parameter of the
error model. This is addressed when those models are ported, and does not change the design.

### 5.6 Threads

The segment tables are immutable once built. The forward-backward state is local to the read.
Events hold no per-read state. The per-thread copy of `Model_Parms`, the `copy()` of every
event and the re-run of `finalize()` on the copies are not needed: threads share the model
read-only and own their read's accumulator.

## 6. Effect on the code base

Line counts on `integration/model-on-tandemD`, `src/igor`, `.cpp` + `.h` + `.tpp`.

| Area | Lines | Under this proposal |
|---|---:|---|
| `Core`: `iterate()` bodies of the four events and `Rec_Event` | 4,665 | removed; replaced by segment tables |
| `Core`: exploration context, spans, offsets, `SafetyMatrix`, `LayeredArray`, `Index_map` | 1,793 | removed |
| `Core`: `GenModel` (EM loop, legacy generation) | 1,337 | replaced by the assembler and a smaller EM loop |
| `Core`: `Model_Parms`, `Model_marginals`, text I/O, migration | 3,524 | replaced by the JSON reader and the importer |
| `Core`: `FastGenerator` | 1,067 | removed |
| `Core`: error models and error counters | 5,367 | unchanged now; become Probe events later |
| `Core`: aligner | 3,490 | unchanged |
| `Core` total | 29,569 | |
| `Model`: handlers, engines, handler factories | 2,214 | kept |
| `Model`: `Topology`, `EventFactory`, `LegacyBridge` | 618 | `Topology` shrinks to a parents table; `EventFactory` kept; `LegacyBridge` becomes the importer |
| `Math` | | kept |

New code: the segment tables, the assembler with its three semirings, the generation
assembler, the JSON reader for `layout` and `segments`.

## 7. Numerical consequences and validation

The current inference is a pruned enumeration. The assembler is exact. For the same model and
the same reads the likelihoods differ, the exact one being larger or equal, and the inferred
parameters differ slightly. The bitwise regression gate on the inference goldens cannot be kept
across this change; this is a consequence of the design, not a defect of the implementation.

Validation replaces it with three checks:

1. The assembler against the current code run exhaustively, on TCR models with
   `SingleErrorRate`: `infer.likelihood_threshold = 0`, `infer.probability_ratio_threshold = 0`,
   a few reads, per-read log-likelihood equal to floating-point precision.
2. Inference recovery: generate from a known model, infer, compare parameters. The
   `test_ModelInference` cases are this kind of check.
3. Generation goldens regenerated once with the indexed sampler; identical across platforms
   from then on.

Comparisons are at tolerance for numbers and order-insensitive for realization lists. The
regression suite then runs on every platform. Today it only passes on Linux, for three reasons:
the `unordered_map` iteration in the sampler and in the text writers, and the 80-bit
`long double` of x86-64.

## 8. Order of work

1. Prototype the chain of §5.3 standalone: human TRB, `SingleErrorRate`, about ten reads,
   compared with the exhaustive run of check 1 above. This decides whether the rest is
   engineering. A few hundred lines.
2. `RecombinationModel` owns the conditioning table, the layout (the current `SeqTypeRegistry`),
   the tensors and the error model; events are resolved in its constructor. This is the
   decision taken on 2026-10-01 (registry in the model, resolution in the constructor, error
   model as a member).
3. Segment tables: gene segments, then junctions.
4. Assembler: the (+, x) pass and the E-step; then (max, x) and best scenarios; then counters.
5. Generation assembler over `SamplingEngine`; retire `FastGenerator` and
   `draw_random_realization()`; regenerate the generation goldens.
6. Model file v3 with `layout`, `segments`, `conditioning`; v1 and v2 readers kept as importers.
7. Probe role: the error model as an event; the boundary window for N-mer models.

## 9. Open points

- Priority. It is the only encoding of exploration order today. Under this design the outer
  order is the conditioning graph's topological order and the inner order is the layout, so
  priority is not needed for correctness. Whether it is kept in the file as data is to be
  decided.
- Parameter tying. The tandem D model shares one allele distribution between D1 and D2. The
  conditioning table needs a way to say that two events share a tensor, and the M-step must
  sum their accumulators before normalizing.
- Per-scenario outputs. Which counters are kept, and in which form the best scenarios are
  reported, decides the extended state of the third semiring.
- Negative deletions. A palindromic extension changes a segment's end content. The segment
  table stores the extended content at that end; to be confirmed when the tables are
  specified.
