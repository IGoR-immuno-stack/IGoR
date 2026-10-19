# Implementation graph: every remaining step, across the plans

**Created**: Oct 19 2026
**Purpose**: one picture to pilot the order of work. It holds every step still to do in:
- [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](REC_EVENT_CAPABILITY_REFACTORING_PLAN.md) (CAP,
  tandem D);
- [ERROR_MODEL_REWRITE_PLAN.md](ERROR_MODEL_REWRITE_PLAN.md) (ERR);
- [ENGINE_MIGRATION_PLAN.md](ENGINE_MIGRATION_PLAN.md) (MIG);
- [GENERATION_REWRITE_PLAN.md](GENERATION_REWRITE_PLAN.md) (GEN): only its §6 is left;
- [ITERATE_GENERIC_REWRITE_PLAN.md](ITERATE_GENERIC_REWRITE_PLAN.md) (ITER): R4, the optional
  S4d, the tier work and its §9;
- the steps of [ARCHITECTURE_SYNTHESIS.md](ARCHITECTURE_SYNTHESIS.md) (SYN) that come after the
  migration;
- what the hand-over of [doc/LAYER_REFACTORING_PROPOSAL.md](../doc/LAYER_REFACTORING_PROPOSAL.md)
  (LAYER) leaves;
- the landing of AA Pgen.

The plans hold the content of each step; this document holds only the order.
**Status**: as of Oct 19 2026, `origin/feature/tandemD` @ `814442b`.

---

## 1. Reading the graph

```mermaid
flowchart LR
  a1[step] ==>|hard block| a2[step]
  b1[step] -.->|soft block| b2[step]
  c1[step] -->|natural order| c2[step]
  m{{milestone}}
  q{decision}
  r[delivered, awaiting review]:::review
  g[moves golden outputs]:::golden
  o[optional]:::optional
  classDef review fill:#dbe8fb,stroke:#1a5fb4,color:#000
  classDef golden fill:#ecdcf7,stroke:#7b2fb5,color:#000
  classDef optional stroke-dasharray:5 5
  classDef milestone fill:#fde7c8,stroke:#c26100,stroke-width:2px,color:#000
  classDef decision fill:#fbdcdc,stroke:#b42318,color:#000
  class m milestone
  class q decision
```

- **Hard block** (thick): the step cannot be done, or cannot stay bitwise, before the other lands.
- **Soft block** (dotted): it can, at the cost of redoing work or of a golden move inside a
  bitwise block.
- **Natural order** (plain): no dependency, only the order the decisions or the file sets make
  sensible. For example, everything after milestone 1 waits for it because tandem D comes first,
  not because of the code.
- **Milestones** (hexagons) are the states the plans name: milestone 1 and 2, and the migration's
  checkpoints into `develop`.
- **Decisions** (diamonds) are the open questions that gate a step.

The subgraphs are the migration's lanes (MIG §4): steps in one lane share files and run in
sequence; steps in different lanes run in parallel.

---

## 2. The graph

```mermaid
flowchart LR

  subgraph REV["Delivered, awaiting review"]
    PRT1["T1 · a tandem-D model loads<br/>feature/tandemD-loading"]
    PRT3["T3 · TraversalPlan, SegmentIndex, FlankCheck<br/>feature/tandemD-T3-plan"]
    PRDOC["Plans: T7, error models, migration, this graph<br/>feature/tandemD-T7-pinning"]
  end

  subgraph TAN["Tandem D, milestone 1 · CAP"]
    T2["T2 · distinct generated names"]
    T7a["T7a · Single over the ordering;<br/>hypermutation and counter guards"]
    T5a["T5a · V-D1-D2-J built from human TRB, in v2"]
    T5b["T5b · reads generated and frozen"]
    T5["T5 · round trip; cost against TRB;<br/>mayfail off the tandem cases"]
    M1{{"Milestone 1<br/>tandem D on the legacy DFS"}}
    DT5{"T5 timings:<br/>is the legacy DFS fast enough?"}
  end

  subgraph LT["Tests first · lane T"]
    L0["L0 · full-precision A/B oracle; VJ track;<br/>driver pins; bridge parity"]
    CT0["CT0 · counter characterization suite"]
    E0["E0 · Single characterization suite"]
    H0["H0 · hypermutation regression track"]
    H1["H1 · hypermutation unit suite"]
  end

  subgraph LE["Error models · lane E · ERR"]
    E1["E1 · SingleErrorModel core;<br/>legacy adapter"]
    ER1["ER1 · pooled error rate, X13"]
    E2["E2 · gating and viterbi out of the error model"]
    E3["E3 · mismatch role; segment-aware bound"]
    E4["E4 · error node; bridge carries r;<br/>JSON reader, X20"]
    T7b["T7b · hypermutation models by name;<br/>the pinning goes"]
    H2["H2 · one context model,<br/>two parameterizations"]
    QE7{"ERR Q7: flank ends<br/>in the TraversalPlan?"}
    H3["H3 · flank segments, B3"]
    H4["H4 · scope sets"]
    HR["HR · sound bounds and repairs"]
    H5["H5 · hypermutation handlers"]
  end

  subgraph LC["Counters · lane C · MIG §2"]
    CR3["CR · CX3, J coverage past its array"]
    CT1["CT1 · dead code out"]
    CR["CR · counter repairs"]
    CT2a["CT2a · snapshot by id; lookups by role"]
    CT3a["CT3a · tandem: best scenarios, errors, Pgen"]
    CT2b["CT2b · scenario record and provenance"]
    CT3b["CT3b · tandem: coverage, D included"]
    CRP["CR-P · palindromes at the mirrored position"]
    CT4["CT4 · counter lifecycle in the driver"]
  end

  subgraph LD["Driver and engine · lane D · one open PR at a time"]
    L1["L1 · EM driver and LegacyDfsEngine"]
    L2["L2 · leaf and dispatch to the engine"]
    CP2{{"Checkpoint into develop:<br/>driver and engine"}}
    QD{"MIG Q3 to Q5: labelled axes,<br/>legacy queue hint, uid order"}
    L3["L3 · uid and scenario record"]
    L4a["L4a · tensor accumulators, ordered merge"]
    Q7{"MIG Q7: hypermutation<br/>through an adapter node?"}
    L4b["L4b · M-step by the handlers"]
    L4c["L4c · the DFS reads the tensors"]
    L5["L5 · the engine reads RecombinationModel"]
    CP5{{"Checkpoint into develop:<br/>engine on RecombinationModel"}}
    L6["L6 · tiers: 0a per run, 0b per iteration"]
    Q1516{"MIG Q15, Q16: events promoted;<br/>where the contexts go"}
    L8["L8 · the legacy types go; LAYER stage 2"]
    CP8{{"Migration complete:<br/>checkpoint into develop"}}
  end

  subgraph LM["Model module · lane M"]
    MAX["Labelled tensor axes"]
    MACC["AccumulatorSet: empty_copy,<br/>ordered merge, chunk slots"]
    MCUT["Handler M-step at the legacy rule;<br/>fixed flag; exact parity tests"]
    MREG["Registry, ordering and layout<br/>owned by RecombinationModel"]
    MIMP["Importer v1/v2 to RecombinationModel;<br/>topology_to_json"]
  end

  subgraph LG["Generation · lane G"]
    GMK["Markov sampler: anchor nucleotide;<br/>from-state in run()"]
    QPR{"Sampling in double<br/>or long double?"}
    QE5{"ERR Q5: errors on<br/>insertions in generation?"}
    L7["L7 · sampling driver and engine;<br/>SequenceGenerator retired"]
    R4["R4 · dinucleotide matrix<br/>in initialize_event"]
    QFG{"FastGenerator:<br/>keep or retire?"}
    T4b["T4b · FastGenerator by SeqTypeId"]
    T6b["T6b · gene_to_seqtype_migr deleted"]
    LEG["Gene_class_legacy in the v1 I/O only;<br/>Utils.h dispersed"]
  end

  subgraph POST["After the migration"]
    PHC["Phase C · validation and admissibility"]
    SYN4b["SYN 4b · sectioned model document;<br/>max_digits10 writer"]
    TOOLS["fmt, validate, overlays;<br/>accumulators section"]
    SYN5["SYN 5 · layout DAG with guards"]
    SYN7["SYN 7 · pure semantics, prepare(read),<br/>evaluate, tier-1a hoists"]
    DFSB["DFS-B prototype,<br/>measured against DFS-A"]
    PHA["Rest of Phase A: footprint reads,<br/>candidates, arity"]
    SYN9["SYN 9, Phase D · ChainDPEngine"]
    DPB["DpBoundedDfs, FFBS,<br/>stochastic EM"]
    KB["k-best scenarios on the DP"]
    CT5["CT5 · counter classes; Observer"]
    TIE["D1/D2 parameter tying;<br/>B9 step 3"]
    COV["Covariate nodes"]
    PPE["Error rate per parent realization"]
    EVM["Events_map keyed by capability"]
    S4d["S4d · tensor-backed no_d_align structure"]
    PERF["Container performance, CAP D4"]
    QB10{"B10 representation:<br/>CAP options b or c"}
    INTU{"int_undefined:<br/>drop the state or keep it?"}
    B10["B10 · optional D2"]
    GD6["GEN D6 · empty chain anchor"]
    M2{{"Milestone 2<br/>tandem D on real repertoires"}}
    Q13{"MIG Q13: AA Pgen<br/>before L2 or after L5?"}
    AA["AA Pgen lands"]
    NC["Non-canonical biology,<br/>VDJ review"]
  end

  %% Review and the milestone-1 path
  PRT1 --> PRDOC
  PRT1 ==> T2
  PRT1 ==> T7a
  PRT3 ==> T7a
  PRT1 ==> T5a
  T2 ==> T5a
  T5a ==> T5b
  T5b ==> T5
  T7a ==> T5
  PRT3 ==> T5
  T5 ==> M1
  T5 ==> DT5

  %% Tandem first: the rest starts after milestone 1
  M1 --> L0
  M1 --> CT0
  M1 --> E0
  M1 --> H0
  M1 --> MCUT
  M1 --> MREG

  %% Error models, Part I
  E0 ==> E1
  E1 ==> ER1
  ER1 --> E2
  E2 --> E3
  PRT3 ==> E3
  E1 ==> E3
  E3 --> E4
  E1 ==> E4
  ER1 ==> E4

  %% Error models, Part II
  H0 --> H1
  H0 ==> T7b
  H1 ==> T7b
  T7a ==> T7b
  H1 ==> H2
  T7b --> H2
  H2 --> H3
  PRT3 ==> H3
  QE7 ==> H3
  H2 ==> H4
  E3 -.-> H4
  H1 ==> HR
  H2 --> HR
  H2 ==> H5
  HR --> H5
  E4 --> H5
  L4b ==> H5

  %% Counters
  CT0 ==> CR3
  CT0 ==> CT1
  CT0 ==> CR
  E2 -. "CX13" .-> CR
  CT0 ==> CT2a
  PRT3 ==> CT2a
  T7a --> CT2a
  CT2a ==> CT3a
  E3 ==> CT3a
  CT2a ==> CT2b
  L3 ==> CT2b
  CT2b ==> CT3b
  CT3a --> CT3b
  CT3b ==> CRP
  L1 ==> CT4
  L4a ==> CT4

  %% Driver and engine
  L0 ==> L1
  CT0 ==> L1
  H0 == "hypermutation models" ==> L1
  PRT3 ==> L1
  T7a ==> L1
  ER1 -.-> L1
  L1 ==> L2
  E2 -.-> L2
  L2 ==> CP2
  L2 ==> L3
  T2 ==> L3
  QD ==> MAX
  MAX ==> L3
  L3 ==> L4a
  L1 --> MACC
  MACC ==> L4a
  L4a ==> L4b
  MCUT ==> L4b
  E4 ==> L4b
  Q7 ==> L4b
  L4b ==> L4c
  L4c ==> L5
  MREG ==> L5
  MREG ==> MIMP
  E4 ==> MIMP
  MIMP ==> L5
  T7b -.-> L5
  L5 ==> CP5
  L5 ==> L6
  L5 ==> L8
  L6 --> L8
  L7 ==> L8
  CT3b ==> L8
  CT4 ==> L8
  H5 -. "or keep the Q7 adapter" .-> L8
  T6b -.-> L8
  Q1516 ==> L8
  L8 ==> CP8

  %% Generation
  GMK ==> L7
  QPR ==> L7
  QE5 -.-> L7
  L0 ==> L7
  E1 ==> L7
  E4 ==> L7
  T2 ==> L7
  L7 ==> R4
  L7 --> QFG
  QFG -- "kept" --> T4b
  T4b ==> T6b
  QFG ==> T6b
  T7b ==> LEG
  T6b ==> LEG
  LEG -.-> L8

  %% After the migration
  L1 -.-> PHC
  L5 ==> PHC
  L5 --> SYN4b
  SYN4b ==> TOOLS
  MACC ==> TOOLS
  L5 -.-> SYN5
  L5 ==> SYN7
  L6 --> SYN7
  SYN7 ==> DFSB
  SYN7 --> PHA
  SYN7 ==> SYN9
  SYN9 ==> DPB
  SYN9 ==> KB
  SYN9 ==> CT5
  L4b ==> TIE
  L4b ==> COV
  L5 ==> COV
  L5 ==> PPE
  E4 ==> PPE
  L5 ==> EVM
  L4c -.-> S4d
  L6 --> PERF
  DT5 --> QB10
  DT5 ==> B10
  QB10 ==> B10
  INTU --> B10
  SYN5 -.-> B10
  L3 -. "if the DFS suffices" .-> B10
  SYN9 -. "if the DFS is too slow" .-> B10
  B10 ==> GD6
  B10 ==> M2
  Q13 ==> AA
  L5 -. "or before L2 opens" .-> AA
  SYN5 -.-> NC
  SYN9 -.-> NC

  classDef review fill:#dbe8fb,stroke:#1a5fb4,color:#000
  classDef golden fill:#ecdcf7,stroke:#7b2fb5,color:#000
  classDef optional stroke-dasharray:5 5
  classDef milestone fill:#fde7c8,stroke:#c26100,stroke-width:2px,color:#000
  classDef decision fill:#fbdcdc,stroke:#b42318,color:#000
  class PRT1,PRT3,PRDOC review
  class ER1,HR,CR3,CR,CRP,L7,SYN4b golden
  class S4d,PERF,NC optional
  class M1,M2,CP2,CP5,CP8 milestone
  class DT5,QE7,QE5,QD,Q7,Q1516,QPR,QFG,QB10,INTU,Q13 decision
```

---

## 3. The paths that matter

**To milestone 1.** The longest chain is T1's review → T2 → T5a → T5b → T5. T3's review and T7a
join it at T5. Nothing in the error-model plan or the migration is on it.

**To the migration's end.** The chain is lane D: L0 → L1 → L2 → L3 → L4a → L4b → L4c → L5 → L8.
Each other lane must arrive in time:
- **lane E**: E0 → E1 → ER1 → E4, before L4b;
- **lane M**: the handler M-step (MCUT) and the `AccumulatorSet` (MACC) before L4; the registry
  and the importer (MREG, MIMP) before L5;
- **lane G**: L7, before L8;
- **lane C**: CT3b and CT4, before L8;
- **hypermutation**: H0 → H1 → T7b, before L5, softly; H5 before L8, or the adapter is kept (Q7).

ER1, the one golden move on the way, lands before L1 so that every lane-D step stays bitwise.

**To real tandem data (milestone 2).** It runs through T5's timings:
- **If the legacy DFS is fast enough:** B10 can follow L3, where adjacency per scenario has a home.
- **If it is too slow, the expected outcome:** B10 waits for `ChainDPEngine`, after L5 and the
  synthesis' step 7.

B10's representation (CAP B10, options (b) and (c)) is decided at the same time.

**Steps that can start right after milestone 1, in parallel.**
- Lane T, the tests: L0, CT0, E0, H0.
- Lane M: MCUT, MREG.
- In lane C, the CX3 repair as soon as CT0 confirms the defect: it can corrupt memory on long reads.

---

## 4. The nodes

Status: **review** — delivered on its branch, awaiting merge · **pending** — ready once its hard
blocks land · **gated** — waits on a decision or on a step of another plan · **optional**.

### Delivered, awaiting review
| Node | Source | Hard blocks | Status |
|---|---|---|---|
| T1 | CAP *T1* | — | review; conflicts with T3 in `Rec_Event.cpp` (both add functions at one spot: keep both) |
| T3 | CAP *T3* | — | review |
| Plans | CAP, ERR, MIG, this document | — | review; the branch contains T1 |

### Tandem D, milestone 1
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| T2 | CAP blocker T2 | T1 | — | pending |
| T7a | CAP *T7*, ERR §4.1 | T1, T3 | — | pending |
| T5a | CAP *T5* | T1, T2 | — | pending; the fixture script, the v2 model, the marginals filled into IGoR's own layout |
| T5b | CAP *T5* | T5a | — | pending |
| T5 | CAP *T5* | T5b, T7a, T3 | — | pending |
| Milestone 1 | CAP | T5 | — | — |
| T5 timings | MIG Q14 | T5 | — | decision |

### Lane T, tests first
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| L0 | MIG §3 | — | milestone 1 (priority) | pending |
| CT0 | MIG §2.6 | — | milestone 1 | pending |
| E0 | ERR §4.2 | — | milestone 1 | pending |
| H0 | ERR §4.3 | — | milestone 1 | pending |
| H1 | ERR §4.3 | — | H0 | pending |

### Lane E, error models
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| E1 | ERR §4.2 | E0 | — | pending |
| ER1 | ERR §4.2, X13 | E1 | — | pending; moves every inference golden once |
| E2 | ERR §4.2 | — | after ER1 | pending |
| E3 | ERR §4.2 | T3, E1 | after E2 | pending |
| E4 | ERR §4.2, X20 | E1, ER1 | after E3 | pending |
| T7b | CAP *T7*, ERR §4.3 | H0, H1, T7a | — | gated on H0, H1 |
| H2 | ERR §4.3 | H1 | after T7b | pending |
| H3 | ERR §4.3, CAP B3 | T3, ERR Q7 | after H2 | gated on Q7 |
| H4 | ERR §4.3 | H2 | E3 | pending |
| HR | ERR §4.3 | H1 | after H2 | pending; moves the hypermutation goldens |
| H5 | ERR §4.3 | H2, L4b | after HR, E4 | gated on L4b |

### Lane C, counters
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| CR, CX3 | MIG §2.3 | CT0 | — | pending; repair early |
| CT1 | MIG §2.6 | CT0 | — | pending |
| CR | MIG §2.6 | CT0 | E2 (CX13) | pending; each repair names the outputs it moves |
| CT2a | MIG §2.6, T7b's counter half | CT0, T3 | after T7a | pending |
| CT3a | MIG §2.5 | CT2a, E3 | — | pending |
| CT2b | MIG §2.6 | CT2a, L3 | — | gated on L3 |
| CT3b | MIG §2.5 | CT2b | after CT3a | gated on L3 |
| CR-P | MIG §2.5 | CT3b | — | moves the coverage goldens once |
| CT4 | MIG §2.6 | L1, L4a | — | gated on L4a |

### Lane D, driver and engine
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| L1 | MIG §3 | L0, CT0, T3, T7a; H0 for hypermutation models | ER1 | pending |
| L2 | MIG §3 | L1 | E2 | pending |
| L3 | MIG §3 | L2, T2, labelled axes | — | gated on Q3–Q5 |
| L4a | MIG §3 | L3, AccumulatorSet | — | pending |
| L4b | MIG §3 | L4a, handler M-step, E4, Q7 | — | gated on Q7 |
| L4c | MIG §3 | L4b | — | pending |
| L5 | MIG §3 | L4c, registry ownership, importer | T7b | pending |
| L6 | MIG §3 | L5 | — | pending; benchmarked |
| L8 | MIG §3 | L5, L7, CT3b, CT4, Q15/Q16 | L6, H5 (or the adapter), T6b, the enums | gated |
| Checkpoints | MIG §4 | L2, L5, L8 | — | — |

### Lane M, Model module
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| Labelled tensor axes | MIG Q3, REV §15.4 | Q3–Q5 | — | gated |
| AccumulatorSet | SYN step 8, MERGE §6.2 | — | after L1, which defines the interface | pending |
| Handler M-step at the legacy rule | MIG §1.1, §3 | — | milestone 1 | pending |
| Registry, ordering, layout | MERGE §4.2, SYN §5 | — | milestone 1 | pending |
| Importer, `topology_to_json` | SYN step 4, LAYER §6 | registry ownership, E4 | — | pending |

### Lane G, generation
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| Markov sampler | GEN §6 | — | — | pending |
| L7 | MIG §3, GEN §6, SYN step 1 | L0, E1, E4, T2, the Markov sampler, the precision choice | ERR Q5 | pending; moves the `generate` goldens once |
| R4 | ITER §6.9 | L7 | — | gated on L7 |
| T4b | CAP *T4* | FastGenerator kept | — | gated on the decision |
| T6b | CAP *T6* | T4b, or FastGenerator retired | — | gated |
| `Gene_class_legacy`, `Utils.h` | CAP, LAYER hand-over | T7b, T6b | — | gated |

### After the migration
| Node | Source | Hard blocks | Soft | Status |
|---|---|---|---|---|
| Phase C | CAP, MIG §7 | L5 | L1 | gated |
| SYN 4b | SYN step 4b | — | after L5, to keep its golden move outside | moves the marginal-file goldens once |
| fmt, validate, overlays | REV §14–15 | SYN 4b, AccumulatorSet | — | gated |
| SYN 5 | SYN step 5 | — | L5 | gated |
| SYN 7 | SYN step 7, ITER §2.10 | L5 | after L6 | gated |
| DFS-B | SYN step 7, REV §10.3 | SYN 7 | — | gated |
| Rest of Phase A | CAP | — | after SYN 7 | gated |
| SYN 9, Phase D | SYN step 9, CAP | SYN 7 | — | gated |
| DpBoundedDfs, FFBS | REV §5.3 | SYN 9 | — | gated |
| k-best on the DP | REV §5.1 | SYN 9 | — | gated |
| CT5 | MIG §2.6 | SYN 9 | — | gated |
| D1/D2 tying, B9 step 3 | CAP, SYN §4 | L4b | — | gated |
| Covariates | SYN, REV §6 | L4b, L5 | — | gated |
| Per-parent error rate | ERR §3.4 | L5, E4 | — | gated |
| `Events_map` by capability | ITER §9 | L5 | — | gated |
| S4d | ITER §6 | — | L4c | optional |
| Container performance | CAP D4 | — | after L6 | optional |
| B10 | CAP B10 | T5's timings, its representation | L3 if the DFS suffices; SYN 9 otherwise; SYN 5; the `int_undefined` decision | gated |
| GEN D6 | GEN §4 | B10 | — | gated |
| Milestone 2 | CAP | B10 | — | — |
| AA Pgen | CAP H1–H9, MIG Q13 | Q13 | after L5, or before L2 opens | gated |
| Non-canonical biology | VDJ review | — | SYN 5, SYN 9 | optional, research |

---

## 5. Keeping it current

- When a step lands, delete its node and its edges, and its row; record the landing in its plan.
- When a plan adds a step, add a node with its hard blocks first, then any soft blocks and its
  natural order, and give it a row.
- A soft block that is overruled — a step done before it, accepting the rework — is recorded in
  the step's plan, and the edge goes.
- A decision node goes once the decision is recorded in its plan; the steps it gated become
  pending.
