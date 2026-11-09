# Glossary

**Status**: reference documentation of the code on `feature/tandemD-summary` @ `4a79e8d` (Nov 9
2026), of the vocabulary of the plans in `docs/` at that commit, and of its equivalents in the AIRR
Standards 2.0.
**Audience**: anyone who reads or writes IGoR's code, its plans or its documentation, and needs to
know what a word means in the code, what it stands for on a receptor, which older names it
replaces, and what AIRR calls it.
**Companions**: [ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) and
[REC_EVENT_CONTRACT.md](REC_EVENT_CONTRACT.md) explain how the pieces named here work together;
[ITERATE_TEST_GUIDE.md](ITERATE_TEST_GUIDE.md) explains how they are tested. The AIRR schemas are
[Rearrangement](https://docs.airr-community.org/en/stable/datarep/rearrangements.html),
[Germline](https://docs.airr-community.org/en/stable/datarep/germline.html),
[Repertoire](https://docs.airr-community.org/en/stable/datarep/metadata.html) and
[Clone](https://docs.airr-community.org/en/stable/datarep/clone.html).

Each table gives, for a term:
- **Definition**: what the term means, in the words of other entries. A word in *italics* has an
  entry of its own. A definition names no identifier, file, link or plan step: a term that is
  itself an identifier is defined by its role. It ends with "Not a …" when another entry is easily
  confused with it.
- **Example**: the term on the human TRB model (`models/human/tcr_beta`) unless stated: a receptor
  sequence, an event and one of its values, or a value the code computes. On that model a receptor
  example and a code example are the same thing. "—" when there is none.
- **In the code**: the class, function, enum, file, JSON key or option that holds the concept, with
  a link to its header or reference section; its conventions (units, index base, layout in memory),
  its lifecycle (when it is built, by whom, how often) and its invariants. "—" when nothing in the
  code holds it yet; for a *(planned)* term, the names its plan gives may follow.
- **Legacy or former name**: the name the concept has in the legacy code or in a plan, or had
  before a rename, and the behaviour a step has since changed. For a term that is not wired yet,
  the code that holds the concept today.
- **AIRR**: the equivalent in the AIRR Standards 2.0 schemas, and where it differs. A bare field
  name is a field of the Rearrangement schema; a field of another object carries the object's name
  (`AlleleDescription.j_cdr3_end`). "—" when AIRR has none. Only the tables of receptors, model,
  segments, scenarios, outputs and generation (§1, §2, §3, §5, §9, §10) have this column: the
  others describe the engine and the target architecture, which AIRR does not cover.

A mark after a term gives its status:
- no mark: code that inference or generation runs today;
- *(unwired)*: merged, but called by neither inference nor generation (MIG §1.1);
- *(planned)*: only in a plan, with the step that brings it;
- *(dormant)*: merged and allocated, but nothing reads or writes it yet.

When a planned step will change a term's code, the term's *In the code* cell ends with a bold note
that names the plan and the step, such as "Subject to change — MIG L5.", and says in one sentence
what that step changes. Reference documents use the same note (§12.2). The step that makes the
change updates the row and drops the note.

Plans are cited by short name and step or section, as listed in §12.1. Where a plan and the code
name one thing differently, the entry uses the code's name. §13 lists the words that mean several
things, §14 the step prefixes that collide, and §15 the AIRR concepts IGoR has no term for.

**AIRR conventions.** AIRR coordinates count from 1 and include both ends: `v_sequence_start` and
`v_sequence_end` on the read, `v_germline_start` and `v_germline_end` on the template. IGoR's
*offsets* count from 0 and also include both ends. A Rearrangement annotates one read with one
scenario. Today only `igor-convert` reads and writes AIRR files, through the Streaming module, and
only reads and alignments. Inference and generation write IGoR's own files; CNT CT7 plans the best
scenarios as Rearrangements.

---

## 1. Receptors, and what IGoR models of them

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **gene** | A germline V, D or J gene. A *model* holds its *templates* in a *library*, and a *gene choice* picks one per *scenario*. A gene is not a *segment*: the two D segments of a *tandem D* model draw from one D library. | TRBV20-1, TRBD1, TRBJ2-7 | — | gene segment (in biology) | `v_call`, `d_call`, `d2_call`, `j_call`, each naming a gene with its allele. `c_call` names a constant gene, which IGoR does not model (§15) |
| **template** | A gene's nucleotide sequence: the value of a *gene choice*'s *realization*, read from the *library*. Each allele is a template of its own. | `TRBD1*01` = `GGGACAGGGGGC` | — | genomic template | `AlleleDescription` (Germline): `label` is the allele's name, `coding_sequence` its coding nucleotides; its `sequence` may also cover the leader, UTRs and RSS. AIRR says germline where IGoR says template. IGoR names a template by its FASTA header: in the shipped TRB files a full IMGT header, whose second field is the allele (`TRBJ1-1*01`) |
| **library** | The *templates* a *gene choice* draws from, and the read's alignments against them, keyed by *gene class*. | `TRBD1*01`, `TRBD2*01`, `TRBD2*02` | The `library` section of a *model document*; the `genomic{V,D,J}s.fasta` files of a model in the *v1 format*. | `Gene_class`, `event_class` | `GermlineSet` (Germline), whose `allele_descriptions` are the templates. `DataProcessing.germline_database` and `germline_set_ref` (Repertoire) name the set a repertoire was annotated with |
| **gene class** | Which *library* and which alignments a *gene choice* reads. It never identifies a *segment*. | `D_gene` for both D gene choices of a tandem model | `Gene_class {V_gene, D_gene, J_gene, Undefined_gene}` ([Utils.h](../../src/igor/Core/Legacy/Utils.h)), held as `Rec_Event::event_class`. **Subject to change — CAP T4b.** After it, the gene class moves off the base class onto `Gene_choice`. | `Gene_class_legacy`, whose junction classes `VD_genes`, `DJ_genes`, `VJ_genes`, `VDJ_genes` survive only in the *v1 format*'s reader and writer | `AlleleDescription.sequence_type`: `V`, `D`, `J`, or `C`, which IGoR does not model |
| **locus** | A chromosomal cluster of V, D, J and C genes. A model describes one. | TRB (`models/human/tcr_beta`) | No identifier of its own: the CLI's `model.chain` (`-chain`) picks a shipped model by its locus, and the model's directory names it. | chain | `locus`, with IMGT's names (IGH, IGK, IGL, TRA, TRB, TRD, TRG); `GermlineSet.locus`; `PCRTarget.pcr_target_locus` (Repertoire) |
| **trimming** | Nucleotides removed from one end of a gene before it is joined: a *deletion* with a positive value. | `d_3_del` = 3 turns `TRBD1*01` into `GGGACAGGG` | — | exonuclease deletion | No field of its own: the template outside `*_germline_start` and `*_germline_end`. A 5′ trim is `d_germline_start` − 1; a 3′ trim is the template's length minus `d_germline_end` |
| **P nucleotides** | A short inverted repeat of a gene's end, left by hairpin opening: a *deletion* with a negative value −k, which appends the k nucleotides nearest the end, reversed and complemented. An end has one deletion event, so it is either trimmed or extended, never both. | `d_3_del` = −2 turns `TRBD1*01` into `GGGACAGGGGGCGC` | — | palindromic insertion, negative deletion | `p3v_length`, `p5d_length`, `p3d_length`, `p5d2_length`, `p3d2_length`, `p5j_length`: a negative deletion's value with its sign changed. AIRR counts P nucleotides in the `np` regions, IGoR in the gene's *segment* |
| **N region** | Non-templated nucleotides added between two joined ends (by TdT). An *insertion* chooses how many, a *Markov chain* which ones. The stretch they occupy is a *segment* of its own. | N1 between TRBV and TRBD (`vd_ins`, `vd_dinucl`); N2 between TRBD and TRBJ (`dj_ins`, `dj_dinucl`) | — | N insertion; the v1 insertion classes `VD_genes`, `DJ_gene`, `VJ_gene` | `n1_length`, `n2_length`, `n3_length` count its nucleotides. `np1`, `np2`, `np3` hold them, with the P nucleotides on either side. A VJ model's N region is `np1` |
| **junction** (biology) | The joined region between two coding ends: their trimmed or extended ends and the N nucleotides between them. The code uses the word more narrowly; see *junction* (code) and §13. | TRBV…N1…TRBD…N2…TRBJ | — | — | `junction`, `junction_aa`, `junction_length`: the CDR3 with its two anchor codons. Not the code's *junction* |
| **CDR3**, **CDR3 anchors** | The CDR3 lies between two conserved codons, the CDR3 anchors: the V gene's Cys (IMGT 104) and the J gene's Phe or Trp (IMGT 118). Not a Markov chain's *anchor*. | junction `CASSLGQAYEQYF`, whose C and F are the anchors | `{V,J}_gene_CDR3_anchors.csv` give each anchor as the 0-based position of its codon's first nucleotide in the template. The anchors serve the CDR3 outputs only, `alignment.extract_cdr3` for reads (`CDR3SeqData`) and `generate.cdr3` for generated sequences: no event models them. | — | `cdr3`, `cdr3_aa`, `cdr3_start`, `cdr3_end` exclude the anchors; `junction` includes them. On the template, the V anchor is the codon before `SequenceDelineationV.cdr3_start`, and the J anchor the codon at `AlleleDescription.j_cdr3_end` (Germline), counted from 1. The `nt_CDR3` column of `generate.cdr3` spans both anchor codons: it is a `junction` |
| **tandem D** | Two D genes in one rearrangement. | the test model V-D1-D2-J built from human TRB (`tst/test_data/tandem_d/`); VDDJ rearrangements in IGH and TRD | A model whose *layout* holds two D *segments*, each with its own *gene choice* and *deletions*, and three *N regions*. Every segment is present in every scenario. **Subject to change — CAP B10.** It makes D2 optional. | V(DD)J, VDDJ | `d2_call` and the other `d2_` fields, `np3`, `n3_length`, `p5d2_length`, `p3d2_length`, `d2_frame`. Two D at most: the planned AIRR output refuses a third (CNT CT7) |
| **D-less rearrangement** | V joined to J with no D nucleotide left. | a TRB read with no D nucleotide | IGoR represents a D trimmed to nothing (an internal segment may be emptied, unlike an *end segment*), not a D absent from the scenario. **Subject to change — CAP B10.** It lets a D be absent from the scenario. | — | A null `d_call`, with `np1` running from V to J |
| **D reading frame** | Which of its three reading frames a D gene is read in. A scenario fixes it, through D's 5′ deletion and D's position relative to V's frame. | IGHD read in frame 1, 2 or 3 | No code computes it. | — | `d_frame`, `d2_frame` |
| **somatic hypermutation** | Point mutations of B-cell receptor genes after recombination. | IGHV mutated at AID hotspots | Modelled by an *error model* whose rate depends on the *N-mer context*: the *hypermutation models*. | — | No field: it shows as differences between `sequence_alignment` and `germline_alignment` inside the genes |
| **productive**, **non-productive** | A productive rearrangement can encode a receptor: in frame, with no stop codon. A non-productive one is not selected and reflects generation alone, so models are usually inferred from non-productive reads. | out-of-frame TRB rearrangements | Inference computes neither. `generate.cdr3` writes `is_inframe`, the junction's length being a multiple of 3, and leaves its `is_productive` column empty. | — | `productive` (AIRR recommends IMGT's definition), `vj_in_frame`, `stop_codon`. `is_inframe` is close to `vj_in_frame` |
| **repertoire** | The reads one inference runs on, or the sequences one generation writes. | the non-productive TRB reads of one donor | The code keeps no metadata about it: donor, sample or protocol. | — | `Repertoire` (Repertoire), with its `Subject`, `Sample` and `DataProcessing`; a generated one has `repertoire_type` `simulated` |
| **Pgen** | The probability that the model generates a nucleotide sequence: the sum of the probabilities of every *scenario* that produces it, without errors. | the Pgen of one TRB CDR3 nucleotide sequence | `Pgen_counter` estimates it for each read during inference. | — | — |
| **AA Pgen** *(planned, AA; merge hazards CAP H1–H9)* | The Pgen of an amino-acid sequence or motif, summed over all its codings. | `CASSLGQAYEQYF`, over every nucleotide sequence that encodes it | Partly merged and dormant: [GeneticCode.h](../../src/igor/Core/GeneticCode.h), `JournaledQuery` and the *pruning mismatch floor*. | — | —; the sequence it is asked for is often a `junction_aa` |

---

## 2. The model and its events

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **model** | The generative model: *events*, their *conditioning graph*, their *marginals*, the *layout* of its *segments* and an *error model*. | `models/human/tcr_beta/models/` | `Model_Parms` holds the structure and `Model_marginals` the parameters. **Subject to change — MIG L5.** The model becomes a `RecombinationModel`. | — | — |
| **event** | One random variable of the model, with its *realizations*, the *segment* it acts on and its *declarations*. | `v_choice`, `d_3_del`, `vd_ins`, `vd_dinucl` | An instance of a class derived from [`Rec_Event`](../../src/igor/Model/Legacy/Rec_Event.h), with `iterate()` for inference and `construct_realization()` for generation. Four types exist. **Subject to change — MIG L8.** The event keeps its definition and its `iterate()` moves to an *engine* (MIG Q15). | recombination event; node (SYN, the JSON schemas) | — |
| **event type** | The kind of an *event*: a *gene choice*, a *deletion*, an *insertion* or a *Markov chain*. | — | `EventType {GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t}` ([Types.h](../../src/igor/Core/Types.h)). In files, the tokens `GeneChoice`, `Deletion`, `Insertion`, `DinucMarkov`: the JSON `type` the *event factory* dispatches on. | `Event_type` | — |
| **gene choice** | The event that picks a *template* and places it on the read from one of its alignments, creating both *ends* and the *content* of its segment. When no alignment survives, the D gene choice tries every position (*exhaustive position scan*). | `j_choice` = `TRBJ2-7*01` | [`Gene_choice`](../../src/igor/Model/Legacy/Genechoice.h). | — | Its realization is `v_call`, `d_call`, `d2_call` or `j_call` (CNT CT7) |
| **deletion** | The event that moves the *end* of a segment that its *side* names, by its value, and trims or extends the content to match (*trimming*, *P nucleotides*). | `v_3_del`, `d_5_del`, `d_3_del`, `j_5_del`; `v_3_del` takes values −4 to 16 | [`Deletion`](../../src/igor/Model/Legacy/Deletion.h). | — | Read off `*_germline_start` and `*_germline_end` (*trimming*), and `p3v_length` and the like (*P nucleotides*) |
| **insertion** | The event that picks the length of an *N region* and places both *ends* of its segment, which stays *unfilled* until the *Markov chain* builds it. | `vd_ins` = 5: five N1 nucleotides | [`Insertion`](../../src/igor/Model/Legacy/Insertion.h). | until ITER R1 and R3, it also allocated the segment's content and recorded no offsets | `n1_length`, `n2_length`, `n3_length` |
| **Markov chain** | The event that creates the *content* of the *unfilled* segment an *insertion* placed, nucleotide by nucleotide, as a first-order Markov chain over A, C, G, T. It starts from the facing nucleotide of its *anchor* and runs away from it. Its marginals are a 4 × 4 transition matrix, normalized per previous nucleotide. | `vd_dinucl` starts from TRBV's last nucleotide and runs 5′ to 3′; `dj_dinucl` starts from TRBJ's first and runs 3′ to 5′ | [`Dinucl_markov`](../../src/igor/Model/Legacy/Dinuclmarkov.h). In inference, it prices the read's nucleotides in that segment. | dinucleotide Markov model; `DinucMarkov` | The N nucleotides of `np1`, `np2`, `np3` |
| **realization** | One value of an event. Its **index** is its coordinate in the *marginals*. | `TRBD1*01`; `d_3_del` = −2; `vd_ins` = 5; a nucleotide, for a Markov chain | `Event_realization`: a name, an integer or string value, and the index, dense from 0. | `label` and `index` in the *model document* | The value of a `*_call`, or of a length field |
| **nickname** | The event's identity, unique in its model: it names the event in model files and outputs. | `d1_5_del` | Required (CAP T2). It keys `Model_Parms::edges`, the marginals' index and stride maps and `processed_events`. **Subject to change — MIG L3.** The *uid* replaces it as the key of those maps. | the *generated name* | — |
| **generated name** | The name IGoR built from an event's type, *gene class*, *side*, *priority* and size before nicknames were required. It is the same for the two D blocks of a tandem model. | `GeneChoice_V_gene_Undefined_side_prio7_size89` | `get_v1_legacy_name()`. Written to v1 files, and to the outputs with `output.event_names = legacy`. | `name`, `get_legacy_name()` | — |
| **event identifier** | A dense integer per event, which indexes the *index map* and the *next-event array*. | — | Set by `Model_Parms` (`get_event_identifier()`). `copy()` keeps it, which is how one *traversal plan* serves every thread's copy of the events. | `event_index` | — |
| **uid** *(unwired)* | A dense integer per event, given in the order the events are added. | — | Set by `Topology`. **Subject to change — MIG L3.** It gives the legacy events one, in nickname order, and keys their maps by it. | *event identifier*; the *nickname* as a key | — |
| **priority** | An event's rank in the *model queue*: higher runs first. | `v_choice`, `j_choice` 7; `d_gene` 6; the deletions 5; `vd_ins` 4; `vd_dinucl` 3; `dj_ins` 2; `dj_dinucl` 1 | In a JSON model it is an engine hint, and `engine_hints.legacy_dfs.event_order` sets the order of the event list that breaks ties. | — | — |
| **model queue** | The order the legacy engine runs the events in: at each step, the highest-priority event whose parents are all queued, ties broken by the order of the event list. The engine's choice rather than the model's. Not the *ordering*. | the tandem model's *joining order* | `Model_Parms::get_model_queue()`. **Subject to change — MIG L5.** Each planned engine chooses its own. | event order | — |
| **joining order** | An event order that follows the order the segments join: D1 to D2, then D1D2 to J, then V, each joint's deletions right before its insertion. | D1 and D2 joined first | The event order of the tandem test model (CAP T5). | — | — |
| **conditioning graph** | The DAG of events in which an edge parent → child conditions the child's distribution on the parent's realization. Each parent adds an axis to the child's *marginals*. | `v_choice` → `j_choice`; `v_choice`, `j_choice` → `d_gene`; `d_gene` → `d_3_del`; `d_5_del` → `d_3_del` | Held in `Model_Parms::edges` and the JSON `structure.conditioning`; `Topology` holds it on the unwired side. | v1 `@Edges`; "topology" in older plans (MERGE §0.1) | — |
| **marginals** | The model's parameters. For each event, one row per configuration of its parents' realizations, holding one probability per realization. | P(TRBD gene given TRBV and TRBJ) | `Model_marginals` keeps them in one flat `long double` array, each event's own realizations last. **Subject to change — MIG L4.** They become one tensor per event. | `model_marginals.txt` | — |
| **index map**, **stride** | The index map holds where each event's current marginal row starts. A parent moves each child's row by its realization's index times the child's stride. | — | `Index_map`, a `LayeredArray<size_t>` keyed by *event identifier*. The strides are `memory_and_offsets`, from `Model_marginals::get_offsets_map()`. | `base_index_map`; the stride is called an "offset" in the code (§13) | — |
| **fixed** | An event whose marginals inference does not update. | — | `Rec_Event::fix()`: `add_to_marginals()` is not called on it. **Subject to change — MIG L4b.** Every node gets such a flag. | — | — |
| **error model** | What prices the mismatches between a scenario's constructed sequence and the read at the *leaf*, and bounds that price for *pruning*. | sequencing errors in TRB reads; hypermutation in IGH | [`Error_rate`](../../src/igor/Model/Legacy/Errorrate.h) and its subclasses (§9). **Subject to change — ERR, SYN step 8.** It becomes a *Probe event*. | error rate (the code's class name, and Single's parameter r) | — |
| **`Model_Parms`** | The legacy model structure: the events, the edges, the error model and the *segment registry*, made ready by *finalize()*. | — | [Model_Parms.h](../../src/igor/Model/Legacy/Model_Parms.h). **Subject to change — MIG L5.** Once it lands, the class is only the v1 importer. | — | — |
| **`finalize()`** | The step that makes a model ready: it runs once per model, after the events are added. | — | `Model_Parms::finalize()`. It refuses two events of the same type on the same segment and side, pins the *legacy ids*, infers a missing *ordering*, freezes the *registry*, resolves each event's *SegmentId* and tells it its *neighbours* ([ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) §2.2). | — | — |
| **`Events_map`** | The legacy lookup of events by (event type, segment name, side), one event per key, which is why *finalize()* refuses a second event on the same segment and side. | `(Deletion_t, "D1_gene_seq", Five_prime)` → `d1_5_del` | **Subject to change — ITER §9.** Queries by capability replace it. | keyed by (`Event_type`, `Gene_class`, `Seq_side`) until CAP B4 | — |
| **model document** | A model as one JSON document, required for any model the v1 format cannot hold. | `tst/test_data/tandem_d/trb_tandem_d_model.json` | A `structure` (`layout`, `events` keyed by nickname, `conditioning`), a `library`, `parameters`, an `error_rate` and `engine_hints`. [ModelSchema.h](../../src/igor/Model/ModelSchema.h) checks it against the schemas of `src/igor/Model/schemas/0.3/`; `ModelJson` reads and writes it. **Subject to change — SYN step 4b.** It adds `provenance` and `accumulators` sections. | the v2 text format, retired at CAP T2; "JSON v3" in REV | — |
| **v1 format** | The legacy pair of text files that holds a model. It cannot hold a tandem model. | `#Deletion;V_gene;Three_prime;5;v_3_del` | `model_parms.txt` (event lines `#Type;GeneClass;Side;priority;nickname`, `@Edges`, `@ErrorRate`) and `model_marginals.txt`. `requires_extended_format()` tells a model it cannot hold. | — | — |

---

## 3. Segments and their geometry

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **segment** | One stretch of the constructed sequence, with two *ends* and a *content*: a gene's retained template, or an *N region*. Not a *gene*. | TRB has five: `V_gene_seq`, `VD_ins_seq`, `D_gene_seq`, `DJ_ins_seq`, `J_gene_seq` | Named by its *segment name* in model documents (`layout.segments`, `targets[].segment`) and messages, identified by its *SegmentId* at run time. The legacy members still say seq type (`get_seq_type()`, `Seq_type_str_p_map`). | seq type, a `Seq_type` value; before that, a *gene class* (`VD_genes` for the VD N region) | No single term: AIRR names the genes by the prefixes `v_`, `d_`, `d2_`, `j_`, and the regions between them `np1`, `np2`, `np3`, from 5′ to 3′ |
| **segment name** | A segment's name. Any name is valid; the six legacy names are a convention. | `D1D2_ins_seq` in the tandem model | A string, `SegmentName` (`Seq_type_String` in legacy code). An event holds the name of its segment (`get_seq_type()`). | seq type, `SeqTypeString`; the `SeqType` enum, alias `Seq_type`, which knows only the six legacy names | — (not `AlleleDescription.sequence_type`, which is a *gene class*) |
| **`SegmentId`** | A segment's dense integer id, used to index every *scenario map*. | — | A `uint16_t`, resolved once by *finalize()*. `kNoSegment` stands for none: past an end of the *ordering*, or not resolved yet. | `SeqTypeId`, `kNoSeqType`; the `SeqType` enum value | — |
| **segment registry** | What converts segment names to ids and back, and holds the *ordering*. | — | [`SegmentRegistry`](../../src/igor/Core/SegmentRegistry.h). Frozen before any map is sized from it; registering a new name then throws. **Subject to change — SYN step 5, MIG L5.** The first moves the ordering into `SegmentLayout` *(planned)*. The second gives the registry to `RecombinationModel`. | `SeqTypeRegistry` | — |
| **layout** | The *segments* a model has and how they follow one another, 5′ to 3′. *Neighbours* and *end segments* are read off it, and so are the segments a *span* holds. The *layout DAG* adds alternatives. Not the *conditioning graph*, nor a storage layout (§13). | TRB: V, VD, D, DJ, J. VJ: V, VJ, J. Tandem: V, VD1, D1, D1D2, D2, D2J, J | A *model document*'s `structure.layout` ([layout.schema.json](../../src/igor/Model/schemas/0.3/layout.schema.json)), tagged by kind. The only kind is `Linear`, one chain given by its `segments`. A *v1 format* file states none, and *finalize()* infers V-VD-D-DJ-J or V-VJ-J. The code holds it as the *segment registry*'s *ordering*, and as `Topology::seqTypeOrder()` (unwired). **Subject to change — SYN step 5.** `SegmentLayout` *(planned)* will hold it apart from the registry. | segment layout (MERGE §0.1); v2's `@Seq_type_order`; "topology" in older plans | No field: a Rearrangement implies one of V-`np1`-J, V-`np1`-D-`np2`-J, V-`np1`-D-`np2`-D2-`np3`-J. CNT CT7 maps each segment to its role, read off the *ordering* |
| **ordering** | The *layout*'s segments as a list, 5′ to 3′; on a *layout DAG*, it will be a topological order of all its segments. *Neighbours* and the segments a *span* holds are read off it. Not the *model queue*. | VDJ: V, VD, D, DJ, J. VJ: V, VJ, J. Tandem: V, VD1, D1, D1D2, D2, D2J, J | A list of `SegmentId`s (`SegmentRegistry::ordering()`). For a `Linear` layout it is the layout itself. The *safety matrix* numbers its rows by position in it, and the *traversal plan* and generation's assembly walk it. **Subject to change — SYN step 5.** It moves into `SegmentLayout` *(planned)*. | — | —, see *layout* |
| **legacy id pinning** | Fixed ids for the six legacy segment names in every model, so that code still reading a segment by its legacy enum value finds it. Code that visits the segments must therefore walk the ordering, not the id range. | a tandem model has an id for `D_gene_seq`, which no event writes | `register_legacy_seq_types()` gives the six legacy names ids 0 to 5; the legacy enum is `SeqType`. **Subject to change — CAP T7b.** The pinning goes. | — | — |
| **neighbour** | The segment just 5′ or 3′ of another in the *layout*. | the neighbours of `VD_ins_seq` are V and D | *finalize()* tells each event its own (`set_adjacent_segments()`); events never look them up. **Subject to change — CAP B10.** Segments become optional, and the neighbours depend on the scenario. | adjacency; code written per event type | — |
| **end** | One side of one segment. An event *places* an end when it sets its *offset*. | TRBV's 3′ end | `SegmentBoundary{id, side}` ([SegmentSpan.h](../../src/igor/Core/SegmentSpan.h)). | boundary | — |
| **side** | Which end of a segment. A *deletion*'s side is the end it moves. A *Markov chain*'s side is its *anchor*'s end facing the chain's *N region*, which sets the direction the chain runs in. Gene choices and insertions have none. | `j_5_del`: 5′. `vd_dinucl`: 3′, its anchor V lying 5′ of it. `dj_dinucl`: 5′ | `SeqSide {Five_prime, Three_prime}`. Gene choices and insertions have `Undefined_side`. | `Seq_side`; CAP's planned `start_side`, which had the opposite convention | The `5` and `3` of `p5d_length`, `p3v_length` |
| **offset** | The read position of the nucleotide at a segment's end. A segment occupies the closed interval from its 5′ to its 3′ offset. | TRBV placed with its 3′ end at read position 40 | `SeqOffset`, an `int`, kept per scenario in `ScenarioContext::segment_ends`, a `SeqOffsetsMap`. `cut_position()` turns an end into a half-open bound (the offset at 5′, the offset + 1 at 3′), so that lengths subtract and adjacent spans add. | `Seq_Offset`; `ScenarioContext::seq_offsets`. "Offset" also names a marginal *stride* (§13) | `*_sequence_start`, `*_sequence_end`: counted from 1 (the offset + 1), and without the P nucleotides, which IGoR's offsets include. An *end segment*'s offset may lie outside the read; an AIRR position may not |
| **span** | The stretch of the constructed sequence between two *ends*, the 5′ one first. It holds the segments strictly between its ends in the *ordering*; these and its two ends decide which events change its length. While events that change it are still to run, a span has no length, only a *profile*; once the scenario has placed both its ends, it has an *observed distance*. | `gap(V, J)` holds N1, D and N2; in the tandem model, `gap(D1, J)` holds D1D2, D2 and D2J | [`SegmentSpan`](../../src/igor/Core/SegmentSpan.h), in read positions: its length is the difference of its ends' half-open bounds (`cut_position()`). `lies_strictly_inside()` tells the segments it holds, and `affects_length_of` the events that change its length. One shape is built, `SegmentSpan::gap(l, r)`: from l's 3′ end to r's 5′ end, neither segment included. **Subject to change — CAP B10, SYN step 5.** On a *layout DAG*, a span will hold the segments of every *layout path* between its ends, and exist in a scenario whose path holds the segments of both its ends. | a junction named by its insertion's seq type, with `VJ_ins_seq` standing for V→J in a VDJ model; "an ordered range of the registry ordering" (ITER §6.10) | — |
| **gap** | A span between two segments whose content is built, as the *traversal plan* sees it when an event runs. | before the D gene choice runs, V→J is one gap; after it, V→D and D→J | Each gap that needs a bound has a *bound key*. | junction | Once every event has run, the gap between two consecutive genes is an `np` region less its P nucleotides: `np1`, `np2`, `np3` |
| **junction** (code) | The code's word for three things: a *gap* that a bound is held for; the *N region* a Markov chain fills; the v1 junction classes. See §13. | the V→D gap; the N1 nucleotides | `JunctionBound`, `DinuclTraversalSpec` and `VD_genes`, in that order. | — | Not AIRR's `junction`, which is the CDR3 with its anchors |
| **content** | A segment's nucleotides. A *gene choice* (its template) or a *Markov chain* creates it; a *deletion* modifies it. | `GGGACAGGG`, TRBD1 trimmed by 3 | Held in the *constructed sequences* map. | — | A gene's content is its piece of `germline_alignment` |
| **placed**, **built**, **unfilled** | Said of a segment when an event runs, in one order of the events, and so at every node of that event. An *end* is placed once the event that sets its offset has run; for an event, a segment is placed when its facing end is. A segment is built once the event that creates its content has run; a *deletion* may still trim it. Unfilled: both ends placed, no content at all yet. Only an *N region* is ever unfilled, between its *insertion* and its *Markov chain*. | when `vd_dinucl` runs, V and D are placed and built, VD is unfilled | Read off the *traversal plan*: an end is placed once the event whose *segment end role* on it is `Creates` has run, a segment built once the event whose *content role* on it is `Creates` has run. For an event: `placed_left`, `placed_right`, `built_left`, `built_right`, `target_built`. `InitialState::placed` lists unfilled segments only, not every placed one. | chosen (`v_chosen`, `d_chosen`); slot, for unfilled; the `int_undefined` placeholders insertions wrote before ITER R1 | — |
| **anchor** | The segment a *Markov chain* takes its first nucleotide from: its neighbour on the chain's *side*. Not a *CDR3 anchor*. | TRBV anchors `vd_dinucl`; TRBJ anchors `dj_dinucl` | `DinuclTraversalSpec::anchor_id`. The code also says that an *end segment* is "anchored by the read", and that an event "anchors" a span when it creates one of the span's two ends (R12). | seed | — (not the anchor codons of `junction`) |
| **end segment** | A segment at an end of the *layout*, with no *neighbour* on one side. It may overhang the read and keeps at least one nucleotide; an internal segment may be deleted to nothing. | TRBV and TRBJ keep a nucleotide; TRBD may vanish | `Deletion::keep_one_nucleotide_` keeps the nucleotide. | anchored on an end of the read | — |
| **observed distance** | The length a scenario has fixed for a gap: the distance between its two placed ends as they stand. | 5 nucleotides between TRBV's 3′ end and TRBD's 5′ end | A *junction bound* is the gap's *profile* read at it. | — | — |
| **absent**, **present**, **not processed** | The three states of a segment in a scenario: written with no nucleotide ("actively absent"), written with some, never written. | TRBD trimmed to nothing is written absent | Held in the scenario's `DynamicSequenceMap`. `occupied()` is true for present segments only. **Subject to change — CAP B10.** It decides how a segment missing from a scenario is written. | — | An absent D gene is a null `d_call` or `d2_call` |
| **flank segment** *(planned, ERR H3, CAP B3)* | A content-only segment holding the template beyond the read, before V or after J, so that a hypermutation context near the read's edge reads it from the scenario. Unrelated to a *flank check*. | the IGHV nucleotides 5′ of the read's first position | — | `friend` access to `Gene_choice` | — |
| **layout DAG**, **existence guard**, **layout path** *(planned, SYN step 5, CAP B10)* | A *layout* with alternatives: a DAG of segments, in which a segment exists when its guard, a predicate on the realizations of *structural events*, holds. A layout path is the chain of segments one scenario follows, from a source to a sink. *Neighbours* then depend on the scenario, and the *ordering* is a topological order of all the segments, not one path. | V-D-J or V-D1-D2-J, chosen per scenario | —. Planned in SYN as `SegmentLayout::paths()`, for the layout paths. | the linear *layout* | — |

---

## 4. What an event declares

Declarations are answered from the event's own data at initialization, never during the walk
([REC_EVENT_CONTRACT.md](REC_EVENT_CONTRACT.md) §2).

| Term | Definition | Example | In the code | Legacy or former name |
|---|---|---|---|---|
| **declaration** | What an event says it does to a segment, so that generic code (*SegmentIndex*, the *traversal plan*, *pending modifier bounds*, the *fold*) never names an event type. | — | A virtual of `Rec_Event`. The first four are the A0 declarations. | capability query; switches on `Event_type` and `Gene_class`; `get_len_min()`, `get_len_max()` |
| **segment end role** | Whether the event sets an *end*'s offset, or shifts one that exists. Each end has exactly one creator. | `v_choice` creates both ends of V; `v_3_del` modifies V's 3′ end; `vd_ins` creates both ends of VD | `get_segment_end_role(id, side)` → `SegmentEndRole {None, Creates, Modifies}`. | offset role (`get_offset_role()`, `OffsetRole`); `Insertion` answered `None` until ITER R3 |
| **content role** | Whether the event allocates a segment's *content*, truncates or extends it, or fills placeholders in it. Kept apart from the *segment end role*: a *gene choice* creates both its segment's ends and its content, but for an *N region* the *insertion* creates the ends and the *Markov chain* the content. | `vd_dinucl` creates VD's content; `d_5_del` modifies D's | `get_seq_construction_role(id)` → `SeqConstructionRole {None, Creates, Modifies, Fills}`. No event answers `Fills`. | sequence construction role; before ITER R1, `Insertion` created and `Dinucl_markov` filled |
| **segment end shift** | How far the event can still shift an *end*. | `v_3_del`, from −4 to 16: `{−16, 4}` | `get_segment_end_shift_bounds(id, side)` → `SegmentEndShift{min, max}`, signed: a 3′ end moves left as nucleotides go. | offset delta (`get_offset_delta_bounds()`, `OffsetDelta`); the eight `*_min_del`, `*_max_del` members of `Gene_choice` and `Deletion` |
| **length contribution** | How many nucleotides the event adds to a segment, or removes from it, over all its realizations. | `vd_ins`, from 0 to 30: `{0, 30}` | `get_length_contribution(id)` → `LengthContribution{min, max}`. | `get_len_min()`, `get_len_max()`; `LengthRange` in ITER |
| **length delta** | The signed length one realization adds to a span it affects: + the template's length, + an insertion's value, − a deletion's value. | `TRBD1*01`: +12 | `length_delta(r)`. | the arithmetic inside four `iterate_initialize_Len_proba()` bodies |
| **`affects_length_of`** | Whether an event's realizations change a span's length. They do when the event creates or modifies an end of a segment strictly inside the span, or modifies one of the span's two ends. Creating one of the two ends does not count: the span is measured from it. | D's deletions change V→J; V's and J's gene choices change none of the gaps between them | One non-virtual rule, read off the *segment end roles* and the ordering (ITER R12). | `has_effect_on()` and four tables, one per event type (ITER R12, R16) |
| **`affects_proba_of`** | Whether an event adds a probability factor to a span without changing its length. | `vd_dinucl`, on V→D and on V→J | False by default; a Markov chain answers true for a span its segment lies strictly inside. | `has_effect_on()` |
| **participant** | An event that takes part in a span's *fold*: one that changes the span's length or adds a probability factor to it. | on V→D: `v_3_del`, `d_5_del`, `vd_ins`, `vd_dinucl` | `participates_in_span()`, which is `affects_length_of() or affects_proba_of()`. | contributor |
| **span factor** | The factor a participant that adds no length brings to a span, given the *unfilled lengths* on the fold's path. For a Markov chain, the *chain bound*. | the best N1 chain of the length `vd_ins` chose | `span_proba_factor(span, unfilled)`. `prepare_span_proba_factor()` computes what it reads, once per iteration. | `p^L` |
| **publisher** | An event that places a segment's *ends* but does not create its *content*, and so publishes the segment's *unfilled length* in the fold. One that creates both publishes 0. | `vd_ins` publishes its value; a gene choice publishes 0 | `places_own_ends()`: a *segment end role* of `Creates` on either end of its own segment. | `creates_own_offsets()` |
| **footprint** *(planned, CAP Phase A, SYN §3)* | Everything an event creates, modifies and reads: the roles and ranges above, plus the context it reads, relative to the nearest built content. It belongs to the model and serves every engine. | a 5-mer hypermutation context reads 2 nucleotides on each side | —. Planned as `context_reads() → {toward, anchor, width, codon_frame}`, for the context. | `get_context_dependency()`, `get_context_seq_types()`, never built |
| **domain**, **arity**, **candidates** *(planned, SYN, CAP Phase A)* | An event's indexed realizations; whether one realization is a scalar or a sequence over a segment's positions (a Markov chain); and the realizations an engine visits, in its order. | — | — | `is_multi_realization()`, `is_branching()`, never built |
| **mismatch role** *(planned, ERR E3)* | Whether an event creates, modifies or leaves alone a segment's *mismatch list*. It will be the test for an *own factor*. | `v_choice` creates V's list; `v_3_del` trims it | — | `SegmentIndex::has_own_factor()`, the test it replaces; `has_mismatch_writer` in CNT |

---

## 5. Reads, scenarios and their state

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **read** | One sequence inference runs on. | one non-productive TRB read covering the CDR3 | `QuerySequenceContext::sequence` and its integer form, with an index (`seq_index`) in IGoR's sequence and alignment files. | indexed sequence | `sequence`, `sequence_id`; `igor-convert` maps `seq_index` to `sequence_id`. `rev_comp` is false: IGoR aligns a read as given |
| **alignment** | One placement of a gene's template on a read, found by the aligner. | `TRBJ2-7*01` aligned at the read's 3′ end | `Alignment_data` ([AlignmentData.h](../../src/igor/Core/AlignmentData.h)), with the template's name, the read position of its first nucleotide (`offset`, which may lie outside the read), a score, mismatches, insertions and deletions. A *gene choice* enumerates its gene class's alignments of the read. | — | `v_call`, `v_score`, `v_cigar`, `v_sequence_start`, `v_germline_start` and the like, for the one alignment a tool keeps. `igor-convert` writes the best alignment of each gene class so, with a `*_alignment_length` column that AIRR does not define, and writes every alignment in the AIRR Alignment schema |
| **scenario** | One *realization* per event: a path from the root of the *walk* to a *leaf*. With the read, it fixes every segment's ends, content and mismatches. Not a class named after it. | `TRBV20-1*01`, `TRBJ2-7*01`, `TRBD1*01`, `v_3_del` = 4, `d_5_del` = 1, `d_3_del` = 3, `vd_ins` = 2, … | The `Scenario` view and `SampledScenario` are other things (§13). | recombination scenario | A `Rearrangement` row annotates a read with one scenario. IGoR keeps every scenario of a read above the threshold; CNT CT7 plans one row per recorded best scenario, with the extra columns `igor_scenario_rank` and `igor_scenario_proba` |
| **scenario probability** | The product of the marginals of the realizations chosen so far; P(s) at the leaf. | — | `ScenarioContext::scenario_proba`. | — | — |
| **error-weighted probability** | The scenario probability times the *error model*'s factor for its mismatches; P(s, read). The marginals accumulate it. | — | `scenario_error_w_proba`, computed at the leaf. | — | —; CNT CT7 writes it divided by the read's likelihood, P(s \| read), as `igor_scenario_proba` |
| **walk** | Inference's depth-first exploration of one read's scenarios: one level per event in *model queue* order, one *node* per realization, with *pruning*. | — | — | DFS; DFS-A in the plans (§11) | — |
| **node**, **hand-off**, **leaf** | A node is one realization of one event in the walk. The hand-off passes on to the next event. The leaf is reached after the last event: the error model, the threshold and the *counters* run there. | — | The hand-off is the call to the next event's `iterate()`, through `iterate_wrap_up()`. `add_to_marginals()` also runs at the leaf. | — | — |
| **`iterate()`** | One event's level of the walk: for each realization the read admits, its probability, its writes, its *flank checks*, its bound factors, the *pruning* test, then the hand-off. | — | Its contract: [REC_EVENT_CONTRACT.md](REC_EVENT_CONTRACT.md) §3.4. **Subject to change — MIG L8.** It becomes a *kernel* of `LegacyDfsEngine`. | a 19-argument signature, replaced by the five *contexts* | — |
| **next-event array** | The event after each one in the queue, none after the last. | — | `next_event_ptr_arr`, indexed by *event identifier*; null after the last event. **Subject to change — MIG L2.** An engine *continuation* replaces it. | — | — |
| **processed events** | The events before the current one in the queue. The others are *pending*. | — | `processed_events`, as `initialize_event()` sees them. | — | — |
| **layer** | One depth of the walk in a layered container: each key holds a stack of layers, so writing at layer n keeps layer n − 1, and going back up the tree costs nothing. Not a module layer (§13). | — | Held in a `LayeredArray`. | memory layer (the `memory_layer_*` members) | — |
| **claim** | The layers of a key granted to an owner. | — | Granted by `request_layer()`, in `initialize_event()`. Every write lands on a claimed layer. `claimed_layer()` is the highest claim, `current_layer()` the layer the key's data stands at (−1 if never written). | a claim made lazily by `set()`, until ITER R3b | — |
| **layer contract** | Claim before you write; write what you claimed on every path that hands off; for a key you also write, read the layer below your claim. | — | [ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) §3.4. | — | — |
| **`LayeredArray`** | The integer-keyed layered container. | — | [LayeredArray.h](../../src/igor/Core/LayeredArray.h), with flat storage `key + layer × count`. | `Enum_fast_memory_map` | — |
| **`DynamicSequenceMap`** | The layered container keyed by segment. For sequences, it tells the *absent*, *present* and *not processed* states apart. | — | [A `LayeredArray` keyed by `SegmentId`](../../src/igor/Core/DynamicSequenceMap.h), sized from a frozen registry. | `Enum_fast_memory_map<Seq_type, …>` | — |
| **scenario maps** | The layered maps of a walk: *constructed sequences*, segment ends (their *offsets*), *mismatch lists*, the *pruning mismatch floor*, the *downstream bound map*, the *safety matrix* and the *index map*. | — | Segment ends are `ScenarioContext::segment_ends`, a `Seq_offsets_map` (`SeqOffsetsMap`) ([ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) §3.3). | maps keyed by `Seq_type`; `seq_offsets` | — |
| **constructed sequence** | Per segment, the nucleotides a scenario has built so far. | TRBV trimmed by `v_3_del` | `Seq_type_str_p_map`. | — | `germline_alignment`, assembled from the genes' kept templates and P nucleotides. Where IGoR's holds the read's N nucleotides, CNT CT7 writes `N` |
| **mismatch list** | Per segment, the read positions where the read differs from the segment's content. N regions have none. | the positions where the read differs from `TRBV20-1*01` | Taken from an alignment by a gene choice, trimmed by a deletion. | — | The positions where `sequence_alignment` and `germline_alignment` differ inside a gene: the `X` operations of its `*_cigar`. `v_identity` and the like summarize them |
| **`int_undefined`** | The nucleotide code for a position allocated but not chosen yet, unlike the code for an ambiguous base. | — | The code 15; `int_N` is 14. No segment holds one at a hand-off. | −1; segments held it at hand-offs until ITER R1 | — |
| **context** | One of the five groups of state that each event's level of the walk receives. Other uses of the word: §13. | — | A struct every `iterate()` receives by reference ([ITERATE_INTERNALS.md](ITERATE_INTERNALS.md) §1). `QuerySequenceContext`: the read and its alignments. `ModelContext`: the marginals, the strides, the events map, the queue. `ScenarioContext`: the scenario probability, constructed sequences, segment ends and mismatch lists. `ExplorationContext`: the bound map, the best scenario so far and the threshold, the index map, the next-event array, the safety matrix and the pruning floor. `AccumulationContext`: the marginals being accumulated, the counters, the error model. | the arguments of the 19-argument `iterate()` | — |
| **`Scenario`** view | A flattened, read-only view of a finished scenario that counters receive at the leaf. | — | [Scenario.h](../../src/igor/Model/Legacy/Scenario.h), indexed by the six legacy ids. **Subject to change — CNT CT2a.** It is sized from the registry. | — | — |
| **scenario record** *(planned, MIG L3)* | Per event, the index of its realization and a payload, for one scenario. Counters and error models read it, and a sampling engine returns it. | `{v_choice: TRBV20-1*01, v_3_del: 4, …}` | —. Planned as `{uid → realization index, payload}`, written by each event at its layer of `ScenarioContext`. | `current_realizations_index_vec` and event members read through pointers; "realization record" in CNT | The `*_call` fields and the lengths of one Rearrangement row (CNT CT7) |
| **segment provenance** *(planned, CNT CT2b)* | Per segment, layered like the offsets: its template, the template interval kept, its palindromic extensions and its interval on the read. | `TRBJ2-7*01` kept from template position 3, placed at read positions 98 to 145 | — | — | `*_germline_start`, `*_germline_end`, `p*_length`, `*_sequence_start`, `*_sequence_end` (CNT CT7) |

---

## 6. The engine's plan

| Term | Definition | Example | In the code | Legacy or former name |
|---|---|---|---|---|
| **engine** | In the plans, one way of computing over a model, which states what it admits. | — | Today, the set-up and read loop of `GenModel::infer_model()` together with the events' `iterate()` bodies: the legacy DFS. The plans' engines are `LegacyDfsEngine`, `ChainDPEngine` and the samplers (§11). | `GenModel` |
| **`SegmentIndex`** | For each segment of the ordering, which events create or modify each end and the content, read from their *declarations*. Independent of the event order, so every engine can share it. | D's content: created by `d_gene`, modified by `d_5_del` and `d_3_del` | [SegmentIndex.h](../../src/igor/Model/SegmentIndex.h). | — |
| **own factor** | The error-bound factor a segment holds because its content can still change against the read: a deletion modifies it, or the event that creates it also places it. Content created in a segment that another event placed (an *N region*'s, by its *Markov chain*) has none: its gap carries the bound. | V, D and J have one; VD has none | `SegmentIndex::has_own_factor(segment)`. The factor sits in the *downstream bound map*. | — |
| **traversal plan** | For one order of the events: what each finds *placed* and *built* when it runs, the *gaps* it bounds and their *bound keys*, and whether the order is admissible. | — | [`TraversalPlan`](../../src/igor/Model/TraversalPlan.h): `EventPlacement` for what each event finds, and `check_admissible()`. Built once per run from the model queue, never modified, shared by every thread. | `kGeneSegments`, `legacy_junction_of()`, partner searches in `events_map` |
| **`EventPlacement`** | One event's entry in the plan: its target, the segments it finds placed or built on either side, and the gaps it lies in or faces. | `vd_ins` finds V's 3′ end and D's 5′ end placed | `target`; `placed_left` and `placed_right`, nearest first, its overlap partners; `built_left`, `built_right` and `target_built`; `enclosing_gap()` and `gap_toward(side)`. | — |
| **`InitialState`** | Segments that exist before the first event runs. | — | A test fixture that writes a segment directly. Its `placed` list holds *unfilled* segments, its `built` list built ones. **Subject to change — CAP Phase D.** A cluster with its boundaries given will be one too. | — |
| **bound key** | Where a factor sits in the *downstream bound map*. | VDJ: V→D and D→J at `VD_ins_seq` and `DJ_ins_seq`, V→J at a key of its own; VJ: V→J at `VJ_ins_seq` | A segment's *own factor* at its id; a gap at the id of the single segment inside it, when that segment has no own factor, and otherwise at a key past the registry's ids (`gap_key(span)`, `bound_key_count()`). | the junction's `Seq_type`, with `VJ_ins_seq` for V→J |
| **admissible** | Said of an order the events can run in. In the plans, each engine has its own admissibility condition. | the tandem joining order is admissible | `TraversalPlan::check_admissible()` refuses any other with `std::invalid_argument`, naming the event and the reason. The error model and the counters have a `check_admissible(const SegmentRegistry &)` of their own. | checks the events made on themselves |
| **tier** | What a computation depends on, which fixes how often it runs: 0a, the model's structure (once per run); 0b, its parameters (once per EM iteration); 1a, the read (once per read); 1b, read and parameters; 2, the scenario (at each node). Other "tiers": §13. | 0b: the junction fold; 2: a trimmed D placed on one branch | Not applied yet. **Subject to change — MIG L6, ITER §2.10.** The first moves the set-up into tiers 0a and 0b; the second's `prepare_for_query()` would be the tier-1a hook. | ITER G10's tiers 0, 1, 2 |

---

## 7. Overlap safety

| Term | Definition | Example | In the code | Legacy or former name |
|---|---|---|---|---|
| **overlap** | Two segments placed on the read sharing a position, which no scenario may contain. | TRBV's 3′ end lying past TRBD's 5′ end | Checked whenever an event places or trims a segment, while its neighbours' ends may still move. | — |
| **pending modifier bounds** | For every segment end, the sum of the *segment end shifts* of the events not processed yet, and for every segment the sum of their *length contributions*. | before `d_5_del` runs, D's 5′ end can still move | [`PendingModifierBounds`](../../src/igor/Model/JunctionGeometry.h): `segment_end_shift(id, side)`, `length(id)`. Rebuilt by gene choices and deletions in `initialize_event()`. | the `*_min_del`, `*_max_del` members and their lookup blocks; `offset_delta()` |
| **reachable interval** | The offsets an end can still reach from where it stands; a single point once nothing pending acts on it. | TRBV's 3′ end at 40, `v_3_del` pending: 24 to 44 | `reachable(id, side, current)`. | `v_3_min_offset`, `v_3_max_offset` and the six like them |
| **overlap verdict** | The answer to whether two segments overlap, one of three: even the best case overlaps; even the worst case does not; it depends on events still to come. | — | `check_overlap()` → `Overlap {Infeasible, Safe, Undetermined}`, in that order. | the `vd_check`, `vj_check`, `dj_check` blocks of each event |
| **safety matrix** | Which pairs of segments are known not to overlap in the scenario being explored. | 7 positions and 21 pairs in the tandem model | [`SafetyMatrix`](../../src/igor/Model/SafetyMatrix.h), by ordering position, one layered 32-bit word per row; hence at most 32 segments. A `Safe` verdict on (A, B) also marks (A, C) for every C 3′ of B (*row propagation*), which moves where a scenario is discarded, never whether. | `Event_safety {VD_safe, VJ_safe, DJ_safe}` in a `Safety_bool_map` |
| **`SafetyCell`** | A pair's address in the matrix. | (V, D) | `{row, column}`, the 5′-most segment first, resolved in `initialize_event()`. | — |
| **flank check** | One overlap check an event makes against its **partner**, its nearest placed neighbour on one side. A gene choice checks both sides, a deletion the side it trims. Unrelated to a *flank segment*. | `j_5_del` checks J against the nearest placed segment 5′ of it | [`FlankCheck`](../../src/igor/Model/FlankCheck.h), resolved from the event's placement by `resolve_flank_checks()`. Each check claims one safety layer. | a check against every placed segment on a side (until CAP T5) |

---

## 8. Bounds and pruning

| Term | Definition | Example | In the code | Legacy or former name |
|---|---|---|---|---|
| **pruning** | Abandoning a node whose *upper bound* lies below the read's best scenario so far times the threshold factor. The leaf makes the same test on the realized probability. | — | `should_prune()`: `upper_bound < seq_max_prob_scenario × proba_threshold_factor`. At the leaf, `is_below_threshold()`. | — |
| **upper bound** | At a node, the scenario probability times every factor of what the scenario can still reach. At the leaf the bound equals the realized probability. | — | `compute_upper_bound()`: `scenario_proba × Π downstream_proba_map[key]`, each key read at the layer the event's snapshot names. At the leaf every key is 1.0. | — |
| **downstream bound map** | The factors of what the scenario can still reach, one per key: an *own factor* (an error-model bound) or a *junction bound*. | — | `Downstream_scenario_proba_bound_map`, a `DynamicSequenceMap<double>` with `bound_key_count()` keys. An event that resolves a factor writes 1.0 there and multiplies its value into `scenario_proba`. | keyed by `Seq_type` |
| **junction bound** | One of an event's three bounds on a gap. Three suffice: a segment has two flanks and lies inside at most one gap. | TRBD's: V→D, D→J, and V→J around it | [`JunctionBound`](../../src/igor/Model/SpanProfile.h), in one of the bound slots `kLeftJunction`, `kRightJunction`, `kEnclosingJunction`, each resolved once from the plan into a span, a key, a layer and a *fold mode*. | six `{vd,dj,vj}_length_best_proba_map` members |
| **fold mode** | How a junction bound treats its gap: the bound holds the gap's *profile*, read at the *observed distance*; or the event splits the gap and the bound is 1.0; or the event splits it and builds the enclosing gap's *retained decomposition*. | TRBD's enclosing slot is `Retain` | `JunctionBound::Fold`: `Yes`, `No` and `Retain`, in that order, the last being `No` with the decomposition built. | — |
| **profile** | For one span and one *reader*: the best probability for each length the span can still reach. | the VD profile's value at 5 nucleotides; about 300 entries at most on the shipped models | [`SpanProfile`](../../src/igor/Model/SpanProfile.h): a dense vector by length. `best_for(distance)` returns nothing when no scenario reaches that length, and the realization is then discarded. | `std::map<int, double>` length maps |
| **reader** | The event that reads a junction table, at its position in the queue. A table depends on its reader and never contains it. When the reader changes the span's length, the table is conditioned on its realization: one profile per realization. | J's gene choice reads V→J | ITER R8 keeps a reader out of its own table; ITER R13 conditions the table on the reader's realization. | consumer, owner |
| **fold** | Building a reader's profiles from the events after it in the queue, state by state rather than path by path: the best probability for each total length. | — | `initialize_Len_proba_bound()`, through `SpanFold::run()`. Equal, bit for bit, to the path-by-path walk it replaced. | `iterate_initialize_Len_proba()`, the path-by-path walk, replaced at ITER R12; "Len_proba machinery" |
| **sweep** | The fold over every event, in reverse queue order, once per EM iteration. Tier 0b. | 0.055 s on the TRB regression model; 24.5 s on the tandem test model (CAP T5) | On one thread (`#pragma omp single`); the other threads copy its result (`adopt_Len_proba_bound()`, ITER S4e). | a fold per thread |
| **`FoldState`**, **`FoldFrontier`** | The part of a partial fold path a later participant can read: its length so far, the parent offsets still open, the *unfilled lengths*. The frontier keeps the best probability per state. | — | [FoldFrontier.h](../../src/igor/Model/Legacy/FoldFrontier.h). | — |
| **conditioning** (in the fold) | A participant whose parent the fold has already chosen reads its marginal at that realization, and maximizes over the other parents' configurations. | the D gene and its two deletions | `SpanConditioning` (ITER R6); `realization_bound()`. | — |
| **unfilled lengths** | Per segment, along one fold path, how many nucleotides its offsets imply that nobody has chosen yet, as its *publisher* set them. A Markov chain's factor depends on them. | VD: the length `vd_ins` chose; D: 0 | [`UnfilledSegmentLengths`](../../src/igor/Model/UnfilledSegmentLengths.h). | `SpanAccumulator`; a dummy `Int_Str` passed through the constructed sequences |
| **chain bound** | A Markov chain's *span factor*: the best chain of L steps over the fifteen nucleotide codes a read position can hold. | — | `chain_bound_`, computed once per sweep (ITER R12/R16 stage 3a). | `p^L`, from the *crude upper bound* |
| **crude upper bound** | An event's largest marginal. | — | `set_crude_upper_bound_proba()`. `Dinucl_markov::span_proba_factor()` uses it as `p^L` beyond the chain lengths it prepared. | the crude scenario bound, deleted at ITER R14 |
| **retained decomposition** | For an event that places itself inside a gap by scanning: for each total length, every (realization, left distance, right distance, bound), so that the scan can stop at the first pruned one. | TRBD's placements between TRBV and TRBJ | [`SpanDecomposition`](../../src/igor/Model/SpanProfile.h), sorted by decreasing bound. | `vj_length_d_position_proba` |
| **exhaustive position scan** | A gene choice with no surviving alignment tries every template at every position between its placed neighbours. | a 3-nucleotide TRBD remnant the aligner missed | `exhaustive_position_fallback_`, `no_alignment_survived`. Enabled for D. | `no_d_align`, still the name of a regression track |
| **surviving core** | The part of a placed gene that no pending deletion can remove. Its mismatches cannot be explained away, which bounds the error factor. | the TRBV positions 5′ of its deepest possible 3′ trim | `Gene_choice::surviving_core_bound()`. | endogenous mismatches; `credited_core_length()` |
| **pruning mismatch floor** *(dormant; CAP H1–H9)* | The mismatches that no reading of an ambiguous (IUPAC or amino-acid) query can remove, which pruning must use for such a query. | a Leu codon read as `YTN` | `Pruning_mismatch_floor_map`. Allocated and passed to the walk; nothing writes or reads it yet. **Subject to change — AA.** The amino-acid Pgen writes and reads it. | — |
| **sound bound**, **relaxation** | A bound is sound when it never lies below the best leaf under its node. A relaxation is a bound that maximizes over some choices independently. | on the uniform TRB model, J's gene choice carries about 7 of the 15 decades of slack at the root | The pruning bound is sound. It remains a relaxation: parents the fold does not choose are maximized over independently. | unsound until ITER R8 and R16 |
| **bound instrument** | A measurement, per depth, of the bound over the best leaf below each node, and of why a node reached no leaf: no child was tested; every child failed the threshold; the waste is further down. | — | [`BoundTightness`](../../src/igor/Model/BoundTightness.h), compiled in by `ENABLE_BOUND_INSTRUMENTATION` only. The causes are `BarrenCause`: `Starved`, `Pruned`, `Hollow`, in that order. | — |

---

## 9. Error models, counters and outputs

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **mismatch** | A read position whose nucleotide differs from the scenario's constructed sequence. | — | Counted per segment in the *mismatch lists*, over genes only. | — | An `X` operation of a `*_cigar`: a difference between `sequence_alignment` and `germline_alignment` |
| **`Single_error_rate`** | The uniform *error model*: a factor (r/3)^k (1 − r)^(n − k) for k mismatches over n positions, with r re-estimated each iteration. | r = 0.000396 in the human TRB model | Charges the positions of every segment with a mismatch list. **Subject to change — ERR E1 to E4.** It becomes a shared `SingleErrorModel` core and an *error node*. | before CAP T7a, ids 0, 2 and 4 only | — |
| **hypermutation models** | The two error models of hypermutation: a mutation probability from a logistic model of the *N-mer context*, or one probability per N-mer. | IGHV, with 5-mer contexts | `Hypermutation_global_errorrate` and `Hypermutation_full_Nmer_errorrate`, in that order. They refuse *layouts* other than V-D-J and V-J. **Subject to change — CAP T7b, ERR H2.** The first lifts the refusal; the second makes them one context model with two parameterizations. | — | — |
| **N-mer context** | The N nucleotides centred on a position, which a hypermutation model reads. | the 5-mer around one IGHV position | — | — | — |
| **scope** | The segments a hypermutation model charges and learns from. | `VDJ_genes` in the shipped files | `apply_to`, `learn_on`; the file key is `apply_on`. **Subject to change — ERR H4.** They become ordered sets of segment names. | — | — |
| **Viterbi run** | An inference that keeps each read's best scenario only. | — | `infer.viterbi`; `add_to_marginals()` then sets rather than adds. **Subject to change — ERR E2.** The option moves from the error model into the engine. | — | — |
| **Probe event**, **error node** *(planned, ERR E4, SYN step 8)* | The error model as a node of the model: it never branches, reads every segment in its scope, and keeps its parameters in a tensor with handlers. | — | — | `Error_rate`; CAP Phase E's `ErrorRate_Rec_Event` | — |
| **counter** | A statistic computed during inference, at each leaf above the threshold and once per read. | — | [`Counter`](../../src/igor/Model/Legacy/Counter.h), called with the `Scenario` view: `Best_scenarios_counter`, `Coverage_err_counter`, `Errors_counter`, `Pgen_counter`. All four refuse *layouts* other than V-D-J and V-J. **Subject to change — CNT CT2a, CT5; MIG L8.** The first lifts the refusal; the other two make the counter an **observer**. | — | — |
| **best scenarios** | Each read's most probable scenarios, with their realizations and mismatches. | — | `Best_scenarios_counter`; ties compared at `output.scenario_tie_digits`. | — | Rearrangement rows, once CNT CT7 lands; IGoR's own file until then |
| **coverage** | Per template position, the posterior weight of the scenarios that cover it, and of those with a mismatch there. | coverage along `TRBV20-1*01` | `Coverage_err_counter`. | — | — |
| **output event names** | How the outputs name the events: by *nickname* or by *generated name*. | — | `output.event_names` (`Event_label`): `nickname`, the default, or `legacy`, the *generated names*. | — | — |
| **per-read summary** *(planned, MIG L1)* | What an engine returns for each read: the likelihood, the mean number of errors, the number of scenarios, the best one. | — | —. Planned as `ReadSummary`. | the error model's per-read getters; `Errors_counter` | — |
| **row sink** *(planned, CNT CT6)* | An output that declares its row schema and is written chunk by chunk, in chunk order, as today's text files or as Parquet. | — | — | — | — |
| **AIRR rearrangement output** *(planned, CNT CT7)* | The best scenarios as AIRR rows, the AIRR V, D and J roles read off the ordering. | `v_call` `TRBV20-1*01`, `junction_aa` `CASSLGQAYEQYF` | — | — | The Rearrangement schema, read off the *scenario record* and the *segment provenance* (CNT §9.3) |

---

## 10. Generation

| Term | Definition | Example | In the code | Legacy or former name | AIRR |
|---|---|---|---|---|---|
| **generated sequence** | A sequence the model generates: the one assembled from a drawn scenario, with errors added if asked. | a synthetic TRB sequence | Written with its index (`seq_index;nt_sequence`), with its scenario's realizations in a second file, and with its CDR3 data in a third when `generate.cdr3` is set. | — | A Rearrangement with `rearrangement_type` `simulated`, in a Repertoire with `repertoire_type` `simulated`. IGoR writes its own files today |
| **draw** | Picking realization indices from an event's marginal row. A Markov chain draws its whole chain. | — | `draw_realization()`, which writes nothing; `pick_realization()`, a walk of the cumulative sum. | `draw_random_realization()`, which also built the sequence | — |
| **construct** | Building what drawn realization indices stand for into the generator's state, with no random draw, marginals or index map: the event's *apply*. It refuses what inference never counts. | `d_3_del` = 3 removes three nucleotides from TRBD's 3′ end | `construct_realization(indices, state)`, into the `GenerationState`. Its refusals follow GEN D2, D3. | `Event::apply` in SYN; `apply_realization` in MERGE | — |
| **propagate** | Moving each child's marginal row by the drawn realization times the child's *stride*. | — | `propagate_realization()`. | — | — |
| **`GenerationState`** | The generator's per-sequence state, per segment. The sequence is assembled in *ordering* order. | — | [GenerationState.h](../../src/igor/Model/Legacy/GenerationState.h), keyed by `SegmentId`. | `unordered_map<Seq_type, string>` | — |
| **`SequenceGenerator`** | The legacy generator: it draws each event in queue order. | — | `draw_random_realization()` on each event. **Subject to change — MIG L7.** A sampling driver and `SamplingEngine` replace it. | — | — |
| **`FastGenerator`** | The fast generator: precomputed CDFs, alias sampling, threads. | — | Run by `--fast`. Still finds segments through `gene_to_seqtype_migr`. **Subject to change — CAP T4b.** Whether it is kept is open. | — | — |
| **anchor nucleotide** | A chain's first state in generation: its *anchor*'s facing nucleotide. | TRBV's last nucleotide starts the N1 chain | It must be constructed before the chain is drawn. | — | — |
| **sampling driver** *(planned, MIG L7)* | The generation counterpart of the *EM driver*: random streams per seed and chunk, error application, writers. | — | — | `SequenceGenerator` | — |

---

## 11. The Model module and the target architecture

The Model module's tensor architecture is merged but unwired (MIG §1.1);
[ENGINE_MIGRATION_PLAN.md](../ENGINE_MIGRATION_PLAN.md) wires it.

| Term | Definition | Example | In the code | Legacy or former name |
|---|---|---|---|---|
| **`RecombinationModel`** *(unwired)* | A model as a *topology* and one *tensor* of weights per event. | — | Tensors indexed by *uid*, with the parents' axes first and the event's own axes last. Handlers borrow the tensors. **Subject to change — MIG L5.** The engine reads it, and it also owns the registry and the *layout*. | `Model_Parms` with `Model_marginals` |
| **`Topology`** *(unwired)* | The events, with their edges and the *ordering*. | — | The events as `Rec_Event`, their edges sorted by uid, the ordering by segment name (`seqTypeOrder()`); built by `topology_from_json()`. | `Model_Parms`; "topology" also meant the *conditioning graph* or the *layout* in older plans |
| **tensor** *(unwired)* | An N-dimensional array of runtime rank, one per event. | `j_choice`'s weights, of shape [V genes, J genes] | `Tensor<T>`, with mdspan views; one per event in `RecombinationModel`. | `Matrix<T>`; one event's rows of `Model_marginals` |
| **handler** *(unwired)* | What holds an event's numbers for one task, inference or sampling. Each is Categorical or Markov. | — | An inference handler is an `InferenceHandler`: a mutable tensor reference and an accumulator, `maximizeLikelihood()`. A sampling handler is a `SamplingHandler`: a const reference and CDFs, `sample()`. | an event's rows of the marginals, with `add_to_marginals()` and `ind_normalize()`; SYN calls the family a parameterization: Categorical, Markov, Bernoulli, Logistic, LinearModel |
| **`InferenceEngine`**, **`SamplingEngine`** *(unwired)* | tk's containers of handlers, for inference and for sampling. | — | `InferenceEngine::run()` resets the accumulators, calls the caller's E-step, normalizes; `SamplingEngine` samples one index per node in topological order. **Subject to change — MIG Q1.** `InferenceEngine` becomes the name of the engine interface. | — |
| **`Navigator`** *(unwired)* | A view over part of a node vector: parents, children or topological order. | — | Non-owning and random-access. | — |
| **`LegacyBridge`** *(unwired)* | The one copy point between the two architectures, both ways. | — | `Model_Parms` to a `Topology`, `Model_marginals` to tensors, and back. **Subject to change — MIG L5.** After it, the bridge becomes the v1 *importer*. | `MarginalsBridge`, proposed then dropped |
| **event factory**, **`Registrar`** | Constructors that register themselves under a type name. | — | [`event_factory::create()`](../../src/igor/Model/EventFactory.h) builds an event from its JSON node, dispatching on its `type`; the handler factories (unwired) build handlers by `EventType`. | — |
| **accumulator**, **`AccumulatorSet`** *(planned, MIG L1)* | What EM accumulates, one per node. | the posterior counts of (`v_choice`, `j_choice`) for one read | —. Planned with `empty_copy()` and an ordered merge over chunk slots. | copies of `Model_marginals` per chunk in `GenModel` |
| **EM driver** *(planned, MIG L1)* | The generic loop around an engine: iterations, data and read policies, chunks, threads, per-read normalization, ordered merge, M-step, counter lifecycle. | — | — | `GenModel::infer_model()` |
| **`LegacyDfsEngine`** *(planned, MIG L1, L8)* | Today's *walk* as one engine. | — | —. Planned with `check_admissible`, `prepare_iteration`, `make_workspace` and `process_read`. Its *kernels* are the `iterate()` bodies, moved to `Inference/Dfs/` at L8. | `GenModel` with `Rec_Event::iterate()` |
| **kernel** *(planned, MIG L8)* | The function that runs one event type's step for one engine. Not a Markov chain's transition matrix (§13). | — | — | an event's `iterate()` body |
| **continuation** *(planned, MIG L2)* | The engine call by which an event hands its probability and bound facts to the engine, which runs the next event or the leaf. | — | — | `iterate_wrap_up()` and the next-event array, which it replaces |
| **DFS-A**, **DFS-B** *(planned, REV §10.3)* | DFS-A: today's walk in engine form, events keeping their loops and writing at their claimed layers. DFS-B: an orchestrator asks each event for its candidates and has it evaluate them, writing detached outcomes. | — | —. DFS-B is planned with `candidates()` and `evaluate()`, as a prototype measured after SYN step 7. | — |
| **`ChainDPEngine`** *(planned, SYN step 9, CAP Phase D)* | A dynamic-programming engine: an outer sum over gene-choice tuples, and a forward-backward pass along the *layout* over *boundary states* inside it. Built on it: a walk pruned by exact DP bounds, exact posterior sampling (FFBS), k-best scenarios. | — | —. The expected engine for real tandem data (MIG Q14). The bounded walk is planned as `DpBoundedDfsEngine`. | CAP Phase D's "DP subgraph decomposition" |
| **boundary state** *(planned)* | What a segment exposes to its neighbour in the DP: an end's offset and its end nucleotides (k of them for an N-mer error model). | TRBV's trimmed 3′ offset and its last nucleotide | — | — |
| **semiring** *(planned)* | The algebra one DP pass runs in: (+, ×) for likelihood and posteriors, (max, ×) for best scenarios and bounds. The *fold* is the (max, ×) case with length as its only boundary variable. | — | — | — |
| **parameter tying** *(planned, CAP B9, SYN §4)* | Two nodes sharing one parameter entry; the M-step sums their accumulators. | D1 and D2 sharing one D gene distribution in a tandem IGH model | — | — |
| **covariate** *(planned, SYN)* | An observed node outside the recombination that conditions events. | the donor's age (AIRR `Subject.age`) conditioning insertion lengths | — | — |
| **structural event** *(planned, SYN step 5)* | An event whose realization chooses a path of the *layout DAG*. | the number of D genes in the rearrangement | — | — |
| **importer** *(planned, SYN step 4, MIG L5)* | What reads v1 text models into the new representations, keeping the legacy indices so that the tensor layout equals the legacy one. | — | —. Planned to read into a *model document* or a `RecombinationModel`. | `LegacyBridge`; the reader of `Model_Parms` |

---

## 12. Plans and process

### 12.1 Plan short names

| Short | Document | Holds |
|---|---|---|
| CAP | [REC_EVENT_CAPABILITY_REFACTORING_PLAN.md](../REC_EVENT_CAPABILITY_REFACTORING_PLAN.md) | tandem D: tasks B0–B11, steps T0–T7, Phases A to E, milestones 1 and 2 |
| ITER | [ITERATE_GENERIC_REWRITE_PLAN.md](../ITERATE_GENERIC_REWRITE_PLAN.md) | the generic `iterate()` bodies: patterns G1–G10, services S1–S5, phase R (R0–R16), decisions O1–O12 |
| PBM | [PROBA_BOUND_MACHINERY.md](../PROBA_BOUND_MACHINERY.md) | the history and measurements of the junction bound |
| MIG | [ENGINE_MIGRATION_PLAN.md](../ENGINE_MIGRATION_PLAN.md) | from `GenModel` to drivers and engines: steps L0–L8, questions Q1–Q16, lanes |
| ERR | [ERROR_MODEL_REWRITE_PLAN.md](../ERROR_MODEL_REWRITE_PLAN.md) | the error models: steps E0–E4 and H0–H5, defects X, repairs ER and HR |
| CNT | [COUNTER_REWRITE_PLAN.md](../COUNTER_REWRITE_PLAN.md) | the counters: steps CT0–CT7, defects CX, repairs CR, questions QC |
| GEN | [GENERATION_REWRITE_PLAN.md](../GENERATION_REWRITE_PLAN.md) | generation: steps G0–G4, defects D1–D6 |
| SYN | [ARCHITECTURE_SYNTHESIS.md](../ARCHITECTURE_SYNTHESIS.md) | the target architecture and its steps |
| REV | [SEGMENT_DECOMPOSITION_REVIEW.md](../SEGMENT_DECOMPOSITION_REVIEW.md) | the design review behind SYN: engines, tiers, DFS-A and DFS-B, the model document |
| LAYER | [doc/LAYER_REFACTORING_PROPOSAL.md](../../doc/LAYER_REFACTORING_PROPOSAL.md) | the module layers and their `Legacy/` directories |
| MERGE | `TK_REFACTORING_MERGE_ANALYSIS.md` | the merge of tk's Model module; not in the repository on this branch |
| AA | [AA_PGEN_IMPLEMENTATION_PLAN.md](../AA_PGEN_IMPLEMENTATION_PLAN.md) | AA Pgen; the hazards of merging it are CAP H1–H9 |
| — | [IMPLEMENTATION_GRAPH.md](../IMPLEMENTATION_GRAPH.md) | the order of every remaining step |
| — | [ITERATE_METHOD_ANALYSIS.md](../ITERATE_METHOD_ANALYSIS.md), [CONTEXT_BASED_ITERATE.md](../CONTEXT_BASED_ITERATE.md) | `iterate()` before the rewrite, and the introduction of the contexts |

### 12.2 Process terms

| Term | Meaning |
|---|---|
| **step** | A unit of work in a plan, named by its plan and label: MIG L3, CAP T7b, ITER R12. §14 lists the labels that collide. |
| **milestone 1**, **milestone 2** | CAP's: tandem D on the legacy DFS with every segment present, delivered at T5; tandem D on real repertoires, which needs an optional D2 (B10). |
| **lane** | MIG §4: steps that share files and run in sequence. T (tests), D (driver and engine), E (error models), C (counters), M (Model module), G (generation). |
| **checkpoint** | A merge of the migration's integration branch into `develop`, after MIG L2, L5 and L8. |
| **regression tracks** | `pixi run test_regression`, on human TRB: `cli`, `align`, `inference`, `no_d_align`, `generate`. |
| **golden**, **golden move** | The stored expected outputs of the regression tracks. A golden move changes them on purpose, in a commit of its own that names the outputs it moves and why. |
| **bitwise** | Said of a step that changes no output of the regression tracks. Marginal files print 6 significant digits, so until MIG L0's A/B oracle compares at full precision, bitwise means equal to 6 digits. |
| **characterization suite** | Tests written against the unmodified code before a refactor, checked by mutating the code: ITER T0, CNT CT0, ERR E0 and H1. |
| **`[!shouldfail]` pin** | A confirmed defect asserted at its intended value under Catch2's `[!shouldfail]`, until the repair that fixes it. |
| **repair** | A behaviour change held back from a refactor, then made once, naming the outputs it moves: ITER's phase R, ERR's ER and HR, CNT's CR. |
| **migration code** | Code that only bridges a migration: a `_migr` file in the `migration` namespace, deleted whole at the step its header names. `vdj_ordering_migr` goes after CAP T7b and CNT CT2a, `gene_to_seqtype_migr` at CAP T6b. |
| **Legacy** | In a path or a namespace (`Model/Legacy/`, `igor::model::legacy`): code from before the layer split that has not reached its final place (LAYER §4). It runs every inference and generation today. |
| **wired**, **unwired** | Whether inference or generation calls a piece of code. |
| **subject to change** | The note that marks a part of a reference document or a README that a planned step changes: "Subject to change — *plan* *step*". |

---

## 13. Words with several meanings

| Word | Meanings | In new text |
|---|---|---|
| segment | a *segment* of the *layout*; a germline *gene* ("gene segment" in biology); a group of events, in SEGMENT_DECOMPOSITION_PROPOSAL; AIRR's `segment` (Alignment schema) and `AlleleDescription.sequence_type`, V, D, J or C: a *gene class* | segment for the layout's pieces, gene for the germline element, gene class for AIRR's `segment` |
| junction | the biological *junction*; a *gap* a bound is held for (`JunctionBound`); the *N region* a Markov chain fills; the v1 junction classes; AIRR's `junction`, the CDR3 with its two anchors | gap for the stretch between built segments, N region for the inserted segment, junction for the biology |
| CDR3 | the region between the two anchor codons, AIRR's `cdr3`; the `nt_CDR3` column of `generate.cdr3`, which includes the anchor codons, as AIRR's `junction` does | CDR3 without the anchors, junction with them |
| anchor | a Markov chain's *anchor*; the *CDR3 anchors*; an *end segment*, "anchored by the read"; an event creating a span's end (R12); the `anchor` of `context_reads()` | qualify it: chain anchor, CDR3 anchor |
| chain | a *Markov chain*; a receptor chain, such as TCR β or the immunoglobulin heavy chain; the *locus* option `model.chain` | Markov chain in full; locus for the option |
| germline | the genes of the germline DNA; AIRR's reference sequences (`GermlineSet`, `germline_alignment`), IGoR's *templates* | template for a sequence of the library |
| layer | a *layer* of a `LayeredArray`; a module layer (Core, Model, Inference…; LAYER); `ScenarioContext::Layer` in REV | layer for the walk's depth, module for the code's division |
| context | the five *contexts*; a handler's conditioning context, its parents' realizations; the sequence an event reads (`context_reads()`); an *N-mer context* | the class names in full; parent configuration; N-mer context |
| tier | the computation *tiers* 0a to 2; ITER's enforcement tiers 1 to 3 (capability, leaf, hand-off), which SYN renames check levels C1 to C3; the test tiers of ITERATE_METHOD_ANALYSIS | tier for computation tiers only |
| scenario | a *scenario*; the counters' `Scenario` view; tk's `SampledScenario` (generation, unwired); `ScenarioContext`; the planned *scenario record* | scenario for the path; the class name otherwise |
| scenario record | MIG's name; CNT says realization record | scenario record |
| slot | an *unfilled* segment (former name); a position of the *layout* (SYN's "a slot in the layout"; "two D slots" in tests); one of an event's three *junction bounds* (`JunctionSlot`); a key of the *downstream bound map*; a *profile* entry (`SpanProfile::slots_`); a storage cell ("a boundary is a slot"; the safety matrix's "one slot per pair"); the chunk slots of the EM reduction (`GenModel`; MIG L4a, CNT CT4); the events' `*_slot()` accessors, which `Coverage_err_counter` reads | slot only where an identifier says it: `JunctionSlot`, `*_slot()`, chunk slots. Unfilled or N region for the segment, segment for a position of the layout, key for the map |
| offset | a segment end's read position; a marginal *stride* (`get_offsets_map()`, `memory_and_offsets`); an alignment's offset on the read | offset for read positions, stride for the marginals |
| layout | the model's *layout* (JSON `structure.layout`); the *layout DAG*; a tensor's or array's storage layout (mdspan's `layout_stride`, the marginals' memory); the numbering of the *bound keys* (`TraversalPlan`'s "layout rule"); a file's sections (REV's "file layout") | layout alone for the segments; tensor or storage layout for memory; key numbering for the bound keys; sections for a file |
| order | event order, the *model queue*; segment order, the *layout*'s *ordering*; the *joining order*, one model queue; `Topology::topologicalOrder()`, the events in an order of the *conditioning graph* | model queue for the events, ordering for the segments |
| engine | an *engine* in the plans' sense; tk's `InferenceEngine`, a handler container, whose alias `LegacyEngine` is not `LegacyDfsEngine`; `GenModel` today; `FastGenerator`, "an approximate engine" | engine for a computation scheme |
| kernel | an engine's per-event function (MIG L8); a Markov chain's transition matrix; REV's "E-step kernel" | transition matrix for the chain |
| bound | a probability *upper bound*; an offset range (`SegmentEndShift`, `PendingModifierBounds`); the error model's upper bound | probability bound; range |
| fold | the junction *fold*; `JunctionBound::Fold`, a mode in which `No` means not folded; ITER's max-product composition ⊗ | — |
| flank | a *flank check*, for overlap; a planned *flank segment*, the template outside the read | — |
| footprint | an event's *footprint* (planned); the remnant of a replaced VH in a CDR3, in biology | — |
| Creates | `SegmentEndRole::Creates`, said of an *end*, and `SeqConstructionRole::Creates`, said of a *content*: two statements. A *gene choice* makes both for its segment; for an *N region* the *insertion* creates the ends and the *Markov chain* the content | — |
| error rate | the *error model* (`Error_rate`, `@ErrorRate`); Single's parameter r | error model for the component, rate for r |
| Pgen | the generation probability; `Pgen_counter`; *AA Pgen* | — |
| importer | `LegacyBridge` after MIG L5; `Model_Parms` as the v1 reader; lane M's importer step | — |
| Probe | SYN's error-model event; SEGMENT_DECOMPOSITION_PROPOSAL's role, which also covered counters; AA Pgen's patched leaf in MIG | — |
| fragment | REV's partial export of a model document; the code's parameters-only document (`ModelJson`) | — |
| domain | an event's realizations (SYN); MERGE's attribute domains, areas of concern | — |
| Markov chain | spelled `Dinucl_markov` (class), `Dinuclmarkov_t` (enum), `DinucMarkov` (file and JSON token), "Dinucl" in the plans | Markov chain in prose, the identifier in code |

---

## 14. Step labels that collide

Cite a step with its plan's short name: the same label means different things in different plans.

| Label | Meanings |
|---|---|
| A | CAP's A0 declarations; CAP Phase A; CNT's access patterns A1–A9 |
| B | CAP's tasks B0–B11 |
| C | CAP Phase C; MIG lane C; check levels C1–C3 (SYN); ITER's standing constraints C1 and C2; CNT's requirement checks C1 and C2 |
| D | the tandem segments D1 and D2; CAP's decisions D1–D5; GEN's defects D1–D6; CAP Phase D; MIG lane D |
| E | ERR's steps E0–E4; CAP Phase E; the E-step |
| G | ITER's patterns G1–G10; GEN's steps G0–G4; MIG lane G |
| H | ERR's hypermutation steps H0–H5; CAP's AA Pgen merge hazards H1–H9 |
| O | ITER's decisions O1–O12, and observations O13–O15 |
| Q | each plan's questions: MIG Q1–Q16, ERR Q1–Q7, CNT QC1–QC10 |
| R | ITER's phase R, rows R0–R16; ERR's repairs ER and HR; CNT's repairs CR |
| S | ITER's services S1–S5 |
| T | CAP's steps T0–T7; MIG lane T; ITER's T0 harness |

---

## 15. AIRR concepts IGoR has no term for yet.

Searching this page for an AIRR name finds its IGoR equivalent. These have none. If one is required further in development used the AIRR term.

| AIRR | What it is | In IGoR |
|---|---|---|
| `c_call` and the other `c_` fields; `AlleleDescription.sequence_type` `C` | the constant region | not modelled: a model ends with J |
| `cdr1`, `cdr2`, `fwr1` to `fwr4`, and their positions | the framework regions and the first two CDRs | not modelled: only the CDR3 anchors are read |
| `Clone`, `Node`, `clone_id`, `Clone.inferred_ancestor` | reads descended from one naive rearrangement, and their lineage | not modelled: inference treats every read as a rearrangement of its own |
| `cell_id`; `SampleProcessing.single_cell`, `physical_linkage` | single-cell data and the pairing of chains | not modelled: a model describes one locus |
| `duplicate_count`, `consensus_count`, `umi_count` | the copies of one sequence | not read: every read of the input counts once |
| `v_rs_start`, `d_rs_5_prime_start`, `d_rs_3_prime_start`, `j_rs_start`, and their ends (Germline) | the recombination signal sequences | not modelled |
| `GenotypeSet`, `Genotype`, `deleted_genes`, `undocumented_alleles` | a donor's germline genotype | not modelled: the *library* stands for it |
| `Study`, `Subject`, `Sample`, `SampleProcessing`, `DataProcessing` | the metadata of a repertoire | not read; the planned *covariates* would read subject fields such as `Subject.age` |
