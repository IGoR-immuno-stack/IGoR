# Tandem D, milestone 1

The model and reads of capability plan T5 (`docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md`). Both are
generated; do not edit them by hand.

- `trb_tandem_d_model.json`: a V-D1-D2-J model built from human TRB by
  `scripts/tandem_d/make_tandem_d_model.py` (T5a). Its header says how.
- `reads/`: 200 reads generated from that model and aligned against TRB's genomic templates, by
  `scripts/tandem_d/make_tandem_d_reads.py` (T5b), seed 42, on one thread:
  - `sequences.csv`: the reads, `seq_index;nt_sequence`;
  - `realizations.csv`: the scenario that generated each, one column per event in the model's
    queue order, then its errors;
  - `V_alignments.csv`, `D_alignments.csv`, `J_alignments.csv`: for each read, the alignments an
    inference reads, those scoring within 55 of its best for V, 35 for D and 10 for J.

Regenerate the model first, then the reads: the reads are drawn from the model's parameters.
