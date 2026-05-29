# Deletion Identity Matrix

This matrix documents every usage of `Gene_class` (or `event_class` / `event_class_str`) within `Deletion.h` and `Deletion.cpp`, classifying each reference to control risk and prevent over-refactoring.

| Code Reference (Line/Scope) | Variable / Expression | Current Type / Usage | Classification | Refactoring Action |
|:---|:---|:---|:---|:---|
| **`Deletion` Constructor** | `Gene_class gene` | Event descriptor input representing genomic origin (V, D, J) | **Genomic & Segment** | Retain genomic origin in an internal `Gene_class` field; map/initialize a target `Seq_type` field. |
| **`Deletion::iterate`** (Line 218) | `switch ((*this).event_class)` | Controls the tree-walk logic branching for V, D, and J segment deletion processing | **Segment Operations** | Switch on target `Seq_type` (e.g. `V_gene_seq`, `D_gene_seq`, `J_gene_seq`) instead of `event_class`. |
| **`Deletion::iterate` (V_gene)** (Lines 201-203, 222-228) | `d_chosen`, `j_chosen` | Queries whether genomic choice events were chosen to check overlaps | **Genomic Identity** | Continue querying genomic choice flags `d_chosen` and `j_chosen` representing structural dependencies. |
| **`Deletion::iterate` (V_gene)** (Lines 205, 225) | `Event_safety::VD_safe`, `Event_safety::VJ_safe` | Enforces sequence overlap constraints at genomic junction boundaries | **Genomic Identity** | Retain current `Event_safety` enum checks (they model structural dependencies, not individual segment indices). |
| **`Deletion::iterate` (V_gene)** (Lines 376, 383, 387) | `V_gene_seq`, `VD_ins_seq`, `VJ_ins_seq` | Updates offsets, mismatches, and sequences on the segment map | **Segment Operations** | Already uses `Seq_type`; no change. Ensure any remaining references use `Seq_type` parameter. |
| **`Deletion::iterate` (D_gene)** (Lines 462, 483, 484) | `v_chosen`, `d_mismatch_list`, `D_gene_seq` | Performs sequence offsets and mismatch adjustments on D segment | **Segment Operations** | Switch hardcoded segment lookups to query target `Seq_type` dynamically. |
| **`Deletion::draw_random_realization`** (Line 1229) | `switch (this->event_class)` | Directs random sampling deletion mutations on constructed sequences | **Segment Operations** | Switch on `Seq_type` and erase bases from the corresponding `Seq_type` segment dynamically. |
| **`Deletion::initialize_event`** (Line 1351) | `switch (this->event_class)` | Requests memory layers for offsets, mismatches, and safety boundaries | **Mixed** | Request memory layers for safety (`VD_safe`, etc.) based on genomic topology; request segment maps (`V_gene_seq`, etc.) based on `Seq_type`. |
| **`Deletion::iterate_initialize_Len_proba`** (Line 1644) | `get_deletion_effective_junctions(this->event_class, ...)` | Maps deletion event to affected sequence junctions (VD, DJ, VJ) | **Segment Operations** | Transition this mapping helper to determine junctions based on the adjacent segment (`Seq_type`) and `Seq_side`. |

---

## Detailed Classification Rationale

### 1. Truly Genomic Identity (Must Stay `Gene_class`)
- **Safety Overlaps (`Event_safety`)**: Checking if a V-deletion overlaps a D-deletion or J-deletion (`VD_safe`, `VJ_safe`) represents topological constraints of the V(D)J recombination process. Even with tandem Ds, the structural boundary checks (e.g. "is there a safety boundary between V and the first D?") are topologic properties of the locus.
- **Gene Presence / Choice flags (`v_chosen`, etc.)**: Checking if a genomic segment of class V, D, or J has been chosen in the scenario.

### 2. Constructed-Segment Operations (Convert to `Seq_type`)
- **Sequence Segment Lookup**: Fetching sequences (`constructed_sequences.at(V_gene_seq)`), offsets, or mismatch arrays. These must query `Seq_type` (like `D1_gene_seq` vs `D2_gene_seq`) because two distinct D segments in tandem D models will have the same `Gene_class` (`D_gene`) but separate constructed sequences.
- **Writing / Serializing Deletion counts**: Identifying which segment is being deleted from in model writes and event maps.
