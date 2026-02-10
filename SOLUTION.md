# Proposal: Resolving Event Name and Stride Mapping for Legacy Import

## Problem Analysis
The user identified an issue where `Model_marginals::get_offsets_map()` keys do not match event names directly, or at least the relationship requires investigation.
- `get_offsets_map()` maps **Parent Event Names** -> **List of (Child, Stride)**.
- We define `MarginalHandler`s by the **Child Event Name**.
- We need to know the strides/offsets of **Parents** to correctly map the linear legacy array into our multi-dimensional `Tensor`.

## Value Discrepancy or Key Confusion?
It is likely that the "mismatch" refers to the fact that we cannot look up an event's structure by querying its own name in `get_offsets_map()` (because that map keys by Parent). Instead, we should use `get_inverse_offset_map()`, which keys by the Event itself.

## Proposed Algorithm

1.  **Use `get_inverse_offset_map()`**:
    - This function returns `UnorderedMap<EventName, List<Pair<ParentEventPtr, Stride>>>`.
    - Key: The name of the event we are building a handler for suitable for `InferenceEngine::register_handler`.
    - Value: The list of conditioning parents and their stride in the flat array.

2.  **Determine Dimension Order**:
    - For a given event `E` (e.g., "d_gene"):
    - Retrieve its parents and strides: `[(v_choice, 1), (other_parent, 50), ...]`.
    - **Sort** these parents by their stride value. This reveals the memory layout (dimension ordering) of the legacy array.
        - *Legacy layout is effectively determined by the order events were pushed to the stack.*

3.  **Construct Shape and Strides**:
    - `Event` itself corresponds to **Dimension 0**, with `Stride = 1`.
    - The sorted parents correspond to **Dimensions 1..N**.
    - Validate that `Stride[i] == Stride[i-1] * Size[Dim[i-1]]`.
        - If this holds, standard row-major (or column-major depending on sort direction) packing is used.
        - If not, we must perform a **permuted copy** or use **strided views**.

4.  **Handling "Mismatch"**:
    - If there is a literal string mismatch (e.g., `GeneChoice_V_gene` vs `V_gene`), we will use the `Rec_Event_name` directly from the `Rec_Event` pointer in the `Model_Parms` graph to ensure consistency.

## Implementation Plan for `LegacyBridge`

1.  **Iterate `Model_Parms::get_edges()`**: To discover all events.
2.  **Call `get_inverse_offset_map(parms)`**: To get the stride data.
3.  **Synthesize `EventDesriptor`**:
    - Name: `event->get_name()`
    - Type: Map `event->get_type()` to our enum.
    - Shape: Derived from sorted strides + event size.
4.  **Data Copy**:
    - Flatten the handler's tensor.
    - Copy from `Model_marginals` giant array using the calculated start index (`get_index_map`) and the strides.

This approach bypasses the confusion of `get_offsets_map` (parents view) by using `get_inverse_offset_map` (childs view), which aligns with our handler-centric architecture.
