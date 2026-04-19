#pragma once

#include <igor/Core/Utils.h>
#include <igor/Core/IntStr.h>

/**
 * @brief Encapsulates per-path mutable state during scenario exploration
 * 
 * ScenarioContext holds state that changes as we recursively explore
 * possible recombination scenarios. Critically, this keeps the hot-path
 * unified: probability computation and sequence building are handled
 * together (accessed in tight loops).
 * 
 * Contains:
 * - Running probability of current scenario path
 * - Partially constructed sequences (gene/junction segments)
 * - Offsets of sequence fragments (5'/3' positions)
 * - Mismatches accumulated for error rate calculation
 * 
 * Design principle: Zero-copy - store references/pointers to minimize
 * overhead in hot path (called millions of times per sequence).
 */
struct ScenarioContext {
    // Scenario probability (multiplied at each recursion level)
    double& scenario_proba;
    
    // Error-weighted probability (computed at leaf nodes)
    // NOTE: This is scenario state computed once we have full sequence + error model evaluation
    // Set by Error_rate::compare_sequences_error_prob() at leaf nodes
    long double scenario_error_w_proba;
    
    // Constructed sequences (built up during recursion)
    // Maps Seq_type → string pointer with memory layer
    Seq_type_str_p_map& constructed_sequences;
    
    // Offsets of each sequence fragment (5'/3' positions)
    Seq_offsets_map& seq_offsets;
    
    // Mismatch positions for error rate calculation
    Mismatch_vectors_map& mismatches_lists;
    
    /**
     * @brief Constructor - binds references to existing state
     * 
     * All parameters are references to existing state that will be
     * modified during scenario exploration.
     */
    ScenarioContext(
        double& scenario_proba_,
        Seq_type_str_p_map& constructed_sequences_,
        Seq_offsets_map& seq_offsets_,
        Mismatch_vectors_map& mismatches_lists_
    ) : scenario_proba(scenario_proba_),
        scenario_error_w_proba(0.0),  // Initialize to 0 (set at leaf nodes)
        constructed_sequences(constructed_sequences_),
        seq_offsets(seq_offsets_),
        mismatches_lists(mismatches_lists_)
    {}
    
    // Prevent copying (would create dangling references)
    ScenarioContext(const ScenarioContext&) = delete;
    ScenarioContext& operator=(const ScenarioContext&) = delete;
    
    // Allow moving (just rebinds references)
    ScenarioContext(ScenarioContext&&) = default;
    ScenarioContext& operator=(ScenarioContext&&) = default;
    
    // ========================================================================
    // Phase 3.5.3: Scenario State Management Abstractions
    // ========================================================================
    // Semantic operations for sequence construction and offset tracking.
    // These inline methods hide memory layer manipulation details and make
    // event code more readable.
    // 
    // PERFORMANCE: All inline - zero overhead compared to direct access.
    // ========================================================================
    
    /**
     * @brief Set sequence offset for a gene/junction
     * 
     * Wrapper for seq_offsets.set_value() with semantic naming.
     * 
     * @param seq_type Type of sequence (V_gene_seq, D_gene_seq, etc.)
     * @param side Which end (Five_prime or Three_prime)
     * @param offset Position in the read sequence
     * @param memory_layer Memory layer for multi-path exploration
     * 
     * PERFORMANCE: Inline, delegates to Seq_offsets_map
     */
    inline void set_offset(
        Seq_type seq_type,
        Seq_side side,
        Seq_Offset offset,
        size_t memory_layer
    ) {
        seq_offsets.set_value(seq_type, side, offset, memory_layer);
    }
    
    /**
     * @brief Get sequence offset for a gene/junction
     * 
     * Wrapper for seq_offsets.at() with semantic naming.
     * 
     * @param seq_type Type of sequence
     * @param side Which end (Five_prime or Three_prime)
     * @param memory_layer Memory layer to read from
     * @return Offset position
     * 
     * PERFORMANCE: Inline memory layer lookup
     */
    inline Seq_Offset get_offset(
        Seq_type seq_type,
        Seq_side side,
        size_t memory_layer
    ) const {
        return seq_offsets.at(seq_type, side, memory_layer);
    }
    
    /**
     * @brief Set constructed sequence segment
     * 
     * Stores a constructed sequence segment (gene or junction)
     * at the specified memory layer.
     * 
     * @param seq_type Type of sequence being stored
     * @param sequence Pointer to the sequence data
     * @param memory_layer Memory layer for storage
     * 
     * PERFORMANCE: Inline, delegates to Seq_type_str_p_map
     */
    inline void set_sequence_segment(
        Seq_type seq_type,
        const Int_Str_ptr sequence,
        size_t memory_layer
    ) {
        constructed_sequences.set_value(seq_type, sequence, memory_layer);
    }
    
    /**
     * @brief Get constructed sequence segment
     * 
     * Retrieves a previously stored sequence segment.
     * 
     * @param seq_type Type of sequence to retrieve
     * @param memory_layer Memory layer to read from
     * @return Pointer to sequence data
     * 
     * PERFORMANCE: Inline memory layer lookup
     */
    inline const Int_Str* get_sequence_segment(
        Seq_type seq_type,
        size_t memory_layer
    ) const {
        return constructed_sequences.at(seq_type, memory_layer);
    }
    
    /**
     * @brief Set mismatch list for a sequence type
     * 
     * Stores alignment mismatches for error rate calculation.
     * 
     * @param seq_type Type of sequence
     * @param mismatches Pointer to vector of mismatch positions
     * @param memory_layer Memory layer for storage
     * 
     * PERFORMANCE: Inline, delegates to Mismatch_vectors_map
     */
    inline void set_mismatches(
        Seq_type seq_type,
        std::vector<int>* mismatches,
        size_t memory_layer
    ) {
        mismatches_lists.set_value(seq_type, mismatches, memory_layer);
    }
    
    /**
     * @brief Get mismatch list for a sequence type
     * 
     * Retrieves stored alignment mismatches.
     * 
     * @param seq_type Type of sequence
     * @param memory_layer Memory layer to read from
     * @return Pointer to vector of mismatch positions
     * 
     * PERFORMANCE: Inline memory layer lookup
     */
    inline std::vector<int>* get_mismatches(
        Seq_type seq_type,
        size_t memory_layer
    ) const {
        return mismatches_lists.at(seq_type, memory_layer);
    }
    
    // ========================================================================
    // Phase 3.5.3: Overloads for implicit (current) memory layer access
    // ========================================================================
    // These overloads delegate to the underlying container's current layer,
    // supporting migration from legacy code that doesn't specify layers.
    // ========================================================================
    
    /**
     * @brief Get offset using current memory layer
     * @param seq_type Type of sequence
     * @param side Which end (Five_prime or Three_prime)
     * @return Offset at current memory layer
     */
    inline Seq_Offset get_offset(
        Seq_type seq_type,
        Seq_side side
    ) const {
        return seq_offsets.at(seq_type, side);
    }
    
    /**
     * @brief Get sequence segment using current memory layer
     * @param seq_type Type of sequence
     * @return Pointer to sequence at current memory layer
     */
    inline const Int_Str* get_sequence_segment(Seq_type seq_type) const {
        return constructed_sequences.at(seq_type);
    }
    
    /**
     * @brief Get mismatches using current memory layer
     * @param seq_type Type of sequence
     * @return Pointer to mismatch vector at current memory layer
     */
    inline std::vector<int>* get_mismatches(Seq_type seq_type) const {
        return mismatches_lists.at(seq_type);
    }
    
    // ========================================================================
    // Phase 3.5.3: Non-const overloads for mutable sequence access
    // ========================================================================
    // Insertion sequences (VD_ins_seq, DJ_ins_seq, VJ_ins_seq) need mutable
    // access because they get filled in during iteration. These non-const
    // overloads allow modifying the retrieved sequences.
    // ========================================================================
    
    /**
     * @brief Get mutable sequence segment (explicit memory layer)
     * @param seq_type Type of sequence
     * @param memory_layer Memory layer to read from
     * @return Mutable pointer to sequence data
     */
    inline Int_Str* get_sequence_segment(
        Seq_type seq_type,
        size_t memory_layer
    ) {
        return const_cast<Int_Str*>(constructed_sequences.at(seq_type, memory_layer));
    }
    
    /**
     * @brief Get mutable sequence segment (current memory layer)
     * @param seq_type Type of sequence
     * @return Mutable pointer to sequence at current memory layer
     */
    inline Int_Str* get_sequence_segment(Seq_type seq_type) {
        return const_cast<Int_Str*>(constructed_sequences.at(seq_type));
    }
};
