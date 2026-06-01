#pragma once

#include <igor/Core/Utils.h>
#include <igor/Core/ScenarioContext.h>
#include <vector>

/**
 * @brief Flattened, read-only view of a completed scenario
 * 
 * This struct provides a simplified interface for Counters to access
 * completed scenario data without exposing memory layer implementation.
 * Created by flattening ScenarioContext at leaf nodes.
 * 
 * Design: Zero-copy via pointers/references to underlying data.
 * Performance: Inline construction, no heap allocation.
 * 
 * ARCHITECTURE NOTE: Counters are passive observers that collect statistics.
 * They don't need to know about memory layer maps - just the final realized
 * scenario. This flattened view hides implementation details and provides
 * a simple, testable interface.
 */
struct Scenario {
    // Scenario probabilities
    const double scenario_proba;
    const long double scenario_error_w_proba;
    
    // Flattened sequence data (indexed by Seq_type)
    // Zero-copy: just stores pointers to existing data
    // Using std::vector (not std::array) for future flexibility if Seq_type becomes runtime-determined
    std::vector<const Int_Str*> sequences;
    
    // Flattened offsets (5'/3' pairs, indexed by Seq_type)
    struct OffsetPair {
        Seq_Offset five_prime;
        Seq_Offset three_prime;
    };
    std::vector<OffsetPair> offsets;
    
    // Flattened mismatches (indexed by Seq_type)
    std::vector<const std::vector<int>*> mismatches;
    
    /**
     * @brief Construct flattened view from ScenarioContext
     * 
     * Flattens memory layer maps to simple arrays indexed by Seq_type.
     * Uses current memory layer (at leaf nodes, all layers are resolved).
     * 
     * @param ctx ScenarioContext at leaf node (all layers resolved)
     * 
     * PERFORMANCE: Inline, zero-copy (stores pointers only)
     * Pre-allocates vectors to V_gene_seq_Max size for cache efficiency
     */
    inline Scenario(const ScenarioContext& ctx)
        : scenario_proba(ctx.scenario_proba),
          scenario_error_w_proba(ctx.scenario_error_w_proba),
          sequences(6, nullptr),  // Pre-reserve vector size with Seq_type cardinality for performance
          offsets(6),
          mismatches(6, nullptr)
    {
        // Flatten constructed_sequences map → vector
        // Uses current memory layer (leaf node has finalized state)
        // NOTE: Not all Seq_type values are always present (e.g., D_gene_seq absent in VJ recombination)
        for (size_t seq_type_idx = 0; seq_type_idx < 6; ++seq_type_idx) {
            Seq_type seq_type = static_cast<Seq_type>(seq_type_idx);
            if (ctx.constructed_sequences.exist(seq_type)) {
                sequences[seq_type_idx] = ctx.get_sequence_segment(seq_type);
            }
            // else: remains nullptr (already initialized)
        }
        
        // Flatten seq_offsets map → vector of pairs
        // NOTE: Offsets only exist for initialized sequences
        for (size_t seq_type_idx = 0; seq_type_idx < 6; ++seq_type_idx) {
            Seq_type seq_type = static_cast<Seq_type>(seq_type_idx);
            if (ctx.seq_offsets.exist(seq_type, Five_prime) &&  ctx.seq_offsets.exist(seq_type, Three_prime)) {
                offsets[seq_type_idx] = {
                    ctx.get_offset(seq_type, Five_prime),
                    ctx.get_offset(seq_type, Three_prime)
                };
            }
            // else: remains default-initialized (zero offsets)
        }
        
        // Flatten mismatches_lists map → vector
        // NOTE: Mismatches only tracked for gene sequences
        for (size_t seq_type_idx = 0; seq_type_idx < 6; ++seq_type_idx) {
            Seq_type seq_type = static_cast<Seq_type>(seq_type_idx);
            if (ctx.mismatches_lists.exist(seq_type)) {
                mismatches[seq_type_idx] = ctx.get_mismatches(seq_type);
            }
            // else: remains nullptr (already initialized)
        }
    }
    
    /**
     * @brief Get sequence for a specific Seq_type
     * 
     * @param seq_type Type of sequence (V_gene_seq, D_gene_seq, etc.)
     * @return Pointer to sequence data, or nullptr if not set
     */
    inline const Int_Str* get_sequence(Seq_type seq_type) const {
        return sequences[seq_type];
    }
    
    /**
     * @brief Get offsets for a specific Seq_type
     * 
     * @param seq_type Type of sequence
     * @return Pair of 5' and 3' offsets
     */
    inline OffsetPair get_offsets(Seq_type seq_type) const {
        return offsets[seq_type];
    }
    
    /**
     * @brief Get mismatches for a specific Seq_type
     * 
     * @param seq_type Type of sequence
     * @return Pointer to vector of mismatch positions, or nullptr if none
     */
    inline const std::vector<int>* get_mismatches(Seq_type seq_type) const {
        return mismatches[seq_type];
    }
    
    /**
     * @brief Build the full scenario sequence by concatenating segments in biological order
     * 
     * Concatenates sequence segments in the correct biological order:
     * - V → VD_ins → D → DJ_ins → J (if D gene present, VDJ recombination)
     * - V → VJ_ins → J (if D gene absent, VJ recombination)
     * 
     * TODO: Refactor Seq_type enum to reflect biological ordering instead of requiring
     * conditional branching logic. Ideal order: V_gene_seq(0), VD_ins_seq(1), D_gene_seq(2),
     * DJ_ins_seq(3), J_gene_seq(4), VJ_ins_seq(5). This would enable simple iteration
     * for sequence concatenation.
     * 
     * @param result Output parameter - cleared and filled with the complete sequence
     * 
     * PERFORMANCE: Inline, zero-copy (modifies output parameter by reference)
     */
    inline void build_full_sequence(Int_Str& result) const {
        result.clear();
        
        // Add V gene segment
        if (sequences[V_gene_seq]) {
            result += *sequences[V_gene_seq];
        }
        
        // Branch based on D gene presence (VDJ vs VJ recombination)
        if (sequences[D_gene_seq]) {
            // VDJ recombination path
            if (sequences[VD_ins_seq]) {
                result += *sequences[VD_ins_seq];
            }
            result += *sequences[D_gene_seq];
            if (sequences[DJ_ins_seq]) {
                result += *sequences[DJ_ins_seq];
            }
        } else {
            // VJ recombination path (no D gene)
            if (sequences[VJ_ins_seq]) {
                result += *sequences[VJ_ins_seq];
            }
        }
        
        // Add J gene segment
        if (sequences[J_gene_seq]) {
            result += *sequences[J_gene_seq];
        }
    }
    
    /**
     * @brief Build the full scenario sequence by concatenating segments in biological order
     * 
     * Convenience overload that returns by value. For performance-critical code,
     * prefer the reference-parameter version to avoid allocation.
     * 
     * @return The complete scenario sequence as Int_Str
     */
    inline Int_Str build_full_sequence() const {
        Int_Str result;
        build_full_sequence(result);
        return result;
    }
};
