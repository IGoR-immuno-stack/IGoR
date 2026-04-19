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
    
    // ===== TEMPORARY ADAPTER INFRASTRUCTURE =====
    // TODO Phase 4.5: Remove after all counters migrated to new interface
    // These references enable NEW→OLD adapter in Counter.cpp for gradual migration
    // Once all counter subclasses override new interface, remove these fields
    Seq_type_str_p_map& constructed_sequences_ref;
    const Seq_offsets_map& seq_offsets_ref;
    Mismatch_vectors_map& mismatches_lists_ref;
    
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
          mismatches(6, nullptr),
          // Store references to original containers (for NEW→OLD adapter)
          constructed_sequences_ref(ctx.constructed_sequences),
          seq_offsets_ref(ctx.seq_offsets),
          mismatches_lists_ref(ctx.mismatches_lists)
    {
        // TODO Phase 4.4: Flatten maps for new counter implementations
        // For now, we skip flattening since legacy counters only use the map references.
        // When counters are migrated to use the flattened vectors (sequences, offsets, mismatches),
        // uncomment and test the flattening code below.
        
        // Flattening disabled to avoid segfaults during Phase 4.2-4.3 (legacy counter compatibility period)
        // The NEW→OLD adapter in Counter.cpp only uses the *_ref fields, not the flattened vectors.
        
        /*
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
        */
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
};
