#pragma once

#include <igor/Core/Utils.h>
#include <igor/Core/IntStr.h>
#include <igor/Core/Aligner.h>
#include <igor/Core/JournaledQuery.h>
#include <optional>
#include <unordered_map>
#include <vector>
#include <string>

/**
 * @brief Encapsulates the input sequence being analyzed
 *
 * QuerySequenceContext holds all input data for a single sequence:
 * - Raw nucleotide sequence
 * - Integer-encoded sequence
 * - Genomic alignment results (currently used for gene choice constraints)
 * - Optional AA-level generation probability journal (patched/motif mode)
 *
 * Separating this enables:
 * - Clear batch processing (iterate over queries with same model)
 * - Parallelization (share ModelContext across threads)
 * - Testing with mock inputs
 *
 * All members are const references - the query is never modified.
 *
 * DESIGN NOTE - Future Extensibility:
 * Currently, gene_alignments specifically stores alignment data from the
 * Aligner and is used to constrain gene choice realizations. This mixes
 * two concepts:
 *   1. Sequence alignment results (Alignment_data structures)
 *   2. Restricting which realizations a RecEvent can explore
 *
 * In the future, we may want to generalize this so that ANY RecEvent
 * (not just gene choices) can have realization constraints. One approach:
 *   - Keep alignment data separate (pure input)
 *   - Add a general "allowed_realizations" map per event type
 *   - Gene choices would derive constraints from alignment data
 *
 * For now, we keep the current structure but name it explicitly as
 * "gene_alignments" to clarify its current scope.
 */
struct QuerySequenceContext {
    // Input sequence (nucleotide string)
    const std::string& sequence;

    // Input sequence (integer-encoded for efficient comparison)
    // - In standard NT mode: exact nucleotide sequence
    // - In patched/motif mode: iupac_union (IUPAC-encoded query)
    const Int_Str& int_sequence;

    // Genomic template alignments per gene class
    // Currently used to constrain gene choice realizations
    // Maps Gene_class → vector of Alignment_data (from Aligner)
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>&
        gene_alignments;

    /**
     * @brief AA-level generation probability journal (optional)
     *
     * Present only in patched/motif mode (when an AA motif is being
     * scored). Contains the JournaledQuery with iupac_union, patches,
     * and intersection information needed for two-track mismatch
     * recording. std::nullopt for standard NT inference.
     */
    const std::optional<JournaledQuery> journaled_query;

    /**
     * @brief Standard NT mode constructor (unchanged behaviour)
     */
    QuerySequenceContext(
        const std::string& sequence_,
        const Int_Str& int_sequence_,
        const std::unordered_map<Gene_class, std::vector<Alignment_data>>& gene_alignments_
    ) : sequence(sequence_),
        int_sequence(int_sequence_),
        gene_alignments(gene_alignments_),
        journaled_query(std::nullopt)
    {}

    /**
     * @brief Patched/motif mode constructor
     *
     * @param sequence_ Amino acid motif string (e.g. "MW[K]X")
     * @param int_sequence_ iupac_union of the motif as Int_Str
     * @param gene_alignments_ Gene alignments
     * @param jq_ JournaledQuery with patches, iupac_union, iupac_intersection
     */
    QuerySequenceContext(
        const std::string& sequence_,
        const Int_Str& int_sequence_,
        const std::unordered_map<Gene_class, std::vector<Alignment_data>>& gene_alignments_,
        JournaledQuery jq_
    ) : sequence(sequence_),
        int_sequence(int_sequence_),
        gene_alignments(gene_alignments_),
        journaled_query(std::move(jq_))
    {}

    // Prevent copying and moving (const references)
    QuerySequenceContext(const QuerySequenceContext&) = delete;
    QuerySequenceContext& operator=(const QuerySequenceContext&) = delete;
    QuerySequenceContext(QuerySequenceContext&&) = delete;
    QuerySequenceContext& operator=(QuerySequenceContext&&) = delete;
};
