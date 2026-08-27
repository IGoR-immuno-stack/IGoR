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
 * - Optional patched/motif query journal (AA-level Pgen mode)
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
    // - In patched/motif mode: journaled_query->iupac_union (IUPAC-encoded)
    const Int_Str& int_sequence;

    // Genomic template alignments per gene class
    // Currently used to constrain gene choice realizations
    // Maps Gene_class → vector of Alignment_data (from Aligner)
    // Only V_gene, D_gene, J_gene keys are valid (not junction classes).
    //
    // FUTURE: May generalize to support realization constraints on
    // any RecEvent type, not just gene choices. Consider separating
    // alignment data from realization constraints in future refactoring.
    const std::unordered_map<Gene_class, std::vector<Alignment_data>>&
        gene_alignments;

    /**
     * @brief Patched/motif query journal (optional)
     *
     * Present only in patched/motif mode, i.e. when scoring an AA motif rather
     * than a concrete nucleotide sequence. Carries iupac_union,
     * iupac_intersection, empty_isect and the patch list.
     * std::nullopt for standard NT inference, which is unaffected by its presence.
     *
     * NOTE: this is the only owning member of an otherwise all-by-reference
     * struct. It is potentially large (three Int_Str of receptor length plus up
     * to 61 alternatives per patch). If query contexts ever become copyable or
     * shared across threads, hold it by shared_ptr<const JournaledQuery> instead.
     */
    const std::optional<JournaledQuery> journaled_query;

    /**
     * @brief Constructor - binds const references to input data (NT mode)
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
     * @param sequence_        Nucleotide rendering of the query (IUPAC union)
     * @param int_sequence_    journaled_query's iupac_union as Int_Str
     * @param gene_alignments_ Gene alignments
     * @param jq_              Journal with patches, iupac_union, iupac_intersection
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
