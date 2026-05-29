/**
 * @file JournaledQuery.h
 * @brief JournaledQuery and Patch structs for AA-level Pgen computation.
 *
 * The JournaledQuery struct unifies exact NT queries, IUPAC NT queries,
 * AA motifs, and arbitrary NT-patch queries under one representation.
 * It is used by the iterate loop to track both floor (conservative) and
 * upper-bound mismatch semantics simultaneously.
 *
 * Patch describes a span of positions where the query has multiple valid
 * NT alternatives (not just the reference).
 */

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <igor/Core/IntStr.h>
#include <igor/Core/Utils.h>
#include <igorCoreExport.h>

/// One span of positions with multiple valid NT alternatives.
/// Covers arbitrary-length patches (e.g. one codon = 3 positions).
struct Patch {
    /// First position within the full receptor nt sequence
    int start;

    /// Number of positions in this patch
    int length;

    /// Valid NT subsequences (each of length `length`).
    /// Does NOT include the reference subsequence.
    /// Used at leaf nodes for exact sequence verification.
    std::vector<Int_Str> alternatives;
};

/// Query representation for Pgen computation — exact NT, IUPAC NT, AA motif, or patched NT.
///
/// Field semantics:
/// - reference: a single concrete valid NT sequence. Used for alignment and
///   for leaf-node display. For AA motifs, this is an arbitrary representative
///   encoding (e.g., first valid codon per position).
/// - iupac_union: per-position IUPAC union across all valid NT sequences.
///   Used for alignment and for the floor track (pruning_mismatch_floor).
/// - iupac_intersection: per-position IUPAC intersection across all valid NT
///   sequences. Used for the upper-bound track (mismatches_lists).
/// - empty_isect[p]: true when the intersection at position p is empty
///   (no IUPAC code can represent it). Position is unconditionally flagged
///   as an upper-bound mismatch.
/// - patches: list of spans where the query has alternatives. Empty for exact
///   NT queries (all fields are then trivially identical).
/// - display_string: the original human-readable input, for diagnostics.
struct CORE_EXPORT JournaledQuery {
    /// One concrete valid NT sequence (leaf check + display)
    Int_Str reference;

    /// Per-pos IUPAC union of all alternatives (alignment + floor)
    Int_Str iupac_union;

    /// Per-pos IUPAC intersection of all alternatives (upper bound)
    Int_Str iupac_intersection;

    /// True where intersection is empty → always flag upper-bound mismatch
    std::vector<bool> empty_isect;

    /// Spans where the query has alternatives (empty for exact NT queries)
    std::vector<Patch> patches;

    /// Original input (AA motif, NT, etc.) for diagnostics
    std::string display_string;
};

namespace EventUtils {

/**
 * @brief Build a JournaledQuery from an OLGA-style AA motif string.
 *
 * @param motif        AA motif, e.g. "CAVX[KSM]DS" (single letters, 'X', or '[...]' groups)
 * @param frame_offset Position of the first codon start within the receptor nt sequence
 * @param receptor_len Total expected nt length of the receptor
 * @return JournaledQuery with reference, iupac_union, iupac_intersection,
 *         empty_isect, and patches populated.
 */
CORE_EXPORT JournaledQuery motif_to_journaled_query(
    const std::string& motif,
    int frame_offset,
    int receptor_len
);

/**
 * @brief Convenience overload for a plain AA sequence (no brackets).
 *
 * Equivalent to motif_to_journaled_query(aa_seq, frame_offset, receptor_len).
 */
inline JournaledQuery aa_to_journaled_query(
    const std::string& aa_seq,
    int frame_offset,
    int receptor_len)
{
    return motif_to_journaled_query(aa_seq, frame_offset, receptor_len);
}

} // namespace EventUtils
