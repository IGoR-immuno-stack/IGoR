#pragma once

#include <array>
#include <cstdint>
#include <forward_list>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Utils.h>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <igor/Core/IntStr.h>
#include <igorCoreExport.h>

namespace EventUtils {

// ============================================================================
// Phase 1: Genetic Code Utilities
// ============================================================================

/**
 * @brief Codon index encoding: codon_index(n0, n1, n2) = n0*16 + n1*4 + n2
 *
 * Maps each of the 64 possible nucleotide triplets to a unique index in [0, 63].
 * Nucleotide values must be in {int_A=0, int_C=1, int_G=2, int_T=3}.
 *
 * @param n0 First nucleotide position (5')
 * @param n1 Second nucleotide position (middle)
 * @param n2 Third nucleotide position (3')
 * @return Integer index in [0, 63]
 */
constexpr int codon_index(int n0, int n1, int n2) {
    return n0 * 16 + n1 * 4 + n2;
}

/**
 * @brief 64-bit bitmask representing a set of codons.
 *
 * Bit k is set iff codon k (encoded by codon_index) is in the set.
 * Used internally for AA motif parsing and JournaledQuery construction.
 */
using CodonMask = uint64_t;

/**
 * @brief Translate a codon (3 nucleotides) to its amino acid.
 *
 * Uses a constexpr lookup table based on the standard genetic code.
 *
 * @param n0 First nucleotide (int_A..int_T, values 0..3)
 * @param n1 Second nucleotide (int_A..int_T, values 0..3)
 * @param n2 Third nucleotide (int_A..int_T, values 0..3)
 * @return Amino acid character ('A'..'Y', '*' for stop), or '?' for invalid input
 */
CORE_EXPORT char translate_codon(int n0, int n1, int n2);

/**
 * @brief Get a CodonMask for a given amino acid.
 *
 * @param aa Single-letter amino acid code ('A'..'Y') or '*' for stop.
 * @return CodonMask with bits set for all codons encoding the target AA.
 *         Returns 0 for unrecognized characters.
 */
CORE_EXPORT CodonMask codon_mask_for_aa(char aa);

/**
 * @brief Translate an Int_Str nucleotide sequence to amino acids.
 *
 * Reads codons starting at frame_offset. Returns empty string if the
 * remaining length past frame_offset is not a multiple of 3.
 *
 * @param seq Nucleotide sequence (values must be int_A..int_T)
 * @param frame_offset Starting position for codon reading
 * @return Amino acid string, or empty if incomplete codons at end
 */
CORE_EXPORT std::string translate_int_seq(const Int_Str& seq, int frame_offset);

/**
 * @brief Get the conservative IUPAC nucleotide triplet for an amino acid.
 *
 * For most AAs, returns a single IUPAC triplet that exactly represents
 * the synonymous codon set (e.g., Phe = TTY, Gly = GGN).
 * For Leu/Arg/Ser, this over-approximates (e.g., Leu = YTN, which
 * also includes Phe codons). Exact handling uses CodonMask via
 * mask_to_iupac_codon.
 *
 * @param aa Single-letter amino acid code
 * @return Array of 3 IUPAC nucleotide codes [pos0, pos1, pos2]
 */
CORE_EXPORT std::array<int, 3> aa_to_iupac_codon(char aa);

/**
 * @brief Get the IUPAC nucleotide triplet for a CodonMask.
 *
 * For each codon position, collects the union of nucleotides that appear
 * at that position across all codons whose bit is set in the mask.
 * Encodes the union as an IUPAC code.
 *
 * @param mask Bitmask of codons
 * @return Array of 3 IUPAC nucleotide codes [pos0, pos1, pos2]
 */
CORE_EXPORT std::array<int, 3> mask_to_iupac_codon(CodonMask mask);

/**
 * @brief Convert a motif character to its CodonMask.
 *
 * - Single-letter AA code (e.g., 'M'): codon_mask_for_aa(ch)
 * - Wildcard 'X': OR of all non-stop codon masks (61 bits)
 * - Stop '*' : 3 stop codons
 * - Unrecognized: throws std::invalid_argument
 *
 * @param ch Motif character
 * @return CodonMask for the character
 */
CORE_EXPORT CodonMask motif_char_to_mask(char ch);

/**
 * @brief Parse an OLGA-style AA motif string into CodonMasks.
 *
 * Supports:
 * - Single-letter AA codes (e.g., "ML")
 * - Wildcard 'X' (any amino acid)
 * - Bracket groups (e.g., "[KSM]" = any of K, S, or M)
 *
 * Example: "CAV[KSM]DS" → {mask_C, mask_A, mask_V, mask_K|S|M, mask_D, mask_S}
 *
 * @param motif AA motif string
 * @return Vector of CodonMasks, one per position
 * @throws std::invalid_argument on malformed input
 */
CORE_EXPORT std::vector<CodonMask> parse_aa_motif(const std::string& motif);

/**
 * @brief Get the IUPAC code for a set of nucleotides (bitmask of 4 bits).
 *
 * Internal helper used by mask_to_iupac_codon and JournaledQuery construction.
 *
 * @param bits Bitmask: bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T
 * @return IUPAC Int_nt code, or int_N if no clean mapping exists
 */
CORE_EXPORT int iupac_from_bits(int bits);

// ============================================================================
// Existing EventUtils API (unchanged)
// ============================================================================

struct GeneChoiceStatus {
  bool exists;
  bool chosen;
  std::shared_ptr<const Rec_Event> event_ptr;
};

CORE_EXPORT GeneChoiceStatus check_gene_choice(
    Gene_class gene,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events);

CORE_EXPORT Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins);

CORE_EXPORT void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets);

CORE_EXPORT int get_insertion_len_max(
    Gene_class gene_pair,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map);

} // namespace EventUtils
