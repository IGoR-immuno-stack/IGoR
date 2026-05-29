/**
 * @file JournaledQuery.cpp
 * @brief Implementation of motif_to_journaled_query factory function.
 *
 * Converts an AA motif string (e.g. "CAVX[KSM]DS") into a JournaledQuery
 * with properly computed reference, iupac_union, iupac_intersection,
 * empty_isect, and patches fields.
 */

#include <igor/Core/JournaledQuery.h>
#include <igor/Core/EventUtils.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

using namespace EventUtils;
using namespace std;

namespace {

/**
 * @brief Enumerate all valid concrete codons (NT triplets) for a CodonMask.
 *
 * @param mask CodonMask (bit k set iff codon_index(k) is valid)
 * @return Vector of Int_Str, each of length 3, containing all valid codons
 */
vector<Int_Str> enumerate_codons(CodonMask mask) {
    vector<Int_Str> result;
    result.reserve(64);
    for (int i = 0; i < 64; ++i) {
        if ((mask >> i) & 1) {
            // Decode codon index to nucleotides
            int n0 = i / 16;
            int n1 = (i / 4) % 4;
            int n2 = i % 4;
            result.push_back({n0, n1, n2});
        }
    }
    return result;
}

/**
 * @brief Compute per-position IUPAC union and intersection for a list of codons.
 *
 * For each codon position (0, 1, 2):
 * - Union: nucleotides that appear at that position in ANY codon
 * - Intersection: nucleotides that appear at that position in ALL codons
 *
 * @param codons Vector of Int_Str (each of length 3)
 * @param union_result Output: per-position IUPAC union codes
 * @param inter_result Output: per-position IUPAC intersection codes
 * @param empty_isect Output: true where intersection is empty
 */
void compute_iupac_from_codons(
    const vector<Int_Str>& codons,
    array<int, 3>& union_result,
    array<int, 3>& inter_result,
    array<bool, 3>& empty_isect)
{
    // Initialize position-wise nucleotide sets
    array<int, 3> union_bits = {0, 0, 0};  // bitmask: bit j set if nucleotide j appears
    array<int, 3> inter_bits = {0b1111, 0b1111, 0b1111}; // start with all nucleotides

    for (const auto& codon : codons) {
        for (int pos = 0; pos < 3; ++pos) {
            union_bits[pos] |= (1 << codon[pos]);
            inter_bits[pos] &= (1 << codon[pos]);
        }
    }

    // Encode as IUPAC codes
    for (int pos = 0; pos < 3; ++pos) {
        union_result[pos] = EventUtils::iupac_from_bits(union_bits[pos]);
        if (inter_bits[pos] == 0) {
            inter_result[pos] = static_cast<int>(int_N); // placeholder
            empty_isect[pos] = true;
        } else {
            inter_result[pos] = EventUtils::iupac_from_bits(inter_bits[pos]);
            empty_isect[pos] = false;
        }
    }
}

} // anonymous namespace

namespace EventUtils {

JournaledQuery motif_to_journaled_query(
    const string& motif,
    int frame_offset,
    int receptor_len)
{
    JournaledQuery jq;
    jq.display_string = motif;

    // Allocate fields to receptor length, initialized to defaults
    jq.reference.resize(receptor_len, static_cast<int>(int_N));
    jq.iupac_union.resize(receptor_len, static_cast<int>(int_N));
    jq.iupac_intersection.resize(receptor_len, static_cast<int>(int_N));
    jq.empty_isect.resize(receptor_len, false);

    // Parse motif into CodonMasks (one per codon position)
    vector<CodonMask> masks = parse_aa_motif(motif);

    // Process each codon position
    for (size_t codon_idx = 0; codon_idx < masks.size(); ++codon_idx) {
        int nt_start = frame_offset + static_cast<int>(codon_idx) * 3;

        // Enumerate all valid concrete codons for this position
        vector<Int_Str> codons = enumerate_codons(masks[codon_idx]);

        if (codons.empty()) {
            throw runtime_error(
                "motif_to_journaled_query: no codons found for motif position " +
                to_string(codon_idx));
        }

        // Pick the first valid codon as reference
        Int_Str ref_codon = codons[0];
        for (int j = 0; j < 3; ++j) {
            jq.reference[nt_start + j] = ref_codon[j];
        }

        // Compute per-position IUPAC union and intersection
        array<int, 3> iupac_union, iupac_intersection;
        array<bool, 3> isect_empty;
        compute_iupac_from_codons(codons, iupac_union, iupac_intersection, isect_empty);

        for (int j = 0; j < 3; ++j) {
            jq.iupac_union[nt_start + j] = iupac_union[j];
            jq.iupac_intersection[nt_start + j] = iupac_intersection[j];
            jq.empty_isect[nt_start + j] = isect_empty[j];
        }

        // Build patch with alternatives (exclude reference)
        if (codons.size() > 1) {
            Patch patch;
            patch.start = nt_start;
            patch.length = 3;
            for (size_t c = 1; c < codons.size(); ++c) {
                patch.alternatives.push_back(codons[c]);
            }
            jq.patches.push_back(patch);
        }
        // If codons.size() == 1, no patch needed (trivial codon)
    }

    return jq;
}

} // namespace EventUtils
