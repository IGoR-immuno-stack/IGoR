/**
 * @file EventUtils.cpp
 * @brief Implementation of EventUtils functions including genetic code utilities.
 */

#include <igor/Core/EventUtils.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Aligner.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace EventUtils {

// ============================================================================
// Phase 1: Genetic Code Utilities — Internal helpers
// ============================================================================

namespace {

/**
 * @brief Compile-time IUPAC code from a nucleotide bitmask.
 *
 * @param bits Bitmask: bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T
 * @return IUPAC Int_nt code, or int_N if no clean mapping exists
 */
constexpr Int_nt iupac_from_bits_impl(int bits) {
    switch (bits) {
    case 0b0001: return int_A; // A
    case 0b0010: return int_C; // C
    case 0b0100: return int_G; // G
    case 0b1000: return int_T; // T
    case 0b0011: return int_M; // A|C
    case 0b0101: return int_R; // A|G
    case 0b1001: return int_W; // A|T
    case 0b0110: return int_S; // C|G
    case 0b1010: return int_Y; // C|T
    case 0b1100: return int_K; // G|T
    case 0b1101: return int_V; // A|C|G
    case 0b1011: return int_H; // A|C|T
    case 0b0111: return int_D; // A|G|T
    case 0b1110: return int_B; // C|G|T
    case 0b1111: return int_N; // A|C|G|T
    default:     return int_N; // fallback
    }
}

/**
 * @brief Standard genetic code lookup table.
 *
 * Index = codon_index(n0, n1, n2) = n0*16 + n1*4 + n2
 * where n0,n1,n2 ∈ {0,1,2,3} mapping to {A,C,G,T}.
 * Codon is read 5'→3': n0=first, n1=second, n2=third.
 */
constexpr std::array<char, 64> genetic_code_table() {
    std::array<char, 64> table{};
    table.fill('?');

    // AAA K  AAC N  AAG K  AAT N
    table[0] = 'K'; table[1] = 'N'; table[2] = 'K'; table[3] = 'N';
    // ACA T  ACC T  ACG T  ACT T
    table[4] = 'T'; table[5] = 'T'; table[6] = 'T'; table[7] = 'T';
    // AGA R  AGC S  AGG R  AGT S
    table[8] = 'R'; table[9] = 'S'; table[10] = 'R'; table[11] = 'S';
    // ATA I  ATC I  ATG M  ATT I
    table[12] = 'I'; table[13] = 'I'; table[14] = 'M'; table[15] = 'I';
    // CAA Q  CAC H  CAG Q  CAT H
    table[16] = 'Q'; table[17] = 'H'; table[18] = 'Q'; table[19] = 'H';
    // CCA P  CCC P  CCG P  CCT P
    table[20] = 'P'; table[21] = 'P'; table[22] = 'P'; table[23] = 'P';
    // CGA R  CGC R  CGG R  CGT R
    table[24] = 'R'; table[25] = 'R'; table[26] = 'R'; table[27] = 'R';
    // CTA L  CTC L  CTG L  CTT L
    table[28] = 'L'; table[29] = 'L'; table[30] = 'L'; table[31] = 'L';
    // GAA E  GAC D  GAG E  GAT D
    table[32] = 'E'; table[33] = 'D'; table[34] = 'E'; table[35] = 'D';
    // GCA A  GCC A  GCG A  GCT A
    table[36] = 'A'; table[37] = 'A'; table[38] = 'A'; table[39] = 'A';
    // GGA G  GGC G  GGG G  GGT G
    table[40] = 'G'; table[41] = 'G'; table[42] = 'G'; table[43] = 'G';
    // GTA V  GTC V  GTG V  GTT V
    table[44] = 'V'; table[45] = 'V'; table[46] = 'V'; table[47] = 'V';
    // TAA *  TAC Y  TAG *  TAT Y
    table[48] = '*'; table[49] = 'Y'; table[50] = '*'; table[51] = 'Y';
    // TCA S  TCC S  TCG S  TCT S
    table[52] = 'S'; table[53] = 'S'; table[54] = 'S'; table[55] = 'S';
    // TGA *  TGC C  TGG W  TGT C
    table[56] = '*'; table[57] = 'C'; table[58] = 'W'; table[59] = 'C';
    // TTA L  TTC F  TTG L  TTT F
    table[60] = 'L'; table[61] = 'F'; table[62] = 'L'; table[63] = 'F';

    return table;
}

// Compile-time initialized lookup table
constexpr auto g_genetic_code = genetic_code_table();

/**
 * @brief Get all nucleotide bits present at a given codon position.
 */
int nucleotides_at_position(CodonMask mask, int pos) {
    int bits = 0;
    for (int i = 0; i < 64; ++i) {
        if ((mask >> i) & 1) {
            // Decode codon index: n0 = i/16, n1 = (i/4)%4, n2 = i%4
            int nt;
            if (pos == 0) nt = i / 16;
            else if (pos == 1) nt = (i / 4) % 4;
            else nt = i % 4;
            bits |= (1 << nt);
        }
    }
    return bits;
}

} // anonymous namespace

// ============================================================================
// Phase 1: Genetic Code Utilities — Public API
// ============================================================================

int iupac_from_bits(int bits) {
    return static_cast<int>(iupac_from_bits_impl(bits));
}

char translate_codon(int n0, int n1, int n2) {
    if (n0 < 0 || n0 > 3 || n1 < 0 || n1 > 3 || n2 < 0 || n2 > 3) {
        return '?';
    }
    return g_genetic_code[codon_index(n0, n1, n2)];
}

CodonMask codon_mask_for_aa(char aa) {
    CodonMask mask = 0;
    for (int i = 0; i < 64; ++i) {
        if (g_genetic_code[i] == aa) {
            mask |= (1ULL << i);
        }
    }
    return mask;
}

std::string translate_int_seq(const Int_Str& seq, int frame_offset) {
    if (frame_offset < 0 || frame_offset >= static_cast<int>(seq.size())) {
        return {};
    }
    size_t remaining = seq.size() - static_cast<size_t>(frame_offset);
    if (remaining % 3 != 0) {
        return {}; // incomplete codon at end
    }
    std::string result;
    result.reserve(remaining / 3);
    for (size_t i = 0; i < remaining; i += 3) {
        result.push_back(translate_codon(seq[frame_offset + i],
                                         seq[frame_offset + i + 1],
                                         seq[frame_offset + i + 2]));
    }
    return result;
}

std::array<int, 3> aa_to_iupac_codon(char aa) {
    CodonMask mask = codon_mask_for_aa(aa);
    return mask_to_iupac_codon(mask);
}

std::array<int, 3> mask_to_iupac_codon(CodonMask mask) {
    std::array<int, 3> result{};
    for (int pos = 0; pos < 3; ++pos) {
        int bits = nucleotides_at_position(mask, pos);
        result[pos] = static_cast<int>(iupac_from_bits_impl(bits));
    }
    return result;
}

CodonMask motif_char_to_mask(char ch) {
    if (ch == 'X') {
        // Wildcard: all 61 non-stop codons
        CodonMask mask = 0;
        for (int i = 0; i < 64; ++i) {
            if (g_genetic_code[i] != '*') {
                mask |= (1ULL << i);
            }
        }
        return mask;
    }
    if (ch == '*') {
        return codon_mask_for_aa('*');
    }
    // Check if it's a valid amino acid letter
    if ((ch >= 'A' && ch <= 'Y') || ch == '*') {
        CodonMask m = codon_mask_for_aa(ch);
        if (m == 0) {
            throw std::invalid_argument(
                std::string("Unrecognized motif character: '") + ch + "'");
        }
        return m;
    }
    throw std::invalid_argument(
        std::string("Unrecognized motif character: '") + ch + "'");
}

std::vector<CodonMask> parse_aa_motif(const std::string& motif) {
    std::vector<CodonMask> result;
    result.reserve(motif.size());

    size_t i = 0;
    while (i < motif.size()) {
        char ch = motif[i];

        if (ch == '[') {
            // Bracket group: find closing bracket
            size_t close = motif.find(']', i + 1);
            if (close == std::string::npos) {
                throw std::invalid_argument(
                    "Unclosed bracket in motif: missing ']'");
            }
            std::string group = motif.substr(i + 1, close - i - 1);
            if (group.empty()) {
                throw std::invalid_argument(
                    "Empty bracket group in motif: '[]'");
            }

            // Check for duplicate characters
            std::vector<bool> seen(256, false);
            CodonMask group_mask = 0;
            for (char gc : group) {
                if (seen[static_cast<unsigned char>(gc)]) {
                    throw std::invalid_argument(
                        std::string("Duplicate character in bracket group: '") +
                        gc + "'");
                }
                seen[static_cast<unsigned char>(gc)] = true;
                group_mask |= motif_char_to_mask(gc);
            }
            result.push_back(group_mask);
            i = close + 1;
        } else {
            // Single character: AA code, wildcard, or stop
            result.push_back(motif_char_to_mask(ch));
            ++i;
        }
    }

    return result;
}

// ============================================================================
// Existing EventUtils API (unchanged)
// ============================================================================

GeneChoiceStatus check_gene_choice(
    Gene_class gene,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map,
    const std::unordered_set<Rec_Event_name> &processed_events) {

  GeneChoiceStatus status;
  status.exists = false;
  status.chosen = false;
  status.event_ptr = nullptr;

  auto key = std::make_tuple(GeneChoice_t, gene, Undefined_side);
  if (events_map.count(key) != 0) {
    status.exists = true;
    status.event_ptr = events_map.at(key);
    if (processed_events.count(status.event_ptr->get_name()) != 0) {
      status.chosen = true;
    }
  }
  return status;
}

Int_Str build_scenario_sequence(Seq_type_str_p_map &constructed_sequences,
                                bool has_v, bool has_d, bool has_j,
                                bool has_vd_ins, bool has_dj_ins,
                                bool has_vj_ins) {

  Int_Str scenario_resulting_sequence;

  if (has_v) {
    scenario_resulting_sequence += (*constructed_sequences[V_gene_seq]);
  }

  if (has_d) {
    if (has_vd_ins) {
      scenario_resulting_sequence += (*constructed_sequences[VD_ins_seq]);
    }
    scenario_resulting_sequence += (*constructed_sequences[D_gene_seq]);
    if (has_dj_ins) {
      scenario_resulting_sequence += (*constructed_sequences[DJ_ins_seq]);
    }
  } else {
    if (has_vj_ins) {
      scenario_resulting_sequence += (*constructed_sequences[VJ_ins_seq]);
    }
  }

  if (has_j) {
    scenario_resulting_sequence += (*constructed_sequences[J_gene_seq]);
  }

  return scenario_resulting_sequence;
}

int get_insertion_len_max(
    Gene_class gene_pair,
    const std::unordered_map<std::tuple<Event_type, Gene_class, Seq_side>,
                             std::shared_ptr<Rec_Event>> &events_map) {
  auto key = std::make_tuple(Insertion_t, gene_pair, Undefined_side);
  if (events_map.count(key) != 0) {
    return events_map.at(key)->get_len_max();
  }
  return 0;
}

void initialize_offset_memory(
    const std::vector<std::pair<std::shared_ptr<const Rec_Event>, int>>
        &offset_vector,
    Index_map &index_map,
    std::forward_list<std::tuple<int, int, int>> &memory_and_offsets) {

  for (auto iter = offset_vector.begin(); iter != offset_vector.end(); ++iter) {
    int event_identitfier = (*iter).first->get_event_identifier();
    index_map.request_memory_layer(event_identitfier);
    memory_and_offsets.emplace_front(
        event_identitfier,
        index_map.get_current_memory_layer(event_identitfier), (*iter).second);
  }
}

} // namespace EventUtils
