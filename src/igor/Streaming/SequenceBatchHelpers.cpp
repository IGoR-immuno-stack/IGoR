/*
 * SequenceBatchHelpers.cpp
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#include <igor/Streaming/SequenceBatchHelpers.h>

#include <sparrow.hpp>
#include <sparrow/builder.hpp>
#include <sparrow/layout/array_access.hpp>

#include <sstream>
#include <algorithm>

namespace igor {

// Helper function implementations

bool has_column(const sparrow::record_batch &batch, const std::string &column_name)
{
    try {
        // Try to get the column - if it doesn't exist, this will throw or return invalid
        const auto &names = batch.names();
        return std::find(names.begin(), names.end(), column_name) != names.end();
    } catch (...) {
        return false;
    }
}

std::string get_string_value(const sparrow::record_batch &batch, const std::string &column_name,
                             size_t row_index, const std::string &default_value)
{

    if (!has_column(batch, column_name)) {
        return default_value;
    }

    try {
        const auto &column = batch.get_column(column_name);
        auto value = column[row_index];

        // Use visitor pattern to extract string from nullable_variant
        return std::visit([&default_value](auto&& arg) -> std::string {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, sparrow::nullable<std::string_view>>) {
                if (arg.has_value()) {
                    return std::string(arg.get());
                }
            }
            return default_value;
        }, value);
    } catch (...) {
        return default_value;
    }
}

int get_int_value(const sparrow::record_batch &batch, const std::string &column_name,
                  size_t row_index, int default_value)
{

    if (!has_column(batch, column_name)) {
        return default_value;
    }

    try {
        const auto &column = batch.get_column(column_name);
        auto value = column[row_index];

        // Use visitor pattern to extract int from nullable_variant
        return std::visit([&default_value](auto&& arg) -> int {
            using T = std::decay_t<decltype(arg)>;
            // Handle various integer types that sparrow might use
            if constexpr (std::is_same_v<T, sparrow::nullable<const int&>> ||
                          std::is_same_v<T, sparrow::nullable<int>> ||
                          std::is_same_v<T, sparrow::nullable<const int32_t&>> ||
                          std::is_same_v<T, sparrow::nullable<int32_t>> ||
                          std::is_same_v<T, sparrow::nullable<const long long&>> ||
                          std::is_same_v<T, sparrow::nullable<long long>>) {
                if (arg.has_value()) {
                    return static_cast<int>(arg.get());
                }
            }
            return default_value;
        }, value);
    } catch (...) {
        return default_value;
    }
}

double get_double_value(const sparrow::record_batch &batch, const std::string &column_name,
                        size_t row_index, double default_value)
{

    if (!has_column(batch, column_name)) {
        return default_value;
    }

    try {
        const auto &column = batch.get_column(column_name);
        auto value = column[row_index];

        // Use visitor pattern to extract double from nullable_variant
        return std::visit([&default_value](auto&& arg) -> double {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, sparrow::nullable<const double&>> ||
                          std::is_same_v<T, sparrow::nullable<double>> ||
                          std::is_same_v<T, sparrow::nullable<const float&>> ||
                          std::is_same_v<T, sparrow::nullable<float>>) {
                if (arg.has_value()) {
                    return static_cast<double>(arg.get());
                }
            }
            return default_value;
        }, value);
    } catch (...) {
        return default_value;
    }
}

size_t get_size_t_value(const sparrow::record_batch &batch, const std::string &column_name,
                        size_t row_index, size_t default_value)
{

    if (!has_column(batch, column_name)) {
        return default_value;
    }

    try {
        const auto &column = batch.get_column(column_name);
        auto value = column[row_index];

        // Use visitor pattern to extract size_t from nullable_variant
        return std::visit([&default_value](auto&& arg) -> size_t {
            using T = std::decay_t<decltype(arg)>;
            // Handle various unsigned integer types that sparrow might use
            if constexpr (std::is_same_v<T, sparrow::nullable<const uint64_t&>> ||
                          std::is_same_v<T, sparrow::nullable<uint64_t>> ||
                          std::is_same_v<T, sparrow::nullable<const size_t&>> ||
                          std::is_same_v<T, sparrow::nullable<size_t>> ||
                          std::is_same_v<T, sparrow::nullable<const unsigned long long&>> ||
                          std::is_same_v<T, sparrow::nullable<unsigned long long>> ||
                          std::is_same_v<T, sparrow::nullable<const uint32_t&>> ||
                          std::is_same_v<T, sparrow::nullable<uint32_t>>) {
                if (arg.has_value()) {
                    return static_cast<size_t>(arg.get());
                }
            }
            return default_value;
        }, value);
    } catch (...) {
        return default_value;
    }
}

std::vector<int> get_int_list_value(const sparrow::record_batch &batch,
                                     const std::string &column_name, size_t row_index)
{
    // Extract value from native Arrow list array
    if (!has_column(batch, column_name)) {
        return std::vector<int>{};
    }

    try {
        const auto& col = batch.get_column(column_name);

        // Access the array element at row_index
        auto list_element = col[row_index];

        std::vector<int> result;

        // Check if the value is null
        if (!list_element.has_value()) {
            return result;
        }

        // Get list_value from the nullable_variant
        // We know it's a list array, so we can try to get the list_value variant
        using list_nullable = sparrow::nullable<sparrow::list_value>;
        auto* list_ptr = std::get_if<list_nullable>(&list_element);

        if (list_ptr && list_ptr->has_value()) {
            const auto& list_val = list_ptr->value();
            // Iterate over the list elements
            for (const auto& item : list_val) {
                if (!item.has_value()) {
                    continue;  // Skip null elements
                }

                // Try to extract as int32 - the element should be in one of the integer variants
                using int32_ref = sparrow::nullable<const int&>;
                auto* int_ptr = std::get_if<int32_ref>(&item);
                if (int_ptr && int_ptr->has_value()) {
                    result.push_back(static_cast<int>(int_ptr->value()));
                }
            }
        }

        return result;
    } catch (...) {
        return std::vector<int>{};
    }
}

std::unordered_map<Gene_class, std::vector<Alignment_data>>
parse_alignments_from_columns(const sparrow::record_batch &batch, size_t row_index)
{

    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

    // Look for common alignment column patterns for V, D, and J genes:
    // - v_gene_name, v_gene_offset, v_gene_score, etc.
    // - d_gene_name, d_gene_offset, d_gene_score, etc.
    // - j_gene_name, j_gene_offset, j_gene_score, etc.

    const std::vector<std::pair<std::string, Gene_class>> gene_prefixes = { { "v_gene", V_gene },
                                                                            { "d_gene", D_gene },
                                                                            { "j_gene", J_gene } };

    for (const auto &[prefix, gene_class] : gene_prefixes) {
        std::string name_col = prefix + "_name";
        std::string offset_col = prefix + "_offset";
        std::string score_col = prefix + "_score";
        std::string five_p_offset_col = prefix + "_five_p_offset";
        std::string three_p_offset_col = prefix + "_three_p_offset";
        std::string insertions_col = prefix + "_insertions";
        std::string deletions_col = prefix + "_deletions";
        std::string align_length_col = prefix + "_align_length";
        std::string mismatches_col = prefix + "_mismatches";

        if (has_column(batch, name_col)) {
            std::string gene_name = get_string_value(batch, name_col, row_index);

            if (!gene_name.empty()) {
                // Extract all 9 fields of Alignment_data
                int offset = get_int_value(batch, offset_col, row_index, 0);
                size_t five_p_offset = get_size_t_value(batch, five_p_offset_col, row_index, 0);
                size_t three_p_offset = get_size_t_value(batch, three_p_offset_col, row_index, 0);
                size_t align_length = get_size_t_value(batch, align_length_col, row_index, 0);
                double score = get_double_value(batch, score_col, row_index, 0.0);

                // Extract nested list structures
                std::vector<int> insertions_vec = get_int_list_value(batch, insertions_col, row_index);
                std::vector<int> deletions_vec = get_int_list_value(batch, deletions_col, row_index);
                std::vector<int> mismatches_vec = get_int_list_value(batch, mismatches_col, row_index);

                // Convert vectors to forward_list for insertions and deletions
                std::forward_list<int> insertions(insertions_vec.begin(), insertions_vec.end());
                std::forward_list<int> deletions(deletions_vec.begin(), deletions_vec.end());

                // Create complete alignment data structure with all 9 fields
                Alignment_data align(gene_name, offset, five_p_offset, three_p_offset,
                                   align_length, insertions, deletions, mismatches_vec, score);

                alignments[gene_class].push_back(align);
            }
        }
    }

    return alignments;
}

SequenceData row_to_sequence_data(const sparrow::record_batch &batch, size_t row_index)
{

    // Validate row index
    if (row_index >= batch.nb_rows()) {
        throw std::out_of_range("row_index " + std::to_string(row_index)
                                + " out of range (batch has " + std::to_string(batch.nb_rows())
                                + " rows)");
    }

    SequenceData seq_data;

    // Get sequence ID - try multiple possible column names
    if (has_column(batch, "sequence_id")) {
        seq_data.index =
                get_int_value(batch, "sequence_id", row_index, static_cast<int>(row_index));
    } else if (has_column(batch, "seq_index")) {
        seq_data.index = get_int_value(batch, "seq_index", row_index, static_cast<int>(row_index));
    } else {
        // Use row index as fallback
        seq_data.index = static_cast<int>(row_index);
    }

    // Get sequence - this is required
    if (!has_column(batch, "sequence")) {
        throw std::runtime_error("Required column 'sequence' not found in record_batch");
    }

    seq_data.sequence = get_string_value(batch, "sequence", row_index);

    if (seq_data.sequence.empty()) {
        throw std::runtime_error("Empty sequence at row " + std::to_string(row_index));
    }

    // Parse alignments from columns
    seq_data.alignments = parse_alignments_from_columns(batch, row_index);

    return seq_data;
}

sparrow::record_batch vector_to_batch(
        const std::vector<std::tuple<int, std::string,
                                     std::unordered_map<Gene_class, std::vector<Alignment_data>>>>
                &sequences)
{

    if (sequences.empty()) {
        // Return empty batch with proper schema
        std::vector<std::string> names = { "sequence_id", "sequence" };
        std::vector<sparrow::array> arrays;

        // Create empty arrays using the builder and wrap in sparrow::array
        arrays.emplace_back(sparrow::build(std::vector<int32_t>{}));
        arrays.emplace_back(sparrow::build(std::vector<std::string>{}));

        return sparrow::record_batch(std::move(names), std::move(arrays));
    }

    // Build column arrays
    std::vector<int32_t> ids;
    std::vector<std::string> seqs;

    ids.reserve(sequences.size());
    seqs.reserve(sequences.size());

    // Prepare alignment data containers for each gene class
    // List fields (insertions, deletions, mismatches) use native Arrow list<int32> arrays
    std::unordered_map<Gene_class, std::vector<std::string>> v_gene_names;
    std::unordered_map<Gene_class, std::vector<int32_t>> v_offsets;
    std::unordered_map<Gene_class, std::vector<uint64_t>> v_five_p_offsets;
    std::unordered_map<Gene_class, std::vector<uint64_t>> v_three_p_offsets;
    std::unordered_map<Gene_class, std::vector<uint64_t>> v_align_lengths;
    std::unordered_map<Gene_class, std::vector<double>> v_scores;
    std::unordered_map<Gene_class, std::vector<std::vector<int32_t>>> v_insertions;  // Native list arrays
    std::unordered_map<Gene_class, std::vector<std::vector<int32_t>>> v_deletions;   // Native list arrays
    std::unordered_map<Gene_class, std::vector<std::vector<int32_t>>> v_mismatches;  // Native list arrays

    // Preallocate for gene classes we expect
    for (auto gc : { V_gene, D_gene, J_gene }) {
        v_gene_names[gc].reserve(sequences.size());
        v_offsets[gc].reserve(sequences.size());
        v_five_p_offsets[gc].reserve(sequences.size());
        v_three_p_offsets[gc].reserve(sequences.size());
        v_align_lengths[gc].reserve(sequences.size());
        v_scores[gc].reserve(sequences.size());
        v_insertions[gc].reserve(sequences.size());
        v_deletions[gc].reserve(sequences.size());
        v_mismatches[gc].reserve(sequences.size());
    }

    for (const auto &seq_tuple : sequences) {
        ids.push_back(static_cast<int32_t>(std::get<0>(seq_tuple)));
        seqs.push_back(std::get<1>(seq_tuple));

        // Extract alignments
        const auto &alignments = std::get<2>(seq_tuple);

        // For each gene class, add alignment data (or empty values if no alignments)
        for (auto gc : { V_gene, D_gene, J_gene }) {
            auto it = alignments.find(gc);
            if (it != alignments.end() && !it->second.empty()) {
                // Take the first (best) alignment for this gene class
                const auto &align = it->second[0];

                v_gene_names[gc].push_back(align.gene_name);
                v_offsets[gc].push_back(static_cast<int32_t>(align.offset));
                v_five_p_offsets[gc].push_back(static_cast<uint64_t>(align.five_p_offset));
                v_three_p_offsets[gc].push_back(static_cast<uint64_t>(align.three_p_offset));
                v_align_lengths[gc].push_back(static_cast<uint64_t>(align.align_length));
                v_scores[gc].push_back(align.score);

                // Convert forward_list to vector for list arrays
                std::vector<int32_t> ins_vec(align.insertions.begin(), align.insertions.end());
                std::vector<int32_t> del_vec(align.deletions.begin(), align.deletions.end());
                std::vector<int32_t> mis_vec(align.mismatches.begin(), align.mismatches.end());

                v_insertions[gc].push_back(std::move(ins_vec));
                v_deletions[gc].push_back(std::move(del_vec));
                v_mismatches[gc].push_back(std::move(mis_vec));
            } else {
                // No alignment for this gene class - add empty/default values
                v_gene_names[gc].push_back("");
                v_offsets[gc].push_back(0);
                v_five_p_offsets[gc].push_back(0);
                v_three_p_offsets[gc].push_back(0);
                v_align_lengths[gc].push_back(0);
                v_scores[gc].push_back(0.0);
                v_insertions[gc].push_back(std::vector<int32_t>{});
                v_deletions[gc].push_back(std::vector<int32_t>{});
                v_mismatches[gc].push_back(std::vector<int32_t>{});
            }
        }
    }

    // Create column names and arrays
    std::vector<std::string> column_names = { "sequence_id", "sequence" };
    std::vector<sparrow::array> arrays;

    // Add basic sequence data using public API
    sparrow::array id_array(sparrow::build(std::move(ids)));
    id_array.set_name("sequence_id");
    arrays.push_back(std::move(id_array));

    sparrow::array seq_array(sparrow::build(std::move(seqs)));
    seq_array.set_name("sequence");
    arrays.push_back(std::move(seq_array));

    // Add alignment columns for each gene class
    for (auto gc : { V_gene, D_gene, J_gene }) {
        std::string prefix;
        switch (gc) {
        case V_gene:
            prefix = "v_gene";
            break;
        case D_gene:
            prefix = "d_gene";
            break;
        case J_gene:
            prefix = "j_gene";
            break;
        default:
            continue;
        }

        // Add columns for this gene class
        std::string name_col = prefix + "_name";
        column_names.push_back(name_col);
        sparrow::array arr_name(sparrow::build(std::move(v_gene_names[gc])));
        arr_name.set_name(name_col);
        arrays.push_back(std::move(arr_name));

        std::string offset_col = prefix + "_offset";
        column_names.push_back(offset_col);
        sparrow::array arr_offset(sparrow::build(std::move(v_offsets[gc])));
        arr_offset.set_name(offset_col);
        arrays.push_back(std::move(arr_offset));

        std::string five_p_offset_col = prefix + "_five_p_offset";
        column_names.push_back(five_p_offset_col);
        sparrow::array arr_five_p(sparrow::build(std::move(v_five_p_offsets[gc])));
        arr_five_p.set_name(five_p_offset_col);
        arrays.push_back(std::move(arr_five_p));

        std::string three_p_offset_col = prefix + "_three_p_offset";
        column_names.push_back(three_p_offset_col);
        sparrow::array arr_three_p(sparrow::build(std::move(v_three_p_offsets[gc])));
        arr_three_p.set_name(three_p_offset_col);
        arrays.push_back(std::move(arr_three_p));

        std::string align_length_col = prefix + "_align_length";
        column_names.push_back(align_length_col);
        sparrow::array arr_align_len(sparrow::build(std::move(v_align_lengths[gc])));
        arr_align_len.set_name(align_length_col);
        arrays.push_back(std::move(arr_align_len));

        std::string score_col = prefix + "_score";
        column_names.push_back(score_col);
        sparrow::array arr_score(sparrow::build(std::move(v_scores[gc])));
        arr_score.set_name(score_col);
        arrays.push_back(std::move(arr_score));

        // Native list arrays with proper child naming
        // Note: We use targeted detail::array_access only to name list children,
        // as this is required by Arrow C interface but not exposed in public API
        std::string insertions_col = prefix + "_insertions";
        column_names.push_back(insertions_col);
        sparrow::array arr_ins(sparrow::build(std::move(v_insertions[gc])));
        arr_ins.set_name(insertions_col);
        // Name the list child array "item" (Arrow standard)
        auto& ins_proxy = sparrow::detail::array_access::get_arrow_proxy(arr_ins);
        if (!ins_proxy.children().empty()) {
            ins_proxy.children()[0].set_name("item");
        }
        arrays.push_back(std::move(arr_ins));

        std::string deletions_col = prefix + "_deletions";
        column_names.push_back(deletions_col);
        sparrow::array arr_del(sparrow::build(std::move(v_deletions[gc])));
        arr_del.set_name(deletions_col);
        auto& del_proxy = sparrow::detail::array_access::get_arrow_proxy(arr_del);
        if (!del_proxy.children().empty()) {
            del_proxy.children()[0].set_name("item");
        }
        arrays.push_back(std::move(arr_del));

        std::string mismatches_col = prefix + "_mismatches";
        column_names.push_back(mismatches_col);
        sparrow::array arr_mis(sparrow::build(std::move(v_mismatches[gc])));
        arr_mis.set_name(mismatches_col);
        auto& mis_proxy = sparrow::detail::array_access::get_arrow_proxy(arr_mis);
        if (!mis_proxy.children().empty()) {
            mis_proxy.children()[0].set_name("item");
        }
        arrays.push_back(std::move(arr_mis));
    }

    return sparrow::record_batch(std::move(column_names), std::move(arrays));
}

} // namespace igor
