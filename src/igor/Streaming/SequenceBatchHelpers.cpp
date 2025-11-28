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
                          std::is_same_v<T, sparrow::nullable<const int32_t&>> ||
                          std::is_same_v<T, sparrow::nullable<const long long&>>) {
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
                          std::is_same_v<T, sparrow::nullable<const float&>>) {
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

std::unordered_map<Gene_class, std::vector<Alignment_data>>
parse_alignments_from_columns(const sparrow::record_batch &batch, size_t row_index)
{

    std::unordered_map<Gene_class, std::vector<Alignment_data>> alignments;

    // For now, return empty alignments
    // TODO: Implement full alignment parsing when format is defined
    // This will be expanded in later tasks to handle actual alignment columns

    // Look for common alignment column patterns:
    // - v_gene_name, v_gene_offset, v_gene_score
    // - d_gene_name, d_gene_offset, d_gene_score
    // - j_gene_name, j_gene_offset, j_gene_score

    const std::vector<std::pair<std::string, Gene_class>> gene_prefixes = { { "v_gene", V_gene },
                                                                            { "d_gene", D_gene },
                                                                            { "j_gene", J_gene } };

    for (const auto &[prefix, gene_class] : gene_prefixes) {
        std::string name_col = prefix + "_name";
        std::string offset_col = prefix + "_offset";
        std::string score_col = prefix + "_score";

        if (has_column(batch, name_col)) {
            std::string gene_name = get_string_value(batch, name_col, row_index);

            if (!gene_name.empty()) {
                int offset = get_int_value(batch, offset_col, row_index, 0);
                double score = get_double_value(batch, score_col, row_index, 0.0);

                // Create a basic alignment data structure
                Alignment_data align(gene_name, offset);
                align.score = score;
                align.align_length = 0; // Will be populated when actual alignment info is available

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

    for (const auto &seq_tuple : sequences) {
        ids.push_back(static_cast<int32_t>(std::get<0>(seq_tuple)));
        seqs.push_back(std::get<1>(seq_tuple));
    }

    // Create Sparrow arrays using the builder
    std::vector<std::string> column_names = { "sequence_id", "sequence" };
    std::vector<sparrow::array> arrays;

    // Use sparrow::build and emplace_back for implicit conversion to sparrow::array
    arrays.emplace_back(sparrow::build(std::move(ids)));
    arrays.emplace_back(sparrow::build(std::move(seqs)));

    // TODO: Add alignment columns in future tasks
    // For now, we only convert the basic sequence data

    return sparrow::record_batch(std::move(column_names), std::move(arrays));
}

} // namespace igor
