/*
 * StreamingTestUtils.cpp
 *
 *  Created on: Feb 3, 2026
 *      Author: IGoR Development Team
 *
 *  Implementation of shared test utilities for Streaming module tests
 */

#include "StreamingTestUtils.h"

namespace igor::test {

std::vector<SequenceTuple> load_murugan_dataset()
{
    std::string test_data_dir = std::string(IGOR_TEST_DATA_DIR) + "/";

    // Parse indexed sequences
    std::unordered_map<int, std::string> sequences;
    std::ifstream seq_file(test_data_dir + "murugan_naive1_noncoding_demo_seqs_indexed_seq.csv");
    std::string line;
    std::getline(seq_file, line); // Skip header

    while (std::getline(seq_file, line)) {
        size_t pos = line.find(';');
        if (pos != std::string::npos) {
            int seq_index = std::stoi(line.substr(0, pos));
            std::string sequence = line.substr(pos + 1);
            sequences[seq_index] = sequence;
        }
    }

    // Parse alignments for each gene class
    std::unordered_map<int, std::unordered_map<Gene_class, std::vector<Alignment_data>>> all_alignments;

    auto parse_alignment_file = [&](const std::string& filename, Gene_class gene_class) {
        std::ifstream file(test_data_dir + filename);
        std::string header;
        std::getline(file, header); // Skip header

        while (std::getline(file, line)) {
            std::vector<std::string> fields;
            size_t start = 0;
            size_t end;

            // Split by semicolon
            while ((end = line.find(';', start)) != std::string::npos) {
                fields.push_back(line.substr(start, end - start));
                start = end + 1;
            }
            fields.push_back(line.substr(start));

            if (fields.size() < 10) continue;

            int seq_index = std::stoi(fields[0]);
            std::string gene_name = fields[1];
            double score = std::stod(fields[2]);
            int offset = std::stoi(fields[3]);

            // Parse insertions list {a,b,c}
            auto parse_list = [](const std::string& str) -> std::vector<int> {
                std::vector<int> result;
                if (str.size() > 2 && str[0] == '{' && str.back() == '}') {
                    std::string content = str.substr(1, str.size() - 2);
                    if (!content.empty()) {
                        size_t pos = 0;
                        while (pos < content.size()) {
                            size_t comma = content.find(',', pos);
                            if (comma == std::string::npos) {
                                result.push_back(std::stoi(content.substr(pos)));
                                break;
                            }
                            result.push_back(std::stoi(content.substr(pos, comma - pos)));
                            pos = comma + 1;
                        }
                    }
                }
                return result;
            };

            auto insertions_vec = parse_list(fields[4]);
            auto deletions_vec = parse_list(fields[5]);
            auto mismatches_vec = parse_list(fields[6]);

            std::forward_list<int> insertions(insertions_vec.begin(), insertions_vec.end());
            std::forward_list<int> deletions(deletions_vec.begin(), deletions_vec.end());
            std::vector<int> mismatches(mismatches_vec.begin(), mismatches_vec.end());

            size_t align_length = std::stoull(fields[7]);
            size_t five_p_offset = std::stoull(fields[8]);
            size_t three_p_offset = std::stoull(fields[9]);

            // Use 9-argument constructor
            Alignment_data align(gene_name, offset, five_p_offset, three_p_offset,
                               align_length, insertions, deletions, mismatches, score);

            all_alignments[seq_index][gene_class].push_back(align);
        }
    };

    parse_alignment_file("murugan_naive1_noncoding_demo_seqs_alignments_V.csv", V_gene);
    parse_alignment_file("murugan_naive1_noncoding_demo_seqs_alignments_D.csv", D_gene);
    parse_alignment_file("murugan_naive1_noncoding_demo_seqs_alignments_J.csv", J_gene);

    // Combine into result vector
    std::vector<SequenceTuple> result;
    for (const auto& [seq_index, sequence] : sequences) {
        result.emplace_back(seq_index, sequence, all_alignments[seq_index]);
    }

    return result;
}

} // namespace igor::test
