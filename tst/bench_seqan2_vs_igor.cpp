#include <igor/Core/Aligner.h>
#include <igor/Core/AlignmentPreset.h>
#ifdef IGOR_WITH_SEQAN2
#include <igor/Core/SeqAn2Aligner.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <forward_list>
#include <map>

static Matrix<double> load_matrix(const std::string& path)
{
    try { return read_substitution_matrix(path); }
    catch (...) {
        std::vector<double> v(15 * 15, -1.0);
        for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
        return Matrix<double>(15, 15, v);
    }
}

static std::vector<std::pair<const int, const std::string>> load_reads(const char* path)
{
    try { return read_txt(path); }
    catch (...) { return {{0, "ACGTACGTACGT"}, {1, "TGCATGCATGCA"}}; }
}

int main(int argc, char** argv)
{
    std::string seq_path = argc > 1 ? argv[1] : std::string(IGOR_SOURCE_DIR) + "/demo/murugan_naive1_noncoding_demo_seqs.txt";
    int repeat = argc > 2 ? (std::max)(1, std::atoi(argv[2])) : 100;
    std::string germline_path = argc > 3 ? argv[3] : std::string(IGOR_SOURCE_DIR) + "/demo/genomicVs_with_primers.fasta";
    std::string matrix_path = argc > 4 ? argv[4] : std::string(IGOR_SOURCE_DIR) + "/models/heavy_pen_substitution_matrix.csv";
    auto base_reads = load_reads(seq_path.c_str());
    std::vector<std::pair<const int, const std::string>> reads_const;
    reads_const.reserve(base_reads.size() * repeat);
    for (int r = 0; r < repeat; ++r) {
        for (const auto& read : base_reads) {
            reads_const.emplace_back(static_cast<int>(reads_const.size()), read.second);
        }
    }
    auto germlines = read_genomic_fasta(germline_path);
    Matrix<double> subst = load_matrix(matrix_path);

    auto avg_len_reads = [&]() {
        size_t total = 0;
        for (const auto& r : base_reads) total += r.second.size();
        return base_reads.empty() ? 0.0 : static_cast<double>(total) / base_reads.size();
    };
    auto avg_len_germlines = [&]() {
        size_t total = 0;
        for (const auto& g : germlines) total += g.second.size();
        return germlines.empty() ? 0.0 : static_cast<double>(total) / germlines.size();
    };

    Aligner legacy(subst, 2, V_gene);
    legacy.set_genomic_sequences(germlines);
    auto t0 = std::chrono::high_resolution_clock::now();
    const double threshold = -1000000.0;
    auto legacy_res = legacy.align_seqs(reads_const, threshold, true);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "{\n";
    std::cout << "  \"input_sequences\": " << base_reads.size() << ",\n";
    std::cout << "  \"repeat\": " << repeat << ",\n";
    std::cout << "  \"benchmark_sequences\": " << reads_const.size() << ",\n";
    std::cout << "  \"germline_templates\": " << germlines.size() << ",\n";
    std::cout << "  \"avg_read_length\": " << avg_len_reads() << ",\n";
    std::cout << "  \"avg_germline_length\": " << avg_len_germlines() << ",\n";
    std::cout << "  \"germline_fasta\": \"" << germline_path << "\",\n";
    auto best_score = [](const auto& xs) {
        double best = -1e300;
        for (const auto& kv : xs) for (const auto& a : kv.second) best = (std::max)(best, a.score);
        return best;
    };
    auto alignment_count = [](const auto& xs) {
        size_t count = 0;
        for (const auto& kv : xs) for (const auto& a : kv.second) { (void)a; ++count; }
        return count;
    };
    struct ScoreComparison {
        size_t common = 0;
        size_t missing_in_candidate = 0;
        size_t missing_in_legacy = 0;
        double sum_abs_delta = 0.0;
        double max_abs_delta = 0.0;
    };
    auto list_best_score = [](const auto& alignments, double& score) {
        bool found = false;
        for (const auto& a : alignments) {
            if (!found || a.score > score) score = a.score;
            found = true;
        }
        return found;
    };
    auto compare_best_scores = [&](const auto& legacy_alignments, const auto& candidate_alignments) {
        ScoreComparison comparison;
        for (const auto& kv : legacy_alignments) {
            double legacy_score = 0.0;
            if (!list_best_score(kv.second, legacy_score)) continue;
            auto candidate_iter = candidate_alignments.find(kv.first);
            if (candidate_iter == candidate_alignments.end()) {
                ++comparison.missing_in_candidate;
                continue;
            }
            double candidate_score = 0.0;
            if (!list_best_score(candidate_iter->second, candidate_score)) {
                ++comparison.missing_in_candidate;
                continue;
            }
            const double abs_delta = std::abs(candidate_score - legacy_score);
            ++comparison.common;
            comparison.sum_abs_delta += abs_delta;
            comparison.max_abs_delta = (std::max)(comparison.max_abs_delta, abs_delta);
        }
        for (const auto& kv : candidate_alignments) {
            double candidate_score = 0.0;
            if (!list_best_score(kv.second, candidate_score)) continue;
            auto legacy_iter = legacy_alignments.find(kv.first);
            if (legacy_iter == legacy_alignments.end()) {
                ++comparison.missing_in_legacy;
                continue;
            }
            double legacy_score = 0.0;
            if (!list_best_score(legacy_iter->second, legacy_score)) ++comparison.missing_in_legacy;
        }
        return comparison;
    };
    auto print_score_comparison = [](const char* prefix, const ScoreComparison& comparison) {
        std::cout << "  \"" << prefix << "_score_common_sequences\": " << comparison.common << ",\n";
        std::cout << "  \"" << prefix << "_score_missing_in_seqan2\": " << comparison.missing_in_candidate << ",\n";
        std::cout << "  \"" << prefix << "_score_missing_in_legacy\": " << comparison.missing_in_legacy << ",\n";
        std::cout << "  \"" << prefix << "_score_mean_abs_delta\": "
                  << (comparison.common ? comparison.sum_abs_delta / comparison.common : 0.0) << ",\n";
        std::cout << "  \"" << prefix << "_score_max_abs_delta\": " << comparison.max_abs_delta;
    };
    std::cout << "  \"legacy_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << ",\n";
    std::cout << "  \"legacy_sequences\": " << legacy_res.size() << ",\n";
    std::cout << "  \"legacy_alignments\": " << alignment_count(legacy_res) << ",\n";
    std::cout << "  \"legacy_best_score\": " << best_score(legacy_res);

#ifdef IGOR_WITH_SEQAN2
    std::vector<std::pair<int, std::string>> reads;
    for (const auto& r : reads_const) reads.push_back({r.first, r.second});

    const int band_half_width = 400;
    SeqAn2AlignConfig legacy_like_banded_config;
    legacy_like_banded_config.band_lower_diag = -band_half_width;
    legacy_like_banded_config.band_upper_diag = band_half_width;
    legacy_like_banded_config.seq2_leading_free = true;
    legacy_like_banded_config.seq2_trailing_free = true;
    legacy_like_banded_config.seq1_leading_free = true;
    legacy_like_banded_config.seq1_trailing_free = true;
    SeqAn2Aligner seqan2_banded(subst, 2, V_gene, legacy_like_banded_config);
    seqan2_banded.set_genomic_sequences(germlines);
    auto b0 = std::chrono::high_resolution_clock::now();
    auto seqan_banded_res = seqan2_banded.align_seqs(reads, threshold, true);
    auto b1 = std::chrono::high_resolution_clock::now();

    SeqAn2AlignConfig unbanded_config = legacy_like_banded_config;
    unbanded_config.band_lower_diag = 1;
    unbanded_config.band_upper_diag = 0; // lower > upper means full-matrix unbanded alignment.
    SeqAn2Aligner seqan2_unbanded(subst, 2, V_gene, unbanded_config);
    seqan2_unbanded.set_genomic_sequences(germlines);
    auto u0 = std::chrono::high_resolution_clock::now();
    auto seqan_unbanded_res = seqan2_unbanded.align_seqs(reads, threshold, true);
    auto u1 = std::chrono::high_resolution_clock::now();

    std::cout << ",\n  \"seqan2_banded_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(b1 - b0).count() << ",\n";
    std::cout << "  \"seqan2_banded_sequences\": " << seqan_banded_res.size() << ",\n";
    std::cout << "  \"seqan2_banded_alignments\": " << alignment_count(seqan_banded_res) << ",\n";
    std::cout << "  \"seqan2_banded_band_lower_diag\": " << -band_half_width << ",\n";
    std::cout << "  \"seqan2_banded_band_upper_diag\": " << band_half_width << ",\n";
    std::cout << "  \"seqan2_v_gene_legacy_like_free_all_ends\": true,\n";
    std::cout << "  \"seqan2_banded_best_score\": " << best_score(seqan_banded_res) << ",\n";
    print_score_comparison("seqan2_banded_vs_legacy", compare_best_scores(legacy_res, seqan_banded_res));
    std::cout << ",\n  \"seqan2_unbanded_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(u1 - u0).count() << ",\n";
    std::cout << "  \"seqan2_unbanded_sequences\": " << seqan_unbanded_res.size() << ",\n";
    std::cout << "  \"seqan2_unbanded_alignments\": " << alignment_count(seqan_unbanded_res) << ",\n";
    std::cout << "  \"seqan2_unbanded_best_score\": " << best_score(seqan_unbanded_res) << ",\n";
    print_score_comparison("seqan2_unbanded_vs_legacy", compare_best_scores(legacy_res, seqan_unbanded_res));
    std::cout << ",\n  \"score_delta_seqan2_unbanded_minus_legacy\": "
              << (best_score(seqan_unbanded_res) - best_score(legacy_res));
#endif
    std::cout << "\n}\n";
    return 0;
}
