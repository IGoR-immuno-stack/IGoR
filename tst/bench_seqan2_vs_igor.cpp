#include <igor/Core/Aligner.h>
#include <igor/Core/AlignmentPreset.h>
#ifdef IGOR_WITH_SEQAN2
#include <igor/Core/SeqAn2Aligner.h>
#endif

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <forward_list>
#include <map>

static Matrix<double> simple_matrix()
{
    std::vector<double> v(15 * 15, -1.0);
    for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
    return Matrix<double>(15, 15, v);
}

static std::vector<std::pair<const int, const std::string>> load_reads(const char* path)
{
    try { return read_txt(path); }
    catch (...) { return {{0, "ACGTACGTACGT"}, {1, "TGCATGCATGCA"}}; }
}

int main(int argc, char** argv)
{
    std::string seq_path = argc > 1 ? argv[1] : std::string(IGOR_SOURCE_DIR) + "/demo/murugan_naive1_noncoding_demo_seqs.txt";
    int repeat = argc > 2 ? std::max(1, std::atoi(argv[2])) : 100;
    auto base_reads = load_reads(seq_path.c_str());
    std::vector<std::pair<const int, const std::string>> reads_const;
    reads_const.reserve(base_reads.size() * repeat);
    for (int r = 0; r < repeat; ++r) {
        for (const auto& read : base_reads) {
            reads_const.emplace_back(static_cast<int>(reads_const.size()), read.second);
        }
    }
    std::vector<std::pair<std::string, std::string>> germlines{{"synthetic", "ACGTACGTACGT"}};
    Matrix<double> subst = simple_matrix();

    Aligner legacy(subst, 2, V_gene);
    legacy.set_genomic_sequences(germlines);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto legacy_res = legacy.align_seqs(reads_const, -100.0, true);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "{\n";
    std::cout << "  \"input_sequences\": " << base_reads.size() << ",\n";
    std::cout << "  \"repeat\": " << repeat << ",\n";
    std::cout << "  \"benchmark_sequences\": " << reads_const.size() << ",\n";
    auto best_score = [](const auto& xs) {
        double best = -1e300;
        for (const auto& kv : xs) for (const auto& a : kv.second) best = std::max(best, a.score);
        return best;
    };
    std::cout << "  \"legacy_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << ",\n";
    std::cout << "  \"legacy_sequences\": " << legacy_res.size() << ",\n";
    std::cout << "  \"legacy_best_score\": " << best_score(legacy_res);

#ifdef IGOR_WITH_SEQAN2
    std::vector<std::pair<int, std::string>> reads;
    for (const auto& r : reads_const) reads.push_back({r.first, r.second});

    const int band_half_width = 50;
    SeqAn2Aligner seqan2_banded(subst, 2, V_gene, AlignmentPreset::with_band(band_half_width));
    seqan2_banded.set_genomic_sequences(germlines);
    auto b0 = std::chrono::high_resolution_clock::now();
    auto seqan_banded_res = seqan2_banded.align_seqs(reads, -100.0, true);
    auto b1 = std::chrono::high_resolution_clock::now();

    SeqAn2AlignConfig unbanded_config;
    unbanded_config.band_lower_diag = 1;
    unbanded_config.band_upper_diag = 0; // lower > upper means full-matrix unbanded alignment.
    unbanded_config.seq1_trailing_free = false;
    unbanded_config.seq2_trailing_free = false;
    SeqAn2Aligner seqan2_unbanded(subst, 2, V_gene, unbanded_config);
    seqan2_unbanded.set_genomic_sequences(germlines);
    auto u0 = std::chrono::high_resolution_clock::now();
    auto seqan_unbanded_res = seqan2_unbanded.align_seqs(reads, -100.0, true);
    auto u1 = std::chrono::high_resolution_clock::now();

    std::cout << ",\n  \"seqan2_banded_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(b1 - b0).count() << ",\n";
    std::cout << "  \"seqan2_banded_sequences\": " << seqan_banded_res.size() << ",\n";
    std::cout << "  \"seqan2_banded_band_lower_diag\": " << -band_half_width << ",\n";
    std::cout << "  \"seqan2_banded_band_upper_diag\": " << band_half_width << ",\n";
    std::cout << "  \"seqan2_banded_best_score\": " << best_score(seqan_banded_res) << ",\n";
    std::cout << "  \"seqan2_unbanded_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(u1 - u0).count() << ",\n";
    std::cout << "  \"seqan2_unbanded_sequences\": " << seqan_unbanded_res.size() << ",\n";
    std::cout << "  \"seqan2_unbanded_best_score\": " << best_score(seqan_unbanded_res) << ",\n";
    std::cout << "  \"score_delta_seqan2_unbanded_minus_legacy\": "
              << (best_score(seqan_unbanded_res) - best_score(legacy_res));
#endif
    std::cout << "\n}\n";
    return 0;
}
