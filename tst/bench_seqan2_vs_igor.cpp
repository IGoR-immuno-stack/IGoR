#include <igor/Core/Aligner.h>
#include <igor/Core/AlignmentPreset.h>
#ifdef IGOR_WITH_SEQAN2
#include <igor/Core/SeqAn2Aligner.h>
#include <seqan/align.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <forward_list>
#include <map>
#include <sstream>
#include <unordered_map>

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
    std::unordered_map<int, size_t> read_lengths;
    std::unordered_map<int, std::string> read_sequences;
    for (const auto& read : reads_const) {
        read_lengths[read.first] = read.second.size();
        read_sequences[read.first] = read.second;
    }
    std::unordered_map<std::string, size_t> germline_lengths;
    std::unordered_map<std::string, std::string> germline_sequences;
    for (const auto& germline : germlines) {
        germline_lengths[germline.first] = germline.second.size();
        germline_sequences[germline.first] = germline.second;
    }
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
    struct BestCigar {
        bool found = false;
        std::string gene_name;
        std::string cigar;
        double score = 0.0;
    };
    struct CigarComparison {
        size_t common = 0;
        size_t equal = 0;
        size_t mismatch = 0;
        size_t missing_in_candidate = 0;
        size_t missing_in_legacy = 0;
        size_t conversion_errors = 0;
        std::vector<std::string> mismatch_samples;
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
    auto best_cigar = [&](int read_id, const auto& alignments, BestCigar& best) {
        for (const auto& a : alignments) {
            std::string cigar;
            try {
                auto read_len = read_lengths.find(read_id);
                auto germline_len = germline_lengths.find(a.gene_name);
                if (read_len == read_lengths.end() || germline_len == germline_lengths.end()) return false;
                cigar = alignment_data_to_cigar_full_span(a, read_len->second, germline_len->second);
            }
            catch (...) { return false; }
            if (!best.found || a.score > best.score ||
                (a.score == best.score && std::make_pair(a.gene_name, cigar) < std::make_pair(best.gene_name, best.cigar))) {
                best.found = true;
                best.gene_name = a.gene_name;
                best.cigar = cigar;
                best.score = a.score;
            }
        }
        return best.found;
    };
    auto short_gene_name = [](const std::string& gene_name) {
        const size_t pipe = gene_name.find('|');
        const std::string compact = pipe == std::string::npos ? gene_name : gene_name.substr(0, pipe);
        return compact.size() <= 32 ? compact : compact.substr(0, 29) + "...";
    };
    auto cigar_visual = [&](const std::string& cigar, const std::string& read, const std::string& germline) {
#ifdef IGOR_WITH_SEQAN2
        using TSeq = seqan2::String<seqan2::Iupac>;
        using TAlign = seqan2::Align<TSeq, seqan2::ArrayGaps>;
        TSeq read_seq;
        TSeq germ_seq;
        seqan2::assign(read_seq, read);
        seqan2::assign(germ_seq, germline);
        TAlign ali;
        seqan2::resize(seqan2::rows(ali), 2);
        seqan2::assignSource(seqan2::row(ali, 0), read_seq);
        seqan2::assignSource(seqan2::row(ali, 1), germ_seq);
        auto& read_row = seqan2::row(ali, 0);
        auto& germ_row = seqan2::row(ali, 1);
        size_t col = 0;
        for (const auto& entry : parse_cigar(cigar)) {
            switch (entry.second) {
            case 'I':
                seqan2::insertGaps(germ_row, col, entry.first);
                col += entry.first;
                break;
            case 'D':
            case 'N':
                seqan2::insertGaps(read_row, col, entry.first);
                col += entry.first;
                break;
            case '=':
            case 'X':
            case 'M':
            case 'S':
                col += entry.first;
                break;
            default:
                break;
            }
        }
        std::string read_gapped;
        std::string match_line;
        std::string germ_gapped;
        const size_t len = seqan2::length(read_row);
        for (size_t i = 0; i < len; ++i) {
            std::ostringstream r;
            std::ostringstream g;
            if (seqan2::isGap(read_row, i)) r << '-'; else r << read_row[i];
            if (seqan2::isGap(germ_row, i)) g << '-'; else g << germ_row[i];
            const char rc = r.str()[0];
            const char gc = g.str()[0];
            read_gapped += rc;
            germ_gapped += gc;
            match_line += (rc == gc && rc != '-') ? '|' : (rc != '-' && gc != '-' ? '*' : ' ');
        }
        std::ostringstream out;
        const size_t width = 100;
        for (size_t i = 0; i < read_gapped.size(); i += width) {
            out << "\n      read " << read_gapped.substr(i, width)
                << "\n           " << match_line.substr(i, width)
                << "\n      germ " << germ_gapped.substr(i, width);
        }
        return out.str();
#else
        (void)cigar; (void)read; (void)germline;
        return std::string("\n      visual alignment unavailable; built without IGOR_WITH_SEQAN2");
#endif
    };
    auto compare_best_cigars = [&](const auto& legacy_alignments, const auto& candidate_alignments) {
        CigarComparison comparison;
        for (const auto& kv : legacy_alignments) {
            BestCigar legacy_cigar;
            if (!best_cigar(kv.first, kv.second, legacy_cigar)) {
                ++comparison.conversion_errors;
                continue;
            }
            auto candidate_iter = candidate_alignments.find(kv.first);
            if (candidate_iter == candidate_alignments.end()) {
                ++comparison.missing_in_candidate;
                continue;
            }
            BestCigar candidate_cigar;
            if (!best_cigar(kv.first, candidate_iter->second, candidate_cigar)) {
                ++comparison.conversion_errors;
                continue;
            }
            ++comparison.common;
            if (legacy_cigar.gene_name == candidate_cigar.gene_name && legacy_cigar.cigar == candidate_cigar.cigar) {
                ++comparison.equal;
            } else {
                ++comparison.mismatch;
                if (comparison.mismatch_samples.size() < 5) {
                    std::ostringstream sample;
                    const auto read_it = read_sequences.find(kv.first);
                    const auto legacy_germline = germline_sequences.find(legacy_cigar.gene_name);
                    const auto candidate_germline = germline_sequences.find(candidate_cigar.gene_name);
                    sample << "read " << kv.first
                           << ": legacy_gene=" << short_gene_name(legacy_cigar.gene_name)
                           << " score=" << legacy_cigar.score
                           << " cigar=" << legacy_cigar.cigar
                           << " | seqan2_gene=" << short_gene_name(candidate_cigar.gene_name)
                           << " score=" << candidate_cigar.score
                           << " cigar=" << candidate_cigar.cigar;
                    if (read_it != read_sequences.end() && legacy_germline != germline_sequences.end()) {
                        sample << "\n      legacy visual:"
                               << cigar_visual(legacy_cigar.cigar, read_it->second, legacy_germline->second);
                    }
                    if (read_it != read_sequences.end() && candidate_germline != germline_sequences.end()) {
                        sample << "\n      seqan2 visual:"
                               << cigar_visual(candidate_cigar.cigar, read_it->second, candidate_germline->second);
                    }
                    comparison.mismatch_samples.push_back(sample.str());
                }
            }
        }
        for (const auto& kv : candidate_alignments) {
            BestCigar candidate_cigar;
            if (!best_cigar(kv.first, kv.second, candidate_cigar)) {
                ++comparison.conversion_errors;
                continue;
            }
            auto legacy_iter = legacy_alignments.find(kv.first);
            if (legacy_iter == legacy_alignments.end()) {
                ++comparison.missing_in_legacy;
                continue;
            }
            BestCigar legacy_cigar;
            if (!best_cigar(kv.first, legacy_iter->second, legacy_cigar)) ++comparison.missing_in_legacy;
        }
        return comparison;
    };
    auto print_section = [](const std::string& title) {
        std::cout << "\n" << title << "\n" << std::string(title.size(), '-') << "\n";
    };
    auto print_metric = [](const std::string& label, const auto& value) {
        std::cout << "  " << std::left << std::setw(28) << label << " : " << value << "\n";
    };
    auto ms = [](auto start, auto end) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    };
    auto print_score_comparison = [&](const std::string& title, const ScoreComparison& comparison) {
        print_section(title + " scores");
        print_metric("common sequences", comparison.common);
        print_metric("missing in SeqAn2", comparison.missing_in_candidate);
        print_metric("missing in legacy", comparison.missing_in_legacy);
        print_metric("mean abs delta", comparison.common ? comparison.sum_abs_delta / comparison.common : 0.0);
        print_metric("max abs delta", comparison.max_abs_delta);
    };
    auto print_cigar_comparison = [&](const std::string& title, const CigarComparison& comparison) {
        print_section(title + " CIGARs");
        print_metric("common sequences", comparison.common);
        print_metric("equal", comparison.equal);
        print_metric("mismatched", comparison.mismatch);
        print_metric("missing in SeqAn2", comparison.missing_in_candidate);
        print_metric("missing in legacy", comparison.missing_in_legacy);
        print_metric("conversion errors", comparison.conversion_errors);
        if (!comparison.mismatch_samples.empty()) {
            std::cout << "  mismatch samples\n";
            for (const auto& sample : comparison.mismatch_samples) std::cout << "    - " << sample << "\n";
        }
    };

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "IGoR / SeqAn2 benchmark\n";
    print_section("Input");
    print_metric("input sequences", base_reads.size());
    print_metric("repeat", repeat);
    print_metric("benchmark sequences", reads_const.size());
    print_metric("germline templates", germlines.size());
    print_metric("avg read length", avg_len_reads());
    print_metric("avg germline length", avg_len_germlines());
    print_metric("germline fasta", germline_path);

    print_section("Legacy IGoR");
    print_metric("time ms", ms(t0, t1));
    print_metric("sequences with hits", legacy_res.size());
    print_metric("alignments", alignment_count(legacy_res));
    print_metric("best score", best_score(legacy_res));

#ifdef IGOR_WITH_SEQAN2
    std::vector<std::pair<const int, const std::string>> reads;
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

    print_section("SeqAn2 banded");
    print_metric("time ms", ms(b0, b1));
    print_metric("sequences with hits", seqan_banded_res.size());
    print_metric("alignments", alignment_count(seqan_banded_res));
    print_metric("band lower diag", -band_half_width);
    print_metric("band upper diag", band_half_width);
    print_metric("free all ends", "yes");
    print_metric("best score", best_score(seqan_banded_res));
    if (ms(b0, b1) > 0) print_metric("speedup vs legacy", static_cast<double>(ms(t0, t1)) / ms(b0, b1));
    print_score_comparison("SeqAn2 banded vs legacy", compare_best_scores(legacy_res, seqan_banded_res));
    print_cigar_comparison("SeqAn2 banded vs legacy", compare_best_cigars(legacy_res, seqan_banded_res));

    print_section("SeqAn2 unbanded");
    print_metric("time ms", ms(u0, u1));
    print_metric("sequences with hits", seqan_unbanded_res.size());
    print_metric("alignments", alignment_count(seqan_unbanded_res));
    print_metric("best score", best_score(seqan_unbanded_res));
    if (ms(u0, u1) > 0) print_metric("speedup vs legacy", static_cast<double>(ms(t0, t1)) / ms(u0, u1));
    print_metric("best score delta", best_score(seqan_unbanded_res) - best_score(legacy_res));
    print_score_comparison("SeqAn2 unbanded vs legacy", compare_best_scores(legacy_res, seqan_unbanded_res));
    print_cigar_comparison("SeqAn2 unbanded vs legacy", compare_best_cigars(legacy_res, seqan_unbanded_res));
#else
    print_section("SeqAn2");
    print_metric("status", "unavailable; built without IGOR_WITH_SEQAN2");
#endif
    return 0;
}
