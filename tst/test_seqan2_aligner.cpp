#include <catch2/catch_test_macros.hpp>

#include <igor/Core/Aligner.h>
#include <igor/Core/SeqAn2Aligner.h>

#ifdef IGOR_WITH_SEQAN2
#include <seqan/align.h>
#include <sstream>

static std::string cigar_visual(const std::string& cigar, const std::string& read, const std::string& germline)
{
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
        if (seqan2::isGap(read_row, i)) r << '-';
        else r << read_row[i];
        if (seqan2::isGap(germ_row, i)) g << '-';
        else g << germ_row[i];
        const char rc = r.str()[0];
        const char gc = g.str()[0];
        read_gapped += rc;
        germ_gapped += gc;
        match_line += (rc == gc && rc != '-') ? '|' : (rc != '-' && gc != '-') ? '*' : ' ';
    }
    std::ostringstream out;
    const size_t width = 100;
    for (size_t i = 0; i < read_gapped.size(); i += width) {
        out << "\n      read " << read_gapped.substr(i, width)
            << "\n           " << match_line.substr(i, width)
            << "\n      germ " << germ_gapped.substr(i, width);
    }
    return out.str();
}

static Matrix<double> simple_matrix()
{
    std::vector<double> v(15 * 15, -1.0);
    for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
    return Matrix<double>(15, 15, v);
}

TEST_CASE("SeqAn2Aligner basic global alignment", "[seqan2]")
{
    SeqAn2Aligner aligner(simple_matrix(), 2, V_gene, AlignmentPreset::with_band(20));
    aligner.set_genomic_sequences({{"g1", "ACGTACGT"}});
    auto alns = aligner.align_seq("ACGTACGT", 1.0, true, true);
    REQUIRE(!alns.empty());
    auto aln = alns.front();
    REQUIRE(aln.gene_name == "g1");
    REQUIRE(aln.five_p_offset == 0);
    REQUIRE(aln.three_p_offset == 7);
    REQUIRE(alignment_data_to_cigar(aln) == "8=");
}

TEST_CASE("SeqAn2Aligner allow_in_dels filters gaps", "[seqan2]")
{
    SeqAn2Aligner aligner(simple_matrix(), 2, V_gene, AlignmentPreset::with_band(20));
    aligner.set_genomic_sequences({{"g1", "ACGTACGT"}});
    auto with_gaps = aligner.align_seq("ACGTTACGT", -100.0, true, false);
    auto no_gaps = aligner.align_seq("ACGTTACGT", -100.0, false, false);
    REQUIRE(!with_gaps.empty());
    REQUIRE(no_gaps.empty());
}

TEST_CASE("known diagonal calibration updates band", "[seqan2][calibration]")
{
    auto preset = AlignmentPreset::with_band(1);
    BandCalibrationParams params;
    params.min_band_slack = 5;
    auto calibrated = calibrate_preset(preset, std::vector<int>{60, 70, 80}, params);
    REQUIRE(calibrated.band_lower_diag <= 55);
    REQUIRE(calibrated.band_upper_diag >= 85);
}

static const std::string READ_136_SEQUENCE = "CTATGTACCTCTGTGCCAGCCCGCGCGCTAAGGGGACAGGGAGGTCTGAAATCAGCCCCA";

TEST_CASE("SeqAn2Aligner read 136 alignment matches legacy aligner", "[seqan2][read136]")
{
    std::string germline_path = std::string(IGOR_SOURCE_DIR) + "/demo/genomicVs_with_primers.fasta";
    auto germlines = read_genomic_fasta(germline_path);
    REQUIRE(!germlines.empty());

    Matrix<double> subst;
    try {
        subst = read_substitution_matrix(std::string(IGOR_SOURCE_DIR) + "/models/heavy_pen_substitution_matrix.csv");
    }
    catch (...) {
        std::vector<double> v(15 * 15, -1.0);
        for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
        subst = Matrix<double>(15, 15, v);
    }

    std::cerr << "DP matrix "<< std::endl;
    std::cerr << subst << std::endl;
    
    
    const int band_half_width = 400;
    SeqAn2AlignConfig config;
    config.band_lower_diag = -band_half_width;
    config.band_upper_diag = band_half_width;
    config.seq2_leading_free = false;
    config.seq1_leading_free = true;
    config.seq1_trailing_free = true;
    config.seq2_trailing_free = true;
    //config.local_alignment = true;
    SeqAn2Aligner seqan_aligner(subst, 30, V_gene, config);
    seqan_aligner.set_genomic_sequences(germlines);

    
    Aligner legacy_aligner(subst, 30, V_gene);
    legacy_aligner.set_genomic_sequences(germlines);

    std::vector<std::pair<const int, const std::string>> reads = {{136, READ_136_SEQUENCE}};
    const double threshold = -1000000.0;

    auto seqan_results = seqan_aligner.align_seqs(reads, threshold, true);
    auto legacy_results = legacy_aligner.align_seqs(reads, threshold, true);

    REQUIRE(seqan_results.find(136) != seqan_results.end());
    REQUIRE(legacy_results.find(136) != legacy_results.end());

    const auto& seqan_alns = seqan_results[136];
    const auto& legacy_alns = legacy_results[136];

    REQUIRE(!seqan_alns.empty());
    REQUIRE(!legacy_alns.empty());

    
    for (const auto& aln : seqan_alns) {
        std::string cigar = alignment_data_to_cigar(aln);
        REQUIRE(!cigar.empty());
        REQUIRE(aln.score > -1e6);
        REQUIRE(!aln.gene_name.empty());
        
        REQUIRE(cigar.find_first_not_of("0123456789=XMIDSN") == std::string::npos);
    }

    
    auto seqan_best = std::max_element(
        seqan_alns.begin(), seqan_alns.end(),
        [](const Alignment_data& a, const Alignment_data& b) { return a.score < b.score; });

    auto legacy_best = std::max_element(
        legacy_alns.begin(), legacy_alns.end(),
        [](const Alignment_data& a, const Alignment_data& b) { return a.score < b.score; });

    const auto seqan_cigar = alignment_data_to_cigar(*seqan_best);

    const size_t read_len = READ_136_SEQUENCE.size();
    size_t germline_len = 0;
    for (const auto& gl : germlines) {
        if (gl.first == legacy_best->gene_name) {
            germline_len = gl.second.size();
            break;
        }
    }
    const auto legacy_cigar = alignment_data_to_cigar_full_span(*legacy_best, read_len, germline_len);

    std::string seqan_germline;
    for (const auto& gl : germlines) {
        if (gl.first == seqan_best->gene_name) {
            seqan_germline = gl.second;
            break;
        }
    }
    std::string legacy_germline;
    for (const auto& gl : germlines) {
        if (gl.first == legacy_best->gene_name) {
            legacy_germline = gl.second;
            break;
        }
    }
    std::cerr << "\n[Read 136 Comparison]\n"
              << "  SeqAn2: gene=" << seqan_best->gene_name
              << " score=" << seqan_best->score
              << " CIGAR=" << seqan_cigar
              << cigar_visual(seqan_cigar, READ_136_SEQUENCE, seqan_germline) << "\n"
              << "  Legacy: gene=" << legacy_best->gene_name
              << " score=" << legacy_best->score
              << " CIGAR=" << legacy_cigar
              << cigar_visual(legacy_cigar, READ_136_SEQUENCE, legacy_germline) << "\n" << std::flush;

    
    REQUIRE(seqan_best->gene_name == legacy_best->gene_name);

    
    REQUIRE(seqan_best->score == legacy_best->score);
    REQUIRE(seqan_cigar == legacy_cigar);
}
#endif