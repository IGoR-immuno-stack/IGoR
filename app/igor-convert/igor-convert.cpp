/*
 * igor-convert.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  Command-line tool for converting between IGoR data formats:
 *  - Parquet (high-performance binary)
 *  - AIRR Rearrangement (TSV/CSV, one row per sequence)
 *  - AIRR Alignment (TSV/CSV, one row per alignment)
 */

#include <igor/Streaming/ParquetReader.h>
#include <igor/Streaming/ParquetWriter.h>
#include <igor/Streaming/AIRRRearrangementReader.h>
#include <igor/Streaming/AIRRRearrangementWriter.h>
#include <igor/Streaming/AIRRAlignmentReader.h>
#include <igor/Streaming/AIRRAlignmentWriter.h>
#include <igor/Streaming/SequenceBatchHelpers.h>

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <tuple>
#include <cstring>

namespace fs = std::filesystem;

// Bring igor types into scope for convenience
using namespace igor;

// Type alias for the legacy sequence tuple format
using SequenceTuple = std::tuple<
    int,
    std::string,
    std::unordered_map<Gene_class, std::vector<Alignment_data>>
>;



//==============================================================================
// Format detection
//==============================================================================

enum class Format
{
    Parquet,
    AIRRRearrangement,
    AIRRAlignment,
    IGoRLegacy,  // Legacy indexed_seq format
    Unknown
};

Format detect_format(const std::string& filepath, bool is_output = false)
{
    fs::path path(filepath);
    auto ext = path.extension().string();

    // Extension-based detection
    if (ext == ".parquet") {
        return Format::Parquet;
    }

    if (ext == ".tsv" || ext == ".csv") {
        // For output files that don't exist yet, try to infer from filename
        if (is_output && !fs::exists(filepath)) {
            std::string filename = path.filename().string();
            std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);

            // Check if filename suggests alignment format
            if (filename.find("align") != std::string::npos) {
                return Format::AIRRAlignment;
            }
            // Default to Rearrangement
            return Format::AIRRRearrangement;
        }

        // Check if file exists and peek at header
        if (!fs::exists(filepath)) {
            return Format::Unknown;
        }

        std::ifstream file(filepath);
        std::string header;
        if (std::getline(file, header)) {
            // Legacy IGoR format has 'seq_index' (semicolon-delimited)
            if (header.find("seq_index") != std::string::npos && header.find(";") != std::string::npos) {
                return Format::IGoRLegacy;
            }
            // AIRR Alignment has 'segment' column
            if (header.find("segment") != std::string::npos) {
                return Format::AIRRAlignment;
            }
            // AIRR Rearrangement has v_call, d_call, j_call
            if (header.find("v_call") != std::string::npos ||
                header.find("d_call") != std::string::npos ||
                header.find("j_call") != std::string::npos) {
                return Format::AIRRRearrangement;
            }
        }
    }

    return Format::Unknown;
}


const char* format_name(Format fmt)
{
    switch (fmt) {
        case Format::Parquet: return "Parquet";
        case Format::AIRRRearrangement: return "AIRR Rearrangement";
        case Format::AIRRAlignment: return "AIRR Alignment";
        case Format::IGoRLegacy: return "IGoR Legacy";
        case Format::Unknown: return "Unknown";
    }
    return "Unknown";
}

//==============================================================================
// Conversion functions
//==============================================================================

void convert_parquet_to_rearrangement(
    const std::string& input,
    const std::string& output,
    igor::airr::Delimiter delim)
{
    std::cout << "Reading Parquet file: " << input << std::endl;
    auto sequences = igor::ParquetReader::read_sequences(input);

    std::cout << "Writing AIRR Rearrangement file: " << output << std::endl;
    std::cout << "  Sequences: " << sequences.size() << std::endl;

    // Convert to SequenceData vector
    std::vector<igor::SequenceData> seq_data;
    seq_data.reserve(sequences.size());
    for (const auto& [id, seq, alignments] : sequences) {
        seq_data.emplace_back(id, seq, alignments);
    }

    igor::airr::rearrangement::write(output, seq_data, delim);
    std::cout << "✓ Conversion complete" << std::endl;
}

void convert_rearrangement_to_parquet(
    const std::string& input,
    const std::string& output,
    igor::CompressionType compression)
{
    std::cout << "Reading AIRR Rearrangement file: " << input << std::endl;
    auto sequences = igor::airr::rearrangement::read_sequences(input);

    std::cout << "Writing Parquet file: " << output << std::endl;
    std::cout << "  Sequences: " << sequences.size() << std::endl;
    std::cout << "  Compression: " << igor::ParquetWriter::compression_name(compression) << std::endl;

    // Convert SequenceData to tuple format
    std::vector<SequenceTuple> tuples;
    tuples.reserve(sequences.size());
    for (const auto& seq : sequences) {
        tuples.emplace_back(seq.index, seq.sequence, seq.alignments);
    }

    igor::ParquetWriter::write_sequences(output, tuples, compression);
    std::cout << "✓ Conversion complete" << std::endl;
}


void convert_parquet_to_alignment(
    const std::string& input,
    const std::string& output,
    igor::airr::Delimiter delim)
{
    std::cout << "Reading Parquet file: " << input << std::endl;
    auto sequences = igor::ParquetReader::read_sequences(input);

    std::cout << "Writing AIRR Alignment file: " << output << std::endl;

    // Convert to SequenceData vector
    std::vector<igor::SequenceData> seq_data;
    seq_data.reserve(sequences.size());
    for (const auto& [id, seq, alignments] : sequences) {
        seq_data.emplace_back(id, seq, alignments);
    }

    igor::airr::alignment::write_sequences(output, seq_data, delim);
    std::cout << "✓ Conversion complete" << std::endl;
}

void convert_alignment_to_parquet(
    const std::string& input,
    const std::string& output,
    igor::CompressionType compression)
{
    std::cout << "Reading AIRR Alignment file: " << input << std::endl;
    auto sequences = igor::airr::alignment::read_sequences(input);

    std::cout << "Writing Parquet file: " << output << std::endl;
    std::cout << "  Sequences: " << sequences.size() << std::endl;
    std::cout << "  Compression: " << igor::ParquetWriter::compression_name(compression) << std::endl;

    // Convert SequenceData to tuple format
    std::vector<SequenceTuple> tuples;
    tuples.reserve(sequences.size());
    for (const auto& seq : sequences) {
        tuples.emplace_back(seq.index, seq.sequence, seq.alignments);
    }

    igor::ParquetWriter::write_sequences(output, tuples, compression);
    std::cout << "✓ Conversion complete" << std::endl;
}


void convert_rearrangement_to_alignment(
    const std::string& input,
    const std::string& output,
    igor::airr::Delimiter out_delim)
{
    std::cout << "Reading AIRR Rearrangement file: " << input << std::endl;
    auto sequences = igor::airr::rearrangement::read_sequences(input);

    std::cout << "Writing AIRR Alignment file: " << output << std::endl;
    igor::airr::alignment::write_sequences(output, sequences, out_delim);
    std::cout << "✓ Conversion complete" << std::endl;
}

void convert_alignment_to_rearrangement(
    const std::string& input,
    const std::string& output,
    igor::airr::Delimiter out_delim)
{
    std::cout << "Reading AIRR Alignment file: " << input << std::endl;
    auto sequences = igor::airr::alignment::read_sequences(input);

    std::cout << "Writing AIRR Rearrangement file: " << output << std::endl;
    igor::airr::rearrangement::write(output, sequences, out_delim);
    std::cout << "✓ Conversion complete" << std::endl;
}

void convert_legacy_to_rearrangement(
    const std::string& input,
    const std::string& output,
    igor::airr::Delimiter out_delim)
{
    std::cout << "Reading IGoR Legacy file: " << input << std::endl;

    // Read legacy indexed_seq format (seq_index;sequence)
    std::ifstream infile(input);
    if (!infile) {
        throw std::runtime_error("Cannot open file: " + input);
    }

    std::vector<SequenceData> sequences;
    std::string line;

    // Skip header
    std::getline(infile, line);

    // Read sequences
    while (std::getline(infile, line)) {
        auto pos = line.find(';');
        if (pos == std::string::npos) continue;

        int seq_id = std::stoi(line.substr(0, pos));
        std::string seq = line.substr(pos + 1);

        // Create minimal SequenceData (no alignment info in indexed_seq)
        SequenceData seq_data;
        seq_data.index = seq_id;
        seq_data.sequence = seq;
        // alignments will be empty

        sequences.push_back(seq_data);
    }

    std::cout << "  Sequences: " << sequences.size() << std::endl;
    std::cout << "Writing AIRR Rearrangement file: " << output << std::endl;
    igor::airr::rearrangement::write(output, sequences, out_delim);
    std::cout << "✓ Conversion complete" << std::endl;
}

void convert_legacy_to_parquet(
    const std::string& input,
    const std::string& output,
    igor::CompressionType compression)
{
    std::cout << "Reading IGoR Legacy file: " << input << std::endl;

    // Read legacy indexed_seq format
    std::ifstream infile(input);
    if (!infile) {
        throw std::runtime_error("Cannot open file: " + input);
    }

    std::vector<SequenceData> sequences;
    std::string line;

    // Skip header
    std::getline(infile, line);

    // Read sequences
    while (std::getline(infile, line)) {
        auto pos = line.find(';');
        if (pos == std::string::npos) continue;

        int seq_id = std::stoi(line.substr(0, pos));
        std::string seq = line.substr(pos + 1);

        SequenceData seq_data;
        seq_data.index = seq_id;
        seq_data.sequence = seq;

        sequences.push_back(seq_data);
    }

    std::cout << "  Sequences: " << sequences.size() << std::endl;
    std::cout << "Writing Parquet file: " << output << std::endl;
    
    // Convert SequenceData to tuples for Parquet writing
    std::vector<SequenceTuple> tuples;
    tuples.reserve(sequences.size());
    for (const auto& seq : sequences) {
        tuples.emplace_back(seq.index, seq.sequence, seq.alignments);
    }
    
    ParquetWriter::write_sequences(output, tuples, compression);
    std::cout << "✓ Conversion complete" << std::endl;
}

//==============================================================================
// Help and usage
//==============================================================================

void print_usage(const char* program)
{
    std::cout << R"(
IGoR Format Converter
=====================

Usage: )" << program << R"( [options] <input> <output>

Convert between IGoR data formats:
  - Parquet (high-performance binary)
  - AIRR Rearrangement (TSV/CSV, one row per sequence)
  - AIRR Alignment (TSV/CSV, one row per alignment)

Options:
  --compression <type>    Parquet compression (NONE, SNAPPY, GZIP, ZSTD, LZ4)
                          Default: SNAPPY

  --delimiter <type>      AIRR delimiter (TAB, COMMA)
                          Default: TAB for .tsv, COMMA for .csv

  --help                  Show this help message

Examples:
  # Parquet to AIRR Rearrangement
  )" << program << R"( sequences.parquet output.tsv

  # AIRR Rearrangement to Parquet with ZSTD compression
  )" << program << R"( --compression ZSTD input.tsv output.parquet

  # Parquet to AIRR Alignment
  )" << program << R"( sequences.parquet alignments.tsv

  # AIRR Rearrangement to AIRR Alignment
  )" << program << R"( rearrangement.tsv alignment.tsv

Supported Conversions:
  Parquet          ⇄  AIRR Rearrangement
  Parquet          ⇄  AIRR Alignment
  AIRR Rearrangement ⇄  AIRR Alignment

Format auto-detection:
  .parquet          → Parquet
  .tsv/.csv         → Auto-detect AIRR schema from header
)";
}

//==============================================================================
// Main
//==============================================================================

int main(int argc, char* argv[])
{
    // Parse arguments
    std::string input;
    std::string output;
    igor::CompressionType compression = igor::CompressionType::SNAPPY;
    igor::airr::Delimiter delimiter = igor::airr::Delimiter::AUTO;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
        else if (arg == "--compression") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --compression requires an argument" << std::endl;
                return 1;
            }
            std::string comp = argv[++i];
            if (comp == "NONE") compression = igor::CompressionType::NONE;
            else if (comp == "SNAPPY") compression = igor::CompressionType::SNAPPY;
            else if (comp == "GZIP") compression = igor::CompressionType::GZIP;
            else if (comp == "ZSTD") compression = igor::CompressionType::ZSTD;
            else if (comp == "LZ4") compression = igor::CompressionType::LZ4;
            else {
                std::cerr << "Error: Unknown compression type: " << comp << std::endl;
                return 1;
            }
        }
        else if (arg == "--delimiter") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --delimiter requires an argument" << std::endl;
                return 1;
            }
            std::string delim = argv[++i];
            if (delim == "TAB") delimiter = igor::airr::Delimiter::TAB;
            else if (delim == "COMMA") delimiter = igor::airr::Delimiter::COMMA;
            else {
                std::cerr << "Error: Unknown delimiter type: " << delim << std::endl;
                return 1;
            }
        }
        else if (input.empty()) {
            input = arg;
        }
        else if (output.empty()) {
            output = arg;
        }
        else {
            std::cerr << "Error: Too many arguments" << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    // Validate arguments
    if (input.empty() || output.empty()) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    if (!fs::exists(input)) {
        std::cerr << "Error: Input file does not exist: " << input << std::endl;
        return 1;
    }

    // Detect formats
    auto input_fmt = detect_format(input);
    auto output_fmt = detect_format(output, true);  // true = is_output

    if (input_fmt == Format::Unknown) {
        std::cerr << "Error: Cannot determine input format for: " << input << std::endl;
        return 1;
    }

    if (output_fmt == Format::Unknown) {
        std::cerr << "Error: Cannot determine output format for: " << output << std::endl;
        return 1;
    }

    std::cout << "Converting: " << format_name(input_fmt)
              << " → " << format_name(output_fmt) << std::endl;
    std::cout << std::endl;

    // Auto-detect delimiter from output extension if not specified
    if (delimiter == igor::airr::Delimiter::AUTO) {
        fs::path out_path(output);
        auto ext = out_path.extension().string();
        if (ext == ".csv") {
            delimiter = igor::airr::Delimiter::COMMA;
        } else {
            delimiter = igor::airr::Delimiter::TAB;
        }
    }

    // Perform conversion
    try {
        if (input_fmt == Format::Parquet && output_fmt == Format::AIRRRearrangement) {
            convert_parquet_to_rearrangement(input, output, delimiter);
        }
        else if (input_fmt == Format::AIRRRearrangement && output_fmt == Format::Parquet) {
            convert_rearrangement_to_parquet(input, output, compression);
        }
        else if (input_fmt == Format::Parquet && output_fmt == Format::AIRRAlignment) {
            convert_parquet_to_alignment(input, output, delimiter);
        }
        else if (input_fmt == Format::AIRRAlignment && output_fmt == Format::Parquet) {
            convert_alignment_to_parquet(input, output, compression);
        }
        else if (input_fmt == Format::AIRRRearrangement && output_fmt == Format::AIRRAlignment) {
            convert_rearrangement_to_alignment(input, output, delimiter);
        }
        else if (input_fmt == Format::AIRRAlignment && output_fmt == Format::AIRRRearrangement) {
            convert_alignment_to_rearrangement(input, output, delimiter);
        }
        // Legacy IGoR format conversions
        else if (input_fmt == Format::IGoRLegacy && output_fmt == Format::AIRRRearrangement) {
            convert_legacy_to_rearrangement(input, output, delimiter);
        }
        else if (input_fmt == Format::IGoRLegacy && output_fmt == Format::Parquet) {
            convert_legacy_to_parquet(input, output, compression);
        }
        else if (input_fmt == Format::IGoRLegacy && output_fmt == Format::AIRRAlignment) {
            // Legacy → Rearrangement → Alignment
            std::string temp = output + ".tmp.tsv";
            convert_legacy_to_rearrangement(input, temp, igor::airr::Delimiter::TAB);
            convert_rearrangement_to_alignment(temp, output, delimiter);
            fs::remove(temp);
        }
        else {
            std::cerr << "Error: Unsupported conversion: "
                      << format_name(input_fmt) << " → " << format_name(output_fmt) << std::endl;
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error during conversion: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
