/*
 * FastGenerator.cpp
 *
 *  Created on: Feb 4, 2026
 *      Fast parallel sequence generation
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <igor/Core/FastGenerator.h>
#include <igor/Core/EventUtils.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Dinuclmarkov.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace igor {
namespace fast {

void FastGenerator::initialize(const Model_Parms &model_parms, const Model_marginals &model_marginals)
{
    event_samplers_.clear();
    event_name_to_index_.clear();

    // Get model queue and index map
    std::queue<std::shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
    std::unordered_map<Rec_Event_name, int> index_map = model_marginals.get_index_map(model_parms, model_queue);

    // Get inverse offset map for dependencies (tells us for each event, which parents affect it)
    model_queue = model_parms.get_model_queue();
    auto inverse_offset_map = model_marginals.get_inverse_offset_map(model_parms, model_queue);

    // Reset queue and count events
    model_queue = model_parms.get_model_queue();
    size_t num_events = 0;
    {
        auto temp_queue = model_queue;
        while (!temp_queue.empty()) {
            num_events++;
            temp_queue.pop();
        }
    }

    // Check for D gene
    has_d_gene_ = false;
    auto events_map = model_parms.get_events_map();
    if (events_map.count(std::make_tuple(GeneChoice_t, std::string("D_gene_seq"), Undefined_side)) > 0) {
        has_d_gene_ = true;
    }

    // First pass: create event samplers and build name-to-index map
    event_samplers_.reserve(num_events);
    size_t event_idx = 0;

    model_queue = model_parms.get_model_queue();
    while (!model_queue.empty()) {
        auto event = model_queue.front();
        model_queue.pop();

        FastEventSampler sampler;
        sampler.type = event->get_type();
        sampler.gene_class = event->get_class();
        sampler.side = event->get_side();
        sampler.name = event->get_name();
        sampler.seq_type = event->get_seq_type();
        sampler.event_index = event_idx;
        sampler.num_realizations = event->size();

        event_name_to_index_[sampler.name] = event_idx;
        event_samplers_.push_back(std::move(sampler));
        event_idx++;
    }

    // Second pass: set up dependencies and initialize samplers
    model_queue = model_parms.get_model_queue();
    event_idx = 0;

    while (!model_queue.empty()) {
        auto event = model_queue.front();
        model_queue.pop();

        FastEventSampler &sampler = event_samplers_[event_idx];
        int base_index = index_map.at(event->get_name());

        // Build dependencies from inverse_offset_map
        // The inverse_offset_map gives (parent_event, stride) pairs
        // where stride is: this_event_size * product_of_earlier_parent_sizes
        // The condition_index = sum(parent_choice * stride / this_event_size)
        if (inverse_offset_map.count(event->get_name()) > 0) {
            const auto &parents = inverse_offset_map.at(event->get_name());

            // Calculate total number of conditions from the strides
            // The largest stride / num_realizations gives us total conditions
            size_t max_stride = sampler.num_realizations;

            for (const auto &parent_stride : parents) {
                const std::string &parent_name = parent_stride.first->get_name();
                size_t stride = static_cast<size_t>(parent_stride.second);

                if (event_name_to_index_.count(parent_name) > 0) {
                    EventDependency dep;
                    dep.parent_event_idx = event_name_to_index_.at(parent_name);
                    // Normalize stride: divide by num_realizations to get the condition multiplier
                    dep.stride = stride / sampler.num_realizations;
                    sampler.dependencies.push_back(dep);

                    // Track max stride to compute total conditions
                    size_t parent_size = event_samplers_[dep.parent_event_idx].num_realizations;
                    if (stride * parent_size > max_stride) {
                        max_stride = stride * parent_size;
                    }
                }
            }

            sampler.num_conditions = max_stride / sampler.num_realizations;
            if (sampler.num_conditions == 0)
                sampler.num_conditions = 1;
        }

        // Initialize based on event type
        if (sampler.type == Dinuclmarkov_t) {
            initialize_dinucl_sampler(sampler, event, model_marginals.marginal_array_smart_p, base_index);
        } else {
            initialize_event_sampler(sampler, event, model_marginals.marginal_array_smart_p, base_index);
        }

        sampler.is_initialized = true;
        event_idx++;
    }

    initialized_ = true;
}

void FastGenerator::initialize_event_sampler(FastEventSampler &sampler, const std::shared_ptr<Rec_Event> &event,
                                             const Marginal_array_p &marginals, int base_index)
{
    auto realizations = event->get_realizations_map();
    size_t num_realizations = realizations.size();

    // For GeneChoice: store gene sequences
    if (sampler.type == GeneChoice_t) {
        sampler.gene_sequences.resize(num_realizations);
        sampler.gene_sequences_int.resize(num_realizations);

        for (const auto &pair : realizations) {
            const Event_realization &real = pair.second;
            sampler.gene_sequences[real.index] = real.value_str;
            sampler.gene_sequences_int[real.index] = real.value_str_int;
        }
    }

    // For Insertion: find max insertion length
    if (sampler.type == Insertion_t) {
        sampler.max_insertions = 0;
        for (const auto &pair : realizations) {
            if (pair.second.value_int > sampler.max_insertions) {
                sampler.max_insertions = pair.second.value_int;
            }
        }
    }

    // For Deletion: build index-to-value mapping
    if (sampler.type == Deletion_t) {
        sampler.min_deletion = INT_MAX;
        sampler.max_deletion = INT_MIN;
        sampler.deletion_idx_to_value.resize(num_realizations);

        for (const auto &pair : realizations) {
            const Event_realization &real = pair.second;
            int val = real.value_int;
            if (val < sampler.min_deletion)
                sampler.min_deletion = val;
            if (val > sampler.max_deletion)
                sampler.max_deletion = val;
            if (real.index >= 0 && static_cast<size_t>(real.index) < num_realizations) {
                sampler.deletion_idx_to_value[real.index] = val;
            }
        }
    }

    // Build conditional probability matrix
    // Layout: marginals[base_index + condition_idx * num_realizations + realization_idx]
    size_t num_conditions = sampler.num_conditions;
    if (num_conditions == 0)
        num_conditions = 1;

    std::vector<std::vector<double>> cond_probs(num_conditions);

    for (size_t c = 0; c < num_conditions; ++c) {
        cond_probs[c].resize(num_realizations, 0.0);
        double sum = 0.0;

        for (const auto &pair : realizations) {
            const Event_realization &real = pair.second;
            double p = static_cast<double>(marginals[base_index + c * num_realizations + real.index]);
            cond_probs[c][real.index] = p;
            sum += p;
        }

        // Normalize
        if (sum > 0) {
            for (size_t i = 0; i < num_realizations; ++i) {
                cond_probs[c][i] /= sum;
            }
        }
    }

    sampler.sampler.initialize(cond_probs, true);
}

void FastGenerator::initialize_dinucl_sampler(FastEventSampler &sampler, const std::shared_ptr<Rec_Event> & /* event */,
                                              const Marginal_array_p &marginals, int base_index)
{
    // The dinucleotide model has a 15x15 or 4x4 transition matrix
    // We need to extract it from the marginals

    // For a 4x4 standard nucleotide matrix:
    // P(next=j | prev=i) is stored at marginals[base_index + i*4 + j]
    // or using IGoR's convention: marginals[base_index + i + 4*j]

    const size_t n = 4; // Standard nucleotides only for generation
    std::vector<double> dinuc_probs(n * n);

    // Extract the transition probabilities
    // IGoR stores these as: marginals[base_index + prev_nt + n_nucleotides * next_nt]
    for (size_t prev = 0; prev < n; ++prev) {
        double row_sum = 0.0;
        for (size_t next = 0; next < n; ++next) {
            double p = static_cast<double>(marginals[base_index + prev + n * next]);
            dinuc_probs[prev * n + next] = p;
            row_sum += p;
        }
        // Normalize each row
        if (row_sum > 0) {
            for (size_t next = 0; next < n; ++next) {
                dinuc_probs[prev * n + next] /= row_sum;
            }
        }
    }

    sampler.dinuc_sampler.initialize(dinuc_probs.data(), n);
}

size_t FastGenerator::compute_condition_index(const FastEventSampler &sampler, const SamplingContext &context) const
{
    if (sampler.dependencies.empty()) {
        return 0; // No dependencies = unconditional
    }

    size_t condition_idx = 0;
    for (const auto &dep : sampler.dependencies) {
        size_t parent_choice = context.sampled_indices[dep.parent_event_idx];
        condition_idx += parent_choice * dep.stride;
    }

    return condition_idx;
}

void FastGenerator::apply_gene_choice(const FastEventSampler &sampler, size_t choice_idx,
                                      std::unordered_map<Seq_type, std::string> &sequences) const
{
    Seq_type seq_type = V_gene_seq;
    if (!EventUtils::try_gene_class_to_gene_seq_type(sampler.gene_class, seq_type)) {
        return;
    }

    if (choice_idx < sampler.gene_sequences.size()) {
        sequences[seq_type] = sampler.gene_sequences[choice_idx];
    }
}

void FastGenerator::apply_deletion(const FastEventSampler &sampler, size_t del_idx,
                                   std::unordered_map<Seq_type, std::string> &sequences) const
{
    if (del_idx >= sampler.deletion_idx_to_value.size()) {
        return;
    }

    int num_del = sampler.deletion_idx_to_value[del_idx];

    Seq_type seq_type = V_gene_seq;
    if (!EventUtils::try_gene_class_to_gene_seq_type(sampler.gene_class, seq_type)) {
        return;
    }

    if (sequences.count(seq_type) == 0) {
        return;
    }

    std::string &seq = sequences[seq_type];

    if (num_del >= 0) {
        // Positive deletion: remove bases
        if (sampler.side == Three_prime) {
            if (static_cast<size_t>(num_del) < seq.size()) {
                seq.erase(seq.size() - num_del);
            }
        } else if (sampler.side == Five_prime) {
            if (static_cast<size_t>(num_del) < seq.size()) {
                seq.erase(0, num_del);
            }
        }
    } else {
        // Negative deletion = palindromic insertion
        int num_palindrome = -num_del;
        if (sampler.side == Three_prime) {
            apply_palindrome(seq, num_palindrome, true);
        } else if (sampler.side == Five_prime) {
            apply_palindrome(seq, num_palindrome, false);
        }
    }
}

void FastGenerator::sample_event(const FastEventSampler &sampler, std::mt19937_64 &rng,
                                 std::unordered_map<Seq_type, std::string> &sequences, std::vector<int> &realization,
                                 SamplingContext &context) const
{
    // Compute condition index from all parent choices
    size_t condition_idx = compute_condition_index(sampler, context);

    // Sample from conditional distribution
    size_t choice_idx = sampler.sampler.sample(condition_idx, rng);

    // Store sampled index for downstream conditionals
    context.sampled_indices[sampler.event_index] = choice_idx;

    // Record realization
    realization.push_back(static_cast<int>(choice_idx));

    // Apply the sampled choice to sequences based on event type
    switch (sampler.type) {
    case GeneChoice_t:
        apply_gene_choice(sampler, choice_idx, sequences);
        break;

    case Deletion_t:
        apply_deletion(sampler, choice_idx, sequences);
        break;

    case Insertion_t: {
        // For insertion, sampled index directly gives insertion length
        // Create placeholder for dinucleotide model
        Seq_type seq_type = VD_ins_seq;
        if (!EventUtils::try_insertion_gene_class_to_seq_type(sampler.gene_class, seq_type)) {
            return;
        }
        sequences[seq_type] = std::string(choice_idx, 'I');
        break;
    }

    default:
        break;
    }
}

void FastGenerator::sample_dinucl_markov(const FastEventSampler &sampler, std::mt19937_64 &rng,
                                         std::unordered_map<Seq_type, std::string> &sequences,
                                         std::vector<int> &realization, SamplingContext &context) const
{
    static const char NT_CHARS[] = { 'A', 'C', 'G', 'T' };
    static const int CHAR_TO_INT[] = {
        0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3
    }; // A=0, C=2, G=6, T=19 offset from 'A'

    auto get_last_nt = [&](const std::string &s) -> int {
        if (s.empty())
            return 0; // Default to A
        char c = s.back();
        int idx = c - 'A';
        if (idx >= 0 && idx < 20 && CHAR_TO_INT[idx] >= 0) {
            return CHAR_TO_INT[idx];
        }
        return 0;
    };

    auto get_first_nt = [&](const std::string &s) -> int {
        if (s.empty())
            return 0; // Default to A
        char c = s.front();
        int idx = c - 'A';
        if (idx >= 0 && idx < 20 && CHAR_TO_INT[idx] >= 0) {
            return CHAR_TO_INT[idx];
        }
        return 0;
    };

    // Process based on seq_type string
    const std::string &st = sampler.seq_type;
    if (st == "VD_ins_seq") {
        std::string &ins_seq = sequences[VD_ins_seq];
        if (!ins_seq.empty()) {
            int prev_nt = get_last_nt(sequences[V_gene_seq]);
            for (size_t i = 0; i < ins_seq.size(); ++i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }

    if (st == "DJ_ins_seq") {
        std::string &ins_seq = sequences[DJ_ins_seq];
        if (!ins_seq.empty()) {
            // DJ insertions are generated from J side going towards D
            int prev_nt = get_first_nt(sequences[J_gene_seq]);
            // Generate in reverse order
            for (int i = static_cast<int>(ins_seq.size()) - 1; i >= 0; --i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }

    if (st == "VJ_ins_seq") {
        std::string &ins_seq = sequences[VJ_ins_seq];
        if (!ins_seq.empty()) {
            int prev_nt = get_last_nt(sequences[V_gene_seq]);
            for (size_t i = 0; i < ins_seq.size(); ++i) {
                if (ins_seq[i] == 'I') {
                    int next_nt = sampler.dinuc_sampler.sample_next(prev_nt, rng);
                    ins_seq[i] = NT_CHARS[next_nt];
                    realization.push_back(next_nt);
                    prev_nt = next_nt;
                }
            }
        }
    }
}

void FastGenerator::apply_palindrome(std::string &seq, int num_bases, bool is_3prime)
{
    auto complement_char = [](char c) -> char {
        switch (c) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return c;
        }
    };

    if (num_bases <= 0 || seq.empty())
        return;

    std::string palindrome;
    if (is_3prime) {
        // Take last num_bases from seq, reverse complement
        size_t start = (num_bases <= static_cast<int>(seq.size())) ? seq.size() - num_bases : 0;
        size_t len = (std::min)(static_cast<size_t>(num_bases), seq.size());
        palindrome = seq.substr(start, len);
        std::reverse(palindrome.begin(), palindrome.end());
        for (char &c : palindrome) {
            c = complement_char(c);
        }
        seq += palindrome;
    } else {
        // Take first num_bases from seq, reverse complement, prepend
        size_t len = (std::min)(static_cast<size_t>(num_bases), seq.size());
        palindrome = seq.substr(0, len);
        std::reverse(palindrome.begin(), palindrome.end());
        for (char &c : palindrome) {
            c = complement_char(c);
        }
        seq = palindrome + seq;
    }
}

std::string FastGenerator::assemble_sequence(const std::unordered_map<Seq_type, std::string> &sequences) const
{
    std::string result;
    result.reserve(500); // Typical sequence length

    // Assemble in order: V + VD_ins + D + DJ_ins + J (for VDJ)
    // or: V + VJ_ins + J (for VJ)

    auto get_seq = [&sequences](Seq_type t) -> const std::string & {
        static const std::string empty;
        auto it = sequences.find(t);
        return (it != sequences.end()) ? it->second : empty;
    };

    result += get_seq(V_gene_seq);

    if (has_d_gene_) {
        result += get_seq(VD_ins_seq);
        result += get_seq(D_gene_seq);
        result += get_seq(DJ_ins_seq);
    } else {
        result += get_seq(VJ_ins_seq);
    }

    result += get_seq(J_gene_seq);

    return result;
}

void FastGenerator::generate_single(std::mt19937_64 &rng, GeneratedSequence &result) const
{
    result.clear();

    std::unordered_map<Seq_type, std::string> sequences;
    sequences.reserve(6);

    result.realizations.reserve(event_samplers_.size());

    // Create sampling context for conditional dependencies
    SamplingContext context;
    context.resize(event_samplers_.size());

    // Process each event in order
    for (const auto &sampler : event_samplers_) {
        std::vector<int> event_realization;
        event_realization.reserve(4);

        if (sampler.type == Dinuclmarkov_t) {
            // Dinucleotide Markov has special handling (generates sequences)
            sample_dinucl_markov(sampler, rng, sequences, event_realization, context);
        } else {
            // All other events use generic sampler
            sample_event(sampler, rng, sequences, event_realization, context);
        }

        result.realizations.push_back(std::move(event_realization));
    }

    result.sequence = assemble_sequence(sequences);
}

std::vector<GeneratedSequence> FastGenerator::generate(size_t num_sequences, const FastGeneratorConfig &config,
                                                       SequenceCallback callback, ProgressCallback progress)
{
    if (!initialized_) {
        throw std::runtime_error("FastGenerator not initialized. Call initialize() first.");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Determine thread count
    size_t num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = get_optimal_thread_count();
    }

    // Get base seed
    uint64_t base_seed = config.base_seed;
    if (base_seed == 0) {
        base_seed = draw_random_seed();
    }

    std::vector<GeneratedSequence> results;
    if (!callback) {
        results.resize(num_sequences);
    }

    std::atomic<size_t> completed{ 0 };
    std::atomic<size_t> next_progress_report{ 0 };
    const size_t progress_interval = (std::max)(size_t(1000), num_sequences / 100);

#pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::mt19937_64 rng(base_seed + thread_id * 1000000ULL);
        GeneratedSequence local_result;

#pragma omp for schedule(dynamic, config.batch_size)
        for (auto i = 0; i < num_sequences; ++i) {
            generate_single(rng, local_result);

            if (callback) {
                callback(i, local_result.sequence, local_result.realizations);
            } else {
                results[i] = local_result;
            }

            size_t done = ++completed;

            // Progress reporting
            if (progress && config.show_progress) {
                size_t expected = next_progress_report.load();
                if (done >= expected) {
                    if (next_progress_report.compare_exchange_strong(expected, expected + progress_interval)) {
                        progress(done, num_sequences);
                    }
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();

    stats_.sequences_generated = num_sequences;
    stats_.total_time_seconds = elapsed;
    stats_.sequences_per_second = num_sequences / elapsed;

    return results;
}

void FastGenerator::generate_to_files(size_t num_sequences, const std::string &seq_filename,
                                      const std::string &real_filename, const FastGeneratorConfig &config,
                                      ProgressCallback progress)
{
    if (!initialized_) {
        throw std::runtime_error("FastGenerator not initialized. Call initialize() first.");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Determine thread count
    size_t num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = get_optimal_thread_count();
    }

    // Get base seed
    uint64_t base_seed = config.base_seed;
    if (base_seed == 0) {
        base_seed = draw_random_seed();
    }

    // Phase 1: Generate all sequences in parallel WITHOUT any I/O
    std::vector<GeneratedSequence> results(num_sequences);

    std::atomic<size_t> completed{ 0 };
    std::atomic<size_t> next_progress_report{ 0 };
    const size_t progress_interval = (std::max)(size_t(1000), num_sequences / 100);

    auto gen_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::mt19937_64 rng(base_seed + thread_id * 1000000ULL);
        GeneratedSequence local_result;

#pragma omp for schedule(static)
        for (auto i = 0; i < num_sequences; ++i) {
            generate_single(rng, local_result);
            results[i] = std::move(local_result);

            size_t done = ++completed;

            // Progress reporting (generation phase: 0-80%)
            if (progress && config.show_progress) {
                size_t expected = next_progress_report.load();
                if (done >= expected) {
                    if (next_progress_report.compare_exchange_strong(expected, expected + progress_interval)) {
                        progress(done * 80 / 100, num_sequences);
                    }
                }
            }
        }
    }

    auto gen_end = std::chrono::high_resolution_clock::now();
    double gen_time = std::chrono::duration<double>(gen_end - gen_start).count();

    // Phase 2: Write to files sequentially with large buffers
    auto io_start = std::chrono::high_resolution_clock::now();

    const size_t write_buffer_size = 4 * 1024 * 1024; // 4MB buffer

    std::ofstream seq_file(seq_filename, std::ios::binary);
    std::ofstream real_file(real_filename, std::ios::binary);

    if (!seq_file || !real_file) {
        throw std::runtime_error("Failed to open output files");
    }

    // Write headers
    seq_file << "seq_index;nt_sequence\n";

    real_file << "seq_index";
    for (const auto &sampler : event_samplers_) {
        real_file << ";" << sampler.name;
    }
    real_file << "\n";

    // Write data with buffering
    std::string seq_buffer, real_buffer;
    seq_buffer.reserve(write_buffer_size);
    real_buffer.reserve(write_buffer_size);

    size_t total_bytes = 0;

    for (size_t i = 0; i < num_sequences; ++i) {
        const auto &result = results[i];

        // Format sequence line
        seq_buffer += std::to_string(i);
        seq_buffer += ';';
        seq_buffer += result.sequence;
        seq_buffer += '\n';

        // Format realizations line
        real_buffer += std::to_string(i);
        for (const auto &event_real : result.realizations) {
            real_buffer += ";(";
            for (size_t j = 0; j < event_real.size(); ++j) {
                if (j > 0)
                    real_buffer += ',';
                real_buffer += std::to_string(event_real[j]);
            }
            real_buffer += ')';
        }
        real_buffer += '\n';

        // Flush buffers when they're large enough
        if (seq_buffer.size() >= write_buffer_size) {
            seq_file.write(seq_buffer.data(), seq_buffer.size());
            total_bytes += seq_buffer.size();
            seq_buffer.clear();
        }
        if (real_buffer.size() >= write_buffer_size) {
            real_file.write(real_buffer.data(), real_buffer.size());
            total_bytes += real_buffer.size();
            real_buffer.clear();
        }

        // Progress reporting (write phase: 80-100%)
        if (progress && config.show_progress && (i % progress_interval == 0)) {
            progress(num_sequences * 80 / 100 + i * 20 / 100, num_sequences);
        }
    }

    // Final flush
    if (!seq_buffer.empty()) {
        seq_file.write(seq_buffer.data(), seq_buffer.size());
        total_bytes += seq_buffer.size();
    }
    if (!real_buffer.empty()) {
        real_file.write(real_buffer.data(), real_buffer.size());
        total_bytes += real_buffer.size();
    }

    auto io_end = std::chrono::high_resolution_clock::now();
    double io_time = std::chrono::duration<double>(io_end - io_start).count();

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();

    // Print timing breakdown
    std::cerr << "  [Timing] Generation: " << gen_time << "s, I/O: " << io_time << "s, Total: " << elapsed << "s"
              << std::endl;

    stats_.sequences_generated = num_sequences;
    stats_.total_time_seconds = elapsed;
    stats_.sequences_per_second = num_sequences / elapsed;
    stats_.bytes_written = total_bytes;
}

} // namespace fast
} // namespace igor
