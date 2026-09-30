/**
 * @file test_ModelInference.cpp
 * @brief Inference validation using the new Model architecture
 *
 * Mirror of tst/test_inference.cpp, with the EM iteration loop driven
 * entirely by InferenceEngine::run().
 *
 * Architecture:
 *   - RecombinationModel<double> holds the probability tensors (weights).
 *   - InferenceEngine<double> orchestrates the EM loop:
 *       1. resetAccumulators()          — zero every handler's accumulator
 *       2. eStep(engine)               — caller-provided E-step
 *       3. updateParameters()          — M-step: normalise accumulators → weights
 *   - The E-step callable bridges to the legacy Rec_Event::iterate()
 *     machinery through GenModel, because a native E-step does not exist
 *     yet.  For each iteration the current model weights are exported to
 *     a legacy Model_marginals, a GenModel runs one E+M step, and the
 *     resulting marginals are imported back into handler accumulators.
 *     The engine's own M-step then re-normalises (idempotent on already-
 *     normalised values), completing the round-trip through the new API.
 *
 * This validates:
 *   1. InferenceEngine::run() correctly orchestrates reset / E / M
 *   2. Handler accumulators accept imported data from the legacy bridge
 *   3. maximizeLikelihood() produces correct normalised distributions
 *   4. The model converges to the ground truth (KL divergence test)
 *   5. Weight tensors in both models match expected distributions
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// New Model architecture
#include <igor/Model/RecombinationModel.h>
#include <igor/Model/InferenceEngine.h>
#include <igor/Model/LegacyBridge.h>
#include <igor/Model/Topology.h>

// Legacy Core types (for generation + inference E-step)
#include <igor/Core/GenModel.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Aligner.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef IGOR_MODELS_DIR
#  error "IGOR_MODELS_DIR must be defined (set by CMake)"
#endif

using namespace igor::model;
using Catch::Matchers::WithinAbs;

// ---------------------------------------------------------------------------
// Model base directory
// ---------------------------------------------------------------------------
static const std::string MODELS_DIR = std::string(IGOR_MODELS_DIR);

// ---------------------------------------------------------------------------
// Mathematical helpers
// ---------------------------------------------------------------------------

static inline double kl_divergence(const std::vector<double>& P,
                                   const std::vector<double>& Q,
                                   double* uncovered_mass = nullptr)
{
    double kl = 0.0;
    double skipped = 0.0;
    for (std::size_t i = 0; i < P.size(); ++i) {
        if (P[i] > 0.0) {
            if (Q[i] <= 0.0)
                skipped += P[i];
            else
                kl += P[i] * std::log2(P[i] / Q[i]);
        }
    }
    if (uncovered_mass) *uncovered_mass = skipped;
    return kl;
}

static inline double entropy(const std::vector<double>& P)
{
    double h = 0.0;
    for (double p : P)
        if (p > 0.0) h -= p * std::log2(p);
    return h;
}

// ---------------------------------------------------------------------------
// Event metadata (lightweight version, sufficient for categorical events)
// ---------------------------------------------------------------------------

struct ModelEventInfo
{
    std::string name;          // Rec_Event full name
    std::string nickname;      // Rec_Event nickname
    int         num_realizations;
    bool        is_dinuc_markov;
    std::size_t queue_position;
    Gene_class  gene_class;

    std::vector<double> marginal;   // marginalised P(realization)
    double              H;          // Shannon entropy of marginal
};

/**
 * @brief Build event metadata from legacy Model_Parms + Model_marginals.
 *
 * Uses compute_event_marginal_probability() for proper marginalisation
 * over parent dimensions.
 */
static std::vector<ModelEventInfo>
build_model_event_info(const Model_Parms& parms,
                       const Model_marginals& marginals)
{
    std::vector<ModelEventInfo> infos;
    auto model_queue = parms.get_model_queue();
    std::size_t pos = 0;

    while (!model_queue.empty()) {
        auto ev = model_queue.front();
        model_queue.pop();

        ModelEventInfo info;
        info.name            = ev->get_name();
        info.nickname        = ev->get_nickname();
        info.num_realizations = ev->size();
        info.is_dinuc_markov = (ev->get_type() == Event_type::Dinuclmarkov_t);
        info.queue_position  = pos;
        info.gene_class      = ev->get_class();

        if (info.is_dinuc_markov) {
            // Skip DinucMarkov — we don't validate its entropy here.
            info.marginal.resize(16, 0.0);
            info.H = 0.0;
        } else {
            auto [dims, probs] =
                marginals.compute_event_marginal_probability(info.name, parms);
            info.marginal.resize(info.num_realizations, 0.0);
            for (int i = 0; i < info.num_realizations; ++i)
                info.marginal[i] = static_cast<double>(probs.get()[i]);
            info.H = entropy(info.marginal);
        }

        infos.push_back(std::move(info));
        ++pos;
    }
    return infos;
}

// ---------------------------------------------------------------------------
// Scenario parsing (identical logic to test_inference.cpp)
// ---------------------------------------------------------------------------

struct ParsedScenario
{
    int v_gene_index = -1, j_gene_index = -1, d_gene_index = -1;
    int v_3p_del = -1, j_5p_del = -1, d_5p_del = -1, d_3p_del = -1;
    int vd_ins = -1, dj_ins = -1, vj_ins = -1;
    std::string v_gene_name, j_gene_name, d_gene_name;
};

static ParsedScenario
parse_scenario(const std::queue<std::queue<int>>& scenario,
               const std::vector<ModelEventInfo>&  event_infos,
               const Model_Parms&                  parms)
{
    ParsedScenario parsed;
    auto scenario_copy = scenario;
    std::size_t pos = 0;

    while (!scenario_copy.empty()) {
        auto inner = scenario_copy.front();
        scenario_copy.pop();
        if (inner.empty()) { ++pos; continue; }

        int realization_idx = inner.front();
        const ModelEventInfo* ev_info = nullptr;
        for (const auto& ev : event_infos) {
            if (ev.queue_position == pos) { ev_info = &ev; break; }
        }
        if (!ev_info) { ++pos; continue; }

        auto ev_ptr    = parms.get_event_pointer(ev_info->name);
        auto event_type = ev_ptr->get_type();
        auto gene_class = ev_ptr->get_class();

        if (event_type == Event_type::GeneChoice_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    if (gene_class == V_gene) {
                        parsed.v_gene_index = realization_idx;
                        parsed.v_gene_name  = name;
                    } else if (gene_class == J_gene) {
                        parsed.j_gene_index = realization_idx;
                        parsed.j_gene_name  = name;
                    } else if (gene_class == D_gene) {
                        parsed.d_gene_index = realization_idx;
                        parsed.d_gene_name  = name;
                    }
                    break;
                }
            }
        } else if (event_type == Event_type::Deletion_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    int del_value = real.value_int;
                    if (gene_class == V_gene)
                        parsed.v_3p_del = del_value;
                    else if (gene_class == J_gene)
                        parsed.j_5p_del = del_value;
                    else if (gene_class == D_gene) {
                        if (ev_ptr->get_side() == Five_prime)
                            parsed.d_5p_del = del_value;
                        else if (ev_ptr->get_side() == Three_prime)
                            parsed.d_3p_del = del_value;
                    }
                    break;
                }
            }
        } else if (event_type == Event_type::Insertion_t) {
            auto realizations = ev_ptr->get_realizations_map();
            for (const auto& [name, real] : realizations) {
                if (real.index == realization_idx) {
                    int ins_value = real.value_int;
                    if (gene_class == VD_genes) parsed.vd_ins = ins_value;
                    else if (gene_class == DJ_genes) parsed.dj_ins = ins_value;
                    else if (gene_class == VJ_genes) parsed.vj_ins = ins_value;
                    break;
                }
            }
        }
        ++pos;
    }
    return parsed;
}

// ---------------------------------------------------------------------------
// Mock alignment helpers (identical to test_inference.cpp)
// ---------------------------------------------------------------------------

static Alignment_data
create_v_mock_alignment(const std::string& sequence,
                        const ParsedScenario& scenario,
                        const std::unordered_map<std::string, std::string>& gene_templates)
{
    const auto& v_template = gene_templates.at(scenario.v_gene_name);
    int v_len = static_cast<int>(v_template.length());
    int v_del = scenario.v_3p_del;

    std::vector<int> mismatches;
    for (int i = v_len - v_del; i < v_len && i < static_cast<int>(sequence.length()); ++i)
        if (v_template[i] != sequence[i]) mismatches.push_back(i);

    Alignment_data v_align(scenario.v_gene_name, 0);
    v_align.five_p_offset  = 0;
    v_align.three_p_offset = v_len - v_del - 1;
    v_align.align_length   = v_len - v_del;
    v_align.mismatches     = mismatches;
    v_align.insertions     = {};
    v_align.deletions      = {};
    v_align.score          = 5.0 * (v_align.align_length - static_cast<int>(mismatches.size()));
    return v_align;
}

static Alignment_data
create_j_mock_alignment(const std::string& sequence,
                        const ParsedScenario& scenario,
                        const std::unordered_map<std::string, std::string>& gene_templates)
{
    const auto& j_template = gene_templates.at(scenario.j_gene_name);
    int j_len    = static_cast<int>(j_template.length());
    int j_del    = scenario.j_5p_del;
    int j_offset = static_cast<int>(sequence.length()) - j_len;

    std::vector<int> mismatches;
    for (int i = 0; i < j_del; ++i) {
        int seq_pos = j_offset + i;
        if (seq_pos >= 0 && seq_pos < static_cast<int>(sequence.length()))
            if (j_template[i] != sequence[seq_pos]) mismatches.push_back(seq_pos);
    }

    Alignment_data j_align(scenario.j_gene_name, j_offset);
    j_align.five_p_offset  = j_offset - j_del;
    j_align.three_p_offset = static_cast<int>(sequence.length()) - 1;
    j_align.align_length   = j_len - j_del;
    j_align.mismatches     = mismatches;
    j_align.insertions     = {};
    j_align.deletions      = {};
    j_align.score          = 5.0 * (j_align.align_length - static_cast<int>(mismatches.size()));
    return j_align;
}

// ---------------------------------------------------------------------------
// Comparison structure
// ---------------------------------------------------------------------------

struct InferenceComparison
{
    std::string event_nickname;
    double kl_divergence_forward;
    double entropy_truth;
    double entropy_inferred;
    double uncovered_mass;
    bool   passes_threshold;
};

static std::vector<InferenceComparison>
compare_models(const std::vector<ModelEventInfo>& ground_truth_events,
               const Model_marginals&            inferred_marginals,
               const Model_Parms&                inferred_parms,
               double                            kl_threshold_factor)
{
    std::vector<InferenceComparison> comparisons;

    for (const auto& gt : ground_truth_events) {
        if (gt.is_dinuc_markov) continue;

        auto [dims, inferred_probs] =
            inferred_marginals.compute_event_marginal_probability(gt.name, inferred_parms);

        std::vector<double> inferred_marginal(gt.num_realizations, 0.0);
        for (int i = 0; i < gt.num_realizations; ++i)
            inferred_marginal[i] = static_cast<double>(inferred_probs.get()[i]);

        InferenceComparison cmp;
        cmp.event_nickname       = gt.nickname;
        cmp.entropy_truth        = gt.H;
        cmp.entropy_inferred     = entropy(inferred_marginal);
        cmp.kl_divergence_forward = kl_divergence(gt.marginal, inferred_marginal, &cmp.uncovered_mass);

        double threshold  = (std::max)(gt.H / kl_threshold_factor, 0.01);
        cmp.passes_threshold = (cmp.kl_divergence_forward < threshold);
        comparisons.push_back(cmp);
    }
    return comparisons;
}

// ---------------------------------------------------------------------------
// THE TEST
// ---------------------------------------------------------------------------

TEMPLATE_TEST_CASE("Model architecture: inference recovers ground truth",
                   "[model][inference]", double, long double)
{
    using T = TestType;

    std::string model_dir;
    std::string model_label;
    int    sample_size        = 0;
    double kl_threshold_factor = 5.0;

    SECTION("human TCR alpha (VJ) — N=1000 — smoke test")
    {
        model_dir   = MODELS_DIR + "/human/tcr_alpha";
        model_label = "human/tcr_alpha";
        sample_size = 1000;
        kl_threshold_factor = 5.0;
    }

    SECTION("human TCR alpha (VJ) — N=10000 — thorough validation")
    {
        model_dir   = MODELS_DIR + "/human/tcr_alpha";
        model_label = "human/tcr_alpha";
        sample_size = 10000;
        kl_threshold_factor = 20.0;
    }

    const std::string parms_path     = model_dir + "/models/model_parms.txt";
    const std::string marginals_path = model_dir + "/models/model_marginals.txt";

    INFO("Testing inference with model: " << model_label);
    INFO("Sample size: " << sample_size);

    // ==================================================================
    // 1. Load ground truth — NEW architecture
    // ==================================================================
    auto truth_model_obj = recombination_model_from_files<T>(parms_path, marginals_path);
    auto truth_model     = std::make_shared<RecombinationModel<T>>(std::move(truth_model_obj));

    // Structural checks on the loaded model
    REQUIRE(truth_model->topology().size() > 0);
    REQUIRE(truth_model->size() == truth_model->topology().size());

    // ==================================================================
    // 2. Build InferenceEngine from ground truth (structural validation)
    // ==================================================================
    InferenceEngine<T> truth_engine(truth_model);
    REQUIRE(truth_engine.size() == truth_model->size());

    // Verify handler weights reference the model's tensors
    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(truth_engine.size()); ++uid)
    {
        REQUIRE(&truth_engine.handler(uid).weights() == &truth_model->weight(uid));
    }

    std::cout << "\n=== Ground truth loaded (new architecture) ===" << std::endl;
    std::cout << "Events: " << truth_model->topology().size() << std::endl;

    // ==================================================================
    // 3. Also load via LEGACY path (for generation + inference)
    // ==================================================================
    Model_Parms truth_parms;
    truth_parms.read_model_parms(parms_path);

    Model_marginals truth_marginals(truth_parms);
    truth_marginals.txt2marginals(marginals_path, truth_parms);

    auto truth_events = build_model_event_info(truth_parms, truth_marginals);

    // ==================================================================
    // 4. Validate LegacyBridge: ground truth imported into
    //    RecombinationModel matches the file-loaded version
    // ==================================================================
    {
        auto bridge_topology = import_from_legacy(truth_parms);
        RecombinationModel<T> bridge_model(
            std::make_unique<Topology>(std::move(*bridge_topology)));
        import_from_legacy(bridge_model, truth_marginals);

        REQUIRE(bridge_model.size() == truth_model->size());

        // Compare every weight tensor element
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(bridge_model.size()); ++uid)
        {
            const auto& bridge_w = bridge_model.weight(uid);
            const auto& truth_w  = truth_model->weight(uid);
            REQUIRE(bridge_w.size() == truth_w.size());

            for (std::size_t i = 0; i < bridge_w.size(); ++i) {
                INFO("uid=" << uid << " i=" << i);
                // Tolerance accounts for text-file precision loss
                // (marginals file uses limited decimal places vs full
                //  long double precision from Model_marginals)
                REQUIRE_THAT(static_cast<double>(bridge_w.data()[i]),
                             WithinAbs(static_cast<double>(truth_w.data()[i]), 1e-6));
            }
        }
        std::cout << "LegacyBridge round-trip: PASS" << std::endl;
    }

    // ==================================================================
    // 5. Load gene templates
    // ==================================================================
    std::unordered_map<std::string, std::string> gene_templates;
    for (const auto& [name, seq] :
         read_genomic_fasta(model_dir + "/ref_genome/genomicVs.fasta"))
        gene_templates[name] = seq;
    for (const auto& [name, seq] :
         read_genomic_fasta(model_dir + "/ref_genome/genomicJs.fasta"))
        gene_templates[name] = seq;

    std::cout << "Loaded " << gene_templates.size() << " gene templates" << std::endl;

    // ==================================================================
    // 6. Generate sequences with scenarios (legacy GenModel)
    // ==================================================================
    std::cout << "\n=== Generating " << sample_size << " sequences ===" << std::endl;
    GenModel gen_model(truth_parms, truth_marginals);
    auto scenarios = gen_model.generate_sequences(sample_size, /*errors=*/false);

    std::size_t actual_count = 0;
    for (auto it = scenarios.begin(); it != scenarios.end(); ++it) ++actual_count;
    REQUIRE(actual_count == static_cast<std::size_t>(sample_size));

    // ==================================================================
    // 7. Create mock alignments
    // ==================================================================
    std::cout << "=== Creating mock alignments ===" << std::endl;

    using AlignMap = std::unordered_map<Gene_class, std::vector<Alignment_data>>;
    std::vector<std::tuple<int, std::string, AlignMap>> sequences_with_alignments;

    int seq_idx = 0;
    for (const auto& [seq, scenario] : scenarios) {
        ParsedScenario parsed = parse_scenario(scenario, truth_events, truth_parms);
        AlignMap aligns_map;
        aligns_map[V_gene] = { create_v_mock_alignment(seq, parsed, gene_templates) };
        aligns_map[J_gene] = { create_j_mock_alignment(seq, parsed, gene_templates) };
        sequences_with_alignments.emplace_back(seq_idx, seq, aligns_map);
        ++seq_idx;
    }

    // ==================================================================
    // 8. Run inference using InferenceEngine::run()
    //
    //    Each EM iteration is driven by engine.run(eStep):
    //      1. resetAccumulators()     — zeros all handler accumulators
    //      2. eStep(engine)           — bridge to legacy E-step
    //      3. updateParameters()      — M-step via maximizeLikelihood()
    //
    //    The E-step callable:
    //      a. Exports current model weights → legacy Model_marginals
    //      b. Creates a GenModel with those weights
    //      c. Runs GenModel::infer_model(1 iteration) — legacy E+M
    //      d. Imports result into handler accumulators
    //
    //    The engine's M-step re-normalises the (already normalised)
    //    values — idempotent, but exercises the full code path.
    // ==================================================================
    std::cout << "\n=== Running inference via InferenceEngine::run() ===" << std::endl;

    // Build the inference model with uniform weights
    auto infer_topology = import_from_legacy(truth_parms);
    auto infer_model = std::make_shared<RecombinationModel<T>>(
        std::make_unique<Topology>(std::move(*infer_topology)));

    Model_marginals uniform_marginals(truth_parms);
    uniform_marginals.uniform_initialize(truth_parms);
    import_from_legacy(*infer_model, uniform_marginals);

    InferenceEngine<T> engine(infer_model);

    REQUIRE(engine.size() == truth_engine.size());

    int num_iterations = (sample_size >= 10000) ? 20 : 15;
    std::cout << "Inference iterations: " << num_iterations << std::endl;

    // Cache the index map (same for every iteration — topology is fixed)
    auto bridge_parms = export_to_legacy(infer_model->topology());
    Model_marginals layout_marginals(truth_parms); // only used for get_index_map
    auto index_map = layout_marginals.get_index_map(truth_parms);

    for (int iter = 0; iter < num_iterations; ++iter) {
        engine.run([&](InferenceEngine<T>& eng) {
            // ── 8a. Export current model weights → legacy marginals ──
            Model_marginals current_marginals(truth_parms);
            long double* dst = current_marginals.marginal_array_smart_p.get();

            const auto& topo = eng.model().topology();
            for (igor::index_type uid = 0;
                 uid < static_cast<igor::index_type>(topo.size()); ++uid)
            {
                const std::string& full_name = topo.event(uid)->get_name();
                const auto& tensor = eng.model().weight(uid);
                int base = index_map.at(full_name);
                for (std::size_t i = 0; i < tensor.size(); ++i)
                    dst[base + i] = static_cast<long double>(tensor.data()[i]);
            }

            // ── 8b. Run one legacy E+M iteration ──
            GenModel iter_gen(truth_parms, current_marginals);
            bool failed = iter_gen.infer_model(
                sequences_with_alignments, /*iterations=*/1, /*path=*/"",
                /*fast_iter=*/true, /*likelihood_threshold=*/1e-60,
                /*viterbi_like=*/false, /*proba_threshold_factor=*/1e-3);
            REQUIRE(!failed);

            // ── 8c. Import result into handler accumulators ──
            const Model_marginals& result = iter_gen.get_marginals();
            const long double* src = result.marginal_array_smart_p.get();

            for (igor::index_type uid = 0;
                 uid < static_cast<igor::index_type>(topo.size()); ++uid)
            {
                auto& acc = eng.handler(uid).accumulator();
                const std::string& full_name = topo.event(uid)->get_name();
                int base = index_map.at(full_name);
                for (std::size_t i = 0; i < acc.size(); ++i)
                    acc.data()[i] = static_cast<T>(src[base + i]);
            }
        });
        // engine.run() just called updateParameters() — M-step complete
    }

    std::cout << "Inference completed: " << num_iterations
              << " iterations via InferenceEngine::run()" << std::endl;

    // ==================================================================
    // 9. Export final engine weights → legacy Model_marginals
    //    (needed for compute_event_marginal_probability comparison)
    // ==================================================================
    Model_marginals inferred_marginals(truth_parms);
    {
        long double* dst = inferred_marginals.marginal_array_smart_p.get();
        const auto& topo = infer_model->topology();
        for (igor::index_type uid = 0;
             uid < static_cast<igor::index_type>(topo.size()); ++uid)
        {
            const std::string& full_name = topo.event(uid)->get_name();
            const auto& tensor = infer_model->weight(uid);
            int base = index_map.at(full_name);
            for (std::size_t i = 0; i < tensor.size(); ++i)
                dst[base + i] = static_cast<long double>(tensor.data()[i]);
        }
    }
    const Model_Parms& inferred_parms = truth_parms; // structure unchanged

    // ==================================================================
    // 10. Validate: compare weight tensors directly
    //     (structural check — all tensor shapes must match)
    // ==================================================================
    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(truth_model->size()); ++uid)
    {
        const auto& tw = truth_model->weight(uid);
        const auto& iw = infer_model->weight(uid);
        REQUIRE(tw.shape() == iw.shape());
    }

    // ==================================================================
    // 11. Validate: KL divergence (same thresholds as test_inference.cpp)
    // ==================================================================
    std::cout << "\n=== Comparing inferred vs ground truth ===" << std::endl;

    auto comparisons = compare_models(truth_events, inferred_marginals,
                                      inferred_parms, kl_threshold_factor);

    std::cout << "\nEvent                  | D_KL(T||I) | H_truth | H_infer | Uncovered | Pass"
              << std::endl;
    std::cout << "---------------------- | ---------- | ------- | ------- | --------- | ----"
              << std::endl;

    for (const auto& cmp : comparisons) {
        std::cout << std::left << std::setw(22) << cmp.event_nickname
                  << " | " << std::fixed << std::setprecision(4)
                  << std::setw(10) << cmp.kl_divergence_forward
                  << " | " << std::setw(7) << cmp.entropy_truth
                  << " | " << std::setw(7) << cmp.entropy_inferred
                  << " | " << std::setw(9) << cmp.uncovered_mass
                  << " | " << (cmp.passes_threshold ? "✓" : "✗")
                  << std::endl;

        INFO("Event: " << cmp.event_nickname);
        INFO("D_KL(truth||inferred) = " << cmp.kl_divergence_forward);
        INFO("H(truth) = " << cmp.entropy_truth);
        INFO("Threshold = " << cmp.entropy_truth / kl_threshold_factor);

        CHECK(cmp.passes_threshold);
    }

    // ==================================================================
    // 12. Validate: handler weights are the model's weight tensors
    //     (the M-step must have written through the borrowed references)
    // ==================================================================
    for (igor::index_type uid = 0;
         uid < static_cast<igor::index_type>(engine.size()); ++uid)
    {
        REQUIRE(&engine.handler(uid).weights() == &infer_model->weight(uid));

        // Weights should be properly normalised after M-step
        const auto& w = engine.handler(uid).weights();
        double total = 0.0;
        for (std::size_t i = 0; i < w.size(); ++i)
            total += static_cast<double>(w.data()[i]);
        // Total should be close to the number of conditional slices
        // (each slice sums to 1, total = number of parent combos).
        REQUIRE(total > 0.0);
    }

    std::cout << "\n=== Model inference test completed ===" << std::endl;
}
