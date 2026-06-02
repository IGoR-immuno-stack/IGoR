#include "sha256.h"

#include <CLI/CLI.hpp>
#include <toml++/toml.h>

#include <igor/Core/Config.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/utsname.h>
#endif

int igor_legacy_main(int argc, char *argv[]);

namespace fs = std::filesystem;

namespace {

#ifndef IGOR_COMPILER_ID
#define IGOR_COMPILER_ID "unknown"
#endif
#ifndef IGOR_COMPILER_VERSION
#define IGOR_COMPILER_VERSION "unknown"
#endif
#ifndef IGOR_COMPILER_PATH
#define IGOR_COMPILER_PATH "unknown"
#endif
#ifndef IGOR_TARGET_SYSTEM_NAME
#define IGOR_TARGET_SYSTEM_NAME "unknown"
#endif
#ifndef IGOR_TARGET_SYSTEM_PROCESSOR
#define IGOR_TARGET_SYSTEM_PROCESSOR "unknown"
#endif

std::string clean_path(std::string s)
{
    auto pos0 = s.find('\0');
    if (pos0 != std::string::npos) {
        s.resize(pos0);
    }
    return s;
}

std::string data_dir()
{
    return clean_path(IGOR_DATA_DIR_STR);
}

std::string trim(std::string s)
{
    auto is_ws = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), is_ws));
    s.erase(std::find_if_not(s.rbegin(), s.rend(), is_ws).base(), s.end());
    return s;
}

std::vector<std::string> split_csv(const std::string &value)
{
    std::vector<std::string> parts;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item = trim(item);
        if (!item.empty()) {
            parts.push_back(item);
        }
    }
    return parts;
}

struct GlobalOptions
{
    fs::path workdir = fs::current_path();
    std::string batch;
    std::optional<int> threads;
    std::optional<fs::path> stdout_file;
};

enum class ModelSource {
    Custom,
    Builtin,
    LastInferred,
};

struct ModelOptions
{
    ModelSource source = ModelSource::Custom;
    std::string species;
    std::string chain;
    fs::path parms;
    fs::path marginals;
    bool load_last_inferred = false;
};

struct SequenceOptions
{
    fs::path input;
};

struct AlignGeneOptions
{
    double threshold = 0.0;
    double gap_penalty = 0.0;
    bool best_align_only = false;
    bool best_gene_only = false;
    int left_offset = 0;
    int right_offset = 0;
};

struct AlignOptions
{
    fs::path genomic_v;
    fs::path genomic_d;
    fs::path genomic_j;
    fs::path anchors_v;
    fs::path anchors_j;
    bool cdr3_input = false;
    bool extract_cdr3 = false;
    AlignGeneOptions v;
    AlignGeneOptions d;
    AlignGeneOptions j;
};

struct InferOptions
{
    int iterations = 0;
    double likelihood_threshold = 0.0;
    double probability_ratio_threshold = 0.0;
    bool viterbi = false;
    bool fix_error_rate = false;
    std::vector<std::string> only;
    std::vector<std::string> not_infer;
};

struct EvaluateOptions
{
    double likelihood_threshold = 0.0;
    double probability_ratio_threshold = 0.0;
    bool viterbi = false;
};

struct GenerateOptions
{
    int seed = -1;
    bool fast = false;
    bool error = true;
    bool cdr3 = false;
    std::string filename_prefix;
    int threads = 0;
};

struct OutputOptions
{
    bool pgen = false;
    int scenarios = 0;
    std::string coverage;
};

struct PipelineOptions
{
    std::vector<std::string> steps;
};

struct ResolvedOptions
{
    ModelOptions model;
    SequenceOptions sequences;
    AlignOptions alignment;
    InferOptions infer;
    EvaluateOptions evaluate;
    GenerateOptions generate;
    OutputOptions output;
    PipelineOptions pipeline;
};

struct ConfigKey
{
    std::string key;
    std::string type;
    std::string default_value;
    std::string description;
};

const std::vector<ConfigKey> &config_schema()
{
    static const std::vector<ConfigKey> schema = [] {
        const std::string dd = data_dir();
        return std::vector<ConfigKey>{
            { "model.source", "enum(custom,builtin,last_inferred)", "custom",
              "Model source used by infer, evaluate, and generate." },
            { "model.species", "string", "human", "Species for built-in model loading." },
            { "model.chain", "string", "beta", "Chain for built-in model loading." },
            { "model.parms", "path", dd + "/scripts/tests/data/input/TRB_model_parms.txt",
              "Model parameters file for custom model loading." },
            { "model.marginals", "path", dd + "/scripts/tests/data/input/TRB_uniform_model_marginals.txt",
              "Model marginals file for custom model loading." },
            { "model.load_last_inferred", "bool", "false", "Load the last inferred model from the workdir." },
            { "genomic.V", "path", dd + "/demo/genomicVs_with_primers.fasta", "V genomic templates FASTA." },
            { "genomic.D", "path", dd + "/demo/genomicDs.fasta", "D genomic templates FASTA." },
            { "genomic.J", "path", dd + "/demo/genomicJs_all_curated.fasta", "J genomic templates FASTA." },
            { "anchors.V", "path", "", "Optional V CDR3 anchors CSV." },
            { "anchors.J", "path", "", "Optional J CDR3 anchors CSV." },
            { "sequences.input", "path", "", "Input sequence file used by `igor run` import step." },
            { "alignment.cdr3_input", "bool", "false", "Treat input sequences as nucleotide CDR3s." },
            { "alignment.extract_cdr3", "bool", "false", "Extract CDR3 features after V/J alignment." },
            { "alignment.V.threshold", "double", "50", "V alignment score threshold." },
            { "alignment.V.gap_penalty", "double", "50", "V alignment gap penalty." },
            { "alignment.V.best_align_only", "bool", "true", "Keep only best V alignments." },
            { "alignment.V.best_gene_only", "bool", "false", "Keep only best V gene candidate." },
            { "alignment.V.left_offset", "int", "-32768", "V minimum alignment offset." },
            { "alignment.V.right_offset", "int", "32767", "V maximum alignment offset." },
            { "alignment.D.threshold", "double", "15", "D alignment score threshold." },
            { "alignment.D.gap_penalty", "double", "50", "D alignment gap penalty." },
            { "alignment.D.best_align_only", "bool", "false", "Keep only best D alignments." },
            { "alignment.D.best_gene_only", "bool", "false", "Keep only best D gene candidate." },
            { "alignment.D.left_offset", "int", "-32768", "D minimum alignment offset." },
            { "alignment.D.right_offset", "int", "32767", "D maximum alignment offset." },
            { "alignment.J.threshold", "double", "15", "J alignment score threshold." },
            { "alignment.J.gap_penalty", "double", "50", "J alignment gap penalty." },
            { "alignment.J.best_align_only", "bool", "true", "Keep only best J alignments." },
            { "alignment.J.best_gene_only", "bool", "false", "Keep only best J gene candidate." },
            { "alignment.J.left_offset", "int", "-32768", "J minimum alignment offset." },
            { "alignment.J.right_offset", "int", "32767", "J maximum alignment offset." },
            { "infer.iterations", "uint", "5", "Number of EM inference iterations." },
            { "infer.likelihood_threshold", "double", "1e-60", "Inference likelihood threshold." },
            { "infer.probability_ratio_threshold", "double", "1e-5", "Inference probability ratio threshold." },
            { "infer.viterbi", "bool", "false", "Use most-likely scenario only during inference." },
            { "infer.fix_error_rate", "bool", "false", "Keep error rate fixed during inference." },
            { "infer.only", "csv", "", "Comma-separated event nicknames to infer exclusively." },
            { "infer.not", "csv", "", "Comma-separated event nicknames to keep fixed." },
            { "evaluate.likelihood_threshold", "double", "1e-60", "Evaluation likelihood threshold." },
            { "evaluate.probability_ratio_threshold", "double", "1e-5", "Evaluation probability ratio threshold." },
            { "evaluate.viterbi", "bool", "false", "Use most-likely scenario only during evaluation." },
            { "generate.seed", "int", "-1", "Generation seed; -1 means time-based random seed." },
            { "generate.fast", "bool", "false", "Use fast parallel generation." },
            { "generate.error", "bool", "true", "Generate sequences with errors." },
            { "generate.cdr3", "bool", "false", "Write generated CDR3 information when supported." },
            { "generate.filename_prefix", "string", "", "Prefix for generated sequence output filenames." },
            { "generate.threads", "uint", "0", "Generation thread count; 0 means auto." },
            { "output.Pgen", "bool", "false", "Enable Pgen counter output." },
            { "output.scenarios", "uint", "0", "Number of best scenarios to record; 0 disables it." },
            { "output.coverage", "string", "", "Coverage counter target such as VJ_gene." },
            { "pipeline.steps", "csv", "import-seqs,align,infer", "Comma-separated steps for `igor run`." },
        };
    }();
    return schema;
}

const ConfigKey *find_config_key(const std::string &key)
{
    for (const auto &item : config_schema()) {
        if (item.key == key) {
            return &item;
        }
    }
    return nullptr;
}

bool parse_bool_text(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); });
    if (value == "true" || value == "1" || value == "yes") {
        return true;
    }
    if (value == "false" || value == "0" || value == "no") {
        return false;
    }
    throw std::runtime_error("expected boolean value, received: " + value);
}

void validate_config_value(const ConfigKey &key, const std::string &value)
{
    try {
        if (key.type == "bool") {
            (void)parse_bool_text(value);
        } else if (key.type == "int") {
            std::size_t pos = 0;
            (void)std::stoi(value, &pos);
            if (pos != value.size()) {
                throw std::invalid_argument("trailing characters");
            }
        } else if (key.type == "uint") {
            std::size_t pos = 0;
            int parsed = std::stoi(value, &pos);
            if (pos != value.size() || parsed < 0) {
                throw std::invalid_argument("not unsigned");
            }
        } else if (key.type == "double") {
            std::size_t pos = 0;
            (void)std::stod(value, &pos);
            if (pos != value.size()) {
                throw std::invalid_argument("trailing characters");
            }
        } else if (key.key == "model.source") {
            if (value != "custom" && value != "builtin" && value != "last_inferred") {
                throw std::invalid_argument("unsupported model source");
            }
        }
    } catch (const std::exception &) {
        throw std::runtime_error("invalid value for config key `" + key.key + "` (" + key.type + "): " + value);
    }
}

struct Config
{
    std::map<std::string, std::string> values;

    static Config defaults()
    {
        Config c;
        for (const auto &item : config_schema()) {
            c.values.emplace(item.key, item.default_value);
        }
        return c;
    }

    std::string get(const std::string &key) const
    {
        auto it = values.find(key);
        if (it == values.end()) {
            throw std::runtime_error("unknown config key: " + key);
        }
        return it->second;
    }

    bool get_bool(const std::string &key) const { return parse_bool_text(get(key)); }

    void set(const std::string &key, const std::string &value)
    {
        const ConfigKey *schema_key = find_config_key(key);
        if (schema_key == nullptr || !values.contains(key)) {
            throw std::runtime_error("unknown config key: " + key);
        }
        validate_config_value(*schema_key, value);
        values[key] = value;
    }
};

ModelSource parse_model_source(const std::string &value)
{
    if (value == "custom") {
        return ModelSource::Custom;
    }
    if (value == "builtin") {
        return ModelSource::Builtin;
    }
    if (value == "last_inferred") {
        return ModelSource::LastInferred;
    }
    throw std::runtime_error("unsupported model source: " + value);
}

std::string model_source_name(ModelSource source)
{
    switch (source) {
    case ModelSource::Custom:
        return "custom";
    case ModelSource::Builtin:
        return "builtin";
    case ModelSource::LastInferred:
        return "last_inferred";
    }
    throw std::runtime_error("unsupported model source");
}

int config_int(const Config &config, const std::string &key)
{
    return std::stoi(config.get(key));
}

double config_double(const Config &config, const std::string &key)
{
    return std::stod(config.get(key));
}

AlignGeneOptions resolve_align_gene(const Config &config, const std::string &gene)
{
    const std::string prefix = "alignment." + gene + ".";
    return AlignGeneOptions{
        config_double(config, prefix + "threshold"), config_double(config, prefix + "gap_penalty"),
        config.get_bool(prefix + "best_align_only"), config.get_bool(prefix + "best_gene_only"),
        config_int(config, prefix + "left_offset"),  config_int(config, prefix + "right_offset"),
    };
}

ResolvedOptions resolve_config(const Config &config)
{
    ResolvedOptions options;
    options.model = ModelOptions{
        parse_model_source(config.get("model.source")),
        config.get("model.species"),
        config.get("model.chain"),
        config.get("model.parms"),
        config.get("model.marginals"),
        config.get_bool("model.load_last_inferred"),
    };
    options.sequences = SequenceOptions{ config.get("sequences.input") };
    options.alignment = AlignOptions{
        config.get("genomic.V"),
        config.get("genomic.D"),
        config.get("genomic.J"),
        config.get("anchors.V"),
        config.get("anchors.J"),
        config.get_bool("alignment.cdr3_input"),
        config.get_bool("alignment.extract_cdr3"),
        resolve_align_gene(config, "V"),
        resolve_align_gene(config, "D"),
        resolve_align_gene(config, "J"),
    };
    options.infer = InferOptions{
        config_int(config, "infer.iterations"),
        config_double(config, "infer.likelihood_threshold"),
        config_double(config, "infer.probability_ratio_threshold"),
        config.get_bool("infer.viterbi"),
        config.get_bool("infer.fix_error_rate"),
        split_csv(config.get("infer.only")),
        split_csv(config.get("infer.not")),
    };
    options.evaluate = EvaluateOptions{
        config_double(config, "evaluate.likelihood_threshold"),
        config_double(config, "evaluate.probability_ratio_threshold"),
        config.get_bool("evaluate.viterbi"),
    };
    options.generate = GenerateOptions{
        config_int(config, "generate.seed"),    config.get_bool("generate.fast"),
        config.get_bool("generate.error"),      config.get_bool("generate.cdr3"),
        config.get("generate.filename_prefix"), config_int(config, "generate.threads"),
    };
    options.output = OutputOptions{
        config.get_bool("output.Pgen"),
        config_int(config, "output.scenarios"),
        config.get("output.coverage"),
    };
    options.pipeline = PipelineOptions{ split_csv(config.get("pipeline.steps")) };
    return options;
}

fs::path config_path(const fs::path &workdir)
{
    return workdir / ".igor" / "config.toml";
}

std::vector<std::string> split_key(const std::string &key)
{
    std::vector<std::string> parts;
    std::stringstream ss(key);
    std::string part;
    while (std::getline(ss, part, '.')) {
        if (!part.empty()) {
            parts.push_back(part);
        }
    }
    return parts;
}

toml::node_view<const toml::node> toml_at_path(const toml::table &table, const std::string &key)
{
    toml::node_view<const toml::node> node{ table };
    for (const auto &part : split_key(key)) {
        node = node[part];
        if (!node) {
            return {};
        }
    }
    return node;
}

std::string number_text(double value)
{
    std::ostringstream out;
    out << std::setprecision(15) << value;
    return out.str();
}

std::string toml_quote(const std::string &value)
{
    std::string quoted = "\"";
    for (char c : value) {
        switch (c) {
        case '\\':
            quoted += "\\\\";
            break;
        case '"':
            quoted += "\\\"";
            break;
        case '\n':
            quoted += "\\n";
            break;
        case '\r':
            quoted += "\\r";
            break;
        case '\t':
            quoted += "\\t";
            break;
        default:
            quoted += c;
            break;
        }
    }
    quoted += "\"";
    return quoted;
}

std::string toml_config_value_text(const std::string &type, const std::string &value)
{
    if (type == "bool") {
        return parse_bool_text(value) ? "true" : "false";
    }
    if (type == "int" || type == "uint") {
        return std::to_string(std::stoll(value));
    }
    if (type == "double") {
        validate_config_value(ConfigKey{ "", "double", "", "" }, value);
        return value;
    }
    return toml_quote(value);
}

std::string toml_value_to_string(toml::node_view<const toml::node> node, const std::string &key)
{
    const toml::node *raw = node.node();
    if (raw == nullptr) {
        throw std::runtime_error("missing TOML value for config key: " + key);
    }
    if (raw->is_string()) {
        return *raw->value<std::string>();
    }
    if (raw->is_boolean()) {
        return *raw->value<bool>() ? "true" : "false";
    }
    if (raw->is_integer()) {
        return std::to_string(*raw->value<int64_t>());
    }
    if (raw->is_floating_point()) {
        return number_text(*raw->value<double>());
    }
    throw std::runtime_error("unsupported TOML value type for config key: " + key);
}

void collect_toml_leaf_keys(const toml::table &table, const std::string &prefix, std::set<std::string> &keys)
{
    for (const auto &[raw_key, node] : table) {
        const std::string key = prefix.empty() ? std::string(raw_key.str()) : prefix + "." + std::string(raw_key.str());
        if (node.is_table()) {
            collect_toml_leaf_keys(*node.as_table(), key, keys);
        } else {
            keys.insert(key);
        }
    }
}

std::map<std::string, std::string> raw_toml_leaf_values(const fs::path &path)
{
    std::map<std::string, std::string> values;
    std::ifstream in(path);
    if (!in) {
        return values;
    }

    std::vector<std::string> table_path;
    std::string line;
    while (std::getline(in, line)) {
        std::string cleaned;
        bool in_string = false;
        bool escaped = false;
        for (char c : line) {
            if (!in_string && c == '#') {
                break;
            }
            cleaned += c;
            if (escaped) {
                escaped = false;
            } else if (c == '\\' && in_string) {
                escaped = true;
            } else if (c == '"') {
                in_string = !in_string;
            }
        }
        cleaned = trim(cleaned);
        if (cleaned.empty()) {
            continue;
        }
        if (cleaned.front() == '[' && cleaned.back() == ']') {
            table_path = split_key(trim(cleaned.substr(1, cleaned.size() - 2)));
            continue;
        }
        const auto equals = cleaned.find('=');
        if (equals == std::string::npos) {
            continue;
        }
        std::string key = trim(cleaned.substr(0, equals));
        std::string value = trim(cleaned.substr(equals + 1));
        std::string full_key;
        for (const auto &part : table_path) {
            if (!full_key.empty()) {
                full_key += ".";
            }
            full_key += part;
        }
        if (!full_key.empty()) {
            full_key += ".";
        }
        full_key += key;
        values[full_key] = value;
    }
    return values;
}

Config read_config(const fs::path &workdir, bool require = true)
{
    Config config = Config::defaults();
    const fs::path path = config_path(workdir);
    if (!fs::exists(path)) {
        if (require) {
            throw std::runtime_error("missing workdir config: " + path.string() + " (run `igor init -w DIR`)");
        }
        return config;
    }

    try {
        const auto raw_values = raw_toml_leaf_values(path);
        const toml::table parsed = toml::parse_file(path.string());
        std::set<std::string> parsed_keys;
        collect_toml_leaf_keys(parsed, "", parsed_keys);
        for (const auto &key : parsed_keys) {
            if (!config.values.contains(key)) {
                throw std::runtime_error("unknown config key in " + path.string() + ": " + key);
            }
        }
        for (const auto &key : parsed_keys) {
            const ConfigKey *schema_key = find_config_key(key);
            const auto raw = raw_values.find(key);
            if (schema_key != nullptr && schema_key->type == "double" && raw != raw_values.end() &&
                !raw->second.empty() && raw->second.front() != '"') {
                config.set(key, raw->second);
            } else {
                config.set(key, toml_value_to_string(toml_at_path(parsed, key), key));
            }
        }
    } catch (const toml::parse_error &e) {
        throw std::runtime_error("invalid TOML config " + path.string() + ": " + std::string(e.description()));
    }
    return config;
}

void write_config(const fs::path &workdir, const Config &config)
{
    fs::create_directories(workdir / ".igor");
    std::ofstream out(config_path(workdir));
    if (!out) {
        throw std::runtime_error("cannot write config: " + config_path(workdir).string());
    }

    std::vector<std::string> current_table;
    for (const auto &item : config_schema()) {
        const auto parts = split_key(item.key);
        const std::vector<std::string> table_parts(parts.begin(), parts.end() - 1);
        if (table_parts != current_table) {
            if (!current_table.empty() || out.tellp() > 0) {
                out << "\n";
            }
            if (!table_parts.empty()) {
                out << "[";
                for (std::size_t i = 0; i < table_parts.size(); ++i) {
                    if (i) {
                        out << ".";
                    }
                    out << table_parts[i];
                }
                out << "]\n";
            }
            current_table = table_parts;
        }
        out << parts.back() << " = " << toml_config_value_text(item.type, config.get(item.key)) << "\n";
    }
}

void ensure_workdir(const fs::path &workdir)
{
    fs::create_directories(workdir / ".igor" / "runs");
    fs::create_directories(workdir / "aligns");
    fs::create_directories(workdir / "output");
    fs::create_directories(workdir / "inference");
    fs::create_directories(workdir / "evaluate");
    fs::create_directories(workdir / "generated");
}

std::string batch_arg(const GlobalOptions &global)
{
    return global.batch.empty() ? "" : global.batch;
}

std::vector<std::string> base_legacy_args(const GlobalOptions &global, const ResolvedOptions &options,
                                          bool include_model, bool include_genomics)
{
    std::vector<std::string> args = { "igor", "-set_wd", global.workdir.string() };
    if (!global.batch.empty()) {
        args.insert(args.end(), { "-batch", global.batch });
    }
    if (global.threads) {
        args.insert(args.end(), { "-threads", std::to_string(*global.threads) });
    }
    if (global.stdout_file) {
        args.insert(args.end(), { "-stdout_f", global.stdout_file->string() });
    }

    const auto add_if = [&](const std::string &flag, const std::string &gene, const std::string &path,
                            std::vector<std::string> &target) {
        if (!path.empty()) {
            target.insert(target.end(), { flag, gene, path });
        }
    };
    if (include_genomics) {
        std::vector<std::string> genomic;
        add_if("-set_genomic", "--V", options.alignment.genomic_v.string(), genomic);
        add_if("-set_genomic", "--D", options.alignment.genomic_d.string(), genomic);
        add_if("-set_genomic", "--J", options.alignment.genomic_j.string(), genomic);
        args.insert(args.end(), genomic.begin(), genomic.end());

        std::vector<std::string> anchors;
        add_if("-set_CDR3_anchors", "--V", options.alignment.anchors_v.string(), anchors);
        add_if("-set_CDR3_anchors", "--J", options.alignment.anchors_j.string(), anchors);
        args.insert(args.end(), anchors.begin(), anchors.end());
    }

    if (include_model) {
        if (options.model.load_last_inferred || options.model.source == ModelSource::LastInferred) {
            args.push_back("-load_last_inferred");
        } else if (options.model.source == ModelSource::Custom) {
            args.insert(args.end(), { "-set_custom_model", options.model.parms.string() });
            if (!options.model.marginals.empty()) {
                args.push_back(options.model.marginals.string());
            }
        } else {
            args.insert(args.end(), { "-species", options.model.species, "-chain", options.model.chain });
        }
    } else if (options.model.source == ModelSource::Builtin) {
        args.insert(args.end(), { "-species", options.model.species, "-chain", options.model.chain });
    }

    if (options.output.pgen) {
        args.insert(args.end(), { "-output", "--Pgen" });
    }
    if (options.output.scenarios > 0) {
        args.insert(args.end(), { "-output", "--scenarios", std::to_string(options.output.scenarios) });
    }
    if (!options.output.coverage.empty()) {
        args.insert(args.end(), { "-output", "--coverage", options.output.coverage });
    }
    return args;
}

int run_legacy(std::vector<std::string> args)
{
    std::vector<char *> raw;
    raw.reserve(args.size());
    for (auto &arg : args) {
        raw.push_back(arg.data());
    }
    return igor_legacy_main(static_cast<int>(raw.size()), raw.data());
}

void append_gene_align_args(std::vector<std::string> &args, const AlignGeneOptions &options, const std::string &gene)
{
    args.push_back("--" + gene);
    args.insert(args.end(), { "---thresh", number_text(options.threshold) });
    args.insert(args.end(), { "---gap_penalty", number_text(options.gap_penalty) });
    args.insert(args.end(), { "---best_align_only", options.best_align_only ? "true" : "false" });
    args.insert(args.end(), { "---best_gene_only", options.best_gene_only ? "true" : "false" });
    args.insert(args.end(),
                { "---offset_bounds", std::to_string(options.left_offset), std::to_string(options.right_offset) });
}

std::vector<std::string> command_args(const std::string &command, const GlobalOptions &global,
                                      const ResolvedOptions &options, const std::vector<std::string> &positional)
{
    if (command == "import-seqs") {
        auto args = base_legacy_args(global, options, false, false);
        args.insert(args.end(), { "-read_seqs", positional.at(0) });
        return args;
    }
    if (command == "align") {
        auto args = base_legacy_args(global, options, false, true);
        args.push_back("-align");
        const std::string gene = positional.at(0);
        if (gene == "all") {
            append_gene_align_args(args, options.alignment.v, "V");
            append_gene_align_args(args, options.alignment.d, "D");
            append_gene_align_args(args, options.alignment.j, "J");
        } else if (gene == "V") {
            append_gene_align_args(args, options.alignment.v, gene);
        } else if (gene == "D") {
            append_gene_align_args(args, options.alignment.d, gene);
        } else {
            append_gene_align_args(args, options.alignment.j, gene);
        }
        if (options.alignment.cdr3_input) {
            args.push_back("--ntCDR3");
        }
        if (options.alignment.extract_cdr3) {
            args.insert(args.end(), { "--feature", "---ntCDR3" });
        }
        return args;
    }
    if (command == "infer") {
        auto args = base_legacy_args(global, options, true, false);
        args.push_back("-infer");
        args.insert(args.end(), { "--L_thresh", number_text(options.infer.likelihood_threshold) });
        args.insert(args.end(), { "--P_ratio_thresh", number_text(options.infer.probability_ratio_threshold) });
        args.insert(args.end(), { "--N_iter", std::to_string(options.infer.iterations) });
        if (options.infer.viterbi) {
            args.push_back("--MLSO");
        }
        if (options.infer.fix_error_rate) {
            args.push_back("--fix_err");
        }
        if (!options.infer.only.empty()) {
            args.push_back("--infer_only");
            for (const auto &item : options.infer.only) {
                args.push_back(item);
            }
        }
        if (!options.infer.not_infer.empty()) {
            args.push_back("--not_infer");
            for (const auto &item : options.infer.not_infer) {
                args.push_back(item);
            }
        }
        return args;
    }
    if (command == "evaluate") {
        auto args = base_legacy_args(global, options, true, false);
        args.push_back("-evaluate");
        args.insert(args.end(), { "--L_thresh", number_text(options.evaluate.likelihood_threshold) });
        args.insert(args.end(), { "--P_ratio_thresh", number_text(options.evaluate.probability_ratio_threshold) });
        if (options.evaluate.viterbi) {
            args.push_back("--MLSO");
        }
        return args;
    }
    if (command == "generate") {
        auto args = base_legacy_args(global, options, true, false);
        args.insert(args.end(), { "-generate", positional.at(0) });
        if (!options.generate.error) {
            args.push_back("--noerr");
        }
        if (options.generate.cdr3) {
            args.push_back("--CDR3");
        }
        if (!options.generate.filename_prefix.empty()) {
            args.insert(args.end(), { "--name", options.generate.filename_prefix });
        }
        if (options.generate.seed != -1) {
            args.insert(args.end(), { "--seed", std::to_string(options.generate.seed) });
        }
        if (options.generate.fast) {
            args.push_back("--fast");
        }
        if (global.threads) {
            args.insert(args.end(), { "--threads", std::to_string(*global.threads) });
        } else {
            args.insert(args.end(), { "--threads", std::to_string(options.generate.threads) });
        }
        return args;
    }
    throw std::runtime_error("unsupported command: " + command);
}

std::string timestamp()
{
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000;
    std::ostringstream out;
    out << std::put_time(std::localtime(&t), "%Y%m%d-%H%M%S") << "-" << std::setw(3) << std::setfill('0') << ms;
    return out.str();
}

std::string batch_prefix(const GlobalOptions &global)
{
    if (global.batch.empty()) {
        return "";
    }
    return global.batch.back() == '_' ? global.batch : global.batch + "_";
}

void add_existing_input(std::vector<fs::path> &inputs, const fs::path &path)
{
    if (!path.empty() && fs::exists(path)) {
        inputs.push_back(fs::absolute(path));
    }
}

bool pipeline_has_step(const ResolvedOptions &options, const std::string &step)
{
    return std::find(options.pipeline.steps.begin(), options.pipeline.steps.end(), step) != options.pipeline.steps.end();
}

std::vector<fs::path> manifest_inputs(const std::string &command, const GlobalOptions &global,
                                      const ResolvedOptions &options, const std::vector<std::string> &positional)
{
    std::vector<fs::path> inputs;
    if (command == "import-seqs") {
        add_existing_input(inputs, positional.at(0));
    }
    if (command == "align") {
        add_existing_input(inputs, options.alignment.genomic_v);
        add_existing_input(inputs, options.alignment.genomic_d);
        add_existing_input(inputs, options.alignment.genomic_j);
        add_existing_input(inputs, options.alignment.anchors_v);
        add_existing_input(inputs, options.alignment.anchors_j);
    }
    if (command == "infer" || command == "evaluate" || command == "generate") {
        add_existing_input(inputs, options.model.parms);
        add_existing_input(inputs, options.model.marginals);
    }
    if (command == "align" || command == "infer" || command == "evaluate") {
        add_existing_input(inputs, global.workdir / "aligns" / (batch_prefix(global) + "indexed_sequences.csv"));
    }
    if (command == "infer" || command == "evaluate") {
        add_existing_input(inputs, global.workdir / "aligns" / (batch_prefix(global) + "V_alignments.csv"));
        add_existing_input(inputs, global.workdir / "aligns" / (batch_prefix(global) + "D_alignments.csv"));
        add_existing_input(inputs, global.workdir / "aligns" / (batch_prefix(global) + "J_alignments.csv"));
    }
    return inputs;
}

std::vector<fs::path> pipeline_manifest_inputs(const ResolvedOptions &options)
{
    std::vector<fs::path> inputs;
    if (pipeline_has_step(options, "import-seqs")) {
        add_existing_input(inputs, options.sequences.input);
    }
    if (pipeline_has_step(options, "align")) {
        add_existing_input(inputs, options.alignment.genomic_v);
        add_existing_input(inputs, options.alignment.genomic_d);
        add_existing_input(inputs, options.alignment.genomic_j);
        add_existing_input(inputs, options.alignment.anchors_v);
        add_existing_input(inputs, options.alignment.anchors_j);
    }
    if ((pipeline_has_step(options, "infer") || pipeline_has_step(options, "evaluate") ||
         pipeline_has_step(options, "generate")) &&
        options.model.source == ModelSource::Custom) {
        add_existing_input(inputs, options.model.parms);
        add_existing_input(inputs, options.model.marginals);
    }
    return inputs;
}

std::vector<fs::path> manifest_inputs_for(const std::string &command, const GlobalOptions &global,
                                          const ResolvedOptions &options,
                                          const std::vector<std::string> &positional)
{
    if (command == "run") {
        return pipeline_manifest_inputs(options);
    }
    return manifest_inputs(command, global, options, positional);
}

struct RuntimeMetadata
{
    std::string os_name = "unknown";
    std::string os_release = "unknown";
    std::string os_version = "unknown";
    std::string machine = "unknown";
};

RuntimeMetadata runtime_metadata()
{
    RuntimeMetadata metadata;
#if defined(__unix__) || defined(__APPLE__)
    struct utsname info;
    if (uname(&info) == 0) {
        metadata.os_name = info.sysname;
        metadata.os_release = info.release;
        metadata.os_version = info.version;
        metadata.machine = info.machine;
    }
#endif
    return metadata;
}

void write_config_tables(std::ostream &out, const std::string &prefix, const Config &config)
{
    std::vector<std::string> current_table;
    for (const auto &item : config_schema()) {
        const auto parts = split_key(item.key);
        std::vector<std::string> table_parts;
        if (!prefix.empty()) {
            table_parts.push_back(prefix);
        }
        table_parts.insert(table_parts.end(), parts.begin(), parts.end() - 1);
        if (table_parts != current_table) {
            out << "\n[";
            for (std::size_t i = 0; i < table_parts.size(); ++i) {
                if (i) {
                    out << ".";
                }
                out << table_parts[i];
            }
            out << "]\n";
            current_table = table_parts;
        }
        out << parts.back() << " = " << toml_config_value_text(item.type, config.get(item.key)) << "\n";
    }
}

std::string joined_argv(const std::vector<std::string> &argv)
{
    std::string joined;
    for (std::size_t i = 0; i < argv.size(); ++i) {
        if (i) {
            joined += ' ';
        }
        joined += argv[i];
    }
    return joined;
}

fs::path write_manifest(const std::string &command, const std::vector<std::string> &argv, const GlobalOptions &global,
                        const Config &config, const ResolvedOptions &options,
                        const std::vector<std::string> &positional, int status)
{
    fs::create_directories(global.workdir / ".igor" / "runs");
    fs::path path = global.workdir / ".igor" / "runs" / (timestamp() + "-" + command + ".toml");

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write run manifest: " + path.string());
    }
    const RuntimeMetadata runtime = runtime_metadata();

    out << "[run]\n";
    out << "command = " << toml_quote(command) << "\n";
    out << "status = " << static_cast<int64_t>(status) << "\n";
    out << "version = " << toml_quote(IGOR_VERSION) << "\n";
    out << "workdir = " << toml_quote(fs::absolute(global.workdir).string()) << "\n";
    out << "batch = " << toml_quote(global.batch) << "\n";
    out << "threads = " << static_cast<int64_t>(global.threads ? *global.threads : options.generate.threads) << "\n";
    out << "seed = " << toml_quote(std::to_string(options.generate.seed)) << "\n";
    out << "argv = " << toml_quote(joined_argv(argv)) << "\n";
    if (!positional.empty()) {
        out << "arg0 = " << toml_quote(positional[0]) << "\n";
    }

    out << "\n[build]\n";
    out << "git_commit = " << toml_quote(IGOR_GIT_COMMIT) << "\n";
    out << "compiled_at = " << toml_quote(std::string(__DATE__) + " " + __TIME__) << "\n";
    out << "cxx_standard = " << static_cast<int64_t>(__cplusplus) << "\n";
    out << "compiler_id = " << toml_quote(IGOR_COMPILER_ID) << "\n";
    out << "compiler_version = " << toml_quote(IGOR_COMPILER_VERSION) << "\n";
    out << "compiler_path = " << toml_quote(IGOR_COMPILER_PATH) << "\n";
    out << "target_system = " << toml_quote(IGOR_TARGET_SYSTEM_NAME) << "\n";
    out << "target_processor = " << toml_quote(IGOR_TARGET_SYSTEM_PROCESSOR) << "\n";

    out << "\n[runtime]\n";
    out << "os_name = " << toml_quote(runtime.os_name) << "\n";
    out << "os_release = " << toml_quote(runtime.os_release) << "\n";
    out << "os_version = " << toml_quote(runtime.os_version) << "\n";
    out << "machine = " << toml_quote(runtime.machine) << "\n";

    write_config_tables(out, "resolved_config", config);

    out << "\n[inputs]\n";
    int i = 0;
    for (const auto &input : manifest_inputs_for(command, global, options, positional)) {
        out << "path" << i << " = " << toml_quote(input.string()) << "\n";
        out << "sha256_" << i << " = " << toml_quote(igor_cli::sha256_file(input)) << "\n";
        ++i;
    }
    out << "count = " << static_cast<int64_t>(i) << "\n";

    out << "\n[artifacts]\n";
    out << "aligns = " << toml_quote((global.workdir / "aligns").string()) << "\n";
    out << "output = " << toml_quote((global.workdir / (batch_prefix(global) + "output")).string()) << "\n";
    out << "inference = " << toml_quote((global.workdir / (batch_prefix(global) + "inference")).string()) << "\n";
    out << "evaluate = " << toml_quote((global.workdir / (batch_prefix(global) + "evaluate")).string()) << "\n";
    out << "generated = " << toml_quote((global.workdir / (batch_prefix(global) + "generated")).string()) << "\n";

    return path;
}

std::string manifest_string(const toml::table &manifest, const std::string &key)
{
    return toml_value_to_string(toml_at_path(manifest, key), key);
}

int manifest_int(const toml::table &manifest, const std::string &key)
{
    return std::stoi(manifest_string(manifest, key));
}

bool manifest_has(const toml::table &manifest, const std::string &key)
{
    return static_cast<bool>(toml_at_path(manifest, key));
}

Config read_manifest_config(const fs::path &manifest, std::string &command, GlobalOptions &global,
                            std::vector<std::string> &positional)
{
    Config config = Config::defaults();
    if (!fs::exists(manifest)) {
        throw std::runtime_error("cannot open replay manifest: " + manifest.string());
    }

    toml::table parsed;
    try {
        parsed = toml::parse_file(manifest.string());
    } catch (const toml::parse_error &e) {
        throw std::runtime_error("invalid replay manifest " + manifest.string() + ": " + std::string(e.description()));
    }

    command = manifest_string(parsed, "run.command");
    global.workdir = manifest_string(parsed, "run.workdir");
    global.batch = manifest_string(parsed, "run.batch");
    if (manifest_has(parsed, "run.threads")) {
        int threads = manifest_int(parsed, "run.threads");
        if (threads > 0) {
            global.threads = threads;
        }
    }
    if (manifest_has(parsed, "run.arg0")) {
        positional.push_back(manifest_string(parsed, "run.arg0"));
    }

    const toml::table *resolved_config = toml_at_path(parsed, "resolved_config").as_table();
    if (resolved_config == nullptr) {
        throw std::runtime_error("replay manifest is missing resolved_config table");
    }
    std::set<std::string> resolved_keys;
    collect_toml_leaf_keys(*resolved_config, "", resolved_keys);
    const auto raw_values = raw_toml_leaf_values(manifest);
    for (const auto &key : resolved_keys) {
        if (!config.values.contains(key)) {
            throw std::runtime_error("unknown resolved config key in replay manifest: " + key);
        }
        const ConfigKey *schema_key = find_config_key(key);
        const auto raw = raw_values.find("resolved_config." + key);
        if (schema_key != nullptr && schema_key->type == "double" && raw != raw_values.end() && !raw->second.empty() &&
            raw->second.front() != '"') {
            config.set(key, raw->second);
        } else {
            config.set(key, toml_value_to_string(toml_at_path(*resolved_config, key), key));
        }
    }

    int count = manifest_int(parsed, "inputs.count");
    for (int i = 0; i < count; ++i) {
        fs::path path = manifest_string(parsed, "inputs.path" + std::to_string(i));
        if (!fs::exists(path)) {
            throw std::runtime_error("replay input is missing: " + path.string());
        }
        std::string actual = igor_cli::sha256_file(path);
        std::string expected = manifest_string(parsed, "inputs.sha256_" + std::to_string(i));
        if (actual != expected) {
            throw std::runtime_error("replay input hash mismatch: " + path.string());
        }
    }
    const ResolvedOptions options = resolve_config(config);
    if (command == "generate") {
        if (options.generate.seed == -1) {
            throw std::runtime_error("replay manifest recorded non-deterministic generation seed -1");
        }
        if (!global.threads && options.generate.threads == 0) {
            throw std::runtime_error("replay manifest recorded auto generation thread count 0");
        }
    }
    return config;
}

void print_config_list()
{
    std::cout << "KEY\tTYPE\tDEFAULT\tDESCRIPTION\n";
    for (const auto &item : config_schema()) {
        std::cout << item.key << "\t" << item.type << "\t" << item.default_value << "\t" << item.description << "\n";
    }
}

int run_command(const std::string &command, const GlobalOptions &global, const Config &config,
                const std::vector<std::string> &positional, const std::vector<std::string> &argv)
{
    const ResolvedOptions options = resolve_config(config);
    auto legacy = command_args(command, global, options, positional);
    int status = run_legacy(legacy);
    fs::path manifest = write_manifest(command, argv, global, config, options, positional, status);
    std::clog << "Run manifest written to: " << manifest << "\n";
    return status;
}

int run_pipeline(const GlobalOptions &global, const Config &config, const std::vector<std::string> &argv)
{
    int status = EXIT_SUCCESS;
    const ResolvedOptions options = resolve_config(config);
    for (const auto &step : options.pipeline.steps) {
        if (step == "import-seqs") {
            if (options.sequences.input.empty()) {
                throw std::runtime_error("pipeline import-seqs requires config key sequences.input");
            }
            status = run_command("import-seqs", global, config, { options.sequences.input.string() }, argv);
        } else if (step == "align") {
            status = run_command("align", global, config, { "all" }, argv);
        } else if (step == "infer" || step == "evaluate") {
            status = run_command(step, global, config, {}, argv);
        } else if (step == "generate") {
            throw std::runtime_error("pipeline generate requires an explicit sequence count and is not supported by `igor run`");
        } else {
            throw std::runtime_error("unknown pipeline step: " + step);
        }
        if (status != EXIT_SUCCESS) {
            break;
        }
    }
    fs::path manifest = write_manifest("run", argv, global, config, options, {}, status);
    std::clog << "Run manifest written to: " << manifest << "\n";
    return status;
}

struct ParsedCli
{
    GlobalOptions global;
    std::string command;
    std::string config_action;
    std::string config_key;
    std::string config_value;
    std::string input;
    std::string gene;
    std::string count;
    std::optional<fs::path> replay;
    std::vector<std::string> argv;
};

CLI::App *add_command(CLI::App &app, const std::string &name, const std::string &description, ParsedCli &parsed)
{
    CLI::App *command = app.add_subcommand(name, description);
    command->fallthrough();
    command->callback([&parsed, name] { parsed.command = name; });
    return command;
}

int parse_cli(int argc, char *argv[], ParsedCli &parsed)
{
    parsed.argv.assign(argv, argv + argc);

    CLI::App app{ "Config-backed IGoR command line interface.", "igor" };
    app.require_subcommand(1, 1);
    app.fallthrough();
    app.set_version_flag("-v,--version", std::string("IGoR version ") + IGOR_VERSION);
    app.add_option("-w,--workdir", parsed.global.workdir, "Use DIR as the IGoR workdir.")->option_text("DIR");
    app.add_option("-b,--batch", parsed.global.batch, "Prefix batch-scoped outputs with NAME.")->option_text("NAME");
    app.add_option("-j,--threads", parsed.global.threads, "Use N worker threads where supported.")
            ->check(CLI::NonNegativeNumber)
            ->option_text("N");
    app.add_option("--stdout-file", parsed.global.stdout_file, "Append standard output to PATH.")->option_text("PATH");

    add_command(app, "datadir", "Print the compiled IGoR data directory.", parsed);
    add_command(app, "init", "Create .igor/config.toml and standard workdir directories.", parsed);

    CLI::App *config = add_command(app, "config", "Read, write, list, or edit the persistent workdir config.", parsed);
    config->require_subcommand(1, 1);
    auto config_action = [&](const std::string &action) {
        return [&parsed, action] { parsed.config_action = action; };
    };
    CLI::App *config_get = config->add_subcommand("get", "Print a config value.");
    config_get->add_option("KEY", parsed.config_key)->required();
    config_get->callback(config_action("get"));
    CLI::App *config_set = config->add_subcommand("set", "Persist a config value.");
    config_set->add_option("KEY", parsed.config_key)->required();
    config_set->add_option("VALUE", parsed.config_value)->required();
    config_set->callback(config_action("set"));
    config->add_subcommand("list", "List supported config keys.")->callback(config_action("list"));
    config->add_subcommand("schema", "List supported config keys.")->callback(config_action("schema"));
    config->add_subcommand("edit", "Open the persistent config in EDITOR.")->callback(config_action("edit"));

    CLI::App *import = add_command(app, "import-seqs", "Import sequences into workdir alignments.", parsed);
    import->add_option("INPUT", parsed.input, "Input .fasta, .csv, .txt, or extensionless text file.")->required();

    CLI::App *align = add_command(app, "align", "Align imported sequences.", parsed);
    align->add_option("--gene", parsed.gene, "Gene class to align.")
            ->required()
            ->check(CLI::IsMember({ "V", "D", "J", "all" }))
            ->option_text("V|D|J|all");

    add_command(app, "infer", "Infer a model from workdir alignments.", parsed);
    add_command(app, "evaluate", "Run one evaluation pass.", parsed);

    CLI::App *generate = add_command(app, "generate", "Generate sequences.", parsed);
    generate->add_option("N", parsed.count, "Number of sequences to generate.")->required()->check(CLI::PositiveNumber);

    CLI::App *run = add_command(app, "run", "Run configured pipeline.steps or replay a run manifest.", parsed);
    run->add_option("--replay", parsed.replay, "Replay a previous run manifest.")->option_text("MANIFEST");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }
    return -1;
}

} // namespace

int main(int argc, char *argv[])
{
    try {
        ParsedCli cli;
        const int parse_status = parse_cli(argc, argv, cli);
        if (parse_status >= 0) {
            return parse_status;
        }

        if (cli.command == "datadir") {
            std::cout << data_dir() << "\n";
            return EXIT_SUCCESS;
        }

        if (cli.command == "init") {
            ensure_workdir(cli.global.workdir);
            if (!fs::exists(config_path(cli.global.workdir))) {
                write_config(cli.global.workdir, Config::defaults());
            }
            std::cout << "Initialized IGoR workdir: " << fs::absolute(cli.global.workdir) << "\n";
            return EXIT_SUCCESS;
        }

        if (cli.command == "config") {
            if (cli.config_action == "list" || cli.config_action == "schema") {
                print_config_list();
                return EXIT_SUCCESS;
            }
            Config config = read_config(cli.global.workdir);
            if (cli.config_action == "get") {
                std::cout << config.get(cli.config_key) << "\n";
                return EXIT_SUCCESS;
            }
            if (cli.config_action == "set") {
                config.set(cli.config_key, cli.config_value);
                write_config(cli.global.workdir, config);
                return EXIT_SUCCESS;
            }
            if (cli.config_action == "edit") {
                const char *editor = std::getenv("EDITOR");
                if (editor == nullptr || std::string(editor).empty()) {
                    editor = "vi";
                }
                std::string cmd = std::string(editor) + " " + config_path(cli.global.workdir).string();
                return std::system(cmd.c_str());
            }
            throw std::runtime_error("unknown config action: " + cli.config_action);
        }

        if (cli.command == "run" && cli.replay) {
            std::string replay_command;
            std::vector<std::string> positional;
            Config config = read_manifest_config(*cli.replay, replay_command, cli.global, positional);
            auto replay_argv = std::vector<std::string>{ "igor", "run", "--replay", cli.replay->string() };
            if (replay_command == "run") {
                return run_pipeline(cli.global, config, replay_argv);
            }
            return run_command(replay_command, cli.global, config, positional, replay_argv);
        }

        Config config = read_config(cli.global.workdir);
        if (cli.command == "import-seqs") {
            config.set("sequences.input", cli.input);
            return run_command(cli.command, cli.global, config, { cli.input }, cli.argv);
        }
        if (cli.command == "align") {
            return run_command(cli.command, cli.global, config, { cli.gene }, cli.argv);
        }
        if (cli.command == "infer" || cli.command == "evaluate") {
            return run_command(cli.command, cli.global, config, {}, cli.argv);
        }
        if (cli.command == "generate") {
            return run_command(cli.command, cli.global, config, { cli.count }, cli.argv);
        }
        if (cli.command == "run") {
            return run_pipeline(cli.global, config, cli.argv);
        }

        throw std::runtime_error("unknown command: " + cli.command);
    } catch (const std::exception &e) {
        std::cerr << "[IGoR] ERROR: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
