#!/usr/bin/env python3
"""Step 1b of doc/LAYER_REFACTORING_PROPOSAL.md: namespaces.

Run from the repository root, after step 1a. Idempotent on the files it fully handles (it
checks for its own markers), so it can be re-run on another branch.

Rules:
  1. Every file under src/igor/<Layer>/Legacy/ is wrapped in `namespace igor::<layer>::legacy`.
     Includes found after the first declaration are hoisted above the namespace; `namespace std`
     blocks (hash specializations) are kept outside it and their legacy type names qualified.
  2. Legacy code of a layer sees the legacy names of the layers below it through
     using-directives placed right after the namespace opening (transitional, removed in 1c).
  3. Legacy files that already opened `namespace igor` (Typedef, ModelJson, JsonDetail,
     gene_to_seqtype_migr, FastGenerator, FastSampling) get that namespace renamed instead.
  4. The three `EventUtils` namespaces would collide through the using-directives: the Core one
     becomes `genetic_code`, the Alignment one `journaled_query`, Model keeps `EventUtils`.
  5. The tk engines moved by 1a leave `igor::model` for `igor::inference` / `igor::generation`.
  6. Streaming: `igor` -> `igor::streaming`, `igor::airr` -> `igor::streaming::airr`.
  7. Qualified names that moved: igor::index_type, igor::model_parms_to_json, igor::migration,
     igor::json_detail, igor::fast, plus the engine factories.
  8. New code (tk Model) qualifies legacy types as `legacy::X`; the engines as `model::legacy::X`.
  9. Consumers (tests, apps) get `using namespace igor::<layer>::legacy;` after their includes.
"""

import os
import re
import subprocess
import sys

ROOT = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
os.chdir(ROOT)

LAYERS = ["Core", "Alignment", "Model", "Generation", "Inference"]
# Legacy namespaces visible from each layer's legacy code (the DAG, bottom-up).
VISIBLE = {
    "Core": [],
    "Alignment": ["Core"],
    "Model": ["Core", "Alignment"],
    "Generation": ["Core", "Alignment", "Model"],
    "Inference": ["Core", "Alignment", "Model", "Generation"],
}
CORE_TYPES_IN_STD = ["Seq_type", "Gene_class_legacy", "Gene_class", "Event_type", "Seq_type_String",
                     "Seq_side", "Int_Str"]
# Legacy names that the new code (tk Model, engines, Streaming) uses unqualified, by owner.
CORE_LEGACY_NAMES = ["index_type", "Event_type", "Gene_class_legacy", "Gene_class", "Seq_side", "Seq_type",
                     "Int_Str", "Rec_Event_name", "SeqTypeRegistry", "Seq_type_String", "SeqTypeId",
                     "GeneChoice_t", "Deletion_t", "Insertion_t", "Dinuclmarkov_t", "Undefined_t",
                     "V_gene", "D_gene", "J_gene", "Undefined_gene"]
ALIGNMENT_LEGACY_NAMES = ["Alignment_data", "Aligner"]
MODEL_LEGACY_NAMES = ["Rec_Event", "Model_Parms", "Model_marginals", "Error_rate", "Gene_choice", "Deletion",
                      "Insertion", "Dinucl_markov", "Event_realization", "Single_error_rate", "Events_map",
                      "Counter"]
MODEL_TK_TYPES = ["RecombinationModel", "Navigator", "Topology", "SampledScenario", "SampledEvent"]


def ns(layer):
    return f"igor::{layer.lower()}::legacy"


def read(path):
    with open(path, encoding="utf-8", errors="surrogateescape") as handle:
        return handle.read()


def write(path, text):
    with open(path, "w", encoding="utf-8", errors="surrogateescape") as handle:
        handle.write(text)


def directives(layer):
    return "".join(f"using namespace {ns(l)};\n" for l in VISIBLE[layer])



# ----------------------------------------------------------------------------- include closure

INCLUDE_LINE = re.compile(r'^\s*#\s*include\s*[<"]igor/([^>"]+)[>"]', re.M)
_closure_cache = {}


def include_closure(path):
    """Every repository header reachable from `path` through <igor/...> includes."""
    if path in _closure_cache:
        return _closure_cache[path]
    seen = set()
    stack = [path]
    while stack:
        current = stack.pop()
        if current in seen or not os.path.exists(current):
            continue
        seen.add(current)
        for rel in INCLUDE_LINE.findall(read(current)):
            stack.append(f"src/igor/{rel}")
    _closure_cache[path] = seen
    return seen


def declared_layers(path, layers):
    """The layers whose legacy namespace a header of the closure declares. A using-directive
    for any other layer would be the only thing in the file that fails to compile, so it is
    not emitted."""
    closure = include_closure(path)
    return [l for l in layers if any(f"src/igor/{l}/Legacy/" in h for h in closure)]


def directives_for(path, layers):
    return "".join(f"using namespace {ns(l)};\n" for l in declared_layers(path, layers))

# ----------------------------------------------------------------------------- rule 1 and 2

def first_code_line(lines):
    """Index of the first line that is neither blank, comment nor preprocessor, outside any
    #if block. Code inside a conditional block that precedes it (Utils.h's Windows shims, which
    include winsock2.h) stays where it is: a system header must not be included inside the
    namespace."""
    in_block = False
    depth = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if in_block:
            if "*/" in s:
                in_block = False
            continue
        if s.startswith("#"):
            directive = s[1:].strip()
            if directive.startswith(("if", "ifdef", "ifndef")):
                depth += 1
            elif directive.startswith("endif"):
                depth = max(0, depth - 1)
            continue
        if not s or s.startswith("//"):
            continue
        if s.startswith("/*"):
            if "*/" not in s:
                in_block = True
            continue
        if s.startswith("*"):
            continue
        if depth == 0:
            return i
    return len(lines)


def preprocessor_depths(lines):
    depths = []
    depth = 0
    for line in lines:
        s = line.strip()
        depths.append(depth)
        if s.startswith("#"):
            directive = s[1:].strip()
            if directive.startswith(("if", "ifdef", "ifndef")):
                depth += 1
            elif directive.startswith("endif"):
                depth = max(0, depth - 1)
    return depths


def wrap_generic(path, layer):
    text = read(path)
    # `::f(` meant "the global f, not the member"; f is now in this layer's legacy namespace.
    # The portable_* shims of Utils.h stay global (they sit in a conditional block, see
    # first_code_line), so their calls keep the leading `::`.
    text = re.sub(r"(?<![\w:>])::(?!portable_)([a-z_][a-z0-9_]*\()", r"legacy::\1", text)
    marker = f"namespace {ns(layer)} {{"
    if marker in text:
        return False
    lines = text.split("\n")
    start = first_code_line(lines)
    head, body = lines[:start], lines[start:]
    # hoist includes found in the body, outside conditional blocks only
    depths = preprocessor_depths(body)
    hoisted = [l for l, d in zip(body, depths) if l.lstrip().startswith("#include") and d == 0]
    body = [l for l, d in zip(body, depths) if not (l.lstrip().startswith("#include") and d == 0)]
    if hoisted:
        head = head + hoisted
    usings = directives_for(path, VISIBLE[layer]).rstrip("\n")
    opening = [marker] + (usings.split("\n") if usings else [])
    closing = f"}} // namespace {ns(layer)}"
    # split around top-level `namespace std {` blocks
    out = []
    i = 0
    while i < len(body):
        line = body[i]
        if re.match(r"\s*namespace std\s*\{", line):
            depth = 0
            j = i
            block = []
            while j < len(body):
                block.append(body[j])
                depth += body[j].count("{") - body[j].count("}")
                j += 1
                if depth == 0:
                    break
            block_text = "\n".join(block)
            if layer == "Core":
                for t in CORE_TYPES_IN_STD:
                    block_text = re.sub(r"(?<![\w:])" + t + r"\b", f"{ns('Core')}::{t}", block_text)
            out.append(closing)
            out.append(block_text)
            out.append(marker)
            i = j
            continue
        out.append(line)
        i += 1
    # trailing blank lines before the closing brace
    while out and out[-1].strip() == "":
        out.pop()
    if path.endswith("Core/Legacy/Utils.h"):
        opening = opening + ["", "// The legacy to_string(Gene_class) overloads below would otherwise hide std::to_string",
                             "// from every legacy namespace that nominates this one.", "using std::to_string;"]
    new = "\n".join(head + [""] + opening + [""] + out + ["", closing, ""])
    write(path, new)
    return True


# ----------------------------------------------------------------------------- rule 3

RENAMED_IGOR_FILES = {
    "src/igor/Core/Legacy/Typedef.h": ("Core", None),
    "src/igor/Model/Legacy/ModelJson.h": ("Model", None),
    "src/igor/Model/Legacy/ModelJson.cpp": ("Model", None),
    "src/igor/Model/Legacy/JsonDetail.h": ("Model", "json_detail"),
    "src/igor/Model/Legacy/gene_to_seqtype_migr.h": ("Model", None),
    "src/igor/Model/Legacy/gene_to_seqtype_migr.cpp": ("Model", None),
    "src/igor/Generation/Legacy/FastGenerator.h": ("Generation", None),
    "src/igor/Generation/Legacy/FastGenerator.cpp": ("Generation", None),
    "src/igor/Generation/Legacy/FastSampling.h": ("Generation", None),
    "src/igor/Generation/Legacy/FastSampling.cpp": ("Generation", None),
}


def rename_igor_namespace(path, layer, sub):
    text = read(path)
    target = ns(layer)
    if f"namespace {target}" in text:
        return False
    if sub:
        text = text.replace(f"namespace igor::{sub} {{", f"namespace {target}::{sub} {{\n{directives(layer)}".rstrip("\n") + "\n" if False else f"namespace {target}::{sub} {{")
        text = text.replace(f"// namespace igor::{sub}", f"// namespace {target}::{sub}")
    else:
        if path.endswith("ModelJson.h"):
            text = text.replace("class Model_Parms;\n\nnamespace igor {\n", f"namespace {target} {{\n\nclass Model_Parms;\n", 1)
        text = text.replace("namespace igor {", f"namespace {target} {{", 1)
        text = re.sub(r"\}\s*//\s*namespace igor\s*$", f"}} // namespace {target}", text, count=1, flags=re.M)
        if path.endswith("Typedef.h"):
            text = text.replace("\n}\n", f"\n}} // namespace {target}\n", 1)
    # using-directives after the opening
    usings = directives_for(path, VISIBLE[layer]).rstrip("\n")
    if usings:
        opening = f"namespace {target}::{sub} {{" if sub else f"namespace {target} {{"
        text = text.replace(opening, opening + "\n" + usings, 1)
    write(path, text)
    return True


# ----------------------------------------------------------------------------- rule 4

EVENTUTILS_RENAMES = [
    ("src/igor/Core/Legacy/GeneticCode.h", "genetic_code"),
    ("src/igor/Core/Legacy/GeneticCode.cpp", "genetic_code"),
    ("src/igor/Alignment/Legacy/JournaledQuery.h", "journaled_query"),
    ("src/igor/Alignment/Legacy/JournaledQuery.cpp", "journaled_query"),
]
GENETIC_CODE_FUNCTIONS = ["codon_index", "translate_codon", "codon_mask_for_aa", "translate_int_seq",
                          "aa_to_iupac_codon", "mask_to_iupac_codon", "motif_char_to_mask",
                          "parse_aa_motif", "iupac_from_bits", "CodonMask"]
JOURNALED_QUERY_FUNCTIONS = ["motif_to_journaled_query", "aa_to_journaled_query"]


def rename_eventutils():
    for path, new in EVENTUTILS_RENAMES:
        text = read(path)
        text = text.replace("namespace EventUtils {", f"namespace {new} {{")
        text = text.replace("// namespace EventUtils", f"// namespace {new}")
        if path.endswith("JournaledQuery.cpp"):
            text = text.replace("using namespace EventUtils;", "using namespace genetic_code;")
        write(path, text)


def requalify_eventutils(text):
    for f in GENETIC_CODE_FUNCTIONS:
        text = re.sub(r"\bEventUtils::" + f + r"\b", "genetic_code::" + f, text)
    for f in JOURNALED_QUERY_FUNCTIONS:
        text = re.sub(r"\bEventUtils::" + f + r"\b", "journaled_query::" + f, text)
    return text


# ----------------------------------------------------------------------------- rule 5 and 6

def rename_engine_namespaces():
    for layer in ("Inference", "Generation"):
        for name in os.listdir(f"src/igor/{layer}"):
            path = f"src/igor/{layer}/{name}"
            if not name.endswith((".h", ".tpp", ".cpp")):
                continue
            text = read(path)
            new = text.replace("igor::model::", f"igor::{layer.lower()}::")
            new = re.sub(r"\bigor::model\b", f"igor::{layer.lower()}", new)
            if new != text:
                write(path, new)


def rename_streaming_namespace():
    for d in ("src/igor/Streaming", "tst/igor/Streaming"):
        for name in os.listdir(d):
            path = f"{d}/{name}"
            if not name.endswith((".h", ".cpp")):
                continue
            text = read(path)
            new = text.replace("namespace igor {", "namespace igor::streaming {")
            new = new.replace("// namespace igor\n", "// namespace igor::streaming\n")
            new = new.replace("namespace igor::airr", "namespace igor::streaming::airr")
            new = new.replace("// namespace igor::airr", "// namespace igor::streaming::airr")
            new = new.replace("namespace igor::test", "namespace igor::streaming::test")
            new = new.replace("// namespace igor::test", "// namespace igor::streaming::test")
            if new != text:
                write(path, new)


STREAMING_TOP_LEVEL = ["CompressionType", "Delimiter", "FileInfo", "ParquetReader", "ParquetWriter",
                       "SequenceData"]


def requalify_moved_names(text, path):
    text = text.replace("igor::index_type", f"{ns('Core')}::index_type")
    text = text.replace("igor::model_parms_to_json", f"{ns('Model')}::model_parms_to_json")
    text = text.replace("igor::kModelJsonSchemaVersion", f"{ns('Model')}::kModelJsonSchemaVersion")
    text = text.replace("igor::migration::", f"{ns('Model')}::migration::")
    text = re.sub(r"\bigor::json_detail\b", f"{ns('Model')}::json_detail", text)
    text = re.sub(r"\bigor::fast\b", f"{ns('Generation')}::fast", text)
    text = text.replace("igor::model::inference_handler_factory", "igor::inference::inference_handler_factory")
    text = text.replace("igor::model::sampling_handler_factory", "igor::generation::sampling_handler_factory")
    for t in ("InferenceEngine", "InferenceHandler", "CategoricalInferenceHandler", "MarkovInferenceHandler"):
        text = re.sub(r"\bigor::model::" + t + r"\b", "igor::inference::" + t, text)
    for t in ("SamplingEngine", "SamplingHandler", "CategoricalSamplingHandler", "MarkovSamplingHandler"):
        text = re.sub(r"\bigor::model::" + t + r"\b", "igor::generation::" + t, text)
    if "/Streaming/" in path or "igor-convert" in path:
        text = re.sub(r"\bigor::airr\b", "igor::streaming::airr", text)
        text = re.sub(r"\bigor::test\b", "igor::streaming::test", text)
        for t in STREAMING_TOP_LEVEL:
            text = re.sub(r"\bigor::" + t + r"\b", "igor::streaming::" + t, text)
        text = text.replace("using namespace igor;", "using namespace igor::streaming;")
    return text


# ----------------------------------------------------------------------------- rule 8

def qualify_outside_strings_and_comments(text, names, prefix):
    """Prefix each whole-word occurrence of `names` that is not already qualified, skipping
    string literals, comments and preprocessor lines."""
    pattern = re.compile(r"(?<![\w:])(" + "|".join(names) + r")\b")
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(text[i:j]); i = j; continue
        if text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(text[i:j]); i = j; continue
        if c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                if text[j] == "\\":
                    j += 1
                j += 1
            out.append(text[i:j + 1]); i = j + 1; continue
        if c == "#" and (i == 0 or text[i - 1] == "\n"):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(text[i:j]); i = j; continue
        # a run of ordinary code up to the next special character
        j = i
        while j < n and text[j] not in '"#/':
            j += 1
        if j == i:
            out.append(c); i += 1; continue
        out.append(pattern.sub(prefix + r"\1", text[i:j]))
        i = j
    return "".join(out)


def qualify_new_code():
    for layer, model_prefix in (("Model", "legacy::"), ("Inference", "model::legacy::"), ("Generation", "model::legacy::")):
        for name in os.listdir(f"src/igor/{layer}"):
            path = f"src/igor/{layer}/{name}"
            if not name.endswith((".h", ".tpp", ".cpp")):
                continue
            text = read(path)
            new = qualify_outside_strings_and_comments(text, CORE_LEGACY_NAMES, "core::legacy::")
            new = qualify_outside_strings_and_comments(new, ALIGNMENT_LEGACY_NAMES, "alignment::legacy::")
            new = qualify_outside_strings_and_comments(new, MODEL_LEGACY_NAMES, model_prefix)
            if layer != "Model":
                # The engines name Model's types unqualified. Forward-declare them in igor::model
                # and bring them in with using-declarations; a stray forward declaration left
                # inside the engine namespace would declare an unrelated type, so drop it.
                new = re.sub(r"^\s*template *<[^>]*> *class RecombinationModel;.*\n", "", new, flags=re.M)
                new = re.sub(r"^(class|struct) (Topology|Navigator|SampledScenario|SampledEvent);.*\n", "", new, flags=re.M)
                opening = f"namespace igor::{layer.lower()} {{\n"
                if opening in new and "using igor::model::RecombinationModel;" not in new:
                    forward = "#include <igor/Model/Forward.h>\n\n"
                    decl = "".join(f"using igor::model::{t};\n" for t in MODEL_TK_TYPES)
                    new = new.replace(opening, forward + opening + decl, 1)
            if new != text:
                write(path, new)


# ----------------------------------------------------------------------------- rule 9

CONSUMER_DIRECTIVES = {
    "tst/igor/Core": ["Core"],
    "tst/igor/Alignment": ["Core", "Alignment"],
    "tst/igor/Model": ["Core", "Alignment", "Model"],
    "tst/igor/Generation": ["Core", "Alignment", "Model", "Generation"],
    "tst/igor/Inference": ["Core", "Alignment", "Model", "Generation", "Inference"],
    "tst/igor/entropy_test_helpers.h": ["Core", "Alignment", "Model", "Generation"],
    "app/igor": ["Core", "Alignment", "Model", "Generation", "Inference"],
    "app/igor-demo": ["Core", "Alignment", "Model", "Generation", "Inference"],
    "tst/igor/Streaming": ["Core", "Alignment"],
    "app/igor-convert": ["Core", "Alignment"],
}


def add_consumer_directives(path, layers):
    text = read(path)
    if "using namespace igor::core::legacy;" in text or "#include <igor/" not in text:
        return False
    lines = text.split("\n")
    stop = first_code_line(lines)
    last_include = max((i for i, l in enumerate(lines[:stop]) if l.lstrip().startswith("#include")), default=-1)
    if last_include < 0:
        return False
    closure = include_closure(path)
    nominated = [ns(l) for l in declared_layers(path, layers)]
    # tests of the engines also open the engine namespace itself
    if path.startswith("tst/igor/Inference/") and any("src/igor/Inference/" in h for h in closure):
        nominated.append("igor::inference")
    if path.startswith("tst/igor/Generation/") and any("src/igor/Generation/" in h for h in closure):
        nominated.append("igor::generation")
    if not nominated:
        return False
    block = [f"using namespace {n};" for n in nominated]
    lines[last_include + 1:last_include + 1] = [""] + block
    text = "\n".join(lines)
    if "using namespace EventUtils;" in text:
        owners = [("GeneticCode.h", "igor::core::legacy::genetic_code"),
                  ("JournaledQuery.h", "igor::alignment::legacy::journaled_query"),
                  ("EventUtils.h", "igor::model::legacy::EventUtils")]
        repl = "".join(f"using namespace {n};\n" for h, n in owners if any(c.endswith(h) for c in closure))
        text = text.replace("using namespace EventUtils;", repl.rstrip("\n"))
    write(path, text)
    return True


STREAMING_CORE_NAMES = ["Gene_class", "Gene_class_legacy", "Int_Str", "Seq_type", "Seq_side",
                        "V_gene", "D_gene", "J_gene", "Undefined_gene"]
STREAMING_ALIGNMENT_NAMES = ["Alignment_data", "Aligner"]


def qualify_streaming():
    """Streaming is new code in igor::streaming: it names the legacy types it consumes in full."""
    paths = [f"src/igor/Streaming/{n}" for n in os.listdir("src/igor/Streaming") if n.endswith((".h", ".cpp"))]
    paths.append("tst/igor/Streaming/StreamingTestUtils.h")
    for path in paths:
        text = read(path)
        # Fully qualified: inside igor::streaming::airr::alignment, `alignment::` would name the
        # enclosing namespace, not igor::alignment.
        new = qualify_outside_strings_and_comments(text, STREAMING_CORE_NAMES, "igor::core::legacy::")
        new = qualify_outside_strings_and_comments(new, STREAMING_ALIGNMENT_NAMES, "igor::alignment::legacy::")
        if new != text:
            write(path, new)


FRIEND_FORWARD = "namespace igor::inference::legacy { class Coverage_err_counter; }\n"


def fix_cross_layer_friends():
    """Deletion and Gene_choice befriend Coverage_err_counter, which 1a moved to Inference. Inside
    igor::model::legacy the bare name would declare a new, unrelated class, so the friend is
    named in full behind a forward declaration. This is an upward reference from Model to
    Inference, to be removed with the observer rework of step 1c."""
    for path in ("src/igor/Model/Legacy/Deletion.h", "src/igor/Model/Legacy/Genechoice.h"):
        text = read(path)
        if "igor::inference::legacy::Coverage_err_counter" in text:
            continue
        text = text.replace("friend class Coverage_err_counter;",
                            "friend class igor::inference::legacy::Coverage_err_counter;")
        marker = "namespace igor::model::legacy {"
        text = text.replace(marker, "// Friend of the events below; see the friend declarations.\n" + FRIEND_FORWARD + marker, 1)
        write(path, text)


def fix_deletion_test():
    path = "tst/igor/Model/Legacy/test_EventUtils.cpp"
    text = read(path)
    if "namespace igor::model::legacy {\nclass DeletionTest" in text:
        return
    start = text.index("class DeletionTest {")
    end = text.index("\n};\n", start) + len("\n};\n")
    text = text[:start] + "namespace igor::model::legacy {\n" + text[start:end] + "} // namespace igor::model::legacy\n" + text[end:]
    write(path, text)


def main():
    # legacy sources
    rename_eventutils()
    for layer in LAYERS:
        d = f"src/igor/{layer}/Legacy"
        for name in sorted(os.listdir(d)):
            path = f"{d}/{name}"
            if not name.endswith((".h", ".cpp")) or name == "TestingExport.h":
                continue
            if path in RENAMED_IGOR_FILES:
                rename_igor_namespace(path, *RENAMED_IGOR_FILES[path])
            else:
                wrap_generic(path, layer)
    fix_cross_layer_friends()
    rename_engine_namespaces()
    rename_streaming_namespace()
    qualify_new_code()
    qualify_streaming()
    # consumers
    fix_deletion_test()
    for top, layers in CONSUMER_DIRECTIVES.items():
        paths = [top] if os.path.isfile(top) else [
            os.path.join(dp, f) for dp, _, fs in os.walk(top) for f in fs if f.endswith((".cpp", ".h"))]
        for path in paths:
            add_consumer_directives(path, layers)
    # qualified names that moved, everywhere
    for top in ("src", "app", "tst"):
        for dp, _, fs in os.walk(top):
            for f in fs:
                if f.endswith((".h", ".cpp", ".tpp", ".md")):
                    path = os.path.join(dp, f)
                    text = read(path)
                    new = requalify_moved_names(requalify_eventutils(text), path)
                    if new != text:
                        write(path, new)
    print("step1b_namespaces: done")


if __name__ == "__main__":
    main()
