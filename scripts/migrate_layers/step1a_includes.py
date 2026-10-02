#!/usr/bin/env python3
"""Step 1a of doc/LAYER_REFACTORING_PROPOSAL.md, text part: rewrite every include and every
export macro that the file moves of step1a_moves.sh break. Run from the repository root,
after step1a_moves.sh. Idempotent: a line already in its final form is left alone.

Rules, in this order:
  1. <igor/Core/X.h>  -> <igor/<Layer>/Legacy/X.h> for every ex-Core header (Config.h stays).
  2. <igor/Model/X>   -> <igor/Inference/X> or <igor/Generation/X> for the moved engines.
  3. The handful of quoted includes become full angle-bracket paths.
  4. Files that left the Core DLL take their new library's export header and macro:
     <igorCoreExport.h> -> <igor/<Layer>/Export.h>, CORE_EXPORT -> <LAYER>_EXPORT.
     Core/Legacy files keep <igorCoreExport.h> and CORE_EXPORT until step 1b.
  5. CORE_TESTING_EXPORT / CORE_TESTING_ENABLED follow the Aligner into Alignment; the macro
     block leaves Utils.h for Alignment/Legacy/TestingExport.h.
  6. The tk factories moved out of Model take INFERENCE_EXPORT / GENERATION_EXPORT.
"""

import os
import re
import subprocess
import sys

ROOT = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
os.chdir(ROOT)

LEGACY = {
    "Core": "CoreEnums StdTypedefs Typedef IntStr GeneticCode SeqTypeRegistry SegmentSpan "
            "LayeredArray DynamicSequenceMap SeqOffsetsMap Utils",
    "Alignment": "Aligner AlignerInternal ExtractFeatures CDR3SeqData JournaledQuery",
    "Model": "Rec_Event Genechoice Deletion Insertion Dinuclmarkov EventUtils gene_to_seqtype_migr "
             "JsonDetail ModelJson Errorrate Singleerrorrate Hypermutationglobalerrorrate "
             "HypermutationfullNmererrorrate Model_Parms Model_marginals Counter "
             "QuerySequenceContext ModelContext ScenarioContext ExplorationContext "
             "AccumulationContext Scenario JunctionGeometry SpanProfile SafetyMatrix "
             "UnfilledSegmentLengths BoundTightness",
    "Inference": "GenModel Pgencounter Bestscenarioscounter Coverageerrcounter Errorscounter",
    "Generation": "FastGenerator FastSampling",
}
ENGINES = {
    "Inference": "InferenceEngine InferenceHandler CategoricalInferenceHandler "
                 "MarkovInferenceHandler InferenceHandlerFactory",
    "Generation": "SamplingEngine SamplingHandler CategoricalSamplingHandler "
                  "MarkovSamplingHandler SamplingHandlerFactory",
}

legacy_layer = {stem: layer for layer, stems in LEGACY.items() for stem in stems.split()}
engine_layer = {stem: layer for layer, stems in ENGINES.items() for stem in stems.split()}

INCLUDE_RE = re.compile(r'(#\s*include\s*)[<"]igor/(Core|Model)/([A-Za-z0-9_]+)\.(h|tpp)[>"]')
QUOTED = {
    '#include "Utils.h"': "#include <igor/Core/Legacy/Utils.h>",
    '#include "Typedef.h"': "#include <igor/Core/Legacy/Typedef.h>",
    '#include "igor/Model/Navigator.h"': "#include <igor/Model/Navigator.h>",
    '#include "igor/Model/MarkovSamplingHandler.h"': "#include <igor/Generation/MarkovSamplingHandler.h>",
}


def rewrite_include(match):
    prefix, old_layer, stem, ext = match.groups()
    if old_layer == "Core":
        if stem == "Config":
            return f"{prefix}<igor/Core/Config.h>"
        layer = legacy_layer.get(stem)
        if layer is None:
            sys.exit(f"unassigned Core header: {stem}.{ext}")
        return f"{prefix}<igor/{layer}/Legacy/{stem}.{ext}>"
    layer = engine_layer.get(stem, "Model")
    return f"{prefix}<igor/{layer}/{stem}.{ext}>"


def layer_of_source(path):
    """The library a source file now belongs to, or None outside src/igor."""
    m = re.match(r"src/igor/([A-Za-z]+)/", path)
    return m.group(1) if m else None


def process(path):
    with open(path, encoding="utf-8", errors="surrogateescape") as handle:
        text = handle.read()
    new = INCLUDE_RE.sub(rewrite_include, text)
    for old, repl in QUOTED.items():
        new = new.replace(old, repl)

    layer = layer_of_source(path)
    in_legacy = "/Legacy/" in path
    if layer in ("Alignment", "Model", "Inference", "Generation") and in_legacy:
        macro = layer.upper()
        new = new.replace("#include <igorCoreExport.h>", f"#include <igor/{layer}/Export.h>")
        new = re.sub(r"\bCORE_NO_EXPORT\b", f"{macro}_NO_EXPORT", new)
        new = re.sub(r"\bCORE_EXPORT\b", f"{macro}_EXPORT", new)
    if layer in ("Inference", "Generation") and not in_legacy:
        # The tk factories left the Model DLL: MODEL_EXPORT would be dllimport on a definition.
        macro = layer.upper()
        new = new.replace("#include <igor/Model/Export.h>", f"#include <igor/{layer}/Export.h>")
        new = re.sub(r"\bMODEL_NO_EXPORT\b", f"{macro}_NO_EXPORT", new)
        new = re.sub(r"\bMODEL_EXPORT\b", f"{macro}_EXPORT", new)
    if layer == "Alignment" or path.startswith("tst/igor/Alignment/"):
        new = new.replace("CORE_TESTING_EXPORT", "ALIGNMENT_TESTING_EXPORT")
        new = new.replace("CORE_TESTING_ENABLED", "ALIGNMENT_TESTING_ENABLED")

    if new != text:
        with open(path, "w", encoding="utf-8", errors="surrogateescape") as handle:
            handle.write(new)
        return True
    return False


TESTING_EXPORT_HEADER = "src/igor/Alignment/Legacy/TestingExport.h"
TESTING_EXPORT_TEXT = """#pragma once

#include <igor/Alignment/Export.h>

// Marks a symbol that is internal to Alignment (not part of its installed public API) but still
// needs to cross the shared library boundary for whitebox tests to link against it directly
// (see AlignerInternal.h). Resolves to a real export/import only when ALIGNMENT_TESTING_ENABLED
// is defined (see tst/igor/Alignment/CMakeLists.txt), so production builds keep these symbols
// hidden.
#ifdef ALIGNMENT_TESTING_ENABLED
#  define ALIGNMENT_TESTING_EXPORT ALIGNMENT_EXPORT
#else
#  define ALIGNMENT_TESTING_EXPORT ALIGNMENT_NO_EXPORT
#endif
"""
UTILS_BLOCK_RE = re.compile(
    r"// Marks a symbol that is internal to Core.*?#ifdef CORE_TESTING_ENABLED\n.*?#endif\n\n",
    re.S)


def relocate_testing_export():
    """Rule 5b: the CORE_TESTING_EXPORT block leaves Utils.h for Alignment/Legacy/TestingExport.h."""
    utils = "src/igor/Core/Legacy/Utils.h"
    if os.path.exists(utils):
        with open(utils, encoding="utf-8") as handle:
            text = handle.read()
        new = UTILS_BLOCK_RE.sub("", text, count=1)
        if new != text:
            with open(utils, "w", encoding="utf-8") as handle:
                handle.write(new)
    if os.path.isdir(os.path.dirname(TESTING_EXPORT_HEADER)) and not os.path.exists(TESTING_EXPORT_HEADER):
        with open(TESTING_EXPORT_HEADER, "w", encoding="utf-8") as handle:
            handle.write(TESTING_EXPORT_TEXT)
    aligner = "src/igor/Alignment/Legacy/Aligner.h"
    if os.path.exists(aligner):
        with open(aligner, encoding="utf-8") as handle:
            text = handle.read()
        anchor = "#include <igor/Alignment/Export.h>\n"
        wanted = anchor + "#include <igor/Alignment/Legacy/TestingExport.h>\n"
        if anchor in text and "TestingExport.h" not in text:
            with open(aligner, "w", encoding="utf-8") as handle:
                handle.write(text.replace(anchor, wanted, 1))


# Rule 7: free functions that used to be called inside the Core DLL and are now called across a
# DLL boundary. On Linux and macOS everything is visible; on Windows an unexported function is
# LNK2019 in the consumer. Keyed by header, the names whose declaration line gets the macro.
CROSS_DLL_FUNCTIONS = {
    "src/igor/Core/Legacy/Utils.h": ("CORE_EXPORT", [
        "extract_string_fields", "show_progress_bar", "close_progress_bar",
        "draw_random_64bits_seed", "translate"]),
    "src/igor/Model/Legacy/Deletion.h": ("MODEL_EXPORT", ["make_transversions", "del_numb_compare"]),
    "src/igor/Model/Legacy/Errorrate.h": ("MODEL_EXPORT", ["add_to_err_rate"]),
}


def export_cross_dll_functions():
    for path, (macro, names) in CROSS_DLL_FUNCTIONS.items():
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().split("\n")
        changed = False
        for i, line in enumerate(lines):
            if "_EXPORT" in line or line.startswith((" ", "\t", "//", "/*", "*")):
                continue
            for name in names:
                if re.match(r"[A-Za-z_][\w:<>,\* &]*\b" + re.escape(name) + r"\s*\(", line):
                    lines[i] = f"{macro} {line}"
                    changed = True
                    break
        if changed:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(lines))


def main():
    relocate_testing_export()
    export_cross_dll_functions()
    changed = 0
    for top in ("src", "app", "tst"):
        for dirpath, _, files in os.walk(top):
            for name in files:
                if name.endswith((".h", ".cpp", ".tpp", ".h.in", ".md")):
                    if process(os.path.join(dirpath, name)):
                        changed += 1
    print(f"step1a_includes: {changed} files rewritten")


if __name__ == "__main__":
    main()
