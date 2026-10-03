#!/usr/bin/env python3
"""Step 1c of doc/LAYER_REFACTORING_PROPOSAL.md: promote a header out of <Layer>/Legacy/.

A promotion, for a concept that lives on with its design unchanged:
  1. `git mv src/igor/<Layer>/Legacy/X.h src/igor/<Layer>/X.h` (and X.cpp if there is one);
  2. the file's namespace becomes igor::<layer>; the transitional using-directives go;
  3. the legacy names it still uses are qualified (`legacy::`, `core::legacy::`, ...), and the
     promoted names of the layers below (`core::`);
  4. its includes of already promoted headers lose their `Legacy/`;
  5. a stub is left at the old path: it includes the promoted header and re-declares its names
     in igor::<layer>::legacy, so legacy consumers keep both their include line and their
     unqualified names. The stub goes when the last consumer switches.

    promote.py <group> [<group> ...]      groups are the keys of PROMOTIONS below

Run from the repository root. Idempotent: an already promoted header is skipped. The CMake
file lists are updated in place; build and run check_directives.py afterwards.
"""

import os
import re
import subprocess
import sys

ROOT = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
os.chdir(ROOT)

# Names still owned by the legacy namespaces, by owner.
CORE_LEGACY = ["Event_type", "GeneChoice_t", "Deletion_t", "Insertion_t", "Dinuclmarkov_t", "Undefined_t",
               "Seq_side", "Five_prime", "Three_prime", "Undefined_side",
               "Seq_type", "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq", "VJ_ins_seq",
               "Seq_type_String", "Rec_Event_name", "Marginal_array_p", "Seq_Offset", "UMCodonTable",
               "Int_Str", "Int_Str_ptr", "Gene_class", "Gene_class_legacy", "V_gene", "D_gene", "J_gene",
               "Undefined_gene", "Int_nt", "Matrix", "Seq_type_str_p_map", "Mismatch_vectors_map",
               "Pruning_mismatch_floor_map", "Index_map", "Downstream_scenario_proba_bound_map",
               "index_type", "first_unfilled_segment", "first_unplaced_segment_end", "kIntNtCount",
               "int_A", "int_C", "int_G", "int_T", "int_R", "int_Y", "int_K", "int_M", "int_S", "int_W",
               "int_B", "int_D", "int_H", "int_V", "int_N", "int_undefined"]
MODEL_LEGACY = ["Rec_Event", "Events_map", "Next_event_ptr", "Gene_choice", "Deletion", "Insertion",
                "Dinucl_markov", "Event_realization", "Error_rate", "OffsetDelta", "LengthContribution",
                "SeqConstructionRole", "OffsetRole", "Model_Parms", "Model_marginals", "Counter",
                "EventUtils"]


# Core vocabulary promoted under a new name (Core/Types.h, Core/IntStr.h): old -> new. Promoted
# and new code uses the new name; the stubs alias the old one for legacy code.
CORE_RENAMED = {"Event_type": "EventType", "Seq_side": "SeqSide", "Seq_type": "SeqType",
                "Seq_type_String": "SeqTypeString", "Rec_Event_name": "EventName",
                "Marginal_array_p": "MarginalArrayPtr", "Seq_Offset": "SeqOffset",
                "UMCodonTable": "CodonTable", "Int_Str": "IntStr"}
# Core vocabulary promoted under its own name.
CORE_PROMOTED_AS_IS = ["GeneChoice_t", "Deletion_t", "Insertion_t", "Dinuclmarkov_t", "Undefined_t",
                       "Five_prime", "Three_prime", "Undefined_side",
                       "V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq", "VJ_ins_seq",
                       "index_type", "genetic_code"]


def vocabulary_promoted():
    return os.path.exists("src/igor/Core/Types.h")


def modernize_core_names(text, in_core):
    """In promoted or new code, replace the legacy spelling of promoted Core vocabulary."""
    prefix = "" if in_core else "core::"
    for old, new in CORE_RENAMED.items():
        text = re.sub(r"(?:(?:igor::)?core::)?legacy::" + old + r"\b", (prefix if not in_core else "") + new, text) \
            if in_core else re.sub(r"((?:igor::)?)core::legacy::" + old + r"\b", r"\1core::" + new, text)
    for name in CORE_PROMOTED_AS_IS:
        text = re.sub(r"(?:(?:igor::)?core::)?legacy::" + name + r"\b", name, text) \
            if in_core else re.sub(r"((?:igor::)?)core::legacy::" + name + r"\b", r"\1core::" + name, text)
    for old in ("CoreEnums", "StdTypedefs", "Typedef"):
        text = text.replace(f"#include <igor/Core/Legacy/{old}.h>", "#include <igor/Core/Types.h>")
    for stem in ("IntStr", "GeneticCode"):
        text = text.replace(f"#include <igor/Core/Legacy/{stem}.h>", f"#include <igor/Core/{stem}.h>")
    # the three legacy headers collapse into one include
    lines, seen = [], False
    for line in text.split("\n"):
        if line == "#include <igor/Core/Types.h>":
            if seen:
                continue
            seen = True
        lines.append(line)
    return "\n".join(lines)


class Symbol:
    """A top-level name of a promoted header and how the legacy stub re-declares it."""

    def __init__(self, name, kind="using", old=None):
        self.name, self.kind, self.old = name, kind, old or name


# group -> layer, list of (stem, [Symbol, ...], {old class name: new class name})
PROMOTIONS = {
    "core-containers": ("Core", [
        ("LayeredArray", [Symbol("LayeredArray")], {}),
        ("SeqTypeRegistry", [Symbol("SeqTypeRegistry"), Symbol("SeqTypeId"), Symbol("kNoSeqType"),
                             Symbol("kLegacySeqTypeCount"), Symbol("legacy_seq_type_registry")], {}),
        ("DynamicSequenceMap", [Symbol("SeqSegmentEmptiness"), Symbol("DynamicSequenceMap")], {}),
        ("SeqOffsetsMap", [Symbol("SeqOffsetsMap", "alias", "Seq_offsets_map")],
         {"Seq_offsets_map": "SeqOffsetsMap"}),
        ("SegmentSpan", [Symbol("SegmentBoundary"), Symbol("SegmentSpan"), Symbol("cut_position"),
                         Symbol("legacy_span_of"), Symbol("legacy_junction_of")], {}),
    ]),
    "core-vocabulary": ("Core", [
        ("IntStr", [Symbol("IntStr", "alias", "Int_Str")], {"Int_Str": "IntStr"}),
        ("GeneticCode", [Symbol("genetic_code", "namespace")], {}),
    ]),
    "model-geometry": ("Model", [
        ("SpanProfile", [Symbol("SpanProfile"), Symbol("SpanDecomposition"), Symbol("JunctionBound")], {}),
        ("SafetyMatrix", [Symbol("kNoOrderingPosition"), Symbol("SafetyCell"), Symbol("SafetyMatrix")], {}),
        ("UnfilledSegmentLengths", [Symbol("UnfilledSegmentLengths")], {}),
        ("BoundTightness", [Symbol("BoundTightness", "namespace")], {}),
        ("JunctionGeometry", [Symbol("JunctionGeometry", "namespace")], {}),
    ]),
}


def read(path):
    with open(path, encoding="utf-8", errors="surrogateescape") as handle:
        return handle.read()


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", errors="surrogateescape") as handle:
        handle.write(text)


def qualify(text, names, prefix):
    """Prefix each whole-word, not yet qualified occurrence of `names`, skipping string literals,
    comments and preprocessor lines."""
    if not names:
        return text
    pattern = re.compile(r"(?<![\w:])(" + "|".join(sorted(names, key=len, reverse=True)) + r")\b")
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i); j = n if j < 0 else j
            out.append(text[i:j]); i = j; continue
        if text.startswith("/*", i):
            j = text.find("*/", i + 2); j = n if j < 0 else j + 2
            out.append(text[i:j]); i = j; continue
        if c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append(text[i:j + 1]); i = j + 1; continue
        if c == "#" and (i == 0 or text[i - 1] == "\n"):
            j = text.find("\n", i); j = n if j < 0 else j
            out.append(text[i:j]); i = j; continue
        j = i
        while j < n and text[j] not in '"#/':
            j += 1
        if j == i:
            out.append(c); i += 1; continue
        out.append(pattern.sub(prefix + r"\1", text[i:j]))
        i = j
    return "".join(out)


def promoted_names():
    """Promoted symbols by layer, from the table and from what is already on disk."""
    names = {}
    for layer, entries in PROMOTIONS.values():
        for stem, symbols, _ in entries:
            if os.path.exists(f"src/igor/{layer}/{stem}.h"):
                names.setdefault(layer, []).extend(s.name for s in symbols)
    return names


def license_of(text):
    end = text.find("#pragma once")
    return text[:end] if end > 0 else ""


def promoted_includes(layer, stem):
    """The promoted headers that the promoted header includes: the stub includes their stubs, so
    a legacy consumer keeps the legacy names it used to get transitively."""
    found = re.findall(r"#include <igor/([A-Za-z]+)/([A-Za-z_0-9]+)\.h>", read(f"src/igor/{layer}/{stem}.h"))
    return [(l, s) for l, s in found if os.path.exists(f"src/igor/{l}/Legacy/{s}.h") and s != stem]


def stub_text(layer, stem, symbols, license_text):
    lines = [license_text + "#pragma once", "",
             f"#include <igor/{layer}/{stem}.h>"]
    lines += [f"#include <igor/{l}/Legacy/{s}.h>" for l, s in promoted_includes(layer, stem)]
    lines += ["",
             f"// {stem}.h was promoted out of Legacy/ (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).",
             "// This stub keeps the legacy include path and the legacy names for the code that has not",
             "// been promoted yet; it goes when its last consumer switches.",
             f"namespace igor::{layer.lower()}::legacy {{", ""]
    for s in symbols:
        target = f"igor::{layer.lower()}::{s.name}"
        if s.kind == "namespace":
            lines.append(f"namespace {s.old} = {target};")
        elif s.kind == "alias":
            lines.append(f"using {s.old} = {target};")
        else:
            lines.append(f"using {target};")
    lines += ["", f"}} // namespace igor::{layer.lower()}::legacy", ""]
    return "\n".join(lines)


def promote(layer, stem, symbols, renames):
    low = layer.lower()
    old_h = f"src/igor/{layer}/Legacy/{stem}.h"
    new_h = f"src/igor/{layer}/{stem}.h"
    if os.path.exists(new_h):
        return False
    for ext in ("h", "cpp"):
        src = f"src/igor/{layer}/Legacy/{stem}.{ext}"
        if os.path.exists(src):
            subprocess.check_call(["git", "mv", src, f"src/igor/{layer}/{stem}.{ext}"])
    for ext in ("h", "cpp"):
        path = f"src/igor/{layer}/{stem}.{ext}"
        if not os.path.exists(path):
            continue
        text = read(path)
        text = text.replace(f"namespace igor::{low}::legacy {{", f"namespace igor::{low} {{")
        text = text.replace(f"}} // namespace igor::{low}::legacy", f"}} // namespace igor::{low}")
        text = re.sub(r"^using namespace igor::[a-z]+::legacy;\n", "", text, flags=re.M)
        text = text.replace(f"igor::{low}::legacy::", f"igor::{low}::")
        for old, new in renames.items():
            text = re.sub(r"\b" + old + r"\b", new, text)
        own = {s.name for s in symbols}
        promoted = promoted_names()
        core_legacy = CORE_LEGACY
        if vocabulary_promoted():
            gone = set(CORE_RENAMED) | set(CORE_PROMOTED_AS_IS)
            core_legacy = [n for n in CORE_LEGACY if n not in gone]
            for old, new in CORE_RENAMED.items():
                if new not in own:
                    text = qualify(text, [old], "\0")
                    text = text.replace("\0" + old, new if layer == "Core" else "core::" + new)
            if layer != "Core":
                text = qualify(text, CORE_PROMOTED_AS_IS, "core::")
        if layer == "Core":
            text = qualify(text, [n for n in core_legacy if n not in own], "legacy::")
        else:
            text = qualify(text, [n for n in promoted.get("Core", []) if n not in own], "core::")
            text = qualify(text, core_legacy, "core::legacy::")
            text = qualify(text, [n for n in MODEL_LEGACY if n not in own], "legacy::")
        # includes of promoted headers lose their Legacy/
        for other_layer, entries in PROMOTIONS.values():
            for other_stem, _, _ in entries:
                if os.path.exists(f"src/igor/{other_layer}/{other_stem}.h") or other_stem == stem:
                    text = text.replace(f"#include <igor/{other_layer}/Legacy/{other_stem}.h>",
                                        f"#include <igor/{other_layer}/{other_stem}.h>")
        write(path, text)
    write(old_h, stub_text(layer, stem, symbols, license_of(read(new_h))))
    subprocess.check_call(["git", "add", old_h])
    # CMake: the promoted header joins the list; the stub keeps the Legacy/ entry.
    cmake = f"src/igor/{layer}/CMakeLists.txt"
    text = read(cmake)
    entry = f"      ${{LAYER_DIR}}/Legacy/{stem}.h\n"
    if entry in text and f"      ${{LAYER_DIR}}/{stem}.h\n" not in text:
        text = text.replace(entry, f"      ${{LAYER_DIR}}/{stem}.h\n" + entry, 1)
    text = text.replace(f"    Legacy/{stem}.cpp\n", f"    {stem}.cpp\n")
    write(cmake, text)
    return True


def restub():
    """Rewrite the stub of every promoted header (their includes depend on what is promoted)."""
    for layer, entries in PROMOTIONS.values():
        for stem, symbols, _ in entries:
            new_h = f"src/igor/{layer}/{stem}.h"
            if os.path.exists(new_h):
                write(f"src/igor/{layer}/Legacy/{stem}.h", stub_text(layer, stem, symbols, license_of(read(new_h))))


def main():
    groups = sys.argv[1:]
    if not groups or any(g not in PROMOTIONS for g in groups):
        sys.exit(f"usage: promote.py <group>...   groups: {', '.join(PROMOTIONS)}")
    for g in groups:
        layer, entries = PROMOTIONS[g]
        for stem, symbols, renames in entries:
            done = promote(layer, stem, symbols, renames)
            print(f"{'promoted' if done else 'skipped '} {layer}/{stem}")
    restub()


if __name__ == "__main__":
    main()
