#!/usr/bin/env python3
"""Build the V-D1-D2-J model of milestone 1 from human TRB (capability plan T5, step T5a).

TRB's D block is duplicated, and the result is written as one JSON model document: its
structure, its template library, its parameters and its error rate.

  - The D gene and its two deletions are copied twice, as D1 and D2. Both gene choices draw
    from TRB's D library and are conditioned on V and J, as TRB's is; each 3' deletion is
    conditioned on its own gene and 5' deletion, as in TRB.
  - TRB's two junctions become three: VD1 copies VD, D2J copies DJ, and D1D2 copies VD (see
    COPIES below).
  - V and J, their deletions and the error rate are TRB's as they are.

Every event's realizations, labels and parameters are those of the TRB event it copies. The
parameters are written by label, as model_parms_to_json() and model_marginals_to_json() lay them
out, so IGoR reads them whatever order it gives the events. The values are TRB's as its
marginals file writes them; IGoR normalizes them on reading, as it does the text file's.

The model is a computation test, biologically implausible by construction: a fixed V-D1-D2-J
layout cannot explain the V-D-J reads a real repertoire also holds.

    make_tandem_d_model.py [--source DIR] [--output FILE]

Run from anywhere: the defaults are relative to the repository. The output is deterministic, so
running the script again on the same source leaves the committed document unchanged.
"""

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "models" / "human" / "tcr_beta" / "models"
DEFAULT_OUTPUT = ROOT / "tst" / "test_data" / "tandem_d" / "trb_tandem_d_model.json"

# The version of the document the reader knows, kModelJsonSchemaVersion in ModelJson.h.
SCHEMA_VERSION = "0.3"
SCHEMA_URI = ("https://raw.githubusercontent.com/IGoR-immuno-stack/IGoR/master/src/igor/Model/schemas/"
              + SCHEMA_VERSION + "/model.schema.json")

SEGMENTS = ["V_gene_seq", "VD1_ins_seq", "D1_gene_seq", "D1D2_ins_seq", "D2_gene_seq", "D2J_ins_seq",
            "J_gene_seq"]

# The tandem model's events, by group: each maps the TRB events it copies to the copies'
# nicknames, all on the group's segment. A parent of a copy is renamed by the copy's group, so
# D1's 3' deletion is conditioned on D1's gene and 5' deletion; a parent outside the group keeps
# its name, so both Ds are conditioned on V and J.
#
# D1D2 copies TRB's VD junction, not DJ: the chain then runs 5' to 3' from D1, as VD's runs from
# V, and as the tandem fixture of the unit cases has it (tst/test_data/model_formats/
# tandem_d_model.json). Either junction would do for a computation test.
COPIES = [
    ("V_gene_seq", {"v_choice": "v_choice", "v_3_del": "v_3_del"}),
    ("VD1_ins_seq", {"vd_ins": "vd1_ins", "vd_dinucl": "vd1_dinucl"}),
    ("D1_gene_seq", {"d_gene": "d1_choice", "d_5_del": "d1_5_del", "d_3_del": "d1_3_del"}),
    ("D1D2_ins_seq", {"vd_ins": "d1d2_ins", "vd_dinucl": "d1d2_dinucl"}),
    ("D2_gene_seq", {"d_gene": "d2_choice", "d_5_del": "d2_5_del", "d_3_del": "d2_3_del"}),
    ("D2J_ins_seq", {"dj_ins": "d2j_ins", "dj_dinucl": "d2j_dinucl"}),
    ("J_gene_seq", {"j_choice": "j_choice", "j_5_del": "j_5_del"}),
]

HEADER = """\
/*
 * The V-D1-D2-J model of milestone 1 (capability plan T5), built from human TRB
 * (models/human/tcr_beta/models) by scripts/tandem_d/make_tandem_d_model.py. Do not edit it:
 * change the script and run it again.
 *
 * TRB's D gene and its deletions are copied twice, as D1 and D2, both conditioned on V and J.
 * VD1 and D1D2 copy TRB's VD junction, D2J its DJ junction. The parameters are TRB's, event by
 * event. A computation test, biologically implausible by construction.
 */
"""


class SourceError(Exception):
    """The TRB files are not what the duplication expects."""


def side_of(event):
    """The end a v1 event acts on, as a document writes it, or None.

    A deletion's line names it. A Markov chain's line says Undefined_side, and its gene class
    gives the direction it runs in, as Model_Parms::read_model_parms() reads it.
    """
    if event["type"] == "Deletion":
        return {"Five_prime": "5'", "Three_prime": "3'"}[event["side"]]
    if event["type"] == "DinucMarkov":
        return {"VD_genes": "3'", "VJ_genes": "3'", "DJ_gene": "5'", "DJ_genes": "5'"}[event["class"]]
    return None


def read_parms(path):
    """The events of a v1 model file by nickname, its edges as (parent, child), and its error rate."""
    lines = path.read_text().splitlines()
    if not lines or lines[0] != "@Event_list":
        raise SourceError(f"{path}: not a v1 model file")

    events = {}
    at = 1
    while at < len(lines) and lines[at].startswith("#"):
        kind, gene_class, side, priority, nickname = lines[at][1:].split(";")
        event = {"type": kind, "class": gene_class, "side": side, "priority": int(priority),
                 "realizations": []}
        at += 1
        while at < len(lines) and lines[at].startswith("%"):
            fields = lines[at][1:].split(";")
            if kind == "GeneChoice":
                label, sequence, index = fields
                event["realizations"].append({"index": int(index), "label": label, "sequence": sequence})
            elif kind in ("Deletion", "Insertion"):
                value, index = fields
                event["realizations"].append({"index": int(index), "label": int(value)})
            elif kind == "DinucMarkov":
                label, index = fields
                event["realizations"].append({"index": int(index), "label": label})
            else:
                raise SourceError(f"{path}: event {nickname} has the unknown type {kind}")
            at += 1
        event["realizations"].sort(key=lambda realization: realization["index"])
        if [r["index"] for r in event["realizations"]] != list(range(len(event["realizations"]))):
            raise SourceError(f"{path}: the indices of {nickname} do not run from 0 without gaps")
        if kind == "DinucMarkov" and [r["label"] for r in event["realizations"]] != ["A", "C", "G", "T"]:
            raise SourceError(f"{path}: the chain {nickname} is not over A, C, G, T")
        events[nickname] = event

    # An edge names its events by their v1 names, which the header fields and the realization
    # count make up.
    v1_names = {}
    for nickname, event in events.items():
        name = (f"{event['type']}_{event['class']}_{event['side']}_prio{event['priority']}"
                f"_size{len(event['realizations'])}")
        v1_names.setdefault(name, []).append(nickname)

    def named(token):
        found = v1_names.get(token, [])
        if len(found) != 1:
            raise SourceError(f"{path}: the edge token {token} names {len(found)} events")
        return found[0]

    if at >= len(lines) or lines[at] != "@Edges":
        raise SourceError(f"{path}: no @Edges after the events")
    at += 1
    edges = []
    while at < len(lines) and lines[at].startswith("%"):
        parent, child = lines[at][1:].split(";")
        edges.append((named(parent), named(child)))
        at += 1

    if lines[at:at + 2] != ["@ErrorRate", "#SingleErrorRate"]:
        raise SourceError(f"{path}: the error rate is not a SingleErrorRate")
    return events, edges, float(lines[at + 2])


def read_marginals(path):
    """The blocks of a v1 marginals file by nickname: the parents a block names, and its rows.

    A row is the event's values, in index order, at one combination of its parents' indices, as
    a "#[parent,index],..." line names it.
    """
    blocks = {}
    nickname = None
    coordinates = None
    for line in path.read_text().splitlines():
        if line.startswith("@"):
            nickname = line[1:]
            blocks[nickname] = {"parents": None, "rows": {}}
        elif line.startswith("$"):
            continue
        elif line.startswith("#"):
            pairs = [pair.strip("[]").split(",") for pair in line[1:].split("],[")] if line != "#" else []
            parents = [parent for parent, _ in pairs]
            block = blocks[nickname]
            if block["parents"] is None:
                block["parents"] = parents
            elif block["parents"] != parents:
                raise SourceError(f"{path}: the rows of {nickname} name different parents")
            coordinates = tuple(int(index) for _, index in pairs)
        elif line.startswith("%"):
            blocks[nickname]["rows"][coordinates] = [float(value) for value in line[1:].split(",")]
    return blocks


def nested(rows, sizes, own, prefix=()):
    """A block's rows nested one level per parent axis, each row one level per own axis."""
    if len(prefix) == len(sizes):
        row = rows[prefix]
        if len(own) == 1:
            return row
        return [row[i * own[1]:(i + 1) * own[1]] for i in range(own[0])]
    return [nested(rows, sizes, own, prefix + (index,)) for index in range(sizes[len(prefix)])]


def build_document(events, edges, error_rate, marginals):
    copied = {source for _, group in COPIES for source in group}
    if copied != set(events):
        raise SourceError("the TRB events are not the ones the duplication copies: "
                          + ", ".join(sorted(copied ^ set(events))))
    parents_of = {nickname: [] for nickname in events}
    for parent, child in edges:
        parents_of[child].append(parent)

    kept = {source for _, group in COPIES for source, copy in group.items() if source == copy}
    structure_events = {}
    conditioning = {}
    library = {}
    nodes = {}
    for segment, group in COPIES:
        for source, copy in group.items():
            event = events[source]

            def renamed(parent):
                if parent in group:
                    return group[parent]
                if parent in kept:
                    return parent
                raise SourceError(f"{copy} would be conditioned on {parent}, which the tandem model does not keep")

            target = {"segment": segment}
            side = side_of(event)
            if side is not None:
                target["side"] = side
            node = {
                "type": event["type"],
                "priority": event["priority"],
                "targets": [target],
                "realizations": [{"index": r["index"], "label": r["label"]} for r in event["realizations"]],
            }
            if event["type"] == "GeneChoice":
                node["library"] = event["class"]
                templates = library.setdefault(event["class"], {})
                for realization in event["realizations"]:
                    templates[realization["label"]] = {"sequence": realization["sequence"]}
            structure_events[copy] = node
            if parents_of[source]:
                conditioning[copy] = sorted(renamed(parent) for parent in parents_of[source])

            # The parameters: the parents in the order the marginals file names them, then the
            # event's own axis, or a Markov chain's two, from then to.
            block = marginals[source]
            if sorted(block["parents"]) != sorted(parents_of[source]):
                raise SourceError(f"the marginals of {source} name other parents than its edges")
            labels = [r["label"] for r in event["realizations"]]
            axes = [{"event": renamed(parent), "labels": [r["label"] for r in events[parent]["realizations"]]}
                    for parent in block["parents"]]
            own = [len(labels)] * (2 if event["type"] == "DinucMarkov" else 1)
            axes += [{"event": copy, "labels": labels} for _ in own]
            sizes = [len(events[parent]["realizations"]) for parent in block["parents"]]
            expected_rows = 1
            for size in sizes:
                expected_rows *= size
            row_length = own[0] * (own[1] if len(own) == 2 else 1)
            if len(block["rows"]) != expected_rows or any(len(row) != row_length for row in block["rows"].values()):
                raise SourceError(f"the marginals of {source} do not fit its realizations")
            nodes[copy] = {"axes": axes, "values": nested(block["rows"], sizes, own)}

    return {
        "$schema": SCHEMA_URI,
        "schema_version": SCHEMA_VERSION,
        "structure": {
            "layout": {"type": "Linear", "segments": SEGMENTS},
            "events": structure_events,
            "conditioning": conditioning,
        },
        "library": library,
        "parameters": {"nodes": nodes},
        "error_rate": {"type": "SingleErrorRate", "rate": error_rate},
    }


def dumps(value, indent=0):
    """The document as json.dumps(indent=2, sort_keys=True) writes it, but for a list of numbers,
    a row of the parameters, which stays on one line: one value per line is hard to read."""
    pad = " " * indent
    inner = " " * (indent + 2)
    if isinstance(value, dict):
        if not value:
            return "{}"
        members = [f"{inner}{json.dumps(key)}: {dumps(value[key], indent + 2)}" for key in sorted(value)]
        return "{\n" + ",\n".join(members) + "\n" + pad + "}"
    if isinstance(value, list):
        if all(isinstance(item, (int, float)) for item in value):
            return json.dumps(value)
        return "[\n" + ",\n".join(inner + dumps(item, indent + 2) for item in value) + "\n" + pad + "]"
    return json.dumps(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source", type=pathlib.Path, default=DEFAULT_SOURCE,
                        help="the directory of TRB's model_parms.txt and model_marginals.txt")
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT, help="the document to write")
    args = parser.parse_args()

    try:
        events, edges, error_rate = read_parms(args.source / "model_parms.txt")
        marginals = read_marginals(args.source / "model_marginals.txt")
        document = build_document(events, edges, error_rate, marginals)
    except SourceError as error:
        sys.exit(f"make_tandem_d_model.py: {error}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(HEADER + dumps(document) + "\n")


if __name__ == "__main__":
    main()
