#!/usr/bin/env python3
"""Per-function coverage for the events being migrated.

The generic rewrite replaces `iterate()` and `initialize_event()` one event at a time, and
each step's characterization suite (the *a* commits of docs/ITERATE_GENERIC_REWRITE_PLAN.md
section 6) is only evidence for the branches it actually reaches. This reports, for a chosen
set of tests, how much of those two functions they execute -- so an uncovered branch is a
visible gap before the collapse rather than a surprise after it.

Usage:
    pixi run build_coverage
    python3 scripts/tests/coverage_report.py -f '[insertion][iterate]'
    python3 scripts/tests/coverage_report.py                # every non-hidden test case

Counters are reset before the run, so the numbers describe exactly the tests selected by -f.
"""

import argparse
import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUILD = os.path.join(ROOT, "build_cov")
OBJ_DIR = os.path.join(BUILD, "src", "igor", "Core", "CMakeFiles", "Core.dir")
DEFAULT_SOURCES = ["Genechoice.cpp", "Deletion.cpp", "Insertion.cpp", "Dinuclmarkov.cpp"]
DEFAULT_PATTERN = r"::(iterate|initialize_event)\("


def find_gcov():
    """The toolchain's gcov, not the host's: the .gcno format is version-locked."""
    candidates = glob.glob(os.path.join(ROOT, ".pixi", "envs", "*", "bin", "*-gcov"))
    if candidates:
        return sorted(candidates)[0]
    found = shutil.which("gcov")
    if not found:
        sys.exit("No gcov found. Install one, or build the pixi environment.")
    return found


def run_tests(test_filter):
    binary = os.path.join(BUILD, "bin", "igor_tests")
    if not os.path.exists(binary):
        sys.exit("No coverage build found. Run: pixi run build_coverage")
    for gcda in glob.glob(os.path.join(BUILD, "**", "*.gcda"), recursive=True):
        os.remove(gcda)
    command = [binary] + ([test_filter] if test_filter else [])
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode not in (0, 1):   # 1 = some test failed, still worth reporting
        sys.exit("Test binary exited with %d" % result.returncode)


def gcov_json(gcov, source, work_dir):
    """Run gcov on one translation unit and return its parsed report, or None."""
    obj = os.path.join(OBJ_DIR, source + ".o")
    if not os.path.exists(obj):
        return None
    subprocess.run([gcov, "--json-format", "-b", "-m", "-o", obj,
                    os.path.join(ROOT, "src", "igor", "Core", source)],
                   cwd=work_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    reports = glob.glob(os.path.join(work_dir, "*.json.gz"))
    if not reports:
        return None
    with gzip.open(reports[0]) as handle:
        parsed = json.load(handle)
    for report in reports:
        os.remove(report)
    return parsed


def summarize(report, source, pattern):
    """One row per matching function: line, branch and block coverage."""
    rows = []
    for entry in report.get("files", []):
        if not entry["file"].endswith(source):
            continue
        lines = entry.get("lines", [])
        for function in entry.get("functions", []):
            name = function.get("demangled_name") or function["name"]
            if not pattern.search(name):
                continue
            start, end = function["start_line"], function["end_line"]
            body = [line for line in lines if start <= line["line_number"] <= end]
            covered = sum(1 for line in body if line["count"] > 0)
            branches = [b for line in body for b in line.get("branches", [])]
            taken = sum(1 for b in branches if b.get("count", 0) > 0)
            rows.append({
                "name": name.split("(")[0],
                "lines": (covered, len(body)),
                "branches": (taken, len(branches)),
                "blocks": (function["blocks_executed"], function["blocks"]),
                "calls": function.get("execution_count", 0),
            })
    return rows


def percent(pair):
    covered, total = pair
    return "  n/a" if total == 0 else "%5.1f%%" % (100.0 * covered / total)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--filter", default="", help="Catch2 test specification to run")
    parser.add_argument("-s", "--source", action="append", default=[],
                        help="Core source to report on (repeatable)")
    parser.add_argument("--pattern", default=DEFAULT_PATTERN,
                        help="regex selecting functions (default: iterate / initialize_event)")
    args = parser.parse_args()

    sources = args.source or DEFAULT_SOURCES
    pattern = re.compile(args.pattern)
    gcov = find_gcov()
    run_tests(args.filter)

    print("Tests run: %s" % (args.filter or "(all non-hidden test cases)"))
    print()
    print("%-44s %8s %8s %8s %8s" % ("FUNCTION", "LINES", "BRANCH", "BLOCKS", "CALLS"))
    print("%-44s %8s %8s %8s %8s" % ("-" * 44, "-" * 8, "-" * 8, "-" * 8, "-" * 8))

    work_dir = tempfile.mkdtemp()
    try:
        for source in sources:
            report = gcov_json(gcov, source, work_dir)
            if report is None:
                continue
            for row in sorted(summarize(report, source, pattern), key=lambda r: r["name"]):
                print("%-44s %8s %8s %8s %8d"
                      % (row["name"][:44], percent(row["lines"]), percent(row["branches"]),
                         percent(row["blocks"]), row["calls"]))
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
