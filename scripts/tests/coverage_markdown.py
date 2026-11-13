#!/usr/bin/env python3
"""The unit tests' coverage figures as Markdown, for CI's run summary and pull request comment.

Reads, in the coverage directory (default build/coverage):
  - summary.json, gcovr's summary of src/ (pixi run coverage_summary);
  - diff.json when it exists, diff-cover's report on the lines a branch changes against its base
    (pixi run coverage_diff).
Writes comment.md there. Both tasks run this script last, so comment.md always holds what the
directory has.

The colours follow the thresholds the pull request comment was configured with (80 % green,
60 % orange), which gcovr's HTML report uses too.

Usage:
    python3 scripts/tests/coverage_markdown.py [coverage directory]
"""

import json
import os
import sys

GREEN = 80.0
ORANGE = 60.0
# Rows of the per-file table and characters of a line list: a GitHub comment holds at most
# 65536 characters.
MAX_FILES = 40
MAX_LINES_TEXT = 200


def mark(percent):
    if percent >= GREEN:
        return "🟢"
    if percent >= ORANGE:
        return "🟠"
    return "🔴"


def line_ranges(lines):
    """[3, 4, 5, 9] -> "3-5, 9", cut at MAX_LINES_TEXT characters."""
    parts = []
    start = previous = None
    for line in sorted(lines):
        if start is None:
            start = previous = line
        elif line == previous + 1:
            previous = line
        else:
            parts.append(f"{start}-{previous}" if previous > start else f"{start}")
            start = previous = line
    if start is not None:
        parts.append(f"{start}-{previous}" if previous > start else f"{start}")
    text = ", ".join(parts)
    return text if len(text) <= MAX_LINES_TEXT else text[:MAX_LINES_TEXT].rsplit(", ", 1)[0] + ", …"


def summary_section(summary):
    rows = [
        "## Coverage of src/ by the unit tests",
        "",
        "| | Covered | Total | |",
        "|---|---:|---:|---:|",
    ]
    for name, key in (("Lines", "line"), ("Functions", "function"), ("Branches", "branch")):
        percent = summary[f"{key}_percent"]
        rows.append(
            f"| {name} | {summary[f'{key}_covered']} | {summary[f'{key}_total']} "
            f"| {mark(percent)} {percent:.1f}% |"
        )
    return rows


def diff_section(diff):
    rows = ["", f"### The lines this branch changes ({diff['diff_name']})", ""]
    total = diff["total_num_lines"]
    if total == 0:
        rows.append("No changed line in src/ is executable code.")
        return rows
    covered = total - diff["total_num_violations"]
    percent = 100.0 * covered / total
    rows.append(
        f"{mark(percent)} **{percent:.1f}%**: the unit tests run {covered} of the {total} "
        "changed executable lines in src/."
    )
    missed = sorted(
        ((path, stats) for path, stats in diff["src_stats"].items() if stats["violation_lines"]),
        key=lambda item: (-len(item[1]["violation_lines"]), item[0]),
    )
    if missed:
        rows += [
            "",
            "<details><summary>Changed lines the unit tests do not run, by file</summary>",
            "",
            "| File | Changed lines run | Lines not run |",
            "|---|---:|---|",
        ]
        for path, stats in missed[:MAX_FILES]:
            percent = stats["percent_covered"]
            rows.append(
                f"| `{path}` | {mark(percent)} {percent:.1f}% "
                f"| {line_ranges(stats['violation_lines'])} |"
            )
        if len(missed) > MAX_FILES:
            rows.append(f"| and {len(missed) - MAX_FILES} more files | | |")
        rows += ["", "</details>"]
    return rows


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else os.path.join("build", "coverage")
    with open(os.path.join(directory, "summary.json"), encoding="utf-8") as handle:
        rows = summary_section(json.load(handle))
    diff_path = os.path.join(directory, "diff.json")
    if os.path.exists(diff_path):
        with open(diff_path, encoding="utf-8") as handle:
            rows += diff_section(json.load(handle))
    with open(os.path.join(directory, "comment.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
