#!/usr/bin/env python3
"""The badges of README.md, drawn here as SVG files rather than by an external service.

Reads build/coverage/summary.json (pixi run coverage_summary) and counts the unit tests ctest
registers in build/ (ctest -L unit -N), then writes into build/badges/:
  - coverage.svg: the line coverage of src/ by the unit tests, coloured with the thresholds of
    the pull request comment (coverage_markdown.py);
  - unit-tests.svg: the number of unit tests.
CI draws them once the unit tests have passed on a push to develop, and publishes them on the
badges branch, from which README.md shows them.

Usage:
    python3 scripts/tests/make_badges.py    (pixi run badges)
"""

import json
import os
import re
import subprocess

from coverage_markdown import GREEN, ORANGE

COLOURS = {"green": "#4c1", "orange": "#fe7d37", "red": "#e05d44"}

# Approximate advance widths, in pixels, of 11 px Verdana: enough to size the two boxes.
NARROW = set("fijlrt.,:;!|' ")
WIDE = set("mwMW%")


def text_width(text):
    return sum(3.8 if c in NARROW else 10.0 if c in WIDE else 7.0 for c in text)


def badge(label, value, colour):
    left = round(text_width(label) + 12)
    right = round(text_width(value) + 12)
    width = left + right
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="20" role="img" aria-label="{label}: {value}">
  <title>{label}: {value}</title>
  <linearGradient id="s" x2="0" y2="100%">
    <stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
    <stop offset="1" stop-opacity=".1"/>
  </linearGradient>
  <clipPath id="r"><rect width="{width}" height="20" rx="3" fill="#fff"/></clipPath>
  <g clip-path="url(#r)">
    <rect width="{left}" height="20" fill="#555"/>
    <rect x="{left}" width="{right}" height="20" fill="{colour}"/>
    <rect width="{width}" height="20" fill="url(#s)"/>
  </g>
  <g fill="#fff" text-anchor="middle" font-family="Verdana,Geneva,DejaVu Sans,sans-serif" font-size="11">
    <text x="{left / 2}" y="15" fill="#010101" fill-opacity=".3">{label}</text>
    <text x="{left / 2}" y="14">{label}</text>
    <text x="{left + right / 2}" y="15" fill="#010101" fill-opacity=".3">{value}</text>
    <text x="{left + right / 2}" y="14">{value}</text>
  </g>
</svg>
"""


def unit_test_count():
    listing = subprocess.run(
        ["ctest", "--test-dir", "build", "-L", "unit", "-N"],
        check=True, capture_output=True, text=True,
    ).stdout
    match = re.search(r"^Total Tests: (\d+)$", listing, re.MULTILINE)
    if not match:
        raise SystemExit("ctest listed no unit test in build/")
    return int(match.group(1))


def main():
    with open(os.path.join("build", "coverage", "summary.json"), encoding="utf-8") as handle:
        percent = json.load(handle)["line_percent"]
    colour = "green" if percent >= GREEN else "orange" if percent >= ORANGE else "red"
    out = os.path.join("build", "badges")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "coverage.svg"), "w", encoding="utf-8") as handle:
        handle.write(badge("coverage", f"{percent:.1f}%", COLOURS[colour]))
    with open(os.path.join(out, "unit-tests.svg"), "w", encoding="utf-8") as handle:
        handle.write(badge("unit tests", f"{unit_test_count()} passing", COLOURS["green"]))


if __name__ == "__main__":
    main()
