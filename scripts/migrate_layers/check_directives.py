#!/usr/bin/env python3
"""Keep the transitional using-directives honest (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).

A `using namespace igor::<layer>::legacy;` is legitimate only if the file reaches a header of
`src/igor/<Layer>/Legacy/` through its <igor/...> includes (and `using namespace igor::<layer>;`
a header of `src/igor/<Layer>/`): otherwise the namespace may not even
be declared, and the directive is the one line that fails to compile. Each promotion of 1c
shrinks include closures, so run this after every one.

    check_directives.py          report the directives that no longer reach their layer
    check_directives.py --fix    remove them

Run from the repository root. Exit status 1 when something is reported and --fix is not given.
"""

import os
import re
import subprocess
import sys

ROOT = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
os.chdir(ROOT)

INCLUDE_LINE = re.compile(r'^\s*#\s*include\s*[<"]igor/([^>"]+)[>"]', re.M)
DIRECTIVE = re.compile(r"^using namespace igor::([a-z]+)(::legacy)?;\n", re.M)
_cache = {}


def read(path):
    with open(path, encoding="utf-8", errors="surrogateescape") as handle:
        return handle.read()


def closure(path):
    if path in _cache:
        return _cache[path]
    seen, stack = set(), [path]
    while stack:
        current = stack.pop()
        if current in seen or not os.path.exists(current):
            continue
        seen.add(current)
        stack.extend(f"src/igor/{rel}" for rel in INCLUDE_LINE.findall(read(current)))
    _cache[path] = seen
    return seen


def own_layer(path):
    m = re.match(r"src/igor/([A-Za-z]+)/Legacy/", path)
    return m.group(1).lower() if m else None


def main():
    fix = "--fix" in sys.argv
    stale = 0
    for top in ("src", "app", "tst"):
        for dirpath, _, files in os.walk(top):
            for name in files:
                if not name.endswith((".h", ".cpp", ".tpp")):
                    continue
                path = os.path.join(dirpath, name)
                text = read(path)
                reached = set()
                for h in closure(path):
                    m = re.match(r"src/igor/([A-Za-z]+)/(Legacy/)?", h)
                    if m and h != path:
                        reached.add((m.group(1).lower(), "::legacy" if m.group(2) else ""))
                # a file under src/igor/<Layer>/ declares its own namespace
                if own_layer(path):
                    reached.add((own_layer(path), "::legacy"))
                # a legacy namespace nests in its layer's namespace, which it therefore declares
                reached |= {(layer, "") for layer, _ in list(reached)}
                dead = [(layer, leg) for layer, leg in DIRECTIVE.findall(text)
                        if layer not in ("math", "streaming") and (layer, leg) not in reached]
                if not dead:
                    continue
                stale += len(dead)
                for layer, leg in dead:
                    print(f"{path}: igor::{layer}{leg} is not reachable")
                    if fix:
                        text = text.replace(f"using namespace igor::{layer}{leg};\n", "")
                if fix:
                    with open(path, "w", encoding="utf-8", errors="surrogateescape") as handle:
                        handle.write(text)
    if stale and not fix:
        sys.exit(1)
    print(f"check_directives: {stale} stale directive(s){' removed' if fix else ''}")


if __name__ == "__main__":
    main()
