#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def load_event_labels(event_path: str) -> list[str]:
    path = Path(event_path)
    if not path.exists():
        raise FileNotFoundError(f"GITHUB_EVENT_PATH not found: {event_path}")

    payload = json.loads(path.read_text(encoding="utf-8"))
    pr = payload.get("pull_request", {})
    labels = pr.get("labels", [])
    return [entry["name"] for entry in labels if isinstance(entry, dict) and "name" in entry]


def compute_tasks(labels: list[str]) -> list[str]:
    selected = sorted(
        {
            label.removeprefix("test::")
            for label in labels
            if label.startswith("test::") and label != "test::"
        }
    )
    return [f"test_{name}" for name in selected]


def print_plan(labels: list[str], tasks: list[str], default_task: str) -> None:
    print("=== PR labels ===")
    if labels:
        for label in labels:
            print(f"- {label}")
    else:
        print("(none)")

    print("\n=== Decision ===")
    if tasks:
        print("Selected Pixi tasks:")
        for task in tasks:
            print(f"- pixi run {task}")
    else:
        print("- no test::XXX label found")
        print(f"- pixi run {default_task}")


def execute_plan(tasks: list[str], default_task: str) -> int:
    commands = tasks if tasks else [default_task]

    for task in commands:
        print(f"::group::pixi run {task}")
        result = subprocess.run(["pixi", "run", task], check=False)
        print("::endgroup::")
        if result.returncode != 0:
            return result.returncode
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolve PR labels into Pixi test tasks.")
    parser.add_argument("--dry-run", action="store_true", help="Only print the planned commands.")
    parser.add_argument("--execute", action="store_true", help="Execute the planned commands.")
    parser.add_argument("--default-task", default="test", help="Fallback Pixi task if no label matches.")
    args = parser.parse_args()

    if args.dry_run == args.execute:
        print("Choose exactly one of --dry-run or --execute", file=sys.stderr)
        return 2

    event_path = os.environ.get("GITHUB_EVENT_PATH")
    if not event_path:
        print("GITHUB_EVENT_PATH is not defined", file=sys.stderr)
        return 2

    labels = load_event_labels(event_path)
    tasks = compute_tasks(labels)

    print_plan(labels, tasks, args.default_task)

    if args.dry_run:
        print("\n=== Mode ===")
        print("- dry-run only, nothing executed")
        return 0

    print("\n=== Mode ===")
    print("- execute")
    return execute_plan(tasks, args.default_task)


if __name__ == "__main__":
    raise SystemExit(main())
