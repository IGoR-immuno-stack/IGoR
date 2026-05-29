#!/usr/bin/env python3


import argparse
import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path


def load_event_labels(event_path: str) -> list[str]:
    path = Path(event_path)
    if not path.exists():
        raise FileNotFoundError(f"GITHUB_EVENT_PATH not found: {event_path}")

    payload = json.loads(path.read_text(encoding="utf-8"))
    pr = payload.get("pull_request", {})
    labels = pr.get("labels", [])
    return [entry["name"] for entry in labels if isinstance(entry, dict) and "name" in entry]


def load_available_test_tasks(pixi_toml_path: str) -> set[str]:
    path = Path(pixi_toml_path)
    if not path.exists():
        raise FileNotFoundError(f"pixi.toml not found: {pixi_toml_path}")

    data = tomllib.loads(path.read_text(encoding="utf-8"))
    tasks = data.get("tasks", {})

    available = set()
    for task_name in tasks:
        if not task_name.startswith("test_"):
            continue
        if task_name == "test":
            continue
        available.add(task_name)

    return available


def resolve_requested_tasks(labels: list[str], available_tasks: set[str]) -> tuple[list[str], list[str]]:
    requested = sorted(
        {
            label.removeprefix("test::")
            for label in labels
            if label.startswith("test::") and label != "test::"
        }
    )

    valid_tasks = []
    unknown_tasks = []

    for name in requested:
        task = f"test_{name}"
        if task in available_tasks:
            valid_tasks.append(task)
        else:
            unknown_tasks.append(task)

    return valid_tasks, unknown_tasks


def print_plan(
    labels: list[str],
    available_tasks: set[str],
    valid_tasks: list[str],
    unknown_tasks: list[str],
    default_task: str,
) -> None:
    print("=== PR labels ===")
    if labels:
        for label in labels:
            print(f"- {label}")
    else:
        print("(none)")

    print("\n=== Available label-selectable Pixi tasks ===")
    if available_tasks:
        for task in sorted(available_tasks):
            print(f"- {task}")
    else:
        print("(none)")

    print("\n=== Decision ===")
    if valid_tasks:
        print("Selected Pixi tasks:")
        for task in valid_tasks:
            print(f"- pixi run {task}")
    else:
        print("- no valid test::XXX label found")
        print(f"- pixi run {default_task}")

    if unknown_tasks:
        print("\n=== Unknown requested tasks ===")
        for task in unknown_tasks:
            print(f"- ERROR: requested task '{task}' does not exist in [tasks] of pixi.toml")


def execute_plan(valid_tasks: list[str], default_task: str) -> int:
    commands = valid_tasks if valid_tasks else [default_task]

    for task in commands:
        print(f"::group::pixi run {task}")
        try:
            result = subprocess.run(["pixi", "run", task], check=False)
        except FileNotFoundError:
            print("ERROR: `pixi` executable not found in PATH.", file=sys.stderr)
            print("::endgroup::")
            return 2

        print("::endgroup::")

        if result.returncode != 0:
            print(
                f"ERROR: pixi task '{task}' failed with exit code {result.returncode}.",
                file=sys.stderr,
            )
            return result.returncode

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolve PR labels into Pixi test tasks.")
    parser.add_argument("--dry-run", action="store_true", help="Only print the planned commands.")
    parser.add_argument("--execute", action="store_true", help="Execute the planned commands.")
    parser.add_argument(
        "--default-task",
        default="test",
        help="Fallback Pixi task if no valid label matches.",
    )
    parser.add_argument(
        "--pixi-toml",
        default="pixi.toml",
        help="Path to pixi.toml",
    )
    args = parser.parse_args()

    if args.dry_run == args.execute:
        print("Choose exactly one of --dry-run or --execute", file=sys.stderr)
        return 2

    event_path = os.environ.get("GITHUB_EVENT_PATH")
    if not event_path:
        print("GITHUB_EVENT_PATH is not defined", file=sys.stderr)
        return 2

    labels = load_event_labels(event_path)
    available_tasks = load_available_test_tasks(args.pixi_toml)
    valid_tasks, unknown_tasks = resolve_requested_tasks(labels, available_tasks)

    print_plan(labels, available_tasks, valid_tasks, unknown_tasks, args.default_task)

    if args.dry_run:
        print("\n=== Mode ===")
        print("- dry-run only, nothing executed")
        return 0

    print("\n=== Mode ===")
    print("- execute")
    return execute_plan(valid_tasks, args.default_task)


if __name__ == "__main__":
    raise SystemExit(main())
