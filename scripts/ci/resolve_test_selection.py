#!/usr/bin/env python3
"""Resolve GitHub PR labels into Pixi test tasks.

This script reads pull request labels from the GitHub event payload and maps
labels of the form ``test::XXX`` to Pixi tasks named ``test_XXX``. Available
tasks are discovered from ``pixi.toml`` so that unknown labels can be reported
cleanly without crashing the workflow.

The script supports two modes:
- ``--dry-run``: print the execution plan without running anything.
- ``--execute``: run the selected Pixi tasks, or a default fallback task if
  no valid label is found.
"""

import argparse
import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path


def load_event_labels(event_path: str) -> list[str]:
    """Load pull request labels from the GitHub event payload.

    Parameters
    ----------
    event_path:
        Path to the JSON file pointed to by ``GITHUB_EVENT_PATH``.

    Returns
    -------
    list[str]
        The list of PR label names found in the payload.

    Raises
    ------
    FileNotFoundError
        If the event payload file does not exist.
    """
    path = Path(event_path)
    if not path.exists():
        raise FileNotFoundError(f"GITHUB_EVENT_PATH not found: {event_path}")

    payload = json.loads(path.read_text(encoding="utf-8"))
    pr = payload.get("pull_request", {})
    labels = pr.get("labels", [])
    return [entry["name"] for entry in labels if isinstance(entry, dict) and "name" in entry]


def load_available_test_tasks(pixi_toml_path: str) -> set[str]:
    """Read label-selectable test tasks from ``pixi.toml``.

    Only tasks defined in the top-level ``[tasks]`` section and whose names
    start with ``test_`` are considered selectable. The generic ``test`` task
    is excluded because it is used as the default fallback task rather than as
    a label-mapped task.

    Parameters
    ----------
    pixi_toml_path:
        Path to the ``pixi.toml`` file.

    Returns
    -------
    set[str]
        The set of available Pixi task names matching the ``test_XXX`` pattern.

    Raises
    ------
    FileNotFoundError
        If the Pixi manifest file does not exist.
    """
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
    """Split requested label-mapped tasks into valid and unknown tasks.

    Labels of the form ``test::XXX`` are converted into Pixi task names
    ``test_XXX``. The function then separates tasks that exist in ``pixi.toml``
    from those that do not.

    Parameters
    ----------
    labels:
        Labels attached to the pull request.
    available_tasks:
        Set of task names discovered in ``pixi.toml``.

    Returns
    -------
    tuple[list[str], list[str]]
        A pair ``(valid_tasks, unknown_tasks)`` where both lists contain Pixi
        task names.
    """
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
    """Print a human-readable execution plan.

    The output includes the raw PR labels, the available Pixi test tasks, the
    tasks selected from labels, and any unknown tasks requested by labels.

    Parameters
    ----------
    labels:
        Raw labels found on the pull request.
    available_tasks:
        Tasks discovered in ``pixi.toml`` and eligible for label selection.
    valid_tasks:
        Tasks that were requested and found in ``pixi.toml``.
    unknown_tasks:
        Tasks requested by labels but not found in ``pixi.toml``.
    default_task:
        Fallback Pixi task to use when no valid label-mapped task is found.
    """
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
    """Execute the selected Pixi tasks or the default fallback task.

    Unknown tasks are not executed because they are filtered out earlier when
    the script compares requested labels against the tasks defined in
    ``pixi.toml``. If a valid Pixi task fails, the function returns its exit
    code. If Pixi itself cannot be found, the function returns ``2``.

    Parameters
    ----------
    valid_tasks:
        The list of validated Pixi tasks to run.
    default_task:
        The fallback task to run when no valid label-mapped task is available.

    Returns
    -------
    int
        Exit code ``0`` on success, ``2`` if ``pixi`` is missing, or the failing
        task's exit code otherwise.
    """
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
    """Parse arguments, build the execution plan, and run it if requested.

    The script requires exactly one of ``--dry-run`` or ``--execute``. Pull
    request labels are read from ``GITHUB_EVENT_PATH`` and available test tasks
    are discovered from ``pixi.toml``.

    Returns
    -------
    int
        Process exit code suitable for use in CI.
    """
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
