#!/usr/bin/env python3
"""The cost of an inference on the tandem TRB model, against the TRB model it is built from
(capability plan T5).

Both models infer one iteration, through the CLI, on matched reads: the tandem model on its frozen
reads (tst/test_data/tandem_d/reads), TRB on as many reads generated and aligned the same way
(make_tandem_d_reads.py --model trb). Each starts from the model that generated its reads, with the
CLI's default inference settings. For each, the script reports:

  - the junction fold's build time, from the CLI's "Initializing probability bounds" to its
    "Initialization of probability bounds over";
  - the iteration's wall time and the CPU time of the whole run;
  - per read, the scenarios explored and the time, from the inference log.

    compare_cost.py [--threads N] [--igor BIN]

Timings depend on the machine; the ratios are what the plan records.
"""

import argparse
import os
import pathlib
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
FROZEN = ROOT / "tst" / "test_data" / "tandem_d" / "reads"
MODELS = {
    "tandem": (ROOT / "tst" / "test_data" / "tandem_d" / "trb_tandem_d_model.json",) * 2,
    "trb": (ROOT / "models" / "human" / "tcr_beta" / "models" / "model_parms.txt",
            ROOT / "models" / "human" / "tcr_beta" / "models" / "model_marginals.txt"),
}
BATCH = "cost"


def igor(binary, workdir, *args, **kwargs):
    return subprocess.run([str(binary), "-w", str(workdir), *args], check=True, **kwargs)


def timed_inference(binary, workdir, threads):
    """Run one iteration, timing the fold from the CLI's own messages. Returns (fold, wall, cpu)."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    fold_start = fold_end = None
    proc = subprocess.Popen([str(binary), "-w", str(workdir), "-b", BATCH, "-j", str(threads), "infer"],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    seen = ""
    while True:
        chunk = proc.stdout.read(1)
        if not chunk:
            break
        seen = (seen + chunk)[-200:]
        if fold_start is None and "Initializing probability bounds" in seen:
            fold_start = time.monotonic()
        if fold_end is None and "Initialization of probability bounds over" in seen:
            fold_end = time.monotonic()
    if proc.wait() != 0:
        sys.exit(f"compare_cost.py: the inference failed in {workdir}")
    wall = time.monotonic() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return fold_end - fold_start, wall, cpu


def per_read(log):
    """Scenarios explored and time, per read, from an inference log."""
    rows = [line.split(";") for line in log.read_text().splitlines()[1:]]
    return [int(row[8]) for row in rows], [float(row[10]) for row in rows]


def measure(name, reads, binary, threads):
    parms, marginals = MODELS[name]
    with tempfile.TemporaryDirectory() as workdir:
        work = pathlib.Path(workdir)
        igor(binary, work, "init", stdout=subprocess.DEVNULL)
        for key, value in [("model.source", "custom"), ("model.parms", parms), ("model.marginals", marginals),
                           ("infer.iterations", 1)]:
            igor(binary, work, "config", "set", key, str(value), stdout=subprocess.DEVNULL)
        aligns = work / "aligns"
        aligns.mkdir(exist_ok=True)
        shutil.copyfile(reads / "sequences.csv", aligns / f"{BATCH}_indexed_sequences.csv")
        for gene in "VDJ":
            shutil.copyfile(reads / f"{gene}_alignments.csv", aligns / f"{BATCH}_{gene}_alignments.csv")
        fold, wall, cpu = timed_inference(binary, work, threads)
        scenarios, seconds = per_read(work / f"{BATCH}_inference" / "inference_logs.txt")
    return {"reads": len(scenarios), "fold": fold, "wall": wall, "cpu": cpu,
            "scenarios mean": statistics.mean(scenarios), "scenarios median": statistics.median(scenarios),
            "time per read": statistics.mean(seconds)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--igor", type=pathlib.Path, default=ROOT / "build" / "bin" / "igor")
    args = parser.parse_args()

    count = len(FROZEN.joinpath("sequences.csv").read_text().splitlines()) - 1
    with tempfile.TemporaryDirectory() as trb_reads:
        subprocess.run([sys.executable, "-I", str(ROOT / "scripts" / "tandem_d" / "make_tandem_d_reads.py"),
                        "--model", "trb", "--count", str(count), "--output", trb_reads, "--igor", str(args.igor)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        results = {"trb": measure("trb", pathlib.Path(trb_reads), args.igor, args.threads),
                   "tandem": measure("tandem", FROZEN, args.igor, args.threads)}

    print(f"One iteration on {count} reads each, {args.threads} threads, {os.cpu_count()} cores")
    print(f"{'':22}{'TRB':>14}{'tandem':>14}{'ratio':>10}")
    for key, unit in [("fold", "s"), ("wall", "s"), ("cpu", "s"), ("time per read", "s"),
                      ("scenarios mean", ""), ("scenarios median", "")]:
        trb, tandem = results["trb"][key], results["tandem"][key]
        print(f"{key + (' (' + unit + ')' if unit else ''):22}{trb:>14.4g}{tandem:>14.4g}{tandem / trb:>10.1f}")


if __name__ == "__main__":
    main()
