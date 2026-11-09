#!/usr/bin/env python3
"""Generate reads from the tandem TRB model and freeze them with their alignments (capability plan
T5, step T5b).

The reads milestone 1 infers from, and the matched TRB reads its cost is compared on. For the
model given, through the IGoR CLI, on one thread:

  1. generate `--count` reads with the legacy generator, seeded, errors included;
  2. align them against TRB's genomic templates, with the CLI's default alignment settings;
  3. keep, for each read, the alignments an inference reads: those scoring within 55 of the read's
     best for V, 35 for D and 10 for J, as igor_legacy_main() reads them
     (read_alignments_seq_csv_score_range()). An inference on the frozen files reads the same
     alignments as one on the full ones, and the D alignments, most of them short and many, are
     what makes the full files large;
  4. write the reads, their realizations (the scenario that generated each) and the alignments.

    make_tandem_d_reads.py [--model tandem|trb] [--count N] [--seed S] [--output DIR] [--igor BIN]

The tandem reads are frozen in tst/test_data/tandem_d/reads, the default output; the TRB ones are
generated where the cost comparison needs them. Generation and alignment are deterministic on one
thread, so running the script again on the same build gives the same files.
"""

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
GENOME = ROOT / "models" / "human" / "tcr_beta" / "ref_genome"
MODELS = {
    "tandem": (ROOT / "tst" / "test_data" / "tandem_d" / "trb_tandem_d_model.json",) * 2,
    "trb": (ROOT / "models" / "human" / "tcr_beta" / "models" / "model_parms.txt",
            ROOT / "models" / "human" / "tcr_beta" / "models" / "model_marginals.txt"),
}
DEFAULT_OUTPUT = ROOT / "tst" / "test_data" / "tandem_d" / "reads"

# The score range igor_legacy_main() reads each gene's alignments with, below the read's best.
SCORE_RANGES = {"V": 55.0, "D": 35.0, "J": 10.0}
BATCH = "reads"


def igor(binary, workdir, *args):
    subprocess.run([str(binary), "-w", str(workdir), *args], check=True, stdout=subprocess.DEVNULL)


def kept_alignments(path, score_range):
    """The header and the rows of an alignments file that score within `score_range` of their
    read's best, in the file's order."""
    lines = path.read_text().splitlines()
    header, rows = lines[0], [line.split(";") for line in lines[1:]]
    best = {}
    for row in rows:
        best[row[0]] = max(best.get(row[0], -1.0), float(row[2]))
    return header, [";".join(row) for row in rows if float(row[2]) >= best[row[0]] - score_range]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--model", choices=sorted(MODELS), default="tandem")
    parser.add_argument("--count", type=int, default=200, help="the number of reads")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=pathlib.Path, help="the directory to write (default: the frozen "
                        "tandem reads, for the tandem model)")
    parser.add_argument("--igor", type=pathlib.Path, default=ROOT / "build" / "bin" / "igor")
    args = parser.parse_args()
    if args.output is None:
        if args.model != "tandem":
            sys.exit("make_tandem_d_reads.py: only the tandem reads are frozen; give --output for the others")
        args.output = DEFAULT_OUTPUT

    parms, marginals = MODELS[args.model]
    with tempfile.TemporaryDirectory() as workdir:
        work = pathlib.Path(workdir)
        igor(args.igor, work, "init")
        for key, value in [("model.source", "custom"), ("model.parms", parms), ("model.marginals", marginals),
                           ("generate.seed", args.seed), ("generate.fast", "false"), ("generate.error", "true"),
                           ("genomic.V", GENOME / "genomicVs.fasta"), ("genomic.D", GENOME / "genomicDs.fasta"),
                           ("genomic.J", GENOME / "genomicJs.fasta")]:
            igor(args.igor, work, "config", "set", key, str(value))
        igor(args.igor, work, "-b", BATCH, "-j", "1", "generate", str(args.count))
        generated = work / f"{BATCH}_generated"
        igor(args.igor, work, "-b", BATCH, "import-seqs", str(generated / "generated_seqs_werr.csv"))
        for gene in SCORE_RANGES:
            igor(args.igor, work, "-b", BATCH, "-j", "1", "align", "--gene", gene)

        args.output.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(generated / "generated_seqs_werr.csv", args.output / "sequences.csv")
        shutil.copyfile(generated / "generated_realizations_werr.csv", args.output / "realizations.csv")
        for gene, score_range in SCORE_RANGES.items():
            header, rows = kept_alignments(work / "aligns" / f"{BATCH}_{gene}_alignments.csv", score_range)
            (args.output / f"{gene}_alignments.csv").write_text("\n".join([header, *rows]) + "\n")
            print(f"{gene}: {len(rows)} alignments kept for {args.count} reads")


if __name__ == "__main__":
    main()
