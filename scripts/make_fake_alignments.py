#!/usr/bin/env python3
"""Generate fake V/J alignment files using model parameters.

The model parameters file is the authoritative source for gene names and
template lengths. We parse it to ensure the alignment metadata matches the
model used for generation and inference, without hard-coding lengths or names.
"""
import argparse
import csv
from pathlib import Path
from typing import List, Tuple


def parse_model_genes(model_params_path: Path) -> Tuple[Tuple[str, int], Tuple[str, int]]:
    """Extract the first V and J gene names and lengths from model params.

    The model_parms file encodes the genomic templates for each gene choice.
    We rely on it to obtain the exact gene identifiers and sequence lengths
    so the fake alignment records are consistent with the model structure
    passed to IGoR during inference.
    """
    current_gene_class = None
    v_gene = None
    j_gene = None

    with model_params_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#GeneChoice;"):
                parts = line.split(";")
                if len(parts) >= 2:
                    current_gene_class = parts[1]
                else:
                    current_gene_class = None
                continue
            if not line.startswith("%"):
                continue

            if current_gene_class not in {"V_gene", "J_gene"}:
                continue

            parts = line.split(";")
            if len(parts) < 3:
                continue
            gene_name = parts[0][1:]
            gene_seq = parts[1]

            if gene_name == "None" or gene_seq == "":
                continue

            if current_gene_class == "V_gene" and v_gene is None:
                v_gene = (gene_name, len(gene_seq))
            elif current_gene_class == "J_gene" and j_gene is None:
                j_gene = (gene_name, len(gene_seq))

    if v_gene is None or j_gene is None:
        raise ValueError("Could not find V or J gene in model params file.")

    return v_gene, j_gene


def read_indexed_sequences(indexed_csv: Path) -> List[Tuple[int, str]]:
    sequences: List[Tuple[int, str]] = []
    with indexed_csv.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, None)
        if not header or header[0] != "seq_index":
            raise ValueError(f"Unexpected header in {indexed_csv}")
        for row in reader:
            if not row:
                continue
            seq_index = int(row[0])
            seq = row[1].strip().upper()
            sequences.append((seq_index, seq))
    return sequences


def write_alignments(path: Path, rows: List[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow([
            "seq_index",
            "gene_name",
            "score",
            "offset",
            "insertions",
            "deletions",
            "mismatches",
            "length",
            "5_p_align_offset",
            "3_p_align_offset",
        ])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build fake IGoR alignment files.")
    parser.add_argument("--model-params", required=True, type=Path)
    parser.add_argument("--indexed-seqs", required=True, type=Path)
    parser.add_argument("--out-v", required=True, type=Path)
    parser.add_argument("--out-j", required=True, type=Path)
    parser.add_argument("--v-score", default=100, type=int)
    parser.add_argument("--j-score", default=100, type=int)
    args = parser.parse_args()

    v_gene, j_gene = parse_model_genes(args.model_params)
    sequences = read_indexed_sequences(args.indexed_seqs)

    v_name, v_len = v_gene
    j_name, j_len = j_gene

    v_rows: List[List[str]] = []
    j_rows: List[List[str]] = []

    for seq_index, seq in sequences:
        seq_len = len(seq)

        v_rows.append([
            str(seq_index),
            v_name,
            str(args.v_score),
            "0",
            "{}",
            "{}",
            "{}",
            str(v_len),
            "0",
            str(v_len - 1),
        ])

        j_offset = seq_len - j_len
        j_rows.append([
            str(seq_index),
            j_name,
            str(args.j_score),
            str(j_offset),
            "{}",
            "{}",
            "{}",
            str(j_len),
            str(j_offset),
            str(seq_len - 1),
        ])


    write_alignments(args.out_v, v_rows)
    write_alignments(args.out_j, j_rows)


if __name__ == "__main__":
    main()
