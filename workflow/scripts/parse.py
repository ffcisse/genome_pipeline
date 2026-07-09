#!/usr/bin/env python
"""Stage 2 -- Parse one genome's protein + CDS FASTA into canonical tables.

Mirrors the exploratory notebooks' "parse" step exactly: protein sequences
get their terminal stop codon(s) stripped and length recorded (same as
01_proteome_overview_FINAL.ipynb section 2.1-2.2); CDS sequences get
uppercased and length recorded (same as dna_seq_analysis.ipynb). Neither
applies the full seq_utils.clean_sequence (stripping X/non-standard
residues) -- that happens per-call in property calculations, a later stage,
same as in the notebooks.
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_utils import iter_fasta_records


def parse_protein(path, genome):
    rows = []
    for rec in iter_fasta_records(path):
        seq = str(rec.seq).rstrip("*")
        rows.append(
            {
                "genome": genome,
                "protein_id": rec.id,
                "sequence": seq,
                "length": len(seq),
            }
        )
    return pd.DataFrame(rows, columns=["genome", "protein_id", "sequence", "length"])


def parse_cds(path, genome):
    rows = []
    for rec in iter_fasta_records(path):
        seq = str(rec.seq).upper()
        rows.append(
            {
                "genome": genome,
                "cds_id": rec.id,
                "sequence": seq,
                "length": len(seq),
            }
        )
    return pd.DataFrame(rows, columns=["genome", "cds_id", "sequence", "length"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--protein-out", required=True)
    parser.add_argument("--cds-out", required=True)
    parser.add_argument("--done", required=True)
    args = parser.parse_args()

    protein_df = parse_protein(args.protein, args.genome)
    cds_df = parse_cds(args.cds, args.genome)

    os.makedirs(os.path.dirname(args.protein_out), exist_ok=True)
    os.makedirs(os.path.dirname(args.cds_out), exist_ok=True)
    protein_df.to_csv(args.protein_out, index=False)
    cds_df.to_csv(args.cds_out, index=False)

    print(f"{args.genome}: wrote {len(protein_df)} proteins -> {args.protein_out}")
    print(f"{args.genome}: wrote {len(cds_df)} CDS -> {args.cds_out}")

    os.makedirs(os.path.dirname(args.done), exist_ok=True)
    with open(args.done, "w") as fh:
        fh.write(f"{args.genome}: {len(protein_df)} proteins, {len(cds_df)} CDS\n")


if __name__ == "__main__":
    main()
