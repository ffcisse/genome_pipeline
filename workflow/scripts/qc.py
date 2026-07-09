#!/usr/bin/env python
"""Stage 1 -- QC one genome's protein + CDS FASTA.

Reports issues the same way the exploratory notebook's "flag issues before
cleaning" step does (count of X residues, stop codons, etc.) rather than
silently fixing them -- cleaning happens downstream, per-property, via
seq_utils.clean_sequence. Only hard-fails on problems that mean a file is
unusable (nothing parsed, duplicate IDs), since a handful of X's or internal
stops is normal in a real proteome.
"""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_utils import VALID_AAS, iter_fasta_records

DNA_BASES = set("ATGC")


def qc_protein(path):
    n = 0
    n_with_x = 0
    n_with_internal_stop = 0
    n_empty = 0
    ids = Counter()
    for rec in iter_fasta_records(path):
        n += 1
        ids[rec.id] += 1
        seq = str(rec.seq)
        if "X" in seq:
            n_with_x += 1
        if "*" in seq.rstrip("*"):
            n_with_internal_stop += 1
        if len(seq.rstrip("*")) == 0:
            n_empty += 1
    dupes = [seq_id for seq_id, count in ids.items() if count > 1]
    return {
        "n_records": n,
        "n_with_X": n_with_x,
        "n_with_internal_stop": n_with_internal_stop,
        "n_empty": n_empty,
        "n_duplicate_ids": len(dupes),
    }


def qc_cds(path):
    n = 0
    n_non_multiple_of_3 = 0
    n_non_atgc = 0
    n_empty = 0
    ids = Counter()
    for rec in iter_fasta_records(path):
        n += 1
        ids[rec.id] += 1
        seq = str(rec.seq).upper()
        if len(seq) == 0:
            n_empty += 1
            continue
        if len(seq) % 3 != 0:
            n_non_multiple_of_3 += 1
        if any(b not in DNA_BASES for b in seq):
            n_non_atgc += 1
    dupes = [seq_id for seq_id, count in ids.items() if count > 1]
    return {
        "n_records": n,
        "n_non_multiple_of_3": n_non_multiple_of_3,
        "n_non_ATGC": n_non_atgc,
        "n_empty": n_empty,
        "n_duplicate_ids": len(dupes),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protein_stats = qc_protein(args.protein)
    cds_stats = qc_cds(args.cds)

    errors = []
    if protein_stats["n_records"] == 0:
        errors.append(f"protein FASTA has zero parseable records: {args.protein}")
    if cds_stats["n_records"] == 0:
        errors.append(f"CDS FASTA has zero parseable records: {args.cds}")
    if protein_stats["n_duplicate_ids"] > 0:
        errors.append(f"protein FASTA has {protein_stats['n_duplicate_ids']} duplicate IDs")
    if cds_stats["n_duplicate_ids"] > 0:
        errors.append(f"CDS FASTA has {cds_stats['n_duplicate_ids']} duplicate IDs")

    lines = [f"QC report: {args.genome}", ""]
    lines.append(f"protein file: {args.protein}")
    for k, v in protein_stats.items():
        lines.append(f"  {k}: {v}")
    lines.append(f"cds file: {args.cds}")
    for k, v in cds_stats.items():
        lines.append(f"  {k}: {v}")
    if protein_stats["n_records"] and cds_stats["n_records"]:
        diff = protein_stats["n_records"] - cds_stats["n_records"]
        lines.append(f"protein/CDS record count difference: {diff}")
    lines.append("")
    if errors:
        lines.append("FAILED:")
        lines.extend(f"  - {e}" for e in errors)
    else:
        lines.append("PASSED")

    report = "\n".join(lines)
    print(report)

    if errors:
        sys.exit(f"QC failed for genome {args.genome}: {'; '.join(errors)}")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    main()
