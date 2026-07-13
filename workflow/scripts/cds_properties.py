#!/usr/bin/env python
"""Stage 4 (Phase 3) -- CDS/codon-usage properties for one genome, via codonW.

codonW (Peden, bioconda package `codonw`) is the authoritative source for
ENC (Nc), GC, and GC3s here -- these are NOT reimplemented in Python, only
codonW's own output is parsed. Two things codonW does not report are
computed directly from the sequences in Python instead: start codon and
stop codon (per gene, categorical -- not a genome-wide tally, to keep the
same one-row-per-CDS granularity as everything else in this pipeline).

Why a synthetic FASTA header (">seq_0", ">seq_1", ...) instead of the real
cds_id: codonW's machine-readable index output truncates the sequence title
to a fixed 25-character column. Real IDs here (e.g.
"jgi|CyamerSoos_1_1|163|CmerSOOS_G163.1") can collide after truncation once
the genome_id prefix is long enough to push the distinguishing suffix past
character 25 -- this is exactly what forced the lossy "rename everything to
gene_N" workaround in the earlier manual analysis
(~/codonw_test/clean_fasta_for_codonw.py), which threw away the ability to
map codonW's output back to real gene IDs. Short synthetic headers avoid
truncation entirely, and results are re-attached to the real cds_id by
*row position*, not by title string -- confirmed safe because codonW
preserves input order 1:1 in both its index and bulk output (verified: fed
it a deliberately malformed/reordered-looking mix of sequences and the
output order exactly matched input order), including for sequences with a
partial trailing codon or internal N's (codonW warns and silently
truncates/excludes rather than dropping the record, so row count is always
preserved).

Conventions (verified empirically against codonW 1.4.4/bioconda, matching
the values already in ~/codonw_results/all_species_codonw.csv from the
earlier manual run on this same data):
  gc, gc3s  -- proportions in [0, 1] (NOT percentages)
  enc       -- codonW's "Nc", natural scale ~20-61 (NOT normalized to [0,1]).
               "*****" (codonW's sentinel for "too few synonymous codons to
               compute Nc", e.g. very short CDS) becomes NaN here.
  gc3s specifically (not naive "GC3") is GC at *synonymous* 3rd-codon-position
  sites only -- the standard population-genetics metric, and literally the
  name of the codonW flag/column this pipeline's outline asked for.

Codon usage: codonW's `-cu` bulk output, with `-machine`, is a fixed-order
64-number table (raw counts, no header) -- CODON_ORDER below is that fixed
order, reverse-engineered by diffing against codonW's own human-readable
labeled output for the same input (see the rule's docstring). Written as a
separate results/cds_properties/<genome>/codon_usage.csv (genome + cds_id +
64 count columns) rather than folded into cds_properties.csv, same
one-file-per-concern reasoning as Phase 2b's disorder.csv: keeps the
scalar-properties file narrow/readable, and any rule should only ever have
one thing responsible for writing a given file.
"""

import argparse
import os
import subprocess
import tempfile

import pandas as pd

STOP_CODONS = {"TAA", "TAG", "TGA"}

# Fixed output order of codonW's `-cu -machine` bulk table: 4 lines of 16
# codons, grouped by 1st base in T,C,A,G order, each group a 4x4 grid over
# (2nd base x 3rd base), also T,C,A,G. Reverse-engineered by running codonW
# with -cu (no -machine) on the same input, which labels every codon, and
# matching value-for-value against the -machine numbers in the same
# position -- see the rule/script docstrings for the verification.
CODON_ORDER = (
    "TTT TCT TAT TGT TTC TCC TAC TGC TTA TCA TAA TGA TTG TCG TAG TGG "
    "CTT CCT CAT CGT CTC CCC CAC CGC CTA CCA CAA CGA CTG CCG CAG CGG "
    "ATT ACT AAT AGT ATC ACC AAC AGC ATA ACA AAA AGA ATG ACG AAG AGG "
    "GTT GCT GAT GGT GTC GCC GAC GGC GTA GCA GAA GGA GTG GCG GAG GGG"
).split()


def write_codonw_fasta(sequences: pd.Series, path: str) -> None:
    with open(path, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq_{i}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j : j + 60] + "\n")


def run_codonw(fasta_path: str, index_path: str, bulk_path: str) -> None:
    subprocess.run(
        [
            "codonw",
            fasta_path,
            index_path,
            bulk_path,
            "-nomenu",
            "-silent",
            "-nowarn",
            "-machine",
            "-enc",
            "-gc",
            "-gc3s",
            "-cu",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def parse_index_file(path: str, n: int) -> pd.DataFrame:
    """codonW's index file: header line, then one row per gene --
    "seq_<i> Nc GC3s GC", whitespace-separated. Nc is "*****" when codonW
    couldn't compute it (too few synonymous-codon-eligible amino acids)."""
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = [c.strip().lower() for c in df.columns]
    df = df.rename(columns={"nc": "enc"})
    assert len(df) == n, f"codonW index file has {len(df)} rows, expected {n}"
    assert list(df["title"]) == [f"seq_{i}" for i in range(n)], (
        "codonW index file row order doesn't match input order -- positional join is unsafe"
    )
    for col in ("enc", "gc3s", "gc"):
        df[col] = pd.to_numeric(df[col], errors="coerce")  # "*****" -> NaN
    return df[["enc", "gc3s", "gc"]].reset_index(drop=True)


def parse_bulk_file(path: str, n: int) -> pd.DataFrame:
    """codonW's `-cu -machine` bulk file: exactly 4 lines per gene, 16
    whitespace-separated counts per line (any trailing "Codons=N"/"Universal
    Genetic code"/gene-name text on some lines is ignored -- always past the
    16th token)."""
    with open(path) as f:
        lines = [line for line in f if line.strip()]
    assert len(lines) == 4 * n, f"codonW bulk file has {len(lines)} lines, expected {4 * n} (4 per gene)"

    rows = []
    for i in range(n):
        chunk = lines[4 * i : 4 * i + 4]
        counts = []
        for line in chunk:
            counts.extend(int(tok) for tok in line.split()[:16])
        assert len(counts) == 64
        gene_name = chunk[3].split()[16]
        assert gene_name == f"seq_{i}", (
            f"codonW bulk file gene {i} is labeled {gene_name!r}, expected seq_{i} -- positional join is unsafe"
        )
        rows.append(counts)
    return pd.DataFrame(rows, columns=[f"codon_{c}" for c in CODON_ORDER])


def start_stop_codons(seq: str) -> tuple[str, str]:
    trimmed_len = len(seq) - (len(seq) % 3)
    start = seq[0:3] if len(seq) >= 3 else ""
    stop = seq[trimmed_len - 3 : trimmed_len] if trimmed_len >= 3 else ""
    return start, stop


def compute_cds_properties(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    n = len(df)
    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "input.fasta")
        index_path = os.path.join(tmp, "index.out")
        bulk_path = os.path.join(tmp, "bulk.out")
        write_codonw_fasta(df["sequence"], fasta_path)
        run_codonw(fasta_path, index_path, bulk_path)
        codonw_df = parse_index_file(index_path, n)
        codon_usage_df = parse_bulk_file(bulk_path, n)

    starts, stops = zip(*(start_stop_codons(seq) for seq in df["sequence"])) if n else ((), ())

    properties = pd.DataFrame(
        {
            "genome": df["genome"].reset_index(drop=True),
            "cds_id": df["cds_id"].reset_index(drop=True),
            "length": df["length"].reset_index(drop=True),
            "is_multiple_of_3": (df["length"] % 3 == 0).reset_index(drop=True),
            "enc": codonw_df["enc"],
            "gc": codonw_df["gc"],
            "gc3s": codonw_df["gc3s"],
            "start_codon": starts,
            "stop_codon": stops,
            "has_canonical_start": [s == "ATG" for s in starts],
            "has_canonical_stop": [s in STOP_CODONS for s in stops],
        }
    )
    codon_usage = pd.concat(
        [
            df[["genome", "cds_id"]].reset_index(drop=True),
            codon_usage_df,
        ],
        axis=1,
    )
    return properties, codon_usage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--cds-table", required=True)
    parser.add_argument("--properties-output", required=True)
    parser.add_argument("--codon-usage-output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.cds_table)
    properties, codon_usage = compute_cds_properties(df)

    os.makedirs(os.path.dirname(args.properties_output), exist_ok=True)
    os.makedirs(os.path.dirname(args.codon_usage_output), exist_ok=True)
    properties.to_csv(args.properties_output, index=False)
    codon_usage.to_csv(args.codon_usage_output, index=False)
    print(
        f"{args.genome}: wrote {len(properties)} CDS x {properties.shape[1]} cols -> {args.properties_output}, "
        f"{codon_usage.shape[1]} cols -> {args.codon_usage_output}"
    )


if __name__ == "__main__":
    main()
