#!/usr/bin/env python
"""Stage 5 (Phase 4) -- One row per genome: median/mean/std of every
numeric property from both master tables (protein-side: composition,
physicochemical, disorder; CDS-side: ENC/GC/GC3s, codon usage), plus gene
counts and the configured group labels (config.yaml's
sensitivity.primary_grouping/subgroup_column, e.g. lifestyle/lineage).

Column-name-agnostic over properties (see summaries_utils.numeric_property_
columns) -- if a future phase adds a new property column to either master
table, it shows up here automatically as <property>_median/_mean/_std,
no edit needed. --group-columns excludes the master tables' own group-label
columns from that treatment (they're strings, not properties).
"""

import argparse
import os

import pandas as pd
from summaries_utils import load_genome_labels, numeric_property_columns


def per_genome_stats(df: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    props = numeric_property_columns(df, label_columns=group_columns)
    stats = df.groupby("genome")[props].agg(["median", "mean", "std"])
    stats.columns = [f"{prop}_{stat}" for prop, stat in stats.columns]
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--group-columns", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)

    protein_stats = per_genome_stats(protein_df, args.group_columns)
    cds_stats = per_genome_stats(cds_df, args.group_columns)

    n_proteins = protein_df.groupby("genome").size().rename("n_proteins")
    n_cds = cds_df.groupby("genome").size().rename("n_cds")

    summary = pd.concat([n_proteins, n_cds, protein_stats, cds_stats], axis=1).reset_index()

    labels = load_genome_labels(args.genomes_tsv, args.group_columns)
    summary = labels.merge(summary, on="genome", how="left", validate="one_to_one")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    summary.to_csv(args.output, index=False)
    print(f"wrote {len(summary)} genomes x {summary.shape[1]} cols -> {args.output}")


if __name__ == "__main__":
    main()
