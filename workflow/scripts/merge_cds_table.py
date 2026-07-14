#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Merge every genome's cds_properties.csv +
codon_usage.csv into one master table, with group labels joined in from
config/genomes.tsv. --group-columns says which genome_table columns to
attach (config.yaml's sensitivity.primary_grouping/subgroup_column) -- not
hardcoded here. Same row-count discipline as merge_protein_table.py:
validated against the concatenated input's row count (not a hardcoded
dataset-specific number), both via pandas' `validate="one_to_one"` merge
check and an explicit count comparison.
"""

import argparse
import os

import pandas as pd
from summaries_utils import assert_row_count, load_genome_labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cds-tables", nargs="+", required=True)
    parser.add_argument("--codon-usage-tables", nargs="+", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--group-columns", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    cds_df = pd.concat([pd.read_csv(f) for f in args.cds_tables], ignore_index=True)
    codon_df = pd.concat([pd.read_csv(f) for f in args.codon_usage_tables], ignore_index=True)
    n_input = len(cds_df)
    assert_row_count(codon_df, n_input, "merge_cds_table (cds vs codon_usage input row count)")

    merged = cds_df.merge(codon_df, on=["genome", "cds_id"], how="inner", validate="one_to_one")
    assert_row_count(merged, n_input, "merge_cds_table (cds+codon_usage join)")

    labels = load_genome_labels(args.genomes_tsv, args.group_columns)
    merged = merged.merge(labels, on="genome", how="left", validate="many_to_one")
    assert_row_count(merged, n_input, "merge_cds_table (group-label join)")
    missing = sorted(merged.loc[merged[args.group_columns].isna().any(axis=1), "genome"].unique())
    if missing:
        raise AssertionError(f"No {args.group_columns} label in {args.genomes_tsv} for genomes: {missing}")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"wrote {len(merged)} CDS x {merged.shape[1]} cols -> {args.output}")


if __name__ == "__main__":
    main()
