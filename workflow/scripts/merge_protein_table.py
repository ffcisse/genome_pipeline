#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Merge every genome's protein_properties.csv +
disorder.csv into one master table, with lifestyle/lineage group labels
joined in from config/genomes.tsv (kept separate from the per-genome result
tables by design -- see config/genomes.tsv and the Snakefile's docstring --
so this merge is the one place group labels actually get attached).

Row-count discipline: rather than hardcode the current dataset's 61,349
figure (which would silently stop being checked the moment a genome is
added/removed), this validates that the merge produced exactly as many rows
as went in -- both via pandas' own `validate="one_to_one"` merge check
(catches an unexpected many-to-many join immediately) and an explicit count
comparison against the concatenated input (catches silent row loss from a
mismatched join key that validate= wouldn't necessarily flag). Both fail
loudly (raise), not warn.
"""

import argparse
import os

import pandas as pd
from summaries_utils import assert_row_count, load_genome_labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-tables", nargs="+", required=True)
    parser.add_argument("--disorder-tables", nargs="+", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protein_df = pd.concat([pd.read_csv(f) for f in args.protein_tables], ignore_index=True)
    disorder_df = pd.concat([pd.read_csv(f) for f in args.disorder_tables], ignore_index=True)
    n_input = len(protein_df)
    assert_row_count(disorder_df, n_input, "merge_protein_table (protein vs disorder input row count)")

    merged = protein_df.merge(disorder_df, on=["genome", "protein_id"], how="inner", validate="one_to_one")
    assert_row_count(merged, n_input, "merge_protein_table (protein+disorder join)")

    labels = load_genome_labels(args.genomes_tsv)
    merged = merged.merge(labels, on="genome", how="left", validate="many_to_one")
    assert_row_count(merged, n_input, "merge_protein_table (group-label join)")
    missing = sorted(merged.loc[merged["lifestyle"].isna(), "genome"].unique())
    if missing:
        raise AssertionError(f"No lifestyle/lineage label in {args.genomes_tsv} for genomes: {missing}")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"wrote {len(merged)} proteins x {merged.shape[1]} cols -> {args.output}")


if __name__ == "__main__":
    main()
