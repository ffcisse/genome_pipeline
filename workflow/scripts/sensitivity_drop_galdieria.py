#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Sensitivity analysis: how much of each property's
lifestyle effect size is actually a Galdieria-lineage artifact.

For every numeric property (both master tables), computes the
extremophile-vs-mesophile rank_biserial twice:
  rank_biserial_full         -- all 6 extremophiles (Cyanidiales + Galdieria)
                                 vs the 3 mesophiles (identical to the
                                 lifestyle row in effect_sizes_lifestyle.csv)
  rank_biserial_no_galdieria -- Cyanidiales-only (the 3 non-Galdieria
                                 extremophiles) vs the same 3 mesophiles,
                                 i.e. the 3 Galdieria genomes dropped
                                 entirely rather than merely relabeled

shrinkage = |rank_biserial_full| - |rank_biserial_no_galdieria|. Positive
shrinkage means the lifestyle grouping's apparent effect size was inflated
by Galdieria specifically -- once Galdieria is removed, Cyanidiales alone
doesn't separate from the mesophiles nearly as much, i.e. "lifestyle" was
substantially standing in for "is this genome Galdieria". A shrinkage near
zero means the property genuinely separates by lifestyle independent of
Galdieria.
"""

import argparse
import os

import pandas as pd
from summaries_utils import compare_groups, numeric_property_columns


def sensitivity_for_table(df: pd.DataFrame, table_name: str) -> list[dict]:
    mesophile = df["lifestyle"] == "mesophile"
    extremophile_full = df["lifestyle"] == "extremophile"
    cyanidiales_only = df["lineage"] == "Cyanidiales"

    rows = []
    for prop in numeric_property_columns(df):
        full = compare_groups(df.loc[extremophile_full, prop], df.loc[mesophile, prop])
        no_gal = compare_groups(df.loc[cyanidiales_only, prop], df.loc[mesophile, prop])
        rows.append(
            {
                "property": prop,
                "table": table_name,
                "rank_biserial_full": full["rank_biserial"],
                "rank_biserial_no_galdieria": no_gal["rank_biserial"],
                "shrinkage": abs(full["rank_biserial"]) - abs(no_gal["rank_biserial"]),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)

    rows = sensitivity_for_table(protein_df, "protein") + sensitivity_for_table(cds_df, "cds")

    result = pd.DataFrame(rows)
    result = result.reindex(
        result["shrinkage"].abs().sort_values(ascending=False, kind="stable").index
    ).reset_index(drop=True)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"wrote {len(result)} rows -> {args.output}")


if __name__ == "__main__":
    main()
