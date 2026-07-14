#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Group-comparison effect sizes for every numeric
property, across BOTH master tables (protein-side + CDS-side), for a
CONFIGURABLE grouping variable (--grouping lifestyle | lineage -- see
rules/summaries.smk for how this becomes a Snakemake wildcard).

With n~61,000 proteins/CDS, Mann-Whitney p-values collapse to ~0 for nearly
every property regardless of whether the difference is biologically
meaningful -- they're reported for reference, but the columns that actually
answer "how big is this difference" are:
  cles          -- common-language effect size = U / (n_a * n_b): the
                   probability a random member of group_a exceeds a random
                   member of group_b (0.5 = no difference).
  rank_biserial -- 2*cles - 1, rescaled to -1..+1 (0 = no difference).

lifestyle has exactly 2 groups -> one comparison per property. lineage has
3 -> pairwise, in the specific order (Galdieria vs Cyanidiales, Galdieria vs
Mesophile_lineage, Cyanidiales vs Mesophile_lineage) that matches how this
analysis is actually framed (Galdieria as the group of interest), not
alphabetical order -- see GROUP_ORDER. This is exactly why grouping is
configurable rather than hardcoded to lifestyle: the earlier manual
analysis's central finding is that several properties (carbon_oxidation_state,
cysteine_fraction, disorder_fraction) separate by Galdieria lineage, not by
lifestyle -- the Cyanidiales extremophiles actually overlap the mesophiles.
"""

import argparse
import itertools
import os

import pandas as pd
from summaries_utils import compare_groups, numeric_property_columns

# Canonical group order per grouping variable -- determines which group is
# "A" (vs "B", i.e. the direction cles/rank_biserial are signed) in each
# pairwise comparison.
GROUP_ORDER = {
    "lifestyle": ["extremophile", "mesophile"],
    "lineage": ["Galdieria", "Cyanidiales", "Mesophile_lineage"],
}


def effect_sizes_for_table(df: pd.DataFrame, table_name: str, group_col: str) -> list[dict]:
    observed = set(df[group_col].dropna().unique())
    expected = GROUP_ORDER[group_col]
    unexpected = observed - set(expected)
    if unexpected:
        raise AssertionError(
            f"Unexpected {group_col} value(s) {sorted(unexpected)} in {table_name} table -- "
            f"update GROUP_ORDER in effect_sizes.py"
        )
    groups_present = [g for g in expected if g in observed]

    rows = []
    for prop in numeric_property_columns(df):
        for group_a, group_b in itertools.combinations(groups_present, 2):
            result = compare_groups(
                df.loc[df[group_col] == group_a, prop],
                df.loc[df[group_col] == group_b, prop],
            )
            result.update(property=prop, table=table_name, group_a=group_a, group_b=group_b)
            rows.append(result)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grouping", required=True, choices=sorted(GROUP_ORDER))
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)

    rows = effect_sizes_for_table(protein_df, "protein", args.grouping)
    rows += effect_sizes_for_table(cds_df, "cds", args.grouping)

    result = pd.DataFrame(rows)
    result = result.reindex(
        result["rank_biserial"].abs().sort_values(ascending=False, kind="stable").index
    ).reset_index(drop=True)
    result = result[
        ["property", "table", "group_a", "group_b", "n_a", "n_b", "median_a", "median_b", "p_value", "cles", "rank_biserial"]
    ]

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"wrote {len(result)} rows -> {args.output}")


if __name__ == "__main__":
    main()
