#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Group-comparison effect sizes for every numeric
property, across BOTH master tables (protein-side + CDS-side), for a
CONFIGURABLE grouping variable (--grouping, e.g. lifestyle | lineage -- see
rules/summaries.smk for how this becomes a Snakemake wildcard, constrained
to config.yaml's sensitivity.primary_grouping/subgroup_column).

With n~61,000 proteins/CDS, Mann-Whitney p-values collapse to ~0 for nearly
every property regardless of whether the difference is biologically
meaningful -- they're reported for reference, but the columns that actually
answer "how big is this difference" are:
  cles          -- common-language effect size = U / (n_a * n_b): the
                   probability a random member of group_a exceeds a random
                   member of group_b (0.5 = no difference).
  rank_biserial -- 2*cles - 1, rescaled to -1..+1 (0 = no difference).

A 2-value grouping produces one comparison per property; an N-value
grouping produces all pairwise comparisons. Which value is "A" (vs "B",
i.e. the sign convention for cles/rank_biserial) and the order pairwise
comparisons are generated in both come from --group-order if given (sourced
from config.yaml's optional group_value_order block), else alphabetical --
see summaries_utils.resolve_group_order. Nothing about specific group names
or values (e.g. "Galdieria", "extremophile") is hardcoded here: this used to
have a GROUP_ORDER dict pinned to this deployment's exact groups, which
meant running on a different genome set required editing this file, not
just config -- see git history for that version if you're curious what
changed and why.
"""

import argparse
import itertools
import os

import pandas as pd
from summaries_utils import compare_groups, genome_table_columns, numeric_property_columns, resolve_group_order


def effect_sizes_for_table(df: pd.DataFrame, table_name: str, group_col: str, label_columns, group_order) -> list[dict]:
    observed = df[group_col].dropna().unique()
    groups_present = resolve_group_order(observed, group_order)

    rows = []
    for prop in numeric_property_columns(df, label_columns=label_columns):
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
    parser.add_argument("--grouping", required=True, help="Column in the master tables, e.g. lifestyle")
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument(
        "--group-order",
        nargs="*",
        default=[],
        help="Preferred order of --grouping's distinct values (which is 'A' in cles/rank_biserial); "
        "omit for alphabetical order",
    )
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    label_columns = genome_table_columns(args.genomes_tsv)
    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)

    rows = effect_sizes_for_table(protein_df, "protein", args.grouping, label_columns, args.group_order)
    rows += effect_sizes_for_table(cds_df, "cds", args.grouping, label_columns, args.group_order)

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
