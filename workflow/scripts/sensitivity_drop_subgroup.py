#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Leave-one-subgroup-out sensitivity analysis: a
standard technique for nested/confounded group designs, testing whether an
apparent effect on a PRIMARY grouping variable (e.g. lifestyle) is actually
being driven by one value of a finer-grained SUBGROUP column (e.g. one
lineage) -- by excluding that subgroup's genomes entirely and recomputing
the primary-grouping effect size. Fully configurable
(--primary-grouping/--subgroup-column/--exclude-subgroup); nothing here is
specific to any one dataset's group names -- this deployment's defaults
(lifestyle/lineage) live in config.yaml's `sensitivity:` block and
workflow/Snakefile's SUBGROUPS, not in this script.

For every numeric property (both master tables), computes the
primary-grouping effect size (rank_biserial) twice:
  rank_biserial_full     -- every genome
  rank_biserial_excluded -- every genome EXCEPT those where
                            subgroup_column == exclude_subgroup (removed
                            entirely from the whole dataset -- not merely
                            relabeled -- before splitting by primary_grouping
                            and comparing again)
  shrinkage = |rank_biserial_full| - |rank_biserial_excluded|

Positive shrinkage means excluding that subgroup shrank the apparent
primary-grouping effect -- i.e. the difference was partly (or wholly) an
artifact of that one subgroup rather than a genuine primary_grouping
effect. A shrinkage near zero means the property separates by
primary_grouping regardless of that subgroup.

primary_grouping must have exactly 2 distinct values present (both before
and after exclusion) -- rank_biserial is inherently a 2-group effect size;
for >2 groups, use effect_sizes.py's pairwise comparisons instead. Group
order (which value is "A") comes from --group-order if given (same
optional config.yaml override effect_sizes.py uses -- see
summaries_utils.resolve_group_order), else alphabetical. shrinkage itself is
computed from absolute values so this ordering choice doesn't change it
either way, but honoring the same configured order keeps the signed
rank_biserial_full/rank_biserial_excluded columns here consistent with the
corresponding row in effect_sizes_<primary_grouping>.csv.

Edge case (flagged, not crashed): if removing exclude_subgroup leaves either
primary_grouping value with fewer than 2 genomes backing it, a Mann-Whitney
comparison would be resting on pseudo-replication from a single genome (or
none) -- rank_biserial_excluded/shrinkage are set to NaN for every property
in that table, with a warning printed, rather than silently reporting a
statistic that looks normal but isn't trustworthy.
"""

import argparse
import os

import pandas as pd
from summaries_utils import compare_groups, genome_table_columns, numeric_property_columns, resolve_group_order


def sensitivity_for_table(
    df: pd.DataFrame,
    table_name: str,
    primary_grouping: str,
    subgroup_column: str,
    exclude_subgroup: str,
    label_columns,
    group_order,
) -> list[dict]:
    observed = df[primary_grouping].dropna().unique()
    if len(set(observed)) != 2:
        raise AssertionError(
            f"--primary-grouping {primary_grouping!r} has {len(set(observed))} distinct value(s) "
            f"in the {table_name} table ({sorted(set(observed))}), expected exactly 2"
        )
    group_a, group_b = resolve_group_order(observed, group_order)

    excluded_df = df[df[subgroup_column] != exclude_subgroup]
    n_genomes_a = excluded_df.loc[excluded_df[primary_grouping] == group_a, "genome"].nunique()
    n_genomes_b = excluded_df.loc[excluded_df[primary_grouping] == group_b, "genome"].nunique()
    safe_to_compare = n_genomes_a >= 2 and n_genomes_b >= 2
    if not safe_to_compare:
        print(
            f"WARNING: excluding {subgroup_column}={exclude_subgroup!r} from the {table_name} table leaves only "
            f"{n_genomes_a} genome(s) for {primary_grouping}={group_a!r} and {n_genomes_b} for "
            f"{primary_grouping}={group_b!r} -- too few for a trustworthy comparison (need >=2 each). "
            f"rank_biserial_excluded/shrinkage will be NaN for every {table_name} property."
        )

    rows = []
    for prop in numeric_property_columns(df, label_columns=label_columns):
        full = compare_groups(df.loc[df[primary_grouping] == group_a, prop], df.loc[df[primary_grouping] == group_b, prop])
        if safe_to_compare:
            excl = compare_groups(
                excluded_df.loc[excluded_df[primary_grouping] == group_a, prop],
                excluded_df.loc[excluded_df[primary_grouping] == group_b, prop],
            )
            rank_biserial_excluded = excl["rank_biserial"]
            shrinkage = abs(full["rank_biserial"]) - abs(rank_biserial_excluded)
        else:
            rank_biserial_excluded = float("nan")
            shrinkage = float("nan")
        rows.append(
            {
                "property": prop,
                "table": table_name,
                "rank_biserial_full": full["rank_biserial"],
                "rank_biserial_excluded": rank_biserial_excluded,
                "shrinkage": shrinkage,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-grouping", required=True, help="Column in the master tables, e.g. lifestyle")
    parser.add_argument("--subgroup-column", required=True, help="Column in the master tables, e.g. lineage")
    parser.add_argument("--exclude-subgroup", required=True, help="Value of --subgroup-column to drop entirely")
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument(
        "--group-order",
        nargs="*",
        default=[],
        help="Preferred order of --primary-grouping's distinct values (which is 'A'); omit for alphabetical order",
    )
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    label_columns = genome_table_columns(args.genomes_tsv)
    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)

    rows = sensitivity_for_table(
        protein_df, "protein", args.primary_grouping, args.subgroup_column, args.exclude_subgroup,
        label_columns, args.group_order,
    ) + sensitivity_for_table(
        cds_df, "cds", args.primary_grouping, args.subgroup_column, args.exclude_subgroup,
        label_columns, args.group_order,
    )

    result = pd.DataFrame(rows)
    result = result.reindex(
        result["shrinkage"].abs().sort_values(ascending=False, kind="stable").index
    ).reset_index(drop=True)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"wrote {len(result)} rows -> {args.output}")


if __name__ == "__main__":
    main()
