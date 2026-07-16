#!/usr/bin/env python
"""Stage 6 (Phase 5) -- Histograms with overlaid KDE for every protein
property, split by the primary grouping and by the subgroup grouping.

Distribution OVERLAP is the point of this figure type, as opposed to
boxplots' summary-statistic view: two groups whose medians/boxplots look
different can still turn out to overlap substantially once you see the
full distributions -- exactly the comparison a lineage-vs-lifestyle
confound shows up in.

Genome-agnostic by construction -- see plot_boxplots.py's docstring; same
design here (property list, group values/order/colors all derived from the
data + config.yaml, nothing hardcoded).
"""

import argparse

import pandas as pd
from visuals_utils import group_order_and_colors, hist_kde_by_group, property_columns, setup_style


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--primary-order", nargs="*", default=[])
    parser.add_argument("--subgroup-order", nargs="*", default=[])
    parser.add_argument("--exclude-properties", nargs="*", default=[])
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()
    df = pd.read_csv(args.master_protein_table)
    label_columns = [args.primary_grouping, args.subgroup_column]
    props = property_columns(df, label_columns=label_columns, exclude=args.exclude_properties)

    primary_order, primary_colors = group_order_and_colors(df, args.primary_grouping, args.primary_order, "primary")
    subgroup_order, subgroup_colors = group_order_and_colors(df, args.subgroup_column, args.subgroup_order, "subgroup")

    for prop in props:
        hist_kde_by_group(
            df, prop, args.primary_grouping, primary_order, primary_colors,
            f"{prop} distribution by {args.primary_grouping}", args.output_dir,
            f"dist_{prop}_by_{args.primary_grouping}",
        )
        hist_kde_by_group(
            df, prop, args.subgroup_column, subgroup_order, subgroup_colors,
            f"{prop} distribution by {args.subgroup_column}", args.output_dir,
            f"dist_{prop}_by_{args.subgroup_column}",
        )

    print(f"wrote 2 figures x {len(props)} properties -> {args.output_dir}")


if __name__ == "__main__":
    main()
