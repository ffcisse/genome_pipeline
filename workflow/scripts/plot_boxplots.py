#!/usr/bin/env python
"""Stage 6 (Phase 5) -- Boxplots and violin plots of every protein
property: split by the primary grouping, by the subgroup grouping, and by
individual species (species colored by their subgroup value, with a median
label on each box).

Genome-agnostic by construction: the property list, group values/order/
colors, and species order are all derived from the data + config.yaml (see
visuals_utils.py). Nothing here names a specific property, genome, or group
value in control flow -- only plot titles/labels are human-readable
strings built FROM the data at runtime.
"""

import argparse

import pandas as pd
from visuals_utils import (
    boxplot_by_group,
    boxplot_by_species,
    group_order_and_colors,
    property_columns,
    setup_style,
    species_order,
    violin_by_group,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--genomes-tsv", required=True)
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

    genome_labels = pd.read_csv(args.genomes_tsv, sep="\t").rename(columns={"genome_id": "genome"})
    genome_order = species_order(genome_labels, args.subgroup_column, args.subgroup_order)
    genome_to_subgroup = dict(zip(genome_labels["genome"], genome_labels[args.subgroup_column]))
    genome_colors = {g: subgroup_colors[genome_to_subgroup[g]] for g in genome_order}

    for prop in props:
        boxplot_by_group(
            df, prop, args.primary_grouping, primary_order, primary_colors,
            f"{prop} by {args.primary_grouping}", args.output_dir,
            f"boxplot_{prop}_by_{args.primary_grouping}",
        )
        violin_by_group(
            df, prop, args.primary_grouping, primary_order, primary_colors,
            f"{prop} by {args.primary_grouping} (violin)", args.output_dir,
            f"violin_{prop}_by_{args.primary_grouping}",
        )
        boxplot_by_group(
            df, prop, args.subgroup_column, subgroup_order, subgroup_colors,
            f"{prop} by {args.subgroup_column}", args.output_dir,
            f"boxplot_{prop}_by_{args.subgroup_column}",
        )
        violin_by_group(
            df, prop, args.subgroup_column, subgroup_order, subgroup_colors,
            f"{prop} by {args.subgroup_column} (violin)", args.output_dir,
            f"violin_{prop}_by_{args.subgroup_column}",
        )
        boxplot_by_species(
            df, prop, genome_order, genome_colors, args.subgroup_column, args.output_dir,
            f"boxplot_{prop}_by_species",
        )

    print(f"wrote 5 figures x {len(props)} properties -> {args.output_dir}")


if __name__ == "__main__":
    main()
