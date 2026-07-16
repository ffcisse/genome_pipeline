#!/usr/bin/env python
"""Stage 6 (Phase 5) -- CDS-side counterpart to plot_boxplots.py/
plot_distributions.py: boxplots, violin plots, per-species boxplots, and
histogram+KDE distributions for every CDS-level scalar property (length,
ENC, GC, GC3s, canonical start/stop -- whatever numeric columns
cds_properties.py actually produced; the 64 raw per-codon count columns are
excluded by visuals_utils.property_columns' default "codon_" prefix rule,
since they're a matrix better shown by plot_clustering.py's correlation
heatmap/PCA, not 64 near-identical individual boxplots).

Genome-agnostic by construction -- see plot_boxplots.py's docstring; same
design here, just reading master_cds_table.csv instead of
master_protein_table.csv.
"""

import argparse

import pandas as pd
from visuals_utils import (
    boxplot_by_group,
    boxplot_by_species,
    group_order_and_colors,
    hist_kde_by_group,
    property_columns,
    setup_style,
    species_order,
    violin_by_group,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--primary-order", nargs="*", default=[])
    parser.add_argument("--subgroup-order", nargs="*", default=[])
    parser.add_argument("--exclude-properties", nargs="*", default=[])
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()
    df = pd.read_csv(args.master_cds_table)
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

    print(f"wrote 6 figures x {len(props)} CDS properties -> {args.output_dir}")


if __name__ == "__main__":
    main()
