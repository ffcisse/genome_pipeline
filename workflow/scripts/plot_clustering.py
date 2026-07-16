#!/usr/bin/env python
"""Stage 6 (Phase 5) -- Two clustering views of the protein property
matrix:
  1. Hierarchical clustering heatmap of genomes x properties (z-scored per
     property across genomes), with a dendrogram over genomes -- genomes
     with similar overall property profiles end up adjacent/merged early,
     regardless of what their primary/subgroup labels say.
  2. Property-property correlation heatmap (Spearman), also hierarchically
     clustered -- reveals which properties move together (redundant/
     capturing the same underlying axis) versus independent.

Genome-agnostic by construction -- see plot_boxplots.py's docstring; row
color annotations use the same programmatic per-value colors as every
other figure, for whatever grouping values are actually present.
"""

import argparse

import pandas as pd
import seaborn as sns
from visuals_utils import group_order_and_colors, property_columns, save_figure, setup_style


def hierarchical_heatmap(species_df, props, primary_grouping, subgroup_column, primary_colors, subgroup_colors, output_dir, name):
    indexed = species_df.set_index("genome")
    matrix = indexed[[f"{p}_median" for p in props]].copy()
    matrix.columns = props
    z = (matrix - matrix.mean()) / matrix.std()

    row_colors = pd.DataFrame(
        {
            subgroup_column: indexed[subgroup_column].map(subgroup_colors),
            primary_grouping: indexed[primary_grouping].map(primary_colors),
        }
    )

    g = sns.clustermap(
        z,
        row_colors=row_colors,
        cmap="vlag",
        center=0,
        figsize=(max(10, len(props) * 0.35), max(6, len(matrix) * 0.5 + 3)),
        dendrogram_ratio=(0.15, 0.12),
        cbar_kws={"label": "z-score"},
    )
    g.ax_heatmap.set_xlabel("property")
    g.ax_heatmap.set_ylabel("genome")
    g.fig.suptitle("Genome x property hierarchical clustering (z-scored medians)", y=1.02)
    save_figure(g.fig, output_dir, name)


def correlation_heatmap(df, props, output_dir, name):
    corr = df[props].corr(method="spearman")
    side = max(8, len(props) * 0.4)
    g = sns.clustermap(
        corr,
        cmap="vlag",
        center=0,
        figsize=(side, side),
        cbar_kws={"label": "Spearman r"},
    )
    g.fig.suptitle("Property-property correlation (Spearman, clustered)", y=1.02)
    save_figure(g.fig, output_dir, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--species-summary", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--primary-order", nargs="*", default=[])
    parser.add_argument("--subgroup-order", nargs="*", default=[])
    parser.add_argument("--exclude-properties", nargs="*", default=[])
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()
    protein_df = pd.read_csv(args.master_protein_table)
    species_df = pd.read_csv(args.species_summary)
    label_columns = [args.primary_grouping, args.subgroup_column]
    props = property_columns(protein_df, label_columns=label_columns, exclude=args.exclude_properties)
    props = [p for p in props if f"{p}_median" in species_df.columns]

    _, primary_colors = group_order_and_colors(protein_df, args.primary_grouping, args.primary_order, "primary")
    _, subgroup_colors = group_order_and_colors(protein_df, args.subgroup_column, args.subgroup_order, "subgroup")

    hierarchical_heatmap(
        species_df, props, args.primary_grouping, args.subgroup_column, primary_colors, subgroup_colors,
        args.output_dir, "clustermap_species_properties",
    )
    correlation_heatmap(protein_df, props, args.output_dir, "corr_heatmap_properties")

    print(f"wrote clustering figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
