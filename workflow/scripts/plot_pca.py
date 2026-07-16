#!/usr/bin/env python
"""Stage 6 (Phase 5) -- PCA of the protein property matrix, at both
per-protein (all 61k-ish rows) and per-species (genome medians) resolution,
each colored by the primary grouping and separately by the subgroup
grouping, plus a loadings plot showing which properties drive PC1/PC2.

Genome-agnostic by construction -- see plot_boxplots.py's docstring;
property list and group values/order/colors come from data + config.yaml,
never hardcoded. The species-median PCA's points are labeled with genome
IDs from genome_table, so an outlier genome (whichever one, for whatever
dataset this runs on) is identifiable directly on the plot rather than
needing a legend cross-reference.

Standardizes (z-score) each property before PCA, since the property matrix
mixes wildly different scales (e.g. pI ~0-14 vs fractions ~0-1) --
without this, PCA would just reflect whichever property happens to have
the largest raw variance, not real covariance structure.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from visuals_utils import group_order_and_colors, property_columns, save_figure, setup_style


def run_pca(matrix: pd.DataFrame):
    scaled = StandardScaler().fit_transform(matrix.to_numpy())
    n_components = min(matrix.shape[0], matrix.shape[1], 10)
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(scaled)
    return pca, scores


def scatter_pca(scores, group_values, order, colors, group_col, pca, title, output_dir, name, alpha, size, annotate=None):
    fig, ax = plt.subplots(figsize=(7, 6))
    for group in order:
        mask = group_values == group
        ax.scatter(
            scores[mask, 0], scores[mask, 1], color=colors[group], label=group, alpha=alpha, s=size, edgecolor="none"
        )
    if annotate is not None:
        for i, txt in enumerate(annotate):
            ax.annotate(txt, (scores[i, 0], scores[i, 1]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% var)")
    ax.set_title(title)
    ax.legend(title=group_col, frameon=False)
    save_figure(fig, output_dir, name)


def loadings_plot(pca, prop_names, output_dir, name, top_n):
    loadings = pd.DataFrame(pca.components_[:2].T, index=prop_names, columns=["PC1", "PC2"])
    for pc in ["PC1", "PC2"]:
        top = loadings[pc].abs().sort_values(ascending=False).head(top_n).index
        vals = loadings.loc[top, pc].sort_values()
        fig, ax = plt.subplots(figsize=(6, max(4, 0.35 * len(vals))))
        bar_colors = ["#d62728" if v < 0 else "#1f77b4" for v in vals]
        ax.barh(vals.index, vals.to_numpy(), color=bar_colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel(f"{pc} loading")
        ax.set_title(f"Top {len(vals)} property loadings on {pc}")
        save_figure(fig, output_dir, f"{name}_{pc}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--species-summary", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--primary-order", nargs="*", default=[])
    parser.add_argument("--subgroup-order", nargs="*", default=[])
    parser.add_argument("--exclude-properties", nargs="*", default=[])
    parser.add_argument("--top-n-loadings", type=int, default=15)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()
    df = pd.read_csv(args.master_protein_table)
    label_columns = [args.primary_grouping, args.subgroup_column]
    props = property_columns(df, label_columns=label_columns, exclude=args.exclude_properties)

    primary_order, primary_colors = group_order_and_colors(df, args.primary_grouping, args.primary_order, "primary")
    subgroup_order, subgroup_colors = group_order_and_colors(df, args.subgroup_column, args.subgroup_order, "subgroup")
    groupings = [
        (args.primary_grouping, primary_order, primary_colors),
        (args.subgroup_column, subgroup_order, subgroup_colors),
    ]

    # --- per-protein PCA ---
    protein_matrix = df[props].dropna()
    pca_p, scores_p = run_pca(protein_matrix)
    idx = protein_matrix.index
    for group_col, order, colors in groupings:
        scatter_pca(
            scores_p, df.loc[idx, group_col].to_numpy(), order, colors, group_col, pca_p,
            f"Per-protein PCA, colored by {group_col}", args.output_dir, f"pca_protein_by_{group_col}",
            alpha=0.35, size=10,
        )
    loadings_plot(pca_p, props, args.output_dir, "pca_protein_loadings", args.top_n_loadings)

    # --- per-species (median) PCA ---
    species_df = pd.read_csv(args.species_summary)
    median_cols = [f"{p}_median" for p in props if f"{p}_median" in species_df.columns]
    species_matrix = species_df[median_cols].dropna()
    species_matrix.columns = [c[: -len("_median")] for c in species_matrix.columns]
    pca_s, scores_s = run_pca(species_matrix)
    sidx = species_matrix.index
    genome_names = species_df.loc[sidx, "genome"].tolist()
    for group_col, order, colors in groupings:
        scatter_pca(
            scores_s, species_df.loc[sidx, group_col].to_numpy(), order, colors, group_col, pca_s,
            f"Per-species (median) PCA, colored by {group_col}", args.output_dir,
            f"pca_species_by_{group_col}", alpha=0.9, size=140, annotate=genome_names,
        )
    loadings_plot(pca_s, species_matrix.columns.tolist(), args.output_dir, "pca_species_loadings", args.top_n_loadings)

    print(f"wrote PCA figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
