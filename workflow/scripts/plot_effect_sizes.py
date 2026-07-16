#!/usr/bin/env python
"""Stage 6 (Phase 5) -- Forest/bar plots of rank-biserial effect sizes per
property, sorted by magnitude: one figure for the primary grouping, and
one per pairwise comparison for the subgroup grouping (a 3-value subgroup
column produces 3 figures, a 5-value one would produce 10 -- however many
pairs Phase 4's effect_sizes_<subgroup>.csv actually contains).

Reads Phase 4's already-computed effect_sizes_<grouping>.csv files directly
-- no statistics are recomputed here, only visualized. Genome-agnostic:
which pairs exist, and how many, come entirely from the CSV's own
group_a/group_b values; nothing here assumes a specific pair count or
group names.

Excludes cds_properties.py's "codon_"-prefixed rows (64 raw per-codon
counts) from the top-N ranking, same DEFAULT_EXCLUDE_PREFIXES rule
visuals_utils.property_columns applies elsewhere: those 64 columns are
strongly GC-content/lineage-correlated and numerically dominate a
combined ranking purely by sheer count, crowding out the individual
protein/CDS-scalar properties this figure exists to surface. The
underlying effect_sizes_<grouping>.csv still has every property, codon
columns included -- this filter only affects which rows this one figure
highlights.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from visuals_utils import DEFAULT_EXCLUDE_PREFIXES, save_figure, setup_style


def forest_plot(df_pair: pd.DataFrame, title: str, output_dir: str, name: str, top_n: int) -> None:
    df_pair = df_pair[~df_pair["property"].str.startswith(DEFAULT_EXCLUDE_PREFIXES)]
    ranked = df_pair.reindex(df_pair["rank_biserial"].abs().sort_values(ascending=False).index).head(top_n)
    ranked = ranked.iloc[::-1]  # smallest-of-the-top at bottom, largest at top for barh
    labels = ranked["property"] + " (" + ranked["table"] + ")"
    colors = ["#d62728" if v < 0 else "#1f77b4" for v in ranked["rank_biserial"]]

    fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(ranked))))
    ax.barh(labels, ranked["rank_biserial"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("rank-biserial correlation")
    ax.set_title(title)
    save_figure(fig, output_dir, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effect-sizes-primary", required=True)
    parser.add_argument("--effect-sizes-subgroup", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()

    primary_df = pd.read_csv(args.effect_sizes_primary)
    forest_plot(
        primary_df,
        f"Top {args.top_n} effect sizes by {args.primary_grouping}",
        args.output_dir,
        f"effect_sizes_{args.primary_grouping}",
        args.top_n,
    )

    subgroup_df = pd.read_csv(args.effect_sizes_subgroup)
    pairs = subgroup_df[["group_a", "group_b"]].drop_duplicates().itertuples(index=False)
    for group_a, group_b in pairs:
        pair_df = subgroup_df[(subgroup_df["group_a"] == group_a) & (subgroup_df["group_b"] == group_b)]
        forest_plot(
            pair_df,
            f"Top {args.top_n} effect sizes: {group_a} vs {group_b} ({args.subgroup_column})",
            args.output_dir,
            f"effect_sizes_{args.subgroup_column}_{group_a}_vs_{group_b}",
            args.top_n,
        )

    print(f"wrote effect-size figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
