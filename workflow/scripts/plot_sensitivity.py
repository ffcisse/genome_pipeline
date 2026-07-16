#!/usr/bin/env python
"""Stage 6 (Phase 5) -- The leave-one-out sensitivity "money figure": a
heatmap of shrinkage (property x excluded_subgroup) for the top properties
by |rank_biserial_full|, showing at a glance which subgroup(s) drive each
property's apparent primary-grouping effect.

Positive shrinkage (see Phase 4's sensitivity_drop_subgroup.py docstring)
means dropping that subgroup shrank the apparent effect -- i.e. that
subgroup was driving it. Negative means dropping it strengthened the
effect (that subgroup had been diluting it, usually because it overlaps
the other primary-grouping side).

Reads Phase 4's already-computed sensitivity_leave_one_out.csv directly --
no statistics recomputed here. Genome-agnostic: which subgroups/properties
appear, and how many of each, come entirely from that CSV; nothing
hardcoded.

Excludes "codon_"-prefixed properties from the top-N ranking, same
DEFAULT_EXCLUDE_PREFIXES rule as plot_effect_sizes.py/
visuals_utils.property_columns -- see plot_effect_sizes.py's docstring for
why (64 raw codon-count columns are strongly GC/lineage-correlated and
would otherwise dominate the ranking purely by count).
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from visuals_utils import DEFAULT_EXCLUDE_PREFIXES, save_figure, setup_style


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sensitivity-leave-one-out", required=True)
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    setup_style()
    df = pd.read_csv(args.sensitivity_leave_one_out)
    df = df[~df["property"].str.startswith(DEFAULT_EXCLUDE_PREFIXES)]

    ranked_props = (
        df[["property", "table", "rank_biserial_full"]]
        .drop_duplicates()
        .assign(abs_full=lambda d: d["rank_biserial_full"].abs())
        .sort_values("abs_full", ascending=False)
        .head(args.top_n)
    )
    top = df.merge(ranked_props[["property", "table"]], on=["property", "table"])
    top = top.assign(label=top["property"] + " (" + top["table"] + ")")

    pivot = top.pivot_table(index="label", columns="excluded_subgroup", values="shrinkage")
    ordered_labels = (ranked_props["property"] + " (" + ranked_props["table"] + ")").tolist()
    pivot = pivot.reindex(ordered_labels)

    vmax = pivot.abs().to_numpy()
    vmax = vmax[~pd.isna(vmax)].max() if vmax.size and not pd.isna(vmax).all() else 1.0

    fig, ax = plt.subplots(figsize=(max(6, 1.4 * pivot.shape[1]), max(5, 0.35 * len(pivot))))
    sns.heatmap(
        pivot,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        ax=ax,
        cbar_kws={"label": "shrinkage"},
    )
    ax.set_xlabel("excluded subgroup")
    ax.set_ylabel("property")
    ax.set_title(f"Leave-one-out sensitivity: shrinkage for top {len(pivot)} properties")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    save_figure(fig, args.output_dir, "sensitivity_leave_one_out_heatmap")

    print(f"wrote sensitivity heatmap -> {args.output_dir}")


if __name__ == "__main__":
    main()
