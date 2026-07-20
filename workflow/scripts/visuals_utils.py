"""Shared helpers for the Phase 5 visuals scripts (plot_*.py) -- mirrors
summaries_utils.py's role for Phase 4: one place for logic every plotting
script needs, instead of copy-pasted per script.

Same discipline as summaries_utils.py, and just as load-bearing here:
nothing in this module (or any plot_*.py script) hardcodes a grouping
column NAME or VALUE, a genome ID, a group count, or a color tied to a
specific group value. Every caller passes grouping columns in explicitly
(config.yaml's sensitivity.primary_grouping/subgroup_column, the same
config Phase 4 uses), and colors/order are derived from whatever group
values are actually present in the data at runtime.
"""

import os

import matplotlib

matplotlib.use("Agg")  # headless -- no display on a login/compute node
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from summaries_utils import DEFAULT_EXCLUDE_PREFIXES, property_columns, resolve_group_order  # noqa: F401

# Two different qualitative palette FAMILIES, keyed by grouping ROLE (which
# config.yaml key the column came from), not by column name or value -- so
# a primary-grouping figure and a subgroup-grouping figure never end up
# reusing the same color for two unrelated things if shown side by side,
# without ever hardcoding e.g. "lifestyle" or "extremophile" to a color.
QUALITATIVE_PALETTES = {
    "primary": "tab10",
    "subgroup": "Dark2",
}


def setup_style() -> None:
    """One consistent, presentation-quality look for every figure this
    phase produces. Call once at the top of every plotting script."""
    plt.rcParams.update(
        {
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.titlesize": 16,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig, output_dir: str, name: str) -> None:
    """Write fig as both <name>.png (300dpi) and <name>.pdf into
    output_dir (created if needed) -- one call, both formats, so no script
    can forget one format for one figure."""
    os.makedirs(output_dir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(output_dir, f"{name}.{ext}"))
    plt.close(fig)


def group_order_and_colors(df: pd.DataFrame, group_col: str, configured_order, role: str):
    """Display order (resolve_group_order -- config override else
    alphabetical, same as Phase 4's effect_sizes/sensitivity scripts) and a
    programmatic color per value actually observed in df[group_col]. Using
    the same order Phase 4 used for this column keeps a figure's left-to-
    right/color order consistent with that column's effect_sizes_<col>.csv
    (e.g. group_a always plotted first)."""
    order = resolve_group_order(df[group_col].dropna().unique(), configured_order)
    cmap = plt.get_cmap(QUALITATIVE_PALETTES.get(role, "tab10"))
    colors = {v: cmap(i % cmap.N) for i, v in enumerate(order)}
    return order, colors


def species_order(genome_labels: pd.DataFrame, subgroup_col: str, configured_order) -> list[str]:
    """Genomes ordered by their subgroup's display order, then
    alphabetically by genome within a subgroup -- "ordered sensibly"
    without hardcoding any genome ID or subgroup value: it's purely a
    function of resolve_group_order's result and the genome_id/subgroup
    columns actually present."""
    subgroup_order = resolve_group_order(genome_labels[subgroup_col].dropna().unique(), configured_order)
    rank = {v: i for i, v in enumerate(subgroup_order)}
    ordered = genome_labels.copy()
    ordered["_rank"] = ordered[subgroup_col].map(rank)
    ordered = ordered.sort_values(["_rank", "genome"])
    return ordered["genome"].tolist()


# ---------------------------------------------------------------------------
# Reusable plot primitives -- shared by plot_boxplots.py, plot_distributions.py,
# and plot_cds_distributions.py (protein-side and CDS-side properties are
# plotted identically; only which master table/columns feed them differs).
# ---------------------------------------------------------------------------


def _boxplot_axes(data, order, colors, ylabel, xlabel, title):
    fig, ax = plt.subplots(figsize=(max(4, len(order) * 0.9), 5))
    bp = ax.boxplot(data, patch_artist=True, showfliers=False)
    for patch, group in zip(bp["boxes"], order):
        patch.set_facecolor(colors[group])
        patch.set_alpha(0.7)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    return fig, ax, bp


def boxplot_by_group(df, prop, group_col, order, colors, title, output_dir, name):
    values = df[prop].astype(float)
    data = [values[df[group_col] == g].dropna() for g in order]
    fig, ax, _ = _boxplot_axes(data, order, colors, prop, group_col, title)
    save_figure(fig, output_dir, name)


def violin_by_group(df, prop, group_col, order, colors, title, output_dir, name):
    values = df[prop].astype(float)
    data = [values[df[group_col] == g].dropna() for g in order]
    fig, ax = plt.subplots(figsize=(max(4, len(order) * 0.9), 5))
    parts = ax.violinplot(data, showmedians=True)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[order[i]])
        body.set_alpha(0.7)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel(prop)
    ax.set_xlabel(group_col)
    ax.set_title(title)
    save_figure(fig, output_dir, name)


def boxplot_by_species(df, prop, genome_order, genome_colors, subgroup_col, output_dir, name):
    values = df[prop].astype(float)
    data = [values[df["genome"] == g].dropna() for g in genome_order]
    fig, ax, bp = _boxplot_axes(
        data, genome_order, genome_colors, prop, "genome", f"{prop} by genome (colored by {subgroup_col})"
    )
    for i, median_line in enumerate(bp["medians"], start=1):
        med = median_line.get_ydata()[0]
        ax.annotate(f"{med:.3g}", (i, med), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save_figure(fig, output_dir, name)


def hist_kde_by_group(df, prop, group_col, order, colors, title, output_dir, name):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    plot_df = df[[prop, group_col]].assign(**{prop: df[prop].astype(float)})
    sns.histplot(
        data=plot_df,
        x=prop,
        hue=group_col,
        hue_order=order,
        palette=colors,
        stat="density",
        common_norm=False,
        kde=True,
        alpha=0.35,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel(prop)
    ax.set_ylabel("density")
    save_figure(fig, output_dir, name)
