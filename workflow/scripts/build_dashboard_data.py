#!/usr/bin/env python
"""Stage 7 (Phase 6a/6b) -- Build the JSON data payload the dashboard HTML
embeds inline.

Two different kinds of numbers go into this payload, and they must never
be confused with each other (see dashboard.smk's docstring and every
section's on-screen caption):

  * EXACT statistics -- property_stats/cds_property_stats (per-genome and
    per-grouping-value min/q1/median/q3/max/mean/std/n), effect_sizes/
    sensitivity pass-throughs, codon usage frequencies, PCA fit on the
    FULL data (loadings, explained variance, per-species scores), and
    property-property correlation matrices -- computed here from the FULL
    master_protein_table/master_cds_table (all ~61k rows for this
    deployment's data), never from the sample below. These drive the
    "exact" box-plot trace and every summary number shown anywhere in the
    dashboard.

  * A DOWNSAMPLED per-protein/per-CDS sample (samples.*/cds_samples.*) --
    a fixed-seed random sample of up to --sample-per-genome rows per
    genome, columnar (structure-of-arrays, categorical columns encoded as
    integer indices into genomes/grouping_values) to keep the standalone
    HTML small. This is the only data behind violin/histogram/KDE/
    cross-property-scatter views and the per-protein/per-CDS PCA scatter,
    which is why the dashboard UI must always caption those views as
    sampled -- even though the PCA *basis* (loadings/explained variance)
    those sampled points are projected onto is itself fit on the full
    data, exactly like Phase 5's static PCA figures.

Genome-agnostic by construction, same discipline as Phase 4/5: property
lists come from summaries_utils.property_columns() (whatever numeric
columns the master tables actually have, minus codon_-prefixed and
config.yaml's visuals.exclude_properties); grouping column NAMES and VALUE
order come from params (config.yaml's sensitivity:/group_value_order:
blocks, same as every other phase); nothing here hardcodes a genome ID,
group value, or property name.

A NOTE ON "COMBINED" PROPERTY SETS (PCA/clustering): master_protein_table
and master_cds_table are only guaranteed to join on `genome` (Phase 4's
own schema -- see merge_protein_table.py/merge_cds_table.py) -- there is
no verified, pipeline-guaranteed row-level key joining an individual
protein row to "its" CDS row (protein_id/cds_id are independently-assigned
identifiers, not a shared key). So anything that needs ROW-level values
(per-protein/per-CDS PCA, property-property correlation, cross-property
scatter) is offered per-table only (protein OR cds, never combined);
anything that operates on GENOME-level medians (per-species PCA, the
genome x property clustermap) CAN safely combine protein+CDS properties,
since genome is a real, guaranteed join key. Combined-property column
labels use "<property> (protein)"/"<property> (cds)" to disambiguate
names that exist in both tables with different meanings (e.g. `length`).
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summaries_utils import DEFAULT_EXCLUDE_PREFIXES, load_genome_labels, property_columns, resolve_group_order

QC_INT_FIELDS = {
    "n_records",
    "n_with_X",
    "n_with_internal_stop",
    "n_non_multiple_of_3",
    "n_non_ATGC",
    "n_empty",
    "n_duplicate_ids",
}


def sigfig(x, n=4):
    """Round x to n significant digits; None/NaN/Inf pass through as None
    (JSON null) -- Plotly/JS skip nulls in a trace instead of choking on a
    non-standard `NaN` token (Python's json module happily emits one, but
    it isn't valid JSON and isn't guaranteed to survive every JS
    JSON.parse)."""
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if x == 0:
        return 0.0
    digits = n - 1 - math.floor(math.log10(abs(x)))
    return round(float(x), max(digits, 0) if digits > 0 else int(digits))


def exact_stats(values: pd.Series) -> dict:
    v = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if v.size == 0:
        return dict(n=0, min=None, q1=None, median=None, q3=None, max=None, mean=None, std=None)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return dict(
        n=int(v.size),
        min=sigfig(v.min()),
        q1=sigfig(q1),
        median=sigfig(med),
        q3=sigfig(q3),
        max=sigfig(v.max()),
        mean=sigfig(v.mean()),
        std=sigfig(v.std(ddof=1)) if v.size > 1 else 0.0,
    )


def property_stats_block(df: pd.DataFrame, props: list, primary_col: str, subgroup_col: str) -> dict:
    """by_genome/by_grouping exact stats for every prop in props -- shared
    shape for both the protein and CDS property_stats blocks."""
    stats = {}
    for prop in props:
        by_genome = {g: exact_stats(sub[prop]) for g, sub in df.groupby("genome")}
        by_primary = {v: exact_stats(sub[prop]) for v, sub in df.groupby(primary_col)}
        by_subgroup = {v: exact_stats(sub[prop]) for v, sub in df.groupby(subgroup_col)}
        stats[prop] = {"by_genome": by_genome, "by_grouping": {primary_col: by_primary, subgroup_col: by_subgroup}}
    return stats


def parse_qc_report(path: str) -> dict:
    """Parse qc.py's fixed plain-text schema (see workflow/scripts/qc.py)
    into the same {protein: {...}, cds: {...}} shape as qc_protein()/
    qc_cds()'s return values, plus record_count_diff/passed. A .qc.done
    file only ever exists on disk if qc.py's own rule succeeded (it
    sys.exit()s before writing the file on any failure), so passed is
    expected True for every file this script is given -- parsed anyway, to
    fail loudly rather than assume, if that ever changes."""
    with open(path) as fh:
        lines = [line.rstrip("\n") for line in fh]

    section = None
    protein, cds = {}, {}
    record_count_diff = None
    passed = None
    for line in lines:
        if line.startswith("protein file:"):
            section = protein
        elif line.startswith("cds file:"):
            section = cds
        elif line.startswith("protein/CDS record count difference:"):
            record_count_diff = int(line.split(":", 1)[1].strip())
        elif line.strip() == "PASSED":
            passed = True
        elif line.startswith("FAILED"):
            passed = False
        elif line.startswith("  ") and ":" in line and section is not None:
            key, val = line.strip().split(":", 1)
            key, val = key.strip(), val.strip()
            if key in QC_INT_FIELDS:
                section[key] = int(val)

    if passed is None:
        raise AssertionError(f"{path}: couldn't find a PASSED/FAILED status line")
    return {"protein": protein, "cds": cds, "record_count_diff": record_count_diff, "passed": passed}


def sample_table(df: pd.DataFrame, props: list, genome_index_of: dict, primary_col, primary_order, subgroup_col, subgroup_order, seed: int, per_genome: int) -> dict:
    """Fixed-seed per-genome random sample, columnar, categorical columns
    as integer indices -- shared shape for both the protein and CDS
    sample blocks. A fresh RandomState per call (rather than a shared/
    threaded-through one) so the CDS sample's draw sequence never depends
    on how many draws the protein sample already made -- both are
    independently reproducible from the same seed."""
    rng = np.random.RandomState(seed)
    parts = []
    for genome, sub in df.groupby("genome"):
        k = min(per_genome, len(sub))
        idx = rng.choice(sub.index.to_numpy(), size=k, replace=False)
        parts.append(sub.loc[idx])
    sampled = pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0]
    return {
        "genome_index": [genome_index_of[g] for g in sampled["genome"]],
        f"{primary_col}_index": [primary_order.index(v) for v in sampled[primary_col]],
        f"{subgroup_col}_index": [subgroup_order.index(v) for v in sampled[subgroup_col]],
        "properties": {prop: [sigfig(v) for v in sampled[prop]] for prop in props},
    }, sampled


def pca_fit(matrix: pd.DataFrame, n_components=None):
    """PCA via standardize (z-score, population std -- matches sklearn's
    StandardScaler default ddof=0) + SVD -- mathematically identical to
    Phase 5's StandardScaler+sklearn.decomposition.PCA (plot_pca.py), but
    implemented in plain numpy so this script doesn't need to add
    scikit-learn as a new dependency just for the dashboard. Returns
    (mean, std, components (k x p), explained_variance_ratio (k,)),
    dropping any row with a NaN in `matrix` first (same as plot_pca.py's
    `.dropna()`).
    """
    clean = matrix.dropna()
    x = clean.to_numpy(dtype=float)
    n_samples, n_features = x.shape
    mean = x.mean(axis=0)
    std = x.std(axis=0, ddof=0)
    std_safe = np.where(std == 0, 1.0, std)
    xs = (x - mean) / std_safe
    k = n_components or min(n_samples, n_features, 10)
    k = max(min(k, n_samples, n_features), 1)
    u, s, vt = np.linalg.svd(xs, full_matrices=False)
    u, s, vt = u[:, :k], s[:k], vt[:k, :]
    scores = u * s
    explained_variance = (s**2) / max(n_samples - 1, 1)
    total_var = (xs.var(axis=0, ddof=1)).sum() if n_samples > 1 else explained_variance.sum()
    explained_variance_ratio = explained_variance / total_var if total_var else np.zeros_like(explained_variance)
    return dict(mean=mean, std=std_safe, components=vt, explained_variance_ratio=explained_variance_ratio, index=clean.index)


def pca_transform(x: np.ndarray, mean, std, components) -> np.ndarray:
    xs = (x - mean) / std
    return xs @ components.T


def pca_species_block(median_matrix: pd.DataFrame, prop_labels: list) -> dict:
    """Per-species (genome-median) PCA for one property set -- always
    exact, since it's only ever 9-ish rows (one per genome)."""
    fit = pca_fit(median_matrix)
    genomes_kept = fit["index"]
    scores = pca_transform(median_matrix.loc[fit["index"]].to_numpy(dtype=float), fit["mean"], fit["std"], fit["components"])
    return {
        "properties": prop_labels,
        "explained_variance_ratio": [sigfig(v) for v in fit["explained_variance_ratio"]],
        "loadings": {prop_labels[i]: [sigfig(v) for v in fit["components"][:, i]] for i in range(len(prop_labels))},
        "scores": {g: [sigfig(v) for v in scores[j]] for j, g in enumerate(genomes_kept)},
    }


def pca_sample_block(full_df: pd.DataFrame, props: list, sampled_df: pd.DataFrame, genome_index_of, primary_col, primary_order, subgroup_col, subgroup_order) -> dict:
    """Per-record PCA (protein-only or CDS-only): fit on the FULL table
    (exact loadings/explained variance, matching Phase 5's per-protein
    PCA exactly), then project only the already-downsampled rows for the
    scatter -- reusing the existing sample rather than shipping every
    point. Rows with a NaN in any prop are dropped from both the fit and
    the projection (same `.dropna()` discipline as Phase 5/pca_fit)."""
    fit = pca_fit(full_df[props])
    sample_clean = sampled_df.dropna(subset=props)
    scores = pca_transform(sample_clean[props].to_numpy(dtype=float), fit["mean"], fit["std"], fit["components"])
    return {
        "properties": props,
        "explained_variance_ratio": [sigfig(v) for v in fit["explained_variance_ratio"]],
        "loadings": {props[i]: [sigfig(v) for v in fit["components"][:, i]] for i in range(len(props))},
        "genome_index": [genome_index_of[g] for g in sample_clean["genome"]],
        f"{primary_col}_index": [primary_order.index(v) for v in sample_clean[primary_col]],
        f"{subgroup_col}_index": [subgroup_order.index(v) for v in sample_clean[subgroup_col]],
        "scores": [[sigfig(v) for v in row] for row in scores],
    }


def spearman_matrix(df: pd.DataFrame, props: list) -> pd.DataFrame:
    return df[props].corr(method="spearman")


def cluster_order(matrix: np.ndarray) -> list:
    """Hierarchical clustering leaf order (average linkage, euclidean --
    seaborn.clustermap's defaults, matching Phase 5's plot_clustering.py)
    over the ROWS of `matrix`. Falls back to the given order for <=2 rows
    (scipy's linkage needs at least 2 observations, and a 2-leaf
    dendrogram carries no ordering information worth computing)."""
    n = matrix.shape[0]
    if n <= 2:
        return list(range(n))
    z = linkage(matrix, method="average", metric="euclidean")
    return dendrogram(z, no_plot=True)["leaves"]


def dendrogram_coords(matrix: np.ndarray) -> dict:
    n = matrix.shape[0]
    if n <= 2:
        return {"icoord": [], "dcoord": [], "leaves": list(range(n))}
    z = linkage(matrix, method="average", metric="euclidean")
    d = dendrogram(z, no_plot=True)
    return {"icoord": d["icoord"], "dcoord": d["dcoord"], "leaves": d["leaves"]}


def codon_usage_block(cds_df: pd.DataFrame, genome_order: list) -> dict:
    """Relative codon-usage frequency (% of all codons in that genome) per
    genome, from the raw per-CDS codon_<TRIPLET> counts (see
    cds_properties.py) summed across every CDS in the genome. This is
    plain frequency, not RSCU (relative synonymous codon usage) -- RSCU
    needs a codon->amino-acid table this pipeline doesn't currently carry
    anywhere, so plain frequency is the honest thing to compute rather
    than approximating RSCU with a hardcoded genetic-code table. Codon
    names are sorted alphabetically (not by a hardcoded biological codon
    table) to stay column-name-agnostic -- this only assumes the
    `codon_` prefix convention already established elsewhere in this
    pipeline (DEFAULT_EXCLUDE_PREFIXES)."""
    codon_cols = sorted(c for c in cds_df.columns if c.startswith(DEFAULT_EXCLUDE_PREFIXES))
    codon_names = [c[len("codon_"):] for c in codon_cols]
    totals = cds_df.groupby("genome")[codon_cols].sum()
    frequencies = {}
    for genome in genome_order:
        if genome not in totals.index:
            frequencies[genome] = [None] * len(codon_names)
            continue
        row = totals.loc[genome]
        grand_total = row.sum()
        frequencies[genome] = [sigfig(100 * v / grand_total) if grand_total else None for v in row]
    return {"codons": codon_names, "frequencies": frequencies}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
    parser.add_argument("--master-cds-table", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--qc-reports", nargs="+", required=True)
    parser.add_argument("--effect-sizes-primary", required=True)
    parser.add_argument("--effect-sizes-subgroup", required=True)
    parser.add_argument("--sensitivity-leave-one-out", required=True)
    parser.add_argument("--primary-grouping", required=True)
    parser.add_argument("--subgroup-column", required=True)
    parser.add_argument("--primary-order", nargs="*", default=[])
    parser.add_argument("--subgroup-order", nargs="*", default=[])
    parser.add_argument("--exclude-properties", nargs="*", default=[])
    parser.add_argument("--sample-per-genome", type=int, required=True)
    parser.add_argument("--sample-seed", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    primary_col, subgroup_col = args.primary_grouping, args.subgroup_column
    group_columns = [primary_col, subgroup_col]

    protein_df = pd.read_csv(args.master_protein_table)
    cds_df = pd.read_csv(args.master_cds_table)
    genome_labels = load_genome_labels(args.genomes_tsv, group_columns)

    props = property_columns(protein_df, label_columns=group_columns, exclude=args.exclude_properties)
    cds_props = property_columns(cds_df, label_columns=group_columns, exclude=args.exclude_properties)

    primary_order = resolve_group_order(genome_labels[primary_col].dropna().unique(), args.primary_order or None)
    subgroup_order = resolve_group_order(genome_labels[subgroup_col].dropna().unique(), args.subgroup_order or None)

    # ---- genomes (ordered by subgroup order, then genome id) ----
    subgroup_rank = {v: i for i, v in enumerate(subgroup_order)}
    genome_labels = genome_labels.copy()
    genome_labels["_rank"] = genome_labels[subgroup_col].map(subgroup_rank)
    genome_labels = genome_labels.sort_values(["_rank", "genome"])
    genome_order = genome_labels["genome"].tolist()
    genome_index_of = {g: i for i, g in enumerate(genome_order)}

    n_proteins_by_genome = protein_df.groupby("genome").size()
    genomes = [
        {
            "genome": row["genome"],
            primary_col: row[primary_col],
            subgroup_col: row[subgroup_col],
            "n_proteins": int(n_proteins_by_genome.get(row["genome"], 0)),
        }
        for _, row in genome_labels.iterrows()
    ]

    # ---- QC ----
    qc = {}
    for path in args.qc_reports:
        genome = os.path.basename(path).split(".qc.done")[0]
        report = parse_qc_report(path)
        header_genome = None
        with open(path) as fh:
            first_line = fh.readline().strip()
        if first_line.startswith("QC report:"):
            header_genome = first_line.split(":", 1)[1].strip()
        if header_genome and header_genome != genome:
            raise AssertionError(f"{path}: filename implies genome '{genome}' but report header says '{header_genome}'")
        qc[genome] = report
    n_cds_total = sum(r["cds"].get("n_records", 0) for r in qc.values())

    # ---- property_stats: exact, from the FULL data ----
    property_stats = property_stats_block(protein_df, props, primary_col, subgroup_col)
    cds_property_stats = property_stats_block(cds_df, cds_props, primary_col, subgroup_col)

    # ---- effect sizes: pass through both groupings, all non-codon properties ----
    def load_effect_sizes(path):
        df = pd.read_csv(path)
        df = df[~df["property"].isin(args.exclude_properties) & ~df["property"].str.startswith(DEFAULT_EXCLUDE_PREFIXES)]
        rows = []
        for _, r in df.iterrows():
            rows.append(
                {
                    "property": r["property"],
                    "table": r["table"],
                    "group_a": r["group_a"],
                    "group_b": r["group_b"],
                    "n_a": int(r["n_a"]),
                    "n_b": int(r["n_b"]),
                    "median_a": sigfig(r["median_a"]),
                    "median_b": sigfig(r["median_b"]),
                    "p_value": sigfig(r["p_value"], 3),
                    "cles": sigfig(r["cles"]),
                    "rank_biserial": sigfig(r["rank_biserial"]),
                }
            )
        return rows

    effect_sizes = {
        primary_col: load_effect_sizes(args.effect_sizes_primary),
        subgroup_col: load_effect_sizes(args.effect_sizes_subgroup),
    }

    default_property = props[0] if props else None
    protein_primary_rows = [r for r in effect_sizes[primary_col] if r["table"] == "protein" and r["rank_biserial"] is not None]
    if protein_primary_rows:
        default_property = max(protein_primary_rows, key=lambda r: abs(r["rank_biserial"]))["property"]

    # ---- sensitivity: pass through, all non-codon properties ----
    sens_df = pd.read_csv(args.sensitivity_leave_one_out)
    sens_df = sens_df[
        ~sens_df["property"].isin(args.exclude_properties) & ~sens_df["property"].str.startswith(DEFAULT_EXCLUDE_PREFIXES)
    ]
    sensitivity_rows = [
        {
            "property": r["property"],
            "table": r["table"],
            "excluded_subgroup": r["excluded_subgroup"],
            "rank_biserial_full": sigfig(r["rank_biserial_full"]),
            "rank_biserial_excluded": sigfig(r["rank_biserial_excluded"]),
            "shrinkage": sigfig(r["shrinkage"]),
        }
        for _, r in sens_df.iterrows()
    ]

    # ---- downsampled per-protein / per-CDS sample, columnar ----
    samples, sampled_protein = sample_table(
        protein_df, props, genome_index_of, primary_col, primary_order, subgroup_col, subgroup_order,
        args.sample_seed, args.sample_per_genome,
    )
    cds_samples, sampled_cds = sample_table(
        cds_df, cds_props, genome_index_of, primary_col, primary_order, subgroup_col, subgroup_order,
        args.sample_seed, args.sample_per_genome,
    )

    # ---- codon usage: exact, from the FULL data ----
    codon_usage = codon_usage_block(cds_df, genome_order)

    # ---- PCA: per-species (exact, protein/cds/combined) + per-record sample (protein-only, cds-only) ----
    protein_median = pd.DataFrame({p: {g: property_stats[p]["by_genome"][g]["median"] for g in genome_order} for p in props}).reindex(genome_order)
    cds_median = pd.DataFrame({p: {g: cds_property_stats[p]["by_genome"][g]["median"] for g in genome_order} for p in cds_props}).reindex(genome_order)
    combined_labels = [f"{p} (protein)" for p in props] + [f"{p} (cds)" for p in cds_props]
    combined_median = pd.concat([protein_median.set_axis([f"{p} (protein)" for p in props], axis=1),
                                 cds_median.set_axis([f"{p} (cds)" for p in cds_props], axis=1)], axis=1)

    pca = {
        "species": {
            "protein": pca_species_block(protein_median, props),
            "cds": pca_species_block(cds_median, cds_props),
            "combined": pca_species_block(combined_median, combined_labels),
        },
        "protein_sample": pca_sample_block(protein_df, props, sampled_protein, genome_index_of, primary_col, primary_order, subgroup_col, subgroup_order),
        "cds_sample": pca_sample_block(cds_df, cds_props, sampled_cds, genome_index_of, primary_col, primary_order, subgroup_col, subgroup_order),
    }

    # ---- clustering: genome x property (combined, z-scored) + property-property correlation (protein, cds) ----
    combined_z = (combined_median - combined_median.mean()) / combined_median.std(ddof=1).replace(0, 1)
    combined_z = combined_z.reindex(genome_order)
    genome_dendro = dendrogram_coords(combined_z.to_numpy())
    genome_leaf_order = [genome_order[i] for i in genome_dendro["leaves"]] if genome_dendro["leaves"] else genome_order

    def correlation_block(df, cols):
        corr = spearman_matrix(df, cols)
        order_idx = cluster_order(corr.to_numpy())
        ordered_cols = [cols[i] for i in order_idx]
        ordered = corr.loc[ordered_cols, ordered_cols]
        return {
            "properties": ordered_cols,
            "matrix": [[sigfig(v) for v in row] for row in ordered.to_numpy()],
        }

    clustering = {
        "genome_property": {
            "properties": combined_labels,
            "genome_order": genome_leaf_order,
            "z": {g: [sigfig(v) for v in combined_z.loc[g]] for g in genome_order},
            "dendrogram": {"icoord": genome_dendro["icoord"], "dcoord": genome_dendro["dcoord"]},
        },
        "correlation": {
            "protein": correlation_block(protein_df, props) if len(props) > 1 else {"properties": props, "matrix": []},
            "cds": correlation_block(cds_df, cds_props) if len(cds_props) > 1 else {"properties": cds_props, "matrix": []},
        },
    }

    payload = {
        "meta": {
            "primary_grouping": primary_col,
            "subgroup_column": subgroup_col,
            "sample_per_genome": args.sample_per_genome,
            "sample_seed": args.sample_seed,
            "n_genomes": len(genomes),
            "n_proteins_total": int(len(protein_df)),
            "n_cds_total": int(n_cds_total),
            "n_sampled_total": int(len(sampled_protein)),
            "n_cds_sampled_total": int(len(sampled_cds)),
            "properties": props,
            "cds_properties": cds_props,
            "default_property": default_property,
            "design_note": (
                "PCA/clustering combine protein+CDS properties only at per-species (genome-median) "
                "resolution -- per-protein/per-CDS views and property-property correlation stay "
                "within one table, since there's no verified per-gene join between the two master "
                "tables in this pipeline's schema (see build_dashboard_data.py's docstring)."
            ),
        },
        "genomes": genomes,
        "grouping_values": {primary_col: primary_order, subgroup_col: subgroup_order},
        "qc": qc,
        "property_stats": property_stats,
        "cds_property_stats": cds_property_stats,
        "effect_sizes": effect_sizes,
        "sensitivity": {"excluded_subgroups": subgroup_order, "rows": sensitivity_rows},
        "samples": samples,
        "cds_samples": cds_samples,
        "codon_usage": codon_usage,
        "pca": pca,
        "clustering": clustering,
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"), allow_nan=False)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(
        f"wrote dashboard data payload -> {args.output} "
        f"({len(genomes)} genomes, {len(props)} protein properties, {len(cds_props)} cds properties, "
        f"{len(sampled_protein)} sampled proteins, {len(sampled_cds)} sampled cds, {size_mb:.2f} MB)"
    )


if __name__ == "__main__":
    main()
