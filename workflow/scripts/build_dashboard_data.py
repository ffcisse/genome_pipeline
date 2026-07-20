#!/usr/bin/env python
"""Stage 7 (Phase 6a) -- Build the JSON data payload the dashboard HTML
embeds inline.

Two different kinds of numbers go into this payload, and they must never
be confused with each other (see dashboard.smk's docstring and the
Property Explorer's on-screen caption):

  * EXACT statistics -- property_stats (per-genome and per-grouping-value
    min/q1/median/q3/max/mean/std/n) and the effect_sizes/sensitivity
    pass-throughs -- computed here from the FULL master_protein_table (all
    ~61k rows for this deployment's data), never from the sample below.
    These drive the "exact" box-plot trace (Plotly's precomputed
    q1/median/q3/lowerfence/upperfence fields) and every summary number
    shown anywhere in the dashboard.

  * A DOWNSAMPLED per-protein sample (samples.*) -- a fixed-seed random
    sample of up to --sample-per-genome proteins per genome, columnar
    (structure-of-arrays, categorical columns encoded as integer indices
    into genomes/grouping_values) to keep the standalone HTML small. This
    is the only data behind violin/histogram/KDE/scatter views, which is
    why the dashboard UI must always caption those views as sampled.

Genome-agnostic by construction, same discipline as Phase 4/5: property
list comes from summaries_utils.property_columns() (whatever numeric
columns master_protein_table actually has, minus codon_-prefixed and
config.yaml's visuals.exclude_properties); grouping column NAMES and VALUE
order come from params (config.yaml's sensitivity:/group_value_order:
blocks, same as every other phase); nothing here hardcodes a genome ID,
group value, or property name.

Effect Sizes/Sensitivity include BOTH protein- and CDS-side properties
(minus codon_-prefixed), same rows Phase 5's plot_effect_sizes.py/
plot_sensitivity.py already plot -- CDS properties aren't sampled
per-protein (no CDS view yet, that's Phase 6b), so a click on a CDS-side
property in the dashboard's Effect Sizes view can't jump to a Property
Explorer plot; the dashboard JS handles that case with a short message
rather than a dead click.
"""

import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-protein-table", required=True)
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
    genome_labels = load_genome_labels(args.genomes_tsv, group_columns)

    props = property_columns(protein_df, label_columns=group_columns, exclude=args.exclude_properties)

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
    property_stats = {}
    for prop in props:
        by_genome = {g: exact_stats(sub[prop]) for g, sub in protein_df.groupby("genome")}
        by_primary = {v: exact_stats(sub[prop]) for v, sub in protein_df.groupby(primary_col)}
        by_subgroup = {v: exact_stats(sub[prop]) for v, sub in protein_df.groupby(subgroup_col)}
        property_stats[prop] = {
            "by_genome": by_genome,
            "by_grouping": {primary_col: by_primary, subgroup_col: by_subgroup},
        }

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

    # ---- downsampled per-protein sample, columnar ----
    rng = np.random.RandomState(args.sample_seed)
    sampled_parts = []
    for genome, sub in protein_df.groupby("genome"):
        k = min(args.sample_per_genome, len(sub))
        idx = rng.choice(sub.index.to_numpy(), size=k, replace=False)
        sampled_parts.append(sub.loc[idx])
    sampled = pd.concat(sampled_parts, ignore_index=True) if sampled_parts else protein_df.iloc[0:0]

    samples = {
        "genome_index": [genome_index_of[g] for g in sampled["genome"]],
        f"{primary_col}_index": [primary_order.index(v) for v in sampled[primary_col]],
        f"{subgroup_col}_index": [subgroup_order.index(v) for v in sampled[subgroup_col]],
        "properties": {prop: [sigfig(v) for v in sampled[prop]] for prop in props},
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
            "n_sampled_total": int(len(sampled)),
            "properties": props,
            "phase6b_note": (
                "CDS/codon-level property views, PCA/clustering, cross-property scatter, "
                "and export/download are planned for Phase 6b."
            ),
        },
        "genomes": genomes,
        "grouping_values": {primary_col: primary_order, subgroup_col: subgroup_order},
        "qc": qc,
        "property_stats": property_stats,
        "effect_sizes": effect_sizes,
        "sensitivity": {"excluded_subgroups": subgroup_order, "rows": sensitivity_rows},
        "samples": samples,
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"), allow_nan=False)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(
        f"wrote dashboard data payload -> {args.output} "
        f"({len(genomes)} genomes, {len(props)} properties, {len(sampled)} sampled proteins, {size_mb:.2f} MB)"
    )


if __name__ == "__main__":
    main()
