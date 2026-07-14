"""Shared helpers for the Phase 4 summaries scripts (merge_*_table.py,
species_summary.py, effect_sizes.py, sensitivity_drop_galdieria.py) --
mirrors seq_utils.py's role for the earlier phases: one place for logic
every summaries script needs, instead of copy-pasted per script.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# Columns that are identifiers/free text, never analysis properties --
# every summaries script excludes these when picking numeric properties.
NON_PROPERTY_COLUMNS = {
    "genome",
    "protein_id",
    "cds_id",
    "sequence",
    "start_codon",
    "stop_codon",
    "name",
    "lifestyle",
    "lineage",
}


def load_genome_labels(genomes_tsv_path: str) -> pd.DataFrame:
    """genome_id/lifestyle/lineage from config/genomes.tsv, renamed to
    `genome` to match the join key every per-genome results table uses."""
    labels = pd.read_csv(genomes_tsv_path, sep="\t")
    return labels.rename(columns={"genome_id": "genome"})[["genome", "lifestyle", "lineage"]]


def assert_row_count(df: pd.DataFrame, expected: int, label: str) -> None:
    """Fail loudly (not a warning) if a merge didn't produce the expected
    row count -- e.g. silently dropping rows on a bad join key, or
    duplicating rows on an unexpected many-to-many match."""
    if len(df) != expected:
        raise AssertionError(f"{label}: expected {expected} rows after merge, got {len(df)}")


def numeric_property_columns(df: pd.DataFrame) -> list[str]:
    """Every column that's a real analysis property: numeric or boolean,
    excluding NON_PROPERTY_COLUMNS. Column-name-agnostic on purpose -- new
    properties added to upstream phases show up here automatically."""
    cols = []
    for col in df.columns:
        if col in NON_PROPERTY_COLUMNS:
            continue
        if pd.api.types.is_numeric_dtype(df[col]) or pd.api.types.is_bool_dtype(df[col]):
            cols.append(col)
    return cols


def compare_groups(a: pd.Series, b: pd.Series) -> dict:
    """Mann-Whitney U, CLES, and rank-biserial correlation for one property
    between two groups (a, b). CLES = U / (n_a * n_b) is the probability a
    random draw from `a` exceeds a random draw from `b` (ties count as
    half) -- this is exactly scipy's U statistic for the first sample
    passed to mannwhitneyu, so `a` must always be "group A" consistently by
    caller convention (documented per-caller, e.g. lifestyle: a=extremophile,
    b=mesophile).

    With n in the tens of thousands, p-values collapse to ~0 for nearly
    every property and stop being informative about *how big* a difference
    is -- rank_biserial (2*CLES - 1, range -1..+1) is the actual effect-size
    signal this pipeline cares about.
    """
    a = pd.to_numeric(a, errors="coerce").dropna().to_numpy(dtype=float)
    b = pd.to_numeric(b, errors="coerce").dropna().to_numpy(dtype=float)
    n_a, n_b = len(a), len(b)
    if n_a == 0 or n_b == 0:
        return dict(
            n_a=n_a, n_b=n_b, median_a=np.nan, median_b=np.nan, p_value=np.nan, cles=np.nan, rank_biserial=np.nan
        )
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    cles = stat / (n_a * n_b)
    return dict(
        n_a=n_a,
        n_b=n_b,
        median_a=float(np.median(a)),
        median_b=float(np.median(b)),
        p_value=p,
        cles=cles,
        rank_biserial=2 * cles - 1,
    )
