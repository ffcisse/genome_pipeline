"""Shared helpers for the Phase 4 summaries scripts (merge_*_table.py,
species_summary.py, effect_sizes.py, sensitivity_drop_subgroup.py) --
mirrors seq_utils.py's role for the earlier phases: one place for logic
every summaries script needs, instead of copy-pasted per script.

Nothing here hardcodes a grouping column NAME or VALUE (e.g. "lifestyle",
"Galdieria") -- every caller passes those in explicitly, sourced from
config.yaml/genomes.tsv via Snakemake `params:`/CLI args. This is what makes
the pipeline's "genome-agnostic" claim actually true for the summaries
stage, not just the per-genome stages.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# Structural/identifier columns from THIS pipeline's own per-gene property
# scripts (protein_properties.py, disorder.py, cds_properties.py) -- fixed
# regardless of genome_table content, since these come from our own known
# output schemas, not user configuration. Group/label columns (whatever
# genome_table actually has -- lifestyle, lineage, your own names) are
# NOT listed here; callers pass those in separately (see
# numeric_property_columns's label_columns param), since only the Snakefile
# and config.yaml know what they're called for a given deployment.
STRUCTURAL_COLUMNS = {
    "genome",
    "protein_id",
    "cds_id",
    "sequence",
    "start_codon",
    "stop_codon",
}


def genome_table_columns(genomes_tsv_path: str) -> list[str]:
    """Every column in genome_table except genome_id -- whatever these are
    (name, lifestyle, lineage, or your own columns), they become label
    columns once merged into a master table, never numeric properties.
    Used to build the exclude list for numeric_property_columns."""
    labels = pd.read_csv(genomes_tsv_path, sep="\t")
    return [c for c in labels.columns if c != "genome_id"]


def load_genome_labels(genomes_tsv_path: str, group_columns: list[str]) -> pd.DataFrame:
    """genome_id + group_columns from genome_table, renamed to `genome` to
    match the join key every per-genome results table uses. group_columns
    is exactly the columns actually wired into this pipeline's grouping
    analyses (config.yaml's sensitivity.primary_grouping/subgroup_column) --
    not necessarily every column genome_table happens to have (e.g. a
    free-text `name` column stays out of the master tables)."""
    labels = pd.read_csv(genomes_tsv_path, sep="\t")
    return labels.rename(columns={"genome_id": "genome"})[["genome", *group_columns]]


def assert_row_count(df: pd.DataFrame, expected: int, label: str) -> None:
    """Fail loudly (not a warning) if a merge didn't produce the expected
    row count -- e.g. silently dropping rows on a bad join key, or
    duplicating rows on an unexpected many-to-many match."""
    if len(df) != expected:
        raise AssertionError(f"{label}: expected {expected} rows after merge, got {len(df)}")


def numeric_property_columns(df: pd.DataFrame, label_columns: list[str] = ()) -> list[str]:
    """Every column that's a real analysis property: numeric or boolean,
    excluding STRUCTURAL_COLUMNS and the given label_columns (pass
    genome_table_columns(genomes_tsv_path), or the specific grouping
    columns a script already knows it merged in). Column-name-agnostic on
    purpose -- new properties added to upstream phases show up here
    automatically."""
    exclude = STRUCTURAL_COLUMNS | set(label_columns)
    cols = []
    for col in df.columns:
        if col in exclude:
            continue
        if pd.api.types.is_numeric_dtype(df[col]) or pd.api.types.is_bool_dtype(df[col]):
            cols.append(col)
    return cols


# Column-name PATTERN that's structurally never an individually-plottable/
# -displayable "property," regardless of dataset: cds_properties.py's own
# fixed naming convention for its 64 raw per-codon count columns (a pipeline
# schema fact, like STRUCTURAL_COLUMNS above -- not a project-specific
# exclusion). A matrix of 64 near-identical codon columns isn't a
# meaningful individual-property view (the Phase 5 correlation heatmap/PCA
# already cover it wholesale), and pooling them into a "top effect size"
# ranking lets them dominate purely by count. Extendable/overridable via
# config.yaml's visuals.exclude_properties for anything else.
DEFAULT_EXCLUDE_PREFIXES = ("codon_",)


def property_columns(
    df: pd.DataFrame,
    label_columns: list[str] = (),
    exclude=(),
    exclude_prefixes=DEFAULT_EXCLUDE_PREFIXES,
) -> list[str]:
    """numeric_property_columns()'s selection, further minus any exact
    names in `exclude` (config.yaml's visuals.exclude_properties) and
    anything matching exclude_prefixes. Shared by Phase 5's plot_*.py and
    Phase 6's dashboard data builder -- one property list, used
    everywhere a human picks/views "a property.\""""
    cols = numeric_property_columns(df, label_columns=label_columns)
    exclude_set = set(exclude)
    return [c for c in cols if c not in exclude_set and not any(c.startswith(p) for p in exclude_prefixes)]


def resolve_group_order(observed_values, configured_order=None) -> list:
    """Preferred order for a grouping column's distinct values -- decides
    which value is "A" (vs "B") in cles/rank_biserial, i.e. which direction
    an effect size is signed, and (for >2 values) the order pairwise
    comparisons are generated in. This is an analyst framing choice (which
    group are you testing FOR) rather than something inferable from the
    data alone, so it's an optional override -- config.yaml's
    group_value_order -- falling back to alphabetical order for any column
    (or any value not covered) it doesn't mention.

    Raises if a configured order doesn't exactly match the observed values
    (a typo or a stale config after genome_table changes should fail loudly,
    not silently drop/duplicate a group).
    """
    observed = set(observed_values)
    if configured_order:
        configured = list(configured_order)
        if set(configured) != observed:
            raise AssertionError(
                f"configured group order {configured} doesn't match the observed values {sorted(observed)}"
            )
        return configured
    return sorted(observed)


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
