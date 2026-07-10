#!/usr/bin/env python
"""Stage 3 (Phase 2a) -- Lightweight per-protein physicochemical properties.

Every property function below is ported unchanged (same formulas, same
edge-case handling) from resources/01_proteome_overview_FINAL.ipynb section
1.3 ("Helper functions"), validated there against the real red-algae
proteomes. The application pattern in compute_properties() also mirrors the
notebook's Section 3 exactly (same columns, same order, same
df["sequence"].apply(...) calls), to minimize any chance of silently
diverging from the validated notebook output.

Do not change formulas/thresholds here without updating and re-validating
against the notebook -- see the Phase 2a commit messages for proposed
(NOT applied) accuracy/efficiency improvements to consider separately.

Deliberately excludes intrinsic disorder prediction (heavy, Phase 2b, its
own SLURM rule) and CDS-derived properties (GC/GC3/codon usage/ENC, a
separate cds_properties rule).

Property groups (built incrementally, one per Phase 2a commit):
  - Core: pct_* amino acid composition, pI, gravy (length already exists
    from Phase 1's parse step)
  - Extended: aliphatic_index, thermostable_fraction, instability_index,
    net_charge_pH7, charge_density, cysteine_fraction, carbon_oxidation_state
  - Aggregation: agg_mean_a3v, agg_Na4vSS, agg_hotspot_fraction
"""

import argparse
import os
import sys

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_utils import clean_sequence

# ============================================================================
# Core group
# ============================================================================

# ---------- Amino acid composition ----------
aa_groups = {
    "charged": set("DEKRH"),
    "acidic": set("DE"),
    "basic": set("KRH"),
    "hydrophobic": set("AVILMFWY"),
    "polar": set("STNQ"),
    "aromatic": set("FWY"),
    "special": set("GPC"),
}


def amino_acid_composition(seq: str) -> dict:
    """Return fraction of residues in each chemical group (uses raw sequence, only strips '*')."""
    seq = seq.replace("*", "")
    length = len(seq)
    if length == 0:
        return {f"pct_{g}": 0 for g in aa_groups}
    return {
        f"pct_{group}": sum(seq.count(aa) for aa in residues) / length
        for group, residues in aa_groups.items()
    }


# ---------- Core properties: pI and GRAVY ----------
def calculate_pi_gravy(seq: str) -> pd.Series:
    """Isoelectric point and GRAVY hydrophobicity via Biopython."""
    seq = clean_sequence(seq)
    if len(seq) == 0:
        return pd.Series({"pI": None, "gravy": None})
    analysis = ProteinAnalysis(seq)
    return pd.Series({"pI": analysis.isoelectric_point(), "gravy": analysis.gravy()})


# ============================================================================
# Orchestration -- mirrors the notebook's Section 3 application pattern.
# ============================================================================


def compute_properties(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    aa_comp_df = df["sequence"].apply(amino_acid_composition).apply(pd.Series)
    df = pd.concat([df, aa_comp_df], axis=1)

    df[["pI", "gravy"]] = df["sequence"].apply(calculate_pi_gravy)

    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--protein-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.protein_table)
    result = compute_properties(df)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"{args.genome}: wrote {len(result)} proteins x {result.shape[1]} cols -> {args.output}")


if __name__ == "__main__":
    main()
