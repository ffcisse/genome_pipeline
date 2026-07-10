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
# Extended group
# ============================================================================


def aliphatic_index(seq: str):
    """Ikai (1980) aliphatic index: relative volume from A, V, I, L side chains.
    Higher => more thermally stable. Not a percentage; can exceed 100."""
    seq = clean_sequence(seq)
    n = len(seq)
    if n == 0:
        return None
    a = seq.count("A") / n * 100
    v = seq.count("V") / n * 100
    i = seq.count("I") / n * 100
    l = seq.count("L") / n * 100
    return a + 2.9 * v + 3.9 * (i + l)


THERMOSTABLE_RESIDUES = set("IVYWREL")


def thermostable_fraction(seq: str):
    """Fraction of residues in {I,V,Y,W,R,E,L} -- residues associated with heat tolerance.
    NOTE: chemically heterogeneous; interpret via per-residue decomposition, not alone."""
    seq = clean_sequence(seq)
    n = len(seq)
    if n == 0:
        return None
    return sum(seq.count(aa) for aa in THERMOSTABLE_RESIDUES) / n


def calculate_instability(seq: str):
    """Guruprasad (1990) instability index. <40 predicted stable, >=40 unstable.
    CAUTION: trained on mesophilic proteins; a weak thermostability proxy for extremophiles."""
    seq = clean_sequence(seq)
    if len(seq) == 0:
        return None
    return ProteinAnalysis(seq).instability_index()


def calculate_net_charge(seq: str, pH: float = 7.0):
    """Net charge at a given pH (default 7.0) via Biopython."""
    seq = clean_sequence(seq)
    if len(seq) == 0:
        return None
    return ProteinAnalysis(seq).charge_at_pH(pH)


def cysteine_fraction(seq: str):
    """Fraction of cysteine residues -- a rough proxy for disulfide-bond potential."""
    seq = clean_sequence(seq)
    n = len(seq)
    if n == 0:
        return None
    return seq.count("C") / n


# Elemental composition (C, H, N, O, S) of each standard amino acid's free/neutral form.
# Used for carbon oxidation state (Zc). Peptide-bond formation (loss of H2O) does not
# change Zc, so summing free-residue formulas is valid.
AA_ELEMENTAL = {
    "G": (2, 5, 1, 2, 0), "A": (3, 7, 1, 2, 0), "S": (3, 7, 1, 3, 0),
    "P": (5, 9, 1, 2, 0), "V": (5, 11, 1, 2, 0), "T": (4, 9, 1, 3, 0),
    "C": (3, 7, 1, 2, 1), "L": (6, 13, 1, 2, 0), "I": (6, 13, 1, 2, 0),
    "N": (4, 8, 2, 3, 0), "D": (4, 7, 1, 4, 0), "Q": (5, 10, 2, 3, 0),
    "K": (6, 14, 2, 2, 0), "E": (5, 9, 1, 4, 0), "M": (5, 11, 1, 2, 1),
    "H": (6, 9, 3, 2, 0), "F": (9, 11, 1, 2, 0), "R": (6, 14, 4, 2, 0),
    "Y": (9, 11, 1, 3, 0), "W": (11, 12, 2, 2, 0),
}


def carbon_oxidation_state(seq: str):
    """Average oxidation state of carbon, Zc (Dick 2014).
    Zc = (-nH + 3nN + 2nO + 2nS) / nC. Lower (more negative) => more reduced protein."""
    seq = clean_sequence(seq)
    total_C = 0
    numerator = 0
    for aa in seq:
        nC, nH, nN, nO, nS = AA_ELEMENTAL[aa]
        total_C += nC
        numerator += -nH + 3 * nN + 2 * nO + 2 * nS
    if total_C == 0:
        return None
    return numerator / total_C


# ============================================================================
# Orchestration -- mirrors the notebook's Section 3 application pattern.
# ============================================================================


def compute_properties(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    aa_comp_df = df["sequence"].apply(amino_acid_composition).apply(pd.Series)
    df = pd.concat([df, aa_comp_df], axis=1)

    df[["pI", "gravy"]] = df["sequence"].apply(calculate_pi_gravy)

    # Sequence-composition properties
    df["aliphatic_index"] = df["sequence"].apply(aliphatic_index)
    df["thermostable_fraction"] = df["sequence"].apply(thermostable_fraction)
    df["cysteine_fraction"] = df["sequence"].apply(cysteine_fraction)
    df["carbon_oxidation_state"] = df["sequence"].apply(carbon_oxidation_state)

    # Biopython-based properties
    df["instability_index"] = df["sequence"].apply(calculate_instability)
    df["net_charge_pH7"] = df["sequence"].apply(calculate_net_charge)

    # Derived -- note: divides by `length` (the Phase 1 column, i.e. raw
    # length with only the terminal stop stripped), NOT len(clean_sequence(seq)).
    # This matches the notebook exactly (df["net_charge_pH7"] / df["length"]).
    df["charge_density"] = df["net_charge_pH7"] / df["length"]

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
