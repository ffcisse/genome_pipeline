#!/usr/bin/env python
"""Stage 3 (Phase 2b) -- Intrinsic disorder prediction for one genome.

Unlike Phase 2a's protein_properties (pure sequence/composition math, fast,
single-core), this is a real neural-network prediction (metapredict, a
PyTorch model) -- CPU inference, ~2h/genome single-threaded on the full
61k-protein dataset in earlier manual testing, hence its own rule with its
own conda env and a real SLURM resource block (see rules/disorder.smk).

Reads Phase 2a's protein_properties.csv for the sequences and protein IDs
only -- does NOT recompute or touch any Phase 2a column. Writes a separate
results/disorder/<genome>/disorder.csv (genome + protein_id + 4 new
columns), joinable back onto protein_properties.csv on (genome, protein_id).
Kept separate rather than merged into protein_properties.csv because a
Snakemake output file should have exactly one rule that produces it -- if
this rule rewrote protein_properties.csv, Snakemake would see Phase 2a's
own output change on a run it didn't do, which breaks its ability to reason
about what's up to date.

CRITICAL: metapredict raises KeyError on non-standard residues (e.g. 'X'),
so every sequence is run through seq_utils.clean_sequence first, exactly as
every Phase 2a property function already does.

Per-protein features (definitions carried over from earlier validated
manual run):
  disorder_mean     -- mean per-residue disorder score (0-1)
  disorder_fraction -- fraction of residues with score >= 0.5 (PRIMARY metric)
  longest_idr       -- longest continuous run of residues with score >= 0.5
  n_idrs            -- count of distinct runs of score >= 0.5 that are >= 5
                        residues long
Sequences that are empty after cleaning get disorder_mean/disorder_fraction
= NaN (undefined -- no residues to average) and longest_idr/n_idrs = 0 (no
regions possible).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
from metapredict import predict_disorder

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_utils import clean_sequence

DISORDER_THRESHOLD = 0.5
MIN_IDR_LENGTH = 5


def longest_and_count_runs(mask: np.ndarray) -> tuple[int, int]:
    """Longest continuous True run in `mask`, and how many such runs are
    at least MIN_IDR_LENGTH long. Plain run-length encoding over the
    per-residue >=0.5 boolean array -- no need for metapredict's own
    return_domains machinery since these are simpler, custom thresholds."""
    longest = 0
    n_idrs = 0
    run = 0
    for is_disordered in mask:
        if is_disordered:
            run += 1
            longest = max(longest, run)
        else:
            if run >= MIN_IDR_LENGTH:
                n_idrs += 1
            run = 0
    if run >= MIN_IDR_LENGTH:
        n_idrs += 1
    return longest, n_idrs


def summarize_scores(scores: np.ndarray) -> tuple[float, float, int, int]:
    mask = scores >= DISORDER_THRESHOLD
    longest, n_idrs = longest_and_count_runs(mask)
    return float(scores.mean()), float(mask.mean()), longest, n_idrs


def compute_disorder(df: pd.DataFrame, threads: int) -> pd.DataFrame:
    torch.set_num_threads(threads)

    cleaned = df["sequence"].apply(clean_sequence)

    # Batch through metapredict as a single dict call (protein_id -> cleaned
    # sequence) rather than looping one sequence at a time: metapredict v3
    # batches/pads internally, and keying by protein_id (rather than passing
    # a bare list) means results come back tagged by ID instead of relying
    # on positional order matching the input list.
    non_empty = {pid: seq for pid, seq in zip(df["protein_id"], cleaned) if seq}
    predictions = (
        predict_disorder(non_empty, device="cpu", show_progress_bar=False)
        if non_empty
        else {}
    )

    means, fractions, longest_idrs, n_idrs_list = [], [], [], []
    for pid in df["protein_id"]:
        if pid in predictions:
            _, scores = predictions[pid]
            mean, fraction, longest, n_idrs = summarize_scores(np.asarray(scores))
        else:
            mean, fraction, longest, n_idrs = np.nan, np.nan, 0, 0
        means.append(mean)
        fractions.append(fraction)
        longest_idrs.append(longest)
        n_idrs_list.append(n_idrs)

    return pd.DataFrame(
        {
            "genome": df["genome"],
            "protein_id": df["protein_id"],
            "disorder_mean": means,
            "disorder_fraction": fractions,
            "longest_idr": longest_idrs,
            "n_idrs": n_idrs_list,
        }
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--protein-properties", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="CPU threads for torch to use -- pass the rule's {threads} so metapredict "
        "uses exactly its SLURM allocation, not the whole node.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.protein_properties, usecols=["genome", "protein_id", "sequence"])
    result = compute_disorder(df, args.threads)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"{args.genome}: wrote {len(result)} proteins x {result.shape[1]} cols -> {args.output}")


if __name__ == "__main__":
    main()
