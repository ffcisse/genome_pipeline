#!/usr/bin/env python
"""Stage 5 (Phase 4) -- Combine every subgroup's sensitivity_drop_{subgroup}.csv
(one per value of config.yaml's `sensitivity.subgroup_column`, computed by
the wildcarded sensitivity_drop_subgroup rule) into one leave-one-out sweep
table, tagging each with the subgroup that was excluded to produce it.

This is more informative than checking one hand-picked subgroup: it shows
how much EACH subgroup individually contributes to the apparent
primary-grouping effect, not just whichever one you thought to check.

--inputs and --subgroups are parallel lists (same order Snakemake's
expand(..., subgroup=SUBGROUPS) generated them in) rather than recovering
the excluded subgroup by parsing it back out of each filename -- subgroup
names can contain underscores (e.g. "Mesophile_lineage"), so an explicit
parallel list is unambiguous where filename-parsing could be fragile.
"""

import argparse
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--subgroups", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if len(args.inputs) != len(args.subgroups):
        raise AssertionError(f"got {len(args.inputs)} --inputs but {len(args.subgroups)} --subgroups")

    frames = []
    for path, subgroup in zip(args.inputs, args.subgroups):
        df = pd.read_csv(path)
        df.insert(0, "excluded_subgroup", subgroup)
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    combined.to_csv(args.output, index=False)
    print(f"wrote {len(combined)} rows ({len(args.subgroups)} subgroups) -> {args.output}")


if __name__ == "__main__":
    main()
