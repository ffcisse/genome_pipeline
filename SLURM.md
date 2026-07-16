# Running on DORI (SLURM)

DORI uses the SLURM scheduler. This pipeline submits jobs the simple way: one `sbatch` job per
phase, running `snakemake --cores $SLURM_CPUS_PER_TASK --use-conda` locally inside a single
allocation -- not Snakemake's separate `snakemake-executor-plugin-slurm` cluster-executor model
(that was considered early on but never actually adopted; if you find references to it elsewhere,
they're stale).

This file is inherently deployment-specific -- "DORI" is this deployment's cluster, and running
elsewhere means different scheduler commands/module setup. What *does* travel with the repo is
`config/config.yaml`'s `slurm:` block: account/QOS/partition/mail_user are read from there (not
hardcoded in any script), so pointing at a different allocation on the same cluster -- or the same
allocation on a different one -- is a config edit, not a script edit.

## Submission model

Each phase has a `submit_phaseN.sh` wrapper around a `run_phaseN.sbatch` script:

```bash
workflow/scripts/submit_phase1.sh   # stage_inputs -> qc -> parse
workflow/scripts/submit_phase2a.sh  # protein_properties
workflow/scripts/submit_phase2b.sh  # disorder (heavy -- see below)
workflow/scripts/submit_phase3.sh   # cds_properties
workflow/scripts/submit_phase4.sh   # merge_*/species_summary/effect_sizes/sensitivity_*
```

Submit via the wrapper, not `sbatch run_phaseN.sbatch` directly, unless your site has defaults for
account/QOS/mail-user: each wrapper reads `config.yaml`'s `slurm:` block and passes it as `sbatch`
CLI flags. See the main [README](README.md#running-it) for the full per-phase resource table and
what each phase produces -- this file only covers what's genuinely SLURM/DORI-specific.

Each `run_phaseN.sbatch` script requests one allocation and runs `snakemake --cores
$SLURM_CPUS_PER_TASK --use-conda` (Phase 2b's `disorder` targets are named explicitly, since that
output isn't part of `rule all`; every other phase runs Snakemake with no explicit target, so it
also picks up any earlier stage that isn't done yet). Re-running an earlier phase's script after
later work is already done is harmless -- Snakemake just confirms everything's up to date.

**Phase 2b (`disorder`) needs real resources**, unlike the other phases: a full 64-core exclusive
node (see `run_phase2b.sbatch`'s `#SBATCH` header). It loads a real PyTorch model and runs
inference over every protein -- do not try to run it inline on a login node or with a small
allocation.

Logs land in `logs/phaseN_<jobid>.{out,err}`.

## `--use-conda` inside SLURM jobs: test before you trust it

More than once on this deployment, a conda env that worked fine when tested interactively on a
login node broke *inside an actual SLURM job* -- different root causes each time (a cold/
first-touch environment on a freshly allocated node; a `libstdc++` ABI conflict between
pip-installed PyTorch and conda-forge-built numpy/scipy -- see `workflow/rules/disorder.smk`'s
shell-command comment for the full story on that one). **Whenever you add or modify a conda env,
submit a real (even tiny) SLURM job that activates it and runs something simple before trusting
it in a full run.** A clean login-node test is not sufficient evidence.

## Per-rule `resources:`

Snakemake rules can declare a `resources:` block (`mem_mb`, `runtime`, etc.) -- `disorder` in
`workflow/rules/disorder.smk` does, documenting what it actually needs even though nothing enforces
it under the current local-`--cores` submission model (that's what a workflow profile, below,
would be for). The other rules don't need one; they're light enough that SLURM's/the sbatch
script's own defaults are enough.

## Not set up: a workflow profile / cluster-executor model

If per-rule SLURM submission (one job per rule instance, submitted and tracked by Snakemake
itself via `snakemake --executor slurm --jobs N` and a workflow profile) ever becomes worth the
overhead -- e.g. many genomes where phases should run genuinely in parallel across separate
nodes rather than one allocation at a time -- that's a real alternative to the model above, but
it isn't set up in this repo. Don't assume it works without building and testing it first.
