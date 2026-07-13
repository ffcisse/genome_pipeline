#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase2b
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=logs/phase2b_%j.out
#SBATCH --error=logs/phase2b_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1/2a's jobs read). Submit via the wrapper:
#
#   workflow/scripts/submit_phase2b.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase2b.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Phase 2b: intrinsic disorder prediction (rule `disorder`, see
# workflow/rules/disorder.smk), one metapredict (PyTorch, CPU inference) pass
# per genome over protein_properties.csv's sequences. Fundamentally heavier
# than Phase 2a's protein_properties -- earlier manual testing measured
# ~2h/genome single-threaded, ~30-40 min total across all 9 genomes on a full
# 64-core node -- hence --exclusive/--mem=0 (whole node, all its memory)
# instead of Phase 2a's small --cpus-per-task=4/--mem=16G.
#
# `disorder` is NOT part of `rule all` (the pipeline's default target) --
# it's additive alongside Phase 2a, not wired into the summaries/dashboard
# chain yet, so a plain `snakemake` run doesn't pick it up. This job asks for
# each genome's results/disorder/<genome>/disorder.csv explicitly instead.
# Requires protein_properties.csv to already exist for every genome (i.e.
# Phase 2a already run) -- if any genome is missing it, Snakemake will build
# it first from the existing qc/parse/protein_properties rules, same as any
# other dependency.

set -euo pipefail

cd "$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

# sbatch jobs don't always inherit a login shell's PATH/conda setup. Prefer
# whatever's already on PATH (e.g. if your site propagates the submitting
# shell's environment); only fall back to sourcing conda directly if
# `snakemake` isn't found. CONDA_BASE defaults to this deployment's anaconda3
# install -- override it (env var, or edit below) if your site sets up conda
# differently (modules, miniforge, a shared install elsewhere, etc.).
if ! command -v snakemake >/dev/null 2>&1; then
    CONDA_BASE="${CONDA_BASE:-$HOME/anaconda3}"
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate base
fi

echo "Job $SLURM_JOB_ID starting on $(hostname) at $(date)"
echo "CPUs available: $(nproc)"

# Genome-agnostic, same as the rest of the pipeline: read genome IDs from
# config/genomes.tsv rather than hardcoding them, and ask Snakemake for each
# genome's disorder.csv directly (since it isn't behind `rule all`/dashboard).
GENOMES=$(tail -n +2 config/genomes.tsv | cut -f1)
TARGETS=()
for g in $GENOMES; do
    TARGETS+=("results/disorder/${g}/disorder.csv")
done

# Make sure the env exists before we try to warm or use it (cheap/no-op if
# already built -- this just reuses the cached .snakemake/conda/<hash>/ dir).
snakemake --use-conda --conda-create-envs-only --cores 1 "${TARGETS[@]}"

# Warm the disorder env's filesystem cache BEFORE Snakemake's real,
# time-pressured run touches it -- earlier attempts at this job hit a
# ModuleNotFoundError for torch on a freshly allocated --exclusive node that
# had never read this (much bigger than analysis.yaml's) env before. A
# single, serial, no-time-pressure import here forces that first read to
# happen (and fail loudly, before the real run) instead of racing
# Snakemake's job.
#
# LD_LIBRARY_PATH=$DISORDER_ENV/lib: same fix as disorder.smk's shell
# command (see the comment there) -- pip-installed torch has no RPATH back
# to this conda env, so without this, the dynamic linker can resolve
# libstdc++.so.6 to DORI's older system copy (GLIBCXX_3.4.25 max) instead of
# this env's own newer one (GLIBCXX_3.4.34), and numpy/scipy (pulled in by
# metapredict) then crash needing GLIBCXX_3.4.29+. Set it here too so this
# warm-up check actually exercises the same working configuration the real
# rule uses, rather than giving a false pass/fail on a different one.
DISORDER_ENV=$(
    snakemake --use-conda --list-conda-envs "${TARGETS[@]}" 2>/dev/null \
        | awk -F'\t' '$1 == "workflow/envs/disorder.yaml" {print $3}'
)
echo "Warming disorder conda env: ${DISORDER_ENV}"
LD_LIBRARY_PATH="${DISORDER_ENV}/lib:${LD_LIBRARY_PATH:-}" \
    "${DISORDER_ENV}/bin/python" -c "import torch, metapredict, pandas, numpy; print('disorder env warm-up OK')"

# --retries: same transient-conda-activation insurance added to Phase 2a's
# job (see run_phase2a.sbatch) -- cheap, and doesn't mask a real,
# reproducible failure (which would just fail again on the retry). Kept even
# with the warm-up above, since it's still free insurance against unrelated
# transient blips.
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda --retries 1 "${TARGETS[@]}"

echo "Done at $(date)"
