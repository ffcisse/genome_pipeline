#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase6
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase6_%j.out
#SBATCH --error=logs/phase6_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1/2a/2b/3/4/5's jobs read). Submit via
# the wrapper:
#
#   workflow/scripts/submit_phase6.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase6.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Phase 6 (6a+6b): the interactive dashboard -- rules `dashboard_data`
# (builds the JSON payload from Phase 1/4 outputs) and `dashboard` (inlines
# that payload + the vendored Plotly bundle + this pipeline's HTML/CSS/JS
# into one standalone results/dashboard/proteome_dashboard.html). Light --
# no new heavy conda env, no per-genome compute, just pandas/numpy/scipy
# over already-computed tables -- same small allocation as Phase 2a/3/4/5.
#
# Like Phase 3/4/5's rules, dashboard_data/dashboard ARE part of `rule all`
# (rule all targets the dashboard HTML directly -- see the Snakefile), so a
# plain `snakemake` invocation picks them up naturally, along with any
# earlier stage that isn't done yet.

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

# --retries: same transient-conda-activation insurance added to Phase
# 2a/2b/3/4/5's jobs -- cheap, and doesn't mask a real, reproducible failure
# (which would just fail again on the retry).
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda --retries 1

echo "Done at $(date)"
