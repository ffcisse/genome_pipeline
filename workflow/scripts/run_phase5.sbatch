#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase5
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase5_%j.out
#SBATCH --error=logs/phase5_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1/2a/2b/3/4's jobs read). Submit via the
# wrapper:
#
#   workflow/scripts/submit_phase5.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase5.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Phase 5: static visualization generation -- boxplots, distributions, PCA,
# clustering, effect-size forest plots, and the leave-one-out sensitivity
# heatmap, all read from Phase 4's already-computed summary tables (rules
# `plot_boxplots`/`plot_distributions`/`plot_cds_distributions`/`plot_pca`/
# `plot_clustering`/`plot_effect_sizes`/`plot_sensitivity` -- see
# workflow/rules/visuals.smk). Light -- plotting ~61k rows with matplotlib/
# seaborn, no heavy model like Phase 2b's disorder -- same small allocation
# as Phase 2a/3/4.
#
# Like Phase 3/4's rules, these ARE part of `rule all` (they replaced the
# touch-only `visuals` stub's real work, though `visuals` itself stays a
# thin aggregator -- see visuals.smk), so a plain `snakemake` invocation
# picks them up naturally alongside dashboard (still a touch-only stub).

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
# 2a/2b/3/4's jobs -- cheap, and doesn't mask a real, reproducible failure
# (which would just fail again on the retry).
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda --retries 1

echo "Done at $(date)"
