#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase3
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase3_%j.out
#SBATCH --error=logs/phase3_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1/2a/2b's jobs read). Submit via the
# wrapper:
#
#   workflow/scripts/submit_phase3.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase3.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Phase 3: CDS/codon-usage properties (rule `cds_properties`, see
# workflow/rules/cds_properties.smk), via codonW -- pure C-program codon
# counting on already-parsed sequences, same cost profile as Phase 2a's
# protein_properties (fast, single-core-per-genome, ~seconds/genome
# interactively) -- NOT Phase 2b's disorder rule (PyTorch, its own exclusive
# node). Small allocation here, same as run_phase2a.sbatch.
#
# Unlike Phase 2b's disorder rule, cds_properties IS part of `rule all` (the
# pipeline's default target) -- it replaced a Phase 0 touch-only stub that
# summaries/visuals/dashboard already depended on, so a plain `snakemake`
# invocation naturally picks it up alongside those (still touch-only)
# downstream stubs, with no separate target list needed the way Phase 2b's
# job needs one for `disorder`.

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

# --retries: same transient-conda-activation insurance added to Phase 2a/2b's
# jobs -- cheap, and doesn't mask a real, reproducible failure (which would
# just fail again on the retry).
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda --retries 1

echo "Done at $(date)"
