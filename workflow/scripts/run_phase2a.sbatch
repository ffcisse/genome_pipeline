#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase2a
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase2a_%j.out
#SBATCH --error=logs/phase2a_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1's job reads). Submit via the wrapper:
#
#   workflow/scripts/submit_phase2a.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase2a.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Runs the full pipeline as currently built. Real logic: stage_inputs -> qc
# -> parse -> protein_properties (Phase 2a: lightweight per-protein
# physicochemical properties -- composition, pI/GRAVY, aliphatic index,
# instability, charge, aggregation propensity; no intrinsic disorder, that's
# Phase 2b's own heavier rule). cds_properties/summaries/visuals/dashboard
# are still `touch`-only stubs from the Phase 0 scaffold -- cheap, included
# here only because `rule all` (the default target) pulls in the whole DAG.
#
# All 9 genomes' protein_properties is light and single-core per genome
# (measured ~4-5s/genome interactively on real data -- under a minute total
# for all 9) -- this does NOT need --exclusive or a big allocation.
# qc/parse already ran for all 9 genomes in the Phase 1 job, and the
# `analysis` conda env is already cached from that run, so this job should
# finish in well under the time budget below; the budget is just headroom,
# not an estimate.

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

snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda

echo "Done at $(date)"
