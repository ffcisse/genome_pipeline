#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase1
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/phase1_%j.out
#SBATCH --error=logs/phase1_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead. Submit via the wrapper, which reads that block and passes
# them as sbatch CLI flags (CLI flags override/supply directives the same as
# if they'd been written in this header):
#
#   workflow/scripts/submit_phase1.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase1.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Runs the full pipeline as currently built (stage_inputs -> qc -> parse, real
# logic; protein_properties/cds_properties/summaries/visuals/dashboard are
# still `touch`-only stubs from the Phase 0 scaffold -- cheap, included here
# just because `rule all` (the default target) pulls in the whole DAG).
#
# All 9 genomes' QC + parsing is light (small proteomes, single-pass FASTA
# parsing) -- this does NOT need --exclusive or a big allocation. The
# generous --time is mostly headroom for --use-conda's one-time env solve
# (workflow/envs/analysis.yaml) on first run; it's cached after that.

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
