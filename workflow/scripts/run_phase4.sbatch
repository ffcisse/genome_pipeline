#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase4
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase4_%j.out
#SBATCH --error=logs/phase4_%j.err
#SBATCH --mail-type=BEGIN,END
#
# Account, QOS, partition, and --mail-user are deliberately NOT set here --
# they're cluster/user-specific and live in config/config.yaml's `slurm:`
# block instead (same block Phase 1/2a/2b/3's jobs read). Submit via the
# wrapper:
#
#   workflow/scripts/submit_phase4.sh
#
# Submitting this file directly (`sbatch workflow/scripts/run_phase4.sbatch`)
# still works, but without an account/QOS/mail-user unless your site has
# defaults for them -- use the wrapper unless you know your site doesn't need
# those.

# Phase 4: summaries -- merge every genome's per-gene property tables into
# master_protein_table.csv/master_cds_table.csv, species_summary.csv, and
# the configurable-grouping effect-size tables (rules `merge_protein_table`,
# `merge_cds_table`, `species_summary`, `effect_sizes`,
# `sensitivity_drop_galdieria` -- see workflow/rules/summaries.smk). Pure
# pandas/scipy over already-computed per-gene tables (61,349 rows) -- light,
# same small allocation as Phase 2a/3, measured ~40s end to end in testing.
#
# Like Phase 3's cds_properties, these rules ARE part of `rule all` (they
# replaced the old touch-only `summaries` stub's real work), so a plain
# `snakemake` invocation picks them up naturally alongside
# visuals/dashboard (still touch-only stubs -- Phase 4 is summaries only).

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
# 2a/2b/3's jobs -- cheap, and doesn't mask a real, reproducible failure
# (which would just fail again on the retry).
snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda --retries 1

echo "Done at $(date)"
