#!/bin/bash
#SBATCH --job-name=genome_pipeline_phase1
#SBATCH -q jgi_normal
#SBATCH -A grp-org-sc-fungal
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/phase1_%j.out
#SBATCH --error=logs/phase1_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ffcisse@berkeley.edu

# Runs the full pipeline as currently built (stage_inputs -> qc -> parse, real
# logic; protein_properties/cds_properties/summaries/visuals/dashboard are
# still `touch`-only stubs from the Phase 0 scaffold -- cheap, included here
# just because `rule all` (the default target) pulls in the whole DAG).
#
# All 9 genomes' QC + parsing is light (small proteomes, single-pass FASTA
# parsing) -- this does NOT need --exclusive or a big allocation. The
# generous --time is mostly headroom for --use-conda's one-time env solve
# (workflow/envs/analysis.yaml) on first run; it's cached after that.
#
# Submit from the repo root:
#   sbatch workflow/scripts/run_phase1.sbatch

set -euo pipefail

source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

cd "$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

echo "Job $SLURM_JOB_ID starting on $(hostname) at $(date)"
echo "CPUs available: $(nproc)"

snakemake --cores "${SLURM_CPUS_PER_TASK}" --use-conda

echo "Done at $(date)"
