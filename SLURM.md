# Running on DORI (SLURM)

DORI uses the SLURM scheduler. Snakemake submits jobs to SLURM through the
`snakemake-executor-plugin-slurm` executor plugin (Snakemake 8+ moved cluster support to
plugins; this repo uses Snakemake 9).

## One-time setup (on DORI, inside your snakemake conda env)
```bash
pip install snakemake-executor-plugin-slurm
```

## Per-rule resources
A rule can declare a `resources:` block (e.g. `slurm_partition`, `mem_mb`, `runtime`,
`cpus_per_task`) that the SLURM executor reads when submitting that rule's jobs. None are set
in this scaffold yet -- Phase 1 will add real values once we know actual runtime/memory needs
per stage. Until then, SLURM's own defaults apply.

## Submit command
```bash
snakemake --executor slurm --jobs 20 --use-conda
```
- `--executor slurm` -- submit each job as its own SLURM job instead of running locally.
- `--jobs 20` -- allow up to 20 jobs queued/running at once.
- `--use-conda` -- activate each rule's declared conda env (workflow/envs/) before it runs.

## Optional: a workflow profile
Instead of retyping flags every time, you can put them in `config/slurm_profile/config.yaml`
and run `snakemake --profile config/slurm_profile`. Not set up yet in this scaffold --
worth adding once per-rule resource requirements stabilize in Phase 1+.
