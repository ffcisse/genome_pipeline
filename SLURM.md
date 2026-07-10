# Running on DORI (SLURM)

DORI uses the SLURM scheduler. Snakemake submits jobs to SLURM through the
`snakemake-executor-plugin-slurm` executor plugin (Snakemake 8+ moved cluster support to
plugins; this repo uses Snakemake 9).

This whole file is inherently deployment-specific -- "DORI" is this deployment's cluster, and
running elsewhere means different scheduler commands/module setup throughout. What *does*
travel with the repo is `config/config.yaml`'s `slurm:` block: account/QOS/partition/mail_user
are read from there (not hardcoded in any script), so pointing at a different allocation on the
same cluster -- or the same allocation on a different one -- is a config edit. See
[Portability](README.md#portability) in the README for the full config-vs-environment
breakdown.

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

## Quick option: one sbatch job running Snakemake locally
For light workloads (e.g. Phase 1's qc/parse stages, or Phase 2a's protein_properties -- 9 small
proteomes, all single-core/single-pass work) it's not worth the executor-plugin's per-rule job
submission overhead. Instead, `workflow/scripts/run_phase1.sbatch` (Phase 1) and
`run_phase2a.sbatch` (Phase 2a) each request one modest allocation and run
`snakemake --cores $SLURM_CPUS_PER_TASK --use-conda` locally inside it -- with no explicit
target, so each one also picks up any earlier stage that isn't done yet, not just its own. Submit
via the wrapper, not `sbatch` directly:
```bash
workflow/scripts/submit_phase1.sh
workflow/scripts/submit_phase2a.sh
```
Each wrapper reads `config/config.yaml`'s `slurm:` block (account/QOS/partition/mail_user) and
passes them as `sbatch` CLI flags -- so pointing this at a different allocation or cluster is a
config edit, not a script edit. (`sbatch workflow/scripts/run_phase1.sbatch` still works
directly, just without an account/QOS/mail-user unless your site has defaults for them.)

Logs land in `logs/phase1_<jobid>.{out,err}` / `logs/phase2a_<jobid>.{out,err}`. Switch to the
`--executor slurm` model above once a stage's per-genome runtime/memory actually warrants
separate jobs -- e.g. Phase 2b's planned intrinsic disorder prediction rule, which is heavy
enough (unlike everything in Phase 1/2a) to need its own per-genome SLURM resources rather than
running inline in one shared allocation. At that point, account/QOS/partition for *that* model
belong in a workflow profile (`config/slurm_profile/config.yaml`, see below), which is the
Snakemake-native equivalent of this wrapper for the per-rule submission model.

## Optional: a workflow profile
Instead of retyping flags every time, you can put them in `config/slurm_profile/config.yaml`
and run `snakemake --profile config/slurm_profile`. Not set up yet in this scaffold --
worth adding once per-rule resource requirements stabilize in Phase 1+.
