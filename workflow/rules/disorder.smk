rule disorder:
    """
    Stage 3 (Phase 2b) -- Intrinsic disorder prediction for one genome, via
    metapredict (a PyTorch model -- CPU inference here). Fundamentally
    different workload from Phase 2a's protein_properties: that rule is pure
    sequence/composition math (fast, single-core, safe on the login node);
    this one runs a neural net over every protein and took ~2h/genome
    single-threaded in earlier manual testing, ~30-40 min total across all 9
    genomes on a full 64-core node. That's why it's its own rule, with its
    own conda env (envs/disorder.yaml) and its own `threads`/`resources`
    below, instead of living inside protein_properties.

    Reads protein_properties.csv for the sequences only -- does not
    recompute or modify any Phase 2a column. Writes a separate
    results/disorder/{genome}/disorder.csv (genome + protein_id + 4 columns),
    joinable back onto protein_properties.csv on (genome, protein_id). See
    workflow/scripts/disorder.py's docstring for why this is a separate file
    rather than a merge into protein_properties.csv.

    `threads` here is Snakemake's mechanism for reserving cores from the
    pool given via `--cores` on the snakemake command line -- it's not a
    SLURM directive by itself. Setting it to 64 means: when this job runs,
    Snakemake treats all 64 cores (of whatever `--cores` value the sbatch
    script requests) as claimed, so no *other* job runs concurrently with
    it -- exactly what you want for something that's going to try to use
    the whole node via torch's CPU thread pool. It also flows into the
    shell command as {threads}, which is how the process actually finds out
    how many cores it's allowed to use.

    `resources: mem_mb` is a similar declaration, but for memory instead of
    cores -- it only *constrains* local scheduling if the snakemake
    invocation is given a matching `--resources mem_mb=N` budget (we don't
    do that here, since this pipeline runs snakemake locally inside a single
    exclusive-node sbatch allocation rather than through a cluster executor
    plugin/profile). Kept here anyway as an explicit, documented statement
    of what this rule actually needs -- useful if this ever moves to a
    profile that does enforce it.
    """
    input:
        OUTDIR + "/protein_properties/{genome}/protein_properties.csv",
    output:
        OUTDIR + "/disorder/{genome}/disorder.csv",
    conda:
        "../envs/disorder.yaml"
    threads: 64
    resources:
        mem_mb=64000,
        runtime=90,
    shell:
        "OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} "
        "python workflow/scripts/disorder.py --genome {wildcards.genome} "
        "--protein-properties {input} --output {output} --threads {threads}"
