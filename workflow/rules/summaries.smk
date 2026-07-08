rule summaries:
    """
    Stage 5 -- Aggregate per-genome protein + CDS properties into cross-genome
    summary tables (group/subgroup comparisons, effect sizes -- Phase 1).

    This is a FAN-IN rule: unlike qc/parse/protein_properties/cds_properties (one
    job per genome), this rule has no `{genome}` wildcard of its own -- it needs
    ALL genomes done first. `expand(...)` is how you say that: it takes a path
    template with a `{genome}` placeholder and GENOMES (the full list from the
    Snakefile), and returns one concrete path per genome, e.g.
    expand(OUTDIR + "/protein_properties/{genome}.properties.done", genome=GENOMES)
    becomes ["results/protein_properties/Cyamer1.properties.done",
    "results/protein_properties/Galsul1.properties.done", ...] -- a list of one file
    per genome this rule requires as input, so Snakemake won't run it until every
    genome has finished both protein_properties and cds_properties.
    """
    input:
        protein=expand(
            OUTDIR + "/protein_properties/{genome}.properties.done", genome=GENOMES
        ),
        cds=expand(OUTDIR + "/cds_properties/{genome}.properties.done", genome=GENOMES),
    output:
        OUTDIR + "/summaries/summaries.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
