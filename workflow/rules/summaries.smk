rule summaries:
    """
    Stage 5 -- Aggregate per-genome protein + CDS properties into cross-genome
    summary tables (group/subgroup comparisons, effect sizes -- Phase 1).

    This is a FAN-IN rule: unlike qc/parse/protein_properties/cds_properties (one
    job per genome), this rule has no `{genome}` wildcard of its own -- it needs
    ALL genomes done first. `expand(...)` is how you say that: it takes a path
    template with a `{genome}` placeholder and GENOMES (the full list from the
    Snakefile), and returns one concrete path per genome, e.g.
    expand(OUTDIR + "/protein_properties/{genome}/protein_properties.csv", genome=GENOMES)
    becomes ["results/protein_properties/Cyamer1/protein_properties.csv",
    "results/protein_properties/Galsul1/protein_properties.csv", ...] -- a list of one file
    per genome this rule requires as input, so Snakemake won't run it until every
    genome has finished both protein_properties and cds_properties.

    cds= below points at cds_properties.csv (Phase 3's real output) rather
    than the old cds_properties/{genome}.properties.done marker from the
    Phase 0 stub -- updated to match once cds_properties.smk became a real
    rule, same as protein= already pointed at protein_properties.csv rather
    than a marker.
    """
    input:
        protein=expand(
            OUTDIR + "/protein_properties/{genome}/protein_properties.csv", genome=GENOMES
        ),
        cds=expand(OUTDIR + "/cds_properties/{genome}/cds_properties.csv", genome=GENOMES),
    output:
        OUTDIR + "/summaries/summaries.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
