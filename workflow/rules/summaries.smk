rule summaries:
    """
    Stage 5 -- Aggregate per-species protein + CDS properties into cross-species
    summary tables (lifestyle/lineage comparisons, effect sizes -- Phase 1).

    This is a FAN-IN rule: unlike qc/parse/protein_properties/cds_properties (one
    job per species), this rule has no `{species}` wildcard of its own -- it needs
    ALL species done first. `expand(...)` is how you say that: it takes a path
    template with a `{species}` placeholder and SPECIES (the full list from the
    Snakefile), and returns one concrete path per species, e.g.
    expand(OUTDIR + "/protein_properties/{species}.properties.done", species=SPECIES)
    becomes ["results/protein_properties/Cyamer1.properties.done",
    "results/protein_properties/Galsul1.properties.done", ...] -- a list of 9 files
    this rule requires as input, so Snakemake won't run it until every species has
    finished both protein_properties and cds_properties.
    """
    input:
        protein=expand(
            OUTDIR + "/protein_properties/{species}.properties.done", species=SPECIES
        ),
        cds=expand(OUTDIR + "/cds_properties/{species}.properties.done", species=SPECIES),
    output:
        OUTDIR + "/summaries/summaries.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
