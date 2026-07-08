rule qc:
    """
    Stage 1 -- Quality control on the raw input FASTA for one species.

    SCAFFOLD ONLY: no real QC logic yet. Phase 1 will add real checks (sequence
    counts, empty/duplicate records, non-standard characters, etc.), same checks
    validated in the exploratory notebooks.

    `{species}` is a wildcard: Snakemake fills it in with each value from SPECIES
    (defined in the Snakefile) when it needs to produce a qc output for that species.
    One rule definition, one job per species -- that's the point of wildcards.
    """
    input:
        protein=config["input"]["protein_dir"] + "/{species}.fasta",
        cds=config["input"]["cds_dir"] + "/{species}.fasta",
    output:
        OUTDIR + "/qc/{species}.qc.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
