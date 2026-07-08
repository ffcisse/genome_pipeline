rule qc:
    """
    Stage 1 -- Quality control on the raw input FASTA for one genome.

    SCAFFOLD ONLY: no real QC logic yet. Phase 1 will add real checks (sequence
    counts, empty/duplicate records, non-standard characters, etc.), same checks
    validated in the exploratory notebooks.

    `{genome}` is a wildcard: Snakemake fills it in with each value from GENOMES
    (defined in the Snakefile) when it needs to produce a qc output for that genome.
    One rule definition, one job per genome -- that's the point of wildcards.
    """
    input:
        protein=config["input"]["protein_dir"] + "/{genome}.fasta",
        cds=config["input"]["cds_dir"] + "/{genome}.fasta",
    output:
        OUTDIR + "/qc/{genome}.qc.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
