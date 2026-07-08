rule protein_properties:
    """
    Stage 3 -- Compute per-protein physicochemical properties for one genome
    (pI, GRAVY, aliphatic index, thermostable fraction, etc. -- ported from the
    exploratory notebook in Phase 1).

    Depends on `parse` for the same genome.
    """
    input:
        OUTDIR + "/parsed/{genome}.parsed.done",
    output:
        OUTDIR + "/protein_properties/{genome}.properties.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
