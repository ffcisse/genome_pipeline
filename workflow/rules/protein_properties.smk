rule protein_properties:
    """
    Stage 3 -- Compute per-protein physicochemical properties for one species
    (pI, GRAVY, aliphatic index, thermostable fraction, etc. -- ported from the
    exploratory notebook in Phase 1).

    Depends on `parse` for the same species.
    """
    input:
        OUTDIR + "/parsed/{species}.parsed.done",
    output:
        OUTDIR + "/protein_properties/{species}.properties.done",
    conda:
        "../envs/analysis.yaml"
    shell:
        "touch {output}"
